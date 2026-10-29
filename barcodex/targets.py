"""
S4 targets: cluster the target variants inside every barcode cluster (one process per shard)
and write <out>.meta/target_pairs.tsv.
"""

import glob
import multiprocessing as mp
import os
import shutil

from .common import log, log_sub, meta_path, open_duckdb, sql_str
from .constants import COST_SCALE
from .distance import BASE_TO_DIGIT, lane_mask, mate_cost, popcount, unpack_seq

# ---- clustering the target variants of ONE barcode --------------------

# Per-process state of the S4 worker pool (filled by init_s4_worker).
S4_STATE = {}

def _align_cost(packed_a, low_pos_a, len_a, packed_b, low_pos_b, len_b, params):
    """
    (hard, low) cost of two variants with a banded alignment (indels allowed).  Unpacks the
    sequences, rebuilds the per-base quality arrays from the low-quality positions and calls
    mate_cost.
    Parameters:
            -- packed_a: bytes: packed sequence A
            -- low_pos_a: bytes: low-quality positions of A (one byte per position)
            -- len_a: int: length of A
            -- packed_b: bytes: packed sequence B
            -- low_pos_b: bytes: low-quality positions of B
            -- len_b: int: length of B
            -- params: tuple: (max_hard, max_low, band, ratio, low_ratio, use_quality)
    Returns:
            -- tuple or None: (hard, low), or None if over the limits
    """
    # step1 : unpack the settings
    max_hard, max_low, band, _ratio, _low_ratio, use_quality = params
    # step2 : unpack both sequences to ACGT bytes
    seq_a = unpack_seq(packed_a, len_a)
    seq_b = unpack_seq(packed_b, len_b)
    # step3 : per-base quality arrays rebuilt from the low-quality positions (0xff = confident, 0 = low)
    qual_a = bytearray(b"\xff" * len_a)         # 0xff = confident, 0 = low quality
    qual_b = bytearray(b"\xff" * len_b)
    if use_quality:
        for pos in low_pos_a:
            qual_a[pos] = 0
        for pos in low_pos_b:
            qual_b[pos] = 0
    # step4 : banded alignment cost with the S4 limits
    return mate_cost(seq_a, bytes(qual_a), seq_b, bytes(qual_b),
                     max_hard, max_low, band, 1 if use_quality else 0)


def cluster_group(items, params):
    """
    Cluster the target variants of ONE barcode.
    Variants are sorted by abundance.  A variant X is absorbed into a more abundant parent P when
      * the edit distance is within limits (hard <= max_hard, low <= max_low), and
      * parent_count >= ratio x child_count (low_ratio instead if only low-quality positions differ).
    Among valid parents the one with the smallest (hard, low) is taken.  Chains (X -> P -> Q) end
    at their root, which becomes the cluster consensus.
    Parameters:
            -- items: list: [(count, key = packed_seq + length_byte, low_quality_positions_bytes)]
            -- params: tuple: (max_hard, max_low, band, ratio, low_ratio, use_quality)
    Returns:
            -- list: [(reads, n_unique_variants, consensus_seq, first_cov, second_cov)], one per target
               cluster;
               first_cov / second_cov = reads of the most / second most abundant variant (0 if only one)
    """
    # step1 : unpack the settings and sort the variants by abundance (most abundant first)
    max_hard, max_low, band, ratio, low_ratio, use_quality = params
    max_parents = S4_STATE["max_parents"]
    #      |----   substep sort by count (descending), then by sequence key
    items.sort(key=lambda t: (-t[0], t[1]))
    n = len(items)
    #      |----   substep very large group: no indel search
    if n > S4_STATE["max_dp_group"]:        # huge group: Hamming-only, indel DP is too slow
        band = 0

    # step2 : prepare every variant
    #      |----   substep sequence as an integer (2 bits/base) and low-quality positions as a lane mask,
    #      |----   substep so that a comparison is a few integer operations
    prepared = []
    for count, key, low_positions in items:
        length = key[-1]
        packed = key[:-1]
        padded_len = (length + 3) & ~3
        low_bits = 0
        if use_quality:
            for pos in low_positions:
                low_bits |= 1 << (2 * (padded_len - 1 - pos))
        prepared.append((int.from_bytes(packed, "big"), low_bits, length, padded_len,
                         packed, low_positions))

    # step3 : find the best parent of every variant (parents always have a smaller index)
    parent_of = [-1] * n
    min_ratio = min(ratio, low_ratio)
    for x in range(1, n):                   # variant 0 is the most abundant: always a root
        child_count = items[x][0]
        threshold = min_ratio * child_count
        child_int, child_low, child_len, padded_len, child_packed, child_low_pos = prepared[x]
        best_cost = best_parent = -1
        for p in range(x if x < max_parents else max_parents):   # only the most abundant can be parents
            parent_count = items[p][0]
            if parent_count < threshold:
                break                       # sorted by count: no later parent can qualify
            parent_int, parent_low, parent_len, _, parent_packed, parent_low_pos = prepared[p]
            result = None
            #      |----   substep same length: Hamming on 2-bit lanes; 'low' = mismatches at low-quality bases
            if child_len == parent_len:
                diff = child_int ^ parent_int
                lanes = (diff | (diff >> 1)) & lane_mask(padded_len)     # 1 per mismatching base
                total = popcount(lanes)
                low = popcount(lanes & (child_low | parent_low)) if use_quality else 0
                hard = total - low
                if hard <= max_hard and low <= max_low:
                    result = (hard, low)
                elif band and max_hard:     # Hamming failed: try an alignment with an indel
                    result = _align_cost(child_packed, child_low_pos, child_len,
                                         parent_packed, parent_low_pos, parent_len, params)
            #      |----   substep different length: banded alignment (indel)
            elif band and abs(child_len - parent_len) <= band:
                result = _align_cost(child_packed, child_low_pos, child_len,
                                     parent_packed, parent_low_pos, parent_len, params)
            if result is None:
                continue
            #      |----   substep coverage rule: parent_count >= ratio x child_count (low_ratio when only low-quality bases differ)
            needed_ratio = low_ratio if result[0] == 0 else ratio
            if parent_count < needed_ratio * child_count:
                continue
            #      |----   substep keep the cheapest valid parent (hard edits first, then low)
            cost = result[0] * COST_SCALE + result[1]
            if best_parent < 0 or cost < best_cost:
                best_cost, best_parent = cost, p
        parent_of[x] = best_parent

    # step4 : resolve chains X -> P -> Q to their root (one forward pass is enough)
    root = list(range(n))
    for x in range(n):
        if parent_of[x] >= 0:
            root[x] = root[parent_of[x]]
    # step5 : per cluster: total reads, number of variants and the top-two variant coverages
    #      |----   substep `items` is sorted by count (descending) and the root has the smallest index of its cluster,
    #      |----   substep so the first member seen is the top variant and the second member seen is the second one
    read_sum, unique_count, first_cov, second_cov = {}, {}, {}, {}
    for x in range(n):
        r = root[x]
        count = items[x][0]
        if r not in read_sum:
            read_sum[r], unique_count[r], first_cov[r], second_cov[r] = count, 1, count, 0
        else:
            read_sum[r] += count
            unique_count[r] += 1
            if unique_count[r] == 2:
                second_cov[r] = count
    # step6 : one result tuple per cluster root (consensus = root sequence)
    return [(read_sum[r], unique_count[r], unpack_seq(prepared[r][4], prepared[r][2]),
             first_cov[r], second_cov[r])
            for r in sorted(read_sum)]


# ---- processing one shard ---------------------------------------------

def run_groups(rows, out_clusters, out_barcode_fq, out_target_fq, params, long_group_note=5000):
    """
    Cluster targets barcode by barcode and write the cluster rows.
    Rows must arrive ORDERED BY centre, so one barcode's rows are adjacent and only one barcode is
    held in memory at a time.
    Parameters:
            -- rows: iterable of (centre, target_str, cnt, tm0, tm1, tm2, tm3) ordered by centre
            -- out_clusters: text file: receives one line per target cluster
            -- out_barcode_fq: text file or None: representative barcode FASTQ
            -- out_target_fq: text file or None: representative target FASTQ
            -- params: tuple: (max_hard, max_low, band, ratio, low_ratio, use_quality)
            -- long_group_note: int: log a note for barcodes with more distinct targets than this
    Returns:
            -- tuple: (n_barcodes, n_clusters, n_reads)
    """
    # step1 : counters and the per-barcode emit() helper
    state = {"n_barcodes": 0, "n_clusters": 0, "n_reads": 0}

    def emit(barcode, items):
        """
        Cluster one barcode's targets and write its result lines (and the optional FASTQ records).
        Parameters:
                -- barcode: str: centre barcode
                -- items: list: [(count, key, low_quality_positions)] of this barcode
        Returns:
                -- None: nothing; updates the counters in `state`
        """
        total_reads = sum(t[0] for t in items)
        if len(items) > long_group_note:
            log_sub("NOTE: barcode %s has %d reads and %d distinct targets (junk/low-complexity "
                "barcode?); top-%d parents only, no indel search"
                % (barcode, total_reads, len(items), S4_STATE["max_parents"]))
        #      |----   substep cluster the targets of this barcode
        clusters = cluster_group(items, params)
        clusters.sort(key=lambda t: -t[0])
        #      |----   substep one output line per target cluster (and the optional FASTQ pair)
        for rank, (reads, n_unique, seq, first_cov, second_cov) in enumerate(clusters, 1):
            out_clusters.write("%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\n"
                               % (barcode, total_reads, rank, reads, n_unique, seq.decode(),
                                  first_cov, second_cov))
            if out_barcode_fq:
                header = "@%s_%d reads=%d" % (barcode, rank, reads)
                out_barcode_fq.write("%s/1\n%s\n+\n%s\n" % (header, barcode, "I" * len(barcode)))
                out_target_fq.write("%s/2\n%s\n+\n%s\n" % (header, seq.decode(), "I" * len(seq)))
        state["n_barcodes"] += 1
        state["n_clusters"] += len(clusters)
        state["n_reads"] += total_reads

    # step2 : stream the rows (ordered by centre) and collect one barcode at a time
    current = None
    items = []
    for centre, target, count, tm0, tm1, tm2, tm3 in rows:
        #      |----   substep a new centre starts: cluster and write the previous barcode first
        if centre != current:                       # rows are sorted: new barcode starts here
            if current is not None:
                emit(current, items)
            current = centre
            items = []
        #      |----   substep pack the target: 2 bits per base + one length byte
        target_bytes = target.encode()
        length = len(target_bytes)
        # sequence -> integer digits 0..3 -> packed bytes, plus the length as last byte
        digits = target_bytes.translate(BASE_TO_DIGIT)
        pad = (-length) & 3
        if pad:
            digits += b"0" * pad
        key = int(digits, 4).to_bytes((length + 3) >> 2, "big") + bytes((length,))
        #      |----   substep low-quality mask (4 x 64 bit words) -> list of low-quality positions
        mask = tm0 | (tm1 << 64) | (tm2 << 128) | (tm3 << 192)
        positions = bytearray()
        while mask:
            lowest = mask & -mask
            positions.append(length - lowest.bit_length())
            mask ^= lowest
        items.append((count, key, bytes(positions)))
    # step3 : write the last barcode and return the totals
    if current is not None:
        emit(current, items)
    return state["n_barcodes"], state["n_clusters"], state["n_reads"]


def init_s4_worker(cfg):
    """
    Pool initializer: store the settings in the global S4_STATE of the worker process.
    Parameters:
            -- cfg: dict: settings built in stage_targets
    Returns:
            -- None: nothing
    """
    S4_STATE.update(cfg)


def shard_rows(shard):
    """
    Stream the rows of one shard ordered by barcode centre (DuckDB, bounded memory).
    Parameters:
            -- shard: int: shard number
    Returns:
            -- generator: yields (centre, target, cnt, tm0, tm1, tm2, tm3) tuples
    """
    # step1 : locate the shard and open a one-thread DuckDB
    state = S4_STATE
    shard_dir = "%s/shards/shard=%d" % (state["work"], shard)
    con = open_duckdb(state["work"], 1, state["s4_mem_gb"], "s4_%d" % shard)
    # step2 : ordered scan: all rows of one barcode centre are adjacent
    cursor = con.execute("SELECT centre, target, cnt, tm0, tm1, tm2, tm3 FROM read_parquet(%s) "
                         "ORDER BY centre" % sql_str(shard_dir + "/*.parquet"))
    # step3 : stream the rows in blocks of 50000 (bounded memory)
    while True:
        rows = cursor.fetchmany(50000)
        if not rows:
            break
        for row in rows:
            yield row
    con.close()


def process_shard(shard):
    """
    Cluster all barcodes of one shard; write res<shard>.tsv (and FASTQs).
    Skips shards that are already finished (resume) or empty.
    Parameters:
            -- shard: int: shard number
    Returns:
            -- tuple: (n_barcodes, n_clusters, n_reads); (0, 0, 0) for a skipped shard
    """
    # step1 : locate the shard; skip it when finished (resume) or empty
    state = S4_STATE
    work = state["work"]
    final = "%s/res%05d.tsv" % (work, shard)
    shard_dir = "%s/shards/shard=%d" % (work, shard)
    if os.path.exists(final) or not glob.glob(shard_dir + "/*.parquet"):
        return 0, 0, 0
    # step2 : open the output files (.tmp files are published at the end)
    out_clusters = open(final + ".tmp", "w")
    out_barcode_fq = out_target_fq = None
    if state["fastq"]:
        out_barcode_fq = open("%s/res%05d.bc.fq" % (work, shard), "w")
        out_target_fq = open("%s/res%05d.tg.fq" % (work, shard), "w")
    # step3 : cluster all barcodes of the shard
    result = run_groups(shard_rows(shard), out_clusters, out_barcode_fq, out_target_fq,
                        state["params"])
    # step4 : close the files and publish the result atomically
    out_clusters.close()
    if out_barcode_fq:
        out_barcode_fq.close()
        out_target_fq.close()
    os.replace(final + ".tmp", final)
    # step5 : free disk space as we go
    if not state["keep_tmp"]:
        shutil.rmtree(shard_dir, ignore_errors=True)        # free disk as we go
    return result


# ---- stage driver -----------------------------------------------------

def run_pool(ctx, threads, func, tasks, initializer=None, initargs=()):
    """
    Yield func(task) for every task, in parallel when threads > 1 (results in completion order).
    Parameters:
            -- ctx: multiprocessing context (spawn)
            -- threads: int: processes; 1 = run inline in this process
            -- func: callable: top-level function applied to every task
            -- tasks: list: task arguments
            -- initializer: callable or None: run once in every worker
            -- initargs: tuple: arguments of the initializer
    Returns:
            -- generator: yields the return value of func for every task
    """
    # step1 : single process: run inline (initializer first)
    if threads <= 1:
        if initializer:
            initializer(*initargs)
        for task in tasks:
            yield func(task)
        return
    # step2 : process pool: imap_unordered, results arrive as tasks finish
    with ctx.Pool(threads, initializer=initializer, initargs=initargs) as pool:
        for result in pool.imap_unordered(func, tasks):
            yield result


def stage_targets(args, work, threads, mem_gb, n_shards):
    """
    S4: cluster the targets of every barcode (one process per shard) and concatenate the shard
    results into <out>.meta/target_pairs.tsv (the input of S5).  With --fastq it also writes
    <out>_barcode.fastq and <out>_target.fastq.
    Parameters:
            -- args: argparse.Namespace: target_dist/max_lowq/max_indel/ratio/lowq_ratio, max_parents,
               max_dp_group, fastq, keep_tmp, out
            -- work: str: work directory
            -- threads: int: worker processes
            -- mem_gb: float: RAM budget in GB (70% is shared by the shard workers)
            -- n_shards: int: number of shards
    Returns:
            -- None: nothing; results are files
    """
    log("S4: clustering targets within barcodes ...")
    # step1 : settings that every worker process receives (S4_STATE)
    ctx = mp.get_context("spawn")
    use_quality = args.target_max_lowq > 0
    cfg = dict(work=work, fastq=args.fastq, keep_tmp=args.keep_tmp,
               max_parents=max(1, args.max_parents), max_dp_group=args.max_dp_group,
               s4_mem_gb=max(0.5, mem_gb * 0.7 / threads),
               params=(args.target_dist, args.target_max_lowq,
                       min(args.target_max_indel, args.target_dist),
                       max(args.target_ratio, 1.0), max(args.target_lowq_ratio, 1.0),
                       use_quality))
    # step2 : cluster the shards in parallel (one shard per task) and add up the counters
    n_barcodes = n_clusters = n_reads = n_done = 0
    for barcodes_i, clusters_i, reads_i in run_pool(ctx, threads, process_shard,
                                                    list(range(n_shards)),
                                                    init_s4_worker, (cfg,)):
        n_barcodes += barcodes_i
        n_clusters += clusters_i
        n_reads += reads_i
        n_done += 1
        if n_done % max(1, n_shards // 10) == 0:
            log_sub("%d/%d shards" % (n_done, n_shards))

    # step3 : concatenate the shard files in shard order
    def concatenate(name_format, dest, header=None):
        """
        Append the per-shard files (in shard order) to dest, deleting them.
        Parameters:
                -- name_format: str: '%s/res%05d.<ext>' pattern filled with (work, shard)
                -- dest: str: output file
                -- header: str or None: header line written first
        Returns:
                -- None: nothing
        """
        with open(dest, "wb") as fo:
            if header:
                fo.write(header.encode())
            for shard in range(n_shards):
                path = name_format % (work, shard)
                if os.path.exists(path):
                    with open(path, "rb") as fi:
                        shutil.copyfileobj(fi, fo, 1 << 22)
                    os.remove(path)

    #      |----   substep meta/target_pairs.tsv = input of S5
    os.makedirs(meta_path(args), exist_ok=True)
    concatenate("%s/res%05d.tsv", meta_path(args, "target_pairs.tsv"),
                "barcode\tbarcode_reads\ttarget_rank\ttarget_reads\ttarget_unique\ttarget_seq\t"
                "first_cov\tsecond_cov\n")
    #      |----   substep optional representative FASTQ pair
    if args.fastq:
        concatenate("%s/res%05d.bc.fq", args.out + "_barcode.fastq")
        concatenate("%s/res%05d.tg.fq", args.out + "_target.fastq")
    # step4 : report the totals
    log_sub("S4 done: %d reads, %d barcode clusters with reads, %d (barcode, target) clusters"
        % (n_reads, n_barcodes, n_clusters))
