"""
S1 ingest: FASTQ -> compacted Parquet tuple store  <work>/tuples/t*.parquet.

Reader processes (one per FASTQ role, see fastq.py) cut the files into chunks and deal them
round-robin to parse workers.  Each worker extracts the barcode and the target of every read
pair (see extract.py), aggregates (barcode, target) -> count (+ AND of the low-quality masks of
barcode and target) in a BOUNDED dict (--agg_entries) and flushes it as a small gz part.  A
background compactor (compaction.py) merges the parts into zstd Parquet, de-duplicating as it
goes.  The reject counters of the extraction are saved to <out>.meta/ingest_stats.json.
"""

import glob
import gzip
import json
import multiprocessing as mp
import os
import queue
import time
import traceback

from .common import dir_gb, log, log_sub, meta_path
from .compaction import Compactor
from .constants import (MASK32, MASK64, MAX_QUAL_BARCODE_LEN, MAX_TARGET_LEN, QUEUE_DEPTH,
                        TARGET_MASK_BITS)
from .extract import SequenceExtractor, max_length, reject_reason_order, spec_from_args
from .fastq import peek_read_len, reader_proc


# ---- parsing: one pair of FASTQ chunks -> aggregated (barcode, target) entries ---

def make_sequence_tables():
    """
    Byte translation tables for sequences.
      to_acgt : any base -> upper-case ACGT (everything else becomes 'A')
      nonacgt : '1' where the base is not ACGT/acgt, else '0'  (builds the low-quality mask)
    Parameters:
            -- None: no parameters
    Returns:
            -- tuple: (to_acgt, nonacgt), both bytes of length 256 for bytes.translate
    """
    # step1 : table to_acgt: any byte -> upper-case ACGT (non-ACGT becomes 'A')
    to_acgt = bytearray(b"A" * 256)
    for c in b"ACGT":
        to_acgt[c] = c
    for c in b"acgt":
        to_acgt[c] = c - 32
    # step2 : table nonacgt: '0' for ACGT / acgt, '1' for everything else (builds the low-quality mask)
    nonacgt = bytearray(b"1" * 256)
    for c in b"ACGTacgt":
        nonacgt[c] = 48
    return bytes(to_acgt), bytes(nonacgt)


def make_lowq_table(phred_offset, min_qual):
    """
    Byte translation table for quality strings: '1' where Phred < min_qual, else '0'.
    Parameters:
            -- phred_offset: int: ASCII offset of the quality encoding (33 = Sanger/Illumina 1.8+)
            -- min_qual: int: Phred values below this count as low quality
    Returns:
            -- bytes: table of length 256 for bytes.translate
    """
    table = bytearray(b"0" * 256)
    for c in range(256):
        if c - phred_offset < min_qual:
            table[c] = 49
    return bytes(table)


def ingest_chunk(barcode_chunk, target_chunk, cfg, agg, tables, extractors, rejects):
    """
    Parse one pair of FASTQ chunks into the aggregation dict.
    The barcode and the target are cut out of their reads by the extractors (whole read, fixed
    position or flanking sequences) and filtered (length, template).  A read
    pair is rejected when a check fails; only the first reason is counted in `rejects`, in the
    order target extraction, barcode extraction, target filters, barcode filters.
    Duplicate (barcode, target) pairs merge with count += 1 and mask &= new_mask, so a base stays
    'low quality' only if it was low in every duplicate.  A mask bit is 1 where the base is low
    quality or not ACGT; the first base of a sequence is the most significant bit.
    Parameters:
            -- barcode_chunk: bytes: FASTQ records of the barcode read (4 lines per record)
            -- target_chunk: bytes: FASTQ records of the target read, same record order
            -- cfg: dict: quality switches (built in stage_ingest)
            -- agg: dict: aggregation dict updated in place, key = b'barcode<TAB>target'
            -- tables: tuple: (to_acgt, nonacgt, target_lowq, barcode_lowq) translation tables
            -- extractors: tuple: (barcode SequenceExtractor, target SequenceExtractor)
            -- rejects: dict: reject key -> read pairs, updated in place
    Returns:
            -- int: number of read pairs in this chunk
    """
    # step1 : unpack the settings (quality switches, translation tables, extractors)
    use_target_qual = cfg["use_target_qual"]
    use_barcode_qual = cfg["use_barcode_qual"]
    use_mask = use_target_qual or use_barcode_qual
    to_acgt, nonacgt, target_lowq, barcode_lowq = tables
    span_barcode, check_barcode = extractors[0].span, extractors[0].check
    span_target, check_target = extractors[1].span, extractors[1].check

    # step2 : split both chunks into lines and pick the (upper-case) sequence / quality lines
    # FASTQ record = 4 lines: header, sequence, '+', quality.
    barcode_lines = barcode_chunk.split(b"\n")
    target_lines = target_chunk.split(b"\n")
    barcode_seqs = b"\n".join(barcode_lines[1::4]).upper().split(b"\n")
    barcode_quals = barcode_lines[3::4]
    target_seqs = b"\n".join(target_lines[1::4]).upper().split(b"\n")
    target_quals = target_lines[3::4]

    # step3 : loop over the read pairs
    get = agg.get
    for barcode_full, barcode_qual, target_full, target_qual in zip(
            barcode_seqs, barcode_quals, target_seqs, target_quals):
        #      |----   substep extract the target, then the barcode; count the first reject reason
        target_span = span_target(target_full)
        if target_span.__class__ is str:
            rejects[target_span] = rejects.get(target_span, 0) + 1
            continue
        barcode_span = span_barcode(barcode_full)
        if barcode_span.__class__ is str:
            rejects[barcode_span] = rejects.get(barcode_span, 0) + 1
            continue
        #      |----   substep filter the target, then the barcode (length, template)
        tg_start, tg_end = target_span
        bc_start, bc_end = barcode_span
        reason = (check_target(target_full, tg_start, tg_end)
                  or check_barcode(barcode_full, bc_start, bc_end))
        if reason is not None:
            rejects[reason] = rejects.get(reason, 0) + 1
            continue
        barcode = barcode_full[bc_start:bc_end]
        target = target_full[tg_start:tg_end]
        #      |----   substep quality-aware barcodes > 32 bp do not fit the mask (targets > 255 bp: see check)
        if use_barcode_qual and bc_end - bc_start > MAX_QUAL_BARCODE_LEN:
            raise RuntimeError("barcode longer than %d bp: quality-aware barcode clustering "
                               "is not possible; use --barcode_max_lowq 0" % MAX_QUAL_BARCODE_LEN)

        #      |----   substep aggregation key = barcode TAB target (non-ACGT target bases become 'A')
        key = barcode + b"\t" + target.translate(to_acgt)
        old = get(key)
        if not use_mask:
            agg[key] = 1 if old is None else old + 1
            continue

        #      |----   substep build the low-quality mask: target bits | (barcode bits << TARGET_MASK_BITS)
        # Low-quality mask of this read: quality below threshold OR base is not ACGT.
        mask = 0
        if use_target_qual:
            mask = (int(target_qual[tg_start:tg_end].translate(target_lowq), 2)
                    | int(target.translate(nonacgt), 2))
        if use_barcode_qual:
            mask |= (int(barcode_qual[bc_start:bc_end].translate(barcode_lowq), 2)
                     | int(barcode.translate(nonacgt), 2)) << TARGET_MASK_BITS
        #      |----   substep first time: store (mask << 32) | 1;  duplicate: AND the masks and add one to the count
        if old is None:
            agg[key] = (mask << 32) | 1
        else:                                               # AND the masks, add one to the count
            agg[key] = ((((old >> 32) & mask) << 32) | ((old & MASK32) + 1))
    # step4 : report the number of read pairs
    return len(barcode_seqs)


# ---- part files: bounded dictionary -> gz text part -------------------

def count_parts(parts_dir):
    """
    Number of text part files (*.tsv.gz) that wait for compaction (used for back-pressure).
    Parameters:
            -- parts_dir: str: directory of the part files
    Returns:
            -- int: number of finished part files
    """
    return sum(1 for x in os.listdir(parts_dir) if x.endswith(".tsv.gz"))


def flush_part(agg, parts_dir, worker_id, part_no, max_parts, abort):
    """
    Write the aggregation dict to one gz part file (8 columns) and empty it.
    Columns: barcode, target, cnt, tm0..tm3 (target low-quality mask as four 64-bit words) and
    bm (barcode low-quality mask).  Waits while too many parts are queued (back-pressure on the
    compactor, so the scratch disk cannot explode).  The file is published atomically.
    Parameters:
            -- agg: dict: key = b'barcode<TAB>target', value = (mask << 32) | count; emptied afterwards
            -- parts_dir: str: directory for the part files
            -- worker_id: int: parse worker number (part of the file name)
            -- part_no: int: running part number of this worker (part of the file name)
            -- max_parts: int: wait while this many parts are pending (0 = do not wait)
            -- abort: multiprocessing.Event: set by the supervisor on any error
    Returns:
            -- None: nothing (returns early without writing if abort is set)
    """
    # step1 : back-pressure: wait while too many parts are waiting for the compactor
    while max_parts and count_parts(parts_dir) >= max_parts:
        if abort.is_set():
            return
        time.sleep(2)
    # step2 : write the dictionary as a gz text part, 20000 lines per write
    #      |----   substep value = (mask << 32) | count  ->  count, target mask (four 64-bit words), barcode mask
    name = "%s/w%03d_%06d" % (parts_dir, worker_id, part_no)
    target_mask_all = (1 << TARGET_MASK_BITS) - 1
    with gzip.open(name + ".tmp", "wb", 1) as fo:
        lines = []
        for key, value in agg.items():
            mask = value >> 32
            target_mask = mask & target_mask_all
            barcode_mask = mask >> TARGET_MASK_BITS
            lines.append(b"%s\t%d\t%d\t%d\t%d\t%d\t%d" % (
                key, value & MASK32,
                target_mask & MASK64, (target_mask >> 64) & MASK64,
                (target_mask >> 128) & MASK64, (target_mask >> 192) & MASK64,
                barcode_mask))
            if len(lines) >= 20000:
                fo.write(b"\n".join(lines))
                fo.write(b"\n")
                lines = []
        if lines:
            fo.write(b"\n".join(lines))
            fo.write(b"\n")
    # step3 : publish the part atomically and empty the dictionary
    os.replace(name + ".tmp", name + ".tsv.gz")             # atomic: compactor never sees a half file
    agg.clear()


# ---- parse worker process --------------------------------------------------

def ingest_worker(worker_id, barcode_queue, target_queue, result_queue, progress, cfg, abort):
    """
    Parse worker: pulls chunk pairs (one barcode chunk and one target chunk), aggregates them,
    flushes a part file whenever the dict reaches cfg['cap'] entries, and reports
    (n_reads, rejects) when both queues deliver the end marker.
    Parameters:
            -- worker_id: int: worker number
            -- barcode_queue: multiprocessing queue: chunks of the barcode read
            -- target_queue: multiprocessing queue: chunks of the target read
            -- result_queue: multiprocessing queue: receives ('done', id, n_reads, rejects) or
               ('error', ...)
            -- progress: multiprocessing.Array: reads parsed so far per worker (for the progress log)
            -- cfg: dict: settings built in stage_ingest
            -- abort: multiprocessing.Event: set by the supervisor on any error
    Returns:
            -- None: nothing; results go through result_queue
    """
    try:
        # step1 : set up the extractors, translation tables, the bounded dictionary and the counters
        extractors = (SequenceExtractor(cfg["barcode_spec"]), SequenceExtractor(cfg["target_spec"]))
        tables = make_sequence_tables() + (
            make_lowq_table(cfg["phred_offset"], cfg["target_min_qual"]),
            make_lowq_table(cfg["phred_offset"], cfg["barcode_min_qual"]))
        agg = {}
        rejects = {}
        cap, parts_dir, max_parts = cfg["cap"], cfg["parts_dir"], cfg["max_parts"]
        n_read = part_no = 0
        # step2 : consume chunk pairs until both queues deliver the end marker
        while True:
            #      |----   substep take one chunk from each role (they must stay in lock-step)
            barcode_item = barcode_queue.get()
            target_item = target_queue.get()
            #      |----   substep end marker: both roles must end together
            if barcode_item is None or target_item is None:
                if not (barcode_item is None and target_item is None):
                    raise RuntimeError("barcode/target read files differ in read number")
                break
            if barcode_item[0] != target_item[0]:
                raise RuntimeError("barcode/target read files differ in read number")
            #      |----   substep parse the chunk pair and aggregate it
            n_read += ingest_chunk(barcode_item[1], target_item[1], cfg, agg, tables,
                                   extractors, rejects)
            del barcode_item, target_item
            progress[worker_id] = n_read
            #      |----   substep write a part file when the bounded dictionary is full
            if len(agg) >= cap:
                flush_part(agg, parts_dir, worker_id, part_no, max_parts, abort)
                part_no += 1
        # step3 : flush the remainder and report the totals to the supervisor
        if agg:
            flush_part(agg, parts_dir, worker_id, part_no, 0, abort)
        result_queue.put(("done", worker_id, n_read, rejects))
    except BaseException:
        result_queue.put(("error", "worker %d" % worker_id, traceback.format_exc()))


# ---- stage helpers ---------------------------------------------------------

def decide_quality(args, barcode_spec, read_len):
    """
    Decide whether barcode / target qualities are used, from the longest barcode the extraction
    can return (barcode masks must fit in one 64-bit word).
    Parameters:
            -- args: argparse.Namespace: barcode_max_lowq, target_max_lowq
            -- barcode_spec: ExtractSpec: barcode extraction settings
            -- read_len: int: length of the first barcode read
    Returns:
            -- tuple: (longest barcode in bp, use_barcode_qual, use_target_qual)
    """
    barcode_len_max = max_length(barcode_spec, read_len)
    #      |----   substep warn for very long barcodes (the tuple store then approaches one row per read)
    if barcode_len_max > 64:
        log_sub("NOTE: barcode up to %d bp (read length %d). Whole-read / variant barcodes are "
                "fine, but if most reads are distinct the tuple store approaches one row per "
                "read (disk!)" % (barcode_len_max, read_len))
    #      |----   substep barcode quality needs barcodes <= 32 bp (one 64-bit word) and --barcode_max_lowq > 0
    use_barcode_qual = args.barcode_max_lowq > 0 and barcode_len_max <= MAX_QUAL_BARCODE_LEN
    if args.barcode_max_lowq > 0 and not use_barcode_qual:
        log_sub("NOTE: barcode can be %d bp (> %d): barcode quality is NOT used (sequence-only)"
                % (barcode_len_max, MAX_QUAL_BARCODE_LEN))
    #      |----   substep target quality is on whenever --target_max_lowq > 0
    return barcode_len_max, use_barcode_qual, args.target_max_lowq > 0


def describe_spec(spec):
    """
    One-line description of an extraction setting (for the log).
    Parameters:
            -- spec: ExtractSpec: settings of one role
    Returns:
            -- str: e.g. 'flank: GTCC...ATCT / TGGT...ACCA (<= 2 mismatches)'
    """
    if spec.method == "pos":
        text = "pos: start %d, length %s" % (spec.start, spec.length or "to read end")
    elif spec.method == "flank":
        text = "flank: %s / %s (<= %d mismatches, %s)" % (spec.up, spec.down,
                                                           spec.flank_mismatch,
                                                           spec.flank_distance)
    else:
        text = "none: whole read"
    if spec.min_len:
        text += ", length >= %d" % spec.min_len
    if spec.allowed_lens:
        text += ", lengths %d-%d only" % (min(spec.allowed_lens), max(spec.allowed_lens))
    if spec.template:
        text += ", template %s (<= %d mismatches)" % (spec.template, spec.template_mismatch)
    return text


def write_ingest_stats(args, specs, n_reads, rejects):
    """
    Save the read-pair counts of the extraction to <out>.meta/ingest_stats.json (read by S5).
    Parameters:
            -- args: argparse.Namespace: uses args.out
            -- specs: tuple: (barcode ExtractSpec, target ExtractSpec)
            -- n_reads: int: read pairs parsed
            -- rejects: dict: reject key -> read pairs
    Returns:
            -- None: nothing
    """
    rows = [(key, rejects.get(key, 0), description)
            for key, description in reject_reason_order(*specs)]
    os.makedirs(meta_path(args), exist_ok=True)
    with open(meta_path(args, "ingest_stats.json"), "w") as fo:
        json.dump({"total_read_pairs": n_reads,
                   "extracted_read_pairs": n_reads - sum(rejects.values()),
                   "rejects": rows}, fo, indent=1)


# ---- stage driver ----------------------------------------------------------

def stage_ingest(args, work, threads, mem_gb):
    """
    S1: FASTQ -> compacted Parquet tuple store  <work>/tuples/t*.parquet.
    Starts one reader process per FASTQ role and the parse workers, supervises them, drives the
    background compaction and writes <out>.meta/ingest_stats.json and <work>/ingest.done.
    Uses args: barcode_read, target_read, method, flank_mismatch, flank_distance,
    barcode_* / target_* options, phred_offset, chunk, agg_entries, compact_parts, merge_files,
    max_reads, no_pigz.
    Parameters:
            -- args: argparse.Namespace: parsed command line
            -- work: str: work directory
            -- threads: int: total processes available
            -- mem_gb: float: RAM budget in GB
    Returns:
            -- int: number of read pairs parsed
    """
    ctx = mp.get_context("spawn")
    use_pigz = not args.no_pigz

    # step1 : extraction settings; decide whether barcode / target qualities are used
    #      |----   substep pos method: the barcode must be complete, a short target is kept as it is;
    #      |----   substep targets > MAX_TARGET_LEN do not fit the quality mask and are rejected
    barcode_spec = spec_from_args(args, "barcode", require_full_len=True)
    target_spec = spec_from_args(args, "target", require_full_len=False, max_len=MAX_TARGET_LEN)
    #      |----   substep read the first record of the barcode FASTQ to learn the read length
    read_len = peek_read_len(args.barcode_read[0], use_pigz)
    barcode_len_max, use_barcode_qual, use_target_qual = decide_quality(args, barcode_spec,
                                                                        read_len)

    # step2 : decide the process layout, scratch directories and the per-worker memory cap
    #      |----   substep 2 reader processes (one per FASTQ role); the rest parse
    n_workers = threads - 2 if threads > 3 else threads
    #      |----   substep scratch directories: text parts and compacted Parquet tuples
    parts_dir = work + "/parts"
    tuples_dir = work + "/tuples"
    os.makedirs(parts_dir, exist_ok=True)
    os.makedirs(tuples_dir, exist_ok=True)

    #      |----   substep max aggregated entries per worker: ~half the RAM budget / ~450 B per entry
    cap = args.agg_entries or int(min(4e6, max(5e4, mem_gb * 1e9 * 0.5 / n_workers / 450.0)))
    #      |----   substep settings shipped to every parse worker
    cfg = dict(barcode_spec=barcode_spec, target_spec=target_spec,
               use_target_qual=use_target_qual, use_barcode_qual=use_barcode_qual,
               phred_offset=args.phred_offset,
               target_min_qual=args.target_min_qual, barcode_min_qual=args.barcode_min_qual,
               cap=cap, parts_dir=parts_dir, max_parts=max(args.compact_parts * 4, 8))
    log("S1: ingest  (%d workers, <=%d aggregated entries/worker, barcode <= %d bp, "
        "barcode quality %s, target quality %s)"
        % (n_workers, cap, barcode_len_max, "on" if use_barcode_qual else "off",
           "on" if use_target_qual else "off"))
    log_sub("barcode: %s" % describe_spec(barcode_spec))
    log_sub("target : %s" % describe_spec(target_spec))

    # step3 : start the reader / parse-worker processes and the background compactor
    #      |----   substep abort flag, result queue and bounded chunk queues (QUEUE_DEPTH per worker and role)
    abort = ctx.Event()
    result_queue = ctx.Queue()
    barcode_queues = [ctx.Queue(maxsize=QUEUE_DEPTH) for _ in range(n_workers)]
    target_queues = [ctx.Queue(maxsize=QUEUE_DEPTH) for _ in range(n_workers)]
    progress = ctx.Array("q", n_workers, lock=False)
    #      |----   substep one parse worker per queue pair, one reader per FASTQ role
    workers = [ctx.Process(target=ingest_worker,
                           args=(w, barcode_queues[w], target_queues[w], result_queue,
                                 progress, cfg, abort))
               for w in range(n_workers)]
    readers = [ctx.Process(target=reader_proc,
                           args=(paths, args.chunk, use_pigz, args.max_reads, queues,
                                 result_queue, abort))
               for paths, queues in ((args.barcode_read, barcode_queues),
                                     (args.target_read, target_queues))]
    processes = workers + readers
    for p in processes:
        p.daemon = True
        p.start()
    #      |----   substep background Parquet compactor, driven from the supervisor loop
    compactor = Compactor(work, max(1, threads // 4), mem_gb * 0.25,
                          args.compact_parts, args.merge_files)

    #      |----   substep helper: stop every process and exit (used on any error)
    def die(message):
        """
        Stop every reader / worker process and exit with an error message.
        Parameters:
                -- message: str: error text for SystemExit
        Returns:
                -- None: never returns (raises SystemExit)
        """
        abort.set()
        for p in processes:
            if p.is_alive():
                p.terminate()
        raise SystemExit(message)

    # step4 : supervise: collect results, drive the compaction, watch for crashes, log progress
    n_done = n_reads = 0
    rejects = {}
    ingest_start = time.time()
    last_log = time.time()
    while n_done < n_workers:
        try:
            #      |----   substep poll the result queue (3 s timeout); an 'error' message stops everything
            message = result_queue.get(timeout=3)
        except queue.Empty:
            message = None
        if message is not None:
            if message[0] == "error":
                die("ERROR in %s:\n%s" % (message[1], message[2]))
            n_done += 1
            n_reads += message[2]
            for key, n in message[3].items():
                rejects[key] = rejects.get(key, 0) + n
        try:
            #      |----   substep let the compactor start a new job when it is idle
            compactor.step()
        except SystemExit as e:
            die(str(e))
        for p in processes:
            #      |----   substep stop when a worker / reader died (usually out of memory)
            if p.exitcode not in (None, 0):
                die("ERROR: a worker/reader process died (exit code %s) - out of memory? "
                    "lower --threads / --agg_entries / --chunk" % p.exitcode)
        #      |----   substep progress line once a minute
        if time.time() - last_log > 60:
            last_log = time.time()
            total = sum(progress)
            log_sub("%d reads parsed (%.1f M/s)"
                    % (total, total / 1e6 / max(1.0, time.time() - ingest_start)))
    # step5 : wait for the processes, finish the compaction and close the stage
    for p in processes:
        p.join(60)
    #      |----   substep merge all remaining parts / Parquet files
    try:
        compactor.finish()
    except SystemExit as e:
        die(str(e))
    #      |----   substep report the reject reasons and save them for the stats file (S5)
    n_rejected = sum(rejects.values())
    log_sub("ingest done: %d read pairs, %d rejected; tuple store %.2f GB"
            % (n_reads, n_rejected, dir_gb(tuples_dir, "t*.parquet")))
    for key, n in sorted(rejects.items(), key=lambda kv: -kv[1]):
        log_sub("  rejected %-30s %d" % (key, n))
    write_ingest_stats(args, (barcode_spec, target_spec), n_reads, rejects)
    #      |----   substep sanity check: the tuple store must not be empty
    if not glob.glob(tuples_dir + "/t*.parquet"):
        raise SystemExit("ERROR: no read pair passed the barcode / target extraction "
                         "(see the reject counts above)")
    #      |----   substep marker file used by --resume
    with open(work + "/ingest.done", "w") as fo:
        fo.write("%d\n" % n_reads)
    return n_reads
