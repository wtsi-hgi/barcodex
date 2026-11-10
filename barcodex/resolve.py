"""
S5 resolve: one target per barcode, then the result files (DuckDB).

Inputs  (intermediate files in <out>.meta/, written by S1 / S2 / S4):
  target_pairs.tsv           (barcode cluster, target cluster) rows
  barcode_clusters.parquet   per barcode-cluster statistics
  raw_stats.json             raw barcode / target counts, reads kept by the barcode clustering
  ingest_stats.json          read pairs parsed / rejected by the extraction
Outputs:
  <out>.barcode_target.details.tsv   one row per final (barcode, target) pair
  <out>.stats.tsv                    summary numbers (see stats.py)
--resolve_only re-runs this stage from <out>.meta/.

Memory: the tables hold no sequences.  Barcodes and targets are keyed by md5_number() (128 bit,
UHUGEINT); the sequences are read again from target_pairs.tsv only while the result files are
written.  Only --merge_dist > 0 keeps a table of the distinct target sequences (table tstr).
"""

import json
import os
import shutil

from .common import log, log_sub, meta_path, sql_str
from .stats import collect_stats, write_stats

PAIR_COLUMNS = ("{'barcode':'VARCHAR', 'barcode_reads':'BIGINT', 'target_rank':'INTEGER', "
                "'target_reads':'BIGINT', 'target_unique':'INTEGER', 'target_seq':'VARCHAR', "
                "'first_cov':'BIGINT', 'second_cov':'BIGINT'}")


def pairs_csv(src):
    """
    SQL table expression that reads target_pairs.tsv (S4 output; no quoting, so no sniffing).
    Parameters:
            -- src: str: path of target_pairs.tsv
    Returns:
            -- str: read_csv(...) expression
    """
    return ("read_csv(%s, delim='\t', header=true, quote='', escape='', columns=%s)"
            % (sql_str(src), PAIR_COLUMNS))


# ---- steps -----------------------------------------------------------------

def merge_targets(con, args, src):
    """
    Optional (--merge_dist > 0): merge near-identical targets ACROSS barcodes.  Creates table tstr
    (target key -> sequence) and rewrites table tmap (target key -> root key) so that every child
    points to its most abundant root.
    Parameters:
            -- con: duckdb connection with table tseq (target, reads, n_barcodes)
            -- args: argparse.Namespace: merge_dist, merge_ratio
            -- src: str: path of target_pairs.tsv
    Returns:
            -- None: nothing; creates tstr and replaces table tmap
    """
    # step1 : distinct target sequences (the only table with sequences)
    merge_dist = args.merge_dist
    n_seeds = merge_dist + 1
    log_sub("merging targets across barcodes (Hamming <= %d, ratio %g)"
            % (merge_dist, args.merge_ratio))
    con.execute("""
        CREATE OR REPLACE TABLE tstr AS
        SELECT s.target, s.seq, length(s.seq) AS len, t.reads
        FROM (SELECT DISTINCT md5_number(target_seq) AS target, target_seq AS seq FROM %s) s
        JOIN tseq t USING (target)
    """ % pairs_csv(src))
    # step2 : pigeonhole: merge_dist + 1 seeds, one row per (target, slot)
    con.execute("""
        CREATE OR REPLACE TABLE seeds AS
        SELECT target, seq, reads, len, slot,
               substr(seq, slot * (len // {k}) + 1, len // {k}) AS seed
        FROM tstr CROSS JOIN range({k}) AS r(slot)
        WHERE len // {k} >= 4
    """.format(k=n_seeds))
    # step3 : best table: best more abundant parent of every child (Hamming <= merge_dist, coverage ratio)
    con.execute("""
        CREATE OR REPLACE TABLE best AS
        SELECT child, parent FROM (
            SELECT c.target AS child, p.target AS parent,
                   hamming(c.seq, p.seq) AS dist, p.reads AS preads, p.seq AS pseq
            FROM seeds c JOIN seeds p
              ON c.len = p.len AND c.slot = p.slot AND c.seed = p.seed
            WHERE c.target <> p.target
              AND (p.reads > c.reads OR (p.reads = c.reads AND p.seq < c.seq))
              AND p.reads >= {ratio!r} * c.reads
        )
        WHERE dist <= {d}
        QUALIFY ROW_NUMBER() OVER (PARTITION BY child ORDER BY dist, preads DESC, pseq) = 1
    """.format(ratio=float(args.merge_ratio), d=merge_dist))
    con.execute("DROP TABLE seeds")
    # step4 : tmap: child -> best parent, parentless targets map to themselves
    con.execute("""
        CREATE OR REPLACE TABLE tmap AS
        SELECT t.target, COALESCE(b.parent, t.target) AS root
        FROM tseq t LEFT JOIN best b ON t.target = b.child
    """)
    # step5 : pointer jumping: follow chains to the root
    while True:
        con.execute("""
            CREATE OR REPLACE TABLE tmap_next AS
            SELECT a.target, b.root FROM tmap a JOIN tmap b ON a.root = b.target
        """)
        n_changed = con.execute("""
            SELECT COUNT(*) FROM tmap a JOIN tmap_next n USING (target)
            WHERE a.root <> n.root""").fetchone()[0]
        con.execute("DROP TABLE tmap")
        con.execute("ALTER TABLE tmap_next RENAME TO tmap")
        if not n_changed:
            break
    # step6 : report the number of merged variants
    n_merged = con.execute("SELECT COUNT(*) FROM tmap WHERE target <> root").fetchone()[0]
    log_sub("%d target variants merged into a more abundant target" % n_merged)


def assign_targets(con, args, barcode_meta, merged):
    """
    Choose ONE target per barcode (the one with the most reads) and give every barcode a status:
      low_coverage: top pair has < min_pair_reads reads
      ambiguous   : top target holds < min_purity of the barcode's reads
      resolved    : barcode had several targets, top one passed both filters
      unique      : barcode had a single target
    Creates tables target_ids, barcode_assign and final (status unique / resolved).
    Parameters:
            -- con: duckdb connection with tables pairs, tseq (and tmap if merged)
            -- args: argparse.Namespace: min_pair_reads, min_purity
            -- barcode_meta: str: <out>.meta/barcode_clusters.parquet
            -- merged: bool: True = --merge_dist > 0 (table tmap exists)
    Returns:
            -- None: nothing; creates the tables
    """
    # step1 : target ids 1, 2, ... (most reads first; ties: smaller md5 key)
    #   target_id breaks ties between equally covered targets of a barcode.  Such a barcode has a
    #   purity <= 0.5, so with --min_purity > 0.5 the tie never changes the result files.
    roots = ("SELECT m.root, SUM(t.reads) AS reads FROM tmap m JOIN tseq t USING (target) GROUP BY m.root"
             if merged else "SELECT target AS root, reads FROM tseq")
    con.execute("""
        CREATE OR REPLACE TABLE target_ids AS
        SELECT root AS tk, ROW_NUMBER() OVER (ORDER BY reads DESC, root) AS target_id FROM (%s)
    """ % roots)
    con.execute("DROP TABLE tseq")

    # step2 : reads per (barcode, target id) -- a view, read once by step3
    #   pair_reads = reads of the pair after target clustering, n_variants = number of different
    #   target sequences merged into the pair, first_/second_cov = the two largest variant coverages
    if merged:
        # S4 clusters merged into one root: sum the reads, top two of all variant coverages
        # (the 2nd largest is the 2nd largest first_cov or the largest second_cov)
        con.execute("""
            CREATE OR REPLACE VIEW bt AS
            SELECT p.bk, MAX(p.barcode_reads) AS barcode_reads, i.target_id, i.tk,
                   SUM(p.target_reads) AS pair_reads, SUM(p.target_unique) AS n_variants,
                   MAX(p.first_cov) AS first_cov,
                   greatest(COALESCE(max(p.first_cov, 2)[2], 0), MAX(p.second_cov)) AS second_cov
            FROM pairs p
            JOIN tmap m ON p.tk = m.target
            JOIN target_ids i ON m.root = i.tk
            GROUP BY p.bk, i.target_id, i.tk
        """)
    else:
        # one S4 row per (barcode cluster, target cluster): no aggregation needed
        con.execute("""
            CREATE OR REPLACE VIEW bt AS
            SELECT p.bk, p.barcode_reads, i.target_id, i.tk,
                   p.target_reads AS pair_reads, p.target_unique AS n_variants,
                   p.first_cov, p.second_cov
            FROM pairs p JOIN target_ids i ON p.tk = i.tk
        """)

    # step3 : ONE target per barcode and its status
    #   GROUP BY + arg_min (most reads, then smallest target_id) instead of a window over all pairs
    log_sub("choosing one target per barcode")
    con.execute("""
        CREATE OR REPLACE TABLE barcode_assign AS
        WITH r AS (
            SELECT bk,
                   arg_min({{'target_id': target_id, 'tk': tk, 'n_variants': n_variants,
                             'first_cov': first_cov, 'second_cov': second_cov}},
                           (-pair_reads, target_id))   AS top,
                   MAX(pair_reads)                      AS pair_reads,
                   MAX(barcode_reads)                   AS barcode_reads,
                   COUNT(*)                             AS n_targets
            FROM bt GROUP BY bk)
        SELECT bk, top.target_id AS target_id, top.tk AS tk,
               pair_reads, barcode_reads, top.n_variants AS n_variants,
               top.first_cov AS first_cov, top.second_cov AS second_cov,
               n_targets,
               CASE WHEN pair_reads < {mp}                            THEN 'low_coverage'
                    WHEN pair_reads * 1.0 / barcode_reads < {pur!r} THEN 'ambiguous'
                    WHEN n_targets > 1                                THEN 'resolved'
                    ELSE 'unique' END AS status
        FROM r
    """.format(mp=int(args.min_pair_reads), pur=float(args.min_purity)))
    con.execute("DROP VIEW bt")
    con.execute("DROP TABLE target_ids")
    #      |----   substep report the number of barcodes per status
    for status, n in con.execute(
            "SELECT status, COUNT(*) FROM barcode_assign GROUP BY status ORDER BY 2 DESC").fetchall():
        log_sub("barcodes %-13s %d" % (status, n))

    # step4 : final pairs = barcodes whose status is unique / resolved, with their S2 statistics
    #   barcode_clusters.parquet holds every barcode cluster: the inner join keeps all final rows
    con.execute("""
        CREATE OR REPLACE TABLE final AS
        SELECT a.bk, a.tk, a.target_id, a.pair_reads, a.barcode_reads, a.n_variants,
               a.first_cov, a.second_cov,
               m."count" AS bc_count, m.n_cluster AS bc_n_cluster,
               m.first_cov AS bc_first_cov, m.second_cov AS bc_second_cov
        FROM barcode_assign a
        JOIN (SELECT md5_number(barcode) AS bk, "count", n_cluster, first_cov, second_cov
              FROM read_parquet(%s)) m USING (bk)
        WHERE a.status IN ('unique', 'resolved')
    """ % sql_str(barcode_meta))


def final_rows_sql(src, merged):
    """
    SQL that adds the barcode and target sequences to the final pairs (one pass over
    target_pairs.tsv; the sequences are never stored in a table).
    Parameters:
            -- src: str: path of target_pairs.tsv
            -- merged: bool: True = --merge_dist > 0 (target sequence = root sequence from tstr)
    Returns:
            -- str: SQL query returning barcode, target and the columns of table final
    """
    if merged:
        # one row per barcode cluster: target_rank 1; the target is the root (maybe another barcode's)
        return ("SELECT c.barcode, s.seq AS target, f.* "
                "FROM (SELECT barcode, md5_number(barcode) AS bk FROM %s WHERE target_rank = 1) c "
                "JOIN final f USING (bk) JOIN tstr s ON s.target = f.tk" % pairs_csv(src))
    return ("SELECT c.barcode, c.target_seq AS target, f.* FROM %s c "
            "JOIN final f ON f.bk = md5_number(c.barcode) AND f.tk = md5_number(c.target_seq)"
            % pairs_csv(src))


def write_details_file(con, args, src, merged):
    """
    Write <out>.barcode_target.details.tsv: one row per final (barcode, target) pair.
    Every final barcode has exactly one target, so the barcode-cluster statistics belong to the
    pair as well.
    Parameters:
            -- con: duckdb connection with table final (and tstr if merged)
            -- args: argparse.Namespace: out
            -- src: str: path of target_pairs.tsv
            -- merged: bool: True = --merge_dist > 0
    Returns:
            -- None: nothing; writes the file
    """
    # step1 : columns of the pair
    #   barcode side (S2): the barcode cluster of the pair
    #     barcode_reads               = reads of all raw barcodes clustered into this barcode
    #     n_diff_barcodes             = raw barcodes clustered into this barcode
    #     first_/second_barcode_cov   = reads of the most / second most abundant raw barcode
    #     first_/second_barcode_ratio = cov / barcode_reads
    #   target side (S4, and S5 with --merge_dist): the target cluster under this barcode
    #     target_reads                = reads of all raw target sequences clustered into this target
    #                                   under this barcode
    #     n_diff_targets              = raw target sequences clustered into this target under this barcode
    #     first_/second_target_cov    = reads of the most / second most abundant target sequence
    #     first_/second_target_ratio  = cov / target_reads
    out = args.out + ".barcode_target.details.tsv"
    log_sub("writing %s" % out)
    con.execute("""
        COPY (SELECT barcode, target,
                     bc_count AS barcode_reads,
                     bc_n_cluster AS n_diff_barcodes,
                     bc_first_cov AS first_barcode_cov,
                     round(bc_first_cov * 1.0 / bc_count, 4) AS first_barcode_ratio,
                     bc_second_cov AS second_barcode_cov,
                     round(bc_second_cov * 1.0 / bc_count, 4) AS second_barcode_ratio,
                     pair_reads AS target_reads,
                     n_variants AS n_diff_targets,
                     first_cov AS first_target_cov,
                     round(first_cov * 1.0 / pair_reads, 4) AS first_target_ratio,
                     second_cov AS second_target_cov,
                     round(second_cov * 1.0 / pair_reads, 4) AS second_target_ratio
              FROM (%s)
              ORDER BY bc_count DESC, pair_reads DESC, barcode)
        TO %s (DELIMITER '\t', HEADER)""" % (final_rows_sql(src, merged), sql_str(out)))


# ---- stage driver ----------------------------------------------------------

def resolve_pairs(args, threads):
    """
    S5: decide ONE target per barcode and write the two result files.
    A barcode is final when its best target has >= --min_pair_reads reads and >= --min_purity of
    the barcode's reads (status 'unique' or 'resolved').
    Parameters:
            -- args: argparse.Namespace: out, db, mem_gb, merge_dist, merge_ratio, min_pair_reads,
               min_purity, clean_meta
            -- threads: int: DuckDB threads
    Returns:
            -- None: nothing; results are files
    """
    # step1 : check duckdb and the intermediate files, then open DuckDB
    try:
        import duckdb
    except ImportError:
        raise SystemExit("ERROR: S5 needs duckdb (pip install duckdb). "
                         "The intermediate files are written; rerun later with --resolve_only")
    #      |----   substep locate the intermediate files written by S1 / S2 / S4
    src = meta_path(args, "target_pairs.tsv")
    barcode_meta = meta_path(args, "barcode_clusters.parquet")
    raw_stats_file = meta_path(args, "raw_stats.json")
    for f in (src, barcode_meta, raw_stats_file):
        if not os.path.exists(f):
            raise SystemExit("ERROR: %s not found (run all steps first)" % f)
    with open(raw_stats_file) as fh:
        raw_stats = json.load(fh)
    #      |----   substep ingest_stats.json is missing for work directories of older versions
    ingest_stats_file = meta_path(args, "ingest_stats.json")
    ingest_stats = None
    if os.path.exists(ingest_stats_file):
        with open(ingest_stats_file) as fh:
            ingest_stats = json.load(fh)

    log("S5: resolving barcode -> target with DuckDB ...")
    #      |----   substep DuckDB with the --mem_gb / --threads limits and a spill directory
    con = duckdb.connect(args.db or ":memory:")
    spill = args.out + ".duckdb.tmp"
    con.execute("SET threads=%d" % threads)
    con.execute("SET memory_limit=%s" % sql_str("%dMB" % int(args.mem_gb * 1000)))
    con.execute("SET temp_directory=%s" % sql_str(spill))
    con.execute("SET preserve_insertion_order=false")      # lets DuckDB stream/spill freely

    # step2 : load the (barcode, target cluster) pairs, sequences replaced by their md5 keys
    log_sub("loading (barcode, target cluster) pairs")
    con.execute("""
        CREATE OR REPLACE TABLE pairs AS
        SELECT md5_number(barcode) AS bk, md5_number(target_seq) AS tk,
               barcode_reads, target_reads, target_unique, first_cov, second_cov
        FROM %s
    """ % pairs_csv(src))

    # step3 : distinct targets (tmap: target -> root, only with --merge_dist)
    con.execute("""
        CREATE OR REPLACE TABLE tseq AS
        SELECT tk AS target, SUM(target_reads) AS reads, COUNT(*) AS n_barcodes
        FROM pairs GROUP BY tk
    """)

    # step4 : optional: merge near-identical targets ACROSS barcodes (--merge_dist)
    merged = args.merge_dist > 0
    if merged:
        merge_targets(con, args, src)

    # step5 : one target per barcode (tables barcode_assign, final) and the statistics
    assign_targets(con, args, barcode_meta, merged)
    stats = collect_stats(con, raw_stats, ingest_stats, args)
    con.execute("DROP TABLE pairs")
    for table in ("tmap", "best"):
        con.execute("DROP TABLE IF EXISTS %s" % table)

    # step6 : result files
    write_details_file(con, args, src, merged)
    log_sub("writing %s" % (args.out + ".stats.tsv"))
    write_stats(args.out + ".stats.tsv", stats)

    # step7 : release DuckDB, delete the spill directory and (--clean_meta) the meta directory
    con.close()
    if not args.db:
        shutil.rmtree(spill, ignore_errors=True)
    if args.clean_meta:
        shutil.rmtree(meta_path(args), ignore_errors=True)
