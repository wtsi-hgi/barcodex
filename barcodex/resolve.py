"""
S5 resolve: one target per barcode, then the result files (DuckDB).

Inputs  (intermediate files in <out>.meta/, written by S1 / S2 / S4):
  target_pairs.tsv           (barcode cluster, target cluster) rows
  barcode_clusters.parquet   per barcode-cluster statistics
  raw_stats.json             raw barcode / target counts, reads kept by the barcode clustering
  ingest_stats.json          read pairs parsed / rejected by the extraction
Outputs:
  <out>.barcode_clusters.tsv   one row per final barcode
  <out>.target_clusters.tsv    one row per final (barcode, target) pair
  <out>.stats.tsv              summary numbers (see stats.py)
--resolve_only re-runs this stage from <out>.meta/.
"""

import json
import os
import shutil

from .common import log, log_sub, meta_path, sql_str
from .stats import collect_stats, write_stats


# ---- steps -----------------------------------------------------------------

def merge_targets(con, args):
    """
    Optional (--merge_dist > 0): merge near-identical targets ACROSS barcodes.  Rewrites table
    tmap (target -> root target) so that every child points to its most abundant root.
    Parameters:
            -- con: duckdb connection with table tseq (target, reads, n_barcodes, len)
            -- args: argparse.Namespace: merge_dist, merge_ratio
    Returns:
            -- None: nothing; replaces table tmap
    """
    # step1 : pigeonhole: merge_dist + 1 seeds
    merge_dist = args.merge_dist
    n_seeds = merge_dist + 1
    log_sub("merging targets across barcodes (Hamming <= %d, ratio %g)"
            % (merge_dist, args.merge_ratio))
    #      |----   substep seeds table: one row per (target, slot)
    con.execute("""
        CREATE OR REPLACE TABLE seeds AS
        SELECT target, reads, len, slot,
               substr(target, slot * (len // {k}) + 1, len // {k}) AS seed
        FROM tseq CROSS JOIN range({k}) AS r(slot)
        WHERE len // {k} >= 4
    """.format(k=n_seeds))
    # step2 : best table: best more abundant parent of every child (Hamming <= merge_dist, coverage ratio)
    con.execute("""
        CREATE OR REPLACE TABLE best AS
        SELECT child, parent FROM (
            SELECT c.target AS child, p.target AS parent,
                   hamming(c.target, p.target) AS dist, p.reads AS preads
            FROM seeds c JOIN seeds p
              ON c.len = p.len AND c.slot = p.slot AND c.seed = p.seed
            WHERE c.target <> p.target
              AND (p.reads > c.reads OR (p.reads = c.reads AND p.target < c.target))
              AND p.reads >= {ratio!r} * c.reads
        )
        WHERE dist <= {d}
        QUALIFY ROW_NUMBER() OVER (PARTITION BY child ORDER BY dist, preads DESC, parent) = 1
    """.format(ratio=float(args.merge_ratio), d=merge_dist))
    # step3 : tmap: child -> best parent, parentless targets map to themselves
    con.execute("""
        CREATE OR REPLACE TABLE tmap AS
        SELECT t.target, COALESCE(b.parent, t.target) AS root
        FROM tseq t LEFT JOIN best b ON t.target = b.child
    """)
    # step4 : pointer jumping: follow chains to the root
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
    # step5 : report the number of merged variants
    n_merged = con.execute("SELECT COUNT(*) FROM tmap WHERE target <> root").fetchone()[0]
    log_sub("%d target variants merged into a more abundant target" % n_merged)


def assign_targets(con, args):
    """
    Choose ONE target per barcode (the one with the most reads) and give every barcode a status:
      low_coverage: top pair has < min_pair_reads reads
      ambiguous   : top target holds < min_purity of the barcode's reads
      resolved    : barcode had several targets, top one passed both filters
      unique      : barcode had a single target
    Creates tables target_ids, bt, barcode_assign and final (status unique / resolved).
    Parameters:
            -- con: duckdb connection with tables pairs, tseq, tmap
            -- args: argparse.Namespace: min_pair_reads, min_purity
    Returns:
            -- None: nothing; creates the tables
    """
    # step1 : target ids T00000001, ... (most reads first)
    con.execute("""
        CREATE OR REPLACE TABLE target_ids AS
        SELECT root AS target_seq,
               'T' || lpad(CAST(ROW_NUMBER() OVER (ORDER BY reads DESC, root) AS VARCHAR), 8, '0') AS target_id
        FROM (SELECT m.root, SUM(t.reads) AS reads
              FROM tmap m JOIN tseq t USING (target) GROUP BY m.root)
    """)

    # step2 : reads per (barcode, target id)
    #   pair_reads = reads of the pair after target clustering, n_variants = number of
    #   different target sequences merged into the pair
    con.execute("""
        CREATE OR REPLACE TABLE bt AS
        SELECT p.barcode, MAX(p.barcode_reads) AS barcode_reads, i.target_id,
               SUM(p.target_reads) AS pair_reads, SUM(p.target_unique) AS n_variants
        FROM pairs p
        JOIN tmap m ON p.target_seq = m.target
        JOIN target_ids i ON m.root = i.target_seq
        GROUP BY p.barcode, i.target_id
    """)

    # step3 : ONE target per barcode and its status
    log_sub("choosing one target per barcode")
    con.execute("""
        CREATE OR REPLACE TABLE barcode_assign AS
        WITH r AS (
            SELECT *,
                   ROW_NUMBER() OVER (PARTITION BY barcode ORDER BY pair_reads DESC, target_id) AS rn,
                   COUNT(*)     OVER (PARTITION BY barcode) AS n_targets
            FROM bt)
        SELECT r.barcode, r.target_id, i.target_seq,
               r.pair_reads, r.barcode_reads, r.n_variants,
               r.pair_reads * 1.0 / r.barcode_reads AS purity,
               r.n_targets,
               CASE WHEN r.pair_reads < {mp}                            THEN 'low_coverage'
                    WHEN r.pair_reads * 1.0 / r.barcode_reads < {pur!r} THEN 'ambiguous'
                    WHEN r.n_targets > 1                                THEN 'resolved'
                    ELSE 'unique' END AS status
        FROM r JOIN target_ids i USING (target_id)
        WHERE r.rn = 1
    """.format(mp=int(args.min_pair_reads), pur=float(args.min_purity)))
    #      |----   substep report the number of barcodes per status
    for status, n in con.execute(
            "SELECT status, COUNT(*) FROM barcode_assign GROUP BY status ORDER BY 2 DESC").fetchall():
        log_sub("barcodes %-13s %d" % (status, n))

    # step4 : final pairs = barcodes whose status is unique / resolved
    con.execute("""
        CREATE OR REPLACE TABLE final AS
        SELECT barcode, target_id, target_seq, pair_reads, barcode_reads, n_variants
        FROM barcode_assign WHERE status IN ('unique','resolved')
    """)


def write_cluster_files(con, args, barcode_meta):
    """
    Write <out>.barcode_clusters.tsv and <out>.target_clusters.tsv for the final pairs.
    Parameters:
            -- con: duckdb connection with tables pairs, tmap, target_ids, final
            -- args: argparse.Namespace: out
            -- barcode_meta: str: <out>.meta/barcode_clusters.parquet
    Returns:
            -- None: nothing; writes the files
    """
    # step1 : top-two target-variant coverages of every final pair
    # Every S4 cluster row carries the reads of its top and second variant; if several S4
    # clusters were merged into one target (--merge_dist), the top two of all of them are taken.
    #      |----   substep helper: coverage rows of the S4 clusters that belong to the final pairs
    def variant_rows(column, extra_filter=""):
        """
        SQL that lists, for every final (barcode, target id), the coverage column of the S4 cluster
        rows that were merged into it.
        Parameters:
                -- column: str: 'first_cov' or 'second_cov'
                -- extra_filter: str: optional WHERE clause appended to the query
        Returns:
                -- str: SQL query returning (barcode, target_id, cov)
        """
        return ("SELECT p.barcode, f.target_id, p.%s AS cov FROM pairs p "
                "JOIN tmap m ON p.target_seq = m.target "
                "JOIN target_ids i ON m.root = i.target_seq "
                "JOIN final f ON f.barcode = p.barcode AND f.target_id = i.target_id %s"
                % (column, extra_filter))
    con.execute("""
        CREATE OR REPLACE TABLE final_cov AS
        SELECT barcode, target_id,
               MAX(cov) FILTER (WHERE rn = 1) AS first_cov,
               COALESCE(MAX(cov) FILTER (WHERE rn = 2), 0) AS second_cov
        FROM (SELECT barcode, target_id, cov,
                     ROW_NUMBER() OVER (PARTITION BY barcode, target_id ORDER BY cov DESC) AS rn
              FROM (%s UNION ALL %s))
        GROUP BY barcode, target_id
    """ % (variant_rows("first_cov"), variant_rows("second_cov", "WHERE p.second_cov > 0")))

    # step2 : barcode clustering file: one row per final barcode
    #   barcode_reads    = reads of the barcode cluster (coverage after clustering)
    #   n_diff_barcodes  = number of different raw barcodes merged into it (before clustering)
    #   first_/second_barcode_cov = reads of the most / second most abundant raw barcode,
    #   first_/second_barcode_ratio = cov / barcode_reads
    tsv = "(DELIMITER '\t', HEADER)"
    log_sub("writing %s" % (args.out + ".barcode_clusters.tsv"))
    con.execute("""
        COPY (SELECT f.barcode AS barcode, f.target_seq AS target,
                     m."count" AS barcode_reads,
                     m.n_cluster AS n_diff_barcodes,
                     m.first_cov AS first_barcode_cov,
                     round(m.first_cov * 1.0 / m."count", 4) AS first_barcode_ratio,
                     m.second_cov AS second_barcode_cov,
                     round(m.second_cov * 1.0 / m."count", 4) AS second_barcode_ratio
              FROM final f JOIN read_parquet(%s) m ON m.barcode = f.barcode
              ORDER BY m."count" DESC, f.barcode)
        TO %s %s""" % (sql_str(barcode_meta), sql_str(args.out + ".barcode_clusters.tsv"), tsv))

    # step3 : target clustering file: one row per final (barcode, target) pair
    #   target_reads    = reads of the pair after target clustering
    #   n_diff_targets  = number of different target sequences merged into it (before clustering)
    #   first_/second_target_cov = reads of the most / second most abundant target sequence,
    #   first_/second_target_ratio = cov / target_reads
    log_sub("writing %s" % (args.out + ".target_clusters.tsv"))
    con.execute("""
        COPY (SELECT f.barcode AS barcode, f.target_seq AS target,
                     f.pair_reads AS target_reads,
                     f.n_variants AS n_diff_targets,
                     c.first_cov AS first_target_cov,
                     round(c.first_cov * 1.0 / f.pair_reads, 4) AS first_target_ratio,
                     c.second_cov AS second_target_cov,
                     round(c.second_cov * 1.0 / f.pair_reads, 4) AS second_target_ratio
              FROM final f JOIN final_cov c ON c.barcode = f.barcode AND c.target_id = f.target_id
              ORDER BY f.pair_reads DESC, f.barcode)
        TO %s %s""" % (sql_str(args.out + ".target_clusters.tsv"), tsv))


# ---- stage driver ----------------------------------------------------------

def resolve_pairs(args, threads):
    """
    S5: decide ONE target per barcode and write the three result files.
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

    # step2 : load the (barcode, target cluster) pairs
    log_sub("loading (barcode, target cluster) pairs")
    con.execute("""
        CREATE OR REPLACE TABLE pairs AS
        SELECT * FROM read_csv(%s, delim='\t', header=true, columns={
            'barcode':'VARCHAR', 'barcode_reads':'BIGINT', 'target_rank':'INTEGER',
            'target_reads':'BIGINT', 'target_unique':'INTEGER', 'target_seq':'VARCHAR',
            'first_cov':'BIGINT', 'second_cov':'BIGINT'})
    """ % sql_str(src))

    # step3 : distinct target sequences (tmap: target -> root, identity for now)
    con.execute("""
        CREATE OR REPLACE TABLE tseq AS
        SELECT target_seq AS target, SUM(target_reads) AS reads,
               COUNT(*) AS n_barcodes, length(target_seq) AS len
        FROM pairs GROUP BY target_seq
    """)
    con.execute("CREATE OR REPLACE TABLE tmap AS SELECT target, target AS root FROM tseq")

    # step4 : optional: merge near-identical targets ACROSS barcodes (--merge_dist)
    if args.merge_dist > 0:
        merge_targets(con, args)

    # step5 : one target per barcode (tables target_ids, bt, barcode_assign, final)
    assign_targets(con, args)

    # step6 : result files
    write_cluster_files(con, args, barcode_meta)
    log_sub("writing %s" % (args.out + ".stats.tsv"))
    stats = collect_stats(con, barcode_meta, raw_stats, ingest_stats, args)
    write_stats(args.out + ".stats.tsv", stats)

    # step7 : release DuckDB, delete the spill directory and (--clean_meta) the meta directory
    con.close()
    if not args.db:
        shutil.rmtree(spill, ignore_errors=True)
    if args.clean_meta:
        shutil.rmtree(meta_path(args), ignore_errors=True)
