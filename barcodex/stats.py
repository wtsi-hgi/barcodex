"""
S5 statistics: <out>.stats.tsv, ordered by step.

Columns
  step         the step the number belongs to (1_input ... 6_final_map)
  stat         name of the number
  value        the number
  percent      for read counts: value / total_read_pairs x 100;  NA for other numbers
  description  what the number means

The read counts add up:
  total_read_pairs     = all rejected rows of 2_extraction + extracted_read_pairs
  extracted_read_pairs = reads_dropped_barcode_clustering (3)
                         + reads_low_coverage_barcodes + reads_ambiguous_barcodes
                         + reads_minor_targets (5) + effective_reads (6)
"""

from .common import log_sub, sql_str

STEPS = ("1_input", "2_extraction", "3_barcode_clustering", "4_target_clustering",
         "5_target_assignment", "6_final_map")


class Stat:
    """
    One row of the stats file.
    Parameters:
            -- step: str: one of STEPS
            -- stat: str: name
            -- value: int, float or None: the number
            -- description: str: what the number means
            -- is_reads: bool: True = a read count (gets a percent of total_read_pairs)
    Returns:
            -- Stat: row object
    """

    def __init__(self, step, stat, value, description, is_reads=False):
        """
        Store the fields.
        Parameters:
                -- see the class docstring
        Returns:
                -- None: constructor
        """
        self.step, self.stat, self.value = step, stat, value
        self.description, self.is_reads = description, is_reads


def number(x):
    """
    Format a statistic: integers without decimals, other numbers with 2 decimals, None as 0.
    Parameters:
            -- x: int, float or None: value to format
    Returns:
            -- str: formatted number
    """
    return "0" if x is None else ("%d" % x if float(x).is_integer() else "%.2f" % x)


# ---- one function per step -------------------------------------------------

def extraction_stats(ingest_stats):
    """
    Steps 1 and 2: read pairs in the input, rejected by the extraction and kept.
    Parameters:
            -- ingest_stats: dict or None: content of ingest_stats.json (S1); None = older run
    Returns:
            -- list[Stat]: rows (empty for older runs)
    """
    if ingest_stats is None:
        return []
    rows = [Stat("1_input", "total_read_pairs", ingest_stats["total_read_pairs"],
                 "read pairs in the FASTQ files", True)]
    rows += [Stat("2_extraction", key, n, description, True)
             for key, n, description in ingest_stats["rejects"]]
    rows.append(Stat("2_extraction", "extracted_read_pairs", ingest_stats["extracted_read_pairs"],
                     "read pairs with a valid barcode and target", True))
    return rows


def barcode_clustering_stats(con, raw_stats):
    """
    Step 3: barcodes before and after the barcode clustering, reads lost in it.
    Parameters:
            -- con: duckdb connection with table barcode_assign
            -- raw_stats: dict: content of raw_stats.json (S2)
    Returns:
            -- list[Stat]: rows
    """
    n_clusters = con.execute("SELECT COUNT(*) FROM barcode_assign").fetchone()[0]
    rows = [Stat("3_barcode_clustering", "raw_barcodes", raw_stats["raw_barcodes"],
                 "distinct barcode sequences before clustering"),
            Stat("3_barcode_clustering", "barcode_clusters", n_clusters,
                 "barcode clusters with reads after clustering")]
    if "reads_clustered" in raw_stats:
        rows += [Stat("3_barcode_clustering", "reads_dropped_barcode_clustering",
                      raw_stats["reads_in"] - raw_stats["reads_clustered"],
                      "reads of barcodes without a parent and < --min_barcode_reads reads "
                      "(or < --min_child_reads)", True),
                 Stat("3_barcode_clustering", "reads_in_barcode_clusters",
                      raw_stats["reads_clustered"], "reads in the barcode clusters", True)]
    return rows


def target_clustering_stats(con, raw_stats):
    """
    Step 4: targets before and after the target clustering within each barcode.
    Parameters:
            -- con: duckdb connection with table pairs (S4 output)
            -- raw_stats: dict: content of raw_stats.json (S2)
    Returns:
            -- list[Stat]: rows
    """
    n_before, n_after = con.execute(
        "SELECT SUM(target_unique), COUNT(*) FROM pairs").fetchone()
    return [Stat("4_target_clustering", "raw_targets", raw_stats["raw_targets"],
                 "distinct target sequences before clustering"),
            Stat("4_target_clustering", "barcode_target_variants", n_before,
                 "distinct (barcode cluster, target sequence) pairs before target clustering"),
            Stat("4_target_clustering", "barcode_target_clusters", n_after,
                 "(barcode cluster, target cluster) pairs after target clustering")]


def assignment_stats(con, args):
    """
    Step 5: one target per barcode: barcodes per status and the reads that are not used.
    Parameters:
            -- con: duckdb connection with table barcode_assign
            -- args: argparse.Namespace: min_pair_reads, min_purity
    Returns:
            -- list[Stat]: rows
    """
    (n_unique, n_resolved, n_ambiguous, n_low, n_multi,
     reads_low, reads_ambiguous, reads_minor) = con.execute("""
        SELECT COUNT(*) FILTER (WHERE status = 'unique'),
               COUNT(*) FILTER (WHERE status = 'resolved'),
               COUNT(*) FILTER (WHERE status = 'ambiguous'),
               COUNT(*) FILTER (WHERE status = 'low_coverage'),
               COUNT(*) FILTER (WHERE n_targets > 1),
               COALESCE(SUM(barcode_reads) FILTER (WHERE status = 'low_coverage'), 0),
               COALESCE(SUM(barcode_reads) FILTER (WHERE status = 'ambiguous'), 0),
               COALESCE(SUM(barcode_reads - pair_reads)
                        FILTER (WHERE status IN ('unique', 'resolved')), 0)
        FROM barcode_assign""").fetchone()
    step = "5_target_assignment"
    return [
        Stat(step, "barcodes_unique", n_unique, "barcodes with a single target: kept"),
        Stat(step, "barcodes_multi_target", n_multi, "barcodes with more than one target cluster"),
        Stat(step, "barcodes_resolved", n_resolved,
             "multi-target barcodes whose top target passed the filters: kept"),
        Stat(step, "barcodes_ambiguous", n_ambiguous,
             "top target has < %g of the barcode's reads (--min_purity): dropped"
             % args.min_purity),
        Stat(step, "barcodes_low_coverage", n_low,
             "top target has < %d reads (--min_pair_reads): dropped" % args.min_pair_reads),
        Stat(step, "reads_low_coverage_barcodes", reads_low, "reads of low-coverage barcodes",
             True),
        Stat(step, "reads_ambiguous_barcodes", reads_ambiguous, "reads of ambiguous barcodes",
             True),
        Stat(step, "reads_minor_targets", reads_minor,
             "reads of kept barcodes that support another (minor) target", True),
    ]


def final_map_stats(con, barcode_meta):
    """
    Step 6: the final barcode -> target map.
    Parameters:
            -- con: duckdb connection with table final
            -- barcode_meta: str: <out>.meta/barcode_clusters.parquet
    Returns:
            -- list[Stat]: rows
    """
    # step1 : barcodes, targets, pairs and reads
    n_barcodes, bc_mean, bc_median = con.execute("""
        SELECT COUNT(*), AVG(m."count"), median(m."count")
        FROM final f JOIN read_parquet(%s) m ON m.barcode = f.barcode""" % sql_str(barcode_meta)).fetchone()
    n_targets, tg_mean, tg_median = con.execute("""
        SELECT COUNT(*), AVG(reads), median(reads)
        FROM (SELECT target_id, SUM(pair_reads) AS reads FROM final GROUP BY target_id)""").fetchone()
    n_pairs, effective = con.execute(
        "SELECT COUNT(*), COALESCE(SUM(pair_reads), 0) FROM final").fetchone()
    # step2 : rows
    step = "6_final_map"
    return [
        Stat(step, "final_barcodes", n_barcodes, "barcodes in the final barcode -> target map"),
        Stat(step, "final_targets", n_targets, "distinct targets in the final map"),
        Stat(step, "final_barcode_target_pairs", n_pairs, "final (barcode, target) pairs"),
        Stat(step, "effective_reads", effective,
             "reads supporting the final (barcode, target) pairs", True),
        Stat(step, "barcode_coverage_mean", number(bc_mean),
             "mean reads per final barcode (after clustering)"),
        Stat(step, "barcode_coverage_median", number(bc_median),
             "median reads per final barcode (after clustering)"),
        Stat(step, "target_coverage_mean", number(tg_mean),
             "mean reads per final target (all its barcodes)"),
        Stat(step, "target_coverage_median", number(tg_median),
             "median reads per final target (all its barcodes)"),
        Stat(step, "mean_barcodes_per_target", number(n_pairs / n_targets if n_targets else 0.0),
             "on average a target has this many barcodes"),
    ]


# ---- all rows and the file -------------------------------------------------

def collect_stats(con, barcode_meta, raw_stats, ingest_stats, args):
    """
    All rows of the stats file, ordered by step.
    Parameters:
            -- con: duckdb connection after assign_targets (tables pairs, final, barcode_assign)
            -- barcode_meta: str: <out>.meta/barcode_clusters.parquet
            -- raw_stats: dict: content of raw_stats.json (S2)
            -- ingest_stats: dict or None: content of ingest_stats.json (S1)
            -- args: argparse.Namespace: min_pair_reads, min_purity
    Returns:
            -- list[Stat]: rows
    """
    return (extraction_stats(ingest_stats) + barcode_clustering_stats(con, raw_stats)
            + target_clustering_stats(con, raw_stats) + assignment_stats(con, args)
            + final_map_stats(con, barcode_meta))


def write_stats(path, stats):
    """
    Write the stats file (step / stat / value / percent / description) and echo it to the log,
    one block per step.
    Parameters:
            -- path: str: output file
            -- stats: list[Stat]: rows
    Returns:
            -- None: nothing
    """
    # step1 : percent of total_read_pairs for read counts
    total = next((s.value for s in stats if s.stat == "total_read_pairs"), 0)

    def percent(row):
        """
        Percent column of one row.
        Parameters:
                -- row: Stat: row
        Returns:
                -- str: e.g. '74.69', or 'NA'
        """
        return "%.2f" % (100.0 * row.value / total) if row.is_reads and total else "NA"

    # step2 : the file
    with open(path, "w") as fo:
        fo.write("step\tstat\tvalue\tpercent\tdescription\n")
        for row in stats:
            fo.write("%s\t%s\t%s\t%s\t%s\n" % (row.step, row.stat, row.value, percent(row),
                                               row.description))
    # step3 : the log, one block per step
    step = None
    for row in stats:
        if row.step != step:
            step = row.step
            log_sub(step)
        log_sub("    %-34s %14s %8s" % (row.stat, row.value, percent(row)))
