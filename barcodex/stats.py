"""
S5 statistics: <out>.stats.tsv, ordered by step.

Columns
  step         the step the number belongs to (1_input ... 6_final_map)
  stat         name of the number: reads_* = read (pair) counts, num_* = numbers of sequences,
               mean_* / median_* = averages
  count        the number
  percent      for reads_* rows: count / reads_raw_input x 100;  NA for other numbers
  description  what the number means

The read counts add up:
  reads_raw_input             = all reads_* reject rows of 2_extraction
                                + reads_dropped_barcode_clustering + reads_in_barcode_clustering
  reads_in_barcode_clustering = reads_dropped_target_clustering + reads_in_target_clustering
  reads_in_target_clustering  = reads_barcodes_low_coverage + reads_barcodes_ambiguous
                                + reads_final_effective
                                + reads of kept barcodes that support another (minor) target
                                  (only resolved barcodes have those, so they are
                                   <= reads_barcodes_resolved)
"""

from .common import log_sub

STEPS = ("1_input", "2_extraction", "3_barcode_clustering", "4_target_clustering",
         "5_target_assignment", "6_final_map")

# reject key suffix of ingest_stats.json (extract.reject_reasons) -> stat name suffix
REJECT_NAMES = {
    "upstream_not_found": "upstream_not_found",
    "downstream_not_found": "downstream_not_found",
    "empty": "missing",
    "too_short": "read_too_short",
    "too_long": "too_long",
    "below_min_len": "incorrect_length",
    "length_filtered": "length_not_allowed",
    "template_mismatch": "incorrect_template",
}


class Stat:
    """
    One row of the stats file.
    Parameters:
            -- step: str: one of STEPS
            -- stat: str: name
            -- count: int, float or None: the number
            -- description: str: what the number means
            -- is_reads: bool: True = a read count (gets a percent of reads_raw_input)
    Returns:
            -- Stat: row object
    """

    def __init__(self, step, stat, count, description, is_reads=False):
        """
        Store the fields.
        Parameters:
                -- see the class docstring
        Returns:
                -- None: constructor
        """
        self.step, self.stat, self.count = step, stat, count
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


def reject_stat_name(key):
    """
    Stat name of one extraction reject key, e.g. 'target_empty' -> 'reads_target_missing'.
    Parameters:
            -- key: str: '<role>_<reason>' from ingest_stats.json
    Returns:
            -- str: 'reads_<role>_<name>' (unknown reasons keep their own name)
    """
    role, _, reason = key.partition("_")
    return "reads_%s_%s" % (role, REJECT_NAMES.get(reason, reason))


# ---- one function per step -------------------------------------------------

def extraction_stats(ingest_stats, raw_stats):
    """
    Steps 1 and 2: read pairs in the input, rejected by the extraction (target reasons first,
    then barcode reasons), and the distinct sequences the extraction produced.
    Parameters:
            -- ingest_stats: dict or None: content of ingest_stats.json (S1); None = older run
            -- raw_stats: dict: content of raw_stats.json (S2)
    Returns:
            -- list[Stat]: rows (reads_* rows only when ingest_stats is known)
    """
    rows = []
    # step1 : input and reject rows; target rows first, each role in check order
    if ingest_stats is not None:
        rows.append(Stat("1_input", "reads_raw_input", ingest_stats["total_read_pairs"],
                         "read pairs in the FASTQ files", True))
        rejects = sorted(ingest_stats["rejects"], key=lambda r: not r[0].startswith("target_"))
        rows += [Stat("2_extraction", reject_stat_name(key), n, description, True)
                 for key, n, description in rejects]
    # step2 : distinct raw sequences and pairs (raw_barcode_target_pairs: runs from v>= this one)
    step = "2_extraction"
    rows += [Stat(step, "num_raw_barcodes", raw_stats["raw_barcodes"],
                  "distinct barcode sequences extracted (before barcode clustering)"),
             Stat(step, "num_raw_targets", raw_stats["raw_targets"],
                  "distinct target sequences extracted (before target clustering)")]
    if "raw_barcode_target_pairs" in raw_stats:
        rows.append(Stat(step, "num_raw_barcode_target_pairs",
                         raw_stats["raw_barcode_target_pairs"],
                         "distinct (barcode sequence, target sequence) pairs extracted "
                         "(before any clustering)"))
    return rows


def barcode_clustering_stats(con, raw_stats):
    """
    Step 3: reads lost / kept by the barcode clustering and the barcodes after it.
    Parameters:
            -- con: duckdb connection with table barcode_assign
            -- raw_stats: dict: content of raw_stats.json (S2)
    Returns:
            -- list[Stat]: rows
    """
    step = "3_barcode_clustering"
    n_clusters = con.execute("SELECT COUNT(*) FROM barcode_assign").fetchone()[0]
    rows = []
    if "reads_clustered" in raw_stats:
        rows += [Stat(step, "reads_dropped_barcode_clustering",
                      raw_stats["reads_in"] - raw_stats["reads_clustered"],
                      "reads of barcodes without a parent and < --min_barcode_reads reads "
                      "(or < --min_child_reads): dropped", True),
                 Stat(step, "reads_in_barcode_clustering", raw_stats["reads_clustered"],
                      "reads in the barcode clusters (kept by the barcode clustering)", True)]
    rows.append(Stat(step, "num_barcodes_after_clustering", n_clusters,
                     "barcodes (barcode clusters) after clustering"))
    return rows


def target_clustering_stats(con, raw_stats, merged):
    """
    Step 4: reads lost / kept by the target clustering within each barcode and the targets
    after it.
    Parameters:
            -- con: duckdb connection with table pairs (S4 output) (and tmap if merged)
            -- raw_stats: dict: content of raw_stats.json (S2)
            -- merged: bool: True = --merge_dist > 0 (table tmap exists)
    Returns:
            -- list[Stat]: rows
    """
    step = "4_target_clustering"
    reads_in, n_targets = con.execute(
        "SELECT COALESCE(SUM(target_reads), 0), COUNT(DISTINCT tk) FROM pairs").fetchone()
    text = "distinct targets (target cluster sequences) after clustering within each barcode"
    if merged:
        n_targets = con.execute("SELECT COUNT(DISTINCT root) FROM tmap").fetchone()[0]
        text = ("distinct targets after clustering within each barcode and merging across "
                "barcodes (--merge_dist)")
    rows = []
    if "reads_clustered" in raw_stats:
        rows.append(Stat(step, "reads_dropped_target_clustering",
                         raw_stats["reads_clustered"] - reads_in,
                         "reads whose target has non-ACGT bases: dropped", True))
    rows += [Stat(step, "reads_in_target_clustering", reads_in,
                  "reads in the target clusters (kept by the target clustering)", True),
             Stat(step, "num_target_after_clustering", n_targets, text)]
    return rows


def assignment_stats(con, args):
    """
    Step 5: one target per barcode: multi-target barcodes and the barcodes / reads per status.
    Parameters:
            -- con: duckdb connection with table barcode_assign
            -- args: argparse.Namespace: min_pair_reads, min_purity
    Returns:
            -- list[Stat]: rows
    """
    (n_multi, reads_multi, n_resolved, reads_resolved, n_ambiguous, reads_ambiguous,
     n_low, reads_low) = con.execute("""
        SELECT COUNT(*) FILTER (WHERE n_targets > 1),
               COALESCE(SUM(barcode_reads) FILTER (WHERE n_targets > 1), 0),
               COUNT(*) FILTER (WHERE status = 'resolved'),
               COALESCE(SUM(barcode_reads) FILTER (WHERE status = 'resolved'), 0),
               COUNT(*) FILTER (WHERE status = 'ambiguous'),
               COALESCE(SUM(barcode_reads) FILTER (WHERE status = 'ambiguous'), 0),
               COUNT(*) FILTER (WHERE status = 'low_coverage'),
               COALESCE(SUM(barcode_reads) FILTER (WHERE status = 'low_coverage'), 0)
        FROM barcode_assign""").fetchone()
    step = "5_target_assignment"
    mp, pur = args.min_pair_reads, args.min_purity
    resolved = ("multi-target barcodes resolved to ONE target: their top target has >= %d reads "
                "(--min_pair_reads) and >= %g of the barcode's reads (--min_purity); kept with "
                "that target" % (mp, pur))
    ambiguous = ("multi-target barcodes NOT resolved: their top target has >= %d reads "
                 "(--min_pair_reads) but < %g of the barcode's reads (--min_purity), no "
                 "dominant target; dropped" % (mp, pur))
    low = ("barcodes (single- or multi-target) whose top target has < %d reads "
           "(--min_pair_reads); dropped" % mp)
    return [
        Stat(step, "reads_barcodes_multi_targets", reads_multi,
             "reads of barcodes with more than one target after target clustering", True),
        Stat(step, "num_barcodes_multi_targets", n_multi,
             "barcodes with more than one target after target clustering"),
        Stat(step, "reads_barcodes_resolved", reads_resolved,
             "all reads of the resolved barcodes (only the top target's reads are kept)", True),
        Stat(step, "num_barcodes_resolved", n_resolved, resolved),
        Stat(step, "reads_barcodes_ambiguous", reads_ambiguous,
             "reads of the ambiguous barcodes", True),
        Stat(step, "num_barcodes_ambiguous", n_ambiguous, ambiguous),
        Stat(step, "reads_barcodes_low_coverage", reads_low,
             "reads of the low-coverage barcodes", True),
        Stat(step, "num_barcodes_low_coverage", n_low, low),
    ]


def final_map_stats(con):
    """
    Step 6: the final barcode -> target map.
    Parameters:
            -- con: duckdb connection with table final (bc_count = reads of the barcode cluster, S2)
    Returns:
            -- list[Stat]: rows
    """
    # step1 : barcodes, targets, pairs and reads
    n_barcodes, bc_mean, bc_median = con.execute(
        "SELECT COUNT(*), AVG(bc_count), median(bc_count) FROM final").fetchone()
    n_targets, tg_mean, tg_median = con.execute("""
        SELECT COUNT(*), AVG(reads), median(reads)
        FROM (SELECT target_id, SUM(pair_reads) AS reads FROM final GROUP BY target_id)""").fetchone()
    n_pairs, effective = con.execute(
        "SELECT COUNT(*), COALESCE(SUM(pair_reads), 0) FROM final").fetchone()
    # step2 : rows
    step = "6_final_map"
    return [
        Stat(step, "num_final_barcodes", n_barcodes,
             "barcodes in the final barcode -> target map (single-target + resolved)"),
        Stat(step, "num_final_targets", n_targets, "distinct targets in the final map"),
        Stat(step, "num_final_barcode_target_pairs", n_pairs, "final (barcode, target) pairs"),
        Stat(step, "reads_final_effective", effective,
             "reads supporting the final (barcode, target) pairs", True),
        Stat(step, "mean_barcode_coverage", number(bc_mean),
             "mean reads per final barcode (after clustering)"),
        Stat(step, "median_barcode_coverage", number(bc_median),
             "median reads per final barcode (after clustering)"),
        Stat(step, "mean_target_coverage", number(tg_mean),
             "mean reads per final target (all its barcodes)"),
        Stat(step, "median_target_coverage", number(tg_median),
             "median reads per final target (all its barcodes)"),
        Stat(step, "mean_barcodes_per_target", number(n_pairs / n_targets if n_targets else 0.0),
             "mean number of final barcodes per final target"),
    ]


# ---- all rows and the file -------------------------------------------------

def collect_stats(con, raw_stats, ingest_stats, args):
    """
    All rows of the stats file, ordered by step.
    Parameters:
            -- con: duckdb connection after assign_targets (tables pairs, final, barcode_assign)
            -- raw_stats: dict: content of raw_stats.json (S2)
            -- ingest_stats: dict or None: content of ingest_stats.json (S1)
            -- args: argparse.Namespace: min_pair_reads, min_purity, merge_dist
    Returns:
            -- list[Stat]: rows
    """
    merged = getattr(args, "merge_dist", 0) > 0
    return (extraction_stats(ingest_stats, raw_stats) + barcode_clustering_stats(con, raw_stats)
            + target_clustering_stats(con, raw_stats, merged) + assignment_stats(con, args)
            + final_map_stats(con))


def write_stats(path, stats):
    """
    Write the stats file (step / stat / count / percent / description) and echo it to the log,
    one block per step.
    Parameters:
            -- path: str: output file
            -- stats: list[Stat]: rows
    Returns:
            -- None: nothing
    """
    # step1 : percent of reads_raw_input for read counts
    total = next((s.count for s in stats if s.stat == "reads_raw_input"), 0)

    def percent(row):
        """
        Percent column of one row.
        Parameters:
                -- row: Stat: row
        Returns:
                -- str: e.g. '74.69', or 'NA'
        """
        return "%.2f" % (100.0 * row.count / total) if row.is_reads and total else "NA"

    # step2 : the file
    with open(path, "w") as fo:
        fo.write("step\tstat\tcount\tpercent\tdescription\n")
        for row in stats:
            fo.write("%s\t%s\t%s\t%s\t%s\n" % (row.step, row.stat, row.count, percent(row),
                                               row.description))
    # step3 : the log, one block per step
    step = None
    for row in stats:
        if row.step != step:
            step = row.step
            log_sub(step)
        log_sub("    %-36s %14s %8s" % (row.stat, row.count, percent(row)))
