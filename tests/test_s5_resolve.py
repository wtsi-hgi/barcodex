"""
S5 resolve (resolve.py) and the stats file (stats.py), from hand-made <out>.meta/ files.
"""

import json
import os
import unittest

import duckdb

from helpers import TempDirTestCase, make_args, read_tsv
from barcodex.common import meta_path
from barcodex.resolve import resolve_pairs
from barcodex.stats import number

PAIR_HEADER = ("barcode\tbarcode_reads\ttarget_rank\ttarget_reads\ttarget_unique\ttarget_seq\t"
               "first_cov\tsecond_cov\n")

# barcode cluster -> (barcode_reads, [(target, reads)])
PAIRS = {
    "BC_UNIQUE":   (10, [("TGT_A", 10)]),
    "BC_RESOLVED": (20, [("TGT_A", 18), ("TGT_B", 2)]),      # purity 0.9 -> kept
    "BC_AMBIG":    (10, [("TGT_B", 6), ("TGT_C", 4)]),       # purity 0.6 -> dropped
    "BC_LOWCOV":   (1, [("TGT_C", 1)]),                      # < --min_pair_reads 2 -> dropped
}


def write_meta(args, with_ingest_stats=True):
    """Write the four S5 input files for PAIRS (40 extracted reads, 3 dropped in S2)."""
    os.makedirs(meta_path(args))
    with open(meta_path(args, "target_pairs.tsv"), "w") as fo:
        fo.write(PAIR_HEADER)
        for barcode, (barcode_reads, targets) in PAIRS.items():
            for rank, (target, reads) in enumerate(targets, 1):
                fo.write("%s\t%d\t%d\t%d\t1\t%s\t%d\t0\n"
                         % (barcode, barcode_reads, rank, reads, target, reads))
    con = duckdb.connect()
    values = ", ".join("('%s', %d, 1, %d, 0)" % (b, r, r) for b, (r, _) in PAIRS.items())
    con.execute("COPY (SELECT * FROM (VALUES %s) v(barcode, count, n_cluster, first_cov, "
                "second_cov)) TO '%s' (FORMAT PARQUET)"
                % (values, meta_path(args, "barcode_clusters.parquet")))
    con.close()
    with open(meta_path(args, "raw_stats.json"), "w") as fo:
        json.dump({"raw_barcodes": 9, "raw_targets": 5, "reads_in": 44, "reads_clustered": 41},
                  fo)
    if with_ingest_stats:
        with open(meta_path(args, "ingest_stats.json"), "w") as fo:
            json.dump({"total_read_pairs": 50, "extracted_read_pairs": 44,
                       "rejects": [["barcode_upstream_not_found", 6, "no flank"]]}, fo)


class TestResolve(TempDirTestCase):

    def run_resolve(self, with_ingest_stats=True):
        """Run S5 and return (final barcode rows, stats dict)."""
        args = make_args(os.path.join(self.tmp, "out"))
        write_meta(args, with_ingest_stats)
        resolve_pairs(args, 1)
        barcodes = read_tsv(args.out + ".barcode_clusters.tsv")
        stats = read_tsv(args.out + ".stats.tsv")
        return barcodes, stats

    def test_one_target_per_barcode(self):
        barcodes, _ = self.run_resolve()
        self.assertEqual(sorted((r["barcode"], r["target"]) for r in barcodes),
                         [("BC_RESOLVED", "TGT_A"), ("BC_UNIQUE", "TGT_A")])

    def test_output_columns(self):
        barcodes, _ = self.run_resolve()
        self.assertEqual(list(barcodes[0]), ["barcode", "target", "barcode_reads",
                                             "n_diff_barcodes", "first_barcode_cov",
                                             "first_barcode_ratio", "second_barcode_cov",
                                             "second_barcode_ratio"])
        targets = read_tsv(os.path.join(self.tmp, "out.target_clusters.tsv"))
        self.assertEqual(list(targets[0]), ["barcode", "target", "target_reads",
                                            "n_diff_targets", "first_target_cov",
                                            "first_target_ratio", "second_target_cov",
                                            "second_target_ratio"])

    def test_stats_by_step_and_read_funnel(self):
        _, rows = self.run_resolve()
        self.assertEqual(list(rows[0]), ["step", "stat", "value", "percent", "description"])
        steps = [r["step"] for r in rows]
        self.assertEqual(steps, sorted(steps))                      # ordered by step
        stats = {r["stat"]: (r["value"], r["percent"]) for r in rows}
        expected = {"total_read_pairs": ("50", "100.00"),
                    "barcode_upstream_not_found": ("6", "12.00"),
                    "extracted_read_pairs": ("44", "88.00"),
                    "reads_dropped_barcode_clustering": ("3", "6.00"),
                    "reads_low_coverage_barcodes": ("1", "2.00"),
                    "reads_ambiguous_barcodes": ("10", "20.00"),
                    "reads_minor_targets": ("2", "4.00"),
                    "effective_reads": ("28", "56.00"),
                    "barcode_clusters": ("4", "NA"), "barcodes_unique": ("1", "NA"),
                    "barcodes_multi_target": ("2", "NA"), "barcodes_resolved": ("1", "NA"),
                    "barcodes_ambiguous": ("1", "NA"), "barcodes_low_coverage": ("1", "NA"),
                    "barcode_target_variants": ("6", "NA"),
                    "barcode_target_clusters": ("6", "NA"),
                    "final_barcodes": ("2", "NA"), "final_targets": ("1", "NA"),
                    "mean_barcodes_per_target": ("2", "NA")}
        self.assertEqual({k: stats[k] for k in expected}, expected)

    def test_older_meta_without_ingest_stats(self):
        _, rows = self.run_resolve(with_ingest_stats=False)
        stats = {r["stat"]: (r["value"], r["percent"]) for r in rows}
        self.assertNotIn("total_read_pairs", stats)
        self.assertEqual(stats["effective_reads"], ("28", "NA"))

    def test_number_format(self):
        self.assertEqual((number(None), number(3.0), number(2.345)), ("0", "3", "2.35"))


if __name__ == "__main__":
    unittest.main()
