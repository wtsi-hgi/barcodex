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
        barcodes = read_tsv(args.out + ".barcode_target.details.tsv")
        stats = read_tsv(args.out + ".stats.tsv")
        return barcodes, stats

    def test_one_target_per_barcode(self):
        barcodes, _ = self.run_resolve()
        self.assertEqual(sorted((r["barcode"], r["target"]) for r in barcodes),
                         [("BC_RESOLVED", "TGT_A"), ("BC_UNIQUE", "TGT_A")])

    def test_output_columns(self):
        barcodes, _ = self.run_resolve()
        self.assertEqual(list(barcodes[0]), ["barcode", "target",
                                             "barcode_reads", "n_diff_barcodes", "first_barcode_cov",
                                             "first_barcode_ratio", "second_barcode_cov",
                                             "second_barcode_ratio",
                                             "target_reads", "n_diff_targets", "first_target_cov",
                                             "first_target_ratio", "second_target_cov",
                                             "second_target_ratio"])
        # only the details file and the stats file are written
        self.assertEqual(sorted(f for f in os.listdir(self.tmp) if f.endswith(".tsv")),
                         ["out.barcode_target.details.tsv", "out.stats.tsv"])

    def test_pair_statistics(self):
        """BC_RESOLVED: barcode side from barcode_clusters.parquet, target side from its top target."""
        barcodes, _ = self.run_resolve()
        row = {r["barcode"]: r for r in barcodes}["BC_RESOLVED"]
        self.assertEqual({k: row[k] for k in ("barcode_reads", "n_diff_barcodes", "first_barcode_cov",
                                              "first_barcode_ratio", "second_barcode_cov",
                                              "target_reads", "n_diff_targets", "first_target_cov",
                                              "first_target_ratio", "second_target_cov")},
                         {"barcode_reads": "20", "n_diff_barcodes": "1", "first_barcode_cov": "20",
                          "first_barcode_ratio": "1.0", "second_barcode_cov": "0",
                          "target_reads": "18", "n_diff_targets": "1", "first_target_cov": "18",
                          "first_target_ratio": "1.0", "second_target_cov": "0"})

    def test_ratios_of_cluster_members(self):
        """
        Barcode 1 = raw barA 6, barB 2, barC 1, barD 1 reads; its target = 3 raw variants 7, 2, 1.
        first_/second_barcode_ratio = 6/10, 2/10; first_/second_target_ratio = 7/10, 2/10.
        """
        args = make_args(os.path.join(self.tmp, "out"))
        os.makedirs(meta_path(args))
        with open(meta_path(args, "target_pairs.tsv"), "w") as fo:
            fo.write(PAIR_HEADER)
            fo.write("BC_1\t10\t1\t10\t3\tTGT_A\t7\t2\n")
        con = duckdb.connect()
        con.execute("COPY (SELECT 'BC_1' AS barcode, 10 AS count, 4 AS n_cluster, 6 AS first_cov, "
                    "2 AS second_cov) TO '%s' (FORMAT PARQUET)" % meta_path(args, "barcode_clusters.parquet"))
        con.close()
        with open(meta_path(args, "raw_stats.json"), "w") as fo:
            json.dump({"raw_barcodes": 4, "raw_targets": 3, "reads_in": 10, "reads_clustered": 10}, fo)
        resolve_pairs(args, 1)
        row = read_tsv(args.out + ".barcode_target.details.tsv")[0]
        self.assertEqual(row, {"barcode": "BC_1", "target": "TGT_A",
                               "barcode_reads": "10", "n_diff_barcodes": "4",
                               "first_barcode_cov": "6", "first_barcode_ratio": "0.6",
                               "second_barcode_cov": "2", "second_barcode_ratio": "0.2",
                               "target_reads": "10", "n_diff_targets": "3",
                               "first_target_cov": "7", "first_target_ratio": "0.7",
                               "second_target_cov": "2", "second_target_ratio": "0.2"})

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


class TestResolveMerge(TempDirTestCase):

    def test_merge_targets_across_barcodes(self):
        """--merge_dist 1: the 1-mismatch variant ACGTACGTAA joins the more abundant ACGTACGTAC."""
        args = make_args(os.path.join(self.tmp, "out"), merge_dist=1, merge_ratio=1.0)
        os.makedirs(meta_path(args))
        # barcode, barcode_reads, rank, target_reads, target_unique, target, first_cov, second_cov
        rows = [("BC_1", 13, 1, 10, 1, "ACGTACGTAC", 10, 0),
                ("BC_1", 13, 2, 3, 2, "ACGTACGTAA", 2, 1),
                ("BC_2", 4, 1, 4, 1, "ACGTACGTAA", 4, 0)]
        with open(meta_path(args, "target_pairs.tsv"), "w") as fo:
            fo.write(PAIR_HEADER)
            for r in rows:
                fo.write("%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\n" % r)
        con = duckdb.connect()
        con.execute("COPY (SELECT * FROM (VALUES ('BC_1', 13, 1, 13, 0), ('BC_2', 4, 1, 4, 0)) "
                    "v(barcode, count, n_cluster, first_cov, second_cov)) TO '%s' (FORMAT PARQUET)"
                    % meta_path(args, "barcode_clusters.parquet"))
        con.close()
        with open(meta_path(args, "raw_stats.json"), "w") as fo:
            json.dump({"raw_barcodes": 2, "raw_targets": 2, "reads_in": 17, "reads_clustered": 17}, fo)
        resolve_pairs(args, 1)
        targets = {r["barcode"]: r for r in read_tsv(args.out + ".barcode_target.details.tsv")}
        self.assertEqual(sorted(targets), ["BC_1", "BC_2"])
        # BC_1: both clusters merged -> 13 reads, 3 variants, top two coverages 10 and 2
        self.assertEqual({k: targets["BC_1"][k] for k in ("target", "target_reads", "n_diff_targets",
                                                          "first_target_cov", "first_target_ratio",
                                                          "second_target_cov", "second_target_ratio")},
                         {"target": "ACGTACGTAC", "target_reads": "13", "n_diff_targets": "3",
                          "first_target_cov": "10", "first_target_ratio": "0.7692",
                          "second_target_cov": "2", "second_target_ratio": "0.1538"})
        # BC_2 only saw the variant: it is reported as the root sequence
        self.assertEqual((targets["BC_2"]["target"], targets["BC_2"]["target_reads"]),
                         ("ACGTACGTAC", "4"))


if __name__ == "__main__":
    unittest.main()
