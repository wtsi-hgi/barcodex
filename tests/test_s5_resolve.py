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


def write_meta(args, with_ingest_stats=True, pairs=None):
    """Write the four S5 input files for PAIRS (40 extracted reads, 3 dropped in S2)."""
    pairs = PAIRS if pairs is None else pairs
    os.makedirs(meta_path(args))
    with open(meta_path(args, "target_pairs.tsv"), "w") as fo:
        fo.write(PAIR_HEADER)
        for barcode, (barcode_reads, targets) in pairs.items():
            for rank, (target, reads) in enumerate(targets, 1):
                fo.write("%s\t%d\t%d\t%d\t1\t%s\t%d\t0\n"
                         % (barcode, barcode_reads, rank, reads, target, reads))
    con = duckdb.connect()
    values = ", ".join("('%s', %d, 1, %d, 0)" % (b, r, r) for b, (r, _) in pairs.items())
    con.execute("COPY (SELECT * FROM (VALUES %s) v(barcode, count, n_cluster, first_cov, "
                "second_cov)) TO '%s' (FORMAT PARQUET)"
                % (values, meta_path(args, "barcode_clusters.parquet")))
    con.close()
    with open(meta_path(args, "raw_stats.json"), "w") as fo:
        json.dump({"raw_barcodes": 9, "raw_targets": 5, "raw_barcode_target_pairs": 12,
                   "reads_in": 44, "reads_clustered": 41}, fo)
    if with_ingest_stats:
        with open(meta_path(args, "ingest_stats.json"), "w") as fo:
            json.dump({"total_read_pairs": 50, "extracted_read_pairs": 44,
                       "rejects": [["barcode_upstream_not_found", 6, "no flank"],
                                   ["target_empty", 0, "empty target"]]}, fo)


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
        self.assertEqual(list(rows[0]), ["step", "stat", "count", "percent", "description"])
        steps = [r["step"] for r in rows]
        self.assertEqual(steps, sorted(steps))                      # ordered by step
        stats = {r["stat"]: (r["count"], r["percent"]) for r in rows}
        # target reject rows come before barcode reject rows
        names = [r["stat"] for r in rows]
        self.assertLess(names.index("reads_target_missing"),
                        names.index("reads_barcode_upstream_not_found"))
        expected = {"reads_raw_input": ("50", "100.00"),
                    "reads_target_missing": ("0", "0.00"),
                    "reads_barcode_upstream_not_found": ("6", "12.00"),
                    "num_raw_barcodes": ("9", "NA"), "num_raw_targets": ("5", "NA"),
                    "num_raw_barcode_target_pairs": ("12", "NA"),
                    "reads_dropped_barcode_clustering": ("3", "6.00"),
                    "reads_in_barcode_clustering": ("41", "82.00"),
                    "num_barcodes_after_clustering": ("4", "NA"),
                    "reads_dropped_target_clustering": ("0", "0.00"),
                    "reads_in_target_clustering": ("41", "82.00"),
                    "num_target_after_clustering": ("3", "NA"),
                    "reads_barcodes_multi_targets": ("30", "60.00"),
                    "num_barcodes_multi_targets": ("2", "NA"),
                    "reads_barcodes_resolved": ("20", "40.00"),
                    "num_barcodes_resolved": ("1", "NA"),
                    "reads_barcodes_ambiguous": ("10", "20.00"),
                    "num_barcodes_ambiguous": ("1", "NA"),
                    "reads_barcodes_low_coverage": ("1", "2.00"),
                    "num_barcodes_low_coverage": ("1", "NA"),
                    "num_final_barcodes": ("2", "NA"), "num_final_targets": ("1", "NA"),
                    "num_final_barcode_target_pairs": ("2", "NA"),
                    "reads_final_effective": ("28", "56.00"),
                    "mean_barcodes_per_target": ("2", "NA")}
        self.assertEqual({k: stats[k] for k in expected}, expected)
        for old_name in ("total_read_pairs", "extracted_read_pairs", "effective_reads",
                         "barcodes_unique", "reads_minor_targets"):
            self.assertNotIn(old_name, stats)

    def test_older_meta_without_ingest_stats(self):
        _, rows = self.run_resolve(with_ingest_stats=False)
        stats = {r["stat"]: (r["count"], r["percent"]) for r in rows}
        self.assertNotIn("reads_raw_input", stats)
        self.assertEqual(stats["reads_final_effective"], ("28", "NA"))

    def test_number_format(self):
        self.assertEqual((number(None), number(3.0), number(2.345)), ("0", "3", "2.35"))


class TestMinPurityReads(TempDirTestCase):

    PAIRS = {
        # 8 / 12 = 0.67 with all reads; without the four 1-read targets 8 / 8 = 1.0
        "BC_NOISY": (12, [("TGT_A", 8), ("TGT_B", 1), ("TGT_C", 1), ("TGT_D", 1),
                          ("TGT_E", 1)]),
        # 8 / 12 = 0.67; without the 1-read targets 8 / 10 = 0.8; without < 3 reads 8 / 8
        "BC_TWO":   (12, [("TGT_A", 8), ("TGT_B", 2), ("TGT_C", 1), ("TGT_D", 1)]),
        # 2 / 3 with all reads; 2 / 2 with 2; with 3 the top target (2 reads, below 3) still
        # counts, so 2 / 2 and not 2 / 0
        "BC_TOP2":  (3, [("TGT_A", 2), ("TGT_B", 1)]),
    }

    def final(self, min_purity_reads):
        """Final barcodes for one --min_purity_reads (--min_purity 0.8, --min_pair_reads 2)."""
        args = make_args(os.path.join(self.tmp, "mpr%d" % min_purity_reads, "out"),
                         min_purity_reads=min_purity_reads)
        os.makedirs(os.path.dirname(args.out))
        write_meta(args, pairs=self.PAIRS)
        resolve_pairs(args, 1)
        rows = read_tsv(args.out + ".barcode_target.details.tsv")
        stats = {r["stat"]: r for r in read_tsv(args.out + ".stats.tsv")}
        return sorted(r["barcode"] for r in rows), stats

    def test_default_uses_all_reads(self):
        barcodes, stats = self.final(1)
        self.assertEqual(barcodes, [])
        self.assertEqual(stats["num_barcodes_ambiguous"]["count"], "3")
        self.assertNotIn("--min_purity_reads", stats["num_barcodes_ambiguous"]["description"])

    def test_ignore_one_read_targets(self):
        barcodes, stats = self.final(2)
        self.assertEqual(barcodes, ["BC_NOISY", "BC_TOP2", "BC_TWO"])
        self.assertIn("--min_purity_reads", stats["num_barcodes_resolved"]["description"])

    def test_top_target_always_counts(self):
        barcodes, _ = self.final(3)
        self.assertEqual(barcodes, ["BC_NOISY", "BC_TOP2", "BC_TWO"])


if __name__ == "__main__":
    unittest.main()
