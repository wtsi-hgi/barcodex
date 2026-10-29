"""
S4 target clustering: edit distance (distance.py) and clustering within a barcode (targets.py).
"""

import io
import os
import unittest

from helpers import TempDirTestCase, make_args, read_tsv
from barcodex.common import meta_path
from barcodex.distance import mate_cost, unpack_seq
from barcodex.targets import init_s4_worker, run_groups, stage_targets

TARGET = "ACGTACGTACGTACGTACGT"
VARIANT = "ACGTACGTACCTACGTACGT"        # 1 substitution at position 10
DELETION = "ACGTACGTACTACGTACGT"        # 1 base deleted
OTHER = "TTTTGGGGCCCCAAAATTTT"

# (max_hard, max_low, band, ratio, low_ratio, use_quality)
PARAMS = (2, 5, 1, 5.0, 1.0, True)


def cluster(rows):
    """run_groups on (centre, target, cnt, tm0) rows; returns the output lines as field lists."""
    init_s4_worker({"max_parents": 200, "max_dp_group": 1000})
    out = io.StringIO()
    run_groups([(c, t, n, m, 0, 0, 0) for c, t, n, m in rows], out, None, None, PARAMS)
    return [line.split("\t") for line in out.getvalue().splitlines()]


class TestDistance(unittest.TestCase):

    def test_hamming_hard_and_low(self):
        confident = b"\xff" * len(TARGET)
        self.assertEqual(mate_cost(TARGET.encode(), confident, VARIANT.encode(), confident,
                                   2, 5, 0, 1), (1, 0))
        low = bytearray(confident)
        low[10] = 0
        self.assertEqual(mate_cost(TARGET.encode(), confident, VARIANT.encode(), bytes(low),
                                   2, 5, 0, 1), (0, 1))

    def test_indel_needs_band(self):
        qual_a, qual_b = b"\xff" * len(TARGET), b"\xff" * len(DELETION)
        self.assertIsNone(mate_cost(TARGET.encode(), qual_a, DELETION.encode(), qual_b,
                                    2, 5, 0, 1))
        self.assertEqual(mate_cost(TARGET.encode(), qual_a, DELETION.encode(), qual_b,
                                   2, 5, 1, 1), (1, 0))

    def test_unpack_seq(self):
        packed = int("0123", 4).to_bytes(1, "big")
        self.assertEqual(unpack_seq(packed, 4), b"ACGT")


class TestClustering(unittest.TestCase):

    def test_error_variants_join_the_main_target(self):
        lines = cluster([("BC1", TARGET, 50, 0), ("BC1", VARIANT, 3, 0), ("BC1", DELETION, 2, 0)])
        self.assertEqual(len(lines), 1)
        barcode, barcode_reads, rank, reads, n_unique, seq, first, second = lines[0]
        self.assertEqual((barcode, seq, reads, n_unique, first, second),
                         ("BC1", TARGET, "55", "3", "50", "3"))

    def test_distinct_targets_stay_apart(self):
        lines = cluster([("BC1", TARGET, 50, 0), ("BC1", OTHER, 40, 0), ("BC2", OTHER, 7, 0)])
        self.assertEqual([(l[0], l[2], l[5]) for l in lines],
                         [("BC1", "1", TARGET), ("BC1", "2", OTHER), ("BC2", "1", OTHER)])

    def test_coverage_ratio(self):
        # 1 hard mismatch needs parent >= 5 x child: 10 vs 3 is not enough
        lines = cluster([("BC1", TARGET, 10, 0), ("BC1", VARIANT, 3, 0)])
        self.assertEqual(len(lines), 2)
        # the same mismatch at a low-quality base only needs low_ratio 1
        low_bit = 1 << (len(VARIANT) - 1 - 10)
        lines = cluster([("BC1", TARGET, 10, 0), ("BC1", VARIANT, 3, low_bit)])
        self.assertEqual(len(lines), 1)


class TestStageTargets(TempDirTestCase):

    def test_stage_writes_target_pairs(self):
        work = os.path.join(self.tmp, "work")
        os.makedirs(work + "/shards/shard=0")
        import duckdb
        con = duckdb.connect()
        con.execute("COPY (SELECT col0 AS centre, col1 AS target, col2::UBIGINT AS cnt, "
                    "0::UBIGINT AS tm0, 0::UBIGINT AS tm1, 0::UBIGINT AS tm2, 0::UBIGINT AS tm3 "
                    "FROM (VALUES ('BC1', '%s', 50), ('BC1', '%s', 3))) "
                    "TO '%s/shards/shard=0/data.parquet' (FORMAT PARQUET)"
                    % (TARGET, VARIANT, work))
        con.close()
        args = make_args(os.path.join(self.tmp, "out"))
        stage_targets(args, work, 1, 1.0, n_shards=1)
        rows = read_tsv(meta_path(args, "target_pairs.tsv"))
        self.assertEqual([(r["barcode"], r["target_seq"], r["target_reads"]) for r in rows],
                         [("BC1", TARGET, "53")])


if __name__ == "__main__":
    unittest.main()
