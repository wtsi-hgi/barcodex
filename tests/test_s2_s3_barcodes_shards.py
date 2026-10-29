"""
S2 barcode clustering (barcodes.py) and S3 sharding (shards.py).
"""

import json
import os
import unittest

from helpers import TempDirTestCase, make_args, read_parquet, write_tuples
from barcodex.barcodes import packed_barcode_sql, rank_key_sql
from barcodex.barcodes import stage_barcodes
from barcodex.common import meta_path
from barcodex.shards import stage_shards

HEAVY = "ACGTACGTACGTACGT"
CHILD = "ACGTACGTACGTACGA"          # 1 mismatch to HEAVY
OTHER = "TTTTGGGGCCCCAAAA"
LONELY = "GATCGATCGATCGATC"         # 1 read, no parent -> dropped (< --min_barcode_reads 2)


class TestBarcodeSql(unittest.TestCase):

    def test_packed_barcode_sql_has_one_term_per_base(self):
        code, conf = packed_barcode_sql(4, 4)
        self.assertEqual(code.count("substr(b,"), 4)
        self.assertEqual(conf.count("bm >>"), 4)

    def test_rank_key_prefers_low_cost(self):
        self.assertIn("79228162514264337593543950336", rank_key_sql("c", "n", "h"))


class TestStageBarcodes(TempDirTestCase):

    def test_child_joins_heavy_parent(self):
        work = os.path.join(self.tmp, "work")
        os.makedirs(work + "/tuples")
        write_tuples(work + "/tuples/t1_000001.parquet",
                     [(HEAVY, "TGT1", 20), (CHILD, "TGT1", 2), (OTHER, "TGT2", 10),
                      (LONELY, "TGT3", 1)])
        args = make_args(os.path.join(self.tmp, "out"))
        stage_barcodes(args, work, 1, 1.0)
        bcmap = dict(read_parquet(work + "/bcmap.parquet", "SELECT barcode, centre FROM t"))
        self.assertEqual(bcmap, {HEAVY: HEAVY, CHILD: HEAVY, OTHER: OTHER})
        clusters = read_parquet(meta_path(args, "barcode_clusters.parquet"),
                                "SELECT barcode, count, n_cluster FROM t ORDER BY barcode")
        self.assertEqual(clusters, [(HEAVY, 22, 2), (OTHER, 10, 1)])
        with open(meta_path(args, "raw_stats.json")) as fh:
            raw = json.load(fh)
        self.assertEqual((raw["raw_barcodes"], raw["reads_in"], raw["reads_clustered"]),
                         (4, 33, 32))


class TestStageShards(TempDirTestCase):

    def test_rows_are_rekeyed_by_centre(self):
        work = os.path.join(self.tmp, "work")
        os.makedirs(work + "/tuples")
        write_tuples(work + "/tuples/t1_000001.parquet",
                     [(HEAVY, "TGT1", 20, 0b0110), (CHILD, "TGT1", 2, 0b0011),
                      (OTHER, "TGT2", 10)])
        import duckdb
        con = duckdb.connect()
        con.execute("COPY (SELECT * FROM (VALUES ('%s','%s',20), ('%s','%s',2), ('%s','%s',10)) "
                    "v(barcode, centre, n)) TO '%s/bcmap.parquet' (FORMAT PARQUET)"
                    % (HEAVY, HEAVY, CHILD, HEAVY, OTHER, OTHER, work))
        con.close()
        args = make_args(os.path.join(self.tmp, "out"))
        stage_shards(args, work, 1, 1.0, n_shards=2)
        rows = read_parquet(work + "/shards/*/*.parquet",
                            "SELECT centre, target, cnt, tm0 FROM t ORDER BY centre")
        # CHILD merged into HEAVY: counts add up, target masks AND-ed
        self.assertEqual(rows, [(HEAVY, "TGT1", 22, 0b0010), (OTHER, "TGT2", 10, 0)])
        self.assertFalse(os.path.exists(work + "/tuples"))      # freed without --keep_tmp
        self.assertTrue(os.path.exists(work + "/shards.done"))


if __name__ == "__main__":
    unittest.main()
