"""
S1 ingest: FASTQ streaming (fastq.py), parsing (ingest.py), part files and compaction
(compaction.py), and the whole stage.
"""

import glob
import gzip
import json
import os
import queue
import threading
import unittest

from helpers import TempDirTestCase, make_args, read_parquet, write_fastq
from barcodex.common import meta_path
from barcodex.compaction import Compactor, compaction_sql
from barcodex.constants import TARGET_MASK_BITS
from barcodex.extract import ExtractSpec, SequenceExtractor
from barcodex.fastq import peek_read_len, reader_proc
from barcodex.ingest import (flush_part, ingest_chunk, make_lowq_table, make_sequence_tables,
                             stage_ingest)

UP, DOWN = "GTCCTCTTCCCATCT", "TGGTGGAAGC"


def fastq_bytes(records):
    """FASTQ chunk (bytes) of [(sequence, quality or None)]."""
    return b"".join(b"@r\n%s\n+\n%s\n" % (s.encode(), (q or "I" * len(s)).encode())
                    for s, q in records)


def parse(barcode_records, target_records, barcode_spec, target_spec, use_qual=False):
    """Run ingest_chunk on one chunk pair; returns (n_reads, agg, rejects)."""
    cfg = {"use_barcode_qual": use_qual, "use_target_qual": use_qual}
    tables = make_sequence_tables() + (make_lowq_table(33, 20), make_lowq_table(33, 20))
    agg, rejects = {}, {}
    n = ingest_chunk(fastq_bytes(barcode_records), fastq_bytes(target_records), cfg, agg,
                     tables, (SequenceExtractor(barcode_spec), SequenceExtractor(target_spec)),
                     rejects)
    return n, agg, rejects


class TestFastq(TempDirTestCase):

    def test_peek_read_len_gz(self):
        path = os.path.join(self.tmp, "r.fq.gz")
        with gzip.open(path, "wt") as fo:
            fo.write("@r\nACGTA\n+\nIIIII\n")
        self.assertEqual(peek_read_len(path, use_pigz=False), 5)

    def test_reader_deals_chunks_round_robin(self):
        path = write_fastq(os.path.join(self.tmp, "r.fq"), [("ACGT", None)] * 5)
        queues = [queue.Queue(), queue.Queue()]
        reader_proc([path], 2, False, 0, queues, queue.Queue(), threading.Event())
        sizes = [[], []]
        for i, q in enumerate(queues):
            while True:
                item = q.get_nowait()
                if item is None:
                    break
                sizes[i].append(item[0])
        self.assertEqual(sizes, [[2, 1], [2]])          # chunk 0, 2 -> worker 0; chunk 1 -> 1


class TestIngestChunk(unittest.TestCase):

    def test_counts_identical_pairs(self):
        n, agg, rejects = parse([("AAAA", None)] * 3 + [("CCCC", None)],
                                [("GGGG", None)] * 4, ExtractSpec("barcode"),
                                ExtractSpec("target"))
        self.assertEqual(n, 4)
        self.assertEqual(agg, {b"AAAA\tGGGG": 3, b"CCCC\tGGGG": 1})
        self.assertEqual(rejects, {})

    def test_flank_extraction_and_reject_order(self):
        barcode = ExtractSpec("barcode", method="flank", up=UP, down=DOWN, template="NNATNN",
                              template_mismatch=0)
        target = ExtractSpec("target", min_len=5)
        reads = [
            (UP + "CCATGG" + DOWN, "ACGTAC"),       # kept
            ("T" * 30, "ACG"),                      # target too short AND no barcode flank:
                                                    #   the barcode extraction is checked first
            (UP + "CCATGG" + DOWN, "ACG"),          # target too short (filter stage)
            (UP + "CCGGGG" + DOWN, "ACGTAC"),       # barcode does not match the template
        ]
        n, agg, rejects = parse([(b, None) for b, _ in reads], [(t, None) for _, t in reads],
                                barcode, target)
        self.assertEqual(agg, {b"CCATGG\tACGTAC": 1})
        self.assertEqual(rejects, {"barcode_upstream_not_found": 1, "target_below_min_len": 1,
                                   "barcode_template_mismatch": 1})

    def test_low_quality_masks_are_anded(self):
        # target base 2 is low quality in both reads, base 0 only in the first one
        n, agg, _ = parse([("AC", None)] * 2, [("ACGT", "#I#I"), ("ACGT", "II#I")],
                          ExtractSpec("barcode"), ExtractSpec("target"), use_qual=True)
        value = agg[b"AC\tACGT"]
        self.assertEqual(value & 0xFFFFFFFF, 2)                       # count
        mask = (value >> 32) & ((1 << TARGET_MASK_BITS) - 1)
        self.assertEqual(mask, 0b0010)                                # only base 2 (of 4)


class TestPartsAndCompaction(TempDirTestCase):

    def test_flush_part_then_compact_merges_duplicates(self):
        parts = os.path.join(self.tmp, "parts")
        tuples = os.path.join(self.tmp, "tuples")
        os.makedirs(parts)
        os.makedirs(tuples)
        abort = threading.Event()
        # the same pair in two parts: counts add up, masks are AND-ed
        flush_part({b"AAAA\tGGGG": (0b0110 << 32) | 2}, parts, 0, 0, 0, abort)
        flush_part({b"AAAA\tGGGG": (0b0011 << 32) | 3, b"CCCC\tTTTT": 1}, parts, 1, 0, 0, abort)
        self.assertEqual(len(glob.glob(parts + "/*.tsv.gz")), 2)
        compactor = Compactor(self.tmp, 1, 1, batch=10, fanout=10)
        compactor.finish()
        self.assertEqual(glob.glob(parts + "/*.tsv.gz"), [])
        rows = read_parquet(tuples + "/t*.parquet",
                            "SELECT barcode, target, cnt, tm0 FROM t ORDER BY barcode")
        self.assertEqual(rows, [("AAAA", "GGGG", 5, 0b0010), ("CCCC", "TTTT", 1, 0)])

    def test_compaction_sql_inputs(self):
        sql = compaction_sql(["p.tsv.gz"], ["a.parquet"], "out.parquet")
        self.assertIn("read_csv(['p.tsv.gz']", sql)
        self.assertIn("read_parquet(['a.parquet'])", sql)
        self.assertIn("UNION ALL", sql)


class TestStageIngest(TempDirTestCase):

    def test_stage_writes_tuples_and_stats(self):
        barcodes = ["AAAACCCC", "AAAACCCC", "GGGGTTTT", "AC"]
        targets = ["ACGTACGTAA", "ACGTACGTAA", "TTTTGGGGCC", "ACGTACGTAA"]
        args = make_args(os.path.join(self.tmp, "out"),
                         barcode_read=[write_fastq(os.path.join(self.tmp, "b.fq"),
                                                   [(b, None) for b in barcodes])],
                         target_read=[write_fastq(os.path.join(self.tmp, "t.fq"),
                                                  [(t, None) for t in targets])],
                         barcode_allowed_lens="8", no_pigz=True, chunk=2)
        work = os.path.join(self.tmp, "work")
        os.makedirs(work)
        self.assertEqual(stage_ingest(args, work, 1, 1.0), 4)
        rows = read_parquet(work + "/tuples/t*.parquet",
                            "SELECT barcode, target, cnt FROM t ORDER BY barcode")
        self.assertEqual(rows, [("AAAACCCC", "ACGTACGTAA", 2), ("GGGGTTTT", "TTTTGGGGCC", 1)])
        with open(meta_path(args, "ingest_stats.json")) as fh:
            stats = json.load(fh)
        self.assertEqual((stats["total_read_pairs"], stats["extracted_read_pairs"]), (4, 3))
        self.assertIn(["barcode_length_filtered", 1, stats["rejects"][-1][2]], stats["rejects"])
        self.assertTrue(os.path.exists(work + "/ingest.done"))


if __name__ == "__main__":
    unittest.main()
