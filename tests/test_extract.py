"""
S1 extraction (barcodex/extract.py): none / pos / flank methods and the filters.
"""

import random
import unittest

import edlib

import helpers  # noqa: F401  (puts the package on sys.path)
from barcodex.extract import (ExtractSpec, SequenceExtractor, extract_sequence,
                              match_approximate, match_hamming_numpy, match_levenshtein,
                              matches_template, max_length, parse_allowed_lens,
                              reject_reason_order, validate_spec)

UP, DOWN = "GTCCTCTTCCCATCT", "TGGTGGAAGC"
BARCODE = "ACGTACGTACGTACGTACGT"
TEMPLATE = "NNNNATNNNNATNN"


def locate(spec, seq):
    """Locate in a str read; returns the extracted str or the reject key."""
    span = SequenceExtractor(spec).locate(seq.encode())
    return span if isinstance(span, str) else seq[span[0]:span[1]]


class TestParseAllowedLens(unittest.TestCase):

    def test_list_and_range(self):
        self.assertEqual(parse_allowed_lens("19,20,21"), {19, 20, 21})
        self.assertEqual(parse_allowed_lens("16,18-20"), {16, 18, 19, 20})

    def test_empty(self):
        self.assertEqual(parse_allowed_lens(None), frozenset())

    def test_invalid(self):
        with self.assertRaises(ValueError):
            parse_allowed_lens("0,20")


def first_window(seq, pattern, max_mismatches, distance):
    """Reference: sequence_utils.py window loop (Hamming by hand, Levenshtein by edlib NW)."""
    k = len(pattern)
    for i in range(len(seq) - k + 1):
        window = seq[i:i + k]
        if distance == "hamming":
            if sum(a != b for a, b in zip(window, pattern)) <= max_mismatches:
                return i
        elif edlib.align(window, pattern, mode="NW", k=max_mismatches)["editDistance"] != -1:
            return i
    return -1


class TestMatchHamming(unittest.TestCase):

    def test_exact(self):
        self.assertEqual(match_hamming_numpy(b"AAAA" + UP.encode(), UP.encode(), 2), 4)

    def test_mismatch(self):
        flank = UP[:3] + "A" + UP[4:]
        self.assertEqual(match_hamming_numpy(b"CC" + flank.encode(), UP.encode(), 1), 2)
        self.assertEqual(match_hamming_numpy(b"CC" + flank.encode(), UP.encode(), 0), -1)

    def test_first_window_before_exact(self):
        # the first window within the limit wins, even before an exact match
        flank = UP[:3] + "A" + UP[4:]
        seq = (flank + "CC" + UP).encode()
        self.assertEqual(match_hamming_numpy(seq, UP.encode(), 1), 0)
        self.assertEqual(match_hamming_numpy(seq, UP.encode(), 0), len(UP) + 2)

    def test_short_read(self):
        self.assertEqual(match_hamming_numpy(UP[:5].encode(), UP.encode(), 2), -1)


class TestMatchLevenshtein(unittest.TestCase):

    def test_deletion(self):
        # the window has the flank length: a deletion pulls in one more base = 2 edits
        deleted = UP[:5] + UP[6:]
        seq = "CC" + deleted + "TTTT"
        self.assertEqual(match_levenshtein(seq.encode(), UP.encode(), 1), -1)
        self.assertEqual(match_levenshtein(seq.encode(), UP.encode(), 2),
                         first_window(seq, UP, 2, "levenshtein"))
        self.assertNotEqual(match_levenshtein(seq.encode(), UP.encode(), 2), -1)

    def test_exact_flank_shifted_window(self):
        # same as sequence_utils.py: the window one base before an exact flank differs by 2
        # edits (one deletion, one insertion), so with 2 edits it is the first match
        seq = ("AAAA" + UP + "CC").encode()
        self.assertEqual(match_levenshtein(seq, UP.encode(), 1), 4)
        self.assertEqual(match_levenshtein(seq, UP.encode(), 2), 3)

    def test_short_read(self):
        self.assertEqual(match_levenshtein(UP[:5].encode(), UP.encode(), 2), -1)

    def test_unknown_distance(self):
        with self.assertRaises(ValueError):
            match_approximate(b"ACGT", b"AC", 0, "edit")


class TestSameAsSequenceUtils(unittest.TestCase):

    def test_random_reads(self):
        # random reads with mutated flanks: the same first window as the plain window loop
        rng = random.Random(1)
        for _ in range(400):
            flank = list(UP)
            for _ in range(rng.randint(0, 3)):
                pos = rng.randrange(len(flank))
                op = rng.choice("sdi")
                if op == "s":
                    flank[pos] = rng.choice("ACGT")
                elif op == "d":
                    del flank[pos]
                else:
                    flank.insert(pos, rng.choice("ACGT"))
            seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 80))) + "".join(flank)
            seq += "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 40)))
            for distance in ("hamming", "levenshtein"):
                for mm in (0, 1, 2):
                    self.assertEqual(match_approximate(seq.encode(), UP.encode(), mm, distance),
                                     first_window(seq, UP, mm, distance), (seq, distance, mm))

    def test_extract_sequence(self):
        seq = ("NN" + UP + BARCODE + DOWN + "AC").encode()
        start = 2 + len(UP)
        self.assertEqual(extract_sequence(seq, UP.encode(), DOWN.encode(), 0, "hamming"),
                         (start, start + len(BARCODE)))
        self.assertEqual(extract_sequence(seq, b"TTTTTTTTTTTT", DOWN.encode(), 1, "hamming"),
                         "upstream not found")
        self.assertEqual(extract_sequence(seq, UP.encode(), b"GGGGGGGGGGGG", 1, "levenshtein"),
                         "downstream not found")


class TestMethods(unittest.TestCase):

    def test_none_whole_read(self):
        self.assertEqual(locate(ExtractSpec("barcode"), BARCODE), BARCODE)
        self.assertEqual(locate(ExtractSpec("barcode"), ""), "barcode_empty")

    def test_pos(self):
        self.assertEqual(locate(ExtractSpec("barcode", method="pos", start=2, length=5),
                                BARCODE), BARCODE[2:7])
        self.assertEqual(locate(ExtractSpec("barcode", method="pos", start=4), BARCODE),
                         BARCODE[4:])

    def test_pos_short_read(self):
        barcode = ExtractSpec("barcode", method="pos", start=10, length=20)
        target = ExtractSpec("target", method="pos", start=10, length=20,
                             require_full_len=False)
        self.assertEqual(locate(barcode, BARCODE), "barcode_too_short")
        self.assertEqual(locate(target, BARCODE), BARCODE[10:])

    def test_flank(self):
        spec = ExtractSpec("barcode", method="flank", up=UP, down=DOWN, flank_mismatch=2)
        self.assertEqual(locate(spec, "NN" + UP + BARCODE + DOWN + "AC"), BARCODE)
        self.assertEqual(locate(spec, "A" * 20 + BARCODE + DOWN), "barcode_upstream_not_found")
        self.assertEqual(locate(spec, UP + BARCODE + "C" * 12), "barcode_downstream_not_found")
        self.assertEqual(locate(spec, UP + DOWN), "barcode_empty")

    def test_flank_levenshtein(self):
        spec = ExtractSpec("barcode", method="flank", up=UP, down=DOWN, flank_mismatch=2,
                           flank_distance="levenshtein")
        hamming = ExtractSpec("barcode", method="flank", up=UP, down=DOWN, flank_mismatch=2)
        up_del = UP[:5] + UP[6:]
        read = "NN" + up_del + BARCODE + DOWN + "AC"
        self.assertEqual(locate(hamming, read), "barcode_upstream_not_found")
        start = first_window(read, UP, 2, "levenshtein") + len(UP)
        end = start + first_window(read[start:], DOWN, 2, "levenshtein")
        self.assertEqual(locate(spec, read), read[start:end])
        self.assertEqual(locate(spec, "A" * 20 + BARCODE + DOWN), "barcode_upstream_not_found")
        self.assertEqual(locate(spec, UP + BARCODE + "C" * 12), "barcode_downstream_not_found")

class TestFilters(unittest.TestCase):

    def test_max_len(self):
        spec = ExtractSpec("target", max_len=10)
        self.assertEqual(locate(spec, "A" * 11), "target_too_long")
        self.assertEqual(locate(spec, "A" * 10), "A" * 10)

    def test_min_len(self):
        spec = ExtractSpec("target", min_len=10)
        self.assertEqual(locate(spec, "A" * 9), "target_below_min_len")
        self.assertEqual(locate(spec, "A" * 10), "A" * 10)
        self.assertEqual(locate(spec, "A" * 200), "A" * 200)      # a minimum, not an exact length

    def test_allowed_lens(self):
        spec = ExtractSpec("barcode", allowed_lens=frozenset({19, 21}))
        self.assertEqual(locate(spec, BARCODE), "barcode_length_filtered")
        self.assertEqual(locate(spec, BARCODE[:19]), BARCODE[:19])

    def test_template(self):
        spec = ExtractSpec("barcode", template=TEMPLATE, template_mismatch=1)
        good = "CCGGATTTAAATGG"
        self.assertEqual(locate(spec, good), good)
        one_off = "CCGGCTTTAAATGG"                                  # 1 fixed base wrong: kept
        self.assertEqual(locate(spec, one_off), one_off)
        two_off = "CCGGCTTTAACTGG"                                  # 2 fixed bases wrong
        self.assertEqual(locate(spec, two_off), "barcode_template_mismatch")
        self.assertEqual(locate(spec, good + "A"), "barcode_template_mismatch")   # wrong length

    def test_matches_template_like_reference(self):
        # check_barcode of sequence_utils.py: N positions are free, the others are counted
        fixed = tuple((i, ord(b)) for i, b in enumerate("NAN") if b != "N")
        self.assertTrue(matches_template(b"CAG", 0, 3, fixed, 3, 0))
        self.assertFalse(matches_template(b"CCG", 0, 3, fixed, 3, 0))
        self.assertTrue(matches_template(b"CCG", 0, 3, fixed, 3, 1))


class TestSpecHelpers(unittest.TestCase):

    def test_validate_method_options(self):
        self.assertEqual(validate_spec(ExtractSpec("barcode")), [])
        self.assertTrue(validate_spec(ExtractSpec("barcode", length=16)))
        self.assertTrue(validate_spec(ExtractSpec("barcode", method="flank", up=UP)))
        self.assertTrue(validate_spec(ExtractSpec("target", up=UP, down=DOWN)))

    def test_validate_template(self):
        self.assertEqual(validate_spec(ExtractSpec("barcode", template=TEMPLATE,
                                                   template_mismatch=1)), [])
        for bad in (ExtractSpec("barcode", template="NNXXNN"),                 # bad letters
                    ExtractSpec("barcode", template="ACGT"),                   # no N
                    ExtractSpec("barcode", template="NNAN", template_mismatch=1),
                    ExtractSpec("barcode", template=TEMPLATE,
                                allowed_lens=frozenset({20})),
                    ExtractSpec("barcode", template=TEMPLATE, min_len=20)):
            self.assertTrue(validate_spec(bad), bad)

    def test_validate_min_len(self):
        self.assertTrue(validate_spec(ExtractSpec("target", min_len=-1)))
        self.assertTrue(validate_spec(ExtractSpec("target", min_len=30,
                                                  allowed_lens=frozenset({20}))))

    def test_max_length(self):
        self.assertEqual(max_length(ExtractSpec("barcode"), 150), 150)
        self.assertEqual(max_length(ExtractSpec("barcode", method="pos", start=10), 150), 140)
        flank = ExtractSpec("barcode", method="flank", up=UP, down=DOWN)
        self.assertEqual(max_length(flank, 150), 150 - len(UP) - len(DOWN))
        self.assertEqual(max_length(ExtractSpec("barcode", template=TEMPLATE), 150),
                         len(TEMPLATE))

    def test_reject_reason_order(self):
        barcode = ExtractSpec("barcode", method="flank", up=UP, down=DOWN, template=TEMPLATE)
        target = ExtractSpec("target", method="flank", up=UP, down=DOWN, min_len=170)
        keys = [key for key, _ in reject_reason_order(barcode, target)]
        self.assertEqual(keys, ["target_upstream_not_found", "target_downstream_not_found",
                                "target_empty",
                                "barcode_upstream_not_found", "barcode_downstream_not_found",
                                "barcode_empty",
                                "target_below_min_len", "barcode_template_mismatch"])


if __name__ == "__main__":
    unittest.main()
