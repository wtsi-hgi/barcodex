"""
S1 extraction (barcodex/extract.py): none / pos / flank methods and the filters.
"""

import random
import unittest


import helpers  # noqa: F401  (puts the package on sys.path)
from barcodex.extract import (ExtractSpec, SequenceExtractor, extract_sequence, find_flank,
                              match_edlib, matches_template, max_length,
                              parse_allowed_lens, reject_reason_order, template_regex,
                              validate_spec)

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


def mutate(rng, flank, n_edits, ops, margin=1):
    """Apply n_edits random edits (ops: s = substitution, d = deletion, i = insertion), keeping
    `margin` bases at each end unchanged."""
    flank = list(flank)
    for _ in range(n_edits):
        pos = rng.randrange(margin, len(flank) - margin)
        op = rng.choice(ops)
        if op == "s":
            flank[pos] = rng.choice([b for b in "ACGT" if b != flank[pos]])
        elif op == "d":
            del flank[pos]
        else:
            flank.insert(pos, rng.choice("ACGT"))
    return "".join(flank)


class TestMatchEdlib(unittest.TestCase):

    def test_exact(self):
        seq = ("AAAA" + UP + "CC").encode()
        for k in (0, 1, 2):
            self.assertEqual(match_edlib(seq, UP.encode(), k), (4, 4 + len(UP), 0))

    def test_deletion_and_insertion(self):
        # one indel = one edit; the match ends exactly at the end of the (mutated) flank
        deleted = UP[:5] + UP[6:]
        inserted = UP[:5] + "A" + UP[5:]
        for flank in (deleted, inserted):
            seq = ("CC" + flank + "TTTT").encode()
            self.assertEqual(match_edlib(seq, UP.encode(), 1), (2, 2 + len(flank), 1))
            self.assertIsNone(match_edlib(seq, UP.encode(), 0))

    def test_substitution_preferred_on_ties(self):
        # last flank base substituted: 'insertion + match' ties with 'substitution'; the
        # substitution (flank length) is chosen, so the next base is not swallowed
        flank = UP[:-1] + ("A" if UP[-1] != "A" else "C")
        seq = ("GG" + flank + UP[-1] + "CCCC").encode()
        self.assertEqual(match_edlib(seq, UP.encode(), 2), (2, 2 + len(UP), 1))

    def test_substitution_preferred_when_edlib_reports_the_insertion(self):
        # combi_1: every read of a barcode has ...TTTTTGG ATAAGTTC... with the downstream flank
        # GTAAGTTC...; 'G>A' and 'flank G + inserted A' both cost 1 and end at the same base, and
        # edlib reports only the insertion; the substitution (flank length) must win
        down = b"GTAAGTTCTTGCTTTGTTCA"
        seq = b"TTTTTGG" + b"A" + down[1:] + b"AACTG"
        self.assertEqual(match_edlib(seq, down, 2), (7, 7 + len(down), 1))

    def test_short_read(self):
        self.assertIsNone(match_edlib(UP[:5].encode(), UP.encode(), 2))


class TestFindFlank(unittest.TestCase):

    def test_exact_first_then_edlib(self):
        # exact flank -> found by find; otherwise the edlib alignment
        rng = random.Random(2)
        for _ in range(400):
            seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 60)))
            seq += mutate(rng, UP, rng.randint(0, 3), "sdi")
            seq = (seq + "".join(rng.choice("ACGT") for _ in range(30))).encode()
            for mm in (0, 1, 2):
                found = find_flank(seq, UP.encode(), mm)
                i = seq.find(UP.encode())
                if i >= 0:
                    self.assertEqual(found, (i, i + len(UP)))
                elif mm == 0:
                    self.assertIsNone(found)
                else:
                    edit = match_edlib(seq, UP.encode(), mm)
                    self.assertEqual(found, None if edit is None else edit[:2])

    def test_one_substitution_keeps_the_flank_length(self):
        rng = random.Random(6)
        for _ in range(400):
            seq = ("".join(rng.choice("ACGT") for _ in range(10)) + mutate(rng, UP, 1, "s")
                   + "".join(rng.choice("ACGT") for _ in range(20))).encode()
            self.assertEqual(find_flank(seq, UP.encode(), 2), (10, 10 + len(UP)))

    def test_indel_next_to_the_barcode(self):
        # deletion of the second-to-last base of the upstream flank: a Hamming window would
        # accept it with 2 mismatches and shift the barcode by one base; edlib needs 1 edit
        read = ("GG" + UP[:-2] + UP[-1] + BARCODE + DOWN).encode()
        self.assertEqual(find_flank(read, UP.encode(), 2), (2, 1 + len(UP)))


class TestExtractSequence(unittest.TestCase):

    def test_extract_sequence(self):
        seq = ("NN" + UP + BARCODE + DOWN + "AC").encode()
        start = 2 + len(UP)
        for mm in (0, 1, 2):
            self.assertEqual(extract_sequence(seq, UP.encode(), DOWN.encode(), mm),
                             (start, start + len(BARCODE)), mm)
        self.assertEqual(extract_sequence(seq, b"TTTTTTTTTTTT", DOWN.encode(), 1),
                         "upstream not found")
        self.assertEqual(extract_sequence(seq, UP.encode(), b"GGGGGGGGGGGG", 1),
                         "downstream not found")

    def test_one_indel_in_each_flank(self):
        # flanks with one indel >= 4 bases away from the barcode: the barcode comes out exactly
        # (an indel right at the junction, or 2 edits, can have another explanation with the
        # same number of edits, e.g. insertion = substitution + 1-base shift)
        rng = random.Random(3)
        for _ in range(300):
            barcode = "".join(rng.choice("ACGT") for _ in range(20))
            up = mutate(rng, UP, 1, "di", margin=4)
            down = mutate(rng, DOWN, 1, "di", margin=4)
            read = ("".join(rng.choice("ACGT") for _ in range(5)) + up + barcode + down
                    + "".join(rng.choice("ACGT") for _ in range(5)))
            span = extract_sequence(read.encode(), UP.encode(), DOWN.encode(), 2)
            self.assertNotIsInstance(span, str, read)
            self.assertEqual(read[span[0]:span[1]], barcode, (up, down))


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

    def test_flank_with_indel(self):
        up_del = UP[:5] + UP[6:]
        spec = ExtractSpec("barcode", method="flank", up=UP, down=DOWN, flank_mismatch=2)
        self.assertEqual(locate(spec, "NN" + up_del + BARCODE + DOWN + "AC"), BARCODE)
        strict = ExtractSpec("barcode", method="flank", up=UP, down=DOWN, flank_mismatch=0)
        self.assertEqual(locate(strict, "NN" + up_del + BARCODE + DOWN + "AC"),
                         "barcode_upstream_not_found")


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

    def test_template_fix(self):
        # fixed bases that differ from the template are corrected; N positions are untouched
        spec = ExtractSpec("barcode", template="NNNNATNNNNATNNNNATNN", template_mismatch=1)
        fix = SequenceExtractor(spec).fix
        self.assertEqual(fix(b"CCCCCTGGGGATCCCCATGG"), b"CCCCATGGGGATCCCCATGG")
        self.assertEqual(fix(b"CCCCATGGGGATCCCCATGG"), b"CCCCATGGGGATCCCCATGG")
        self.assertEqual(fix(b"CCCCNTGGGGATCCCCATGG"), b"CCCCATGGGGATCCCCATGG")
        self.assertEqual(SequenceExtractor(ExtractSpec("barcode")).fix(b"ACGT"), b"ACGT")

    def test_template_fix_like_reference(self):
        # same result as the per-base loop of the request, for random sequences
        def reference(seq, temp):
            return "".join(t if t != "N" and s != t else s for s, t in zip(seq, temp))
        temp = "NNATNNNNATNNNNATNNNNATNNNNATNNNNATNNNN"
        fix = SequenceExtractor(ExtractSpec("barcode", template=temp, template_mismatch=2)).fix
        rng = random.Random(4)
        for _ in range(2000):
            seq = "".join(t if t != "N" and rng.random() > 0.05 else rng.choice("ACGTN")
                          for t in temp)
            self.assertEqual(fix(seq.encode()).decode(), reference(seq, temp))

    def test_matches_template_like_reference(self):
        # check_barcode of sequence_utils.py: N positions are free, the others are counted
        fixed = tuple((i, ord(b)) for i, b in enumerate("NAN") if b != "N")
        self.assertTrue(matches_template(b"CAG", 0, 3, fixed, 3, 0))
        self.assertFalse(matches_template(b"CCG", 0, 3, fixed, 3, 0))
        self.assertTrue(matches_template(b"CCG", 0, 3, fixed, 3, 1))
        exact = template_regex("NAN")                       # regex fast path, same answers
        self.assertTrue(matches_template(b"CAG", 0, 3, fixed, 3, 0, exact))
        self.assertFalse(matches_template(b"CCG", 0, 3, fixed, 3, 0, exact))
        self.assertTrue(matches_template(b"CCG", 0, 3, fixed, 3, 1, exact))


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
