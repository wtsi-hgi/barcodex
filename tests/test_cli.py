"""
Command line (barcodex/cli.py): options, validation and the command log.
"""

import contextlib
import io
import os
import unittest

from helpers import TempDirTestCase
from barcodex import __version__
from barcodex.cli import command_text, main, parse_args, validate_args, write_cmd_log
from barcodex.extract import spec_from_args


def args_for(*options):
    """Parse options with dummy input / output names."""
    return parse_args(["--barcode_read", "bc.fq", "--target_read", "tg.fq", "--out", "out"]
                      + list(options))


def error_of(args):
    """Message of the SystemExit raised by validate_args."""
    with unittest.TestCase().assertRaises(SystemExit) as error:
        validate_args(args)
    return str(error.exception)


class TestOptions(unittest.TestCase):

    def test_method_is_required(self):
        self.assertIn("--method is required", error_of(args_for()))

    def test_methods_are_independent(self):
        args = args_for("--method", "pos", "flank", "--barcode_len", "16",
                        "--target_up", "ACGTACGT", "--target_down", "TTGGCCAA")
        validate_args(args)
        self.assertEqual(spec_from_args(args, "barcode", True).method, "pos")
        target = spec_from_args(args, "target", False)
        self.assertEqual((target.method, target.flank_mismatch), ("flank", 2))

    def test_flank_mismatch_is_shared(self):
        args = args_for("--method", "flank", "flank", "--flank_mismatch", "1",
                        "--barcode_up", "AAAA", "--barcode_down", "CCCC",
                        "--target_up", "GGGG", "--target_down", "TTTT")
        self.assertEqual([spec_from_args(args, role, True).flank_mismatch
                          for role in ("barcode", "target")], [1, 1])

    def test_flank_distance(self):
        args = args_for("--method", "flank", "flank", "--flank_distance", "levenshtein",
                        "--barcode_up", "AAAA", "--barcode_down", "CCCC",
                        "--target_up", "GGGG", "--target_down", "TTTT")
        validate_args(args)
        self.assertEqual([spec_from_args(args, role, True).flank_distance
                          for role in ("barcode", "target")], ["levenshtein", "levenshtein"])
        self.assertEqual(args_for("--method", "none", "none").flank_distance, "hamming")
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            args_for("--method", "flank", "none", "--flank_distance", "edit")

    def test_unknown_method(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            args_for("--method", "whole", "none")

    def test_option_of_another_method(self):
        message = error_of(args_for("--method", "none", "none", "--barcode_len", "16"))
        self.assertIn("need the barcode method pos", message)

    def test_flank_without_flanks(self):
        message = error_of(args_for("--method", "flank", "none", "--barcode_up", "ACGT"))
        self.assertIn("barcode method flank needs --barcode_up and --barcode_down", message)

    def test_invalid_template_and_min_len(self):
        message = error_of(args_for("--method", "none", "none", "--barcode_temp", "NNXX",
                                    "--target_min_len", "-5"))
        self.assertIn("--barcode_temp may only contain", message)
        self.assertIn("--target_min_len must be >= 0", message)

    def test_target_min_len_above_mask_limit(self):
        message = error_of(args_for("--method", "none", "none", "--target_min_len", "300"))
        self.assertIn("255 bp limit", message)

    def test_no_short_options(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            parse_args(["-b", "bc.fq"])

    def test_version(self):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), self.assertRaises(SystemExit):
            parse_args(["--version"])
        self.assertEqual(stdout.getvalue().split(), ["barcodex", __version__])
        self.assertEqual(__version__, "0.1")


class TestCommandLog(TempDirTestCase):

    def test_command_text_one_option_per_line(self):
        text = command_text(["--barcode_read", "a.fq", "b.fq", "--method", "flank", "none"])
        self.assertEqual(text.splitlines(), ["barcodex \\", "    --barcode_read a.fq b.fq \\",
                                             "    --method flank none"])

    def test_cmd_log_is_appended_in_output_dir(self):
        args = args_for("--out", os.path.join(self.tmp, "sub", "prefix"))
        for _ in range(2):
            path = write_cmd_log(args, ["--out", args.out, "--method", "none", "none"])
        self.assertEqual(path, os.path.join(self.tmp, "sub", "cmd.log"))
        with open(path) as fh:
            text = fh.read()
        self.assertEqual(text.count("# "), 2)
        self.assertIn("barcodex %s" % __version__, text)
        self.assertIn("--method none none", text)

    def test_main_writes_cmd_log_before_running(self):
        out = os.path.join(self.tmp, "run")
        with self.assertRaises((SystemExit, OSError)):  # stops at the missing FASTQ files
            main(["--barcode_read", "missing.fq", "--target_read", "missing.fq", "--out", out,
                  "--method", "none", "none", "--threads", "1"])
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "cmd.log")))


if __name__ == "__main__":
    unittest.main()
