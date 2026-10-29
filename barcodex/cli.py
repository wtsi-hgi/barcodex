"""
Command line: options, validation, the command log and the entry point main().

NAMING CONVENTION
  Every command-line option is a long snake_case option (--max_reads, --barcode_len, ...),
  identical to the variable name used in the code (args.max_reads).  There are no short options.
"""

import argparse
import os
import shlex
import socket
import sys
from datetime import datetime

from . import __version__
from .common import default_mem_gb, default_threads, log
from .constants import MAX_TARGET_LEN
from .extract import FLANK_DISTANCES, METHODS, spec_from_args, validate_spec

# file (in the output directory) that collects the command of every run
CMD_LOG = "cmd.log"


# ---- options ---------------------------------------------------------------

def add_extraction_options(group, role):
    """
    Add the extraction options of one role (--<role>_start, _len, _up, _down, _min_len,
    _allowed_lens).  Barcode and target get the same set.
    Parameters:
            -- group: argparse argument group
            -- role: str: 'barcode' or 'target'
    Returns:
            -- None: nothing
    """
    group.add_argument("--%s_start" % role, type=int, default=0,
                       help="method pos: 0-based start")
    group.add_argument("--%s_len" % role, type=int, default=0,
                       help="method pos: length (0 = to the end of the read)")
    group.add_argument("--%s_up" % role, default=None,
                       help="method flank: upstream flanking sequence")
    group.add_argument("--%s_down" % role, default=None,
                       help="method flank: downstream flanking sequence")
    group.add_argument("--%s_min_len" % role, type=int, default=0,
                       help="reject %ss shorter than this (0 = off)" % role)
    group.add_argument("--%s_allowed_lens" % role, default=None,
                       help="keep only %ss of these lengths, e.g. '19,20,21' or '18-22' "
                            "(any method; recommended with flank)" % role)


def parse_args(argv=None):
    """
    Define and parse the command line.  All options are long snake_case options: the
    command-line name equals the variable name used in the code (--max_reads -> args.max_reads).
    Parameters:
            -- argv: list[str] or None: argument list (None = sys.argv[1:])
    Returns:
            -- argparse.Namespace: parsed options
    """
    # step1 : create the parser
    ap = argparse.ArgumentParser(
        prog="barcodex",
        description="Out-of-core, quality-aware barcode -> target association for huge "
                    "paired-end datasets with one-target-per-barcode resolution.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # step2 : input FASTQ files, output prefix and the extraction methods
    ap.add_argument("--barcode_read", nargs="+",
                    help="FASTQ(.gz) holding the barcode (R1 or R2, one or more files)")
    ap.add_argument("--target_read", nargs="+",
                    help="FASTQ(.gz) holding the target sequence, same order")
    ap.add_argument("--out", help="output prefix")
    ap.add_argument("--method", nargs=2, choices=METHODS, metavar=("BARCODE", "TARGET"),
                    help="REQUIRED: how the barcode and the target are taken from their reads, "
                         "chosen independently: none = the whole read; pos = by position "
                         "(--*_start, --*_len); flank = between flanking sequences (--*_up, "
                         "--*_down). Example: --method flank none")
    ap.add_argument("--flank_mismatch", type=int, default=2,
                    help="method flank: max mismatches in each flanking sequence, "
                         "barcode and target; base qualities are not used to find the flanks")
    ap.add_argument("--flank_distance", choices=FLANK_DISTANCES, default="hamming",
                    help="method flank: how a flank is matched to a window of its length: "
                         "hamming = substitutions only (numpy); levenshtein = substitutions "
                         "and indels (edlib)")

    # step3 : barcode extraction (S1)
    g = ap.add_argument_group("barcode extraction (in --barcode_read)")
    add_extraction_options(g, "barcode")
    g.add_argument("--barcode_temp", default=None,
                   help="barcode template of A/C/G/T/N, e.g. NNNNATNNNNATNN: the barcode must "
                        "have this length and match the fixed (non-N) bases")
    g.add_argument("--barcode_temp_mismatch", type=int, default=1,
                   help="max mismatches at the fixed bases of --barcode_temp")

    # step4 : barcode clustering (S2)
    g = ap.add_argument_group("barcode clustering")
    g.add_argument("--barcode_dist", type=int, default=1,
                   help="max edits between barcodes at CONFIDENT bases")
    g.add_argument("--barcode_max_lowq", type=int, default=1,
                   help="max extra mismatches at LOW-quality barcode bases (cheap). "
                        "0 = ignore barcode quality. Barcodes > 32 bp: quality is not used")
    g.add_argument("--barcode_min_qual", type=int, default=20,
                   help="Phred below this = low quality (barcode read)")
    g.add_argument("--barcode_indel", type=int, default=0,
                   help="max indels in a barcode (indel matching ignores quality)")
    g.add_argument("--barcode_ratio", type=float, default=5.0,
                   help="barcode parent count >= ratio x child count")
    g.add_argument("--min_child_reads", type=int, default=1,
                   help="S2 speed-up: barcodes with fewer reads than this are skipped "
                        "entirely (neither merged into a cluster nor kept). Default 1 = "
                        "cluster everything, singletons included. 2 = skip singletons: "
                        "much faster, but their reads no longer count towards cluster coverage")
    g.add_argument("--min_barcode_reads", type=int, default=None,
                   help="drop barcodes with fewer reads that match no parent "
                        "(default: --min_pair_reads); keeps the barcode table small")

    # step5 : target extraction (S1)
    g = ap.add_argument_group("target extraction (in --target_read)")
    add_extraction_options(g, "target")

    # step6 : target clustering (S4)
    g = ap.add_argument_group("target clustering within a barcode (PCR / sequencing errors)")
    g.add_argument("--target_dist", type=int, default=2,
                   help="max confident edits (indels + high-quality mismatches)")
    g.add_argument("--target_max_lowq", type=int, default=5,
                   help="max extra mismatches at LOW-quality bases (free). 0 = ignore quality")
    g.add_argument("--target_min_qual", type=int, default=20,
                   help="Phred below this = low quality (target read)")
    g.add_argument("--target_max_indel", type=int, default=1,
                   help="max indels in target (0 = off)")
    g.add_argument("--target_ratio", type=float, default=5.0,
                   help="pair coverage rule: parent reads >= ratio x child reads")
    g.add_argument("--target_lowq_ratio", type=float, default=1.0,
                   help="same, when differences are only at low-quality bases")

    # step7 : pair resolution options (S5) and run-mode switches
    g = ap.add_argument_group("pair resolution (S5, DuckDB)")
    g.add_argument("--min_pair_reads", type=int, default=2,
                   help="min reads supporting the (barcode, target) pair")
    g.add_argument("--min_purity", type=float, default=0.8,
                   help="min fraction of a barcode's reads supporting its top target")
    g.add_argument("--merge_dist", type=int, default=0,
                   help="merge targets within this Hamming distance ACROSS barcodes "
                        "(0 = exact only; keep 0 for short targets)")
    g.add_argument("--merge_ratio", type=float, default=5.0,
                   help="across-barcode merge: parent total reads >= ratio x child total reads")
    g.add_argument("--db", default="",
                   help="DuckDB file to keep (default: in-memory, spills to <out>.duckdb.tmp)")
    g.add_argument("--resolve_only", action="store_true",
                   help="skip S1-S4, re-run S5 on the existing <out>.meta/ files")
    g.add_argument("--no_resolve", action="store_true",
                   help="stop after S4 (writes only <out>.meta/target_pairs.tsv)")
    g.add_argument("--clean_meta", action="store_true",
                   help="delete <out>.meta/ after S5 (then --resolve_only is no longer possible)")

    # step8 : resources and housekeeping options
    g = ap.add_argument_group("resources")
    g.add_argument("--threads", type=int, default=default_threads(),
                   help="worker processes (default: job slots from LSF/SLURM, else affinity)")
    g.add_argument("--mem_gb", type=float, default=default_mem_gb(),
                   help="total RAM budget; set to your bsub -M / --mem value")
    g.add_argument("--shards", type=int, default=64, help="number of S3/S4 shards")
    g.add_argument("--chunk", type=int, default=100000, help="reads per parse chunk")
    g.add_argument("--agg_entries", type=int, default=0,
                   help="max aggregated (barcode,target) entries per worker before a part is "
                        "written (0 = derive from --mem_gb)")
    g.add_argument("--compact_parts", type=int, default=12,
                   help="compact text parts into Parquet when this many are pending")
    g.add_argument("--merge_files", type=int, default=6,
                   help="merge this many same-level Parquet files into the next level")
    g.add_argument("--max_cand", type=int, default=5000,
                   help="seeds shared by more than this many heavy barcodes are ignored")
    g.add_argument("--max_dp_group", type=int, default=1000,
                   help="barcodes with more distinct target variants than this skip the "
                        "(slow) indel search")
    g.add_argument("--max_parents", type=int, default=200,
                   help="within one barcode, each target variant is compared only with the N "
                        "most abundant variants (bounds the O(n^2) cost of huge/junk barcodes)")
    g.add_argument("--tmp", default="",
                   help="parent dir for the work dir (default: next to --out); use fast scratch")
    g.add_argument("--keep_tmp", action="store_true")
    g.add_argument("--resume", action="store_true",
                   help="continue a previous run from the last finished step")
    g.add_argument("--no_pigz", action="store_true", help="do not use pigz even if installed")
    g.add_argument("--max_reads", type=int, default=0,
                   help="pilot run: only use the first ~N reads")

    # step9 : general options
    ap.add_argument("--phred_offset", type=int, default=33)
    ap.add_argument("--barcode_map", action="store_true", help="write barcode -> centre table")
    ap.add_argument("--fastq", action="store_true", help="write representative FASTQ pair")
    ap.add_argument("--selftest", action="store_true", help="run on synthetic data and exit")
    ap.add_argument("--version", action="version", version="%(prog)s " + __version__)
    # step10 : parse
    return ap.parse_args(argv)


# ---- validation ------------------------------------------------------------

def validate_args(args):
    """
    Check the inputs and the extraction options before any work starts.
    Parameters:
            -- args: argparse.Namespace: parsed options
    Returns:
            -- None: nothing (raises SystemExit with all problems found)
    """
    # step1 : input files and the required --method
    errors = []
    if not args.barcode_read or not args.target_read:
        errors.append("--barcode_read and --target_read are required")
    elif len(args.barcode_read) != len(args.target_read):
        errors.append("--barcode_read and --target_read need the same number of files")
    if not args.method:
        errors.append("--method is required, e.g. --method flank none "
                      "(barcode method, target method; each one of %s)" % ", ".join(METHODS))
    # step2 : extraction options of both roles
    for role in ("barcode", "target"):
        try:
            max_len = MAX_TARGET_LEN if role == "target" else 0
            errors += validate_spec(spec_from_args(args, role, True, max_len))
        except ValueError as e:
            errors.append("--%s_allowed_lens: %s" % (role, e))
    if errors:
        raise SystemExit("ERROR:\n  " + "\n  ".join(errors))
    # step3 : flank method without a length filter can produce very short barcodes (S2 seeds)
    barcode_spec = spec_from_args(args, "barcode", require_full_len=True)
    if barcode_spec.method == "flank" and not (barcode_spec.allowed_lens or barcode_spec.min_len
                                               or barcode_spec.template):
        log("NOTE: barcode method flank without a length filter: very short barcodes "
            "can make S2 fail ('barcodes too short'); consider --barcode_allowed_lens, "
            "--barcode_min_len or --barcode_temp")


# ---- command log -----------------------------------------------------------

def command_text(argv):
    """
    The command as a copy-and-paste-able shell text: program, then one option per line.
    Parameters:
            -- argv: list[str] or None: arguments given to main() (None = sys.argv[1:])
    Returns:
            -- str: command text (lines joined with ' \\' line continuations)
    """
    # step1 : program: python <script> / python -m barcodex / installed command
    if argv is None:
        script, argv = sys.argv[0], sys.argv[1:]
        if script.endswith("__main__.py"):
            program = "python -m barcodex"
        elif script.endswith(".py"):
            program = "python " + shlex.quote(os.path.abspath(script))
        else:
            program = shlex.quote(script)
    else:
        program = "barcodex"
    # step2 : one line per option and its values
    lines = [program]
    for arg in argv:
        if arg.startswith("--") or len(lines) == 1:
            lines.append(shlex.quote(arg))
        else:
            lines[-1] += " " + shlex.quote(arg)
    return " \\\n    ".join(lines)


def write_cmd_log(args, argv):
    """
    Append the command of this run to <output dir>/cmd.log (date, version, host, working
    directory and the command), so the run can be documented and repeated.
    Parameters:
            -- args: argparse.Namespace: uses args.out
            -- argv: list[str] or None: arguments given to main()
    Returns:
            -- str: path of the log file
    """
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, CMD_LOG)
    with open(path, "a") as fo:
        fo.write("# %s  barcodex %s  host %s\n" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                                                 __version__, socket.gethostname()))
        fo.write("cd %s\n%s\n\n" % (shlex.quote(os.getcwd()), command_text(argv)))
    return path


# ---- entry point -----------------------------------------------------------

def main(argv=None):
    """
    Parse the command line and run the steps (S1..S5), the self-test or S5 only.
    Parameters:
            -- argv: list[str] or None: argument list (None = sys.argv[1:])
    Returns:
            -- None: nothing (exits with a message on invalid input)
    """
    # step1 : parse; --selftest leaves early
    args = parse_args(argv)
    if args.selftest:
        from .selftest import selftest
        return selftest()
    if not args.out:
        raise SystemExit("ERROR: --out is required")
    threads = max(1, args.threads)
    # step2 : --resolve_only: S5 only
    if args.resolve_only:
        log("command written to %s" % write_cmd_log(args, argv))
        from .resolve import resolve_pairs
        resolve_pairs(args, threads)
        return
    # step3 : validate the options, record the command, then run S1 -> S5
    validate_args(args)
    log("command written to %s" % write_cmd_log(args, argv))
    from .steps import run_steps
    run_steps(args, threads)
