"""
Self-test: run all steps on small synthetic data sets and check the result.
"""

import json
import random
import tempfile

from .common import log

# flanking sequences used by the 'flank' scenario
SELFTEST_UP = "GTCCTCTTCCCATCT"
SELFTEST_DOWN = "TGGTGGAAGCTACTG"

# read-funnel rows of the stats file that must add up to total_read_pairs
FUNNEL_AFTER_EXTRACTION = ("reads_dropped_barcode_clustering", "reads_low_coverage_barcodes",
                           "reads_ambiguous_barcodes", "reads_minor_targets", "effective_reads")


def write_synthetic_fastq(tmp_dir, scenario, rnd, n_barcodes=300, barcode_len=16, target_len=60,
                          reads_per_barcode=40):
    """
    Write a synthetic FASTQ pair.  Barcodes carry ~1% errors (the changed base gets a LOW
    quality), targets ~0.4% errors (also low quality).
      pos      : barcode read = barcode + 30 random bases;            target read = target
      flank    : barcode read = 5 random + UP + barcode + DOWN + 10;  target read = target + 20
                 random bases; 2% of the barcode reads lose their upstream flank
    Parameters:
            -- tmp_dir: str: output directory
            -- scenario: str: 'pos' or 'flank'
            -- rnd: random.Random: random generator
            -- n_barcodes / barcode_len / target_len / reads_per_barcode: int: data set size
    Returns:
            -- tuple: (truth dict barcode -> target, barcode FASTQ path, target FASTQ path,
               number of reads without upstream flank)
    """
    # step1 : helpers
    def random_seq(n):
        """
        Random DNA sequence.
        Parameters:
                -- n: int: length
        Returns:
                -- str: sequence of A/C/G/T
        """
        return "".join(rnd.choice("ACGT") for _ in range(n))

    def mutate(seq, rate):
        """
        Introduce substitution errors; every changed base gets a LOW quality ('#').
        Parameters:
                -- seq: str: original sequence
                -- rate: float: probability that a base is replaced by a random base
        Returns:
                -- tuple: (mutated sequence, quality string with '#' at mutated and 'I' at other
                   positions)
        """
        out, qual = [], []
        for base in seq:
            if rnd.random() < rate:
                out.append(rnd.choice("ACGT"))
                qual.append("#")
            else:
                out.append(base)
                qual.append("I")
        return "".join(out), "".join(qual)

    # step2 : draw the true barcode -> target table
    truth = {}
    while len(truth) < n_barcodes:
        truth[random_seq(barcode_len)] = random_seq(target_len)
    # step3 : write the FASTQ pair
    barcode_fq, target_fq = tmp_dir + "/bc.fq", tmp_dir + "/tg.fq"
    n_no_flank = 0
    with open(barcode_fq, "w") as f_barcode, open(target_fq, "w") as f_target:
        read_no = 0
        for barcode, target in truth.items():
            for _ in range(reads_per_barcode):
                read_no += 1
                bc_seq, bc_qual = mutate(barcode, 0.01)
                tg_seq, tg_qual = mutate(target, 0.004)
                if scenario == "flank":
                    up = SELFTEST_UP
                    if rnd.random() < 0.02:
                        up = random_seq(len(SELFTEST_UP))
                        n_no_flank += 1
                    bc_seq = random_seq(5) + up + bc_seq + SELFTEST_DOWN + random_seq(10)
                    bc_qual = "I" * (5 + len(up)) + bc_qual + "I" * (len(SELFTEST_DOWN) + 10)
                    tg_seq, tg_qual = tg_seq + random_seq(20), tg_qual + "I" * 20
                else:
                    bc_seq, bc_qual = bc_seq + random_seq(30), bc_qual + "I" * 30
                f_barcode.write("@r%d\n%s\n+\n%s\n" % (read_no, bc_seq, bc_qual))
                f_target.write("@r%d\n%s\n+\n%s\n" % (read_no, tg_seq, tg_qual))
    return truth, barcode_fq, target_fq, n_no_flank


def read_stats(path):
    """
    Read a stats file into a dict.
    Parameters:
            -- path: str: <out>.stats.tsv
    Returns:
            -- dict: stat -> value (str)
    """
    with open(path) as fh:
        next(fh)
        return dict(line.rstrip("\n").split("\t")[1:3] for line in fh)


def run_scenario(scenario, seed=7):
    """
    Run one scenario.  PASSES when >= 97% of the true barcodes map to their true target, none
    maps to a wrong target, <= 2% of the output barcodes are spurious and the read funnel of the
    stats file adds up.
    Parameters:
            -- scenario: str: 'pos' or 'flank'
            -- seed: int: random seed
    Returns:
            -- bool: True = PASS
    """
    from .cli import main
    # step1 : synthetic data
    rnd = random.Random(seed)
    tmp_dir = tempfile.mkdtemp(prefix="barcodex_selftest_%s_" % scenario)
    barcode_len, target_len = 16, 60
    truth, barcode_fq, target_fq, n_no_flank = write_synthetic_fastq(
        tmp_dir, scenario, rnd, barcode_len=barcode_len, target_len=target_len)
    # step2 : run all steps (small chunks / parts so every stage is exercised)
    out = tmp_dir + "/st"
    if scenario == "flank":
        extraction = ["--method", "flank", "pos", "--barcode_up", SELFTEST_UP,
                      "--barcode_down", SELFTEST_DOWN, "--barcode_allowed_lens", "15-17",
                      "--target_len", str(target_len)]
    else:
        extraction = ["--method", "pos", "none", "--barcode_len", str(barcode_len)]
    main(["--barcode_read", barcode_fq, "--target_read", target_fq, "--out", out] + extraction
         + ["--threads", "2", "--mem_gb", "2", "--chunk", "1500", "--agg_entries", "2000",
            "--compact_parts", "2", "--merge_files", "2", "--shards", "4"])
    # step3 : compare the output with the truth
    got = {}
    with open(out + ".barcode_clusters.tsv") as fh:
        next(fh)
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            got[fields[0]] = fields[1]
    n_ok = sum(1 for bc, tg in truth.items() if got.get(bc) == tg)
    n_wrong = sum(1 for bc, tg in got.items() if bc in truth and truth[bc] != tg)
    n_extra = sum(1 for bc in got if bc not in truth)
    log("SELFTEST %s: %d/%d true barcodes mapped correctly, %d wrong, %d spurious barcodes"
        % (scenario, n_ok, len(truth), n_wrong, n_extra))
    # step4 : the read funnel must add up
    stats = read_stats(out + ".stats.tsv")
    with open(out + ".meta/ingest_stats.json") as fh:
        rejected = sum(row[1] for row in json.load(fh)["rejects"])
    funnel_ok = (int(stats["total_read_pairs"]) == rejected + int(stats["extracted_read_pairs"])
                 and int(stats["extracted_read_pairs"])
                 == sum(int(stats[k]) for k in FUNNEL_AFTER_EXTRACTION))
    if scenario == "flank":
        # a random replacement flank can still match within 2 mismatches, hence <=
        funnel_ok = funnel_ok and int(stats["barcode_upstream_not_found"]) <= n_no_flank
    log("SELFTEST %s: read funnel %s" % (scenario, "adds up" if funnel_ok else "does NOT add up"))
    # step5 : PASS rule
    passed = (n_ok >= 0.97 * len(truth) and n_wrong == 0 and n_extra <= 0.02 * len(truth)
              and funnel_ok)
    log("SELFTEST %s %s   (files in %s)" % (scenario, "PASS" if passed else "FAIL", tmp_dir))
    return passed


def selftest():
    """
    Run all scenarios: 'pos' (barcode by position, target = whole read) and 'flank'
    (barcode between flanks, target by position).
    Parameters:
            -- None: no parameters
    Returns:
            -- None: nothing (raises SystemExit(1) when a scenario fails)
    """
    results = [run_scenario(scenario) for scenario in ("pos", "flank")]
    if not all(results):
        raise SystemExit(1)
