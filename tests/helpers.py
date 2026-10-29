"""
Shared helpers of the unit tests: package import path, default options, temporary directories
and small input files.
"""

import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from barcodex.cli import parse_args  # noqa: E402
from barcodex.common import open_duckdb, sql_str  # noqa: E402

# column layout of the S1 tuple store (see compaction.PART_COLUMNS)
TUPLE_COLUMNS = ("barcode", "target", "cnt", "tm0", "tm1", "tm2", "tm3", "bm")


def make_args(out, **overrides):
    """
    Default command-line options with an output prefix, plus overrides.
    Parameters:
            -- out: str: output prefix
            -- overrides: option name -> value
    Returns:
            -- argparse.Namespace: options as main() would see them
    """
    args = parse_args(["--out", out, "--threads", "1", "--mem_gb", "1"])
    for name, value in overrides.items():
        if not hasattr(args, name):
            raise AttributeError("unknown option %s" % name)
        setattr(args, name, value)
    return args


def write_fastq(path, records):
    """
    Write a FASTQ file.
    Parameters:
            -- path: str: output file
            -- records: list[tuple]: [(sequence, quality)]; quality None = all 'I'
    Returns:
            -- str: path
    """
    with open(path, "w") as fo:
        for i, (seq, qual) in enumerate(records):
            fo.write("@r%d\n%s\n+\n%s\n" % (i, seq, qual or "I" * len(seq)))
    return path


def write_tuples(path, rows):
    """
    Write a tuple-store Parquet file (as produced by S1).
    Parameters:
            -- path: str: output file
            -- rows: list[tuple]: (barcode, target, cnt) or (barcode, target, cnt, tm0, tm1, tm2,
               tm3, bm); missing masks are 0
    Returns:
            -- str: path
    """
    work = os.path.dirname(path)
    con = open_duckdb(work, 1, 1, "test")
    values = ", ".join(
        "(%s, %s, %d, %d, %d, %d, %d, %d)" % ((sql_str(r[0]), sql_str(r[1])) + tuple(
            (list(r[2:]) + [0] * 6)[:6])) for r in rows)
    con.execute("COPY (SELECT col0::VARCHAR AS barcode, col1::VARCHAR AS target, "
                "col2::UBIGINT AS cnt, col3::UBIGINT AS tm0, col4::UBIGINT AS tm1, "
                "col5::UBIGINT AS tm2, col6::UBIGINT AS tm3, col7::UBIGINT AS bm "
                "FROM (VALUES %s)) TO %s (FORMAT PARQUET)" % (values, sql_str(path)))
    con.close()
    return path


def read_parquet(path, query="SELECT * FROM t ORDER BY ALL"):
    """
    Read a Parquet file (or glob) with DuckDB.
    Parameters:
            -- path: str: file or glob
            -- query: str: SQL on the view t
    Returns:
            -- list[tuple]: rows
    """
    import duckdb
    con = duckdb.connect()
    con.execute("CREATE VIEW t AS SELECT * FROM read_parquet(%s)" % sql_str(path))
    rows = con.execute(query).fetchall()
    con.close()
    return rows


def read_tsv(path):
    """
    Read a tab-separated file with a header line.
    Parameters:
            -- path: str: file
    Returns:
            -- list[dict]: one dict per row
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


class TempDirTestCase(unittest.TestCase):
    """
    TestCase with a fresh temporary directory self.tmp (removed afterwards).
    """

    def setUp(self):
        """Create the temporary directory."""
        self.tmp = tempfile.mkdtemp(prefix="barcodex_test_")

    def tearDown(self):
        """Remove the temporary directory."""
        shutil.rmtree(self.tmp, ignore_errors=True)
