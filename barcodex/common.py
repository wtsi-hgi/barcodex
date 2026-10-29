"""
Common helpers used by every stage: logging, SQL quoting, resources (threads, RAM, DuckDB)
and paths.
"""

import glob
import os
from datetime import datetime

# ---- logging and SQL quoting ------------------------------------------

# printed in front of every sub-step message
SUB_PREFIX = "    |---- "      # printed in front of every sub-step message

def log(message, sub=False):
    """
    Print one progress line:  <date time> [|---- ]<message>.
    Main steps use log(...), their sub-steps use log_sub(...).
    Parameters:
            -- message: str: text to print
            -- sub: bool: True = indent the line with SUB_PREFIX (sub-step line); default False
    Returns:
            -- None: prints to stdout and flushes immediately
    """
    prefix = SUB_PREFIX if sub else ""
    print(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {prefix}{message}", flush=True)


def log_sub(message):
    """
    Print a sub-step line (indented with SUB_PREFIX).  Shortcut for log(message, sub=True).
    Parameters:
            -- message: str: text to print
    Returns:
            -- None: prints to stdout
    """
    log(message, sub=True)


def sql_str(text):
    """
    Quote a Python string as an SQL string literal (single quotes are doubled).
    Parameters:
            -- text: str: raw text, e.g. a file path
    Returns:
            -- str: the text wrapped in single quotes and safe to paste into SQL
    """
    return "'" + text.replace("'", "''") + "'"


def sql_list(files):
    """
    Quote a list of file names as an SQL list literal, e.g. ['a.parquet','b.parquet'].
    Parameters:
            -- files: list[str]: file names
    Returns:
            -- str: SQL list literal
    """
    return "[" + ",".join(sql_str(f) for f in files) + "]"


# ---- resources (threads, RAM, DuckDB) and paths -----------------------

def default_threads():
    """
    Number of CPU slots of THIS job.  os.cpu_count() returns the cores of the whole node on a
    cluster (e.g. 128), so batch-scheduler variables are checked first.
    Parameters:
            -- None: no parameters
    Returns:
            -- int: slots from LSB_DJOB_NUMPROC / SLURM_CPUS_PER_TASK / NSLOTS / PBS_NP, else the CPU
               affinity of the process, else os.cpu_count(), else 1
    """
    # step1 : batch-scheduler variables first (LSF / SLURM / SGE / PBS)
    for var in ("LSB_DJOB_NUMPROC", "SLURM_CPUS_PER_TASK", "NSLOTS", "PBS_NP"):
        value = os.environ.get(var, "")
        if value.isdigit() and int(value) > 0:
            return int(value)
    # step2 : no scheduler variable: CPU affinity of this process, else all cores
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def default_mem_gb():
    """
    Default RAM budget: 85% of the container (cgroup v2 or v1) memory limit, or 32 GB if there is none.
    Parameters:
            -- None: no parameters
    Returns:
            -- float: memory budget in GB
    """
    # step1 : container (cgroup v2, then v1) memory limit -> 85% of it
    #      |----   substep a huge value (>= 2^50) means 'no limit' and is ignored
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            with open(path) as fh:
                text = fh.read().strip()
            if text.isdigit() and int(text) < (1 << 50):
                return round(int(text) / 1e9 * 0.85, 1)
        except OSError:
            pass
    # step2 : no limit found: assume 32 GB
    return 32.0


def open_duckdb(work, threads, mem_gb, name, path=":memory:"):
    """
    Open a DuckDB connection with thread / memory limits and its own spill directory.
    preserve_insertion_order is switched off so DuckDB can stream and spill freely.
    Parameters:
            -- work: str: work directory; the spill directory is <work>/duck_<name>
            -- threads: int: DuckDB threads (at least 1)
            -- mem_gb: float: DuckDB memory_limit in GB (at least 300 MB)
            -- name: str: label of the connection, used in the spill directory name
            -- path: str: database file, or ':memory:' (default)
    Returns:
            -- duckdb.DuckDBPyConnection: ready-to-use connection (SystemExit with a message if duckdb
               is not installed)
    """
    # step1 : import duckdb (clear error message when it is missing)
    try:
        import duckdb
    except ImportError:
        raise SystemExit("ERROR: barcodex needs duckdb  (pip install duckdb)")
    # step2 : connect and apply the resource limits
    #      |----   substep threads, memory_limit and a private spill directory
    #      |----   substep preserve_insertion_order=false lets DuckDB stream / spill freely
    con = duckdb.connect(path)
    con.execute("SET threads=%d" % max(1, threads))
    con.execute("SET memory_limit=%s" % sql_str("%dMB" % max(300, int(mem_gb * 1000))))
    con.execute("SET temp_directory=%s" % sql_str(os.path.join(work, "duck_" + name)))
    con.execute("SET preserve_insertion_order=false")      # lets DuckDB stream/spill freely
    return con


def meta_path(args, name=""):
    """
    Path of the intermediate files that S5 needs (kept next to the outputs so --resolve_only works):
      <out>.meta/barcode_clusters.parquet   per barcode-cluster statistics   (written by S2)
      <out>.meta/raw_stats.json             raw barcode / target counts      (written by S2)
      <out>.meta/target_pairs.tsv           (barcode, target cluster) rows   (written by S4)
    Parameters:
            -- args: argparse.Namespace: uses args.out (output prefix)
            -- name: str: file name inside the meta directory; '' (default) returns the directory itself
    Returns:
            -- str: path of the meta directory or of the file inside it
    """
    return os.path.join(args.out + ".meta", name) if name else args.out + ".meta"


def dir_gb(path, pattern):
    """
    Total size of the files in a directory that match a glob pattern (used for log messages).
    Parameters:
            -- path: str: directory
            -- pattern: str: glob pattern, e.g. 't*.parquet'
    Returns:
            -- float: size in GB
    """
    return sum(os.path.getsize(f) for f in glob.glob(os.path.join(path, pattern))) / 1e9
