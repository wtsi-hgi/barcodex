"""
S3 shards: tuples x bcmap -> Parquet shards keyed by centre barcode.
"""

import os
import shutil

from .common import log, log_sub, open_duckdb, sql_str

# ---- stage driver -----------------------------------------------------

def stage_shards(args, work, threads, mem_gb, n_shards):
    """
    S3: replace every barcode by its centre, merge identical (centre, target) rows (counts add,
    target masks AND) and split the result into n_shards Hive-partitioned Parquet shards by
    hash(centre), so that all rows of one barcode cluster are in the same shard.
    Deletes the tuple store and bcmap.parquet afterwards (unless --keep_tmp) and writes
    <work>/shards.done.
    Parameters:
            -- args: argparse.Namespace: keep_tmp
            -- work: str: work directory
            -- threads: int: DuckDB threads
            -- mem_gb: float: RAM budget in GB (85% is given to DuckDB)
            -- n_shards: int: number of shards (--shards)
    Returns:
            -- None: nothing; results are files in <work>/shards/shard=<n>/
    """
    log("S3: re-keying tuples by corrected barcode into %d shard(s) ..." % n_shards)
    # step1 : remove old shards and open DuckDB
    shards_dir = work + "/shards"
    shutil.rmtree(shards_dir, ignore_errors=True)
    con = open_duckdb(work, threads, mem_gb * 0.85, "shard")
    # step2 : tuples JOIN barcode map -> (centre, target) rows; counts add, masks AND; write shards by hash(centre)
    n_rows = con.execute("""
        COPY (
            SELECT m.centre AS centre, t.target AS target, SUM(t.cnt)::UBIGINT AS cnt,
                   bit_and(t.tm0) AS tm0, bit_and(t.tm1) AS tm1,
                   bit_and(t.tm2) AS tm2, bit_and(t.tm3) AS tm3,
                   (hash(m.centre) %% %d)::INTEGER AS shard
            FROM read_parquet(%s) t JOIN read_parquet(%s) m ON t.barcode = m.barcode
            GROUP BY m.centre, t.target
        ) TO %s (FORMAT PARQUET, PARTITION_BY (shard), COMPRESSION ZSTD)
    """ % (n_shards, sql_str(work + "/tuples/t*.parquet"), sql_str(work + "/bcmap.parquet"),
           sql_str(shards_dir))).fetchone()[0]
    con.close()
    # step3 : report the size of the shards
    shard_gb = sum(os.path.getsize(os.path.join(d, f))
                   for d, _, files in os.walk(shards_dir) for f in files) / 1e9
    log_sub("%d distinct (barcode cluster, target sequence) rows; shards %.2f GB" % (n_rows, shard_gb))
    # step4 : free disk space: delete the tuple store and bcmap.parquet
    if not args.keep_tmp:
        shutil.rmtree(work + "/tuples", ignore_errors=True)
        if os.path.exists(work + "/bcmap.parquet"):
            os.remove(work + "/bcmap.parquet")
    # step5 : marker file used by --resume
    with open(work + "/shards.done", "w") as fo:
        fo.write("%d\n" % n_shards)
