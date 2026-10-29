"""
S1 compaction: gz text parts written by the parse workers -> tiered, de-duplicated Parquet
files  <work>/tuples/t*.parquet.
"""

import glob
import os
import threading
import traceback

from .common import dir_gb, log_sub, open_duckdb, sql_list, sql_str

# Column layout of the text parts (must match flush_part).
PART_COLUMNS = ("{'barcode':'VARCHAR','target':'VARCHAR','cnt':'UBIGINT',"
                "'tm0':'UBIGINT','tm1':'UBIGINT','tm2':'UBIGINT','tm3':'UBIGINT',"
                "'bm':'UBIGINT'}")

def compaction_sql(part_files, parquet_files, out_file):
    """
    SQL that merges text parts and/or Parquet files into one Parquet file.
    Rows with the same (barcode, target) are merged: counts add up, masks are AND-ed.
    Parameters:
            -- part_files: list[str]: gz text parts (may be empty)
            -- parquet_files: list[str]: Parquet files to merge (may be empty)
            -- out_file: str: Parquet file to write
    Returns:
            -- str: COPY ... TO statement
    """
    # step1 : one sub-query per input kind (text parts and / or Parquet files)
    sources = []
    if part_files:
        sources.append("SELECT * FROM read_csv(%s, delim='\t', header=false, columns=%s, "
                       "compression='gzip')" % (sql_list(part_files), PART_COLUMNS))
    if parquet_files:
        sources.append("SELECT * FROM read_parquet(%s)" % sql_list(parquet_files))
    # step2 : merge rows with the same (barcode, target): counts add up, masks are AND-ed; write ZSTD Parquet
    return ("COPY (SELECT barcode, target, SUM(cnt)::UBIGINT AS cnt, "
            "bit_and(tm0) AS tm0, bit_and(tm1) AS tm1, bit_and(tm2) AS tm2, "
            "bit_and(tm3) AS tm3, bit_and(bm) AS bm FROM (%s) "
            "GROUP BY barcode, target) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"
            % (" UNION ALL ".join(sources), sql_str(out_file)))


class Compactor:
    """
    Background tiered compaction:  text parts -> t1_*.parquet -> t2_* -> t3_* -> t4_* ...
    The parent process calls step() regularly; each call starts at most one background DuckDB job.
    Merging early keeps the scratch disk proportional to the number of DISTINCT (barcode, target)
    pairs instead of the number of reads.  finish() runs the remaining jobs after all workers ended.
    Parameters:
            -- work: str: work directory (parts/ and tuples/ are below it)
            -- threads: int: DuckDB threads of one job
            -- mem_gb: float: DuckDB memory of one job
            -- batch: int: text parts merged per level-1 job (--compact_parts)
            -- fanout: int: same-level Parquet files merged per higher-level job (--merge_files)
    Returns:
            -- Compactor: object driven by stage_ingest through step() and finish()
    """

    def __init__(self, work, threads, mem_gb, batch, fanout):
        """
        Store the settings and the directories; no job is started yet.
        Parameters:
                -- work: str: work directory
                -- threads: int: DuckDB threads of one job
                -- mem_gb: float: DuckDB memory of one job
                -- batch: int: text parts per level-1 job
                -- fanout: int: same-level Parquet files per higher-level job
        Returns:
                -- None: constructor
        """
        self.work, self.threads, self.mem_gb = work, threads, mem_gb
        self.parts_dir = work + "/parts"
        self.tuples_dir = work + "/tuples"
        self.batch = batch              # text parts per level-1 job
        self.fanout = fanout            # same-level Parquet files per higher-level job
        self.thread = None
        self.error = None
        self.sequence = 0

    def _pick_job(self, final):
        """
        Choose the next compaction job.  Level 1 takes all pending text parts when at least `batch`
        are waiting (or any at all when final=True).  Levels 2..4 merge `fanout` Parquet files of the
        level below.
        Parameters:
                -- final: bool: True when all workers have finished (compact whatever is left)
        Returns:
                -- tuple or None: (output_level, text_part_files, parquet_files), or None if nothing to
                   do
        """
        # step1 : level 1: compact all pending text parts once `batch` are waiting (or at the end)
        parts = sorted(glob.glob(self.parts_dir + "/*.tsv.gz"))
        if len(parts) >= self.batch or (final and parts):
            return 1, parts, []
        # step2 : levels 2..4: merge `fanout` same-level Parquet files into the next level
        for level in (1, 2, 3):
            files = sorted(glob.glob("%s/t%d_*.parquet" % (self.tuples_dir, level)))
            if len(files) >= self.fanout:
                return level + 1, [], files
        return None

    def _compact(self, part_files, parquet_files, tmp_out):
        """
        Run one compaction in a fresh DuckDB connection (closed afterwards).
        Parameters:
                -- part_files: list[str]: gz text parts to merge
                -- parquet_files: list[str]: Parquet files to merge
                -- tmp_out: str: temporary output file (renamed by the caller)
        Returns:
                -- int: number of rows written
        """
        con = open_duckdb(self.work, self.threads, self.mem_gb, "compact")
        try:
            return con.execute(compaction_sql(part_files, parquet_files, tmp_out)).fetchone()[0]
        finally:
            con.close()

    def _run(self, job):
        """
        Thread body: compact, publish the result atomically, delete the inputs.  Any exception is
        stored in self.error and re-raised by _check() in the supervisor.
        Parameters:
                -- job: tuple: (level, part_files, parquet_files) from _pick_job
        Returns:
                -- None: nothing
        """
        try:
            # step1 : name the output file t<level>_<sequence>.parquet
            level, part_files, parquet_files = job
            self.sequence += 1
            out = "%s/t%d_%06d.parquet" % (self.tuples_dir, level, self.sequence)
            # step2 : compact into <out>.tmp
            n_rows = self._compact(part_files, parquet_files, out + ".tmp")
            # step3 : publish atomically, then delete the inputs
            os.replace(out + ".tmp", out)
            for f in part_files + parquet_files:
                os.remove(f)
            # step4 : report
            log_sub("compacted %d part(s) + %d file(s) -> L%d: %d rows (%.2f GB on disk in tuples/)"
                % (len(part_files), len(parquet_files), level, n_rows,
                   dir_gb(self.tuples_dir, "t*.parquet")))
        except BaseException:
            self.error = traceback.format_exc()

    def _check(self):
        """
        Raise SystemExit if the background job failed.
        Parameters:
                -- None: no parameters
        Returns:
                -- None: nothing (raises SystemExit with the traceback on error)
        """
        if self.error:
            raise SystemExit("ERROR in compaction:\n" + self.error)

    def step(self):
        """
        Start a new background job if the previous one is finished and work is pending.
        Called every few seconds by the supervisor loop of stage_ingest.
        Parameters:
                -- None: no parameters
        Returns:
                -- None: nothing
        """
        self._check()
        if self.thread is not None and self.thread.is_alive():
            return
        job = self._pick_job(False)
        if job:
            self.thread = threading.Thread(target=self._run, args=(job,), daemon=True)
            self.thread.start()

    def finish(self):
        """
        Run jobs until everything is merged (called after all workers are done).
        Parameters:
                -- None: no parameters
        Returns:
                -- None: nothing
        """
        # step1 : loop until no job is left
        while True:
            #      |----   substep wait for the running job and re-raise its error
            if self.thread is not None:
                self.thread.join()
                self.thread = None
            self._check()
            #      |----   substep final=True: also compact the last text parts
            job = self._pick_job(True)
            if job is None:
                return
            self.thread = threading.Thread(target=self._run, args=(job,), daemon=True)
            self.thread.start()
