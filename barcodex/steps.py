"""
Steps of a run: S1 ingest -> S2 barcodes -> S3 shards -> S4 targets -> S5 resolve.
With --resume every step whose marker file exists in the work directory is skipped.
"""

import os
import shutil

from .barcodes import stage_barcodes
from .common import log, log_sub, meta_path
from .ingest import stage_ingest
from .resolve import resolve_pairs
from .shards import stage_shards
from .targets import stage_targets


def prepare_work_dir(args):
    """
    Work directory  <tmp or output dir>/<output name>.work.  A fresh start deletes old work /
    meta directories; --resume keeps them.
    Parameters:
            -- args: argparse.Namespace: out, tmp, resume
    Returns:
            -- str: path of the work directory (created)
    """
    base_dir = args.tmp if args.tmp else os.path.dirname(os.path.abspath(args.out))
    work = os.path.join(base_dir, os.path.basename(args.out) + ".work")
    if args.resume and os.path.isdir(work):
        log("Resuming from existing temporary directory %s." % work)
    else:
        shutil.rmtree(work, ignore_errors=True)
        shutil.rmtree(meta_path(args), ignore_errors=True)
    os.makedirs(work, exist_ok=True)
    return work


def run_steps(args, threads):
    """
    Run S1..S5 (skipping the steps already finished when --resume is given).
    Parameters:
            -- args: argparse.Namespace: validated options
            -- threads: int: worker processes / DuckDB threads
    Returns:
            -- None: nothing; results are files
    """
    # step1 : resources and work directory
    mem_gb = args.mem_gb
    log("resources: %d process(es), %.1f GB RAM budget" % (threads, mem_gb))
    work = prepare_work_dir(args)
    n_shards = max(1, args.shards)

    # step2 : run S1 -> S4; with --resume every step whose marker file exists is skipped
    def stage_done(marker):
        """
        True when --resume is on and the marker file of a step exists in the work directory.
        Parameters:
                -- marker: str: marker file name, e.g. 'ingest.done'
        Returns:
                -- bool: True = the step can be skipped
        """
        return args.resume and os.path.exists(os.path.join(work, marker))

    if stage_done("shards.done"):
        log_sub("resume: shards already built")
    else:
        if stage_done("bcmap.done"):
            log_sub("resume: barcode map already built")
        else:
            if stage_done("ingest.done"):
                log_sub("resume: ingest already done")
            else:
                #      |----   substep S1 ingest   (FASTQ -> tuples/t*.parquet, marker ingest.done)
                shutil.rmtree(work + "/parts", ignore_errors=True)
                shutil.rmtree(work + "/tuples", ignore_errors=True)
                stage_ingest(args, work, threads, mem_gb)
            #      |----   substep S2 barcodes (tuples -> bcmap.parquet + meta files, marker bcmap.done)
            stage_barcodes(args, work, threads, mem_gb)
        #      |----   substep S3 shards   (tuples x bcmap -> shards/, marker shards.done)
        stage_shards(args, work, threads, mem_gb, n_shards)
    #      |----   substep S4 targets  (shards -> <out>.meta/target_pairs.tsv)
    stage_targets(args, work, threads, mem_gb, n_shards)
    # step3 : delete the work directory unless --keep_tmp
    if not args.keep_tmp:
        shutil.rmtree(work, ignore_errors=True)
    # step4 : S5: one target per barcode and the result files (skipped with --no_resolve)
    if not args.no_resolve:
        resolve_pairs(args, threads)
    log("Done")
