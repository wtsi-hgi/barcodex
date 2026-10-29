"""
S1 FASTQ streaming: reader processes cut the FASTQ files into chunks and deal them out to the
parse workers.
"""

import gzip
import queue
import shutil
import subprocess
import traceback
from itertools import islice

def open_stream(path, use_pigz):
    """
    Open a FASTQ file (plain or .gz) for binary reading.  .gz files are decompressed by an external
    pigz process when allowed and installed (much faster than Python gzip), else by gzip.open.
    Parameters:
            -- path: str: FASTQ file name
            -- use_pigz: bool: allow pigz
    Returns:
            -- tuple: (binary file object, subprocess.Popen or None); pass both to close_stream()
    """
    if path.endswith(".gz"):
        if use_pigz and shutil.which("pigz"):
            proc = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE,
                                    bufsize=1 << 22)
            return proc.stdout, proc
        return gzip.open(path, "rb"), None
    return open(path, "rb"), None


def close_stream(stream, proc):
    """
    Close a stream opened by open_stream and stop its pigz process (if any).
    Parameters:
            -- stream: file object returned by open_stream
            -- proc: subprocess.Popen or None: second value returned by open_stream
    Returns:
            -- None: nothing
    """
    try:
        stream.close()
    finally:
        if proc is not None:
            proc.terminate()
            proc.wait()


def peek_read_len(path, use_pigz):
    """
    Length of the first read of a FASTQ file (used to derive the barcode length when --barcode_len is 0).
    Parameters:
            -- path: str: FASTQ file name
            -- use_pigz: bool: allow pigz
    Returns:
            -- int: length of the first sequence line (0 for an empty file)
    """
    stream, proc = open_stream(path, use_pigz)
    try:
        stream.readline()
        return len(stream.readline().strip())
    finally:
        close_stream(stream, proc)


def put_abortable(target_queue, item, abort):
    """
    queue.put that gives up when the abort event is set (avoids deadlocks when a process failed).
    Parameters:
            -- target_queue: multiprocessing queue: destination
            -- item: object: item to send
            -- abort: multiprocessing.Event: set by the supervisor on any error
    Returns:
            -- bool: True = item queued, False = aborted before it could be queued
    """
    while True:
        try:
            target_queue.put(item, timeout=2)
            return True
        except queue.Full:
            if abort.is_set():
                return False


def reader_proc(paths, chunk_records, use_pigz, max_reads, worker_queues, result_queue, abort):
    """
    Reader process (one per FASTQ role).  Cuts the FASTQ files into chunks of chunk_records records
    and deals them round-robin to the worker queues.  The barcode reader and the target reader use
    the same rule, so chunk i of both roles reaches the same worker.  Sends None to every queue at
    the end.  Errors are sent to result_queue.
    Parameters:
            -- paths: list[str]: FASTQ files of this role, in order
            -- chunk_records: int: FASTQ records per chunk (--chunk)
            -- use_pigz: bool: allow pigz
            -- max_reads: int: stop after about this many reads (0 = all; pilot runs)
            -- worker_queues: list of multiprocessing queues: one per parse worker
            -- result_queue: multiprocessing queue: receives ('error', 'reader', traceback) on failure
            -- abort: multiprocessing.Event: set by the supervisor on any error
    Returns:
            -- None: nothing; results go through the queues
    """
    n_workers = len(worker_queues)
    try:
        # step1 : counters for the round-robin dealing and for --max_reads
        chunk_no = reads_done = 0
        stop = False
        # step2 : cut every FASTQ file into chunks of chunk_records records and deal them out
        for path in paths:
            #      |----   substep open the file (pigz when available)
            stream, proc = open_stream(path, use_pigz)
            try:
                while True:
                    #      |----   substep read 4 lines per record
                    lines = list(islice(stream, 4 * chunk_records))
                    if not lines:
                        break
                    #      |----   substep a partial record means a truncated file
                    if len(lines) % 4:
                        raise RuntimeError("truncated FASTQ: %s" % path)
                    n_records = len(lines) // 4
                    #      |----   substep chunk i goes to worker i % n_workers (the same rule is used for the barcode and the target reader)
                    if not put_abortable(worker_queues[chunk_no % n_workers],
                                         (n_records, b"".join(lines)), abort):
                        return
                    del lines
                    chunk_no += 1
                    reads_done += n_records
                    #      |----   substep pilot run: stop after about --max_reads reads
                    if max_reads and reads_done >= max_reads:
                        stop = True
                        break
            finally:
                close_stream(stream, proc)
            if stop:
                break
        # step3 : send the end-of-input marker (None) to every worker queue
        for q in worker_queues:                              # end-of-input marker
            if not put_abortable(q, None, abort):
                return
    except BaseException:
        result_queue.put(("error", "reader", traceback.format_exc()))
