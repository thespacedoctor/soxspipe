#!/usr/bin/env python
"""
*Uncompress ESO fits.Z frames*

Author
: David Young

Date Created
: April 11, 2023
"""

import os
import sys

os.environ["TERM"] = "vt100"


BATCH_SIZE = 25
MISSING_COMMAND_MESSAGE = (
    "The uncompress command was not found. Please install it or manually "
    "uncompress all `.Z` files before running `soxspipe prep` again."
)


def _build_batches(directory):
    """list the .Z archives in a directory and group them into batches

    **Key Arguments:**

    - ``directory`` -- directory containing .Z file to uncompress

    **Return:**

    - ``batches`` -- lists of at most 25 archive paths, in sorted filename order
    - ``count`` -- the total number of archives found
    """
    count = 0
    batches = []
    batch = []
    for d in sorted(os.listdir(directory)):
        filepath = os.path.join(directory, d)
        if (
            os.path.isfile(filepath)
            and "fits" in d
            and os.path.splitext(filepath)[1] == ".Z"
        ):
            batch.append(filepath)
            count += 1
            if len(batch) == BATCH_SIZE:
                batches.append(batch)
                batch = []
    if len(batch) > 0:
        batches.append(batch)
    return batches, count


def _run_batch(log, batch, uncompressedCount, count):
    """run the uncompress command on one batch and report progress

    **Key Arguments:**

    - ``log`` -- logger
    - ``batch`` -- the archive paths to uncompress in this batch
    - ``uncompressedCount`` -- number of archives handled so far, including this batch
    - ``count`` -- the total number of archives found
    """
    from subprocess import PIPE, Popen

    cmd = ["uncompress", "-f", *batch]
    try:
        # THE COMMAND IS A FIXED ARGUMENT LIST, NEVER A SHELL STRING, SO S603 DOES NOT APPLY
        p = Popen(cmd, stdout=PIPE, stderr=PIPE)  # noqa: S603
        stdout, stderr = p.communicate()
        log.debug(f"output: {stdout}")
        if not stderr and p.returncode == 0:
            if uncompressedCount > len(batch):
                # CURSOR UP ONE LINE AND CLEAR LINE
                sys.stdout.flush()
                sys.stdout.write("\x1b[1A\x1b[2K")
            percent = (float(uncompressedCount) / float(count)) * 100.0
            print(
                f"Decompressed {uncompressedCount}/{count} fits.Z files ({percent:.1f}%)"
            )
    except FileNotFoundError:
        print(MISSING_COMMAND_MESSAGE)
        sys.exit(0)
    except OSError as error:
        log.error(f"Could not uncompress .Z files: {error}")
        return

    stderrMessage = stderr.decode("ascii", errors="replace")
    if p.returncode == 127:
        print(stderrMessage)
        print(MISSING_COMMAND_MESSAGE)
        sys.exit(0)
    if p.returncode:
        log.error(
            f"Could not uncompress .Z files (exit code {p.returncode}): "
            f"{stderrMessage}"
        )


def uncompress(log, directory):
    """uncompress ESO fits.Z frames

    **Key Arguments:**

    - ``log`` -- logger
    - ``directory`` -- directory containing .Z file to uncompress

    ```python
    from soxspipe.commonutils import uncompress
    uncompress(
        log=log,
        directory="/path/to/raw_data/"
    )
    ```
    """

    log.debug("starting the ``uncompress`` function")

    # GENERATE A LIST OF FILE PATHS
    batches, count = _build_batches(directory)

    uncompressedCount = 0
    for batch in batches:
        uncompressedCount += len(batch)
        _run_batch(log, batch, uncompressedCount, count)

    log.debug("completed the ``uncompress`` function")
    return
