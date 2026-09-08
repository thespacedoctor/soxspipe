#!/usr/bin/env python
# encoding: utf-8
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

    from subprocess import PIPE, Popen

    # GENERATE A LIST OF FILE PATHS
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
            if len(batch) == 25:
                batches.append(batch)
                batch = []
    if len(batch) > 0:
        batches.append(batch)

    uncompressedCount = 0
    missingCommandMessage = (
        "The uncompress command was not found. Please install it or manually "
        "uncompress all `.Z` files before running `soxspipe prep` again."
    )
    for batch in batches:
        uncompressedCount += len(batch)
        cmd = ["uncompress", "-f", *batch]
        try:
            p = Popen(cmd, stdout=PIPE, stderr=PIPE)
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
            print(missingCommandMessage)
            sys.exit(0)
        except OSError as error:
            log.error(f"Could not uncompress .Z files: {error}")
            continue

        stderrMessage = stderr.decode("ascii", errors="replace")
        if p.returncode == 127:
            print(stderrMessage)
            print(missingCommandMessage)
            sys.exit(0)
        if p.returncode:
            log.error(
                f"Could not uncompress .Z files (exit code {p.returncode}): "
                f"{stderrMessage}"
            )

    log.debug("completed the ``uncompress`` function")
    return None
