#!/usr/bin/env python
# encoding: utf-8
"""
*Uncompress ESO fits.Z frames*

Author
: David Young

Date Created
: April 11, 2023
"""
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


def uncompress(
        log,
        directory):
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

    log.debug('starting the ``uncompress`` function')

    from subprocess import Popen, PIPE, STDOUT

    # GENERATE A LIST OF FILE PATHS
    count = 0
    for d in os.listdir(directory):
        filepath = os.path.join(directory, d)
        if os.path.isfile(filepath) and "fits" in d and os.path.splitext(filepath)[1] == ".Z":
            count += 1

    if count > 0:
        cmd = f"""uncompress {directory}/*fits.Z"""
        try:
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
            stdout, stderr = p.communicate()
            log.debug(f'output: {stdout}')
            if not stderr:
                print(f"Decompressed {count} fits.Z files")
        except Exception as e:
            log.error(f'Could not uncompress .Z files')

        if stderr and "uncompress" in stderr.decode('ascii'):
            print(f'The uncompress command was not found. Please install it or manually uncompress all `.Z` files before running `soxspipe prep` again.')
            sys.exit(0)

    log.debug('completed the ``uncompress`` function')
    return None
