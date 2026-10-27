#!/usr/bin/env python
"""
*write FITS products so that an interrupted run never leaves a partial file at the final path*

A product is written to a temporary name in the directory of the final path and renamed onto it when the write is
complete. A reader of the final path sees the previous file or the whole new one, never half of either.

This module imports nothing from soxspipe, so any module can import it without creating an import cycle.

Author
: David Young

Date Created
: October 27, 2026
"""

import os
import uuid
import warnings

from astropy.io import fits

TEMPORARY_SUFFIX = ".partial"


def write_hdulist_atomically(hduList, filePath, **writetoKwargs):
    """*write an HDUList to a temporary file beside ``filePath`` and rename it into place*

    **Key Arguments:**

    - ``hduList`` -- the astropy HDUList to write
    - ``filePath`` -- the final path of the file
    - ``writetoKwargs`` -- keyword arguments for ``HDUList.writeto``. ``overwrite`` is *False* unless set

    A failed write removes the temporary file and leaves ``filePath`` as it was. A file already at ``filePath``
    raises ``OSError`` unless ``overwrite`` is set.

    **Usage:**

    ```python
    from soxspipe.commonutils.atomic_fits import write_hdulist_atomically
    write_hdulist_atomically(hduList, filePath, overwrite=True, checksum=True)
    ```
    """
    overwrite = writetoKwargs.pop("overwrite", False)
    if os.path.exists(filePath) and not overwrite:
        raise OSError(f"File {filePath!r} already exists.")

    directory, basename = os.path.split(filePath)
    temporaryPath = os.path.join(directory, f".{basename}.{os.getpid()}.{uuid.uuid4().hex[:8]}{TEMPORARY_SUFFIX}")
    try:
        hduList.writeto(temporaryPath, overwrite=True, **writetoKwargs)
        os.replace(temporaryPath, filePath)
    except BaseException:
        if os.path.exists(temporaryPath):
            os.remove(temporaryPath)
        raise


def is_complete_fits(filePath):
    """*return True when the file opens as FITS and holds all the data its headers promise*

    **Key Arguments:**

    - ``filePath`` -- the path of the file to test

    **Return:**

    - ``isComplete`` -- *False* when the file is empty, is not FITS, or is shorter than its headers say

    **Usage:**

    ```python
    from soxspipe.commonutils.atomic_fits import is_complete_fits
    if not is_complete_fits(productPath):
        ...
    ```
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            with fits.open(filePath) as hduList:
                for hdu in hduList:
                    hdu.header  # noqa: B018
        except (OSError, ValueError, TypeError):
            return False

    return not any("truncated" in str(warning.message) for warning in caught)
