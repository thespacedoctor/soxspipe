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

import contextlib
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
    raises ``OSError`` unless ``overwrite`` is set. A symbolic link at ``filePath`` is kept and its target is
    replaced.

    The write protects against a killed process, not against power loss (nothing is synced to disk). A process
    killed with ``SIGKILL`` leaves its ``.partial`` file behind; it is harmless and nothing sweeps it. The
    ``overwrite=False`` check is best effort when two processes write the same path at once, which soxspipe does not
    support (one process per workspace).

    **Usage:**

    ```python
    from soxspipe.commonutils.atomic_fits import write_hdulist_atomically
    write_hdulist_atomically(hduList, filePath, overwrite=True, checksum=True)
    ```
    """
    overwrite = writetoKwargs.pop("overwrite", False)
    if os.path.exists(filePath) and not overwrite:
        raise OSError(f"File {filePath!r} already exists.")

    # WRITE BESIDE THE REAL FILE, SO A SYMLINK AT THE PRODUCT PATH SURVIVES THE RENAME
    filePath = os.path.realpath(filePath)
    directory, basename = os.path.split(filePath)
    temporaryPath = os.path.join(directory, f".{basename}.{os.getpid()}.{uuid.uuid4().hex[:8]}{TEMPORARY_SUFFIX}")
    try:
        hduList.writeto(temporaryPath, overwrite=True, **writetoKwargs)
        os.replace(temporaryPath, filePath)
    except BaseException:
        with contextlib.suppress(OSError):
            os.remove(temporaryPath)
        raise


def is_complete_fits(filePath):
    """*return True when the file opens as FITS and ends exactly where its last HDU says it should*

    **Key Arguments:**

    - ``filePath`` -- the path of the file to test

    **Return:**

    - ``isComplete`` -- *False* when the file is empty, is not FITS, cannot be read, or its size is not the end of
      its last HDU

    The test compares the file size with the end of the last HDU (data start plus data span, padding included), so
    no warning text is involved. A file cut exactly at the end of an HDU looks whole, since FITS records no HDU
    count. ``warnings.catch_warnings`` is not thread-safe, so do not call this from threads that also rely on
    warning filters.

    **Usage:**

    ```python
    from soxspipe.commonutils.atomic_fits import is_complete_fits
    if not is_complete_fits(productPath):
        ...
    ```
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            with fits.open(filePath) as hduList:
                lastHdu = None
                for lastHdu in hduList:
                    lastHdu.header  # noqa: B018
                if lastHdu is None:
                    return False
                info = lastHdu.fileinfo()
                return os.path.getsize(filePath) == info["datLoc"] + info["datSpan"]
        except Exception:
            return False
