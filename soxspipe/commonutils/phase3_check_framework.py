#!/usr/bin/env python
"""
*the records, the registry and the file reader behind the Phase 3 checker*

``phase3_check`` is the public entry point. This module holds what every check module needs, so a check module can
import it without a circular import: the ``Finding`` record, ``fail`` and ``warn``, the ``Phase3Product`` a check
reads, and the ``register_check`` registry.

Author
: David Young

Date Created
: November 24, 2026
"""

from dataclasses import dataclass

FAILURE = "failure"
WARNING = "warning"

PRIMARY = "primary"
EXTENSION = "extension"


@dataclass(frozen=True)
class Finding:
    """*one fault or doubt that a check found in a product*

    **Key Arguments:**

    - ``severity`` -- ``FAILURE`` or ``WARNING``
    - ``message`` -- what is wrong, in one sentence
    - ``check`` -- the name of the check that found it, set by the runner
    """

    severity: str
    message: str
    check: str = ""


def fail(message):
    """*return a failure finding*

    **Key Arguments:**

    - ``message`` -- what is wrong

    **Return:**

    - ``finding`` -- a ``Finding`` with the ``FAILURE`` severity
    """
    return Finding(FAILURE, message)


def warn(message):
    """*return a warning finding*

    **Key Arguments:**

    - ``message`` -- what is missing or doubtful

    **Return:**

    - ``finding`` -- a ``Finding`` with the ``WARNING`` severity
    """
    return Finding(WARNING, message)


@dataclass(frozen=True)
class Phase3Product:
    """*the parts of one opened file that the checks read*

    **Key Arguments:**

    - ``path`` -- the ``Path`` of the file
    - ``hdus`` -- the open ``HDUList``
    - ``primary`` -- the primary header
    - ``extension`` -- the header of the first extension, or *None* when the file has none
    - ``table`` -- the first extension's table data, or *None* when it is not a binary table or cannot be read
    - ``tableError`` -- why a binary table cannot be read, or *None*
    - ``isBinTable`` -- the first extension is a binary table, whether or not its data can be read
    """

    path: object
    hdus: object
    primary: object
    extension: object
    table: object
    tableError: object = None
    isBinTable: bool = False

    def header(self, hdu):
        """*return the primary or the extension header, or None when the file has no extension*

        **Key Arguments:**

        - ``hdu`` -- ``PRIMARY`` or ``EXTENSION``

        **Return:**

        - ``header`` -- the header, or *None* when the file has no extension
        """
        return self.primary if hdu == PRIMARY else self.extension


@dataclass(frozen=True)
class Check:
    """*one registered check*"""

    name: str
    function: object
    needsTable: bool
    needsTableHeader: bool


# EVERY REGISTERED CHECK, IN THE ORDER THEY RUN
CHECKS = {}


def register_check(function=None, *, needsTable=False, needsTableHeader=False):
    """*register a function as a Phase 3 check, as a plain or a called decorator*

    The check is named after the function. It takes a ``Phase3Product`` and returns a list of ``Finding`` records.

    **Key Arguments:**

    - ``function`` -- the check function, when the decorator is used without a call
    - ``needsTable`` -- the check reads the table data, so the runner skips it when the file has no readable table.
      Default *False*
    - ``needsTableHeader`` -- the check reads the table's column cards but not its data, so the runner skips it only
      when the first extension is not a binary table. Default *False*

    **Return:**

    - ``function`` -- the function, unchanged, or a decorator when called with arguments

    **Raises:**

    - `ValueError` when a check of that name is already registered

    **Usage:**

    ```python
    @register_check(needsTable=True)
    def check_flux_is_not_all_zero(product):
        return []
    ```
    """

    def decorator(checkFunction):
        if checkFunction.__name__ in CHECKS:
            raise ValueError(f"A Phase 3 check named {checkFunction.__name__} is already registered")
        CHECKS[checkFunction.__name__] = Check(checkFunction.__name__, checkFunction, needsTable, needsTableHeader)
        return checkFunction

    return decorator if function is None else decorator(function)


# THE SUFFIX OF A PHASE 3 FITS FILE, AND THE STEM THAT MARKS A PRODUCT FOR ESO. ONLY THE `_ph3` PRODUCT IS BOUND BY
# ESO'S 68-CHARACTER LIMITS: THE INTERNAL MERGED SPECTRUM HAS A LONGER FILE NAME
FITS_SUFFIX = ".fits"
PHASE3_STEM_SUFFIX = "_ph3"
MAX_FILE_NAME_LENGTH = 68


def is_phase3_product(product):
    """*return True when the file is named as a Phase 3 product, ending in ``_ph3.fits``*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to look at

    **Return:**

    - ``isPhase3`` -- True when the file name ends in ``_ph3.fits``
    """
    return product.path.name.endswith(PHASE3_STEM_SUFFIX + FITS_SUFFIX)


def open_product(filePath):
    """*open a file for checking*

    **Key Arguments:**

    - ``filePath`` -- the path to the file

    **Return:**

    - ``product`` -- the ``Phase3Product``, with the file left open for the caller to close through ``product.hdus``

    **Raises:**

    - `OSError` when the file does not exist or is not a FITS file, `ValueError` when the path is empty
    """
    from pathlib import Path

    from astropy.io import fits

    hdus = fits.open(filePath, checksum=False)
    extension = hdus[1] if len(hdus) > 1 else None
    table, tableError = None, None
    if isinstance(extension, fits.BinTableHDU):
        try:
            table = extension.data
        except (ValueError, OSError) as error:
            tableError = str(error)
    return Phase3Product(
        path=Path(filePath),
        hdus=hdus,
        primary=hdus[0].header,
        extension=extension.header if extension is not None else None,
        table=table,
        tableError=tableError,
        isBinTable=isinstance(extension, fits.BinTableHDU),
    )
