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

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from soxspipe.commonutils.phase3 import MAX_VALUE_LENGTH

if TYPE_CHECKING:
    from astropy.io.fits import FITS_rec, HDUList, Header

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


def fail(message: str) -> Finding:
    """*return a failure finding*

    **Key Arguments:**

    - ``message`` -- what is wrong

    **Return:**

    - ``finding`` -- a ``Finding`` with the ``FAILURE`` severity

    **Usage:**

    ```python
    findings = [fail("EXPTIME is 0, it must be greater than 0")]
    ```
    """
    return Finding(FAILURE, message)


def warn(message: str) -> Finding:
    """*return a warning finding*

    **Key Arguments:**

    - ``message`` -- what is missing or doubtful

    **Return:**

    - ``finding`` -- a ``Finding`` with the ``WARNING`` severity

    **Usage:**

    ```python
    findings = [warn("recommended keyword SPEC_ERR is missing from the primary header")]
    ```
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

    path: Path
    hdus: HDUList
    primary: Header
    extension: Header | None
    table: FITS_rec | None
    tableError: str | None = None
    isBinTable: bool = False

    def header(self, hdu: str) -> Header | None:
        """*return the primary or the extension header, or None when the file has no extension*

        **Key Arguments:**

        - ``hdu`` -- ``PRIMARY`` or ``EXTENSION``

        **Return:**

        - ``header`` -- the header, or *None* when the file has no extension

        **Usage:**

        ```python
        extensionHeader = product.header(EXTENSION)
        ```
        """
        return self.primary if hdu == PRIMARY else self.extension


@dataclass(frozen=True)
class Check:
    """*one registered check*"""

    name: str
    function: Callable[[Phase3Product], list[Finding]]
    needsTable: bool
    needsTableHeader: bool


# EVERY REGISTERED CHECK, IN THE ORDER THEY RUN
CHECKS: dict[str, Check] = {}


def register_check(
    function: Callable[[Phase3Product], list[Finding]] | None = None,
    *,
    needsTable: bool = False,
    needsTableHeader: bool = False,
) -> Callable:
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

    def decorator(checkFunction: Callable[[Phase3Product], list[Finding]]) -> Callable:
        if checkFunction.__name__ in CHECKS:
            raise ValueError(f"A Phase 3 check named {checkFunction.__name__} is already registered")
        CHECKS[checkFunction.__name__] = Check(checkFunction.__name__, checkFunction, needsTable, needsTableHeader)
        return checkFunction

    return decorator if function is None else decorator(function)


# THE SUFFIX OF A PHASE 3 FITS FILE, AND THE STEM THAT MARKS A PRODUCT FOR ESO. ONLY THE `_ph3` PRODUCT IS BOUND BY
# ESO'S 68-CHARACTER LIMITS: THE INTERNAL MERGED SPECTRUM HAS A LONGER FILE NAME
FITS_SUFFIX = ".fits"
PHASE3_STEM_SUFFIX = "_ph3"
FAILED_SUFFIX = "_FAILED"
# THE FILE NAME IS THE VALUE OF ORIGFILE, SO THE LIMIT ON A FILE NAME (SECTION 3.4.1) IS THE LIMIT ON A VALUE
MAX_FILE_NAME_LENGTH = MAX_VALUE_LENGTH


def is_phase3_name(name: str) -> bool:
    """*return True when a file name is that of a Phase 3 product, ``<stem>_ph3.fits`` or ``<stem>_ph3_FAILED.fits``*

    **Key Arguments:**

    - ``name`` -- the file name, without its folder

    **Return:**

    - ``isPhase3`` -- True when the name ends in ``_ph3.fits`` or ``_ph3_FAILED.fits``

    **Usage:**

    ```python
    isPhase3 = is_phase3_name("SOXS_20260102T030405_STARE_ph3.fits")
    ```
    """
    stem = name[: -len(FITS_SUFFIX)] if name.endswith(FITS_SUFFIX) else None
    return stem is not None and stem.removesuffix(FAILED_SUFFIX).endswith(PHASE3_STEM_SUFFIX)


def is_phase3_product(product: Phase3Product) -> bool:
    """*return True when the file is named as a Phase 3 product, see ``is_phase3_name``*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to look at

    **Return:**

    - ``isPhase3`` -- True when the file name is that of a Phase 3 product

    **Usage:**

    ```python
    if is_phase3_product(product):
        findings += header_card_faults(product.primary, "primary")
    ```
    """
    return is_phase3_name(product.path.name)


def table_columns(header: Header) -> list[tuple[int, str, str]]:
    """*return the columns of a table header as ``(index, TTYPE, TFORM)``, the index counting from 1*

    **Key Arguments:**

    - ``header`` -- the header of a binary table extension

    **Return:**

    - ``columns`` -- one tuple per column, in file order, with an empty string for a missing card

    **Usage:**

    ```python
    names = [name for _, name, _ in table_columns(product.extension)]
    ```
    """
    return [
        (index, str(header.get(f"TTYPE{index}", "")), str(header.get(f"TFORM{index}", "")))
        for index in range(1, header.get("TFIELDS", 0) + 1)
    ]


def open_product(filePath: str | Path) -> Phase3Product:
    """*open a file for checking*

    **Key Arguments:**

    - ``filePath`` -- the path to the file

    **Return:**

    - ``product`` -- the ``Phase3Product``, with the file left open for the caller to close through ``product.hdus``

    **Raises:**

    - `OSError` when the file does not exist or is not a FITS file, `ValueError` when the path is empty

    **Usage:**

    ```python
    product = open_product("SOXS_20260102T030405_STARE_ph3.fits")
    ```
    """
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
