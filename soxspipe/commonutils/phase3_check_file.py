#!/usr/bin/env python
"""
*the checks of the file: its name, its HDU structure, its checksums and the ancillary files it links*

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import TYPE_CHECKING

from soxspipe.commonutils.phase3_check_framework import (
    EXTENSION,
    FITS_SUFFIX,
    MAX_FILE_NAME_LENGTH,
    PRIMARY,
    Finding,
    Phase3Product,
    fail,
    is_phase3_product,
    register_check,
    table_columns,
)

if TYPE_CHECKING:
    from astropy.io.fits import Header

# THE SIZE OF THE PIECES A LINKED FILE IS READ IN WHEN ITS MD5 SUM IS COMPUTED
MD5_CHUNK_BYTES = 1024 * 1024


@register_check
def check_file_name(product: Phase3Product) -> list[Finding]:
    """*report a file name that does not end in ``.fits``, or a ``_ph3`` file name over 68 characters*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_file_name(product)
    ```
    """
    name = product.path.name
    findings = []
    if not name.endswith(FITS_SUFFIX):
        findings.append(fail(f"the file name {name} does not end in {FITS_SUFFIX}"))
    if is_phase3_product(product) and len(name) > MAX_FILE_NAME_LENGTH:
        findings.append(
            fail(f"the file name is {len(name)} characters, over the {MAX_FILE_NAME_LENGTH} that ESO accepts")
        )
    return findings


@register_check
def check_hdu_structure(product: Phase3Product) -> list[Finding]:
    """*report a file that is not an empty primary HDU and one single-row binary table with no variable-length array*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_hdu_structure(product)
    ```
    """
    findings = []
    if product.primary.get("NAXIS", 0) != 0:
        findings.append(fail(f"the primary HDU holds data, NAXIS is {product.primary['NAXIS']} and must be 0"))
    extensionCount = len(product.hdus) - 1
    if extensionCount != 1:
        findings.append(fail(f"the file must have exactly one extension, it has {extensionCount}"))
    if product.extension is None:
        return findings
    if product.tableError:
        findings.append(fail(f"the table cannot be read: {product.tableError}"))
    elif product.table is None:
        findings.append(fail("the first extension is not a BINTABLE"))
    if product.table is None:
        return findings
    header = product.extension
    for keyword, expected in (("NAXIS2", 1), ("PCOUNT", 0), ("GCOUNT", 1)):
        if header.get(keyword) != expected:
            findings.append(fail(f"{keyword} of the table is {header.get(keyword)}, it must be {expected}"))
    variableColumns = [name for _, name, form in table_columns(header) if re.match(r"\d*[PQ]", form)]
    if variableColumns:
        findings.append(
            fail(f"the columns {variableColumns} are variable-length arrays, which Phase 3 does not accept")
        )
    return findings


@register_check
def check_checksums(product: Phase3Product) -> list[Finding]:
    """*report an HDU whose ``CHECKSUM`` or ``DATASUM`` does not match its content*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_checksums(product)
    ```
    """
    findings = []
    for hdu, name in zip(product.hdus, (PRIMARY, EXTENSION), strict=False):
        for keyword, verify in (("CHECKSUM", hdu.verify_checksum), ("DATASUM", hdu.verify_datasum)):
            if keyword in hdu.header and verify() != 1:
                findings.append(fail(f"{keyword} of the {name} header does not match the content"))
    return findings


# ANCILLARY FILES ARE LINKED FROM THE PRIMARY HEADER BY ASSONi (THE FILE NAME), ASSOCi (THE CATEGORY OF A NON-FITS FILE)
# AND ASSOMi (ITS MD5 SUM)
ASSOCIATION_PATTERN = r"ASSO([NCM])(\d+)"
ANCILLARY_PREFIX = "ANCILLARY."


def _associations(header: Header) -> tuple[set[int], set[int], set[int]]:
    """*return the indexes of the ``ASSON``, ``ASSOC`` and ``ASSOM`` cards of a header as three sets*"""
    found: dict[str, set[int]] = {"N": set(), "C": set(), "M": set()}
    for keyword in header:
        match = re.fullmatch(ASSOCIATION_PATTERN, keyword)
        if match:
            found[match.group(1)].add(int(match.group(2)))
    return found["N"], found["C"], found["M"]


def _md5_of_file(path: Path) -> str:
    """*return the md5 sum of a file, read in pieces so a large file is never held in memory*"""
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(MD5_CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _linked_path(product: Phase3Product, index: int, name: object) -> tuple[Path | None, list[Finding]]:
    """*return the path of the file ``ASSON<index>`` names, or None and the failure when it may not be read*

    The value must be a bare file name. A value with a folder in it, an absolute path or a link that leads outside
    the product's folder is refused without opening the file.
    """
    if not isinstance(name, str) or not name or name in (".", "..") or Path(name).name != name:
        return None, [fail(f"ASSON{index} is {name!r}, it must be a bare file name of a file beside the product")]
    folder = product.path.parent.resolve()
    path = (folder / name).resolve()
    if path.parent != folder:
        return None, [fail(f"ASSON{index} names {name!r}, which leads outside the folder of the product")]
    if not path.is_file():
        return None, [fail(f"ASSON{index} names {name!r}, which is not found beside the product")]
    return path, []


def _linked_fits_findings(index: int, path: Path) -> list[Finding]:
    """*return the findings for a linked FITS file: an ancillary ``PRODCATG`` and no links of its own*"""
    from astropy.io import fits

    try:
        header = fits.getheader(path)
    except (OSError, ValueError) as error:
        return [fail(f"ASSON{index} names a FITS file that cannot be opened: {error}")]
    findings = []
    prodcatg = str(header.get("PRODCATG", ""))
    if not prodcatg.startswith(ANCILLARY_PREFIX):
        findings.append(
            fail(f"ASSON{index} is a FITS file with PRODCATG {prodcatg!r}, it must start with {ANCILLARY_PREFIX}")
        )
    if _associations(header)[0]:
        findings.append(
            fail(f"ASSON{index} is a FITS file that links files in turn, nested associations are not allowed")
        )
    return findings


def _linked_file_findings(product: Phase3Product, index: int, name: object) -> list[Finding]:
    """*return the findings for the file that ``ASSON<index>`` names*"""
    path, findings = _linked_path(product, index, name)
    if path is None:
        return findings
    if len(name) > MAX_FILE_NAME_LENGTH:
        findings.append(fail(f"the file name {name!r} of ASSON{index} is over {MAX_FILE_NAME_LENGTH} characters"))
    if path.suffix == FITS_SUFFIX:
        return findings + _linked_fits_findings(index, path)
    category = product.primary.get(f"ASSOC{index}")
    if not str(category).startswith(ANCILLARY_PREFIX):
        findings.append(
            fail(f"ASSOC{index} is {category!r}, a non-FITS file needs a category that starts with {ANCILLARY_PREFIX}")
        )
    md5 = product.primary.get(f"ASSOM{index}")
    if md5 is None:
        findings.append(fail(f"ASSOM{index} is missing, a non-FITS file needs its md5 sum"))
    elif md5 != _md5_of_file(path):
        findings.append(fail(f"ASSOM{index} is not the md5 sum of {name!r}"))
    return findings


@register_check
def check_ancillary_links(product: Phase3Product) -> list[Finding]:
    """*report an ASSON, ASSOC or ASSOM card that does not describe a file beside the product*

    An ``ASSON`` value must be a bare file name of a file in the product's own folder. The file is never read from
    anywhere else, so a header cannot point the checker at an arbitrary file.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_ancillary_links(product)
    ```
    """
    linked, categories, sums = _associations(product.primary)
    findings = [
        fail(f"{prefix}{index} has no matching ASSON{index}")
        for prefix, indexes in (("ASSOC", categories), ("ASSOM", sums))
        for index in sorted(indexes - linked)
    ]
    for index in sorted(linked):
        findings.extend(_linked_file_findings(product, index, product.primary[f"ASSON{index}"]))
    return findings
