#!/usr/bin/env python
"""
*the checks of the file: its name, its HDU structure, its checksums and the ancillary files it links*

Author
: David Young

Date Created
: November 24, 2026
"""

from soxspipe.commonutils.phase3_check_framework import (
    EXTENSION,
    FITS_SUFFIX,
    MAX_FILE_NAME_LENGTH,
    PRIMARY,
    fail,
    is_phase3_product,
    register_check,
)


@register_check
def check_file_name(product):
    """*report a file name that does not end in ``.fits``, or a ``_ph3`` file name over 68 characters*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
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
def check_hdu_structure(product):
    """*report a file that is not an empty primary HDU and one single-row binary table with no variable-length array*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
    """
    import re

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
    variableColumns = [
        header.get(f"TTYPE{index}")
        for index in range(1, header.get("TFIELDS", 0) + 1)
        if re.match(r"\d*[PQ]", str(header.get(f"TFORM{index}", "")))
    ]
    if variableColumns:
        findings.append(
            fail(f"the columns {variableColumns} are variable-length arrays, which Phase 3 does not accept")
        )
    return findings


@register_check
def check_checksums(product):
    """*report an HDU whose ``CHECKSUM`` or ``DATASUM`` does not match its content*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
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


def _associations(header):
    """*return the indexes of the ``ASSON``, ``ASSOC`` and ``ASSOM`` cards of a header as three sets*"""
    import re

    found = {"N": set(), "C": set(), "M": set()}
    for keyword in header:
        match = re.fullmatch(ASSOCIATION_PATTERN, keyword)
        if match:
            found[match.group(1)].add(int(match.group(2)))
    return found["N"], found["C"], found["M"]


def _linked_file_findings(product, index, name):
    """*return the findings for the file that ``ASSON<index>`` names*"""
    import hashlib

    from astropy.io import fits

    path = product.path.parent / name
    if not path.is_file():
        return [fail(f"ASSON{index} names {name!r}, which is not found beside the product")]
    findings = []
    if len(name) > MAX_FILE_NAME_LENGTH:
        findings.append(fail(f"the file name {name!r} of ASSON{index} is over {MAX_FILE_NAME_LENGTH} characters"))
    if path.suffix == FITS_SUFFIX:
        prodcatg = str(fits.getheader(path).get("PRODCATG", ""))
        if not prodcatg.startswith(ANCILLARY_PREFIX):
            findings.append(
                fail(f"ASSON{index} is a FITS file with PRODCATG {prodcatg!r}, it must start with {ANCILLARY_PREFIX}")
            )
        if _associations(fits.getheader(path))[0]:
            findings.append(
                fail(f"ASSON{index} is a FITS file that links files in turn, nested associations are not allowed")
            )
        return findings
    category = product.primary.get(f"ASSOC{index}")
    if not str(category).startswith(ANCILLARY_PREFIX):
        findings.append(
            fail(f"ASSOC{index} is {category!r}, a non-FITS file needs a category that starts with {ANCILLARY_PREFIX}")
        )
    md5 = product.primary.get(f"ASSOM{index}")
    if md5 is None:
        findings.append(fail(f"ASSOM{index} is missing, a non-FITS file needs its md5 sum"))
    elif md5 != hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest():
        findings.append(fail(f"ASSOM{index} is not the md5 sum of {name!r}"))
    return findings


@register_check
def check_ancillary_links(product):
    """*report an ASSON, ASSOC or ASSOM card that does not describe a file beside the product*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
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
