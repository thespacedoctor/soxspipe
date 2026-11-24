#!/usr/bin/env python
"""
*the checks of the spectrum table: its columns, their descriptors and the arrays in them*

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from soxspipe.commonutils.phase3_check_framework import (
    Finding,
    Phase3Product,
    fail,
    register_check,
    table_columns,
)

if TYPE_CHECKING:
    from astropy.io.fits import Header

# THE FIRST THREE COLUMNS OF A SPECTRUM, IN ORDER, AND THE TOKEN THAT MARKS THE MAIN FLUX AND ERROR IN THEIR UCDS
REQUIRED_COLUMNS = ("WAVE", "FLUX", "ERR")
MAIN_UCD_TAG = "meta.main"
VOCLASS_VALUE = "SPECTRUM v2.0"
WAVE_UCDS = ("em.wl;obs.atmos", "em.wl")
EXPOSURE_UCD = "time.duration;obs.exposure"

# THE UTYPES OF THE MAIN COLUMNS FOR VOCLASS `SPECTRUM v2.0`
MAIN_UTYPES = {
    "WAVE": "spec:Data.SpectralAxis.Value",
    "FLUX": "spec:Data.FluxAxis.Value",
    "ERR": "spec:Data.FluxAxis.Accuracy.StatError",
}


def _column_indexes(header: Header) -> dict[str, int]:
    """*return the 1-based index of each column name in a table header, the first one when a name repeats*"""
    indexes: dict[str, int] = {}
    for index, name, _ in table_columns(header):
        indexes.setdefault(name, index)
    return indexes


@register_check(needsTableHeader=True)
def check_column_names(product: Phase3Product) -> list[Finding]:
    """*report a table that does not start with WAVE, FLUX, ERR, or that repeats a column name*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_column_names(product)
    ```
    """
    header = product.extension
    names = [name for _, name, _ in table_columns(header)]
    findings = []
    if tuple(names[:3]) != REQUIRED_COLUMNS:
        findings.append(fail(f"the first three columns must be WAVE, FLUX, ERR in that order, not {names[:3]}"))
    findings.extend(
        fail(f"the column {name} appears more than once") for name in sorted(set(names)) if names.count(name) > 1
    )
    return findings


def _column_descriptor_findings(header: Header, name: str, index: int) -> list[Finding]:
    """*return the findings for the unit, UCD and UType cards of one column*"""
    findings = []
    for prefix in ("TUNIT", "TUCD", "TUTYP"):
        if f"{prefix}{index}" not in header:
            findings.append(fail(f"{prefix}{index} is missing for the column {name}"))
    ucd = header.get(f"TUCD{index}")
    if name == "WAVE" and ucd is not None and ucd not in WAVE_UCDS:
        findings.append(fail(f"TUCD{index} of WAVE is {ucd!r}, it must be em.wl;obs.atmos (air) or em.wl (vacuum)"))
    if name in ("FLUX", "ERR") and ucd is not None and MAIN_UCD_TAG not in ucd.split(";"):
        findings.append(fail(f"TUCD{index} of {name} is {ucd!r}, it must hold the token {MAIN_UCD_TAG}"))
    utype = header.get(f"TUTYP{index}")
    if name in MAIN_UTYPES and utype is not None and utype != MAIN_UTYPES[name]:
        findings.append(fail(f"TUTYP{index} of {name} is {utype!r}, it must be {MAIN_UTYPES[name]!r}"))
    if name == "EXPOSURE":
        if ucd is not None and ucd != EXPOSURE_UCD:
            findings.append(fail(f"TUCD{index} of EXPOSURE is {ucd!r}, it must be {EXPOSURE_UCD!r}"))
        if utype:
            findings.append(fail(f"TUTYP{index} of EXPOSURE is {utype!r}, it must be an empty string"))
    return findings


def _same_as_flux_findings(header: Header, indexes: dict[str, int]) -> list[Finding]:
    """*return the findings for an ERR or BGFLUX column that does not share the unit (and UCD) of FLUX*"""
    findings = []
    if "FLUX" not in indexes:
        return findings
    flux = indexes["FLUX"]
    fluxUnit = header.get(f"TUNIT{flux}")
    fluxUcd = [token for token in str(header.get(f"TUCD{flux}", "")).split(";") if token != MAIN_UCD_TAG]
    for name in ("ERR", "BGFLUX"):
        index = indexes.get(name)
        if index is None:
            continue
        if header.get(f"TUNIT{index}") != fluxUnit:
            findings.append(fail(f"the unit of {name} is {header.get(f'TUNIT{index}')!r}, it must be the unit of FLUX"))
        ucd = [token for token in str(header.get(f"TUCD{index}", "")).split(";") if token != MAIN_UCD_TAG]
        if name == "BGFLUX" and ucd != fluxUcd:
            findings.append(fail("the TUCD of BGFLUX must be the TUCD of FLUX, without the meta.main token"))
    return findings


@register_check(needsTableHeader=True)
def check_column_descriptors(product: Phase3Product) -> list[Finding]:
    """*report a column without TUNIT, TUCD or TUTYP, with a wrong UCD or UType, or a wrong VOCLASS*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_column_descriptors(product)
    ```
    """
    header = product.extension
    indexes = _column_indexes(header)
    findings = []
    if header.get("VOCLASS", VOCLASS_VALUE) != VOCLASS_VALUE:
        findings.append(fail(f"VOCLASS is {header['VOCLASS']!r}, it must be {VOCLASS_VALUE!r}"))
    for name, index in indexes.items():
        findings.extend(_column_descriptor_findings(header, name, index))
    findings.extend(_same_as_flux_findings(header, indexes))
    if "WAVE" in indexes:
        findings.extend(
            fail(f"{prefix}{indexes['WAVE']} is missing, the WAVE column needs its range")
            for prefix in ("TDMIN", "TDMAX")
            if f"{prefix}{indexes['WAVE']}" not in header
        )
    return findings


@register_check(needsTable=True)
def check_table_arrays(product: Phase3Product) -> list[Finding]:
    """*report arrays of unequal length, a NELEM that is not their length, a non-integer QUAL or a non-finite value*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_table_arrays(product)
    ```
    """
    import numpy as np

    if len(product.table) == 0:
        return []
    arrays = {name: np.atleast_1d(product.table[name][0]) for name in product.table.names}
    first = next(iter(arrays.values()), np.array([]))
    findings = [
        fail(f"the column {name} has {array.size} elements, the first column has {first.size}")
        for name, array in arrays.items()
        if array.size != first.size
    ]
    nelem = product.extension.get("NELEM")
    if nelem is not None and nelem != first.size:
        findings.append(fail(f"NELEM is {nelem}, the arrays have {first.size} elements"))
    if "QUAL" in arrays and not np.issubdtype(arrays["QUAL"].dtype, np.integer):
        findings.append(fail("the QUAL column must hold integers"))
    findings.extend(
        fail(f"the column {name} holds NaN or infinite values")
        for name, array in arrays.items()
        if np.issubdtype(array.dtype, np.floating) and not np.isfinite(array).all()
    )
    return findings


@register_check(needsTable=True)
def check_wavelength_axis(product: Phase3Product) -> list[Finding]:
    """*report a WAVE array that does not strictly increase, or that lies outside TDMIN and TDMAX*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_wavelength_axis(product)
    ```
    """
    import numpy as np

    indexes = _column_indexes(product.extension)
    if "WAVE" not in indexes or len(product.table) == 0:
        return []
    wave = np.atleast_1d(product.table["WAVE"][0])
    findings = []
    if wave.size > 1 and not (np.diff(wave) > 0).all():
        findings.append(fail("WAVE does not strictly increase"))
    index = indexes["WAVE"]
    tdmin, tdmax = product.extension.get(f"TDMIN{index}"), product.extension.get(f"TDMAX{index}")
    if tdmin is not None and wave.min() < tdmin:
        findings.append(fail(f"WAVE holds a value below TDMIN{index} ({tdmin})"))
    if tdmax is not None and wave.max() > tdmax:
        findings.append(fail(f"WAVE holds a value above TDMAX{index} ({tdmax})"))
    return findings
