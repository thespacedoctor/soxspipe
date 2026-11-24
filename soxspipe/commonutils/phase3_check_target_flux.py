#!/usr/bin/env python
"""
*the checks of the target position and flux-calibration keywords of a Phase 3 spectrum*

The standard already makes these keywords mandatory, and ``phase3_check_headers`` checks their types, ``RADESYS``,
``EQUINOX``, and ``FLUXERR`` of an uncalibrated spectrum. These checks add the house rules of the spec (DY-2393): an
absolute spectrum has a ``FLUXERR`` of -2, ``TOT_FLUX`` is always false because no slit-loss correction is applied, and
the target position is on the sky.

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

import math

from soxspipe.commonutils.phase3_check_framework import EXTENSION, PRIMARY, Finding, Phase3Product, fail, register_check
from soxspipe.commonutils.phase3_check_headers import HOUSE_RULE

ABSOLUTE_FLUXERR = -2
RA_RANGE = (0.0, 360.0)
DEC_RANGE = (-90.0, 90.0)


def _is_number(value: object) -> bool:
    """*return True for a finite int or float, a bool never passing as a number*"""
    return isinstance(value, int | float) and not isinstance(value, bool) and math.isfinite(value)


@register_check
def check_absolute_spectrum_has_the_absolute_flux_error(product: Phase3Product) -> list[Finding]:
    """*report an absolute spectrum whose FLUXERR is not -2*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_absolute_spectrum_has_the_absolute_flux_error(product)
    ```
    """
    primary = product.primary
    if primary.get("FLUXCAL") != "ABSOLUTE" or "FLUXERR" not in primary:
        return []
    fluxerr = primary["FLUXERR"]
    if _is_number(fluxerr) and fluxerr == ABSOLUTE_FLUXERR:
        return []
    return [fail(f"FLUXERR is {fluxerr}, it must be {ABSOLUTE_FLUXERR} when FLUXCAL is 'ABSOLUTE' ({HOUSE_RULE})")]


@register_check
def check_total_flux_is_false(product: Phase3Product) -> list[Finding]:
    """*report a TOT_FLUX that is true, because no slit-loss correction is applied*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_total_flux_is_false(product)
    ```
    """
    if product.primary.get("TOT_FLUX") is True:
        return [fail(f"TOT_FLUX is True, it must be False: the spectrum has no slit-loss correction ({HOUSE_RULE})")]
    return []


@register_check
def check_target_position_is_on_the_sky(product: Phase3Product) -> list[Finding]:
    """*report an RA outside 0 to 360 deg, or a DEC outside -90 to 90 deg, in the primary or extension header*

    A value that is not a number is left to ``check_keyword_value_types``.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_target_position_is_on_the_sky(product)
    ```
    """
    findings = []
    for hduName in (PRIMARY, EXTENSION):
        header = product.header(hduName)
        if header is None:
            continue
        for keyword, (low, high), isHighIncluded in (("RA", RA_RANGE, False), ("DEC", DEC_RANGE, True)):
            value = header.get(keyword)
            if not _is_number(value):
                continue
            if value < low or value > high or (value == high and not isHighIncluded):
                closing = "]" if isHighIncluded else ")"
                findings.append(
                    fail(
                        f"{keyword} is {value} in the {hduName} header, outside its allowed range (house rule): "
                        f"it must be in [{low:g}, {high:g}{closing} deg"
                    )
                )
    return findings
