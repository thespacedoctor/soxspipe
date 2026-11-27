#!/usr/bin/env python
"""
*the checks of the wavelength-calibration and resolution keywords of a Phase 3 spectrum*

The standard already rates ``SPEC_RES`` as mandatory, ``SPEC_ERR`` and ``SPEC_SYE`` as recommended, and ``LAMRMS`` and
``LAMNLIN`` as optional. ``phase3_check_headers`` checks their presence, their type, and that ``SPEC_RES`` is above 0.
These checks add house rules of this ticket (DY-2358): ``LAMRMS`` and ``SPEC_ERR`` are above 0, ``LAMNLIN`` is at
least 1, and ``SPEC_ERR = LAMRMS / sqrt(LAMNLIN)``. ESO's own X-Shooter products use ``LAMRMS / LAMNLIN``, which the
standard does not say.

``SPEC_SYE`` is recommended, so a product without it (every NIR product until its sky-shift error formulas exist,
DY-2359) gets a warning from the presence check and no failure. The standard lets ``SPEC_SYE`` be 0 when the systematic
error is insignificant, so the house rule of DY-2359 is only that it is not negative.

Author
: David Young

Date Created
: November 27, 2026
"""

from __future__ import annotations

import math

from soxspipe.commonutils.phase3_check_framework import PRIMARY, Finding, Phase3Product, fail, register_check

# WHERE THESE RULES COME FROM, NAMED IN THE MESSAGES: THE STANDARD GIVES THESE KEYWORDS NO RANGE AND NO FORMULA CHECK
WAVELENGTH_RULE = "house rule, DY-2358"
SPEC_SYE_RULE = "house rule, DY-2359"

# SPEC_ERR IS WRITTEN AT FULL PRECISION, SO THIS ONLY ALLOWS FOR A WRITER THAT ROUNDS IT
SPEC_ERR_TOLERANCE = 1e-4
MIN_LINE_COUNT = 1


def _is_number(value: object) -> bool:
    """*return True for a finite int or float, a bool never passing as a number*"""
    return isinstance(value, int | float) and not isinstance(value, bool) and math.isfinite(value)


@register_check
def check_wavelength_keyword_ranges(product: Phase3Product) -> list[Finding]:
    """*report a* ``LAMRMS`` *or* ``SPEC_ERR`` *that is not above 0, or a* ``LAMNLIN`` *below 1*

    A value that is not a number is left to ``check_keyword_value_types``.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_wavelength_keyword_ranges(product)
    ```
    """
    header = product.header(PRIMARY)
    rules = (
        ("LAMRMS", lambda value: value > 0, "greater than 0"),
        ("SPEC_ERR", lambda value: value > 0, "greater than 0"),
        ("LAMNLIN", lambda value: value >= MIN_LINE_COUNT, f"at least {MIN_LINE_COUNT}"),
    )
    return [
        fail(f"{keyword} is {header[keyword]}, outside its allowed range ({WAVELENGTH_RULE}): it must be {text}")
        for keyword, isValid, text in rules
        if keyword in header and _is_number(header[keyword]) and not isValid(header[keyword])
    ]


@register_check
def check_spec_err_is_lamrms_over_the_root_of_lamnlin(product: Phase3Product) -> list[Finding]:
    """*report a* ``SPEC_ERR`` *that is not* ``LAMRMS / sqrt(LAMNLIN)``

    The rule is not applied unless all three keywords are numbers and ``LAMNLIN`` is at least 1: a missing keyword is
    reported by the presence check, and a ``LAMNLIN`` out of range by ``check_wavelength_keyword_ranges``.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_spec_err_is_lamrms_over_the_root_of_lamnlin(product)
    ```
    """
    header = product.header(PRIMARY)
    values = [header.get(keyword) for keyword in ("LAMRMS", "LAMNLIN", "SPEC_ERR")]
    if not all(_is_number(value) for value in values) or values[1] < MIN_LINE_COUNT:
        return []
    lamrms, lamnlin, specErr = values
    expected = lamrms / math.sqrt(lamnlin)
    if math.isclose(specErr, expected, rel_tol=SPEC_ERR_TOLERANCE):
        return []
    return [
        fail(
            f"SPEC_ERR is {specErr:.6g}, it must be LAMRMS / sqrt(LAMNLIN) = {expected:.6g} "
            f"(LAMRMS {lamrms:.6g}, LAMNLIN {lamnlin}) ({WAVELENGTH_RULE})"
        )
    ]


@register_check
def check_spec_sye_is_not_negative(product: Phase3Product) -> list[Finding]:
    """*report a* ``SPEC_SYE`` *below 0*

    Zero is allowed: the standard sets ``SPEC_SYE`` to 0 when the systematic error is insignificant. A value that is not
    a number is left to ``check_keyword_value_types``, and a missing keyword to the presence check.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_spec_sye_is_not_negative(product)
    ```
    """
    value = product.header(PRIMARY).get("SPEC_SYE")
    if not _is_number(value) or value >= 0:
        return []
    return [fail(f"SPEC_SYE is {value}, outside its allowed range ({SPEC_SYE_RULE}): it must be at least 0")]
