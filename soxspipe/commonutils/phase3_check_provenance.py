#!/usr/bin/env python
"""
*the checks of the frame-accounting keywords of a Phase 3 spectrum: NCOMBINE, PROVi, OBIDi, MJD-END, EXPTIME*

The recipes build these keywords from their on-target frame list (DY-2356). The checks look for the faults a wrong list
leaves behind: a frame count that is not the number of provenance cards, a gap in the numbered cards, a provenance name
that still has its ``.fits`` suffix, an observation that ends before it starts, and a per-arm ``EXPTIME`` that is not
the ``TEXPTIME``.

A joined product has the median ``EXPTIME`` of its arms and the maximum ``TEXPTIME``, so the two differ by design
(Dave, 2026-11-24). The joined writer drops the ``ESO SEQ`` cards of the arm and a per-arm product keeps them, so the
card ``ESO SEQ ARM`` tells the two apart.

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from soxspipe.commonutils.phase3_check_framework import Finding, Phase3Product, fail, register_check

if TYPE_CHECKING:
    from astropy.io.fits import Header

PROVENANCE_PATTERN = re.compile(r"PROV(\d+)")
OBSERVATION_BLOCK_PATTERN = re.compile(r"OBID(\d+)")

# A PER-ARM PRODUCT KEEPS THIS CARD OF ITS ARM. THE JOINED WRITER DROPS EVERY `ESO SEQ` CARD
PER_ARM_CARD = "ESO SEQ ARM"

# HEADER TIMES ARE ROUNDED, SO EXPTIME AND TEXPTIME MAY DIFFER BY THIS MANY SECONDS WITHOUT BEING A FAULT
EXPOSURE_TOLERANCE_S = 1e-3


def _indices(header: Header, pattern: re.Pattern[str]) -> list[int]:
    """*return the indices of the numbered cards in a header that match a pattern, in ascending order*"""
    return sorted(int(match.group(1)) for key in header if (match := pattern.fullmatch(key)))


def _gaps(prefix: str, indices: list[int]) -> list[Finding]:
    """*return a failure for each numbered card that is missing below the highest one*"""
    return [
        fail(f"{prefix}{number} is missing but {prefix}{indices[-1]} is there: the {prefix}i cards must run from 1")
        for number in range(1, indices[-1] + 1 if indices else 1)
        if number not in indices
    ]


@register_check
def check_ncombine_matches_the_provenance_cards(product: Phase3Product) -> list[Finding]:
    """*report an NCOMBINE that is not the number of PROVi cards*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_ncombine_matches_the_provenance_cards(product)
    ```
    """
    ncombine = product.primary.get("NCOMBINE")
    count = len(_indices(product.primary, PROVENANCE_PATTERN))
    if not isinstance(ncombine, int) or isinstance(ncombine, bool) or ncombine == count:
        return []
    return [fail(f"NCOMBINE is {ncombine}, but the header has {count} PROVi cards: there is one for each frame")]


@register_check
def check_numbered_cards_run_from_one(product: Phase3Product) -> list[Finding]:
    """*report a gap in the PROVi cards or in the OBIDi cards*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_numbered_cards_run_from_one(product)
    ```
    """
    return _gaps("PROV", _indices(product.primary, PROVENANCE_PATTERN)) + _gaps(
        "OBID", _indices(product.primary, OBSERVATION_BLOCK_PATTERN)
    )


@register_check
def check_provenance_names_have_no_fits_suffix(product: Phase3Product) -> list[Finding]:
    """*report a PROVi card whose value still ends in* ``.fits``

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_provenance_names_have_no_fits_suffix(product)
    ```
    """
    return [
        fail(f"PROV{number} is {product.primary[f'PROV{number}']!r}: the name of a raw frame has no .fits suffix")
        for number in _indices(product.primary, PROVENANCE_PATTERN)
        if str(product.primary[f"PROV{number}"]).endswith(".fits")
    ]


@register_check
def check_observation_ends_after_it_starts(product: Phase3Product) -> list[Finding]:
    """*report an MJD-END that is not after MJD-OBS*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_observation_ends_after_it_starts(product)
    ```
    """
    start, end = product.primary.get("MJD-OBS"), product.primary.get("MJD-END")
    if not all(isinstance(value, int | float) and not isinstance(value, bool) for value in (start, end)):
        return []
    if end > start:
        return []
    return [fail(f"MJD-END ({end}) is not after MJD-OBS ({start})")]


@register_check
def check_per_arm_exposure_equals_total_exposure(product: Phase3Product) -> list[Finding]:
    """*report a per-arm product whose EXPTIME is not its TEXPTIME*

    A product with no ``ESO SEQ ARM`` card is joined, and is not checked.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_per_arm_exposure_equals_total_exposure(product)
    ```
    """
    exptime, texptime = product.primary.get("EXPTIME"), product.primary.get("TEXPTIME")
    if PER_ARM_CARD not in product.primary:
        return []
    if not all(isinstance(value, int | float) and not isinstance(value, bool) for value in (exptime, texptime)):
        return []
    if abs(exptime - texptime) <= EXPOSURE_TOLERANCE_S:
        return []
    return [fail(f"EXPTIME ({exptime} s) is not TEXPTIME ({texptime} s): a per-arm product has one total exposure")]
