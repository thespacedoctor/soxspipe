#!/usr/bin/env python
"""
*The error and flux columns of a merged spectrum, in the SOXS layout and in the X-Shooter layout*

A SOXS merged spectrum holds ``FLUX_COUNTS`` and ``ERR_COUNTS`` (the error in electrons), and, once a response was
applied, ``FLUX`` and ``ERR``. An X-Shooter merged spectrum keeps the older layout: ``FLUX_COUNTS`` and ``VARIANCE``,
and ``FLUX_CALIBRATED`` and ``ERR_CALIBRATED`` once calibrated. The layout follows the columns, so code that reads a
spectrum finds the error and the calibrated flux whichever instrument made it.

Author
: David Young

Date Created
: November 27, 2026
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

ERR_COUNTS_COLUMN = "ERR_COUNTS"
VARIANCE_COLUMN = "VARIANCE"
FLUX_COLUMN = "FLUX"
ERR_COLUMN = "ERR"
# THE X-SHOOTER NAMES OF THE CALIBRATED FLUX AND ERROR
XSH_FLUX_COLUMN = "FLUX_CALIBRATED"
XSH_ERR_COLUMN = "ERR_CALIBRATED"


def error_column_name(columns: Iterable[str]) -> str | None:
    """*return the name of the column that holds the error or the variance of the counts*

    **Key Arguments:**

    - ``columns`` -- the column names of a merged spectrum

    **Return:**

    - ``name`` -- ``ERR_COUNTS`` or ``VARIANCE``, or *None* when the spectrum has neither
    """
    names = list(columns)
    for name in (ERR_COUNTS_COLUMN, VARIANCE_COLUMN):
        if name in names:
            return name
    return None


def error_counts(spectrum: Any) -> Any:
    """*return the 1-sigma error of the counts of a merged spectrum as a float array*

    **Key Arguments:**

    - ``spectrum`` -- a merged spectrum (data frame or table) with an ``ERR_COUNTS`` or ``VARIANCE`` column

    **Return:**

    - ``error`` -- the error at each sample, in electrons, in float64. ``ERR_COUNTS`` is returned as it is, and the
      error of an X-Shooter spectrum is the square root of its ``VARIANCE``
    """
    import numpy as np

    name = error_column_name(spectrum.columns)
    if name is None:
        raise ValueError(f"The spectrum has neither an {ERR_COUNTS_COLUMN} nor a {VARIANCE_COLUMN} column.")
    values = np.asarray(spectrum[name], dtype=np.float64)
    return values if name == ERR_COUNTS_COLUMN else np.sqrt(values)


def variance_counts(spectrum: Any) -> Any:
    """*return the variance of the counts of a merged spectrum as a float array*

    **Key Arguments:**

    - ``spectrum`` -- a merged spectrum (data frame or table) with an ``ERR_COUNTS`` or ``VARIANCE`` column

    **Return:**

    - ``variance`` -- the variance at each sample, in float64: the square of ``ERR_COUNTS``, or ``VARIANCE`` itself
    """
    import numpy as np

    name = error_column_name(spectrum.columns)
    if name == VARIANCE_COLUMN:
        return np.asarray(spectrum[name], dtype=np.float64)
    return error_counts(spectrum) ** 2


def calibrated_flux_column_name(columns: Iterable[str]) -> str:
    """*return the name of the calibrated flux column of a flux-calibrated spectrum*

    **Key Arguments:**

    - ``columns`` -- the column names of a flux-calibrated spectrum

    **Return:**

    - ``name`` -- ``FLUX``, or ``FLUX_CALIBRATED`` for an X-Shooter spectrum
    """
    return FLUX_COLUMN if FLUX_COLUMN in list(columns) else XSH_FLUX_COLUMN
