"""Contracts for the packaged SOXS VIS single- and multi-pinhole predicted-line tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

import soxspipe

pytestmark = pytest.mark.unit

CALIBRATIONS_DIR = Path(soxspipe.__file__).parent / "resources" / "static_calibrations" / "soxs"

# ONE ROW PER (LINE, PINHOLE): 135 LINES x 1 OR 9 PINHOLES
LINE_COUNT = 135
PINHOLES_PER_LINE = {"SPH": 1, "MPH": 9}
TABLES = [(tech, binning) for tech in PINHOLES_PER_LINE for binning in ("1x1", "2x1")]

# DY-1578: ORDER-1 Th LINES BLENDED WITH A STRONGER NEIGHBOUR. EVERY PINHOLE SITS
# +1.1, +2.0 AND -1.5 PX OFF THE SPATIAL SOLUTION, AND THE SINGLE PINHOLE SITS
# +1.8, +2.7 AND -0.8 PX OFF THE DISPERSION SOLUTION, SO NEITHER FIT MAY USE THEM
BLENDED_WAVELENGTHS = [657.765628, 658.461292, 659.346205]
WAVELENGTH_TOLERANCE = 1e-4


def _read_table(tech: str, binning: str) -> Table:
    return Table.read(CALIBRATIONS_DIR / f"SOXS_{tech}_TAB_VIS_ThAr_{binning}.fits")


def _rows_at(table: Table, wavelength: float) -> np.ndarray:
    return np.abs(np.asarray(table["wavelength"]) - wavelength) < WAVELENGTH_TOLERANCE


@pytest.mark.parametrize(("tech", "binning"), TABLES)
@pytest.mark.parametrize("wavelength", BLENDED_WAVELENGTHS)
def test_blended_th_line_is_deleted_in_every_pinhole(tech: str, binning: str, wavelength: float) -> None:
    # ARRANGE
    table = _read_table(tech, binning)

    # ACT
    rows = _rows_at(table, wavelength)

    # ASSERT
    assert rows.sum() == PINHOLES_PER_LINE[tech]
    assert np.all(np.asarray(table["delete"])[rows] == 1)


@pytest.mark.parametrize(("tech", "binning"), TABLES)
def test_only_the_blended_th_lines_are_deleted(tech: str, binning: str) -> None:
    # ARRANGE
    table = _read_table(tech, binning)
    blended = np.zeros(len(table), dtype=bool)
    for wavelength in BLENDED_WAVELENGTHS:
        blended |= _rows_at(table, wavelength)

    # ACT
    deleted = np.asarray(table["delete"]) == 1

    # ASSERT
    assert len(table) == LINE_COUNT * PINHOLES_PER_LINE[tech]
    assert np.array_equal(deleted, blended)
