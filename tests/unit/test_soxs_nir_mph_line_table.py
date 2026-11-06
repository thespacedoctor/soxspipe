"""Contracts for the packaged SOXS NIR multi-pinhole predicted-line table (DY-1989)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.table import Table

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
RESOURCES_DIR = PACKAGE_DIR / "resources"
CALIBRATIONS_DIR = RESOURCES_DIR / "static_calibrations" / "soxs"
MPH_TABLE = "SOXS_MPH_TAB_NIR_ArHgNeXe.fits"
SPH_TABLE = "ArHgNeXe_clean_within_2.0pixel.fits"

# DY-1989: O10 AND O11 ARE REBUILT FROM LINES CHECKED ON THE ARC FRAME
CURATED_ORDERS = [10, 11]
PINHOLES_PER_LINE = 9
# A PINHOLE FURTHER THAN THIS FROM A SMOOTH CURVE ALONG THE SLIT IS A CORRUPT POSITION
MAX_GEOMETRY_RESIDUAL_PX = 0.5
SLIT_GEOMETRY_POLY_DEG = 2
# THE O11 RED HALF HAD NO FULL SET BEFORE DY-1989
O11_RED_HALF_NM = 1800.0
MIN_O11_RED_HALF_LINES = 8
# WITH THE O11 RED HALF NOW CONSTRAINED, THE X FIT NEEDS ONE MORE ORDER AND WAVELENGTH DEGREE
NIR_X_ORDER_DEG = 5
NIR_X_WAVELENGTH_DEG = 6


def _read_table(name: str) -> Table:
    return Table.read(CALIBRATIONS_DIR / name)


def _nir_spatial_solution_settings() -> dict:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return settings["soxs-spat-solution"]["nir"]


def _nir_line_tables() -> dict:
    detectorParams = yaml.safe_load((RESOURCES_DIR / "soxs_detector_parameters.yaml").read_text())
    return detectorParams["NIR"]["predicted pinhole lines"]


def test_nir_multi_pinhole_frames_read_the_curated_table() -> None:
    # ACT
    tables = _nir_line_tables()

    # ASSERT
    assert tables["multi"]["1x1"] == MPH_TABLE
    assert tables["single"]["1x1"] == SPH_TABLE


def test_nir_spatial_solution_x_fit_uses_the_raised_degrees() -> None:
    # ACT
    settings = _nir_spatial_solution_settings()

    # ASSERT
    assert settings["order-deg"][0] == NIR_X_ORDER_DEG
    assert settings["wavelength-deg"][0] == NIR_X_WAVELENGTH_DEG


@pytest.mark.parametrize("order", CURATED_ORDERS)
def test_every_curated_line_is_a_full_pinhole_set(order: int) -> None:
    # ARRANGE
    table = _read_table(MPH_TABLE)
    rows = table[table["order"] == order]

    # ACT
    wavelengths, counts = np.unique(np.asarray(rows["wavelength"]), return_counts=True)

    # ASSERT
    assert len(wavelengths) > 0
    assert np.all(counts == PINHOLES_PER_LINE)
    for wavelength in wavelengths:
        line = rows[rows["wavelength"] == wavelength]
        assert sorted(line["slit_index"]) == list(range(PINHOLES_PER_LINE))


@pytest.mark.parametrize("order", CURATED_ORDERS)
def test_curated_pinhole_positions_follow_a_smooth_slit_geometry(order: int) -> None:
    # ARRANGE
    table = _read_table(MPH_TABLE)
    rows = table[table["order"] == order]

    # ACT
    offenders = []
    for wavelength in np.unique(np.asarray(rows["wavelength"])):
        line = rows[rows["wavelength"] == wavelength]
        slit = np.asarray(line["slit_position"])
        for axis in ("detector_x", "detector_y"):
            position = np.asarray(line[axis])
            residual = position - np.polyval(np.polyfit(slit, position, SLIT_GEOMETRY_POLY_DEG), slit)
            if np.max(np.abs(residual)) >= MAX_GEOMETRY_RESIDUAL_PX:
                offenders.append(f"O{order} {wavelength} {axis}")

    # ASSERT
    assert offenders == []


def test_o11_red_half_has_enough_lines() -> None:
    # ARRANGE
    table = _read_table(MPH_TABLE)
    rows = table[(table["order"] == 11) & (table["wavelength"] > O11_RED_HALF_NM)]

    # ACT
    redLines = np.unique(np.asarray(rows["wavelength"]))

    # ASSERT
    assert len(redLines) >= MIN_O11_RED_HALF_LINES


def test_orders_outside_the_curation_are_unchanged() -> None:
    # ARRANGE
    curated = _read_table(MPH_TABLE)
    shipped = _read_table(SPH_TABLE)

    # ACT
    curatedRest = curated[~np.isin(curated["order"], CURATED_ORDERS)]
    shippedRest = shipped[~np.isin(shipped["order"], CURATED_ORDERS)]
    # THE SHIPPED TABLE REPEATS SOME (ORDER, WAVELENGTH, SLIT_INDEX) KEYS, SO POSITIONS BREAK THE TIES
    sortKeys = ["order", "wavelength", "slit_index", "detector_x", "detector_y"]
    curatedRest.sort(sortKeys)
    shippedRest.sort(sortKeys)

    # ASSERT
    assert curatedRest.colnames == shippedRest.colnames
    assert len(curatedRest) == len(shippedRest)
    for column in curatedRest.colnames:
        assert np.array_equal(np.asarray(curatedRest[column]), np.asarray(shippedRest[column])), column
