"""Contracts for the packaged SOXS NIR multi-pinhole predicted-line table (DY-1989, DY-2084)."""

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
# DY-2084: XeI AND HgI WAVELENGTHS ARE NIST RITZ VALUES, NOT NIST OBSERVED VALUES; NO LINE MOVES MORE THAN 0.13 nm
RITZ_IONS = ["XeI", "HgI"]
MAX_RITZ_SHIFT_NM = 0.13
# (ORDER, ION, NIST OBSERVED nm, NIST RITZ nm). THE ARC MEASURES THE RITZ VALUE FOR THE FIRST SIX
RITZ_LINES = [
    (12, "XeI", 1605.203, 1605.329),
    (13, "XeI", 1465.985, 1466.080),
    (14, "XeI", 1365.649, 1365.706),
    (16, "XeI", 1227.289, 1227.1902),
    (16, "HgI", 1215.9, 1216.008),
    (14, "HgI", 1342.658, 1342.635),
    # THE ARC PREFERS NIST OBSERVED FOR 1719.868 AND 1732.942 (BY 0.3 A), BUT LEFT THERE THEY MAKE THE SET-OFFSET GUARD
    # CLIP 1692.017, AND THE O11 SKY GETS WORSE, SO ALL O11 HgI LINES TAKE RITZ
    (11, "HgI", 1692.017, 1692.027),
    (11, "HgI", 1719.868, 1719.828),
    (11, "HgI", 1732.942, 1732.896),
]
WAVELENGTH_ATOL_NM = 5e-4
# (ORDER, ION, nm) SETS REMOVED IN DY-2084: THE ARC FEATURE IS NOT THE LISTED LINE, WHOSE NIST WAVELENGTH IS RIGHT
DROPPED_SETS = [
    (13, "NeI", 1534.82),
    (13, "Ar", 1534.88),
    (16, "ArII", 1167.79),
    (17, "ArII", 1167.79),
    (16, "ArII", 1172.05),
    (17, "ArII", 1172.05),
    (16, "ArII", 1202.71),
]
DROPPED_SET_ATOL_NM = 0.006


def _read_table(name: str) -> Table:
    return Table.read(CALIBRATIONS_DIR / name)


def _uncurated_rows(name: str) -> Table:
    table = _read_table(name)
    return table[~np.isin(table["order"], CURATED_ORDERS)]


def _stable_sorted(table: Table) -> Table:
    # THE SHIPPED TABLE REPEATS SOME (ORDER, WAVELENGTH, SLIT_INDEX) KEYS, SO POSITIONS BREAK THE TIES
    rows = table.copy()
    rows.sort(["order", "slit_index", "detector_x", "detector_y"], kind="stable")
    return rows


def _ions(table: Table) -> np.ndarray:
    return np.array([i.decode().strip() if isinstance(i, bytes) else str(i).strip() for i in table["ion"]])


def _set_rows(table: Table, order: int, ion: str, wavelength: float) -> np.ndarray:
    return (
        (np.asarray(table["order"]) == order)
        & (_ions(table) == ion)
        & np.isclose(np.asarray(table["wavelength"], dtype=float), wavelength, rtol=0, atol=DROPPED_SET_ATOL_NM)
    )


def _dropped_set_rows(table: Table) -> np.ndarray:
    rows = np.zeros(len(table), dtype=bool)
    for order, ion, wavelength in DROPPED_SETS:
        rows |= _set_rows(table, order, ion, wavelength)
    return rows


def _clean_list_rows_kept() -> Table:
    table = _uncurated_rows(SPH_TABLE)
    return table[~_dropped_set_rows(table)]


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


@pytest.mark.parametrize(("order", "ion", "wavelength"), DROPPED_SETS)
def test_misidentified_set_is_removed(order: int, ion: str, wavelength: float) -> None:
    # ARRANGE
    clean = _uncurated_rows(SPH_TABLE)
    curated = _read_table(MPH_TABLE)

    # ACT
    inClean = _set_rows(clean, order, ion, wavelength).sum()
    inCurated = _set_rows(curated, order, ion, wavelength).sum()

    # ASSERT
    assert inClean == PINHOLES_PER_LINE
    assert inCurated == 0


def test_orders_outside_the_curation_keep_the_clean_list_geometry() -> None:
    # ARRANGE
    columns = [name for name in _read_table(MPH_TABLE).colnames if name != "wavelength"]

    # ACT
    curatedRest = _stable_sorted(_uncurated_rows(MPH_TABLE))
    shippedRest = _stable_sorted(_clean_list_rows_kept())

    # ASSERT
    assert curatedRest.colnames == shippedRest.colnames
    assert len(curatedRest) == len(shippedRest)
    for column in columns:
        assert np.array_equal(np.asarray(curatedRest[column]), np.asarray(shippedRest[column])), column


def test_only_ritz_ion_wavelengths_differ_from_the_clean_list() -> None:
    # ARRANGE
    curated = _stable_sorted(_uncurated_rows(MPH_TABLE))
    shipped = _stable_sorted(_clean_list_rows_kept())

    # ACT
    shift = np.asarray(curated["wavelength"], dtype=float) - np.asarray(shipped["wavelength"], dtype=float)
    moved = shift != 0

    # ASSERT
    assert moved.any()
    assert set(_ions(curated[moved])) <= set(RITZ_IONS)
    assert np.max(np.abs(shift)) <= MAX_RITZ_SHIFT_NM


@pytest.mark.parametrize(("order", "ion", "observed", "ritz"), RITZ_LINES)
def test_line_carries_its_nist_ritz_wavelength(order: int, ion: str, observed: float, ritz: float) -> None:
    # ARRANGE
    table = _read_table(MPH_TABLE)
    rows = table[(table["order"] == order) & (_ions(table) == ion)]

    # ACT
    wavelengths = np.unique(np.asarray(rows["wavelength"], dtype=float))

    # ASSERT
    assert np.any(np.isclose(wavelengths, ritz, rtol=0, atol=WAVELENGTH_ATOL_NM))
    assert not np.any(np.isclose(wavelengths, observed, rtol=0, atol=WAVELENGTH_ATOL_NM))
