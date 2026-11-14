"""Contracts for the packaged SOXS NIR multi-pinhole predicted-line table (DY-1989, DY-2084, DY-1970)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.table import Table

import soxspipe
from tests.unit._nir_line_tables import (
    EDITED_ORDER,
    HGI_NM,
    MPH_O10_ROWS_BEFORE,
    MPH_TABLE,
    NEI_NM,
    PINHOLES_PER_LINE,
    POSITION_ATOL_PX,
    SLIT_INDEX_AT_CENTRE,
    SPH_TABLE,
    WAVELENGTH_ATOL_NM,
    XEI_NM,
    ions,
    order_checksums,
    positions_by_slit,
    read_table,
    set_mask,
    set_rows,
    untouched_rows_checksum,
)

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
RESOURCES_DIR = PACKAGE_DIR / "resources"

# DY-1989: O10 AND O11 ARE REBUILT FROM LINES CHECKED ON THE ARC FRAME
CURATED_ORDERS = [10, 11]
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

# DY-1970: THREE O10 ArII SETS ARE OTHER LINES, 7 TO 9 A BLUER; THE TABLE SWAPS THEM FOR NeI, HgI AND XeI SETS THAT
# CARRY THE NIST WAVELENGTH OF THE LINE SEEN ON THE ARC
O10_REMOVED_ARII_NM = [1948.118, 1982.643, 1995.242]
# (ION, nm, [(X, Y) BY SLIT INDEX 0 TO 8]) OF THE O10 SETS ADDED IN DY-1970
O10_ADDED_SETS = [
    (
        "XeI",
        XEI_NM,
        [
            (238.94428400821, 384.2126726446532),
            (238.0291672711235, 388.3097409456561),
            (237.01387119679896, 392.4341127104791),
            (236.24280510734806, 396.55301744552594),
            (235.24139291562463, 400.64647062444766),
            (234.42978633556734, 404.72281130535845),
            (233.54475099718334, 408.9019194476989),
            (232.7501887397893, 413.02986753969475),
            (231.87786571163622, 417.10651614419385),
        ],
    ),
    (
        "NeI",
        NEI_NM,
        [
            (1133.1216187967545, 377.3162562167707),
            (1131.9770113938387, 381.3074396418454),
            (1130.9558039850574, 385.1798617035703),
            (1129.9770475807022, 389.10916048943017),
            (1128.92364794532, 393.0473836343341),
            (1127.9042316532855, 396.9720631194816),
            (1126.7902833534977, 400.9025266902693),
            (1125.8810051541302, 404.8512953425118),
            (1124.8812355416576, 408.7863302559909),
        ],
    ),
    (
        "HgI",
        HGI_NM,
        [
            (1278.8037203547938, 391.1359111152702),
            (1277.6575418458676, 394.996090481074),
            (1276.6034385248267, 398.85481309895357),
            (1275.5471443022102, 402.7272151588729),
            (1274.551644963742, 406.6691489833562),
            (1273.4727593716086, 410.4867301338521),
            (1272.3909330669321, 414.5750878728325),
            (1271.3225078958224, 418.3978550030975),
            (1270.2998868208856, 422.28452520525985),
        ],
    ),
]
# THE ORDER 10 ROWS THAT ARE NOT IN AN ADDED SET WERE NOT TOUCHED
TOUCHED_SETS = [(ion, wavelength) for ion, wavelength, _ in O10_ADDED_SETS]
# SHA-256 OF THE UNTOUCHED ORDER 10 ROWS (EVERY ROW BUT THE ADDED SETS), IN FILE ORDER
# REGENERATE WITH untouched_rows_checksum(read_table(MPH_TABLE), EDITED_ORDER, TOUCHED_SETS);
# A LATER LEGITIMATE EDIT TO ORDER 10 MUST RE-RECORD IT
UNTOUCHED_O10_CHECKSUM = "ebc444705e14876023920060464918eb387c7c00843186f27700d43f832f0cf2"
# SHA-256 OF EVERY ROW OF EACH UNEDITED ORDER ON DEVELOP b9ab9446, IN FILE ORDER
# REGENERATE WITH order_checksums(read_table(MPH_TABLE));
# A LATER LEGITIMATE EDIT TO ANY OF THESE ORDERS MUST RE-RECORD IT
UNEDITED_ORDER_CHECKSUMS = {
    9: "5145cb61322757b7173fc23a2c9cec44f5f74e53cf6ee617926841fc4fd82f09",
    11: "aa56f3527760a2d88299240d517f3a98de10bfeec9fa9c7a431f36dfab1f1fbd",
    12: "3bf55257e87b67c53720590575c1cd341ac9cead0177f4965ddbfb662420695e",
    13: "9e8f07f8b388191084b55750d6c17ab8f1e91cb6044387a98a5f951701902ed0",
    14: "7d23b5adfc85c58fce36586d0fafb71c7949f4c065f9639681158a6dff2f3b3d",
    15: "3aaef09975e18769be823c5020185c07e8112b81505ecc8cb87e4c702a7f0108",
    16: "3bef626deaf1cd600aa0971b59f31692ec583f2ae165895f98626d562fa47ab3",
    17: "d3da003cd0a07c443fd38578f18e39650f7823d3bcb9ad18ae669a88f8de2776",
    18: "2bdd32cf4e8eabf96e3b61f89906a978d70a18fc39b0960be5efcda5930065cd",
    19: "8f372c54790a7367bcd5eb5e4510e430f2e25e416b9334bf3faa786a7161e707",
    20: "f6e1368155a6253460751c97cf4ab2ccd10987fac658c0612949243065ab77d6",
    21: "a693398eea40ba6d151d8e50373265b29a4d55f349cf9ff370df55bcbb22d8c7",
    22: "a650f53acabff5cb84faf0b4fdf4d13961a8e6ef59d9c8bf8f67ee01bd20c714",
    23: "5e8089ed3e6af14587778d190e51830d4996ec647f02eb914b3ad08e798d36de",
    24: "620216674bf04fdb6053ac473a4ec5eab3c2791470e8a16562b4a3ed30e26c82",
}


def _uncurated_rows(name: str) -> Table:
    table = read_table(name)
    return table[~np.isin(table["order"], CURATED_ORDERS)]


def _stable_sorted(table: Table) -> Table:
    # THE SHIPPED TABLE REPEATS SOME (ORDER, WAVELENGTH, SLIT_INDEX) KEYS, SO POSITIONS BREAK THE TIES
    rows = table.copy()
    rows.sort(["order", "slit_index", "detector_x", "detector_y"], kind="stable")
    return rows


def _dropped_set_rows(table: Table) -> np.ndarray:
    rows = np.zeros(len(table), dtype=bool)
    for order, ion, wavelength in DROPPED_SETS:
        rows |= set_mask(table, order, ion, wavelength, DROPPED_SET_ATOL_NM)
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
    table = read_table(MPH_TABLE)
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
    table = read_table(MPH_TABLE)
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
    table = read_table(MPH_TABLE)
    rows = table[(table["order"] == 11) & (table["wavelength"] > O11_RED_HALF_NM)]

    # ACT
    redLines = np.unique(np.asarray(rows["wavelength"]))

    # ASSERT
    assert len(redLines) >= MIN_O11_RED_HALF_LINES


@pytest.mark.parametrize(("order", "ion", "wavelength"), DROPPED_SETS)
def test_misidentified_set_is_removed(order: int, ion: str, wavelength: float) -> None:
    # ARRANGE
    clean = _uncurated_rows(SPH_TABLE)
    curated = read_table(MPH_TABLE)

    # ACT
    inClean = set_mask(clean, order, ion, wavelength, DROPPED_SET_ATOL_NM).sum()
    inCurated = set_mask(curated, order, ion, wavelength, DROPPED_SET_ATOL_NM).sum()

    # ASSERT
    assert inClean == PINHOLES_PER_LINE
    assert inCurated == 0


def test_orders_outside_the_curation_keep_the_clean_list_geometry() -> None:
    # ARRANGE
    columns = [name for name in read_table(MPH_TABLE).colnames if name != "wavelength"]

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
    assert set(ions(curated[moved])) <= set(RITZ_IONS)
    assert np.max(np.abs(shift)) <= MAX_RITZ_SHIFT_NM


@pytest.mark.parametrize(("order", "ion", "observed", "ritz"), RITZ_LINES)
def test_line_carries_its_nist_ritz_wavelength(order: int, ion: str, observed: float, ritz: float) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    rows = table[(table["order"] == order) & (ions(table) == ion)]

    # ACT
    wavelengths = np.unique(np.asarray(rows["wavelength"], dtype=float))

    # ASSERT
    assert np.any(np.isclose(wavelengths, ritz, rtol=0, atol=WAVELENGTH_ATOL_NM))
    assert not np.any(np.isclose(wavelengths, observed, rtol=0, atol=WAVELENGTH_ATOL_NM))


@pytest.mark.parametrize("wavelength", O10_REMOVED_ARII_NM)
def test_mislabelled_arii_set_is_absent_from_order_10(wavelength: float) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    rows = set_rows(table, EDITED_ORDER, "ArII", wavelength)

    # ASSERT
    assert len(rows) == 0


@pytest.mark.parametrize(("ion", "wavelength", "positions"), O10_ADDED_SETS)
def test_added_order_10_line_is_a_full_set_with_its_measured_centre(
    ion: str, wavelength: float, positions: list[tuple[float, float]]
) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    centre = positions[SLIT_INDEX_AT_CENTRE]

    # ACT
    rows = set_rows(table, EDITED_ORDER, ion, wavelength)
    middle = rows[rows["slit_index"] == SLIT_INDEX_AT_CENTRE]

    # ASSERT
    assert sorted(rows["slit_index"]) == list(range(PINHOLES_PER_LINE))
    assert len(middle) == 1
    assert (float(middle["detector_x"][0]), float(middle["detector_y"][0])) == pytest.approx(
        centre, rel=0, abs=POSITION_ATOL_PX
    )


@pytest.mark.parametrize(("ion", "wavelength", "positions"), O10_ADDED_SETS)
def test_added_order_10_line_has_every_pinhole_at_its_recorded_position(
    ion: str, wavelength: float, positions: list[tuple[float, float]]
) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    actual = positions_by_slit(set_rows(table, EDITED_ORDER, ion, wavelength))

    # ASSERT
    assert actual.keys() == set(range(PINHOLES_PER_LINE))
    for slit, position in enumerate(positions):
        assert actual[slit] == pytest.approx(position, rel=0, abs=POSITION_ATOL_PX)


def test_order_10_row_count_is_unchanged_with_three_sets_out_and_three_in() -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    rowCount = int(np.sum(np.asarray(table["order"]) == EDITED_ORDER))

    # ASSERT
    assert rowCount == MPH_O10_ROWS_BEFORE


def test_order_10_rows_the_edit_did_not_touch_match_the_recorded_checksum() -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    checksum = untouched_rows_checksum(table, EDITED_ORDER, TOUCHED_SETS)

    # ASSERT
    assert checksum == UNTOUCHED_O10_CHECKSUM


@pytest.mark.parametrize("order", sorted(UNEDITED_ORDER_CHECKSUMS))
def test_orders_other_than_10_match_develop_row_for_row(order: int) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    checksums = order_checksums(table)

    # ASSERT
    assert checksums[order] == UNEDITED_ORDER_CHECKSUMS[order]
