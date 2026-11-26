"""Contracts for the packaged SOXS NIR single-pinhole predicted-line table (DY-1970, DY-2448, DY-2172)."""

from __future__ import annotations

import numpy as np
import pytest

from tests.unit._nir_line_tables import (
    DROPPED_SET_ATOL_NM,
    DROPPED_SETS,
    EDITED_ORDER,
    NEI_NM,
    PINHOLES_PER_LINE,
    POSITION_ATOL_PX,
    SLIT_INDEX_AT_CENTRE,
    SPH_O10_ROWS_BEFORE,
    SPH_TABLE,
    WAVELENGTH_ATOL_NM,
    ions,
    order_checksums,
    positions_by_slit,
    read_table,
    set_rows,
    untouched_rows_checksum,
)

pytestmark = pytest.mark.unit

# DY-1970: ArII LABELS WITH NO PEAK AT THE LABEL WAVELENGTH, REMOVED FROM ORDER 10 (9 ROWS EACH)
REMOVED_ARII_NM = [1902.3024, 1948.1178, 1963.88837, 1982.6433, 1995.2419, 2000.1502, 2007.8336]
# DY-2448: THE ArII 2002.881 SET SAT 8.9 px FROM THE NEAREST PEAK (x = 1677.7, 2002.97 nm, 0.87 A FROM THE LABEL, NO
# OTHER NIST LINE NEAR IT) AND 3.3 px FROM A BRIGHTER UNIDENTIFIED PEAK (x = 1672.1, 2002.52 nm), SO IT IS REMOVED
ARII_2002_881_NM = 2002.88075
ARII_2002_881_ATOL_NM = 1e-3
# 7 ArII SETS OUT, 1 HgI SET IN (DY-1970), THEN THE ArII 2002.881 SET OUT (DY-2448)
SPH_O10_ROWS_AFTER = SPH_O10_ROWS_BEFORE - (len(REMOVED_ARII_NM) + 1) * PINHOLES_PER_LINE + PINHOLES_PER_LINE
# DY-2172: THE O10 HgI 1970.018 SET TAKES ITS NIST RITZ WAVELENGTH
HGI_NM = 1970.02
# DY-2172: (ORDER, ION, NIST OBSERVED nm, NIST RITZ nm) OF THE O10 XeI AND HgI SETS THAT TAKE THEIR RITZ WAVELENGTH
O10_RITZ_LINES = [
    (EDITED_ORDER, "HgI", 1970.018, HGI_NM),
    (EDITED_ORDER, "XeI", 1878.8154, 1878.814),
    (EDITED_ORDER, "XeI", 1879.166, 1879.02),
]
NEI_SHIFT_X_PX = 3.86
NEI_SHIFT_Y_PX = 0.30
# (SLIT_INDEX, X, Y) OF THE NeI 1957.7118 ROWS ON DEVELOP b9ab9446, BEFORE THE SHIFT
DEVELOP_NEI_POSITIONS = [
    (0, 1112.3971400588853, 379.5850152341245),
    (1, 1127.2656845023248, 380.80856450969635),
    (2, 1125.6673376823676, 384.704018457522),
    (3, 1122.7816284283977, 388.7188885397159),
    (4, 1124.1473745234518, 392.63836472177053),
    (5, 1126.8774173584497, 396.8803647864288),
    (6, 1121.146428561533, 404.67422274087926),
    (7, 1110.2702447855547, 407.7294245109801),
    (8, 1111.511801728674, 422.0960986487003),
]
# (X, Y) OF THE HgI 1970.018 ROWS, BY SLIT INDEX
# SLIT 4 IS THE MEASURED ARC POSITION MINUS THE DSOL CLUSTER SHIFT (1276.42 - 2.549, 406.92 + 0.171)
HGI_POSITIONS = [
    (1278.1230753910518, 391.55776213191405),
    (1276.9768968821256, 395.4179414977179),
    (1275.9227935610847, 399.2766641155975),
    (1274.8664993384682, 403.14906617551685),
    (1273.871, 407.091),
    (1272.7921144078666, 410.90858115049605),
    (1271.71028810319, 414.9969388894763),
    (1270.6418629320804, 418.8197060197414),
    (1269.6192418571436, 422.7063762219037),
]
# THE ORDER 10 LINES THE DY-1970 EDIT ADDED (HgI) OR MOVED (NeI); THE OTHER ORDER 10 ROWS ARE UNTOUCHED, EXCEPT THAT
# DY-2448 ALSO REMOVED THE ArII 2002.881 SET FROM ORDER 10 (SO NO ROW OF IT IS HASHED)
TOUCHED_SETS = [("NeI", NEI_NM), ("HgI", HGI_NM)]
# SHA-256 OF THE UNTOUCHED ORDER 10 ROWS (EVERY ROW BUT THE TOUCHED SETS), IN FILE ORDER
# REGENERATE WITH untouched_rows_checksum(read_table(SPH_TABLE), EDITED_ORDER, TOUCHED_SETS);
# A LATER LEGITIMATE EDIT TO ORDER 10 MUST RE-RECORD IT; RE-RECORDED IN DY-2448 FOR THE REMOVED ArII 2002.881 SET AND IN
# DY-2172 FOR THE RITZ WAVELENGTHS OF THE O10 XeI SETS
UNTOUCHED_O10_CHECKSUM = "82d2057d1a0e82d77b172874601354847a4780ab83c452fc8bbfe29a041305f1"
# SHA-256 OF EVERY ROW OF EACH ORDER OTHER THAN 10 (UNEDITED ON DEVELOP b9ab9446, EXCEPT THAT DY-2172 RE-RECORDED
# ORDERS 11 TO 21 FOR THE RITZ WAVELENGTHS AND THE REMOVED SETS), IN FILE ORDER
# REGENERATE WITH order_checksums(read_table(SPH_TABLE));
# A LATER LEGITIMATE EDIT TO ANY OF THESE ORDERS MUST RE-RECORD IT
UNEDITED_ORDER_CHECKSUMS = {
    9: "53e1d3777799a1df3500cae5d2f7100e3bfc418d0441eedb23e772dc957b6cd9",
    11: "12c276d2b05fedd54c4101b613faca493fe6b1be5cc817df6985e68e6ff28ec7",
    12: "4dd5aff44121b74c4b815a0853e242ab3902ba0c8bb6fa8cc0a4dbabf71f3d26",
    13: "a4f6e772159a3a94803fe11372d6b0befa192b3a7735b93aed8649a516621236",
    14: "a782bf3fdeaa5616255670046809e7f7281f2d42d8a1da891911ce3ff9bb51d3",
    15: "e7dbdc879f994fc39a74eb54fc855b8de00d1381b54739fd71d27901d06cd560",
    16: "4029bd6a8631d22fa1bb59229e4ad268b571ed7db41236598ac2e37a5a4238db",
    17: "231943821a58df96b918e1b755a40f7613b8e9d01c4839ded9262592fc3b3034",
    18: "202a08c945d5617df655308175a23c175607455a40be50046848e1cdc43137e6",
    19: "b33f3297b5631f0c0308b4ddc86067f7a0165691a50b636cd680103e446d8c0b",
    20: "3e5029ba01f47acc7f8ef4402e373748d30752ad35d20387d621fb32356b7b9f",
    21: "d3f5df8078efeecc9326cf8c08b4ff7022a5de32023bd10dbf9b98d9ff619f6a",
    22: "0f04755b8de378a94b42ebdab9eee1d364dd227206c4443cedde7f28fe4dc8cd",
    23: "cfaa014f4f046a660d1a7b17820ce3d044c777a46e4c31506e148d726acc6efa",
    24: "0df973dad267f1c549cb9fda93f6a428ae15e7ace4e07a7092f3ee271b2c448c",
}


@pytest.mark.parametrize("wavelength", REMOVED_ARII_NM)
def test_mislabelled_arii_set_is_absent_from_order_10(wavelength: float) -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    rows = set_rows(table, EDITED_ORDER, "ArII", wavelength)

    # ASSERT
    assert len(rows) == 0


def test_arii_2002_881_is_absent_from_every_order() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    nearby = np.isclose(
        np.asarray(table["wavelength"], dtype=float), ARII_2002_881_NM, rtol=0, atol=ARII_2002_881_ATOL_NM
    )

    # ASSERT
    assert not nearby.any()


def test_order_10_row_count_drops_by_the_removed_sets_less_the_added_set() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    rowCount = int(np.sum(np.asarray(table["order"]) == EDITED_ORDER))

    # ASSERT
    assert rowCount == SPH_O10_ROWS_AFTER


def test_nei_1957_rows_sit_one_line_spacing_from_their_develop_positions() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)
    expected = {s: (x + NEI_SHIFT_X_PX, y + NEI_SHIFT_Y_PX) for s, x, y in DEVELOP_NEI_POSITIONS}

    # ACT
    actual = positions_by_slit(set_rows(table, EDITED_ORDER, "NeI", NEI_NM))

    # ASSERT
    assert actual.keys() == expected.keys()
    for slit, (x, y) in expected.items():
        assert actual[slit] == pytest.approx((x, y), rel=0, abs=POSITION_ATOL_PX)


def test_hgi_1970_has_a_full_set_with_the_validated_slit_4_position() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    rows = set_rows(table, EDITED_ORDER, "HgI", HGI_NM)
    positions = positions_by_slit(rows)

    # ASSERT
    assert sorted(rows["slit_index"]) == list(range(PINHOLES_PER_LINE))
    assert positions[SLIT_INDEX_AT_CENTRE] == pytest.approx(
        HGI_POSITIONS[SLIT_INDEX_AT_CENTRE], rel=0, abs=POSITION_ATOL_PX
    )


def test_hgi_1970_pinholes_sit_at_their_validated_positions() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    actual = positions_by_slit(set_rows(table, EDITED_ORDER, "HgI", HGI_NM))

    # ASSERT
    assert actual.keys() == set(range(PINHOLES_PER_LINE))
    for slit, position in enumerate(HGI_POSITIONS):
        assert actual[slit] == pytest.approx(position, rel=0, abs=POSITION_ATOL_PX)


def test_order_10_rows_the_edit_did_not_touch_match_the_recorded_checksum() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    checksum = untouched_rows_checksum(table, EDITED_ORDER, TOUCHED_SETS)

    # ASSERT
    assert checksum == UNTOUCHED_O10_CHECKSUM


@pytest.mark.parametrize("order", sorted(UNEDITED_ORDER_CHECKSUMS))
def test_orders_other_than_10_match_develop_row_for_row(order: int) -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    checksums = order_checksums(table)

    # ASSERT
    assert checksums[order] == UNEDITED_ORDER_CHECKSUMS[order]


@pytest.mark.parametrize(("order", "ion", "observed", "ritz"), O10_RITZ_LINES)
def test_order_10_line_carries_its_nist_ritz_wavelength(order: int, ion: str, observed: float, ritz: float) -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)
    rows = table[(np.asarray(table["order"]) == order) & (ions(table) == ion)]

    # ACT
    wavelengths = np.unique(np.asarray(rows["wavelength"], dtype=float))

    # ASSERT
    assert np.any(np.isclose(wavelengths, ritz, rtol=0, atol=WAVELENGTH_ATOL_NM / 10))
    assert not np.any(np.isclose(wavelengths, observed, rtol=0, atol=WAVELENGTH_ATOL_NM / 10))


@pytest.mark.parametrize(("order", "ion", "wavelength"), DROPPED_SETS)
def test_misidentified_set_is_absent(order: int, ion: str, wavelength: float) -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)
    nearby = (
        (np.asarray(table["order"]) == order)
        & (ions(table) == ion)
        & np.isclose(np.asarray(table["wavelength"], dtype=float), wavelength, rtol=0, atol=DROPPED_SET_ATOL_NM)
    )

    # ACT
    rowCount = int(nearby.sum())

    # ASSERT
    assert rowCount == 0
