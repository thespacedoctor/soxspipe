"""Contracts for the packaged SOXS NIR single-pinhole predicted-line table, order 10 edit (DY-1970)."""

from __future__ import annotations

import numpy as np
import pytest

from tests.unit._nir_line_tables import (
    EDITED_ORDER,
    HGI_NM,
    NEI_NM,
    PINHOLES_PER_LINE,
    SPH_O10_ROWS_BEFORE,
    SPH_TABLE,
    order_checksums,
    read_table,
    set_rows,
)

pytestmark = pytest.mark.unit

# DY-1970: ArII LABELS WITH NO PEAK AT THE LABEL WAVELENGTH, REMOVED FROM ORDER 10 (9 ROWS EACH)
REMOVED_ARII_NM = [1902.3024, 1948.1178, 1963.88837, 1982.6433, 1995.2419, 2000.1502, 2007.8336]
KEPT_ARII_NM = 2002.88075
# 7 ArII SETS OUT, 1 HgI SET IN
SPH_O10_ROWS_AFTER = SPH_O10_ROWS_BEFORE - len(REMOVED_ARII_NM) * PINHOLES_PER_LINE + PINHOLES_PER_LINE
NEI_SHIFT_X_PX = 3.86
NEI_SHIFT_Y_PX = 0.30
POSITION_ATOL_PX = 1e-6
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
# HgI 1970.018 AT SLIT 4: THE MEASURED ARC POSITION MINUS THE DSOL CLUSTER SHIFT (1276.42 - 2.549, 406.92 + 0.171)
HGI_SLIT_4_X = 1273.871
HGI_SLIT_4_Y = 407.091
SLIT_INDEX_AT_CENTRE = 4
# SHA-256 OF EVERY ROW OF EACH UNEDITED ORDER ON DEVELOP b9ab9446, IN FILE ORDER
UNEDITED_ORDER_CHECKSUMS = {
    9: "53e1d3777799a1df3500cae5d2f7100e3bfc418d0441eedb23e772dc957b6cd9",
    11: "9c9e370b2c631fb84db83815cce1c0b7ec434a87204fdea812819261bf0cd611",
    12: "128ca8cb44ace676395a00893a57bb1b8f15e05bb08c9ddf5836f3084d01ccc3",
    13: "e148cf5739bcd5fbcb08864744a0f9f341524a292dcd3168698e975b1f3096ab",
    14: "48efdda065475295b40aaae70add22cf61a7f49c86b3b961efa4b64fe26edb89",
    15: "a78a750f3b200f23c733f136209aeacfa44bc5d20ee1331aa339bf0454f48da7",
    16: "1c071f288219bef4cc614ba84bedf74dc89e4373329583ed9e3031f5e2e1d61c",
    17: "7e9b9679671449104981243aa13314d66d44a5b7677e86c348aafe6b724b9fc2",
    18: "9c97831606d2a49e3569d89bf6c612767afb8d531067b5c4f58f24248c302f1b",
    19: "442b9b0cc664d3217be3826674a9ffcd6fc45d77cb7a80357b4beb148600ca4e",
    20: "608fd9d4171e49d446ca81f073315e68c94508b0019ec61524baa49402140474",
    21: "7ea08772bcd49407597b6ddf878ea33e815500484bfc716406a28e23d11920c7",
    22: "0f04755b8de378a94b42ebdab9eee1d364dd227206c4443cedde7f28fe4dc8cd",
    23: "cfaa014f4f046a660d1a7b17820ce3d044c777a46e4c31506e148d726acc6efa",
    24: "0df973dad267f1c549cb9fda93f6a428ae15e7ace4e07a7092f3ee271b2c448c",
}


def _positions_by_slit(rows) -> dict[int, tuple[float, float]]:
    columns = zip(rows["slit_index"], rows["detector_x"], rows["detector_y"], strict=True)
    return {int(s): (float(x), float(y)) for s, x, y in columns}


@pytest.mark.parametrize("wavelength", REMOVED_ARII_NM)
def test_mislabelled_arii_set_is_absent_from_order_10(wavelength: float) -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    rows = set_rows(table, "ArII", wavelength)

    # ASSERT
    assert len(rows) == 0


def test_arii_2002_881_keeps_its_full_set_in_order_10() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    rows = set_rows(table, "ArII", KEPT_ARII_NM)

    # ASSERT
    assert sorted(rows["slit_index"]) == list(range(PINHOLES_PER_LINE))


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
    actual = _positions_by_slit(set_rows(table, "NeI", NEI_NM))

    # ASSERT
    assert actual.keys() == expected.keys()
    for slit, (x, y) in expected.items():
        assert actual[slit] == pytest.approx((x, y), rel=0, abs=POSITION_ATOL_PX)


def test_hgi_1970_has_a_full_set_with_the_validated_slit_4_position() -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    rows = set_rows(table, "HgI", HGI_NM)
    positions = _positions_by_slit(rows)

    # ASSERT
    assert sorted(rows["slit_index"]) == list(range(PINHOLES_PER_LINE))
    assert positions[SLIT_INDEX_AT_CENTRE] == pytest.approx((HGI_SLIT_4_X, HGI_SLIT_4_Y), rel=0, abs=POSITION_ATOL_PX)


@pytest.mark.parametrize("order", sorted(UNEDITED_ORDER_CHECKSUMS))
def test_orders_other_than_10_match_develop_row_for_row(order: int) -> None:
    # ARRANGE
    table = read_table(SPH_TABLE)

    # ACT
    checksums = order_checksums(table)

    # ASSERT
    assert checksums[order] == UNEDITED_ORDER_CHECKSUMS[order]
