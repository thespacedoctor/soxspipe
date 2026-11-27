"""Contracts for the ``SKY SHIFT O<n> ERR`` QCs: the uncertainty (nm) of each VIS order's applied sky shift (DY-2359).

The formula follows the classification of the shift that DY-2461 records (``measured``, ``model`` or ``none``), and
the design of DY-2307 gives each case its own formula. NIR gets no error QC here (DY-2360).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.phase3 import fit_card_comment
from soxspipe.commonutils.sky_shift_error import MEDIAN_ERROR_FACTOR
from soxspipe.commonutils.ug_wavelength_transfer import combine_r_sigma_px
from tests.unit._horne_nir_scene import build_extractor as build_nir_extractor
from tests.unit._horne_nir_scene import build_orders as build_nir_orders
from tests.unit._horne_vis_scene import (
    KEYWORDS,
    RI_DISPERSION,
    UG_DISPERSION,
    build_extractor,
    build_sky_scene,
    night_header,
    qc_by_name,
    vis_parameters,
)

pytestmark = pytest.mark.unit

ERROR_CARD = re.compile(r"^SKY SHIFT O(\d+) ERR$")
SHIFT_CARD = re.compile(r"^SKY SHIFT O(\d+)$")
FAINT_HEIGHT = 2.0
# THE VALUES OF THE SHIPPED VIS DETECTOR PARAMETERS (nm AND px)
UNTUNED_NM = 0.0754
ATTACHED_ARC_NM = 0.0051
R_FROM_I_RMS_PX = 0.364
U_RMS_PX = 0.96
G_RMS_PX = 1.14
ERROR_DECIMALS = 4


def _errors(extractor: horne_extraction) -> dict[int, float]:
    if "qc_name" not in extractor.qc:
        return {}
    return {
        int(ERROR_CARD.match(name).group(1)): value
        for name, value in zip(extractor.qc["qc_name"], extractor.qc["qc_value"], strict=True)
        if ERROR_CARD.match(name)
    }


def _tune(extractor: horne_extraction, scene: pd.DataFrame) -> pd.DataFrame:
    return extractor.tune_wavelength_calibration_to_skylines(scene.copy(), "VIS")


def _measured_stats(log: object, tmp_path: Path, scene: pd.DataFrame) -> dict:
    """The sky-line statistics the tuning measures for each order, from an extractor of its own."""
    (tmp_path / "stats").mkdir()
    extractor = build_extractor(log, tmp_path / "stats", night_header())
    orders = np.sort(scene["order"].unique())
    dispersions = extractor._order_dispersions(scene, orders, True)
    return extractor._measure_orders(scene.copy(), orders, "VIS", True, dispersions)[1]


# MEASURED ORDERS


def test_vis_i_measured_has_the_standard_error_of_the_median_of_its_kept_line_offsets(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = build_sky_scene()
    stats = _measured_stats(log, tmp_path, scene)[4]

    # ACT
    _tune(extractor, scene)

    # ASSERT
    expected = MEDIAN_ERROR_FACTOR * stats["scatterNm"] / np.sqrt(stats["nLines"])
    assert stats["nLines"] >= 2
    assert _errors(extractor)[4] == pytest.approx(expected, abs=10**-ERROR_DECIMALS)
    assert _errors(extractor)[4] > 0


def test_vis_r_measured_alone_has_the_standard_error_of_the_median_of_its_kept_line_offsets(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE: A FAINT i LEAVES r WITH NOTHING TO BE PREDICTED FROM
    extractor = build_extractor(log, tmp_path, night_header())
    scene = build_sky_scene(iHeight=FAINT_HEIGHT)
    stats = _measured_stats(log, tmp_path, scene)[1]

    # ACT
    _tune(extractor, scene)

    # ASSERT
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "measured"
    expected = MEDIAN_ERROR_FACTOR * stats["scatterNm"] / np.sqrt(stats["nLines"])
    assert _errors(extractor)[1] == pytest.approx(expected, abs=10**-ERROR_DECIMALS)


# MODELLED ORDERS


def test_vis_r_predicted_from_i_has_the_r_from_i_rms_times_the_dispersion(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune(extractor, build_sky_scene(rHeight=FAINT_HEIGHT))

    # ASSERT
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "predicted"
    assert _errors(extractor)[1] == pytest.approx(R_FROM_I_RMS_PX * RI_DISPERSION, abs=10**-ERROR_DECIMALS)


def test_vis_r_combined_with_the_prediction_has_the_inverse_variance_sigma_of_the_combination(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = build_sky_scene()
    stats = _measured_stats(log, tmp_path, scene)[1]

    # ACT
    _tune(extractor, scene)

    # ASSERT: THE MEASURED SIGMA IS THE ONE THE COMBINATION WEIGHTS THE SHIFT WITH
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "combined"
    sigmaPx = stats["scatterNm"] / RI_DISPERSION / np.sqrt(stats["nLines"])
    expected = combine_r_sigma_px(sigmaPx, R_FROM_I_RMS_PX) * RI_DISPERSION
    assert _errors(extractor)[1] == pytest.approx(expected, abs=10**-ERROR_DECIMALS)
    assert _errors(extractor)[1] < R_FROM_I_RMS_PX * RI_DISPERSION


def test_vis_u_and_g_from_the_transfer_model_have_their_rms_times_their_dispersion(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune(extractor, build_sky_scene())

    # ASSERT
    assert qc_by_name(extractor, "UG SHIFT SOURCE")["qc_value"] == "model"
    assert _errors(extractor)[2] == pytest.approx(U_RMS_PX * abs(UG_DISPERSION), abs=10**-ERROR_DECIMALS)
    assert _errors(extractor)[3] == pytest.approx(G_RMS_PX * abs(UG_DISPERSION), abs=10**-ERROR_DECIMALS)


# ORDERS WITH NO SHIFT


def test_vis_i_refused_has_the_untuned_value(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune(extractor, build_sky_scene(iHeight=FAINT_HEIGHT))

    # ASSERT
    assert qc_by_name(extractor, "SKY SHIFT O4 SRC")["qc_value"] == "none"
    assert _errors(extractor)[4] == UNTUNED_NM


def test_vis_r_refused_has_the_untuned_value(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune(extractor, build_sky_scene(rHeight=FAINT_HEIGHT, iHeight=FAINT_HEIGHT))

    # ASSERT
    assert qc_by_name(extractor, "SKY SHIFT O1 SRC")["qc_value"] == "none"
    assert _errors(extractor)[1] == UNTUNED_NM


def test_vis_u_and_g_left_unshifted_have_the_untuned_value(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE PLACEHOLDER ROTATOR ANGLE MAKES THE TRANSFER MODEL SKIP u AND g
    placeholder = {KEYWORDS["ABSROT"]: 0.0, KEYWORDS["PAR_ANG_START"]: 0.0}
    extractor = build_extractor(log, tmp_path, placeholder)

    # ACT
    _tune(extractor, build_sky_scene())

    # ASSERT
    assert qc_by_name(extractor, "SKY SHIFT O2 SRC")["qc_value"] == "none"
    assert _errors(extractor)[2] == UNTUNED_NM
    assert _errors(extractor)[3] == UNTUNED_NM


def test_vis_u_and_g_of_an_arm_without_the_transfer_model_have_the_untuned_value(log: object, tmp_path: Path) -> None:
    # ARRANGE
    parameters = {key: value for key, value in vis_parameters().items() if key != "ug-wavelength-transfer"}
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams=parameters)

    # ACT
    _tune(extractor, build_sky_scene())

    # ASSERT
    assert _errors(extractor)[2] == UNTUNED_NM
    assert _errors(extractor)[3] == UNTUNED_NM


# SKIPPED TUNING


def test_every_order_has_the_untuned_value_when_tuning_is_switched_off(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    extractor.recipeSettings = {"tune-wavelength-to-skylines": False}

    # ACT
    extractor._tune_or_skip(build_sky_scene(), "VIS")

    # ASSERT
    assert _errors(extractor) == dict.fromkeys((1, 2, 3, 4), UNTUNED_NM)
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skipped: switched off"


def test_every_order_has_the_attached_arc_value_when_an_arc_is_attached(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    extractor._arc_attachment = lambda: (1.0, 1.0)

    # ACT
    extractor._tune_or_skip(build_sky_scene(), "VIS")

    # ASSERT
    assert _errors(extractor) == dict.fromkeys((1, 2, 3, 4), ATTACHED_ARC_NM)
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skipped: attached arc"


def test_a_skipped_vis_frame_has_error_cards_but_no_shift_cards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    extractor._arc_attachment = lambda: (1.0, 1.0)

    # ACT
    extractor._tune_or_skip(build_sky_scene(), "VIS")

    # ASSERT
    assert not [name for name in extractor.qc["qc_name"] if SHIFT_CARD.match(name)]


# THE CARDS


def test_every_vis_order_with_a_shift_card_has_an_error_card(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune(extractor, build_sky_scene(rHeight=FAINT_HEIGHT))

    # ASSERT
    shiftOrders = {int(SHIFT_CARD.match(n).group(1)) for n in extractor.qc["qc_name"] if SHIFT_CARD.match(n)}
    assert shiftOrders == {1, 2, 3, 4}
    assert set(_errors(extractor)) == shiftOrders


def test_an_error_card_is_a_per_order_header_qc_in_nm_rounded_to_four_decimals(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune(extractor, build_sky_scene())

    # ASSERT
    row = qc_by_name(extractor, "SKY SHIFT O3 ERR")
    assert row["qc_unit"] == "nm"
    assert row["to_header"]
    assert row["qc_order"] == 3
    assert row["qc_value"] == round(row["qc_value"], ERROR_DECIMALS)
    assert row["qc_comment"]


@pytest.mark.parametrize(
    "scene",
    [build_sky_scene(), build_sky_scene(rHeight=FAINT_HEIGHT), build_sky_scene(iHeight=FAINT_HEIGHT)],
    ids=["combined", "predicted", "measured-and-none"],
)
def test_every_error_comment_fits_its_fits_card(log: object, tmp_path: Path, scene: pd.DataFrame) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune(extractor, scene)

    # ASSERT
    rows = extractor.qc.loc[extractor.qc["qc_name"].str.match(ERROR_CARD)]
    assert len(rows) == 4
    for name, value, comment in zip(rows["qc_name"], rows["qc_value"], rows["qc_comment"], strict=True):
        assert fit_card_comment(f"ESO QC {name}", value, comment) == comment, name


def test_a_second_extraction_keeps_the_larger_error_of_an_order_in_one_row(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE FIRST EXTRACTION MEASURES i, THE SECOND REFUSES IT, AND BOTH WRITE TO ONE QC TABLE
    extractor = build_extractor(log, tmp_path, night_header())
    _tune(extractor, build_sky_scene())
    measured = _errors(extractor)[4]

    # ACT
    _tune(extractor, build_sky_scene(iHeight=FAINT_HEIGHT))

    # ASSERT
    assert measured < UNTUNED_NM
    assert qc_by_name(extractor, "SKY SHIFT O4 ERR")["qc_value"] == UNTUNED_NM


def test_a_second_extraction_with_a_smaller_error_does_not_lower_it(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    _tune(extractor, build_sky_scene(iHeight=FAINT_HEIGHT))

    # ACT
    _tune(extractor, build_sky_scene())

    # ASSERT
    assert qc_by_name(extractor, "SKY SHIFT O4 ERR")["qc_value"] == UNTUNED_NM


def test_a_vis_arm_without_the_untuned_value_has_no_card_for_the_orders_that_need_it_and_one_warning(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE: i IS REFUSED, SO u AND g ARE LEFT UNSHIFTED TOO, AND THE THREE NEED THE UNTUNED VALUE. r IS MEASURED
    parameters = {key: value for key, value in vis_parameters().items() if key != "sky-shift-uncertainty-nm"}
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams=parameters)

    # ACT
    _tune(extractor, build_sky_scene(iHeight=FAINT_HEIGHT))

    # ASSERT
    assert set(_errors(extractor)) == {1}
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "orders 2, 3, 4" in warnings[0]


def test_a_skipped_vis_frame_of_an_arm_without_the_value_has_no_error_cards_and_one_warning(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    parameters = {key: value for key, value in vis_parameters().items() if key != "sky-shift-uncertainty-nm"}
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams=parameters)
    extractor._arc_attachment = lambda: (1.0, 1.0)

    # ACT
    extractor._tune_or_skip(build_sky_scene(), "VIS")

    # ASSERT
    assert _errors(extractor) == {}
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "attached-arc" in warnings[0]


# NO ERROR CARDS FOR THE OTHER ARMS (DY-2360 ADDS NIR)


def test_nir_global_fit_writes_no_error_cards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_nir_extractor(log, tmp_path)

    # ACT
    extractor._tune_or_skip(build_nir_orders(), "NIR")

    # ASSERT
    assert _errors(extractor) == {}


def test_nir_with_tuning_switched_off_writes_no_error_cards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_nir_extractor(log, tmp_path)
    extractor.recipeSettings = {"tune-wavelength-to-skylines": False}

    # ACT
    extractor._tune_or_skip(build_nir_orders(), "NIR")

    # ASSERT
    assert _errors(extractor) == {}


def test_x_shooter_vis_writes_no_error_cards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header(), inst="XSH")

    # ACT
    _tune(extractor, build_sky_scene())

    # ASSERT
    assert _errors(extractor) == {}
    assert qc_by_name(extractor, "SKY SHIFT O4 SRC")["qc_value"] == "measured"


def test_tuning_all_orders_together_writes_no_error_cards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    extractor.tune_wavelength_calibration_to_skylines(build_sky_scene(), "VIS", byOrder=False)

    # ASSERT
    assert _errors(extractor) == {}
