"""The pure parts of the sky-shift uncertainty of an order (DY-2359).

The uncertainty of the shift applied to a VIS order depends on where the shift came from. These tests cover the
formulas and the detector-parameter reads behind them. The per-order QCs that use them are tested with the extraction.
"""

from __future__ import annotations

import math

import pytest

from soxspipe.commonutils.sky_shift_error import (
    ATTACHED_ARC_KEY,
    UNTUNED_KEY,
    median_standard_error_nm,
    uncertainty_of_untuned_order_nm,
)
from soxspipe.commonutils.ug_wavelength_transfer import (
    MIN_MEASURED_SIGMA_PX,
    combine_r_shift_px,
    combine_r_sigma_px,
    load_ug_rms_px,
)
from tests.unit._horne_vis_scene import vis_parameters

pytestmark = pytest.mark.unit


def test_the_standard_error_of_the_median_is_the_scatter_times_1_2533_over_root_n() -> None:
    # ACT
    error = median_standard_error_nm(scatterNm=0.01, nLines=25)

    # ASSERT
    assert error == pytest.approx(1.2533 * 0.01 / 5, rel=1e-4)


def test_the_combined_sigma_is_the_inverse_variance_sigma_of_the_measurement_and_the_prediction() -> None:
    # ACT
    sigma = combine_r_sigma_px(measuredSigmaPx=0.2, predictionSigmaPx=0.364)

    # ASSERT
    assert sigma == pytest.approx(1 / math.sqrt(1 / 0.2**2 + 1 / 0.364**2))
    assert sigma < 0.2


def test_the_combined_sigma_uses_the_floored_measurement_sigma_that_weights_the_combined_shift() -> None:
    # ARRANGE: THE COMBINED SHIFT FLOORS A MEASUREMENT SIGMA BELOW MIN_MEASURED_SIGMA_PX, AND THE SIGMA MUST AGREE
    tiny = MIN_MEASURED_SIGMA_PX / 100

    # ACT
    sigma = combine_r_sigma_px(measuredSigmaPx=tiny, predictionSigmaPx=0.364)

    # ASSERT
    assert sigma == pytest.approx(1 / math.sqrt(1 / MIN_MEASURED_SIGMA_PX**2 + 1 / 0.364**2))
    weightMeasured = 1 / MIN_MEASURED_SIGMA_PX**2
    weightPredicted = 1 / 0.364**2
    assert sigma**2 == pytest.approx(1 / (weightMeasured + weightPredicted))
    assert combine_r_shift_px(1.0, tiny, 0.0, 0.364) == pytest.approx(
        weightMeasured / (weightMeasured + weightPredicted)
    )


def test_the_u_and_g_rms_come_from_the_leave_one_night_out_values_of_the_design() -> None:
    # ACT
    rms = load_ug_rms_px(vis_parameters())

    # ASSERT
    assert rms == {2: 0.96, 3: 1.14}


def test_an_arm_without_a_transfer_block_has_no_u_and_g_rms() -> None:
    # ACT
    rms = load_ug_rms_px({})

    # ASSERT
    assert rms is None


@pytest.mark.parametrize("bad", [None, "fast", 0.0, -1.0, float("nan"), float("inf")])
def test_an_unusable_u_and_g_rms_raises(bad: object) -> None:
    # ARRANGE
    orders = {
        2: {"m": 1.0, "cos": 1.0, "sin": 1.0, "rms_px": bad},
        3: {"m": 1.0, "cos": 1.0, "sin": 1.0, "rms_px": 1.0},
    }

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="rms_px"):
        load_ug_rms_px({"ug-wavelength-transfer": {"orders": orders}})


def test_a_missing_u_and_g_rms_raises() -> None:
    # ARRANGE
    orders = {2: {"m": 1.0, "cos": 1.0, "sin": 1.0}, 3: {"m": 1.0, "cos": 1.0, "sin": 1.0, "rms_px": 1.0}}

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="rms_px"):
        load_ug_rms_px({"ug-wavelength-transfer": {"orders": orders}})


@pytest.mark.parametrize(("kind", "expected"), [(UNTUNED_KEY, 0.0754), (ATTACHED_ARC_KEY, 0.0051)])
def test_an_untuned_vis_order_takes_the_per_arm_yaml_value(kind: str, expected: float) -> None:
    # ACT
    error = uncertainty_of_untuned_order_nm(vis_parameters(), kind)

    # ASSERT
    assert error == expected


@pytest.mark.parametrize("block", [{}, {"sky-shift-uncertainty-nm": {"untuned": 0.1}}])
def test_an_arm_without_the_requested_yaml_value_has_none(block: dict) -> None:
    # ACT
    error = uncertainty_of_untuned_order_nm(block, ATTACHED_ARC_KEY)

    # ASSERT
    assert error is None


@pytest.mark.parametrize("bad", [None, "slow", 0.0, -0.1, float("nan")])
def test_an_unusable_yaml_value_raises(bad: object) -> None:
    # ACT AND ASSERT
    with pytest.raises(ValueError, match="untuned"):
        uncertainty_of_untuned_order_nm({"sky-shift-uncertainty-nm": {"untuned": bad}}, UNTUNED_KEY)
