"""Contracts for the VIS u/g wavelength-transfer model (DY-1608)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from astropy.io import fits

import soxspipe
from soxspipe.commonutils.slit_illumination_qc import CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES
from soxspipe.commonutils.ug_wavelength_transfer import (
    MIN_MATCHED_LINES,
    MIN_MEASURED_SIGMA_PX,
    MIN_R_SECOND_LINE_SNR,
    MIN_SKY_LINE_SNR,
    calibration_rotator_angle,
    combine_r_shift_px,
    detector_px_to_nm_shift,
    load_r_transfer_coefficients,
    load_ug_transfer_coefficients,
    nm_shift_to_detector_px,
    predict_r_detector_shift,
    predict_ug_detector_shifts,
    r_prediction_skip_reason,
    r_second_line_snr_reason,
    sky_line_snr_reason,
    ug_skip_reason,
    usable_rotator_angle,
)

pytestmark = pytest.mark.unit

KEYWORDS = {"ABSROT": "ESO ADA ABSROT END", "PAR_ANG_START": "ESO TEL PARANG START"}


def _kw(alias: str) -> str:
    return KEYWORDS[alias]


def _header(**cards: object) -> fits.Header:
    header = fits.Header()
    for alias, value in cards.items():
        header[KEYWORDS[alias]] = value
    return header


# THE 2026-10-16 HUB FIT (output/rotated_mph_study/results/ug_transfer_coefficients.yaml)
HUB_COEFFICIENTS = {
    2: {"m": -2.714, "cos": -4.624, "sin": -5.17},
    3: {"m": -3.489, "cos": -5.233, "sin": -5.935},
}


@pytest.mark.parametrize(
    ("deltaR", "deltaI", "scienceAngle", "calibrationAngle", "expectedU", "expectedG"),
    [
        # NIGHT 20260705T023329 AGAINST DAY 20260705T111638 (results/vis_arcs.csv)
        (-0.93625, -0.45750, 54.13553, 0.0, -0.3835, -0.2112),
        # NIGHT 20260705T052629 AGAINST DAY 20260705T111638
        (0.34000, 0.51500, -19.76611, 0.0, 0.8606, 0.8239),
        # NIGHT 20260706T071809 AGAINST DAY 20260706T112958
        (2.27375, 2.61500, -204.28288, 0.0, 0.0787, -0.9661),
    ],
)
def test_predictor_reproduces_the_hub_night_against_day_pairs(
    deltaR: float,
    deltaI: float,
    scienceAngle: float,
    calibrationAngle: float,
    expectedU: float,
    expectedG: float,
) -> None:
    # ACT
    shifts = predict_ug_detector_shifts(deltaR, deltaI, scienceAngle, calibrationAngle, HUB_COEFFICIENTS)

    # ASSERT
    assert shifts == {2: pytest.approx(expectedU, abs=0.01), 3: pytest.approx(expectedG, abs=0.01)}


def test_predictor_gives_no_rotator_term_when_science_and_calibration_share_an_angle() -> None:
    # ACT
    shifts = predict_ug_detector_shifts(1.0, 3.0, 123.0, 123.0, HUB_COEFFICIENTS)

    # ASSERT
    assert shifts == {2: pytest.approx(-2.714 * 2.0), 3: pytest.approx(-3.489 * 2.0)}


def test_usable_rotator_angle_reads_a_real_angle() -> None:
    # ACT
    angle = usable_rotator_angle(_header(ABSROT=-66.03356, PAR_ANG_START=140.611), _kw)

    # ASSERT
    assert angle == pytest.approx(-66.03356)


def test_usable_rotator_angle_keeps_a_zero_angle_with_a_real_parallactic_angle() -> None:
    # ACT
    angle = usable_rotator_angle(_header(ABSROT=0.0, PAR_ANG_START=21.1), _kw)

    # ASSERT
    assert angle == 0.0


@pytest.mark.parametrize(
    "header",
    [
        pytest.param(_header(ABSROT=0.0, PAR_ANG_START=0.0), id="placeholder-0-0"),
        pytest.param(_header(PAR_ANG_START=12.0), id="missing"),
        pytest.param(_header(ABSROT="n/a", PAR_ANG_START=12.0), id="not-a-number"),
        pytest.param(_header(ABSROT=True, PAR_ANG_START=12.0), id="boolean"),
    ],
)
def test_usable_rotator_angle_rejects_unusable_angles(header: fits.Header) -> None:
    # ACT
    angle = usable_rotator_angle(header, _kw)

    # ASSERT
    assert angle is None


def test_calibration_angle_is_read_from_a_rotated_night_arc() -> None:
    # ACT
    angle = calibration_rotator_angle(_header(ABSROT=54.13553, PAR_ANG_START=21.119), _kw)

    # ASSERT
    assert angle == pytest.approx(54.13553)


@pytest.mark.parametrize(
    "header",
    [
        pytest.param(_header(), id="daytime-arc-without-absrot"),
        pytest.param(_header(ABSROT=0.0, PAR_ANG_START=0.0), id="placeholder-0-0"),
    ],
)
def test_calibration_angle_falls_back_to_the_park_angle(header: fits.Header) -> None:
    # ACT
    angle = calibration_rotator_angle(header, _kw)

    # ASSERT
    assert CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES == 0.0
    assert angle == CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES


# SIGNED DISPERSION: WAVELENGTH RISES WITH PIXEL IN r/i AND FALLS IN u/g
RI_DISPERSION_NM_PER_PX = 0.04
UG_DISPERSION_NM_PER_PX = -0.025


def test_features_moving_to_larger_pixels_in_a_rising_order_need_a_negative_nm_shift() -> None:
    # A FEATURE MOVED +2 PX NOW READS 0.08 NM TOO RED, SO THE TUNING ADDS -0.08 NM

    # ACT
    shiftNm = detector_px_to_nm_shift(2.0, RI_DISPERSION_NM_PER_PX)
    shiftPx = nm_shift_to_detector_px(-0.08, RI_DISPERSION_NM_PER_PX)

    # ASSERT
    assert shiftNm == pytest.approx(-0.08)
    assert shiftPx == pytest.approx(2.0)


def test_features_moving_to_larger_pixels_in_a_falling_order_need_a_positive_nm_shift() -> None:
    # A FEATURE MOVED +2 PX NOW READS 0.05 NM TOO BLUE, SO THE TUNING ADDS +0.05 NM

    # ACT
    shiftNm = detector_px_to_nm_shift(2.0, UG_DISPERSION_NM_PER_PX)
    shiftPx = nm_shift_to_detector_px(0.05, UG_DISPERSION_NM_PER_PX)

    # ASSERT
    assert shiftNm == pytest.approx(0.05)
    assert shiftPx == pytest.approx(2.0)


@pytest.mark.parametrize("dispersion", [0.0, float("nan"), float("inf")])
def test_conversion_refuses_a_degenerate_dispersion(dispersion: float) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match="dispersion"):
        nm_shift_to_detector_px(0.05, dispersion)
    with pytest.raises(ValueError, match="dispersion"):
        detector_px_to_nm_shift(2.0, dispersion)


def _shipped_detector_parameters() -> dict:
    path = Path(soxspipe.__file__).parent / "resources" / "soxs_detector_parameters.yaml"
    return yaml.safe_load(path.read_text())


def test_shipped_vis_parameters_carry_the_hub_coefficients_and_their_provenance() -> None:
    # ARRANGE
    visParameters = _shipped_detector_parameters()["VIS"]

    # ACT
    coeffs = load_ug_transfer_coefficients(visParameters)

    # ASSERT
    assert coeffs == HUB_COEFFICIENTS
    block = visParameters["ug-wavelength-transfer"]
    assert block["fit_date"] == "2026-10-16"
    assert block["n_arcs"] == 72
    assert block["n_nights"] == 18
    assert block["absrot_range_deg"] == [-204.7, 186.9]
    assert "rotated_mph_study" in block["source"]


def test_detector_parameters_without_the_block_give_no_coefficients() -> None:
    # ACT
    coeffs = load_ug_transfer_coefficients(_shipped_detector_parameters()["NIR"])

    # ASSERT
    assert coeffs is None


@pytest.mark.parametrize(
    "orders",
    [
        pytest.param({2: {"m": 1.0, "cos": 1.0, "sin": 1.0}}, id="missing-g"),
        pytest.param({2: {"m": 1.0, "cos": 1.0}, 3: {"m": 1.0, "cos": 1.0, "sin": 1.0}}, id="missing-term"),
        pytest.param({2: {"m": "x", "cos": 1.0, "sin": 1.0}, 3: {"m": 1.0, "cos": 1.0, "sin": 1.0}}, id="text"),
        pytest.param({2: {"m": float("nan"), "cos": 1.0, "sin": 1.0}, 3: {"m": 1.0, "cos": 1.0, "sin": 1.0}}, id="nan"),
        pytest.param({2: {"m": 1.0, "cos": 1.0, "sin": 1.0}, 3: {"m": 1.0, "cos": float("inf"), "sin": 1.0}}, id="inf"),
    ],
)
def test_a_malformed_block_is_refused(orders: dict) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match="ug-wavelength-transfer"):
        load_ug_transfer_coefficients({"ug-wavelength-transfer": {"orders": orders}})


def _stats(nLines: int, lineSnr: float) -> dict:
    return {"nLines": nLines, "lineSnr": lineSnr, "shiftPx": 0.5}


BRIGHT_SNR = 3 * MIN_SKY_LINE_SNR


def _passing_stats() -> dict:
    return {
        1: _stats(MIN_MATCHED_LINES[1], BRIGHT_SNR),
        4: _stats(MIN_MATCHED_LINES[4], BRIGHT_SNR),
    }


def test_guard_thresholds_are_the_real_data_values() -> None:
    # ASSERT
    assert MIN_MATCHED_LINES == {1: 2, 4: 10}


def test_guards_pass_with_enough_well_agreeing_lines_and_a_science_angle() -> None:
    # ACT
    reason = ug_skip_reason(_passing_stats(), scienceAngleDeg=-19.8)

    # ASSERT
    assert reason is None


@pytest.mark.parametrize(
    ("change", "expectedWords"),
    [
        pytest.param({1: None}, "order 1: no sky shift", id="r-not-measured"),
        pytest.param(
            {4: _stats(MIN_MATCHED_LINES[4] - 1, BRIGHT_SNR)}, "order 4: 9 sky lines (min 10)", id="too-few-i-lines"
        ),
        pytest.param(
            {1: _stats(MIN_MATCHED_LINES[1] - 1, BRIGHT_SNR)}, "order 1: 1 sky lines (min 2)", id="too-few-r-lines"
        ),
        pytest.param({1: _stats(10, MIN_SKY_LINE_SNR - 0.1)}, "order 1: sky-line SNR", id="faint-r-sky"),
        pytest.param({4: _stats(10, MIN_SKY_LINE_SNR - 0.1)}, "order 4: sky-line SNR", id="faint-i-sky"),
        pytest.param({4: _stats(10, float("nan"))}, "order 4: sky-line SNR nan", id="nan-snr"),
    ],
)
def test_guards_name_the_failing_sky_order(change: dict, expectedWords: str) -> None:
    # ARRANGE
    stats = _passing_stats()
    for order, value in change.items():
        if value is None:
            del stats[order]
        else:
            stats[order] = value

    # ACT
    reason = ug_skip_reason(stats, scienceAngleDeg=-19.8)

    # ASSERT
    assert reason is not None
    assert expectedWords in reason


@pytest.mark.parametrize(
    ("lineSnr", "expected"),
    [
        pytest.param(MIN_SKY_LINE_SNR, None, id="at-threshold"),
        pytest.param(MIN_SKY_LINE_SNR - 0.05, f"order 4: sky-line SNR {MIN_SKY_LINE_SNR - 0.05:.1f}", id="below"),
        pytest.param(float("nan"), "order 4: sky-line SNR nan", id="nan"),
    ],
)
def test_sky_line_snr_reason_refuses_only_a_faint_or_unknown_sky(lineSnr: float, expected: str | None) -> None:
    # ACT
    reason = sky_line_snr_reason(4, lineSnr)

    # ASSERT
    if expected is None:
        assert reason is None
    else:
        assert reason.startswith(expected)
        assert f"(min {MIN_SKY_LINE_SNR})" in reason


@pytest.mark.parametrize(
    ("secondLineSnr", "expected"),
    [
        pytest.param(MIN_R_SECOND_LINE_SNR, None, id="at-threshold"),
        pytest.param(1.8, f"order 1: 2nd sky line SNR 1.8 (min {MIN_R_SECOND_LINE_SNR})", id="noise-level"),
        pytest.param(float("nan"), f"order 1: 2nd sky line SNR nan (min {MIN_R_SECOND_LINE_SNR})", id="nan"),
    ],
)
def test_r_second_line_reason_refuses_a_noise_level_second_line(secondLineSnr: float, expected: str | None) -> None:
    # ACT
    reason = r_second_line_snr_reason(secondLineSnr)

    # ASSERT
    assert reason == expected


def test_guards_fail_without_a_usable_science_angle() -> None:
    # ACT
    reason = ug_skip_reason(_passing_stats(), scienceAngleDeg=None)

    # ASSERT
    assert reason == "no usable science ABSROT"


# THE 2026-10-19 HUB FIT (output/r_from_i_dy1645/results/r_from_i_coefficients.yaml)
HUB_R_COEFFICIENTS = {"i": 1.104, "cos": 0.554, "sin": -0.12, "rms_px": 0.364}


def test_r_predictor_reproduces_the_hub_formula() -> None:
    # ARRANGE
    # 1.104 * 0.5 + 0.554 * (cos(-19.76611 deg) - 1) - 0.12 * sin(-19.76611 deg)

    # ACT
    shift = predict_r_detector_shift(0.5, -19.76611, 0.0, HUB_R_COEFFICIENTS)

    # ASSERT
    assert shift == pytest.approx(0.5600, abs=0.001)


def test_r_predictor_gives_only_the_scaled_i_shift_when_science_and_calibration_share_an_angle() -> None:
    # ACT
    shift = predict_r_detector_shift(2.0, 123.0, 123.0, HUB_R_COEFFICIENTS)

    # ASSERT
    assert shift == pytest.approx(1.104 * 2.0)


def test_shipped_vis_parameters_carry_the_r_from_i_coefficients_and_their_provenance() -> None:
    # ARRANGE
    visParameters = _shipped_detector_parameters()["VIS"]

    # ACT
    coeffs = load_r_transfer_coefficients(visParameters)

    # ASSERT
    assert coeffs == HUB_R_COEFFICIENTS
    block = visParameters["r-from-i-transfer"]
    assert block["fit_date"] == "2026-10-19"
    assert block["n_arcs"] == 72
    assert block["n_nights"] == 18
    assert "output/r_from_i_dy1645/fit_r_from_i.py" in block["source"]


def test_detector_parameters_without_the_r_block_give_no_r_coefficients() -> None:
    # ACT
    coeffs = load_r_transfer_coefficients(_shipped_detector_parameters()["NIR"])

    # ASSERT
    assert coeffs is None


@pytest.mark.parametrize(
    "block",
    [
        pytest.param({"i": 1.0, "cos": 1.0, "sin": 1.0}, id="missing-rms"),
        pytest.param({"i": "x", "cos": 1.0, "sin": 1.0, "rms_px": 0.3}, id="text"),
        pytest.param({"i": float("nan"), "cos": 1.0, "sin": 1.0, "rms_px": 0.3}, id="nan"),
        pytest.param({"i": 1.0, "cos": float("inf"), "sin": 1.0, "rms_px": 0.3}, id="inf"),
        pytest.param({"i": 1.0, "cos": 1.0, "sin": 1.0, "rms_px": 0.0}, id="zero-rms"),
        pytest.param("not a mapping", id="not-a-mapping"),
    ],
)
def test_a_malformed_r_block_is_refused(block: object) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match="r-from-i-transfer"):
        load_r_transfer_coefficients({"r-from-i-transfer": block})


def test_r_prediction_guards_pass_for_a_well_measured_i_and_a_science_angle() -> None:
    # ACT
    reason = r_prediction_skip_reason(_stats(MIN_MATCHED_LINES[4], BRIGHT_SNR), scienceAngleDeg=-19.8)

    # ASSERT
    assert reason is None


@pytest.mark.parametrize(
    ("iStats", "angle", "expectedWords"),
    [
        pytest.param(None, -19.8, "order 4: no sky shift", id="i-not-measured"),
        pytest.param(
            _stats(MIN_MATCHED_LINES[4] - 1, BRIGHT_SNR), -19.8, "order 4: 9 sky lines (min 10)", id="few-i-lines"
        ),
        pytest.param(_stats(10, MIN_SKY_LINE_SNR - 0.1), -19.8, "order 4: sky-line SNR", id="faint-i-sky"),
        pytest.param(_stats(10, BRIGHT_SNR), None, "no usable science ABSROT", id="no-angle"),
    ],
)
def test_r_prediction_guards_name_the_failure(iStats: object, angle: object, expectedWords: str) -> None:
    # ACT
    reason = r_prediction_skip_reason(iStats, scienceAngleDeg=angle)

    # ASSERT
    assert reason is not None
    assert expectedWords in reason


def test_inverse_variance_combination_weights_the_estimates_by_their_sigmas() -> None:
    # ARRANGE
    # SIGMA 0.2 AGAINST 0.4 GIVES WEIGHTS 25 AND 6.25, SO THE MEAN SITS 80 PERCENT OF THE WAY TO THE MEASUREMENT

    # ACT
    combined = combine_r_shift_px(1.0, 0.2, 2.0, 0.4)

    # ASSERT
    assert combined == pytest.approx(1.2)


def test_inverse_variance_combination_floors_a_vanishing_measured_sigma() -> None:
    # ACT
    combined = combine_r_shift_px(1.0, 0.0, 2.0, 0.364)

    # ASSERT
    weightMeasured = 1 / MIN_MEASURED_SIGMA_PX**2
    weightPredicted = 1 / 0.364**2
    assert combined == pytest.approx(
        (weightMeasured * 1.0 + weightPredicted * 2.0) / (weightMeasured + weightPredicted)
    )
