"""Contracts for the pure NIR global sky-line fit that replaces per-order median shifts (DY-1968)."""

from __future__ import annotations

import copy
import pickle

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.nir_sky_refit import (
    MAX_GLOBAL_RMS_PX,
    MAX_TIGHT_TOLERANCE_NM,
    MIN_GLOBAL_LINES,
    SLOPE_SHRINK_LINES,
    apply_sky_corrections,
    fit_nir_sky_model,
    geometry_from_rows,
    global_seed_px,
    match_peaks_to_catalogue,
    match_tolerance_nm,
    order_geometry,
    read_sky_tuning_settings,
    sky_tuning_settings,
)
from soxspipe.commonutils.ug_wavelength_transfer import detector_px_to_nm_shift

pytestmark = pytest.mark.unit

ORDERS = list(range(11, 19))
ORDER_MID = 14.5
ORDER_HALF = 3.5
HALF_RANGE_NM = 40.0
PIXEL_SCALE_NM = 0.06
# KNOWN MODEL (px): A0 + A1 * m + A2 * t, WITH m THE ORDER NORMALISED TO [-1, 1] AND t THE POSITION ALONG THE ORDER
A0, A1, A2 = -3.5, 0.6, 0.4
SETTINGS = sky_tuning_settings()


def _geometry() -> dict[int, order_geometry]:
    # THE SIGN OF THE DISPERSION ALTERNATES, SO THE SIGN CONVENTION IS EXERCISED IN BOTH DIRECTIONS
    return {
        o: order_geometry(
            centreNm=1000.0 + 100.0 * (o - 11),
            halfRangeNm=HALF_RANGE_NM,
            dispersionNmPerPx=PIXEL_SCALE_NM * (1 if o % 2 else -1),
            medianPixelScaleNm=PIXEL_SCALE_NM,
        )
        for o in ORDERS
    }


def _model_px(order: int, t: np.ndarray, extraSlope: float = 0.0) -> np.ndarray:
    return A0 + A1 * (order - ORDER_MID) / ORDER_HALF + (A2 + extraSlope) * t


def _lines(
    orders=ORDERS,
    perOrder: int = 10,
    tRange=(-0.9, 0.9),
    noisePx: float = 0.02,
    extraSlopes: dict | None = None,
    seed: int = 7,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    geometry = _geometry()
    frames = []
    for o in orders:
        t = np.linspace(tRange[0], tRange[1], perOrder)
        px = _model_px(o, t, (extraSlopes or {}).get(o, 0.0)) + rng.normal(0.0, noisePx, perOrder)
        observed = geometry[o].centreNm + t * HALF_RANGE_NM
        shiftNm = detector_px_to_nm_shift(px, geometry[o].dispersionNmPerPx)
        frames.append(pd.DataFrame({"order": o, "catalogueNm": observed + shiftNm, "observedNm": observed}))
    return pd.concat(frames, ignore_index=True)


# SETTINGS


def test_settings_default_to_the_shipped_values_when_the_block_is_missing():
    # ACT
    settings = read_sky_tuning_settings(None)

    # ASSERT
    assert settings == sky_tuning_settings(globalFit=True, orderDeg=1, alongOrderDeg=1, slopeMinLines=8, clipSigma=3.0)


def test_settings_read_every_key_of_the_block():
    # ARRANGE
    block = {"global-fit": False, "order-deg": 2, "along-order-deg": 0, "slope-min-lines": 12, "clip-sigma": 2.5}

    # ACT
    settings = read_sky_tuning_settings(block)

    # ASSERT
    assert settings == sky_tuning_settings(
        globalFit=False, orderDeg=2, alongOrderDeg=0, slopeMinLines=12, clipSigma=2.5
    )


@pytest.mark.parametrize(
    "block",
    [
        {"clip-sigma": float("nan")},
        {"clip-sigma": float("inf")},
        {"clip-sigma": 0},
        {"clip-sigma": "3"},
        {"global-fit": "False"},
        {"order-deg": -1},
        {"order-deg": 1.5},
        {"order-deg": True},
        {"along-order-deg": 9},
        {"slope-min-lines": 2},
        {"unknown-key": 1},
    ],
)
def test_settings_reject_values_that_cannot_drive_a_fit(block):
    # ACT / ASSERT
    with pytest.raises(ValueError, match="sky-tuning"):
        read_sky_tuning_settings(block)


# FITTING


def test_fit_recovers_the_known_order_and_along_order_model():
    # ARRANGE
    lines = _lines()

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    assert result.fallbackReason is None
    assert set(result.corrections) == set(ORDERS)
    for o, correction in result.corrections.items():
        assert correction.shiftPx == pytest.approx(_model_px(o, np.array(0.0)), abs=0.02)
        assert correction.slope_px == pytest.approx(A2, abs=0.03)
    assert result.coefficients[0] == pytest.approx(A0, abs=0.02)
    assert result.coefficients[1] == pytest.approx(A1, abs=0.02)
    assert result.coefficients[2] == pytest.approx(A2, abs=0.02)
    assert result.rmsPx < 0.05


def test_fit_clips_outlier_lines_and_reports_the_clip_fraction():
    # ARRANGE
    lines = _lines()
    geometry = _geometry()
    badRows = [3, 17, 40, 62]
    for row in badRows:
        o = int(lines.loc[row, "order"])
        lines.loc[row, "catalogueNm"] += detector_px_to_nm_shift(2.0, geometry[o].dispersionNmPerPx)

    # ACT
    result = fit_nir_sky_model(lines, geometry, SETTINGS)

    # ASSERT
    assert result.fallbackReason is None
    assert result.coefficients[0] == pytest.approx(A0, abs=0.03)
    assert result.clipFraction == pytest.approx(len(badRows) / len(lines))
    assert sorted(result.lines.index[result.lines["clipped"]]) == badRows


def test_order_with_two_lines_takes_the_global_shift_and_slope():
    # ARRANGE
    lines = _lines(orders=[o for o in ORDERS if o != 15])
    twoLines = _lines(orders=[15], perOrder=2)
    lines = pd.concat([lines, twoLines], ignore_index=True)

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    correction = result.corrections[15]
    assert correction.nLines == 2
    assert correction.perOrderSlopePx == 0.0
    assert correction.shiftPx == pytest.approx(_model_px(15, np.array(0.0)), abs=0.03)


def test_order_without_lines_still_gets_the_model_shift():
    # ARRANGE
    lines = _lines(orders=[o for o in ORDERS if o != 18])

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    correction = result.corrections[18]
    assert correction.nLines == 0
    assert correction.shiftPx == pytest.approx(_model_px(18, np.array(0.0)), abs=0.05)


def test_well_sampled_order_gets_a_shrunk_per_order_slope():
    # ARRANGE
    extra = 0.8
    lines = _lines(extraSlopes={12: extra}, perOrder=10)
    lines = pd.concat([lines[lines["order"] != 12], _lines(orders=[12], perOrder=20, extraSlopes={12: extra})])

    # ACT
    result = fit_nir_sky_model(lines.reset_index(drop=True), _geometry(), SETTINGS)

    # ASSERT
    correction = result.corrections[12]
    shrink = correction.nLines / (correction.nLines + SLOPE_SHRINK_LINES)
    assert correction.nLines == 20
    assert 0.0 < correction.perOrderSlopePx < extra
    assert correction.perOrderSlopePx == pytest.approx(extra * shrink, rel=0.1)
    assert correction.slope_px == pytest.approx(result.coefficients[2] + correction.perOrderSlopePx, abs=1e-9)


def test_order_whose_lines_cover_less_than_half_the_order_gets_no_per_order_slope():
    # ARRANGE
    lines = _lines(orders=[o for o in ORDERS if o != 12])
    narrow = _lines(orders=[12], perOrder=14, tRange=(0.0, 0.4), extraSlopes={12: 0.8})
    lines = pd.concat([lines, narrow], ignore_index=True)

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    assert result.corrections[12].perOrderSlopePx == 0.0


def test_order_with_fewer_lines_than_the_minimum_gets_no_per_order_slope():
    # ARRANGE
    lines = _lines(orders=[o for o in ORDERS if o != 12])
    few = _lines(orders=[12], perOrder=SETTINGS.slopeMinLines - 1, extraSlopes={12: 0.8})
    lines = pd.concat([lines, few], ignore_index=True)

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    assert result.corrections[12].perOrderSlopePx == 0.0


def test_too_few_lines_falls_back_with_a_reason_and_no_corrections():
    # ARRANGE
    lines = _lines(orders=ORDERS[:2], perOrder=(MIN_GLOBAL_LINES - 1) // 2)

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    assert result.corrections == {}
    assert "sky lines" in result.fallbackReason
    assert str(MIN_GLOBAL_LINES) in result.fallbackReason


def test_scatter_above_the_limit_falls_back_with_a_reason():
    # ARRANGE
    lines = _lines(noisePx=4.0 * MAX_GLOBAL_RMS_PX, perOrder=12)

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    assert result.corrections == {}
    assert "rms" in result.fallbackReason


def test_a_single_order_cannot_fit_an_order_dependence_and_falls_back():
    # ARRANGE
    lines = _lines(orders=[14], perOrder=30)

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    assert result.corrections == {}
    assert "rank" in result.fallbackReason


def test_a_line_from_an_order_without_geometry_is_rejected():
    # ARRANGE
    lines = _lines()
    lines.loc[0, "order"] = 99

    # ACT / ASSERT
    with pytest.raises(ValueError, match="99"):
        fit_nir_sky_model(lines, _geometry(), SETTINGS)


def test_fit_leaves_the_input_lines_unchanged():
    # ARRANGE
    lines = _lines()
    before = lines.copy()

    # ACT
    fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    pd.testing.assert_frame_equal(lines, before)


def test_result_lines_carry_the_residual_before_and_the_model_after():
    # ARRANGE
    lines = _lines()

    # ACT
    result = fit_nir_sky_model(lines, _geometry(), SETTINGS)

    # ASSERT
    assert {"order", "catalogueNm", "observedNm", "residualPx", "modelPx", "clipped"} <= set(result.lines.columns)
    after = result.lines["residualPx"] - result.lines["modelPx"]
    assert np.abs(after).max() < 0.1
    assert np.abs(result.lines["residualPx"]).max() > 1.0


# CORRECTIONS


def test_correction_adds_the_nm_shift_row_by_row_from_the_local_dispersion():
    # ARRANGE
    result = fit_nir_sky_model(_lines(), _geometry(), SETTINGS)
    correction = result.corrections[12]
    geometry = _geometry()[12]
    wavelengths = geometry.centreNm + np.array([-0.5, 0.0, 0.5]) * HALF_RANGE_NM
    pixelScales = np.array([0.054, 0.06, 0.066])

    # ACT
    shifts = correction.correction_nm(wavelengths, pixelScales)

    # ASSERT
    t = np.array([-0.5, 0.0, 0.5])
    expectedPx = correction.shiftPx + correction.slope_px * t
    localDispersion = geometry.dispersionNmPerPx * pixelScales / geometry.medianPixelScaleNm
    assert shifts == pytest.approx(-localDispersion * expectedPx)
    assert correction.shift_nm == pytest.approx(shifts[1] * geometry.medianPixelScaleNm / pixelScales[1])


def test_correction_without_pixel_scales_uses_the_median_dispersion():
    # ARRANGE
    result = fit_nir_sky_model(_lines(), _geometry(), SETTINGS)
    correction = result.corrections[13]

    # ACT
    shift = correction.correction_nm(np.array([_geometry()[13].centreNm]))

    # ASSERT
    assert shift[0] == pytest.approx(correction.shift_nm)


def test_applying_corrections_returns_a_new_frame_and_leaves_unlisted_orders_alone():
    # ARRANGE
    result = fit_nir_sky_model(_lines(), _geometry(), SETTINGS)
    geometry = _geometry()
    frame = pd.DataFrame(
        {
            "order": [11, 11, 12, 99],
            "wavelengthMean": [geometry[11].centreNm, geometry[11].centreNm + 10, geometry[12].centreNm, 1500.0],
            "pixelScaleNm": PIXEL_SCALE_NM,
        }
    )
    before = frame.copy()

    # ACT
    shifted = apply_sky_corrections(frame, result.corrections)

    # ASSERT
    pd.testing.assert_frame_equal(frame, before)
    assert shifted is not frame
    expected = result.corrections[11].correction_nm(frame.loc[:1, "wavelengthMean"].to_numpy(), None)
    assert (shifted["wavelengthMean"] - frame["wavelengthMean"]).iloc[:2].to_numpy() == pytest.approx(expected)
    assert shifted.loc[3, "wavelengthMean"] == 1500.0


def test_corrections_survive_pickling_and_deep_copy():
    # ARRANGE
    result = fit_nir_sky_model(_lines(), _geometry(), SETTINGS)
    correction = result.corrections[14]

    # ACT
    unpickled = pickle.loads(pickle.dumps(correction))  # noqa: S301 - ROUND-TRIPS AN OBJECT THIS TEST BUILT
    copied = copy.deepcopy(correction)

    # ASSERT
    assert unpickled == correction
    assert copied == correction


# GEOMETRY


def test_geometry_is_read_from_an_orders_rows():
    # ARRANGE
    wavelengths = pd.Series([1000.0, 1010.0, np.nan, 1020.0])
    pixelScales = pd.Series([0.05, 0.06, np.nan, 0.07])

    # ACT
    geometry = geometry_from_rows(wavelengths, pixelScales, -0.06)

    # ASSERT
    assert geometry == order_geometry(
        centreNm=1010.0, halfRangeNm=10.0, dispersionNmPerPx=-0.06, medianPixelScaleNm=0.06
    )


def test_geometry_of_an_order_without_wavelengths_is_rejected():
    # ACT / ASSERT
    with pytest.raises(ValueError, match="wavelength"):
        geometry_from_rows(pd.Series([np.nan]), pd.Series([0.06]), 0.06)


# MATCHING


def test_peaks_match_the_nearest_catalogue_line_within_the_tolerance():
    # ARRANGE
    catalogue = np.array([1500.0, 1501.0, 1502.0])
    peaks = np.array([1500.03, 1501.5, 1501.98])

    # ACT
    peakIdx, catalogueIdx = match_peaks_to_catalogue(peaks, catalogue, toleranceNm=0.1)

    # ASSERT
    assert peakIdx.tolist() == [0, 2]
    assert catalogueIdx.tolist() == [0, 2]


def test_a_line_half_a_nanometre_away_is_not_matched_to_a_neighbour_after_seeding():
    # ARRANGE
    catalogue = np.array([1500.0, 1500.47])
    peaks = np.array([1500.0 + 0.02])

    # ACT
    peakIdx, catalogueIdx = match_peaks_to_catalogue(peaks, catalogue, toleranceNm=MAX_TIGHT_TOLERANCE_NM)

    # ASSERT
    assert catalogueIdx.tolist() == [0]


def test_a_peak_between_two_catalogue_lines_is_matched_once():
    # ARRANGE
    catalogue = np.array([1500.00, 1500.08])
    peaks = np.array([1500.05])

    # ACT
    peakIdx, catalogueIdx = match_peaks_to_catalogue(peaks, catalogue, toleranceNm=0.1)

    # ASSERT
    assert peakIdx.tolist() == [0]
    assert catalogueIdx.tolist() == [1]


def test_matching_with_no_peaks_or_no_catalogue_returns_nothing():
    # ACT
    noPeaks = match_peaks_to_catalogue(np.array([]), np.array([1500.0]), toleranceNm=0.1)
    noCatalogue = match_peaks_to_catalogue(np.array([1500.0]), np.array([]), toleranceNm=0.1)

    # ASSERT
    assert noPeaks[0].size == noPeaks[1].size == noCatalogue[0].size == noCatalogue[1].size == 0


def test_tolerance_scales_with_the_arc_rms_and_is_capped():
    # ACT
    tight = match_tolerance_nm(arcRmsPx=0.2, dispersionNmPerPx=-0.06)
    capped = match_tolerance_nm(arcRmsPx=5.0, dispersionNmPerPx=0.06)

    # ASSERT
    assert tight == pytest.approx(2.5 * 0.2 * 0.06)
    assert capped == MAX_TIGHT_TOLERANCE_NM


def test_tolerance_without_an_arc_rms_is_the_cap():
    # ACT / ASSERT
    assert match_tolerance_nm(arcRmsPx=None, dispersionNmPerPx=0.06) == MAX_TIGHT_TOLERANCE_NM
    assert match_tolerance_nm(arcRmsPx=float("nan"), dispersionNmPerPx=0.06) == MAX_TIGHT_TOLERANCE_NM


# SEEDING


def test_seed_is_the_weighted_median_shift_of_the_line_rich_orders():
    # ARRANGE: ORDER 15 IS LINE-RICH BUT ALIASED; ORDERS 20 AND 21 HAVE FEW LINES AND ARE LEFT OUT
    dispersions = dict.fromkeys([11, 12, 13, 15, 20, 21], 0.06)
    shiftsNm = {11: -0.21, 12: -0.22, 13: -0.20, 15: 0.47, 20: -0.60, 21: -0.55}
    lineCounts = {11: 20, 12: 24, 13: 22, 15: 21, 20: 2, 21: 1}

    # ACT
    seedPx = global_seed_px(shiftsNm, lineCounts, dispersions)

    # ASSERT
    assert seedPx == pytest.approx(0.21 / 0.06, abs=0.2)


def test_seed_skips_orders_without_a_dispersion_or_lines():
    # ARRANGE
    shiftsNm = {11: -0.2, 12: 0.0}

    # ACT
    seedPx = global_seed_px(shiftsNm, {11: 10, 12: 0}, {11: 0.05})

    # ASSERT
    assert seedPx == pytest.approx(4.0)


def test_seed_is_none_when_no_order_has_lines():
    # ACT / ASSERT
    assert global_seed_px({11: 0.0}, {11: 0}, {11: 0.05}) is None
