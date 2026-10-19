"""Behaviour contracts for the robust log-space spline fit of the instrument response."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from scipy.interpolate import BSpline, LSQUnivariateSpline

from soxspipe.commonutils.response_fit import (
    _interior_knots,
    _ResponseFitConvergenceError,
    bin_log_response,
    evaluate_log_spline,
    fit_log_spline,
    fit_response_curve,
    grow_telluric_wings,
    mask_response_samples,
    resolve_response_settings,
)

pytestmark = pytest.mark.unit

BIN_WIDTH = 2.0
KNOT_SPACING = 40.0


def _truth(wave: np.ndarray) -> np.ndarray:
    """Return a smooth response that rises tenfold from 350 to 850 nm, steepest in the blue."""
    scale = 1.0 - np.exp(-500.0 / 150.0)
    return np.exp(np.log(10.0) * (1.0 - np.exp(-(wave - 350.0) / 150.0)) / scale)


def _noisy_ratio(wave: np.ndarray, noise: float, seed: int = 7) -> np.ndarray:
    return _truth(wave) * (1.0 + noise * np.random.default_rng(seed).standard_normal(wave.size))


def _fit(
    wave: np.ndarray,
    ratio: np.ndarray,
    telluricRegions: list[list[float]] | None = None,
    log: object = None,
    maxIterations: int = 10,
    penalty: float = 0.0,
) -> tuple[object, np.ndarray, tuple[float, float], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    mask = mask_response_samples(wave, ratio, telluricRegions or [], [], 2.5)
    binWave, binLog, binErr = bin_log_response(wave, ratio, mask, BIN_WIDTH, 5)
    spline, kept = fit_log_spline(
        binWave, binLog, binErr, KNOT_SPACING, 4.0, 3.0, maxIterations, log=log, penalty=penalty
    )
    fitRange = (float(binWave[kept].min()), float(binWave[kept].max()))
    return spline, kept, fitRange, (binWave, binLog, binErr)


class _RecordingLog:
    def __init__(self) -> None:
        self.warnings: list[str] = []

    def warning(self, message: str) -> None:
        self.warnings.append(message)


def test_recovers_a_steep_blue_rise_within_two_percent_at_the_blue_edge() -> None:
    # ARRANGE
    wave = np.linspace(350.0, 850.0, 20000)
    ratio = _noisy_ratio(wave, noise=0.10)

    # ACT
    spline, _, fitRange, _ = _fit(wave, ratio)
    modelled = evaluate_log_spline(spline, np.array([352.0]), fitRange)[0]

    # ASSERT
    assert modelled == pytest.approx(_truth(np.array([352.0]))[0], rel=0.02)


def test_one_sided_absorption_bands_bias_the_continuum_by_less_than_one_percent() -> None:
    # ARRANGE: UNMASKED BANDS RAISE THE RATIO BY 30 TO 100 PERCENT
    wave = np.linspace(350.0, 850.0, 20000)
    ratio = _noisy_ratio(wave, noise=0.03)
    bands = [(520.0, 528.0, 1.8), (600.0, 612.0, 1.6), (680.0, 690.0, 2.0), (760.0, 768.0, 1.3)]
    outsideBands = np.ones(wave.size, dtype=bool)
    for lower, upper, factor in bands:
        inBand = (wave >= lower) & (wave <= upper)
        ratio[inBand] *= factor
        outsideBands &= ~((wave >= lower - 4.0) & (wave <= upper + 4.0))

    # ACT
    spline, _, fitRange, _ = _fit(wave, ratio)
    modelled = evaluate_log_spline(spline, wave, fitRange)

    # ASSERT
    bias = np.abs(modelled[outsideBands] / _truth(wave[outsideBands]) - 1.0)
    assert bias.max() < 0.01


def test_a_masked_gap_of_150_nm_is_bridged_without_overshoot() -> None:
    # ARRANGE
    wave = np.linspace(350.0, 850.0, 20000)
    ratio = _noisy_ratio(wave, noise=0.03)

    # ACT
    spline, _, fitRange, _ = _fit(wave, ratio, telluricRegions=[[550.0, 700.0]])
    gap = np.linspace(555.0, 695.0, 200)
    modelled = evaluate_log_spline(spline, gap, fitRange)

    # ASSERT
    assert np.abs(modelled / _truth(gap) - 1.0).max() < 0.02


@pytest.mark.parametrize("penalty", [0.0, 1.0], ids=["least-squares", "penalised"])
def test_extrapolation_beyond_the_data_is_finite_monotonic_and_close_to_a_log_linear_continuation(
    penalty: float,
) -> None:
    # ARRANGE
    wave = np.linspace(400.0, 850.0, 20000)
    ratio = _noisy_ratio(wave, noise=0.03)
    spline, _, fitRange, _ = _fit(wave, ratio, penalty=penalty)
    blue = np.linspace(370.0, 400.0, 31)
    red = np.linspace(850.0, 880.0, 31)

    # ACT
    blueModel = evaluate_log_spline(spline, blue, fitRange)
    redModel = evaluate_log_spline(spline, red, fitRange)

    # ASSERT
    assert np.isfinite(blueModel).all()
    assert np.isfinite(redModel).all()
    assert np.all(np.diff(blueModel) > 0)
    assert np.all(np.diff(redModel) >= 0) or np.all(np.diff(redModel) <= 0)
    for edge, grid, model in ((fitRange[0], blue, blueModel), (fitRange[1], red, redModel)):
        slope = (np.log(_truth(np.array([edge + 1.0]))) - np.log(_truth(np.array([edge - 1.0])))) / 2.0
        continuation = np.exp(np.log(_truth(np.array([edge]))) + slope * (grid - edge))
        assert np.abs(model / continuation - 1.0).max() < 0.20


def test_too_few_bins_raises_a_convergence_error() -> None:
    # ARRANGE
    binWave = np.array([400.0, 410.0, 420.0])

    # ACT AND ASSERT
    with pytest.raises(_ResponseFitConvergenceError):
        fit_log_spline(binWave, np.zeros(3), np.full(3, 0.01), KNOT_SPACING, 4.0, 3.0, 5)


def test_all_non_finite_input_raises_with_the_no_finite_samples_message() -> None:
    # ARRANGE: EVERY SAMPLE HAD ZERO EXTRACTED FLUX
    wave = np.linspace(1000.0, 2000.0, 400)
    ratio = np.full(400, np.inf)

    # ACT
    mask = mask_response_samples(wave, ratio, [], [], 5.0)
    binWave, binLog, binErr = bin_log_response(wave, ratio, mask, 5.0, 5)

    # ASSERT
    with pytest.raises(_ResponseFitConvergenceError, match="no finite samples"):
        fit_log_spline(binWave, binLog, binErr, 60.0, 4.0, 2.5, 5)


def test_cap_hit_refit_returns_a_spline_that_fits_the_returned_kept_bins() -> None:
    # ARRANGE: A BLOCK OF HIGH BINS KEEPS THE CLIP CHANGING PAST THE SINGLE ALLOWED PASS
    wave = np.linspace(350.0, 850.0, 20000)
    ratio = _noisy_ratio(wave, noise=0.02)
    ratio[(wave > 500.0) & (wave < 540.0)] *= 1.5
    recordingLog = _RecordingLog()

    # ACT
    spline, kept, _, (binWave, binLog, binErr) = _fit(wave, ratio, log=recordingLog, maxIterations=1)

    # ASSERT: THE SPLINE EQUALS A DIRECT WEIGHTED FIT OF THE KEPT BINS ON THE SAME KNOTS
    assert (~kept).any()
    assert len(recordingLog.warnings) == 1
    assert "iteration cap" in recordingLog.warnings[0]
    direct = LSQUnivariateSpline(binWave[kept], binLog[kept], spline.get_knots()[1:-1], w=1.0 / binErr[kept], k=3)
    np.testing.assert_allclose(spline(binWave[kept]), direct(binWave[kept]), atol=1e-10)


def test_converging_before_the_cap_logs_no_warning() -> None:
    # ARRANGE
    wave = np.linspace(350.0, 850.0, 20000)
    ratio = _noisy_ratio(wave, noise=0.02)
    recordingLog = _RecordingLog()

    # ACT
    _fit(wave, ratio, log=recordingLog)

    # ASSERT
    assert recordingLog.warnings == []


def test_mask_drops_non_finite_non_positive_telluric_and_stellar_line_samples() -> None:
    # ARRANGE
    wave = np.array([500.0, 501.0, np.nan, 503.0, 700.0, 656.0, 800.0, 900.0])
    ratio = np.array([1.0, np.inf, 1.0, -1.0, 1.0, 1.0, 0.0, 1.0])

    # ACT
    mask = mask_response_samples(wave, ratio, [[690.0, 710.0]], [656.3], 2.5)

    # ASSERT: A NEGATIVE RATIO (NEGATIVE NOISY COUNTS) IS KEPT; A ZERO RATIO IS NOT
    assert mask.tolist() == [True, False, False, True, False, False, False, True]


def test_a_low_signal_bin_with_negative_count_samples_is_not_biased() -> None:
    # ARRANGE: COUNTS WITH SIGNAL-TO-NOISE OF 1.5 PER SAMPLE, SO MANY SAMPLES ARE NEGATIVE
    sampleCount = 4000
    wave = np.linspace(2000.0, 2001.9, sampleCount)
    trueResponse = 40.0
    absoluteFlux = 100.0
    counts = absoluteFlux / trueResponse * (1.0 + np.random.default_rng(3).standard_normal(sampleCount) / 1.5)
    ratio = absoluteFlux / counts
    mask = mask_response_samples(wave, ratio, [], [], 5.0)

    # ACT
    _, binLog, _ = bin_log_response(wave, ratio, mask, 2.0, 5)

    # ASSERT
    assert np.exp(binLog[0]) == pytest.approx(trueResponse, rel=0.05)


def test_bins_below_the_minimum_sample_count_are_dropped_and_errors_have_a_floor() -> None:
    # ARRANGE: A FULL BIN OF IDENTICAL VALUES THEN A BIN OF TWO SAMPLES
    wave = np.concatenate([np.linspace(400.0, 401.9, 10), [410.0, 410.5]])
    ratio = np.concatenate([np.full(10, 2.0), [3.0, 3.0]])
    mask = np.ones(wave.size, dtype=bool)

    # ACT
    binWave, binLog, binErr = bin_log_response(wave, ratio, mask, 2.0, 5)

    # ASSERT
    assert binWave.size == 1
    assert binLog[0] == pytest.approx(np.log(2.0))
    assert binErr[0] > 0.0


@pytest.mark.parametrize(
    ("badWavelength", "badResponse"),
    [(1500.0, np.inf), (1500.0, np.nan), (np.nan, 2.0)],
    ids=["infinite-response", "nan-response", "nan-wavelength"],
)
def test_non_finite_samples_are_dropped_with_one_warning(badWavelength: float, badResponse: float) -> None:
    # ARRANGE: ONE NON-FINITE SAMPLE IN AN OTHERWISE FLAT RESPONSE
    wave = np.linspace(1000.0, 2000.0, 400)
    ratio = np.full(400, 2.0)
    wave[200] = badWavelength
    ratio[200] = badResponse
    recordingLog = _RecordingLog()

    # ACT
    mask = mask_response_samples(wave, ratio, [], [], 5.0, log=recordingLog)

    # ASSERT
    assert int(mask.sum()) == 399
    assert not mask[200]
    assert len(recordingLog.warnings) == 1
    assert "Dropped 1 raw response samples" in recordingLog.warnings[0]


def test_noisy_bins_at_a_low_signal_end_are_judged_against_their_own_noise() -> None:
    # ARRANGE: 2 PERCENT SAMPLE NOISE UP TO 800 NM, THEN 60 PERCENT NOISE IN THE LOW-SIGNAL RED END
    wave = np.linspace(350.0, 900.0, 22000)
    noise = np.where(wave < 800.0, 0.02, 0.60)
    counts = (1.0 + noise * np.random.default_rng(11).standard_normal(wave.size)) / _truth(wave)
    ratio = 1.0 / counts

    # ACT
    _, kept, _, (binWave, _, _) = _fit(wave, ratio)
    noisyEnd = binWave > 800.0

    # ASSERT: AT MOST A FEW NOISY BINS ARE CLIPPED, AS EXPECTED FROM NOISE ALONE
    assert np.count_nonzero(~kept[noisyEnd]) / np.count_nonzero(noisyEnd) < 0.10


def _synthetic_response_fit() -> tuple[object, np.ndarray, np.ndarray]:
    """Fit the VIS defaults to a smooth synthetic response and return the fit and the caller's input arrays."""
    wave = np.linspace(400.0, 800.0, 2000)
    ratio = np.full(wave.size, 2.0)
    fitSettings = resolve_response_settings({}, "VIS", "soxs")
    return fit_response_curve(wave, ratio, fitSettings), wave, ratio


def test_a_response_fit_does_not_change_when_the_callers_input_arrays_are_mutated() -> None:
    # ARRANGE
    fit, wave, ratio = _synthetic_response_fit()
    waveBefore, ratioBefore = fit.wave.copy(), fit.ratio.copy()

    # ACT
    wave[:] = 0.0
    ratio[:] = 99.0

    # ASSERT
    np.testing.assert_array_equal(fit.wave, waveBefore)
    np.testing.assert_array_equal(fit.ratio, ratioBefore)


@pytest.mark.parametrize("field", ["wave", "ratio", "mask", "binWave", "binLog", "keptMask"])
def test_writing_into_a_response_fit_array_raises_value_error(field: str) -> None:
    # ARRANGE
    fit, _, _ = _synthetic_response_fit()

    # ACT / ASSERT
    with pytest.raises(ValueError, match="read-only"):
        getattr(fit, field)[0] = getattr(fit, field)[1]


def test_unusable_bins_do_not_raise_numpy_warnings_and_leave_usable_bins_unchanged() -> None:
    # ARRANGE: A CLEAN FIT, THEN THE SAME BINS WITH THREE UNUSABLE ONES ADDED (ZERO, NEGATIVE AND NAN ERRORS)
    wave = np.linspace(350.0, 850.0, 20000)
    _, keptClean, _, (binWave, binLog, binErr) = _fit(wave, _noisy_ratio(wave, noise=0.05))
    badWave = np.array([351.1, 352.3, 353.5])
    badErr = np.array([0.0, -1.0, np.nan])
    allWave = np.concatenate([binWave, badWave])
    allLog = np.concatenate([binLog, np.zeros(3)])
    allErr = np.concatenate([binErr, badErr])
    order = np.argsort(allWave)

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _, kept = fit_log_spline(allWave[order], allLog[order], allErr[order], KNOT_SPACING, 4.0, 3.0, 10)

    # ASSERT
    isBad = np.isin(allWave[order], badWave)
    assert not kept[isBad].any()
    np.testing.assert_array_equal(kept[~isBad], keptClean)


WING_MAX_NM = 60.0
WING_REF_NM = 40.0
WING_THRESHOLD = 0.03
WING_BIN_ERROR = 0.005
WING_REGION = [[1340.0, 1460.0]]
WING_DECAY_NM = 15.0
WING_AMPLITUDE = 0.2


class _DebugLog:
    def __init__(self) -> None:
        self.debugs: list[str] = []

    def debug(self, message: str) -> None:
        self.debugs.append(message)


def _wing_bins(
    lowerEdge: float = 1340.0, upperEdge: float = 1460.0, slope: float = 0.0, wingSign: float = 1.0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return 5 nm bins either side of a masked region, on a straight log continuum with exponential wings."""
    blue = np.arange(lowerEdge - 200.0, lowerEdge, 5.0) + 2.5
    red = np.arange(upperEdge, upperEdge + 200.0, 5.0) + 2.5
    binWave = np.concatenate([blue, red])
    distance = np.where(binWave < lowerEdge, lowerEdge - binWave, binWave - upperEdge)
    binLog = 0.5 + slope * (binWave - lowerEdge) + wingSign * WING_AMPLITUDE * np.exp(-distance / WING_DECAY_NM)
    return binWave, binLog, np.full(binWave.size, WING_BIN_ERROR)


def _grow(binWave: np.ndarray, binLog: np.ndarray, binErr: np.ndarray, log: object = None) -> np.ndarray:
    return grow_telluric_wings(binWave, binLog, binErr, WING_REGION, WING_MAX_NM, WING_REF_NM, WING_THRESHOLD, log=log)


def test_wing_bins_raised_above_the_continuum_are_flagged_and_bins_beyond_the_wing_are_not() -> None:
    # ARRANGE: EXPONENTIAL WINGS FALL TO THE 0.03 THRESHOLD ABOUT 28 NM FROM THE REGION EDGE
    binWave, binLog, binErr = _wing_bins()
    distance = np.where(binWave < 1340.0, 1340.0 - binWave, binWave - 1460.0)
    rise = WING_AMPLITUDE * np.exp(-distance / WING_DECAY_NM)

    # ACT
    wingMask = _grow(binWave, binLog, binErr)

    # ASSERT: BINS WELL ABOVE THE THRESHOLD ARE FLAGGED, BINS WELL BELOW IT AND BEYOND THE REACH ARE NOT
    assert wingMask.dtype == bool
    assert wingMask[rise > 1.5 * WING_THRESHOLD].all()
    assert not wingMask[rise < 0.5 * WING_THRESHOLD].any()
    assert not wingMask[distance > WING_MAX_NM].any()


def test_a_genuine_linear_slope_in_the_response_is_not_flagged() -> None:
    # ARRANGE: A STEEP BUT STRAIGHT CONTINUUM AND NO WINGS
    binWave, binLog, binErr = _wing_bins(slope=-0.002)
    binLog = 0.5 - 0.002 * (binWave - 1340.0)

    # ACT
    wingMask = _grow(binWave, binLog, binErr)

    # ASSERT
    assert not wingMask.any()


def test_a_side_with_fewer_than_three_reference_bins_is_skipped_without_error() -> None:
    # ARRANGE: THE BLUE REFERENCE WINDOW [1240, 1280) HOLDS ONLY TWO BINS; THE RED SIDE IS COMPLETE
    binWave, binLog, binErr = _wing_bins()
    inBlueReference = (binWave >= 1240.0) & (binWave < 1280.0)
    dropped = inBlueReference & (np.arange(binWave.size) % 4 != 0)
    keep = ~dropped
    recordingLog = _DebugLog()

    # ACT
    wingMask = _grow(binWave[keep], binLog[keep], binErr[keep], log=recordingLog)

    # ASSERT
    keptWave = binWave[keep]
    assert np.count_nonzero((keptWave >= 1240.0) & (keptWave < 1280.0)) == 2
    assert not wingMask[keptWave < 1340.0].any()
    assert wingMask[keptWave > 1460.0].any()
    assert len(recordingLog.debugs) == 1


def test_a_negative_excursion_next_to_a_region_is_never_flagged() -> None:
    # ARRANGE: EXPONENTIAL WINGS THAT DIP BELOW THE CONTINUUM
    binWave, binLog, binErr = _wing_bins(wingSign=-1.0)

    # ACT
    wingMask = _grow(binWave, binLog, binErr)

    # ASSERT
    assert not wingMask.any()


def test_a_wing_max_of_zero_disables_growth() -> None:
    # ARRANGE
    binWave, binLog, binErr = _wing_bins()

    # ACT
    wingMask = grow_telluric_wings(binWave, binLog, binErr, WING_REGION, 0.0, WING_REF_NM, WING_THRESHOLD)

    # ASSERT
    assert wingMask.shape == binWave.shape
    assert not wingMask.any()


def test_the_walk_stops_at_the_first_bin_that_is_not_raised() -> None:
    # ARRANGE: A SECOND RAISED BIN BEHIND A FLAT BIN BELONGS TO THE CONTINUUM, NOT THE WING
    binWave, binLog, binErr = _wing_bins()
    binLog = np.full(binWave.size, 0.5)
    nearest = np.flatnonzero(binWave < 1340.0)[-1]
    binLog[nearest] += 0.2
    binLog[nearest - 2] += 0.2

    # ACT
    wingMask = _grow(binWave, binLog, binErr)

    # ASSERT
    assert np.flatnonzero(wingMask).tolist() == [nearest]


def test_an_unusable_bin_between_the_edge_and_a_raised_bin_stops_the_walk() -> None:
    # ARRANGE: THE BIN NEAREST THE BLUE EDGE HAS NO ERROR, SO THE RAISED BINS BEHIND IT ARE NOT CONTIGUOUS
    binWave, binLog, binErr = _wing_bins()
    nearestBlue = np.flatnonzero(binWave < 1340.0)[-1]
    binErr[nearestBlue] = np.nan

    # ACT
    wingMask = _grow(binWave, binLog, binErr)

    # ASSERT
    assert not wingMask[binWave < 1340.0].any()
    assert wingMask[binWave > 1460.0].any()


def test_a_missing_bin_hole_in_the_wing_stops_the_walk() -> None:
    # ARRANGE: A STELLAR-LINE-STYLE HOLE REMOVES THE TWO BINS AT 1327.5 AND 1322.5 NM
    binWave, binLog, binErr = _wing_bins()
    keep = ~np.isin(binWave, [1327.5, 1322.5])

    # ACT
    wingMask = _grow(binWave[keep], binLog[keep], binErr[keep])

    # ASSERT: ONLY THE TWO BINS BEFORE THE HOLE ARE FLAGGED ON THE BLUE SIDE
    assert binWave[keep][wingMask & (binWave[keep] < 1340.0)].tolist() == [1332.5, 1337.5]


def test_a_bin_inside_a_telluric_region_is_not_used_as_reference_or_flagged() -> None:
    # ARRANGE: A STRAY BIN INSIDE THE REGION THAT SITS FAR ABOVE THE CONTINUUM
    binWave, binLog, binErr = _wing_bins(wingSign=0.0)
    binWave = np.append(binWave, 1400.0)
    binLog = np.append(binLog, 5.0)
    binErr = np.append(binErr, WING_BIN_ERROR)

    # ACT
    wingMask = _grow(binWave, binLog, binErr)

    # ASSERT
    assert not wingMask.any()


WING_FIT_SETTINGS = {
    "telluric_regions": [[600.0, 720.0]],
    "stellar_lines": [],
    "line_half_width_nm": 2.5,
    "bin_width_nm": 5.0,
    "min_bin_samples": 5,
    "knot_spacing_nm": 60.0,
    "clip_sigma_low": 4.0,
    "clip_sigma_high": 2.5,
    "max_iteration": 10,
    "telluric_wing_max_nm": 60.0,
    "telluric_wing_ref_nm": 40.0,
    "telluric_wing_threshold": 0.03,
}


def _gentle_truth(wave: np.ndarray) -> np.ndarray:
    """Return a response whose log is nearly straight over a wing reference window."""
    offset = wave - 350.0
    return np.exp(1.0 - 0.002 * offset + 1.0e-6 * offset**2)


def _ratio_with_wings() -> tuple[np.ndarray, np.ndarray]:
    """Return a noisy response with a 120 nm telluric gap whose wings raise the ratio by up to 25 percent."""
    wave = np.linspace(350.0, 850.0, 10000)
    ratio = _gentle_truth(wave) * (1.0 + 0.05 * np.random.default_rng(5).standard_normal(wave.size))
    distance = np.where(wave < 600.0, 600.0 - wave, np.where(wave > 720.0, wave - 720.0, np.inf))
    return wave, ratio * (1.0 + 0.25 * np.exp(-distance / 10.0))


def _gap_error(fit: object) -> float:
    gap = np.linspace(605.0, 715.0, 200)
    return float(np.abs(fit.evaluate(gap) / _gentle_truth(gap) - 1.0).max())


def test_a_telluric_gap_with_raised_wings_is_bridged_within_two_percent() -> None:
    # ARRANGE
    wave, ratio = _ratio_with_wings()

    # ACT
    fit = fit_response_curve(wave, ratio, WING_FIT_SETTINGS)

    # ASSERT
    assert _gap_error(fit) < 0.02
    assert fit.wingMask.any()
    assert not fit.keptMask[fit.wingMask].any()


def test_the_same_gap_overshoots_when_wing_growth_is_disabled() -> None:
    # ARRANGE
    wave, ratio = _ratio_with_wings()
    disabled = {**WING_FIT_SETTINGS, "telluric_wing_max_nm": 0}

    # ACT
    fit = fit_response_curve(wave, ratio, disabled)

    # ASSERT
    assert _gap_error(fit) > 0.03
    assert not fit.wingMask.any()


def test_a_fit_settings_block_without_wing_keys_fits_with_growth_disabled() -> None:
    # ARRANGE
    wave, ratio = _ratio_with_wings()
    legacy = {key: value for key, value in WING_FIT_SETTINGS.items() if not key.startswith("telluric_wing")}

    # ACT
    fit = fit_response_curve(wave, ratio, legacy)

    # ASSERT
    assert not fit.wingMask.any()


def test_the_exclude_mask_removes_bins_from_the_fit_from_the_start() -> None:
    # ARRANGE
    wave = np.linspace(350.0, 850.0, 20000)
    _, _, _, (binWave, binLog, binErr) = _fit(wave, _noisy_ratio(wave, noise=0.03))
    excluded = (binWave > 500.0) & (binWave < 540.0)

    # ACT
    _, kept = fit_log_spline(binWave, binLog, binErr, KNOT_SPACING, 4.0, 3.0, 10, excludeMask=excluded)

    # ASSERT
    assert not kept[excluded].any()
    assert kept[~excluded].all()


PENALTY_KNOT_SPACING = 20.0
PENALTY_BIN_ERROR = 0.01


def _penalty_truth(wave: np.ndarray) -> np.ndarray:
    """Return a log response with a gentle curve, in log units."""
    offset = wave - 1000.0
    return 1.5 - 0.0012 * offset + 4.0e-7 * offset**2


def _penalty_bins(
    segments: list[tuple[float, float]], noise: float = 0.0, seed: int = 4
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return 5 nm bins over the given wavelength segments, on the penalty truth plus Gaussian noise."""
    binWave = np.concatenate([np.arange(lower, upper, 5.0) + 2.5 for lower, upper in segments])
    binLog = _penalty_truth(binWave) + noise * np.random.default_rng(seed).standard_normal(binWave.size)
    return binWave, binLog, np.full(binWave.size, PENALTY_BIN_ERROR)


def test_a_penalty_of_zero_gives_the_same_curve_as_a_plain_least_squares_spline() -> None:
    # ARRANGE
    wave = np.linspace(350.0, 850.0, 20000)
    _, keptClean, _, (binWave, binLog, binErr) = _fit(wave, _noisy_ratio(wave, noise=0.03))
    knots = _interior_knots(binWave, KNOT_SPACING)
    plain = LSQUnivariateSpline(binWave, binLog, knots, w=1.0 / binErr, k=3)

    # ACT
    spline, _ = fit_log_spline(binWave, binLog, binErr, KNOT_SPACING, 4.0, 3.0, 1, penalty=0.0)

    # ASSERT: ONE ITERATION IS THE FIRST FIT, WHICH IS THE PLAIN SPLINE ON ALL USABLE BINS
    assert isinstance(spline, LSQUnivariateSpline)
    np.testing.assert_array_equal(spline(binWave), plain(binWave))


def test_a_penalised_curve_crosses_a_200_nm_masked_gap_close_to_a_straight_log_line() -> None:
    # ARRANGE: A NOISY RESPONSE WITH NO BINS BETWEEN 1400 AND 1600 NM AND 20 NM KNOTS
    binWave, binLog, binErr = _penalty_bins([(1000.0, 1400.0), (1600.0, 2000.0)], noise=0.01)
    gap = np.linspace(1410.0, 1590.0, 50)
    edgeLine = np.interp(gap, [1400.0, 1600.0], _penalty_truth(np.array([1400.0, 1600.0])))

    # ACT
    penalised, _ = fit_log_spline(binWave, binLog, binErr, PENALTY_KNOT_SPACING, 4.0, 3.0, 10, penalty=1.0)
    plain, _ = fit_log_spline(binWave, binLog, binErr, PENALTY_KNOT_SPACING, 4.0, 3.0, 10)

    # ASSERT: THE TRUTH CURVES BY ABOUT 0.016 IN THE GAP, SO THE STRAIGHT EDGE LINE IS THE NATURAL BRIDGE
    penalisedError = np.abs(penalised(gap) - edgeLine).max()
    assert penalisedError < 0.03
    assert np.abs(plain(gap) - edgeLine).max() > 2.0 * penalisedError


def test_a_penalised_fit_recovers_a_sparse_end_beyond_a_180_nm_gap() -> None:
    # ARRANGE: A DENSE BODY, FOUR BINS, A 180 NM GAP, THEN TWO BINS AT THE END (THE K-BAND PATTERN)
    binWave, binLog, binErr = _penalty_bins([(1000.0, 1700.0), (1700.0, 1720.0), (1900.0, 1910.0)], noise=0.005)
    endWave = np.array([1902.5, 1907.5])

    # ACT
    spline, kept = fit_log_spline(binWave, binLog, binErr, PENALTY_KNOT_SPACING, 4.0, 3.0, 10, penalty=1.0)

    # ASSERT
    assert kept.all()
    np.testing.assert_allclose(np.exp(spline(endWave)), np.exp(_penalty_truth(endWave)), rtol=0.03)
    assert np.exp(spline(1905.0)) == pytest.approx(np.exp(_penalty_truth(np.array([1905.0]))[0]), rel=0.03)


def test_a_penalised_spline_supports_derivatives_and_log_linear_extrapolation() -> None:
    # ARRANGE
    binWave, binLog, binErr = _penalty_bins([(1000.0, 1400.0)])
    spline, kept = fit_log_spline(binWave, binLog, binErr, PENALTY_KNOT_SPACING, 4.0, 3.0, 10, penalty=1.0)
    fitRange = (float(binWave[kept].min()), float(binWave[kept].max()))
    beyond = np.array([1420.0, 1440.0])

    # ACT
    slope = spline(1200.0, nu=1)
    model = evaluate_log_spline(spline, beyond, fitRange)

    # ASSERT: THE TRUE LOG SLOPE AT 1200 NM IS -0.0012 + 2 * 4E-7 * 200
    assert float(slope) == pytest.approx(-0.0012 + 8.0e-7 * 200.0, abs=2e-4)
    edge = fitRange[1]
    expected = np.exp(spline(edge) + spline(edge, nu=1) * (beyond - edge))
    np.testing.assert_allclose(model, expected)


def test_a_very_large_penalty_tends_to_the_weighted_straight_line_in_log() -> None:
    # ARRANGE
    binWave, binLog, binErr = _penalty_bins([(1000.0, 2000.0)], noise=0.01)
    line = np.polyfit(binWave, binLog, 1, w=1.0 / binErr)

    # ACT
    spline, _ = fit_log_spline(binWave, binLog, binErr, PENALTY_KNOT_SPACING, 4.0, 3.0, 10, penalty=1.0e9)

    # ASSERT
    np.testing.assert_allclose(spline(binWave), np.polyval(line, binWave), atol=2e-3)


def test_a_penalised_fit_with_too_few_bins_raises_a_convergence_error() -> None:
    # ARRANGE
    binWave = np.array([1000.0, 1010.0, 1020.0])

    # ACT AND ASSERT
    with pytest.raises(_ResponseFitConvergenceError, match="too few bins"):
        fit_log_spline(binWave, np.zeros(3), np.full(3, 0.01), PENALTY_KNOT_SPACING, 4.0, 3.0, 5, penalty=1.0)


def test_a_penalised_fit_of_bins_at_one_wavelength_raises_a_convergence_error() -> None:
    # ARRANGE: FOUR BINS AT ONE WAVELENGTH CANNOT FIX THE STRAIGHT LINE THE PENALTY LEAVES FREE
    binWave = np.full(4, 1000.0)

    # ACT AND ASSERT
    with pytest.raises(_ResponseFitConvergenceError, match="one wavelength"):
        fit_log_spline(binWave, np.zeros(4), np.full(4, 0.01), PENALTY_KNOT_SPACING, 4.0, 3.0, 5, penalty=1.0)


def test_fit_response_curve_passes_the_smoothing_penalty_to_the_fit() -> None:
    # ARRANGE
    wave, ratio = _ratio_with_wings()
    penalised = {**WING_FIT_SETTINGS, "smoothing_penalty": 10.0, "knot_spacing_nm": 20.0}

    # ACT
    fit = fit_response_curve(wave, ratio, penalised)

    # ASSERT
    assert isinstance(fit.spline, BSpline)
    assert _gap_error(fit) < 0.02


def test_a_linear_algebra_failure_in_the_penalised_least_squares_becomes_a_convergence_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    binWave, binLog, binErr = _penalty_bins([(1000.0, 1400.0)])

    def _singular(*_: object, **__: object) -> None:
        raise np.linalg.LinAlgError("Singular matrix")

    monkeypatch.setattr(np.linalg, "lstsq", _singular)

    # ACT AND ASSERT
    with pytest.raises(_ResponseFitConvergenceError, match="Singular matrix"):
        fit_log_spline(binWave, binLog, binErr, PENALTY_KNOT_SPACING, 4.0, 3.0, 5, penalty=1.0)


def test_the_cap_hit_refit_keeps_the_penalty() -> None:
    # ARRANGE: A BLOCK OF HIGH BINS KEEPS THE CLIP CHANGING PAST THE SINGLE ALLOWED PASS
    binWave, binLog, binErr = _penalty_bins([(1000.0, 1400.0)], noise=0.002)
    binLog[(binWave > 1100.0) & (binWave < 1130.0)] += 0.3
    recordingLog = _RecordingLog()

    # ACT
    spline, kept = fit_log_spline(
        binWave, binLog, binErr, PENALTY_KNOT_SPACING, 4.0, 3.0, 1, log=recordingLog, penalty=1.0
    )

    # ASSERT
    assert isinstance(spline, BSpline)
    assert (~kept).any()
    assert len(recordingLog.warnings) == 1
