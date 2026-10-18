"""Behaviour contracts for the robust log-space spline fit of the instrument response."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from scipy.interpolate import LSQUnivariateSpline

from soxspipe.commonutils.response_fit import (
    _ResponseFitConvergenceError,
    bin_log_response,
    evaluate_log_spline,
    fit_log_spline,
    fit_response_curve,
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
) -> tuple[object, np.ndarray, tuple[float, float], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    mask = mask_response_samples(wave, ratio, telluricRegions or [], [], 2.5)
    binWave, binLog, binErr = bin_log_response(wave, ratio, mask, BIN_WIDTH, 5)
    spline, kept = fit_log_spline(binWave, binLog, binErr, KNOT_SPACING, 4.0, 3.0, maxIterations, log=log)
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


def test_extrapolation_beyond_the_data_is_finite_monotonic_and_close_to_a_log_linear_continuation() -> None:
    # ARRANGE
    wave = np.linspace(400.0, 850.0, 20000)
    ratio = _noisy_ratio(wave, noise=0.03)
    spline, _, fitRange, _ = _fit(wave, ratio)
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
