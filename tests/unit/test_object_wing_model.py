"""Contracts for the object wing model: the object light that a wavelength-only sky fit absorbs."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from soxspipe.commonutils.object_wing_model import (
    core_spectrum,
    estimate_object_wings,
    faint_end_zero_point,
    profile_peak,
    slit_profile,
    wing_model,
)

pytestmark = pytest.mark.unit

BIN_ARCSEC = 0.25


def test_profile_peak_returns_centre_of_brightest_bin():
    # ARRANGE
    slitPosition = np.repeat(np.arange(-5.0, 5.0, BIN_ARCSEC), 30) + 0.01
    residual = 100.0 * np.exp(-0.5 * ((slitPosition - 1.1) / 0.4) ** 2)

    # ACT
    centre = profile_peak(slitPosition, residual)

    # ASSERT
    assert centre == pytest.approx(1.1, abs=BIN_ARCSEC / 2)


def test_profile_peak_ignores_sparse_bins():
    # ARRANGE
    slitPosition = np.repeat(np.arange(-5.0, 5.0, BIN_ARCSEC), 30) + 0.01
    residual = 100.0 * np.exp(-0.5 * ((slitPosition - 1.1) / 0.4) ** 2)
    slitPosition = np.append(slitPosition, -9.0)
    residual = np.append(residual, 1.0e6)

    # ACT
    centre = profile_peak(slitPosition, residual)

    # ASSERT
    assert centre == pytest.approx(1.1, abs=BIN_ARCSEC / 2)


def _smooth_spectrum(wavelength):
    """Return a smooth, positive spectrum with a few percent of ripple."""
    return 1000.0 * (1.0 + 0.3 * np.sin(wavelength / 40.0))


def _core_scene(nRows=2000, seed=1):
    """Build shuffled pixels whose core columns carry the spectrum and whose wing columns carry a tenth of it."""
    rng = np.random.default_rng(seed)
    columns = np.arange(-5.5, 5.51, BIN_ARCSEC)
    rows = np.arange(nRows)
    slitPosition = np.tile(columns, nRows)
    wavelength = 500.0 + 0.05 * np.repeat(rows, len(columns)) + 0.2 * slitPosition
    level = np.where(np.abs(slitPosition - 0.5) < 0.5, 1.0, 0.1)
    residual = _smooth_spectrum(wavelength) * level + rng.normal(0.0, 5.0, wavelength.size)
    order = rng.permutation(wavelength.size)
    return wavelength[order], slitPosition[order], residual[order]


def test_core_spectrum_follows_object_spectrum_at_every_pixel():
    # ARRANGE
    wavelength, slitPosition, residual = _core_scene()
    isUsable = np.ones(wavelength.size, dtype=bool)

    # ACT
    spectrum, nCore = core_spectrum(wavelength, slitPosition, residual, isUsable, 0.5)

    # ASSERT
    assert spectrum.shape == wavelength.shape
    assert nCore > 0
    np.testing.assert_allclose(spectrum, _smooth_spectrum(wavelength), rtol=0.02)


def test_core_spectrum_counts_only_usable_pixels_inside_the_core():
    # ARRANGE
    wavelength, slitPosition, residual = _core_scene()
    isUsable = np.abs(slitPosition) < 3.0
    expectedCore = np.count_nonzero(isUsable & (np.abs(slitPosition - 0.5) < 0.5))

    # ACT
    _, nCore = core_spectrum(wavelength, slitPosition, residual, isUsable, 0.5)

    # ASSERT
    assert nCore == expectedCore


def test_core_spectrum_ignores_unusable_pixels():
    # ARRANGE
    wavelength, slitPosition, residual = _core_scene()
    isUsable = np.ones(wavelength.size, dtype=bool)
    isUsable[::7] = False
    residual = np.where(isUsable, residual, 1.0e9)

    # ACT
    spectrum, _ = core_spectrum(wavelength, slitPosition, residual, isUsable, 0.5)

    # ASSERT
    np.testing.assert_allclose(spectrum, _smooth_spectrum(wavelength), rtol=0.02)


PIXELS_PER_BIN = 40
CORE_LEVEL = 100.0
PIXEL_ERROR = 2.0


def _profile_scene(binValues):
    """Build pixels, 40 per 0.25 arcsec bin, whose residual is the core level times the given profile value."""
    centres = BIN_ARCSEC * np.arange(len(binValues))
    slitPosition = np.repeat(centres, PIXELS_PER_BIN)
    coreSpectrum = np.full(slitPosition.shape, CORE_LEVEL)
    residual = CORE_LEVEL * np.repeat(binValues, PIXELS_PER_BIN)
    error = np.full(slitPosition.shape, PIXEL_ERROR)
    return slitPosition, residual, error, coreSpectrum


def test_slit_profile_is_residual_over_core_spectrum_per_bin():
    # ARRANGE
    binValues = np.array([-0.02, -0.01, 0.0, 0.5, 1.0, 0.4, 0.0, -0.01])
    slitPosition, residual, error, coreSpectrum = _profile_scene(binValues)

    # ACT
    binCentres, profile, _ = slit_profile(slitPosition, residual, error, coreSpectrum)

    # ASSERT
    np.testing.assert_allclose(binCentres, BIN_ARCSEC * np.arange(len(binValues)))
    np.testing.assert_allclose(profile, binValues, atol=1e-12)


def test_slit_profile_noise_is_summed_error_over_summed_core_spectrum():
    # ARRANGE
    slitPosition, residual, error, coreSpectrum = _profile_scene(np.zeros(6))
    expectedNoise = np.sqrt(PIXELS_PER_BIN) * PIXEL_ERROR / (PIXELS_PER_BIN * CORE_LEVEL)

    # ACT
    _, _, profileNoise = slit_profile(slitPosition, residual, error, coreSpectrum)

    # ASSERT
    np.testing.assert_allclose(profileNoise, expectedNoise)


def test_slit_profile_is_nan_where_a_bin_has_too_few_pixels():
    # ARRANGE
    slitPosition, residual, error, coreSpectrum = _profile_scene(np.zeros(6))
    isKept = (slitPosition != 3 * BIN_ARCSEC) | (np.arange(slitPosition.size) % PIXELS_PER_BIN < 5)

    # ACT
    _, profile, profileNoise = slit_profile(slitPosition[isKept], residual[isKept], error[isKept], coreSpectrum[isKept])

    # ASSERT
    assert np.isnan(profile[3]) and np.isnan(profileNoise[3])
    assert np.isfinite(np.delete(profile, 3)).all()


def test_slit_profile_skips_pixels_where_the_core_spectrum_is_not_positive():
    # ARRANGE
    slitPosition, residual, error, coreSpectrum = _profile_scene(np.zeros(6))
    isBroken = np.arange(slitPosition.size) % 2 == 0
    residual = np.where(isBroken, 1.0e9, residual)
    coreSpectrum = np.where(isBroken, 0.0, coreSpectrum)

    # ACT
    _, profile, _ = slit_profile(slitPosition, residual, error, coreSpectrum)

    # ASSERT
    np.testing.assert_allclose(profile, 0.0, atol=1e-12)


NAN = float("nan")


def test_faint_end_zero_point_is_lowest_mean_over_one_arcsec_of_bins():
    # ARRANGE
    profile = np.array([0.5, 0.3, 0.1, -0.1, -0.1, -0.1, -0.1, 0.0, 0.2, 0.6])
    profileNoise = np.full(profile.shape, 0.02)
    binCentres = BIN_ARCSEC * np.arange(profile.size)

    # ACT
    zeroPoint, zeroNoise, runCentre = faint_end_zero_point(binCentres, profile, profileNoise)

    # ASSERT
    assert zeroPoint == pytest.approx(-0.1)
    assert zeroNoise == pytest.approx(np.sqrt(4) * 0.02 / 4)
    assert runCentre == pytest.approx(binCentres[3:7].mean())


def test_faint_end_zero_point_skips_runs_that_contain_a_nan_bin():
    # ARRANGE
    profile = np.array([5.0, -1.0, -1.0, NAN, -1.0, -1.0, 5.0, 0.1, 0.1, 0.1, 0.1, 0.2])
    profileNoise = np.full(profile.shape, 0.02)
    binCentres = BIN_ARCSEC * np.arange(profile.size)

    # ACT
    zeroPoint, _, runCentre = faint_end_zero_point(binCentres, profile, profileNoise)

    # ASSERT
    assert zeroPoint == pytest.approx(0.1)
    assert runCentre == pytest.approx(binCentres[7:11].mean())


def test_faint_end_zero_point_is_nan_without_a_complete_run_of_finite_bins():
    # ARRANGE
    profile = np.array([-1.0, -1.0, -1.0, NAN, -1.0, -1.0, -1.0, NAN, -1.0])
    profileNoise = np.full(profile.shape, 0.02)
    binCentres = BIN_ARCSEC * np.arange(profile.size)

    # ACT
    zeroPoint, zeroNoise, runCentre = faint_end_zero_point(binCentres, profile, profileNoise)

    # ASSERT
    assert np.isnan([zeroPoint, zeroNoise, runCentre]).all()


N_ROWS = 2000
FIRST_WAVELENGTH = 500.0
WAVELENGTH_PER_ROW = 0.05
WAVELENGTH_PER_ARCSEC = 0.2
SLIT_HALF_LENGTH = 5.5
SLIT_DRIFT_PER_ROW = 0.0137
MOFFAT_BETA = 4.765
MOFFAT_FWHM = 1.0
MOFFAT_ALPHA = MOFFAT_FWHM / (2.0 * np.sqrt(2.0 ** (1.0 / MOFFAT_BETA) - 1.0))
WING_SCALE_ARCSEC = 3.0
UNMASKED_EDGE_ARCSEC = 2.5
SKY_LEVEL = 50.0
BRIGHT_PEAK = 1.0e5
GATE_FRACTION = 0.02


def _object_spectrum(wavelength, peak):
    """Return the injected object spectrum: a smooth ripple around the peak flux."""
    return peak * (1.0 + 0.3 * np.sin(wavelength / 40.0))


def _object_profile(slitPosition, wingAmplitude):
    """Return the injected spatial profile: a unit-peak Moffat plus a one-sided exponential on the +slit side."""
    moffat = (1.0 + (slitPosition / MOFFAT_ALPHA) ** 2) ** -MOFFAT_BETA
    oneSided = wingAmplitude * np.exp(-slitPosition / WING_SCALE_ARCSEC) * (slitPosition > 0)
    return moffat + oneSided


@dataclass(frozen=True)
class synthetic_order:
    """One shuffled order with the truth that was injected into it."""

    wavelength: np.ndarray
    slitPosition: np.ndarray
    flux: np.ndarray
    skyModel: np.ndarray
    error: np.ndarray
    isUsable: np.ndarray
    spectrum: np.ndarray
    profile: np.ndarray
    unmaskedLevel: float
    injectedBias: float

    def estimate(self, isUsable=None, gateFraction=GATE_FRACTION):
        """Run the estimate on this order, optionally with another usable mask."""
        isUsable = self.isUsable if isUsable is None else isUsable
        return estimate_object_wings(
            self.wavelength, self.slitPosition, self.flux, self.skyModel, self.error, isUsable, gateFraction
        )


def _build_order(peak=BRIGHT_PEAK, biasFraction=1.0, extraNoise=0.0, columnStep=BIN_ARCSEC, wingAmplitude=0.0, seed=7):
    """Build a shuffled order: object plus flat sky plus noise, and a sky model that carries a share of the wings.

    The sky model is the sky plus ``biasFraction`` times the object spectrum times the true mean profile level
    outside +-2.5 arcsec, which is what a sky fit without an object mask would absorb.
    """
    rng = np.random.default_rng(seed)
    columns = np.arange(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH + columnStep / 2, columnStep)
    rows = np.repeat(np.arange(N_ROWS), len(columns))
    drift = (rows * SLIT_DRIFT_PER_ROW) % 1.0 - 0.5
    slitPosition = np.tile(columns, N_ROWS) + columnStep * drift
    wavelength = FIRST_WAVELENGTH + WAVELENGTH_PER_ROW * rows + WAVELENGTH_PER_ARCSEC * slitPosition
    spectrum = _object_spectrum(wavelength, peak)
    profile = _object_profile(slitPosition, wingAmplitude)
    unmaskedLevel = profile[np.abs(slitPosition) > UNMASKED_EDGE_ARCSEC].mean()
    objectLight = spectrum * profile
    sigma = np.sqrt(np.maximum(SKY_LEVEL + objectLight, 1.0) + extraNoise**2)
    flux = SKY_LEVEL + objectLight + rng.normal(0.0, 1.0, wavelength.size) * sigma
    skyModel = SKY_LEVEL + biasFraction * unmaskedLevel * spectrum
    shuffle = rng.permutation(wavelength.size)
    return synthetic_order(
        wavelength[shuffle],
        slitPosition[shuffle],
        flux[shuffle],
        skyModel[shuffle],
        sigma[shuffle],
        np.ones(wavelength.size, dtype=bool),
        spectrum[shuffle],
        profile[shuffle],
        float(unmaskedLevel),
        float(biasFraction * unmaskedLevel * np.median(spectrum)),
    )


def test_estimate_recovers_the_injected_bias_of_a_bright_star():
    # ARRANGE
    order = _build_order()

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    assert estimate.reason == ""
    assert estimate.biasLevel == pytest.approx(order.injectedBias, rel=0.05)
    assert estimate.centre == pytest.approx(0.0, abs=BIN_ARCSEC)


FAINT_END_ARCSEC = 4.6


def test_wing_model_matches_the_injected_wings_above_the_faint_end():
    # ARRANGE
    order = _build_order()
    estimate = order.estimate()
    faintLevel = order.profile[np.abs(order.slitPosition) > FAINT_END_ARCSEC].mean()
    truth = order.spectrum * (order.profile - faintLevel)
    isWing = (np.abs(order.slitPosition) > 2.0) & (np.abs(order.slitPosition) < 3.5)
    isOuter = np.abs(order.slitPosition) > 1.5

    # ACT
    wings = wing_model(estimate, order.slitPosition)

    # ASSERT
    # LINEAR INTERPOLATION OF 0.25 ARCSEC BIN MEANS OVERSHOOTS A CONVEX MOFFAT WING BY ABOUT H^2/8 TIMES ITS CURVATURE
    # (5 TO 15 PERCENT AT 2 ARCSEC), SO THE BANDS ARE WIDER THAN THE NOISE
    np.testing.assert_allclose(wings[isWing], truth[isWing], rtol=0.2)
    assert wings[isOuter].sum() == pytest.approx(truth[isOuter].sum(), rel=0.15)


WING_AMPLITUDE = 0.02


def _bin_of_each_pixel(order, estimate):
    """Return the index of the estimate's profile bin that each pixel of the order falls in."""
    steps = (order.slitPosition - estimate.binCentres[0]) / BIN_ARCSEC
    return np.clip(np.rint(steps).astype(int), 0, len(estimate.binCentres) - 1)


def test_asymmetric_wings_put_the_zero_point_on_the_faint_side_and_the_profile_follows_them():
    # ARRANGE
    order = _build_order(wingAmplitude=WING_AMPLITUDE)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    assert estimate.biasLevel == pytest.approx(order.injectedBias, rel=0.05)
    _, _, runCentre = faint_end_zero_point(estimate.binCentres, estimate.profile, np.zeros(estimate.profile.shape))
    assert runCentre < -FAINT_END_ARCSEC
    binIndex = _bin_of_each_pixel(order, estimate)
    binMean = np.bincount(binIndex, weights=order.profile) / np.bincount(binIndex)
    faintLevel = order.profile[order.slitPosition < -FAINT_END_ARCSEC].mean()
    truth = order.spectrum * (binMean[binIndex] - faintLevel)
    wings = wing_model(estimate, estimate.binCentres[binIndex])
    isWing = (truth > 20.0) & (np.abs(order.slitPosition) > 1.0)
    assert np.count_nonzero(isWing & (order.slitPosition > 0)) > 0
    ratio = wings[isWing] / truth[isWing]
    # THE CORE SPECTRUM IS A ROLLING MEDIAN OVER SEVERAL COLUMNS, SO IT WOBBLES BY ABOUT 5 PERCENT AS THE SLIT DRIFTS
    assert np.median(ratio) == pytest.approx(1.0, abs=0.02)
    assert np.abs(ratio - 1.0).max() < 0.15


FAINT_PEAK = 5000.0
NOISY_EXTRA_SIGMA = 1500.0
HALF_BINNED_STEP = 2 * BIN_ARCSEC


def _bin_index_of(order):
    """Return the 0.25 arcsec bin number of each pixel, counted from the lowest slit position."""
    return np.floor((order.slitPosition - order.slitPosition.min()) / BIN_ARCSEC + 0.5).astype(int)


def test_sky_dominated_faint_star_is_not_corrected_and_the_reason_names_the_gate_fraction():
    # ARRANGE
    order = _build_order(peak=FAINT_PEAK)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert not estimate.isCorrected
    assert "gate fraction" in estimate.reason


def test_large_but_noisy_bias_is_not_corrected_and_the_reason_names_the_significance():
    # ARRANGE
    order = _build_order(extraNoise=NOISY_EXTRA_SIGMA)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert not estimate.isCorrected
    assert "significance" in estimate.reason
    assert estimate.biasLevel > GATE_FRACTION * estimate.skyLevel


def test_sky_model_without_object_light_is_not_corrected():
    # ARRANGE
    order = _build_order(biasFraction=-1.0)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert not estimate.isCorrected
    assert estimate.zeroPoint > 0
    assert "not negative" in estimate.reason


def test_an_empty_profile_bin_does_not_stop_the_zero_point_being_found_elsewhere():
    # ARRANGE
    order = _build_order()
    isUsable = _bin_index_of(order) != 2

    # ACT
    estimate = order.estimate(isUsable)

    # ASSERT
    assert np.isnan(estimate.profile[2])
    assert estimate.isCorrected, estimate.reason
    assert estimate.biasLevel == pytest.approx(order.injectedBias, rel=0.05)


def test_slit_without_a_complete_run_of_finite_bins_is_not_corrected():
    # ARRANGE
    order = _build_order()
    isUsable = _bin_index_of(order) % 4 != 2

    # ACT
    estimate = order.estimate(isUsable)

    # ASSERT
    assert not estimate.isCorrected
    assert "run" in estimate.reason


def test_half_as_many_slit_columns_still_give_a_finite_profile_and_the_bias():
    # ARRANGE
    order = _build_order(columnStep=HALF_BINNED_STEP)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    assert np.isfinite(estimate.profile).all()
    assert estimate.biasLevel == pytest.approx(order.injectedBias, rel=0.05)


def test_slit_loss_fraction_is_the_share_of_the_object_light_the_sky_fit_removes():
    # ARRANGE
    order = _build_order()
    binCentres = np.arange(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH + BIN_ARCSEC / 2, BIN_ARCSEC)
    binLight = _object_profile(binCentres, 0.0)
    expected = len(binCentres) * order.unmaskedLevel / binLight.sum()

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.slitLossFraction == pytest.approx(expected, rel=0.1)


def test_estimate_reports_too_few_core_pixels():
    # ARRANGE
    order = _build_order()
    isUsable = order.wavelength < FIRST_WAVELENGTH + 2.0

    # ACT
    estimate = order.estimate(isUsable)

    # ASSERT
    assert not estimate.isCorrected
    assert "core pixels" in estimate.reason


def test_estimate_without_usable_pixels_is_not_corrected():
    # ARRANGE
    order = _build_order(peak=FAINT_PEAK)

    # ACT
    estimate = order.estimate(np.zeros(order.flux.shape, dtype=bool))

    # ASSERT
    assert not estimate.isCorrected
    assert estimate.reason != ""


def test_estimate_rejects_arrays_of_different_shapes():
    # ARRANGE
    order = _build_order(peak=FAINT_PEAK)

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="same shape"):
        estimate_object_wings(
            order.wavelength[:-1],
            order.slitPosition,
            order.flux,
            order.skyModel,
            order.error,
            order.isUsable,
            GATE_FRACTION,
        )


def test_wing_model_of_an_uncorrected_estimate_is_all_zeros():
    # ARRANGE
    order = _build_order(peak=FAINT_PEAK)
    estimate = order.estimate()

    # ACT
    wings = wing_model(estimate, order.slitPosition)

    # ASSERT
    assert not estimate.isCorrected
    assert wings.shape == order.slitPosition.shape
    assert not wings.any()


def test_wing_model_never_returns_a_value_that_is_not_finite():
    # ARRANGE
    order = _build_order()
    wavelength = order.wavelength.copy()
    wavelength[0] = np.nan
    isUsable = np.isfinite(wavelength)
    estimate = estimate_object_wings(
        wavelength, order.slitPosition, order.flux, order.skyModel, order.error, isUsable, GATE_FRACTION
    )

    # ACT
    wings = wing_model(estimate, order.slitPosition)

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    assert np.isnan(estimate.coreSpectrum[0])
    assert np.isfinite(wings).all()
    assert wings[0] == 0.0


def test_wing_model_rejects_slit_positions_of_a_different_length():
    # ARRANGE
    order = _build_order(peak=FAINT_PEAK)
    estimate = order.estimate()

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="slit positions"):
        wing_model(estimate, order.slitPosition[:-1])


def test_profile_peak_is_nan_when_no_bin_has_enough_pixels():
    # ARRANGE
    slitPosition = np.linspace(-5.0, 5.0, 10)

    # ACT
    centre = profile_peak(slitPosition, np.ones(10))

    # ASSERT
    assert np.isnan(centre)


def test_core_spectrum_is_nan_when_no_pixel_is_near_the_centre():
    # ARRANGE
    wavelength, slitPosition, residual = _core_scene(nRows=100)
    isUsable = np.ones(wavelength.size, dtype=bool)

    # ACT
    spectrum, nCore = core_spectrum(wavelength, slitPosition, residual, isUsable, 50.0)

    # ASSERT
    assert nCore == 0
    assert np.isnan(spectrum).all()


def test_faint_end_zero_point_is_nan_for_a_profile_shorter_than_one_run():
    # ARRANGE
    profile = np.array([-0.1, -0.1, -0.1])

    # ACT
    zeroPoint, _, _ = faint_end_zero_point(BIN_ARCSEC * np.arange(3), profile, np.full(3, 0.02))

    # ASSERT
    assert np.isnan(zeroPoint)


def test_estimate_reports_that_no_slit_bin_has_enough_usable_pixels():
    # ARRANGE
    order = _build_order(peak=FAINT_PEAK)
    isUsable = np.zeros(order.flux.shape, dtype=bool)
    isUsable[:10] = True

    # ACT
    estimate = order.estimate(isUsable)

    # ASSERT
    assert not estimate.isCorrected
    assert "enough usable pixels" in estimate.reason
