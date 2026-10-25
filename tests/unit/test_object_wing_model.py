"""Contracts for the object wing model: the object light that a wavelength-only sky fit absorbs."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pytest

from soxspipe.commonutils.object_wing_model import (
    core_spectrum,
    estimate_object_wings,
    far_slit_zero_point,
    profile_peak,
    profile_sources,
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


def test_far_slit_zero_point_is_the_mean_profile_beyond_three_arcsec_from_the_centre():
    # ARRANGE
    binCentres = BIN_ARCSEC * np.arange(-24, 25)
    profile = np.where(np.abs(binCentres - 0.5) > 3.0, -0.01, 1.0)
    profile[binCentres < -4.0] = -0.03
    profileNoise = np.full(profile.shape, 0.02)
    isFar = np.abs(binCentres - 0.5) > 3.0

    # ACT
    zeroPoint, zeroNoise = far_slit_zero_point(binCentres, profile, profileNoise, 0.5)

    # ASSERT
    # 8 BINS FROM -6.00 TO -4.25 AT -0.03, 6 MORE LEFT BINS DOWN TO -2.75 AND 10 RIGHT BINS FROM 3.75 AT -0.01
    assert isFar.sum() == 24
    assert zeroPoint == pytest.approx((8 * -0.03 + 16 * -0.01) / 24)
    assert zeroNoise == pytest.approx(np.sqrt(24) * 0.02 / 24)


def test_far_slit_zero_point_skips_bins_that_are_not_finite():
    # ARRANGE
    binCentres = BIN_ARCSEC * np.arange(-24, 25)
    profile = np.where(np.abs(binCentres) > 3.0, -0.01, 1.0)
    profile[0] = NAN
    profileNoise = np.full(profile.shape, 0.02)

    # ACT
    zeroPoint, _ = far_slit_zero_point(binCentres, profile, profileNoise, 0.0)

    # ASSERT
    assert zeroPoint == pytest.approx(-0.01)


def test_far_slit_zero_point_is_nan_with_fewer_than_one_arcsec_of_far_bins():
    # ARRANGE
    binCentres = BIN_ARCSEC * np.arange(16)
    profile = np.full(binCentres.shape, -0.01)
    profileNoise = np.full(profile.shape, 0.02)

    # ACT
    zeroPoint, zeroNoise = far_slit_zero_point(binCentres, profile, profileNoise, 0.0)

    # ASSERT
    # ONLY THE BINS AT 3.25, 3.50 AND 3.75 ARE MORE THAN 3 ARCSEC FROM THE CENTRE
    assert np.count_nonzero(np.abs(binCentres) > 3.0) == 3
    assert np.isnan(zeroPoint)
    assert np.isnan(zeroNoise)


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
# THE ESTIMATE TAKES THE PROFILE BINS CENTRED MORE THAN 3 ARCSEC FROM THE OBJECT AS SKY (FAR_SLIT_ARCSEC). THE OBJECT
# CENTRE IS ITSELF A BIN CENTRE, SO THE NEAREST FAR BIN STARTS HALF A BIN FURTHER OUT
FAR_ARCSEC = 3.0 + BIN_ARCSEC / 2
SKY_LEVEL = 50.0
BRIGHT_PEAK = 1.0e5
MIN_LEAK_SIGMA = 0.03


def _object_spectrum(wavelength, peak):
    """Return the injected object spectrum: a smooth ripple around the peak flux."""
    return peak * (1.0 + 0.3 * np.sin(wavelength / 40.0))


def _object_profile(slitPosition, wingAmplitude):
    """Return the injected spatial profile: a unit-peak Moffat plus a one-sided exponential on the +slit side."""
    moffat = (1.0 + (slitPosition / MOFFAT_ALPHA) ** 2) ** -MOFFAT_BETA
    oneSided = wingAmplitude * np.exp(-slitPosition / WING_SCALE_ARCSEC) * (slitPosition > 0)
    return moffat + oneSided


@dataclass(frozen=True, eq=False)
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
    farLevel: float
    visibleBias: float

    def estimate(self, isUsable=None, minLeakSigma=MIN_LEAK_SIGMA):
        """Run the estimate on this order, optionally with another usable mask."""
        isUsable = self.isUsable if isUsable is None else isUsable
        return estimate_object_wings(
            self.wavelength, self.slitPosition, self.flux, self.skyModel, self.error, isUsable, minLeakSigma
        )


def _build_order(peak=BRIGHT_PEAK, biasFraction=1.0, extraNoise=0.0, columnStep=BIN_ARCSEC, wingAmplitude=0.0, seed=7):
    """Build a shuffled order: object plus flat sky plus noise, and a sky model that carries a share of the wings.

    The sky model is the sky plus ``biasFraction`` times the object spectrum times the true mean profile level
    outside +-2.5 arcsec, which is what a sky fit without an object mask would absorb.

    The estimate takes the profile more than 3 arcsec from the object as sky, so the object light out there is
    invisible to it. ``farLevel`` is the true mean profile there, and ``visibleBias`` is the injected bias less that
    light: the bias the estimate should find.
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
    farLevel = profile[np.abs(slitPosition) > FAR_ARCSEC].mean()
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
        float(farLevel),
        float((biasFraction * unmaskedLevel - farLevel) * np.median(spectrum)),
    )


def test_estimate_recovers_the_injected_bias_of_a_bright_star_less_its_light_beyond_three_arcsec():
    # ARRANGE
    order = _build_order()

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    assert estimate.reason == ""
    assert estimate.biasLevel == pytest.approx(order.visibleBias, rel=0.05)
    # THE MOFFAT WINGS PUT ABOUT A QUARTER OF THE INJECTED BIAS BEYOND 3 ARCSEC, WHERE THE ESTIMATE CANNOT SEE IT
    assert order.visibleBias < 0.85 * order.injectedBias
    assert estimate.centre == pytest.approx(0.0, abs=BIN_ARCSEC)


def test_wing_model_matches_the_injected_wings_above_the_far_slit_level():
    # ARRANGE
    order = _build_order()
    estimate = order.estimate()
    truth = order.spectrum * (order.profile - order.farLevel)
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


def test_one_sided_wing_light_beyond_three_arcsec_is_taken_as_sky_and_the_profile_follows_the_wings():
    # ARRANGE
    order = _build_order(wingAmplitude=WING_AMPLITUDE)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    binIndex = _bin_of_each_pixel(order, estimate)
    binMean = np.bincount(binIndex, weights=order.profile) / np.bincount(binIndex)
    # A ONE-SIDED WING MAKES THE FAR LEVEL DEPEND ON HOW MANY FAR BINS EACH SIDE HAS, SO IT IS TAKEN FROM THE BINS
    farLevel = binMean[np.abs(estimate.binCentres - estimate.centre) > 3.0].mean()
    # THE FAR SLIT ON THE WING SIDE HOLDS OBJECT LIGHT, SO THE ESTIMATE FINDS LESS THAN THE INJECTED BIAS (DY-1964)
    expectedBias = (order.unmaskedLevel - farLevel) * np.median(order.spectrum)
    assert estimate.biasLevel == pytest.approx(expectedBias, rel=0.05)
    assert expectedBias < 0.5 * order.injectedBias
    truth = order.spectrum * (binMean[binIndex] - farLevel)
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


TINY_LEAK_PEAK = 2000.0
FEW_PERCENT_LEAK_PEAK = 4.0e4


def test_leak_small_against_the_noise_of_a_faint_sky_pixel_is_not_corrected():
    # ARRANGE
    order = _build_order(peak=TINY_LEAK_PEAK)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert not estimate.isCorrected
    assert "noise of a faint-sky pixel" in estimate.reason
    assert estimate.biasLevel < MIN_LEAK_SIGMA * np.median(order.error)


def test_leak_of_a_few_percent_of_the_sky_is_corrected_when_it_stands_above_the_faint_sky_noise():
    # ARRANGE
    order = _build_order(peak=FEW_PERCENT_LEAK_PEAK)

    # ACT
    estimate = order.estimate()

    # ASSERT
    # DY-1964: NIR ORDERS 18-23 LEAKED 4-7 PERCENT OF THE SKY AND THE OLD 10 PERCENT GATE SKIPPED THEM
    assert estimate.isCorrected, estimate.reason
    assert estimate.biasLevel < 0.07 * estimate.skyLevel
    assert estimate.biasLevel == pytest.approx(order.visibleBias, rel=0.1)


def test_large_but_noisy_bias_is_not_corrected_and_the_reason_names_the_significance():
    # ARRANGE
    order = _build_order(extraNoise=NOISY_EXTRA_SIGMA)

    # ACT
    # THE BIAS IS ALSO FAR BELOW THE NOISE OF ONE PIXEL, SO THAT GATE IS OPENED TO REACH THE SIGNIFICANCE TEST
    estimate = order.estimate(minLeakSigma=0.0)

    # ASSERT
    assert not estimate.isCorrected
    assert "significance" in estimate.reason


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
    assert estimate.biasLevel == pytest.approx(order.visibleBias, rel=0.05)


def test_slit_with_too_little_profile_far_from_the_object_is_not_corrected():
    # ARRANGE
    order = _build_order()
    isUsable = np.abs(order.slitPosition) < 3.2

    # ACT
    estimate = order.estimate(isUsable)

    # ASSERT
    assert not estimate.isCorrected
    assert "arcsec from the object" in estimate.reason


SKY_LINE_PEAK = 2000.0
SKY_LINE_SHARE = 0.2
DIM_END_ARCSEC = -4.5
DIM_END_DEFICIT = 0.03


def _with_sky_lines_and_a_dim_slit_end(order, seed=11):
    """Return the order with bright sky lines on a share of the wavelengths and one slit end that is a little dim.

    The dim end lacks ``DIM_END_DEFICIT`` of the sky beyond ``DIM_END_ARCSEC``, as an illumination error would, so its
    residual scales with the sky, not with the object. The sky model holds the lines but not the dim end.
    """
    rng = np.random.default_rng(seed)
    isLine = (np.sin(order.wavelength * 7.0) > np.cos(np.pi * SKY_LINE_SHARE)).astype(float)
    sky = SKY_LEVEL + SKY_LINE_PEAK * isLine
    isDim = order.slitPosition < DIM_END_ARCSEC
    objectLight = order.spectrum * order.profile
    sigma = np.sqrt(np.maximum(sky + objectLight, 1.0))
    flux = sky * (1.0 - DIM_END_DEFICIT * isDim) + objectLight + rng.normal(0.0, 1.0, sigma.size) * sigma
    skyModel = order.skyModel + SKY_LINE_PEAK * isLine
    return replace(order, flux=flux, skyModel=skyModel, error=sigma)


def test_a_dim_slit_end_under_bright_sky_lines_does_not_inflate_the_bias():
    # ARRANGE
    order = _with_sky_lines_and_a_dim_slit_end(_build_order())

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    # THE FAR BINS AVERAGE THE DIM END WITH THE REST OF THE FAR SLIT, AND ITS SKY-LINE PIXELS ARE LEFT OUT
    assert estimate.biasLevel == pytest.approx(order.visibleBias, rel=0.1)


def test_a_first_fit_that_used_only_bright_sky_pixels_is_not_corrected():
    # ARRANGE
    order = _with_sky_lines_and_a_dim_slit_end(_build_order())
    isFitted = order.skyModel > np.median(order.skyModel)

    # ACT
    estimate = estimate_object_wings(
        order.wavelength,
        order.slitPosition,
        order.flux,
        order.skyModel,
        order.error,
        order.isUsable,
        MIN_LEAK_SIGMA,
        isFitted=isFitted,
    )

    # ASSERT
    assert not estimate.isCorrected
    assert "faint sky" in estimate.reason


def test_faint_sky_noise_comes_from_the_far_slit_not_the_object_core():
    # ARRANGE
    order = _build_order()
    farNoise = np.median(order.error[np.abs(order.slitPosition) > FAR_ARCSEC])

    # ACT
    estimate = order.estimate()

    # ASSERT
    # THE BRIGHT CORE PIXELS WOULD RAISE THE MEDIAN ERROR BY ABOUT 20 PERCENT AND MAKE THE LEAK GATE STRICTER
    assert np.median(order.error) > 1.1 * farNoise
    assert estimate.faintSkyNoise == pytest.approx(farNoise, rel=0.02)


def test_half_as_many_slit_columns_still_give_a_finite_profile_and_the_bias():
    # ARRANGE
    order = _build_order(columnStep=HALF_BINNED_STEP)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    assert np.isfinite(estimate.profile).all()
    assert estimate.biasLevel == pytest.approx(order.visibleBias, rel=0.05)


def test_slit_loss_fraction_is_the_share_of_the_object_light_the_sky_fit_removes():
    # ARRANGE
    order = _build_order()
    binCentres = np.arange(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH + BIN_ARCSEC / 2, BIN_ARCSEC)
    binLight = _object_profile(binCentres, 0.0)
    expected = len(binCentres) * (order.unmaskedLevel - order.farLevel) / (binLight - order.farLevel).sum()

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
            MIN_LEAK_SIGMA,
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
        wavelength, order.slitPosition, order.flux, order.skyModel, order.error, isUsable, MIN_LEAK_SIGMA
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


OBJECT_MASK_ARCSEC = 2.0
UPPER_CLIP_SIGMA = 2.0
N_CLIP_PASSES = 5


def _object_share(order, isFitted, flux):
    """Return the least-squares share of the object spectrum in the flux of the fitted pixels, above the sky."""
    spectrum = order.spectrum[isFitted]
    return float(np.sum((flux[isFitted] - SKY_LEVEL) * spectrum) / np.sum(spectrum**2))


def _upper_clipped_first_fit(order):
    """Mimic a first sky fit with an object mask and an upper-only clip.

    Return the mask of the pixels the fit used, the share of the object spectrum the fit absorbed, and the sky model.
    """
    isOutsideMask = np.abs(order.slitPosition) > OBJECT_MASK_ARCSEC
    isFitted = isOutsideMask
    for _ in range(N_CLIP_PASSES):
        share = _object_share(order, isFitted, order.flux)
        skyModel = SKY_LEVEL + share * order.spectrum
        isFitted = isOutsideMask & (order.flux - skyModel < UPPER_CLIP_SIGMA * order.error)
    share = _object_share(order, isFitted, order.flux)
    return isFitted, share, SKY_LEVEL + share * order.spectrum


def _estimate_with(order, skyModel, isFitted=None):
    """Run the estimate on the order against a given sky model, optionally with the mask of the fitted pixels.

    The leak gate is open: the one-sided wing leaves a bias well below the noise of one pixel, and these tests are
    about which pixels build the profile, not about the gate.
    """
    kwargs = {} if isFitted is None else {"isFitted": isFitted}
    return estimate_object_wings(
        order.wavelength, order.slitPosition, order.flux, skyModel, order.error, order.isUsable, 0.0, **kwargs
    )


def test_profile_from_the_fitted_pixels_does_not_overstate_the_wings_an_upper_clipped_fit_saw():
    # ARRANGE
    order = _build_order(wingAmplitude=WING_AMPLITUDE)
    isFitted, share, skyModel = _upper_clipped_first_fit(order)
    assert share > 0

    # ACT
    allPixelWings = wing_model(_estimate_with(order, skyModel), order.slitPosition)
    fittedEstimate = _estimate_with(order, skyModel, isFitted)
    fittedWings = wing_model(fittedEstimate, order.slitPosition)

    # ASSERT
    # THE FITTED-PIXEL PROFILE REMOVES PART OF THE OBJECT LIGHT AND NEVER MORE THAN THE FIT SAW, SO NO DEFICIT IS LEFT.
    # THE REST IS THE WING LIGHT BEYOND 3 ARCSEC, WHICH THE ESTIMATE TAKES AS SKY (DY-1964). THE ALL-PIXEL PROFILE
    # COUNTS THE CLIPPED WING PIXELS BEYOND 3 ARCSEC IN ITS ZERO POINT AND REMOVES NOTHING
    assert fittedEstimate.isCorrected, fittedEstimate.reason
    fittedLeft = _object_share(order, isFitted, order.flux - fittedWings)
    allPixelLeft = _object_share(order, isFitted, order.flux - allPixelWings)
    assert 0.0 <= fittedLeft < 0.9 * share
    assert allPixelLeft == pytest.approx(share)


def test_slit_loss_fraction_still_counts_the_object_core_when_the_fit_did_not_use_it():
    # ARRANGE
    order = _build_order()
    isFitted = np.abs(order.slitPosition) > OBJECT_MASK_ARCSEC
    binCentres = np.arange(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH + BIN_ARCSEC / 2, BIN_ARCSEC)
    binLight = _object_profile(binCentres, 0.0)
    expected = len(binCentres) * (order.unmaskedLevel - order.farLevel) / (binLight - order.farLevel).sum()

    # ACT
    estimate = _estimate_with(order, order.skyModel, isFitted)

    # ASSERT
    assert estimate.slitLossFraction == pytest.approx(expected, rel=0.1)


def test_estimate_without_any_pixel_the_first_fit_used_is_not_corrected():
    # ARRANGE
    order = _build_order()

    # ACT
    estimate = _estimate_with(order, order.skyModel, np.zeros(order.flux.shape, dtype=bool))

    # ASSERT
    assert not estimate.isCorrected
    assert "first sky fit" in estimate.reason


COMPANION_OFFSET_ARCSEC = 3.0
# THE COMPANION'S OWN LIGHT FILLS THE FAR SLIT, SO A CROWDED ORDER ONLY REACHES THE CROWDING TEST WHEN ITS SKY MODEL
# HOLDS EVEN MORE OBJECT LIGHT THAN THAT (DY-1964)
CROWDED_BIAS_FRACTION = 2000.0


def _source_profile(binCentres, sources, noise=0.001, seed=3):
    """Return a binned profile of unit-width Gaussian sources, (centre, height) pairs, with Gaussian noise."""
    rng = np.random.default_rng(seed)
    profile = np.zeros(binCentres.shape)
    for centre, height in sources:
        profile += height * np.exp(-0.5 * ((binCentres - centre) / 0.45) ** 2)
    profileNoise = np.full(binCentres.shape, noise)
    return profile + rng.normal(0.0, noise, binCentres.shape), profileNoise


SLIT_BINS = np.arange(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH + BIN_ARCSEC / 2, BIN_ARCSEC)


def test_profile_sources_finds_one_source_in_a_single_star_profile():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(0.5, 1.0)])

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [0.5], atol=BIN_ARCSEC)


def test_profile_sources_finds_every_bright_source_on_a_crowded_slit_brightest_first():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(-3.5, 0.6), (0.0, 1.0), (3.5, 0.4)])

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [0.0, -3.5, 3.5], atol=BIN_ARCSEC)


def test_profile_sources_ignores_a_faint_bump_that_is_a_small_fraction_of_the_profile():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(0.0, 1.0), (-3.0, 0.08)])

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    assert len(sources) == 1


def test_profile_sources_ignores_a_bump_that_is_not_significant_against_its_noise():
    # ARRANGE
    profile, _ = _source_profile(SLIT_BINS, [(0.0, 1.0), (3.0, 0.5)])
    profileNoise = np.full(SLIT_BINS.shape, 0.2)

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    assert len(sources) == 1


def test_profile_sources_merges_peaks_closer_than_the_minimum_separation():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(0.0, 1.0), (1.0, 0.9)], noise=0.0001)

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    assert len(sources) == 1


def test_profile_sources_skips_empty_bins():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(-3.0, 0.7), (2.0, 1.0)])
    profile[10] = np.nan
    profileNoise[10] = np.nan

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [2.0, -3.0], atol=BIN_ARCSEC)


def test_profile_sources_is_empty_without_finite_bins():
    # ARRANGE
    profile = np.full(SLIT_BINS.shape, np.nan)

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profile.copy())

    # ASSERT
    assert sources.size == 0


def _with_companion(order, fraction):
    """Return the order with a second star, ``fraction`` as bright as the first, added to its flux."""
    companion = fraction * order.spectrum * _object_profile(order.slitPosition - COMPANION_OFFSET_ARCSEC, 0.0)
    return replace(order, flux=order.flux + companion)


def test_bright_star_with_a_companion_on_the_slit_is_not_corrected_and_the_reason_names_the_sources():
    # ARRANGE
    order = _with_companion(_build_order(biasFraction=CROWDED_BIAS_FRACTION), 0.5)

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert not estimate.isCorrected
    assert estimate.isCrowdedSlit
    assert "2 sources" in estimate.reason
    np.testing.assert_allclose(np.sort(estimate.sourceCentres), [0.0, COMPANION_OFFSET_ARCSEC], atol=BIN_ARCSEC)
    assert estimate.biasLevel > 0


def test_corrected_bright_star_reports_its_one_source():
    # ARRANGE
    order = _build_order()

    # ACT
    estimate = order.estimate()

    # ASSERT
    assert estimate.isCorrected, estimate.reason
    assert not estimate.isCrowdedSlit
    np.testing.assert_allclose(estimate.sourceCentres, [0.0], atol=BIN_ARCSEC)


def test_crowded_slit_that_another_gate_already_stops_is_not_flagged_as_crowded():
    # ARRANGE
    order = _with_companion(_build_order(biasFraction=CROWDED_BIAS_FRACTION), 0.5)

    # ACT
    estimate = order.estimate(minLeakSigma=1.0e6)

    # ASSERT
    assert not estimate.isCorrected
    assert not estimate.isCrowdedSlit
    assert "noise of a faint-sky pixel" in estimate.reason
    assert len(estimate.sourceCentres) == 2


def test_profile_sources_still_finds_a_companion_when_an_unlit_slit_end_sits_far_below_the_profile():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(-3.5, 1.0), (0.2, 0.6)])
    profile[:2] = [-3.3, -1.8]

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [-3.5, 0.2], atol=BIN_ARCSEC)


def test_profile_sources_bridges_empty_bins_at_the_peak_of_a_source():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(-3.0, 1.0), (1.5, 0.8)])
    gap = np.abs(SLIT_BINS - 1.5) < 0.3
    profile[gap] = np.nan
    profileNoise[gap] = np.nan

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    assert len(sources) == 2
    assert sources[0] == pytest.approx(-3.0, abs=BIN_ARCSEC)
    # A BRIDGED PEAK CAN LAND ANYWHERE FROM ONE VALID FLANK OF THE GAP TO THE OTHER
    assert abs(sources[1] - 1.5) <= 0.3 + BIN_ARCSEC


def test_profile_sources_does_not_make_a_peak_where_a_falling_profile_has_a_gap():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(-3.0, 1.0)])
    profile = profile - 0.05 * SLIT_BINS
    gap = np.abs(SLIT_BINS - 1.5) < 0.6
    profile[gap] = np.nan
    profileNoise[gap] = np.nan

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [-3.0], atol=BIN_ARCSEC)


def test_profile_sources_skips_bins_without_a_positive_noise():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(0.0, 1.0), (3.0, 0.5)])
    profileNoise[SLIT_BINS > 2.0] = 0.0

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [0.0], atol=BIN_ARCSEC)


def _invalidate(profile, profileNoise, isInvalid, how):
    """Return copies of the profile and noise with the chosen bins made invalid in the chosen way."""
    profile, profileNoise = profile.copy(), profileNoise.copy()
    if how == "nan":
        profile[isInvalid] = np.nan
    else:
        profileNoise[isInvalid] = 0.0
    return profile, profileNoise


@pytest.mark.parametrize("how", ["nan", "zero noise"])
@pytest.mark.parametrize(
    "isInvalid",
    [SLIT_BINS < -4.0, SLIT_BINS > 4.0, np.abs(SLIT_BINS + 1.0) < 0.2],
    ids=["leading", "trailing", "interior"],
)
def test_profile_sources_finds_both_sources_whatever_the_invalid_bins(isInvalid, how):
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(-2.5, 1.0), (1.5, 0.7)])
    profile, profileNoise = _invalidate(profile, profileNoise, isInvalid, how)

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [-2.5, 1.5], atol=BIN_ARCSEC)


@pytest.mark.parametrize("nValid", [0, 1])
def test_profile_sources_is_empty_with_fewer_than_two_valid_bins(nValid):
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(0.0, 1.0)])
    profileNoise[nValid:] = 0.0

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    assert sources.size == 0


def test_profile_sources_does_not_report_a_source_in_the_last_valid_bin():
    # ARRANGE
    profile, profileNoise = _source_profile(SLIT_BINS, [(0.0, 1.0), (SLIT_BINS[-1], 0.8)])

    # ACT
    sources = profile_sources(SLIT_BINS, profile, profileNoise)

    # ASSERT
    np.testing.assert_allclose(sources, [0.0], atol=BIN_ARCSEC)
