"""Contracts for ``flat_slit_gradient.flat_slit_ratio``, the along-slit sky illumination from the flat (DY-2147)."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from soxspipe.commonutils.flat_slit_gradient import (
    MAX_SLOPE_PER_ARCSEC,
    SOXS_NIR_GRADIENT_SCALE_BY_ORDER,
    _profile_extremes,
    flat_slit_ratio,
    gradient_scale,
)

pytestmark = pytest.mark.unit

ORDER = 12
COLUMNS = 120
# THE WAVELENGTH OF COLUMN 0 AND THE WAVELENGTH STEP PER COLUMN (NM)
START_WAVELENGTH = 1700.0
DISPERSION = 0.37
# SLIT POSITIONS FROM -6 TO +6 ARCSEC IN 0.05 ARCSEC STEPS, ON EVERY COLUMN, SO EVERY 0.1 ARCSEC BIN HOLDS 100 PIXELS
SLIT_STEPS = np.round(np.arange(-120, 121) * 0.05, 2)
KEEP_RANGE = (-5.0, 5.0)
# THE MEASURED NIR GRADIENT OF THE FLAT (DY-1990), IN FRACTION PER ARCSEC
FLAT_SLOPE = -0.0035


def _pixels():
    """The slit position and the column of every pixel of one synthetic order."""
    slit = np.tile(SLIT_STEPS, COLUMNS)
    column = np.repeat(np.arange(COLUMNS), SLIT_STEPS.size)
    return slit, column


def _wavelength(column):
    """The wavelength of each pixel, from its column."""
    return START_WAVELENGTH + DISPERSION * column


def _flat(slit, column, profile):
    """A flat that is a spectral shape times the slit profile, so the profile is the same on every column."""
    return (500.0 + 10.0 * column) * profile(slit)


def _log_text(log):
    return " ".join(message for _, message in log.messages)


def test_ratio_is_the_inverse_of_the_flat_profile_normalised_to_mean_one(log):
    # ARRANGE
    slit, column = _pixels()
    profile = lambda s: 1.0 + FLAT_SLOPE * s + 0.0004 * s**2  # noqa: E731
    flat = _flat(slit, column, profile)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    expected = 1.0 / profile(np.clip(slit, *KEEP_RANGE))
    expected = expected / expected[isSky].mean()
    np.testing.assert_allclose(ratio, expected, rtol=1e-3)
    assert ratio[isSky].mean() == pytest.approx(1.0)


def test_edge_roll_off_outside_the_kept_range_does_not_change_the_ratio(log):
    # ARRANGE
    slit, column = _pixels()
    isSky = np.abs(slit) < 4.0
    smooth = flat_slit_ratio(
        log, ORDER, _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s), slit, _wavelength(column), KEEP_RANGE, isSky
    )
    rolloff = lambda s: (1.0 + FLAT_SLOPE * s) * np.where(np.abs(s) > 5.0, 0.3, 1.0)  # noqa: E731

    # ACT
    withRollOff = flat_slit_ratio(
        log, ORDER, _flat(slit, column, rolloff), slit, _wavelength(column), KEEP_RANGE, isSky
    )

    # ASSERT
    np.testing.assert_allclose(withRollOff, smooth, rtol=1e-9)


def test_ratio_is_one_and_the_reason_is_logged_when_too_few_slit_bins_are_populated(log):
    # ARRANGE: ONLY 1.5 ARCSEC OF THE KEPT RANGE HOLDS PIXELS, WHICH IS 15 BINS OF 0.1 ARCSEC
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)
    isNarrow = np.abs(slit) < 0.75
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(
        log,
        ORDER,
        np.where(isNarrow | (np.abs(slit) > 5.0), flat, np.nan),
        slit,
        _wavelength(column),
        KEEP_RANGE,
        isSky,
    )

    # ASSERT
    assert (ratio == 1.0).all()
    assert ratio.shape == slit.shape
    assert f"order {ORDER}" in _log_text(log)
    assert "too few populated flat bins" in _log_text(log)


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_ratio_is_one_and_the_reason_is_logged_when_the_slope_is_too_steep(log, sign):
    # ARRANGE: A SLOPE TWICE THE LIMIT, IN EITHER DIRECTION
    slit, column = _pixels()
    steep = sign * 2 * MAX_SLOPE_PER_ARCSEC
    flat = _flat(slit, column, lambda s: 1.0 + steep * s)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert f"order {ORDER}" in _log_text(log)
    assert "slope" in _log_text(log)
    assert "%/arcsec" in _log_text(log)
    assert "per arcsec" not in _log_text(log)


def test_ratio_is_applied_when_the_slope_is_just_inside_the_limit(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + 0.9 * MAX_SLOPE_PER_ARCSEC * s)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    assert not (ratio == 1.0).all()
    assert "slope" not in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_the_order_has_no_kept_range(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), None, np.ones(slit.shape, dtype=bool))

    # ASSERT
    assert (ratio == 1.0).all()
    assert "no slit range to fit" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_no_pixel_is_sky(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, np.zeros(slit.shape, dtype=bool))

    # ASSERT
    assert (ratio == 1.0).all()
    assert "no sky pixels" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_the_flat_profile_is_not_positive(log):
    # ARRANGE: THE PROFILE FALLS BELOW ZERO AT BOTH ENDS OF THE KEPT RANGE, WITH NO NET SLOPE FOR THE SLOPE GUARD TO SEE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 - 0.05 * s**2)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "not positive" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_no_flat_pixel_is_inside_the_kept_range(log):
    # ARRANGE: THE KEPT RANGE LIES OUTSIDE EVERY SLIT POSITION OF THE ORDER
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), (20.0, 30.0), isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "too few populated flat bins" in _log_text(log)


def test_a_slit_coverage_that_changes_with_wavelength_does_not_tilt_the_profile(log):
    # ARRANGE: THE FLAT BRIGHTENS TO THE RED, AND THE SLIT ABOVE +2 ARCSEC IS ONLY LIT IN THE RED HALF OF THE ORDER, AS
    # A CURVED ORDER IS. A RAW MEDIAN PER SLIT BIN WOULD READ THE RED BRIGHTNESS AS A RISE ALONG THE SLIT
    slit, column = _pixels()
    isCovered = (slit <= 2.0) | (column >= COLUMNS // 2)
    slit, column = slit[isCovered], column[isCovered]
    profile = lambda s: 1.0 + FLAT_SLOPE * s  # noqa: E731
    flat = (500.0 + 400.0 * column / COLUMNS) * profile(slit)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    expected = 1.0 / profile(np.clip(slit, *KEEP_RANGE))
    np.testing.assert_allclose(ratio, expected / expected[isSky].mean(), rtol=3e-3)


def test_ratio_is_one_and_the_reason_is_logged_when_the_profile_bends_steeply_with_no_net_slope(log):
    # ARRANGE: A U-SHAPED PROFILE, ABOUT 6 %/ARCSEC STEEP AT THE ENDS OF THE KEPT RANGE AND LEVEL IN THE MIDDLE, SO ITS
    # END-TO-END SLOPE IS ZERO
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + 0.006 * s**2)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "slope" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_the_profile_dips_below_zero_between_check_points(log):
    # ARRANGE: THE PROFILE'S MINIMUM, JUST BELOW ZERO, SITS BETWEEN TWO POINTS OF A 0.1 ARCSEC CHECK GRID, AT BOTH OF
    # WHICH IT IS POSITIVE. THE PIXELS SIT AT +-0.025 ARCSEC FROM EACH BIN CENTRE, SO EACH BIN MEDIAN IS THE PROFILE
    # AT THE BIN CENTRE, LESS A CONSTANT 1.25e-5 THAT THE DIP MUST BEAT
    slit, column = _pixels()
    slit = slit + 0.025
    flat = _flat(slit, column, lambda s: 0.02 * (s - 0.05) ** 2 - 3e-5)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "not positive" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_once_when_there_is_no_master_flat(log):
    # ARRANGE
    slit, column = _pixels()

    # ACT
    ratio = flat_slit_ratio(log, ORDER, None, slit, _wavelength(column), KEEP_RANGE, np.ones(slit.shape, dtype=bool))

    # ASSERT
    assert (ratio == 1.0).all()
    assert _log_text(log).count("there is no master flat") == 1


def test_a_real_flat_roll_off_toward_the_slit_ends_is_corrected_not_refused(log):
    # ARRANGE: A -0.35 %/ARCSEC SLOPE WITH A CURVATURE THAT MAKES THE GRADIENT 2 %/ARCSEC AT THE ENDS OF THE KEPT RANGE,
    # AS THE SOXS NIR MASTER FLATS HAVE
    slit, column = _pixels()
    profile = lambda s: 1.0 + FLAT_SLOPE * s - 0.002 * s**2  # noqa: E731
    flat = _flat(slit, column, profile)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    expected = 1.0 / profile(np.clip(slit, *KEEP_RANGE))
    np.testing.assert_allclose(ratio, expected / expected[isSky].mean(), rtol=2e-3)


def test_ratio_is_one_and_the_reason_is_logged_when_the_flat_is_not_finite_anywhere(log):
    # ARRANGE
    slit, column = _pixels()
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, np.full(slit.shape, np.nan), slit, _wavelength(column), KEEP_RANGE, isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "too few populated flat bins" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_no_usable_pixel_is_in_the_slit_interior(log):
    # ARRANGE: THE FIT RANGE LIES ENTIRELY BEYOND THE 3 ARCSEC INTERIOR THAT SETS THE LEVEL OF EACH WAVELENGTH CHUNK
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), (3.5, 5.0), np.abs(slit) < 4.0)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "too few populated flat bins" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_every_chunk_level_is_not_positive(log):
    # ARRANGE: A FLAT THAT IS NEGATIVE EVERYWHERE, SO NO CHUNK HAS A LEVEL TO DIVIDE BY
    slit, column = _pixels()
    flat = -_flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, np.abs(slit) < 4.0)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "too few populated flat bins" in _log_text(log)


def test_a_quadratic_term_too_small_to_have_a_vertex_gives_the_extremes_of_a_line():
    # ARRANGE: A SUBNORMAL QUADRATIC TERM, WHOSE VERTEX -b / 2a OVERFLOWS
    coefficients = np.array([5e-324, -0.0035, 1.0])

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        minimum, maxGradient, meanLevel, meanGradient = _profile_extremes(coefficients, -5.0, 5.0)

    # ASSERT
    assert minimum == pytest.approx(1.0 - 0.0035 * 5.0)
    assert maxGradient == pytest.approx(0.0035)
    assert meanLevel == pytest.approx(1.0)
    assert meanGradient == pytest.approx(-0.0035)


def _linear_ratio(log, gradientScale, **kwargs):
    """The ratio for a flat whose slit profile is a straight line of slope ``FLAT_SLOPE``, and the slit positions."""
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)
    isSky = np.abs(slit) < 4.0
    ratio = flat_slit_ratio(
        log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky, gradientScale=gradientScale, **kwargs
    )
    return ratio, slit, isSky


def test_ratio_follows_k_times_the_flat_gradient_when_the_gradient_scale_is_k(log):
    # ARRANGE
    k = 0.6

    # ACT
    ratio, slit, isSky = _linear_ratio(log, k)

    # ASSERT: THE PROFILE IS THE MEAN LEVEL PLUS K TIMES ITS DEVIATION FROM THE MEAN, WHICH IS 1 FOR A SYMMETRIC LINE
    expected = 1.0 / (1.0 + k * FLAT_SLOPE * np.clip(slit, *KEEP_RANGE))
    expected = expected / expected[isSky].mean()
    np.testing.assert_allclose(ratio, expected, rtol=1e-3)
    assert ratio[isSky].mean() == pytest.approx(1.0)


def test_ratio_is_the_same_with_a_gradient_scale_of_one_as_without_one(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s + 0.0004 * s**2)
    isSky = np.abs(slit) < 4.0
    unscaled = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky)

    # ACT
    scaled = flat_slit_ratio(log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, isSky, gradientScale=1.0)

    # ASSERT
    np.testing.assert_array_equal(scaled, unscaled)


def test_the_gradient_scale_does_not_change_the_gradient_the_log_reports(log):
    # ARRANGE
    unscaledLog = type(log)()
    _linear_ratio(unscaledLog, 1.0)

    # ACT
    _linear_ratio(log, 0.6)

    # ASSERT
    gradientText = _log_text(unscaledLog).split(", scaled by")[0]
    assert gradientText in _log_text(log)
    assert "gradient of -0.35 %/arcsec" in gradientText


@pytest.mark.parametrize("k", [0.6, 1.0])
def test_the_log_gives_the_gradient_scale_next_to_the_flat_gradient(log, k):
    # ACT
    _linear_ratio(log, k)

    # ASSERT
    assert "the master flat gives a gradient of -0.35 %/arcsec along the slit" in _log_text(log)
    assert f"scaled by k = {k:.2f}" in _log_text(log)


def test_a_refused_flat_stays_refused_whatever_the_gradient_scale(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + 2 * MAX_SLOPE_PER_ARCSEC * s)

    # ACT
    ratio = flat_slit_ratio(
        log, ORDER, flat, slit, _wavelength(column), KEEP_RANGE, np.abs(slit) < 4.0, gradientScale=0.6
    )

    # ASSERT
    assert (ratio == 1.0).all()
    assert "beyond the limit" in _log_text(log)


@pytest.mark.parametrize(
    ("order", "expected"),
    [(15, 0.80), (16, 0.75), (17, 0.60), (18, 0.70), (19, 0.60), (20, 0.60), (21, 0.60), (22, 0.60)],
)
def test_soxs_nir_orders_in_the_table_have_their_scale(order, expected):
    # ACT
    k = gradient_scale("SOXS", "NIR", order)

    # ASSERT
    assert k == expected


@pytest.mark.parametrize("order", [10, 11, 12, 13, 14, 23, 24])
def test_soxs_nir_orders_outside_the_table_keep_the_flat_gradient(order):
    # ACT
    k = gradient_scale("SOXS", "NIR", order)

    # ASSERT
    assert k == 1.0


@pytest.mark.parametrize(("instrument", "arm"), [("SOXS", "VIS"), ("XSHOOTER", "NIR"), ("XSHOOTER", "VIS")])
def test_vis_and_xshooter_keep_the_flat_gradient_in_every_order(instrument, arm):
    # ACT
    scales = {gradient_scale(instrument, arm, order) for order in range(1, 30)}

    # ASSERT
    assert scales == {1.0}


def test_the_instrument_and_arm_are_read_without_regard_to_case():
    # ACT
    k = gradient_scale("soxs", "nir", 17)

    # ASSERT
    assert k == SOXS_NIR_GRADIENT_SCALE_BY_ORDER[17]
