"""Contracts for ``flat_slit_gradient.flat_slit_ratio``, the along-slit sky illumination from the flat (DY-2147)."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.flat_slit_gradient import MAX_SLOPE_PER_ARCSEC, flat_slit_ratio

pytestmark = pytest.mark.unit

ORDER = 12
COLUMNS = 50
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
    ratio = flat_slit_ratio(log, ORDER, flat, slit, KEEP_RANGE, isSky)

    # ASSERT
    expected = 1.0 / profile(np.clip(slit, *KEEP_RANGE))
    expected = expected / expected[isSky].mean()
    np.testing.assert_allclose(ratio, expected, rtol=1e-3)
    assert ratio[isSky].mean() == pytest.approx(1.0)


def test_edge_roll_off_outside_the_kept_range_does_not_change_the_ratio(log):
    # ARRANGE
    slit, column = _pixels()
    isSky = np.abs(slit) < 4.0
    smooth = flat_slit_ratio(log, ORDER, _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s), slit, KEEP_RANGE, isSky)
    rolloff = lambda s: (1.0 + FLAT_SLOPE * s) * np.where(np.abs(s) > 5.0, 0.3, 1.0)  # noqa: E731

    # ACT
    withRollOff = flat_slit_ratio(log, ORDER, _flat(slit, column, rolloff), slit, KEEP_RANGE, isSky)

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
        log, ORDER, np.where(isNarrow | (np.abs(slit) > 5.0), flat, np.nan), slit, KEEP_RANGE, isSky
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
    ratio = flat_slit_ratio(log, ORDER, flat, slit, KEEP_RANGE, isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert f"order {ORDER}" in _log_text(log)
    assert "slope" in _log_text(log)


def test_ratio_is_applied_when_the_slope_is_just_inside_the_limit(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + 0.9 * MAX_SLOPE_PER_ARCSEC * s)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, KEEP_RANGE, isSky)

    # ASSERT
    assert not (ratio == 1.0).all()
    assert "slope" not in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_the_order_has_no_kept_range(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, None, np.ones(slit.shape, dtype=bool))

    # ASSERT
    assert (ratio == 1.0).all()
    assert "no kept slit range" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_no_pixel_is_sky(log):
    # ARRANGE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, KEEP_RANGE, np.zeros(slit.shape, dtype=bool))

    # ASSERT
    assert (ratio == 1.0).all()
    assert "no sky pixels" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_the_flat_profile_is_not_positive(log):
    # ARRANGE: THE PROFILE FALLS BELOW ZERO AT BOTH ENDS OF THE KEPT RANGE, WITH NO NET SLOPE FOR THE SLOPE GUARD TO SEE
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 - 0.05 * s**2)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, KEEP_RANGE, isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "not positive" in _log_text(log)


def test_ratio_is_one_and_the_reason_is_logged_when_no_flat_pixel_is_inside_the_kept_range(log):
    # ARRANGE: THE KEPT RANGE LIES OUTSIDE EVERY SLIT POSITION OF THE ORDER
    slit, column = _pixels()
    flat = _flat(slit, column, lambda s: 1.0 + FLAT_SLOPE * s)
    isSky = np.abs(slit) < 4.0

    # ACT
    ratio = flat_slit_ratio(log, ORDER, flat, slit, (20.0, 30.0), isSky)

    # ASSERT
    assert (ratio == 1.0).all()
    assert "too few populated flat bins" in _log_text(log)
