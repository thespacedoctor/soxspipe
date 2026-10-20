"""Contracts for the combined photon-noise signal-to-noise of an order's sky lines (DY-1663)."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.sky_line_snr import combine_line_snrs, sky_line_snrs

pytestmark = pytest.mark.unit

PIXEL_NM = 0.05
WAVELENGTH = np.arange(500.0, 520.0, PIXEL_NM)
READ_NOISE = 3.0
# A FLAT SKY OF 16 e- AND A READ NOISE OF 3 e- GIVE A PER-PIXEL PHOTON NOISE OF SQRT(16 + 9) = 5 e-
BACKGROUND = 16.0


def _sky(centres: np.ndarray, heights: np.ndarray) -> np.ndarray:
    sky = np.full_like(WAVELENGTH, BACKGROUND)
    for centre, height in zip(centres, heights, strict=True):
        sky += height * np.exp(-0.5 * ((WAVELENGTH - centre) / 0.08) ** 2)
    return sky


def test_one_line_snr_is_its_height_over_the_photon_noise() -> None:
    # ARRANGE
    centre = np.array([510.0])
    sky = _sky(centre, np.array([20.0]))

    # ACT
    snr = combine_line_snrs(sky_line_snrs(WAVELENGTH, sky, centre, READ_NOISE))

    # ASSERT
    assert snr == pytest.approx(4.0, rel=1e-6)


def test_line_snrs_add_in_quadrature() -> None:
    # ARRANGE
    centres = np.array([505.0, 510.0])
    sky = _sky(centres, np.array([15.0, 20.0]))

    # ACT
    snr = combine_line_snrs(sky_line_snrs(WAVELENGTH, sky, centres, READ_NOISE))

    # ASSERT
    assert snr == pytest.approx(5.0, rel=1e-6)


def test_a_line_offset_by_a_few_pixels_is_still_found() -> None:
    # ARRANGE
    centre = np.array([510.0])
    sky = _sky(centre + 3 * PIXEL_NM, np.array([20.0]))

    # ACT
    snr = combine_line_snrs(sky_line_snrs(WAVELENGTH, sky, centre, READ_NOISE))

    # ASSERT
    assert snr == pytest.approx(4.0, rel=1e-6)


def test_a_catalogue_line_with_no_sky_line_adds_nothing() -> None:
    # ARRANGE
    centres = np.array([505.0, 510.0])
    sky = _sky(centres[1:], np.array([20.0]))

    # ACT
    snr = combine_line_snrs(sky_line_snrs(WAVELENGTH, sky, centres, READ_NOISE))

    # ASSERT
    assert snr == pytest.approx(4.0, rel=1e-6)


def test_a_dip_at_a_catalogue_line_does_not_count_as_negative_signal() -> None:
    # ARRANGE
    centres = np.array([505.0, 510.0])
    sky = _sky(centres, np.array([-10.0, 20.0]))

    # ACT
    snr = combine_line_snrs(sky_line_snrs(WAVELENGTH, sky, centres, READ_NOISE))

    # ASSERT
    assert snr == pytest.approx(4.0, rel=1e-6)


def test_catalogue_lines_outside_the_order_and_nan_sky_are_ignored() -> None:
    # ARRANGE
    centres = np.array([490.0, 510.0, 530.0])
    sky = _sky(centres[1:2], np.array([20.0]))
    sky[:10] = np.nan

    # ACT
    snr = combine_line_snrs(sky_line_snrs(WAVELENGTH, sky, centres, READ_NOISE))

    # ASSERT
    assert snr == pytest.approx(4.0, rel=1e-6)


def test_an_order_with_no_catalogue_lines_has_zero_snr() -> None:
    # ACT
    snr = combine_line_snrs(sky_line_snrs(WAVELENGTH, np.full_like(WAVELENGTH, BACKGROUND), np.array([]), READ_NOISE))

    # ASSERT
    assert snr == 0.0


def test_per_line_snrs_are_sorted_brightest_first() -> None:
    # ARRANGE
    centres = np.array([505.0, 510.0, 515.0])
    sky = _sky(centres, np.array([5.0, 20.0, -10.0]))

    # ACT
    snrs = sky_line_snrs(WAVELENGTH, sky, centres, READ_NOISE)

    # ASSERT
    assert snrs == pytest.approx([4.0, 1.0, 0.0], rel=1e-6)


def test_unsorted_samples_give_the_same_snrs_as_sorted_ones() -> None:
    # ARRANGE
    centres = np.array([505.0, 510.0])
    sky = _sky(centres, np.array([15.0, 20.0]))
    order = np.random.default_rng(7).permutation(WAVELENGTH.size)

    # ACT
    snrs = sky_line_snrs(WAVELENGTH[order], sky[order], centres, READ_NOISE)

    # ASSERT
    assert snrs == pytest.approx([4.0, 3.0], rel=1e-6)
