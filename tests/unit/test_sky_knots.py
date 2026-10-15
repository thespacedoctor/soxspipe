"""Uniform b-spline knot placement for the sky fit (DY-1287)."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.sky_knots import knots_with_samples, native_wavelength_step, uniform_interior_knots

pytestmark = pytest.mark.unit


def test_native_wavelength_step_is_the_median_step_between_dispersion_neighbours_on_one_line():
    # ARRANGE: TWO CROSS-DISPERSION LINES, 0.02 NM PER PIXEL, THE SECOND LINE OFFSET IN WAVELENGTH (A TILTED SLIT)
    dispersionPixel = np.tile(np.arange(10), 2)
    crossPixel = np.repeat([0, 1], 10)
    wavelength = 500.0 + 0.02 * dispersionPixel + 0.007 * crossPixel
    shuffle = np.random.default_rng(1).permutation(wavelength.size)

    # ACT
    step = native_wavelength_step(wavelength[shuffle], dispersionPixel[shuffle], crossPixel[shuffle])

    # ASSERT: THE 0.007 NM OFFSET BETWEEN LINES IS NOT A DISPERSION STEP
    assert step == pytest.approx(0.02)


def test_native_wavelength_step_skips_gaps_in_the_dispersion_direction():
    # ARRANGE: A MISSING PIXEL (A CLIPPED BAD PIXEL) MUST NOT COUNT AS ONE DOUBLE-SIZED STEP
    dispersionPixel = np.array([0, 1, 2, 5, 6, 9, 12])
    crossPixel = np.zeros(dispersionPixel.size, dtype=int)
    wavelength = 700.0 - 0.05 * dispersionPixel

    # ACT
    step = native_wavelength_step(wavelength, dispersionPixel, crossPixel)

    # ASSERT: THE WAVELENGTH FALLS ALONG THE DISPERSION AXIS, BUT THE STEP IS A SIZE
    assert step == pytest.approx(0.05)


def test_native_wavelength_step_raises_without_any_neighbouring_pair():
    with pytest.raises(ValueError, match="neighbouring"):
        native_wavelength_step(np.array([500.0, 501.0]), np.array([0, 2]), np.array([0, 0]))


def test_uniform_interior_knots_are_evenly_spaced_and_strictly_inside_the_range():
    knots = uniform_interior_knots(500.0, 501.0, 0.3)

    np.testing.assert_allclose(knots, [500.3, 500.6, 500.9])


def test_uniform_interior_knots_never_place_a_knot_on_the_last_wavelength():
    knots = uniform_interior_knots(500.0, 501.0, 0.25)

    np.testing.assert_allclose(knots, [500.25, 500.5, 500.75])


@pytest.mark.parametrize("spacing", [0.0, -0.1, np.nan, np.inf])
def test_uniform_interior_knots_reject_a_spacing_that_is_not_finite_and_positive(spacing):
    with pytest.raises(ValueError, match="spacing"):
        uniform_interior_knots(500.0, 501.0, spacing)


def test_knots_with_samples_keeps_one_sample_per_interval_by_default():
    # ARRANGE: NO SAMPLE BETWEEN 502.5 AND 502.7, SO THE REDDER OF THE TWO GOES
    wavelength = np.arange(500.0, 506.0, 1.0)

    # ACT
    kept = knots_with_samples(np.array([502.5, 502.7, 504.5]), wavelength)

    # ASSERT
    np.testing.assert_allclose(kept, [502.5, 504.5])


def test_knots_with_samples_merges_intervals_with_fewer_than_min_samples():
    # ARRANGE: 20 SAMPLES; KNOTS EVERY 2 SAMPLES, MIN 5 SAMPLES PER INTERVAL
    wavelength = np.arange(20, dtype=float) + 0.5
    knots = np.arange(2.0, 20.0, 2.0)

    # ACT
    kept = knots_with_samples(knots, wavelength, minSamples=5)

    # ASSERT: EVERY INTERVAL, THE TWO END INTERVALS INCLUDED, HOLDS AT LEAST 5 SAMPLES
    edges = np.concatenate(([-np.inf], kept, [np.inf]))
    counts = np.histogram(wavelength, bins=edges)[0]
    assert counts.min() >= 5
    np.testing.assert_allclose(kept, [6.0, 12.0])


def test_knots_with_samples_drops_the_last_knots_when_the_red_end_interval_is_too_sparse():
    wavelength = np.arange(10, dtype=float)

    kept = knots_with_samples(np.array([3.5, 7.5, 8.5]), wavelength, minSamples=3)

    np.testing.assert_allclose(kept, [3.5])


def test_knots_with_samples_ignores_a_sample_exactly_on_a_knot():
    wavelength = np.array([0.0, 1.0, 2.0, 3.0])

    kept = knots_with_samples(np.array([1.0, 2.0]), wavelength, minSamples=1)

    # NOTHING LIES STRICTLY BETWEEN 1.0 AND 2.0
    np.testing.assert_allclose(kept, [1.0])


def test_knots_with_samples_returns_no_knots_without_samples():
    assert knots_with_samples(np.array([1.0]), np.array([])).size == 0
