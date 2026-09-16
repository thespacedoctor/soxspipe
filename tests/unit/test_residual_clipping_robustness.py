"""Contracts keeping the sigma-clip decision insensitive to last-bit arithmetic noise.

`fit_polynomials` clips lines whose residual exceeds a sigma threshold, then
re-fits on what survives. That loop amplifies: a residual differing in its last
bits can land on either side of the threshold, which changes the next fit, which
changes the next clip. In the real-data NIR offset reduction this turned a
difference of order 1e-16 into five differently-clipped lines and a twelve-bin
shift in the red end of the merged spectrum.

Residuals are therefore quantised before they are compared against the threshold,
at a precision far finer than the measurement but far coarser than float noise.
"""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.create_dispersion_map import (
    RESIDUAL_QUANTISATION_DECIMALS,
    quantise_residuals,
)

pytestmark = pytest.mark.unit


def test_quantisation_is_far_finer_than_any_real_measurement():
    """The quantum is negligible against centroiding precision, which is ~0.01 px."""
    # ARRANGE / ACT
    quantum = 10.0**-RESIDUAL_QUANTISATION_DECIMALS

    # ASSERT
    assert quantum <= 1e-6, "the quantum must be far below any measurable residual"


def test_residuals_differing_only_in_their_last_bits_quantise_to_one_value():
    """Two solves that differ by float noise yield identical clipping inputs."""
    # ARRANGE: THE SAME RESIDUAL AS TWO SOLVES ON DIFFERENT HARDWARE MIGHT PRODUCE
    residual = 0.4567891234567
    noisy = np.array([residual, np.nextafter(residual, 1.0), np.nextafter(residual, 0.0)])

    # ACT
    quantised = quantise_residuals(noisy)

    # ASSERT
    assert len(np.unique(quantised)) == 1


def test_quantisation_preserves_residuals_that_genuinely_differ():
    """A real difference in the data survives quantisation untouched."""
    # ARRANGE
    residuals = np.array([0.10, 0.25, 1.75, -0.5])

    # ACT
    quantised = quantise_residuals(residuals)

    # ASSERT
    np.testing.assert_allclose(quantised, residuals, rtol=0, atol=1e-9)


def test_quantisation_leaves_nan_alone():
    """Missing residuals stay missing rather than becoming a number."""
    # ARRANGE
    residuals = np.array([0.1, np.nan, 0.3])

    # ACT
    quantised = quantise_residuals(residuals)

    # ASSERT
    assert np.isnan(quantised[1])
    np.testing.assert_allclose(quantised[[0, 2]], [0.1, 0.3], rtol=0, atol=1e-9)


def test_a_threshold_comparison_is_stable_across_last_bit_noise():
    """The clip decision itself does not flip when residuals carry float noise."""
    # ARRANGE: A RESIDUAL SITTING EXACTLY ON THE THRESHOLD, AS THE CLIP LOOP FINDS
    threshold = 0.5
    onEdge = 0.5
    runA = np.array([np.nextafter(onEdge, 1.0)])
    runB = np.array([np.nextafter(onEdge, 0.0)])

    # ACT
    clippedA = quantise_residuals(runA) > threshold
    clippedB = quantise_residuals(runB) > threshold

    # ASSERT
    assert clippedA.tolist() == clippedB.tolist()
