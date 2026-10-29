"""The object mask and the significance gate of the sky fit on a nod A-B difference (DY-2057)."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.difference_sky_fit import gate_insignificant_model, signed_object_slit_ranges

pytestmark = pytest.mark.unit

SLIT_HALF_LENGTH = 5.5
PIXELS_PER_COLUMN = 60
COLUMNS = 400
NOISE = 10.0
TRACE_SIGMA = 0.4
POSITIVE_CENTRE, NEGATIVE_CENTRE = 1.5, -2.0
TRACE_PEAK = 6.0 * NOISE
# THE RANGE MAY BE WIDER THAN THE TRACE (THE MASK TAKES IN THE WINGS) BUT MUST COVER ITS CORE
CORE_HALF_WIDTH = 2.0 * TRACE_SIGMA
# TWO BRIGHT TRACES AND THEIR WINGS MAY TAKE UP TO THIS SHARE OF THE SLIT
MAX_MASKED_SLIT_FRACTION = 0.6


def _difference_pixels(seed: int, tracePeaks: tuple[float, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The slit position, flux and error of every pixel of an order with a positive and a negative trace."""
    rng = np.random.default_rng(seed)
    slitGrid = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, PIXELS_PER_COLUMN)
    slitPosition = np.tile(slitGrid, COLUMNS)
    profile = tracePeaks[0] * np.exp(-0.5 * ((slitPosition - POSITIVE_CENTRE) / TRACE_SIGMA) ** 2)
    profile = profile - tracePeaks[1] * np.exp(-0.5 * ((slitPosition - NEGATIVE_CENTRE) / TRACE_SIGMA) ** 2)
    flux = profile + rng.normal(0.0, NOISE, slitPosition.size)
    return slitPosition, flux, np.full(slitPosition.size, NOISE)


def _covers(objectRanges: list, centre: float) -> bool:
    return any(lower <= centre - CORE_HALF_WIDTH and centre + CORE_HALF_WIDTH <= upper for lower, upper in objectRanges)


def test_the_positive_and_the_negative_trace_each_get_a_slit_range() -> None:
    # ARRANGE
    slitPosition, flux, error = _difference_pixels(1, (TRACE_PEAK, TRACE_PEAK))

    # ACT
    objectRanges = signed_object_slit_ranges(slitPosition, flux, error)

    # ASSERT
    assert len(objectRanges) == 2
    assert _covers(objectRanges, POSITIVE_CENTRE)
    assert _covers(objectRanges, NEGATIVE_CENTRE)


def test_a_half_amplitude_negative_trace_is_still_masked() -> None:
    # ARRANGE
    slitPosition, flux, error = _difference_pixels(2, (TRACE_PEAK, 0.5 * TRACE_PEAK))

    # ACT
    objectRanges = signed_object_slit_ranges(slitPosition, flux, error)

    # ASSERT
    assert _covers(objectRanges, POSITIVE_CENTRE)
    assert _covers(objectRanges, NEGATIVE_CENTRE)


def test_a_uniform_sky_residual_and_noise_alone_give_no_object_range() -> None:
    # ARRANGE
    slitPosition, flux, error = _difference_pixels(3, (0.0, 0.0))

    # ACT
    objectRanges = signed_object_slit_ranges(slitPosition, flux + 4.0 * NOISE, error)

    # ASSERT
    assert objectRanges == []


def test_the_mask_leaves_most_of_the_slit_for_the_sky_fit() -> None:
    # ARRANGE
    slitPosition, flux, error = _difference_pixels(4, (TRACE_PEAK, TRACE_PEAK))

    # ACT
    objectRanges = signed_object_slit_ranges(slitPosition, flux, error)

    # ASSERT
    maskedWidth = sum(upper - lower for lower, upper in objectRanges)
    assert maskedWidth < MAX_MASKED_SLIT_FRACTION * 2 * SLIT_HALF_LENGTH


def _line_model(samples: int, amplitude: float) -> np.ndarray:
    """A flat model with one Gaussian skyline in the middle."""
    index = np.arange(samples)
    return amplitude * np.exp(-0.5 * ((index - samples // 2) / 6.0) ** 2)


def test_a_model_that_is_noise_sized_is_zeroed_everywhere() -> None:
    # ARRANGE
    samples = 3000
    rng = np.random.default_rng(5)
    model = rng.normal(0.0, 0.02 * NOISE, samples)
    weights = np.full(samples, 1.0 / NOISE)

    # ACT
    gated = gate_insignificant_model(model, weights, np.ones(samples, dtype=bool), 3.0)

    # ASSERT
    assert not gated.any()


def test_a_strong_skyline_survives_the_gate_unchanged_at_its_core() -> None:
    # ARRANGE
    samples = 3000
    model = _line_model(samples, 5.0 * NOISE)
    weights = np.full(samples, 1.0 / NOISE)

    # ACT
    gated = gate_insignificant_model(model, weights, np.ones(samples, dtype=bool), 3.0)

    # ASSERT
    core = slice(samples // 2 - 5, samples // 2 + 6)
    assert np.array_equal(gated[core], model[core])


def test_a_negative_skyline_survives_the_gate_as_a_positive_one_does() -> None:
    # ARRANGE
    samples = 3000
    model = -_line_model(samples, 5.0 * NOISE)
    weights = np.full(samples, 1.0 / NOISE)

    # ACT
    gated = gate_insignificant_model(model, weights, np.ones(samples, dtype=bool), 3.0)

    # ASSERT
    core = slice(samples // 2 - 5, samples // 2 + 6)
    assert np.array_equal(gated[core], model[core])


def test_the_gate_returns_a_new_array_and_leaves_the_model_alone() -> None:
    # ARRANGE
    model = np.full(500, 0.001)
    original = model.copy()

    # ACT
    gated = gate_insignificant_model(model, np.full(500, 1.0), np.ones(500, dtype=bool), 3.0)

    # ASSERT
    assert np.array_equal(model, original)
    assert gated is not model


# THE SHIFT-AND-ADD FRAME OF A NOD PAIR HOLDS THE OBJECT AT FULL AMPLITUDE AND TWICE AT HALF AMPLITUDE WITH THE
# OTHER SIGN
# (DY-1674). A BRIGHT STANDARD'S MOFFAT WINGS FILL MOST OF THE SLIT, SO THE SCATTER OF THE SLIT BINS IS NOT THE NOISE
BRIGHT_PEAK = 60.0 * NOISE
MOFFAT_ALPHA = 0.8
MOFFAT_BETA = 1.8
SHIFTED_NEGATIVE_CENTRE = 2.0 * POSITIVE_CENTRE - NEGATIVE_CENTRE


def _moffat(slitPosition: np.ndarray, centre: float) -> np.ndarray:
    return (1.0 + ((slitPosition - centre) / MOFFAT_ALPHA) ** 2) ** -MOFFAT_BETA


def test_a_bright_object_with_broad_wings_does_not_hide_its_two_half_amplitude_negative_traces() -> None:
    # ARRANGE
    rng = np.random.default_rng(6)
    slitGrid = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, PIXELS_PER_COLUMN)
    slitPosition = np.tile(slitGrid, COLUMNS)
    profile = BRIGHT_PEAK * (
        _moffat(slitPosition, POSITIVE_CENTRE)
        - 0.5 * _moffat(slitPosition, NEGATIVE_CENTRE)
        - 0.5 * _moffat(slitPosition, SHIFTED_NEGATIVE_CENTRE)
    )
    flux = profile + rng.normal(0.0, NOISE, slitPosition.size)

    # ACT
    objectRanges = signed_object_slit_ranges(slitPosition, flux, np.full(slitPosition.size, NOISE))

    # ASSERT
    for centre in (POSITIVE_CENTRE, NEGATIVE_CENTRE, SHIFTED_NEGATIVE_CENTRE):
        assert any(lower <= centre <= upper for lower, upper in objectRanges)
