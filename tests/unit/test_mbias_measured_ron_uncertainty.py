"""The master bias `ERRS` follows the read noise measured on the raw frames, not the header value (DY-2438)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.recipes.soxs_mbias import soxs_mbias

pytestmark = pytest.mark.unit

HEADER_RON = 3.3
MEASURED_RON = 3.83
FRAME_COUNT = 5


def _recipe(log: Any) -> soxs_mbias:
    """Return a recipe holding the state the uncertainty step reads."""
    recipe = soxs_mbias.__new__(soxs_mbias)
    recipe.log = log
    recipe.detectorParams = {"ron": HEADER_RON * u.electron}
    return recipe


def _stacked_bias(shape: tuple[int, int] = (4, 4)) -> CCDData:
    """Return a stacked bias whose `ERRS` is the header read noise over the root of the frame count.

    This is what `clip_and_stack` produces from raw bias frames prepared with a constant header-noise map.
    """
    frame = CCDData(np.zeros(shape), unit=u.electron)
    frame.uncertainty = np.full(shape, HEADER_RON / np.sqrt(FRAME_COUNT), dtype=np.float32)
    return frame


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_uncertainty_follows_the_measured_raw_ron_over_root_n(log: Any) -> None:
    """Every pixel of the master bias `ERRS` becomes the measured raw RON over the root of the frame count."""
    # ARRANGE
    recipe = _recipe(log)
    frame = _stacked_bias()

    # ACT
    recipe._scale_uncertainty_to_raw_ron(frame, MEASURED_RON)

    # ASSERT
    expected = MEASURED_RON / np.sqrt(FRAME_COUNT)
    assert np.median(frame.uncertainty.array) == pytest.approx(expected, rel=1e-5)
    assert _warnings(log) == []


def test_uncertainty_keeps_the_per_pixel_frame_count_when_scaling(log: Any) -> None:
    """A pixel that fewer frames contributed to keeps its larger error after the rescale."""
    # ARRANGE
    recipe = _recipe(log)
    frame = _stacked_bias()
    frame.uncertainty = np.where(
        np.arange(16).reshape(4, 4) == 0,
        HEADER_RON / np.sqrt(FRAME_COUNT - 2),
        HEADER_RON / np.sqrt(FRAME_COUNT),
    ).astype(np.float32)

    # ACT
    recipe._scale_uncertainty_to_raw_ron(frame, MEASURED_RON)

    # ASSERT
    assert frame.uncertainty.array[0, 0] == pytest.approx(MEASURED_RON / np.sqrt(FRAME_COUNT - 2), rel=1e-5)
    assert frame.uncertainty.array[1, 1] == pytest.approx(MEASURED_RON / np.sqrt(FRAME_COUNT), rel=1e-5)


@pytest.mark.parametrize("rawRon", [None, 0.0, float("nan"), float("inf"), -1.0])
def test_uncertainty_keeps_header_ron_and_warns_when_the_measurement_failed(log: Any, rawRon: float | None) -> None:
    """Without a usable measured read noise the header-based `ERRS` stays and a warning says so."""
    # ARRANGE
    recipe = _recipe(log)
    frame = _stacked_bias()
    original = frame.uncertainty.array.copy()

    # ACT
    recipe._scale_uncertainty_to_raw_ron(frame, rawRon)

    # ASSERT
    np.testing.assert_array_equal(frame.uncertainty.array, original)
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "unmasked pixels" in warnings[0]
    assert "header read noise" in warnings[0]
    assert f"{HEADER_RON:0.2f}" in warnings[0]
    assert f"RAW RON = {rawRon}" in warnings[0]


@pytest.mark.parametrize("headerRon", [0.0, float("nan"), float("inf")])
def test_uncertainty_is_not_scaled_when_the_header_ron_is_unusable(log: Any, headerRon: float) -> None:
    """A header read noise that cannot be a divisor leaves the uncertainty alone and warns once."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.detectorParams = {"ron": headerRon * u.electron}
    frame = _stacked_bias()
    original = frame.uncertainty.array.copy()

    # ACT
    recipe._scale_uncertainty_to_raw_ron(frame, MEASURED_RON)

    # ASSERT
    np.testing.assert_array_equal(frame.uncertainty.array, original)
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "header read noise" in warnings[0]
