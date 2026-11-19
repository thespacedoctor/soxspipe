"""The second-pass flat normalisation leaves out the pixels where the first-pass master flat is unusable (DY-2279)."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest
from astropy import stats as astropy_stats
from astropy.nddata import CCDData

from tests.factories.unusable_flat import (
    FLAT_SHAPE,
    USABLE_FLAT_VALUE,
    flat_with_unusable_pixels,
)
from tests.unit.test_soxs_mflat_worker_characterization import _frame
from tests.unit.test_soxs_mflat_worker_characterization import _recipe as _worker_recipe

pytestmark = pytest.mark.unit

FRAME_VALUE = 100.0


def _normalise_second_pass(log: Any, monkeypatch: pytest.MonkeyPatch) -> tuple[list[CCDData], list[np.ndarray]]:
    """Run the second pass with no pixel excluded by the mask, and record the sample behind each level."""
    realClippedStats = astropy_stats.sigma_clipped_stats
    samples: list[np.ndarray] = []

    def spy(data: np.ndarray, *args: object, **kwargs: object) -> tuple[float, float, float]:
        samples.append(np.array(data))
        return realClippedStats(data, *args, **kwargs)

    monkeypatch.setattr(astropy_stats, "sigma_clipped_stats", spy)
    frame = _frame(np.full(FLAT_SHAPE, FRAME_VALUE))
    noExclusion = np.zeros(FLAT_SHAPE, dtype=bool)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        normalised = _worker_recipe(log)._normalise_flats_second_pass([frame], noExclusion, flat_with_unusable_pixels())
    return normalised, samples


def test_second_pass_level_uses_only_the_pixels_where_the_flat_is_usable(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE / ACT
    _, samples = _normalise_second_pass(log, monkeypatch)

    # ASSERT
    assert len(samples) == 1
    np.testing.assert_allclose(samples[0], np.full(2, FRAME_VALUE / USABLE_FLAT_VALUE))


def test_second_pass_normalises_the_frame_to_the_level_of_the_usable_pixels(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE / ACT
    normalised, _ = _normalise_second_pass(log, monkeypatch)

    # ASSERT
    np.testing.assert_allclose(normalised[0].data, np.full(FLAT_SHAPE, 2.0))
