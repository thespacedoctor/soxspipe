"""The stare recipe logs a median-offset flat-fit fallback at INFO, and a zero-offset one as a warning."""

from __future__ import annotations

from typing import Any

import pytest

import soxspipe.commonutils.slit_illumination_mask as slit_illumination_mask
from soxspipe.commonutils.slit_illumination_mask import fit_offset
from soxspipe.recipes.soxs_stare import soxs_stare

pytestmark = pytest.mark.unit


def _recipe(log: Any) -> soxs_stare:
    recipe = object.__new__(soxs_stare)
    recipe.log = log
    recipe.recipeSettings = {"sky-subtraction": {"flat-slit-gradient": True}}
    recipe.slitIlluminationFitRanges = {11: (-4.0, 4.0), 12: (-4.0, 4.0), 13: (-4.0, 4.0)}
    return recipe


def test_a_median_offset_fallback_is_logged_at_info_and_a_zero_offset_fallback_as_a_warning(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    offsets = {
        11: fit_offset(0.1, 0.0, "measured"),
        12: fit_offset(-0.056, 0.0, "median"),
        13: fit_offset(0.0, 0.3, "zero"),
    }
    monkeypatch.setattr(slit_illumination_mask, "flat_fit_offsets", lambda result: offsets)

    # ACT
    _recipe(log)._warn_flat_fit_fallbacks(result=object())

    # ASSERT
    levels = {message.split(":")[0]: level for level, message in log.messages}
    assert levels == {"order 12": "info", "order 13": "warning"}
