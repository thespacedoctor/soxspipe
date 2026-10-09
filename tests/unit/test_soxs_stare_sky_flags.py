"""The stare recipe's bad-pixel QC ignores the pixels that sky subtraction flags as missing (DY-1284)."""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes.soxs_stare import soxs_stare

pytestmark = pytest.mark.unit

stareModule = importlib.import_module("soxspipe.recipes.soxs_stare")

FRAME_SHAPE = (3, 3)


def _frame(mask: np.ndarray) -> CCDData:
    return CCDData(
        np.ones(FRAME_SHAPE),
        unit=u.electron,
        mask=mask,
        uncertainty=StdDevUncertainty(np.ones(FRAME_SHAPE)),
    )


def test_the_bad_pixel_qc_does_not_count_pixels_that_only_the_sky_subtraction_flagged(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detectorMask = np.zeros(FRAME_SHAPE, dtype=bool)
    detectorMask[0, 0] = True
    skySubtractionMask = detectorMask.copy()
    skySubtractionMask[1:, :] = True
    objectFrame = _frame(detectorMask)
    calls: dict[str, list] = {"excludeMask": []}

    class FakeSubtractSky:
        def __init__(self, **_: Any) -> None:
            return None

        def subtract(self) -> tuple[Any, Any, Any, str, str]:
            return _frame(skySubtractionMask), _frame(skySubtractionMask), None, "qc", "products"

    def record_generic_checks(**kwargs: Any) -> str:
        calls["excludeMask"].append(kwargs.get("excludeMask"))
        return "qc"

    monkeypatch.setattr(stareModule, "subtract_sky", FakeSubtractSky)
    monkeypatch.setattr(stareModule, "generic_quality_checks", record_generic_checks)
    monkeypatch.setattr(stareModule, "spectroscopic_image_quality_checks", lambda **_: "qc")
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.settings = {}
    recipe.recipeSettings = {"sky-subtraction": {"subtract_sky": True}}
    recipe.subtractSky = True
    recipe.qc = "qc"
    recipe.products = "products"
    recipe.sofName = "stare"
    recipe.recipeName = "soxs-stare"
    recipe.startNightDate = "2024-01-02"
    recipe.debug = False
    recipe._write_sky_products = lambda *_: "product.fits"

    recipe._subtract_sky(objectFrame, "map.fits", "disp.fits", "orders.fits")

    [excludeMask] = calls["excludeMask"]
    expected = skySubtractionMask & ~detectorMask
    np.testing.assert_array_equal(excludeMask, expected)
