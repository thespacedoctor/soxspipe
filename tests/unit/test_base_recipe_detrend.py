"""Calibration-subtraction contracts for the base recipe."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes import base_recipe
from tests.factories import instrument_header, product_table

pytestmark = pytest.mark.unit


def _frame(value: float, exposure: float = 10.0) -> CCDData:
    """Return a uniform, prepared calibration frame."""
    shape = (4, 4)
    return CCDData(
        np.full(shape, value, dtype=np.float32),
        unit=u.electron,
        meta=instrument_header(
            overrides={
                "DPR_TYPE": "OBJECT",
                "DPR_TECH": "ECHELLE,SLIT",
                "DPR_CATG": "SCIENCE",
                "EXPTIME": exposure,
            }
        ),
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(np.ones(shape), unit=u.electron),
    )


def _recipe(log: Any) -> base_recipe:
    """Return the minimum recipe state needed by ``detrend``."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: keyword
    recipe.detectorParams = {}
    recipe.recipeSettings = {"subtract_background": False}
    recipe.darkDetrendWarningIssued1 = False
    recipe.darkDetrendWarningIssued2 = False
    recipe.sofName = "synthetic"
    recipe.recipeName = "soxs-stare"
    recipe.settings = {}
    recipe.products = product_table().iloc[0:0].copy()
    recipe.qc = pd.DataFrame()
    recipe.startNightDate = "2024-01-02"
    return recipe


def test_detrend_requires_a_calibration_frame(log: Any) -> None:
    """Calling detrend without a calibration fails before modifying the input."""
    recipe = _recipe(log)

    with pytest.raises(TypeError, match="needs at least a master-bias"):
        recipe.detrend(_frame(20.0))


def test_detrend_applies_bias_matching_dark_and_flat_in_order(log: Any) -> None:
    """Bias, equal-exposure dark, and normalized flat retain CCDData semantics."""
    recipe = _recipe(log)

    calibrated = recipe.detrend(
        _frame(20.0),
        master_bias=_frame(2.0),
        dark=_frame(3.0),
        master_flat=_frame(0.5),
    )

    np.testing.assert_allclose(calibrated.data, np.full((4, 4), 30.0))
    assert calibrated.data.dtype == np.float32
    assert not recipe.darkDetrendWarningIssued1
    assert not recipe.darkDetrendWarningIssued2


def test_detrend_scales_a_different_exposure_dark_once(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """A different-exposure dark is scaled and emits its warning only once."""
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    quicklookCalls: list[dict[str, object]] = []
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: quicklookCalls.append(kwargs))
    dark = _frame(2.0, exposure=5.0)

    calibrated = recipe.detrend(_frame(20.0, exposure=10.0), dark=dark)

    np.testing.assert_allclose(calibrated.data, np.full((4, 4), 16.0))
    assert recipe.darkDetrendWarningIssued1
    assert recipe.darkDetrendWarningIssued2
    assert len(quicklookCalls) == 2
    assert quicklookCalls[0]["CCDObject"] is dark
    assert quicklookCalls[1]["CCDObject"] is calibrated


def test_detrend_delegates_background_subtraction_when_enabled(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An order table delegates to background subtraction before flat correction."""
    import importlib

    import soxspipe.commonutils.toolkit as toolkit

    recipeModule = importlib.import_module("soxspipe.recipes.base_recipe")
    recipe = _recipe(log)
    recipe.recipeSettings = {"subtract_background": True}
    inputFrame = _frame(20.0)
    expectedFrame = _frame(17.0)
    backgroundFrame = _frame(3.0)
    expectedProducts = product_table().iloc[0:0].copy()
    expectedProducts.loc[0] = ["product"] * len(expectedProducts.columns)
    received: dict[str, object] = {}

    class Background:
        """Deterministic background collaborator."""

        def __init__(self, **kwargs: object) -> None:
            received.update(kwargs)

        def subtract(self) -> tuple[CCDData, CCDData, pd.DataFrame]:
            return backgroundFrame, expectedFrame, expectedProducts

    monkeypatch.setattr(recipeModule, "subtract_background", Background)
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)

    result = recipe.detrend(inputFrame, master_bias=_frame(2.0), order_table="orders.fits")

    assert result is expectedFrame
    assert recipe.products is expectedProducts
    assert received["frame"].data[0, 0] == pytest.approx(18.0)
    assert received["orderTable"] == "orders.fits"
    assert received["productsTable"].empty


def test_clip_and_stack_clips_an_outlier_and_propagates_uncertainty(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Stacking clips rogue samples and retains CCDData metadata and uncertainty."""
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    recipe.kw = lambda keyword: f"ESO {keyword.replace('_', ' ')}" if keyword.startswith("DPR_") else keyword
    recipe.imageType = "OBJECT"
    recipe.verbose = True
    recipe.debug = False
    recipe.recipeSettings.update(
        {
            "stacked-clipping-sigma": 2.0,
            "stacked-clipping-iterations": 3,
            "frame-clipping-sigma": 5.0,
            "frame-clipping-iterations": 2,
        }
    )
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **_: None)
    frames = [_frame(10.0), _frame(10.0), _frame(10.0)]
    frames[2].data[1, 1] = 1000.0

    result = recipe.clip_and_stack(
        frames,
        "soxs_stare",
        ignore_input_masks=False,
        post_stack_clipping=True,
    )

    assert result.header[recipe.kw("DPR_TYPE")] == "OBJECT"
    assert result.data[1, 1] == pytest.approx(10.0)
    assert result.uncertainty.array[2, 2] == pytest.approx(1 / np.sqrt(3))


def test_clip_and_stack_short_circuits_single_frame_and_rejects_empty_input(log: Any) -> None:
    """The stack boundary retains a single frame and rejects an empty collection."""
    recipe = _recipe(log)
    frame = _frame(10.0)

    assert recipe.clip_and_stack([frame], "soxs_stare") is frame
    with pytest.raises(ValueError, match="No frames were sent"):
        recipe.clip_and_stack([], "soxs_stare")
