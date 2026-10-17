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

XSH_INSTRUMENT = "XSHOOTER"


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
    recipe.detectorParams = {"ron": 3.8 * u.electron}
    recipe.recipeSettings = {"subtract_background": False}
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
    # ONLY NIR ACCEPTS A BIAS AND A DARK TOGETHER (DY-1580)
    recipe.arm = "NIR"

    calibrated = recipe.detrend(
        _frame(20.0),
        master_bias=_frame(2.0),
        dark=_frame(3.0),
        master_flat=_frame(0.5),
    )

    np.testing.assert_allclose(calibrated.data, np.full((4, 4), 30.0))
    assert calibrated.data.dtype == np.float32
    assert not recipe.darkDetrendWarningIssued2


@pytest.mark.parametrize("instrument", ["SOXS", XSH_INSTRUMENT])
@pytest.mark.parametrize("arm", ["UVB", "VIS", "vis"])
def test_detrend_refuses_a_uvb_vis_frame_with_both_a_bias_and_a_dark(log: Any, arm: str, instrument: str) -> None:
    """A UVB/VIS master dark still carries the bias pedestal, so bias plus dark would subtract it twice (DY-1580)."""
    recipe = _recipe(log)
    recipe.arm = arm
    recipe.inst = instrument

    with pytest.raises(ValueError, match="still carry the bias pedestal"):
        recipe.detrend(_frame(20.0), master_bias=_frame(2.0), dark=_frame(3.0))


@pytest.mark.parametrize(
    ("arm", "masterBias", "dark", "expectedValue"),
    [
        pytest.param("NIR", 2.0, 3.0, 15.0, id="nir-bias-and-dark"),
        pytest.param("NIR", None, 3.0, 17.0, id="nir-dark-only"),
        pytest.param("VIS", 2.0, None, 18.0, id="vis-bias-only"),
        pytest.param("VIS", None, 3.0, 17.0, id="vis-dark-only"),
    ],
)
def test_detrend_accepts_every_other_bias_and_dark_combination(
    log: Any,
    arm: str,
    masterBias: float | None,
    dark: float | None,
    expectedValue: float,
) -> None:
    """NIR, and UVB/VIS with only one of bias or dark, detrend as before."""
    recipe = _recipe(log)
    recipe.arm = arm

    calibrated = recipe.detrend(
        _frame(20.0),
        master_bias=False if masterBias is None else _frame(masterBias),
        dark=False if dark is None else _frame(dark),
    )

    np.testing.assert_allclose(calibrated.data, np.full((4, 4), expectedValue))


def test_detrend_scales_a_different_exposure_dark_once(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """An X-Shooter different-exposure dark is scaled and emits its warning only once."""
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    recipe.inst = XSH_INSTRUMENT
    quicklookCalls: list[dict[str, object]] = []
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: quicklookCalls.append(kwargs))
    dark = _frame(2.0, exposure=5.0)

    calibrated = recipe.detrend(_frame(20.0, exposure=10.0), dark=dark)

    np.testing.assert_allclose(calibrated.data, np.full((4, 4), 16.0))
    assert recipe.darkDetrendWarningIssued2
    assert len(quicklookCalls) == 2
    assert quicklookCalls[0]["CCDObject"] is dark
    assert quicklookCalls[1]["CCDObject"] is calibrated


def _exptime_warnings(log: Any) -> list[str]:
    """Return every warning the recipe logged about differing exposure times."""
    return [message for level, message in log.messages if level == "warning" and "differing exposure-times" in message]


def test_detrend_does_not_warn_about_a_different_exposure_soxs_dark(log: Any) -> None:
    """SOXS darks are subtracted unscaled by design, so a different exposure time is not worth a warning (DY-1553)."""
    recipe = _recipe(log)
    recipe.arm = "NIR"
    recipe.inst = "SOXS"

    recipe.detrend(_frame(20.0, exposure=3.75), dark=_frame(2.0, exposure=5.0))

    assert _exptime_warnings(log) == []


def test_detrend_warns_once_about_a_different_exposure_xshooter_dark(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An X-Shooter dark of a different exposure time is scaled, which earns exactly one warning per recipe."""
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    recipe.arm = "NIR"
    recipe.inst = XSH_INSTRUMENT
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)

    recipe.detrend(_frame(20.0, exposure=10.0), dark=_frame(2.0, exposure=5.0))
    recipe.detrend(_frame(20.0, exposure=10.0), dark=_frame(2.0, exposure=5.0))

    assert len(_exptime_warnings(log)) == 1


def test_detrend_does_not_warn_about_xshooter_exposure_times_within_the_tolerance(log: Any) -> None:
    """X-Shooter exposure times that differ only by header float noise are equal, so they are not warned about."""
    recipe = _recipe(log)
    recipe.arm = "NIR"
    recipe.inst = XSH_INSTRUMENT

    recipe.detrend(_frame(20.0, exposure=5.0), dark=_frame(2.0, exposure=5.0 + 1e-6))

    assert _exptime_warnings(log) == []


def _spy_on_subtract_dark(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, object]]:
    """Record the keyword arguments of every real ``ccdproc.subtract_dark`` call."""
    import ccdproc

    realSubtractDark = ccdproc.subtract_dark
    calls: list[dict[str, object]] = []

    def spy(*args: object, **kwargs: object) -> CCDData:
        calls.append(kwargs)
        return realSubtractDark(*args, **kwargs)

    monkeypatch.setattr(ccdproc, "subtract_dark", spy)
    return calls


@pytest.mark.parametrize(
    ("instrument", "darkExposure", "frameExposure", "expectedScale", "expectedValue"),
    [
        # SOXS NEVER SCALES, WHETHER OR NOT THE EXPOSURE TIMES MATCH (DY-1263 RULING)
        pytest.param("SOXS", 0.4, 0.9, False, 18.0, id="soxs-sub-second-mismatch-unscaled"),
        pytest.param("SOXS", 2.0, 5.0, False, 18.0, id="soxs-whole-second-mismatch-unscaled"),
        pytest.param("SOXS", 5.0, 5.0, False, 18.0, id="soxs-equal-unscaled"),
        # X-SHOOTER SCALES WHENEVER THE EXPOSURE TIMES DIFFER: 20 - 2 * 0.9 / 0.4 = 15.5
        pytest.param(XSH_INSTRUMENT, 0.4, 0.9, True, 15.5, id="xsh-sub-second-mismatch-scaled"),
        pytest.param(XSH_INSTRUMENT, 5.0, 5.0, False, 18.0, id="xsh-equal-unscaled"),
        # JUST OVER THE HEADER-NOISE TOLERANCE IS A REAL MISMATCH: 20 - 2 * 5.0 / 5.01
        pytest.param(XSH_INSTRUMENT, 5.01, 5.0, True, 20.0 - 2.0 * 5.0 / 5.01, id="xsh-just-over-tolerance-scaled"),
    ],
)
def test_subtract_dark_frame_scales_only_differing_xshooter_exposures(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    instrument: str,
    darkExposure: float,
    frameExposure: float,
    expectedScale: bool,
    expectedValue: float,
) -> None:
    """SOXS darks are never scaled; X-Shooter darks are scaled only when the exposure times differ."""
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    recipe.inst = instrument
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)
    calls = _spy_on_subtract_dark(monkeypatch)
    frame = _frame(20.0, exposure=frameExposure)

    subtracted = recipe._subtract_dark_frame(frame, _frame(2.0, exposure=darkExposure), frame)

    assert len(calls) == 1
    assert calls[0].get("scale", False) is expectedScale
    np.testing.assert_allclose(subtracted.data, np.full((4, 4), expectedValue))


def test_subtract_dark_frame_treats_header_noise_as_equal_xshooter_exposures(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """X-Shooter exposure times that differ only by float header noise are not scaled."""
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)
    recipe.inst = XSH_INSTRUMENT
    calls = _spy_on_subtract_dark(monkeypatch)
    frame = _frame(20.0, exposure=5.0)

    recipe._subtract_dark_frame(frame, _frame(2.0, exposure=5.0 + 1e-6), frame)

    assert len(calls) == 1
    assert calls[0].get("scale", False) is False


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


def test_clip_and_stack_keeps_the_single_frame_flux_scale(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A stack is the mean of its frames with the first frame's EXPTIME, so counts per second match one frame.

    DY-1611: stare SOFs now stack every frame of an OB, and flux calibration divides by the stack's EXPTIME.
    """
    # ARRANGE
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    recipe.kw = lambda keyword: f"ESO {keyword.replace('_', ' ')}" if keyword.startswith("DPR_") else keyword
    recipe.imageType = "OBJECT"
    recipe.verbose = True
    recipe.debug = False
    recipe.recipeSettings.update({"stacked-clipping-sigma": 50.0, "stacked-clipping-iterations": 1})
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **_: None)
    frames = [_frame(10.0), _frame(20.0), _frame(30.0)]
    for frame in frames:
        frame.header["EXPTIME"] = 300.0

    # ACT
    result = recipe.clip_and_stack(frames, "soxs_stare", ignore_input_masks=False, post_stack_clipping=False)

    # ASSERT
    assert result.data[0, 0] == pytest.approx(20.0)
    assert result.header["EXPTIME"] == 300.0


def test_clip_and_stack_short_circuits_single_frame_and_rejects_empty_input(log: Any) -> None:
    """The stack boundary retains a single frame and rejects an empty collection."""
    recipe = _recipe(log)
    frame = _frame(10.0)

    assert recipe.clip_and_stack([frame], "soxs_stare") is frame
    with pytest.raises(ValueError, match="No frames were sent"):
        recipe.clip_and_stack([], "soxs_stare")


BIAS_PEDESTAL = 647.0
SKY_LEVEL = 4.0
READ_NOISE = 3.8
MASTER_BIAS_ERROR = 0.5


def _pedestal_frames(flatValue: float, stackSize: int) -> tuple[CCDData, CCDData, CCDData]:
    """Return a prepared sky frame, its master bias and its master flat.

    The frame sits on a bias pedestal and carries the uncertainty map that the preparation step gives it: the
    deviation of the raw (pedestal included) counts. ``stackSize`` frames of this kind have been mean-combined.
    """
    shape = (4, 4)
    recipe = _recipe(None)
    recipe.kw = lambda keyword: f"ESO {keyword.replace('_', ' ')}" if keyword.startswith("DPR_") else keyword
    recipe.detectorParams = {"ron": READ_NOISE * u.electron}
    rawFrame = _frame(BIAS_PEDESTAL + SKY_LEVEL)
    rawFrame = recipe._add_uncertainty_map(rawFrame)
    rawFrame.uncertainty = StdDevUncertainty(rawFrame.uncertainty.array / np.sqrt(stackSize), unit=u.electron)
    masterBias = _frame(BIAS_PEDESTAL)
    masterBias.uncertainty = StdDevUncertainty(np.full(shape, MASTER_BIAS_ERROR), unit=u.electron)
    masterFlat = _frame(flatValue)
    masterFlat.uncertainty = StdDevUncertainty(np.zeros(shape), unit=u.electron)
    return rawFrame, masterBias, masterFlat


@pytest.mark.parametrize(("flatValue", "stackSize"), [(1.0, 1), (0.5, 1), (1.0, 4)])
def test_detrend_error_excludes_the_bias_pedestal_from_the_photon_term(
    log: Any,
    flatValue: float,
    stackSize: int,
) -> None:
    """The error of a bias-subtracted sky pixel is its own photon noise, read noise and bias error, over the flat."""
    recipe = _recipe(log)
    recipe.detectorParams = {"ron": READ_NOISE * u.electron}
    rawFrame, masterBias, masterFlat = _pedestal_frames(flatValue, stackSize)

    calibrated = recipe.detrend(rawFrame, master_bias=masterBias, master_flat=masterFlat)

    expectedVariance = (SKY_LEVEL + READ_NOISE**2) / stackSize + MASTER_BIAS_ERROR**2
    expectedError = np.sqrt(expectedVariance) / flatValue
    np.testing.assert_allclose(calibrated.uncertainty.array, np.full((4, 4), expectedError), rtol=1e-4)
    np.testing.assert_allclose(calibrated.data, np.full((4, 4), SKY_LEVEL / flatValue), rtol=1e-4)


def test_detrend_leaves_the_input_uncertainty_untouched(log: Any) -> None:
    """Correcting the photon term must not change the frame that was handed in."""
    recipe = _recipe(log)
    recipe.detectorParams = {"ron": READ_NOISE * u.electron}
    rawFrame, masterBias, masterFlat = _pedestal_frames(1.0, 1)
    before = rawFrame.uncertainty.array.copy()

    recipe.detrend(rawFrame, master_bias=masterBias, master_flat=masterFlat)

    np.testing.assert_array_equal(rawFrame.uncertainty.array, before)


def _helper_inputs(rawValue: float, biasValue: float) -> tuple[base_recipe, CCDData, CCDData]:
    """Return a recipe, a 4x4 prepared frame with unit error and a master bias, for the pedestal helper."""
    recipe = _recipe(None)
    recipe.detectorParams = {"ron": READ_NOISE * u.electron}
    return recipe, _frame(rawValue), _frame(biasValue)


@pytest.mark.parametrize("nanIn", ["data", "bias"])
def test_pedestal_helper_leaves_the_error_of_a_non_finite_pixel_unchanged(nanIn: str) -> None:
    """A NaN in the frame or in the master bias keeps that pixel's error as it was."""
    recipe, frame, masterBias = _helper_inputs(BIAS_PEDESTAL + SKY_LEVEL, BIAS_PEDESTAL)
    (frame if nanIn == "data" else masterBias).data[1, 2] = np.nan

    corrected = recipe._remove_bias_pedestal_from_uncertainty(frame, masterBias)

    assert corrected.uncertainty.array[1, 2] == pytest.approx(1.0)
    expectedScale = (SKY_LEVEL + READ_NOISE**2) / (BIAS_PEDESTAL + SKY_LEVEL + READ_NOISE**2)
    assert corrected.uncertainty.array[0, 0] == pytest.approx(np.sqrt(expectedScale), rel=1e-5)


def test_pedestal_helper_keeps_only_read_noise_when_the_frame_is_below_the_bias() -> None:
    """Counts below the bias level carry no photon noise, so only the read-noise share of the variance remains."""
    rawValue = 100.0
    recipe, frame, masterBias = _helper_inputs(rawValue, BIAS_PEDESTAL)

    corrected = recipe._remove_bias_pedestal_from_uncertainty(frame, masterBias)

    expectedScale = READ_NOISE**2 / (rawValue + READ_NOISE**2)
    np.testing.assert_allclose(corrected.uncertainty.array, np.sqrt(expectedScale), rtol=1e-5)


def test_pedestal_helper_leaves_a_non_positive_pixel_unchanged() -> None:
    """A pixel at or below zero counts has a read-noise-only error already, so the scale is one."""
    recipe, frame, masterBias = _helper_inputs(-5.0, BIAS_PEDESTAL)

    corrected = recipe._remove_bias_pedestal_from_uncertainty(frame, masterBias)

    np.testing.assert_allclose(corrected.uncertainty.array, 1.0, rtol=1e-6)


def test_pedestal_helper_keeps_the_mask_and_an_electron_stddev_uncertainty() -> None:
    """The corrected copy keeps the input mask and a standard-deviation uncertainty in electrons."""
    recipe, frame, masterBias = _helper_inputs(BIAS_PEDESTAL + SKY_LEVEL, BIAS_PEDESTAL)
    frame.mask[2, 3] = True

    corrected = recipe._remove_bias_pedestal_from_uncertainty(frame, masterBias)

    assert corrected is not frame
    np.testing.assert_array_equal(corrected.mask, frame.mask)
    assert isinstance(corrected.uncertainty, StdDevUncertainty)
    assert corrected.uncertainty.unit == u.electron


def test_pedestal_helper_returns_a_frame_without_uncertainty_as_is() -> None:
    """A frame with no uncertainty map has nothing to correct."""
    recipe, frame, masterBias = _helper_inputs(BIAS_PEDESTAL + SKY_LEVEL, BIAS_PEDESTAL)
    frame.uncertainty = None

    assert recipe._remove_bias_pedestal_from_uncertainty(frame, masterBias) is frame
