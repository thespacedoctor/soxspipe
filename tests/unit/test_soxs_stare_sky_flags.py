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


def _frame(mask: np.ndarray | None) -> CCDData:
    return CCDData(
        np.ones(FRAME_SHAPE),
        unit=u.electron,
        mask=mask,
        uncertainty=StdDevUncertainty(np.ones(FRAME_SHAPE)),
    )


def _run_subtract_sky_and_capture_exclude_mask(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    objectFrame: CCDData,
    skySubtractionMask: np.ndarray | None = None,
    excludeMask: np.ndarray | None = None,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Run `_subtract_sky`; return the `excludeMask` the bad-pixel QC received and the subtractor's kwargs.

    When `skySubtractionMask` is None the subtractor double flags what the real one does: the object mask plus the
    `excludeMask` it was constructed with.
    """
    calls: dict[str, list] = {"excludeMask": []}
    constructorKwargs: list[dict[str, Any]] = []

    class FakeSubtractSky:
        def __init__(self, **kwargs: Any) -> None:
            constructorKwargs.append(kwargs)

        def subtract(self) -> tuple[Any, Any, Any, str, str]:
            if skySubtractionMask is not None:
                flags = skySubtractionMask
            else:
                received = constructorKwargs[0]["excludeMask"]
                flags = np.zeros(FRAME_SHAPE, dtype=bool) if received is None else received.copy()
                if objectFrame.mask is not None:
                    flags = flags | objectFrame.mask
            return _frame(flags), _frame(flags), None, "qc", "products"

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

    recipe._subtract_sky(objectFrame, "map.fits", "disp.fits", "orders.fits", excludeMask)

    [qcExcludeMask] = calls["excludeMask"]
    return qcExcludeMask, constructorKwargs


def test_the_bad_pixel_qc_does_not_count_pixels_that_only_the_sky_subtraction_flagged(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detectorMask = np.zeros(FRAME_SHAPE, dtype=bool)
    detectorMask[0, 0] = True
    skySubtractionMask = detectorMask.copy()
    skySubtractionMask[1:, :] = True

    excludeMask, _ = _run_subtract_sky_and_capture_exclude_mask(
        log, monkeypatch, _frame(detectorMask), skySubtractionMask
    )

    np.testing.assert_array_equal(excludeMask, skySubtractionMask & ~detectorMask)


def test_the_bad_pixel_qc_excludes_every_sky_flag_when_the_object_frame_has_no_mask(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    skySubtractionMask = np.zeros(FRAME_SHAPE, dtype=bool)
    skySubtractionMask[1:, :] = True

    excludeMask, _ = _run_subtract_sky_and_capture_exclude_mask(log, monkeypatch, _frame(None), skySubtractionMask)

    assert excludeMask.dtype == bool
    np.testing.assert_array_equal(excludeMask, skySubtractionMask)


def _sky_recipe(log: Any, monkeypatch: pytest.MonkeyPatch, *, subtractSky: bool, model: Any = "model") -> tuple:
    """Return a stare recipe whose sky subtraction is a recording double, and the list of its constructor kwargs."""
    constructorKwargs: list[dict[str, Any]] = []

    class RecordingSubtractSky:
        def __init__(self, **kwargs: Any) -> None:
            constructorKwargs.append(kwargs)

        def subtract(self) -> tuple[Any, Any, Any, str, str]:
            if model is None:
                return None, None, None, "qc", "products"
            return _frame(None), _frame(np.zeros(FRAME_SHAPE, dtype=bool)), None, "qc", "products"

    monkeypatch.setattr(stareModule, "subtract_sky", RecordingSubtractSky)
    monkeypatch.setattr(stareModule, "generic_quality_checks", lambda **_: "qc")
    monkeypatch.setattr(stareModule, "spectroscopic_image_quality_checks", lambda **_: "qc")
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.settings = {}
    recipe.recipeSettings = {"sky-subtraction": {"subtract_sky": subtractSky}}
    recipe.subtractSky = subtractSky
    recipe.qc = "qc"
    recipe.products = "products"
    recipe.sofName = "stare"
    recipe.recipeName = "soxs-stare"
    recipe.startNightDate = "2024-01-02"
    recipe.debug = False
    recipe._write_sky_products = lambda *_: "product.fits"
    return recipe, constructorKwargs


def _overlap_mask() -> np.ndarray:
    excludeMask = np.zeros(FRAME_SHAPE, dtype=bool)
    excludeMask[2, :] = True
    return excludeMask


def test_the_overlap_mask_reaches_the_sky_subtraction_constructor(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    recipe, constructorKwargs = _sky_recipe(log, monkeypatch, subtractSky=True)
    excludeMask = _overlap_mask()

    recipe._subtract_sky(_frame(None), "map.fits", "disp.fits", "orders.fits", excludeMask)

    [kwargs] = constructorKwargs
    assert kwargs["excludeMask"] is excludeMask


def test_the_sky_subtraction_constructor_gets_no_mask_by_default(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    recipe, constructorKwargs = _sky_recipe(log, monkeypatch, subtractSky=True)

    recipe._subtract_sky(_frame(None), "map.fits", "disp.fits", "orders.fits")

    [kwargs] = constructorKwargs
    assert kwargs["excludeMask"] is None


def test_the_master_flat_reaches_the_sky_subtraction_constructor(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    recipe, constructorKwargs = _sky_recipe(log, monkeypatch, subtractSky=True)
    masterFlat = _frame(None)

    recipe._subtract_sky(_frame(None), "map.fits", "disp.fits", "orders.fits", None, masterFlat)

    [kwargs] = constructorKwargs
    assert kwargs["masterFlat"] is masterFlat


def test_the_sky_subtraction_constructor_gets_no_master_flat_by_default(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    recipe, constructorKwargs = _sky_recipe(log, monkeypatch, subtractSky=True)

    recipe._subtract_sky(_frame(None), "map.fits", "disp.fits", "orders.fits")

    [kwargs] = constructorKwargs
    assert kwargs["masterFlat"] is None


def test_the_bad_pixel_qc_does_not_count_pixels_excluded_by_the_overlap_mask(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # THE DETECTOR MASK HAS ONE BAD PIXEL OF ITS OWN; THE OVERLAP MASK IS PASSED IN AND FLAGGED BY THE SUBTRACTOR DOUBLE
    detectorMask = np.zeros(FRAME_SHAPE, dtype=bool)
    detectorMask[0, 0] = True
    excludeMask = _overlap_mask()

    qcExcludeMask, constructorKwargs = _run_subtract_sky_and_capture_exclude_mask(
        log, monkeypatch, _frame(detectorMask), excludeMask=excludeMask
    )

    assert constructorKwargs[0]["excludeMask"] is excludeMask
    # EXACTLY THE OVERLAP-MASKED PIXELS DROP OUT OF THE BAD-PIXEL COUNT, AND THE DETECTOR'S BAD PIXEL STAYS IN IT
    np.testing.assert_array_equal(qcExcludeMask, excludeMask)
    assert excludeMask.any()
    assert not qcExcludeMask[0, 0]


@pytest.mark.parametrize("subtractSky", [False, True], ids=["sky-off", "no-model"])
def test_without_a_sky_model_the_returned_frame_is_a_copy_flagged_with_the_overlap_mask(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    subtractSky: bool,
) -> None:
    recipe, _ = _sky_recipe(log, monkeypatch, subtractSky=subtractSky, model=None)
    detectorMask = np.zeros(FRAME_SHAPE, dtype=bool)
    detectorMask[0, 0] = True
    combined = _frame(detectorMask)
    excludeMask = _overlap_mask()

    model, skySubtracted, productPath = recipe._subtract_sky(
        combined, "map.fits", "disp.fits", "orders.fits", excludeMask
    )

    assert model is False
    assert productPath is None
    assert skySubtracted is not combined
    np.testing.assert_array_equal(skySubtracted.mask, detectorMask | excludeMask)
    np.testing.assert_array_equal(skySubtracted.data, combined.data)
    np.testing.assert_array_equal(skySubtracted.uncertainty.array, combined.uncertainty.array)
    # THE INPUT FRAME IS NOT MUTATED
    np.testing.assert_array_equal(combined.mask, detectorMask)


@pytest.mark.parametrize("subtractSky", [False, True], ids=["sky-off", "no-model"])
def test_without_a_sky_model_a_frame_with_no_mask_is_flagged_with_only_the_overlap_mask(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    subtractSky: bool,
) -> None:
    recipe, _ = _sky_recipe(log, monkeypatch, subtractSky=subtractSky, model=None)
    combined = _frame(None)
    excludeMask = _overlap_mask()

    _, skySubtracted, _ = recipe._subtract_sky(combined, "map.fits", "disp.fits", "orders.fits", excludeMask)

    np.testing.assert_array_equal(skySubtracted.mask, excludeMask)
    assert combined.mask is None


@pytest.mark.parametrize("subtractSky", [False, True], ids=["sky-off", "no-model"])
def test_without_a_sky_model_and_without_an_overlap_mask_the_input_frame_is_returned_itself(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    subtractSky: bool,
) -> None:
    recipe, _ = _sky_recipe(log, monkeypatch, subtractSky=subtractSky, model=None)
    combined = _frame(None)

    _, skySubtracted, _ = recipe._subtract_sky(combined, "map.fits", "disp.fits", "orders.fits")

    assert skySubtracted is combined


def _binning_recipe(log: Any, arm: str = "VIS") -> soxs_stare:
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = arm
    recipe.kw = lambda keyword: keyword
    return recipe


def _binned_frame(**header: object) -> CCDData:
    frame = _frame(None)
    frame.header.update(header)
    return frame


def test_the_binning_is_read_when_both_header_keys_are_present(log: Any) -> None:
    frame = _binned_frame(WIN_BINX=2, WIN_BINY=4)

    assert _binning_recipe(log)._read_frame_binning(frame) == (2, 4)


@pytest.mark.parametrize(
    "header",
    [{"WIN_BINX": 2}, {"WIN_BINY": 2}, {"WIN_BINX": 2, "WIN_BINY": "two"}, {}],
    ids=["only-x", "only-y", "not-an-integer", "neither"],
)
def test_the_binning_falls_back_to_one_by_one_with_a_warning_when_a_header_key_is_unusable(
    log: Any, header: dict[str, object]
) -> None:
    result = _binning_recipe(log)._read_frame_binning(_binned_frame(**header))

    assert result == (1, 1)
    assert [text for level, text in log.messages if level == "warning" and "binning" in text]


def test_the_nir_binning_is_one_by_one_without_a_warning(log: Any) -> None:
    result = _binning_recipe(log, arm="NIR")._read_frame_binning(_binned_frame(WIN_BINX=2, WIN_BINY=2))

    assert result == (1, 1)
    assert not [level for level, _ in log.messages if level == "warning"]
