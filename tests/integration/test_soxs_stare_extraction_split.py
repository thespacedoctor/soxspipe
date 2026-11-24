"""The stare tail is split into an extraction half and a finishing half, so a nod position can extract alone (DY-2290).

`_extract_stack` runs the sky fit and the flat-fielded extraction of one stack. `_finish_extraction` flux-calibrates or
builds the response curve, writes the ASCII spectrum and draws the plots. `_reduce_stacked_object` is the two in turn,
and its output is pinned by the stare characterization tests.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_stare
from soxspipe.recipes.nod_stare_stacking import reduction_inputs, response_inputs, stack_extraction
from tests.factories import synthetic_ccd
from tests.integration.test_soxs_stare_reduction_characterization import (
    RECIPE_READ_NOISE,
    _patch_stare,
    _stare_recipe,
)

pytestmark = pytest.mark.integration

HANDOFF = object()
NO_TRACE_MESSAGE = "No object trace was found in the standard-star frame"
SKY_SEEDS = (530, 531, 532)


def _inputs(responseFunctionPath: Any = False) -> reduction_inputs:
    return reduction_inputs(False, False, False, "order.fits", "disp.fits", "map.fits", responseFunctionPath)


class _extractor:
    """An extraction double that records its keyword arguments and reports a spill, a fallback and a handoff."""

    built: list[dict[str, Any]] = []
    slitEdgeOrders: list[int] = []
    traceFallback = False
    extractionNotes = [("WARNING", "an order was skipped")]
    flatPassHandoff = HANDOFF

    def __init__(self, **kwargs: Any) -> None:
        type(self).built.append(kwargs)

    def extract(self) -> tuple[Any, Any, pd.DataFrame, dict[str, int], str]:
        name = "notflat.fits" if self.built[-1].get("notFlattened") else "extracted.fits"
        spectrum = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})
        return self.built[-1]["qcTable"], self.built[-1]["productsTable"], spectrum, {"10-11": 11}, name


def _extractor_class(*, slitEdgeOrders: list[int] | None = None, traceFallback: bool = False) -> type[_extractor]:
    class extractor(_extractor):
        pass

    extractor.built = []
    extractor.slitEdgeOrders = slitEdgeOrders or []
    extractor.traceFallback = traceFallback
    return extractor


def _recipe(
    log: Any,
    tmpPath: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    isStandard: bool = False,
    slitEdgeOrders: list[int] | None = None,
    traceFallback: bool = False,
    skyResult: tuple[Any, ...] | None = None,
) -> tuple[soxs_stare, dict[str, list[dict[str, object]]], type[_extractor]]:
    recipe = _stare_recipe(log, tmpPath)
    recipe.generateReponseCurve = isStandard
    skyResult = skyResult or tuple(synthetic_ccd(seed=seed, prepared=True) for seed in SKY_SEEDS)
    captured = _patch_stare(recipe, monkeypatch, tmpPath, skyResult=skyResult)
    extractor = _extractor_class(slitEdgeOrders=slitEdgeOrders, traceFallback=traceFallback)
    monkeypatch.setattr(import_module("soxspipe.commonutils"), "horne_extraction", extractor)
    recipe.scienceHeaders = []
    recipe.subtractSky = True
    recipe.slitIlluminationOffsets = None
    return recipe, captured, extractor


def _extraction(**overrides: Any) -> stack_extraction:
    fields: dict[str, Any] = {
        "productPath": "product.fits",
        "combinedObject": synthetic_ccd(seed=3, prepared=True),
        "mergedSpectrum": pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]}),
        "orderJoins": {"10-11": 11},
        "extractionPath": "extracted.fits",
        "responseInputs": None,
        "notFlatPath": None,
        "notFlatSpectrum": None,
        "traceFallback": False,
        "slitEdgeSpill": False,
        "extractionNotes": [],
    }
    return stack_extraction(**{**fields, **overrides})


def _spy_finish(recipe: soxs_stare, monkeypatch: pytest.MonkeyPatch) -> dict[str, list[tuple[Any, ...]]]:
    calls: dict[str, list[tuple[Any, ...]]] = {"fluxCal": [], "response": [], "compute": [], "ascii": [], "plot": []}
    monkeypatch.setattr(
        recipe, "_flux_calibrate_spectrum", lambda *args: calls["fluxCal"].append(args) or "fluxcal.fits"
    )
    monkeypatch.setattr(
        recipe, "_generate_response_curve", lambda *args, **kwargs: calls["response"].append((args, kwargs)) or False
    )
    monkeypatch.setattr(recipe, "_compute_response", lambda *args: calls["compute"].append(args) or False)
    monkeypatch.setattr(recipe, "_write_ascii_spectrum", lambda *args: calls["ascii"].append(args))
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: calls["plot"].append(args))
    return calls


def test_the_extraction_half_returns_the_extraction_and_does_not_finish_it(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, extractor = _recipe(log, tmp_path, monkeypatch, slitEdgeOrders=[11])
    calls = _spy_finish(recipe, monkeypatch)
    stack = synthetic_ccd(seed=3, prepared=True)

    # ACT
    extraction = recipe._extract_stack(stack, [stack], False, _inputs(), stack)

    # ASSERT
    assert isinstance(extraction, stack_extraction)
    assert extraction.combinedObject is stack
    assert extraction.productPath == str(tmp_path / "OBJECT_VIS_SKYSUB_RESIDUALS.fits")
    assert extraction.extractionPath == "extracted.fits"
    assert extraction.orderJoins == {"10-11": 11}
    assert extraction.mergedSpectrum["WAVE"].tolist() == [500.0]
    assert extraction.slitEdgeSpill is True
    assert extraction.traceFallback is False
    assert extraction.extractionNotes == [("WARNING", "an order was skipped")]
    assert len(extractor.built) == 1
    assert extractor.built[0]["readNoise"] == RECIPE_READ_NOISE
    assert calls == {"fluxCal": [], "response": [], "compute": [], "ascii": [], "plot": []}
    assert recipe.dateObs == stack.header["DATE-OBS"]


def test_a_science_stack_keeps_no_inputs_for_a_response_curve(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch)
    stack = synthetic_ccd(seed=3, prepared=True)

    # ACT
    extraction = recipe._extract_stack(stack, [stack], False, _inputs(), stack)

    # ASSERT
    assert extraction.responseInputs is None
    assert extraction.notFlatPath is None


def test_a_standard_stack_keeps_the_inputs_of_its_not_flattened_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    skyResult = tuple(synthetic_ccd(seed=seed, prepared=True) for seed in SKY_SEEDS)
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch, isStandard=True, skyResult=skyResult)
    stack = synthetic_ccd(seed=3, prepared=True)

    # ACT
    extraction = recipe._extract_stack(stack, [stack], False, _inputs(), stack)

    # ASSERT
    inputsOfResponse = extraction.responseInputs
    assert isinstance(inputsOfResponse, response_inputs)
    assert (inputsOfResponse.twoDMap, inputsOfResponse.dispMap) == ("map.fits", "disp.fits")
    assert inputsOfResponse.flatPass is HANDOFF
    assert inputsOfResponse.skyModel is skyResult[0]
    assert inputsOfResponse.unflattened is skyResult[1]


def test_a_standard_with_a_response_table_or_without_a_trace_keeps_no_response_inputs(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    (tmp_path / "table").mkdir()
    (tmp_path / "fallback").mkdir()
    withTable, _, _ = _recipe(log, tmp_path / "table", monkeypatch, isStandard=True)
    withFallback, _, _ = _recipe(log, tmp_path / "fallback", monkeypatch, isStandard=True, traceFallback=True)
    stack = synthetic_ccd(seed=3, prepared=True)

    # ACT
    tabled = withTable._extract_stack(stack, [stack], False, _inputs("resp.fits"), stack)
    fallback = withFallback._extract_stack(stack, [stack], False, _inputs(), stack)

    # ASSERT
    assert tabled.responseInputs is None
    assert fallback.responseInputs is None
    assert fallback.traceFallback is True


def test_the_finishing_half_flux_calibrates_with_a_response_table_then_writes_and_plots(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch)
    calls = _spy_finish(recipe, monkeypatch)
    extraction = _extraction(extractionNotes=[("WARNING", "note")])
    frames = [extraction.combinedObject]

    # ACT
    failure = recipe._finish_extraction(extraction, frames, "resp.fits")

    # ASSERT
    assert failure is False
    assert calls["fluxCal"] == [("resp.fits", extraction.mergedSpectrum, extraction.combinedObject, frames)]
    assert calls["ascii"] == [("extracted.fits", True)]
    assert calls["plot"] == [(extraction.mergedSpectrum, {"10-11": 11}, "fluxcal.fits", [("WARNING", "note")], False)]


def test_the_finishing_half_fails_a_standard_that_found_no_trace_and_notes_why(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch, isStandard=True)
    calls = _spy_finish(recipe, monkeypatch)
    extraction = _extraction(traceFallback=True)

    # ACT
    failure = recipe._finish_extraction(extraction, [extraction.combinedObject], False)

    # ASSERT
    assert isinstance(failure, str)
    assert failure.startswith(NO_TRACE_MESSAGE)
    assert calls["response"] == []
    assert calls["compute"] == []
    assert calls["plot"][0][3] == [("ERROR", failure)]
    assert calls["plot"][0][4] is True


def test_the_finishing_half_builds_the_response_from_the_inputs_of_a_stack_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch, isStandard=True)
    calls = _spy_finish(recipe, monkeypatch)
    inputsOfResponse = response_inputs("UNFLAT", "SKY", "map.fits", "disp.fits", HANDOFF)
    extraction = _extraction(responseInputs=inputsOfResponse)
    frames = [extraction.combinedObject]

    # ACT
    failure = recipe._finish_extraction(extraction, frames, False)

    # ASSERT
    assert failure is False
    assert calls["response"] == [
        (("UNFLAT", "SKY", "map.fits", "disp.fits", "extracted.fits", {"10-11": 11}, HANDOFF), {"objectFrames": frames})
    ]


def test_the_finishing_half_builds_the_response_from_a_not_flattened_extraction_already_made(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch, isStandard=True)
    calls = _spy_finish(recipe, monkeypatch)
    extraction = _extraction(notFlatPath="notflat.fits")
    frames = [extraction.combinedObject]

    # ACT
    failure = recipe._finish_extraction(extraction, frames, False)

    # ASSERT
    assert failure is False
    assert calls["response"] == []
    assert calls["compute"] == [("extracted.fits", "notflat.fits", {"10-11": 11}, frames)]


def test_a_standard_that_spills_over_the_slit_edge_fails_even_when_the_response_passes(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch, isStandard=True)
    _spy_finish(recipe, monkeypatch)
    extraction = _extraction(notFlatPath="notflat.fits", slitEdgeSpill=True)

    # ACT
    failure = recipe._finish_extraction(extraction, [extraction.combinedObject], False)

    # ASSERT
    assert failure is True
    assert any(level == "error" and "spills over the slit edge" in message for level, message in log.messages)


def test_the_not_flattened_extraction_is_handed_the_flat_pass_and_returns_its_path_and_spectrum(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, extractor = _recipe(log, tmp_path, monkeypatch, isStandard=True)

    # ACT
    path, spectrum = recipe._extract_not_flattened("UNFLAT", "SKY", "map.fits", "disp.fits", HANDOFF)

    # ASSERT
    assert path == "notflat.fits"
    assert spectrum["WAVE"].tolist() == [500.0]
    built = extractor.built[0]
    assert built["notFlattened"] is True
    assert built["readNoise"] == RECIPE_READ_NOISE
    assert built["flatPass"] is HANDOFF
    assert built["skySubtractedFrame"] == "UNFLAT"
    assert built["subtractedFrame"] == "SKY"


def test_the_response_is_computed_from_the_flat_and_not_flattened_extractions(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _recipe(log, tmp_path, monkeypatch, isStandard=True)
    received: dict[str, Any] = {}

    class response:
        def __init__(self, **kwargs: Any) -> None:
            received.update(kwargs)

        def get(self) -> tuple[Any, Any, str]:
            return recipe.qc, recipe.products, "failed"

    airmass = import_module("soxspipe.commonutils.sequence_airmass")
    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", response)
    monkeypatch.setattr(airmass, "require_equal_exposure_times", lambda *_: None)
    monkeypatch.setattr(airmass, "mean_airmass", lambda *_: 1.2)
    recipe.kw = import_module("soxspipe.commonutils").keyword_lookup(log=log, settings=recipe.settings).get
    frame = synthetic_ccd(seed=1, prepared=True)

    # ACT
    failure = recipe._compute_response("flat.fits", "notflat.fits", {"10-11": 11}, [frame])

    # ASSERT
    assert failure == "failed"
    assert received["stdExtractionPath"] == "flat.fits"
    assert received["stdNotFlatExtractionPath"] == "notflat.fits"
    assert received["orderJoins"] == {"10-11": 11}
    assert received["airmass"] == 1.2
    assert received["sofName"] == recipe.productStem
