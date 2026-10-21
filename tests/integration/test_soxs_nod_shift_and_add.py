"""The stacked nod path extracts the shift-and-add of the A-B and B-A differences once per pass (DY-1674).

A science nod makes one Horne extraction. A standard star makes two: the flat-fielded pass and the not-flattened pass
that feeds the response curve. When the object position cannot be measured at both nod positions, each pass falls
back to the A-B and the B-A extractions. The per-cycle path is not touched.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.nod_shift_and_add import shift_and_add_plan
from soxspipe.commonutils.slit_alignment import SHIFT_COLUMNS, build_two_d_map_images, slit_alignment_result
from tests.integration.test_soxs_nod_cycle_characterization import (
    _ab_frames,
    _cycle_recipe,
    _patch_extraction_collaborators,
)
from tests.integration.test_soxs_nod_reduction_characterization import (
    SINGLE_PAIR,
    TWO_LOCATIONS,
    _extract_each_cycle,
    _nod_recipe,
    _patch,
    _StubResponse,
)

pytestmark = pytest.mark.integration

SHAPE = (32, 32)
NAN = float("nan")
# ERRORS OF THE PREPARED SYNTHETIC FRAMES, AND THE ERROR OF THE MEAN OF TWO INDEPENDENT DIFFERENCES OF THEM
FRAME_SIGMA = 1.5
DIFFERENCE_SIGMA = FRAME_SIGMA * np.sqrt(2)
COMBINED_SIGMA = DIFFERENCE_SIGMA / np.sqrt(2)


def _small_map_images():
    """Return the map of one order whose slit coordinate runs along the columns, 0.2 arcsec per pixel."""
    rows, columns = np.mgrid[0 : SHAPE[0], 0 : SHAPE[1]].astype(float)
    return build_two_d_map_images(500.0 + 0.02 * rows, 0.2 * (columns - 16.0), np.ones(SHAPE))


def _install_alignment(monkeypatch: pytest.MonkeyPatch, positionA: float, positionB: float) -> list[dict[str, Any]]:
    """Make the slit alignment report an object at ``positionA`` for the A group and ``positionB`` for the B group."""
    slitAlignment = import_module("soxspipe.commonutils.slit_alignment")
    received: list[dict[str, Any]] = []
    mapImages = _small_map_images()

    def align(**kwargs: Any) -> list[slit_alignment_result]:
        received.append(kwargs)
        results = []
        for group, position in zip(kwargs["groups"], (positionA, positionB), strict=True):
            shifts = pd.DataFrame(
                {
                    "frame": group.frameNames,
                    "mjd": [frame.header["MJD-OBS"] for frame in group.frames],
                    "order": 1,
                    "position_arcsec": position,
                    "shift_arcsec": 0.0,
                    "is_fallback": False,
                    "is_reference": True,
                }
            )[SHIFT_COLUMNS]
            results.append(slit_alignment_result(tuple(group.frames), shifts, 0, mapImages))
        return results

    monkeypatch.setattr(slitAlignment, "align_frame_groups_along_slit", align)
    return received


def _set_aperture(recipe: Any, slitLength: float = 10) -> None:
    """Give the recipe the settings the throw guard reads: the Horne aperture length and the dispersion axis."""
    recipe.recipeSettings = {**recipe.recipeSettings, "horne-extraction-slit-length": slitLength}
    recipe.detectorParams = {"dispersion-axis": "y"}


def _real_cycle(recipe: Any, monkeypatch: pytest.MonkeyPatch) -> list[dict[str, object]]:
    """Run the real cycle method with a recording Horne extractor, and return the extractor constructions."""
    monkeypatch.delattr(recipe, "process_single_ab_nodding_cycle")
    recipe.recipeSettings = {**recipe.recipeSettings, "save_single_frame_extractions": False}
    recipe.turnOffMP = True
    return _patch_extraction_collaborators(monkeypatch, recipe)[0]


def _stacked_run(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    isStandard: bool,
    positions: tuple[float, float] | None = (1.5, -1.5),
    debug: bool = False,
    slitLength: float = 10,
):
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    recipe.debug = debug
    # THE STUB MAP HAS 0.2 ARCSEC PIXELS, SO THE DEFAULT 10 PIXEL APERTURE IS 2 ARCSEC AND THE 3 ARCSEC THROW CLEARS IT
    _set_aperture(recipe, slitLength)
    recipe.generateReponseCurve = isStandard
    _, captured = _patch(recipe, monkeypatch, tmp_path)
    if positions is not None:
        _install_alignment(monkeypatch, *positions)
    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", _StubResponse)
    constructions = _real_cycle(recipe, monkeypatch)
    recipe.produce_product()
    return recipe, constructions, captured


def test_a_science_nod_with_a_measured_position_is_extracted_once(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    _, constructions, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT: ONE EXTRACTION, AND ONE SPECTRUM TO STACK
    assert len(constructions) == 1
    assert len(captured["stack_extractions"]) == 1
    spectra = captured["stack_extractions"][0]["args"][0]
    assert len(spectra) == 1
    assert captured["stack_extractions"][0]["orderJoins"] == {"10-11": 1}


def test_the_one_extraction_reads_the_combined_frame_with_the_error_of_a_mean(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    _, constructions, _ = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT: THE FRAME IS NEITHER THE A-B DIFFERENCE NOR THE B-A DIFFERENCE. WHERE THE SHIFT STAYS ON THE DETECTOR
    # ITS ERROR IS THE ERROR OF THE MEAN OF TWO DIFFERENCES
    extraction = constructions[0]
    combined = extraction["skySubtractedFrame"]
    assert combined.uncertainty.array[16, 8] == pytest.approx(COMBINED_SIGMA, rel=1e-5)
    assert extraction["unflattenedFrame"].uncertainty.array[16, 8] == pytest.approx(COMBINED_SIGMA, rel=1e-5)
    assert combined.header["ESO SEQ CUMOFF Y"] == 3.0
    assert extraction["notFlattened"] is False


def test_the_sky_frame_is_the_b_frame_averaged_with_the_a_frame_shifted_the_same_way(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    _, constructions, _ = _stacked_run(log, tmp_path, monkeypatch, isStandard=True)

    # ASSERT: BOTH PASSES, WITH THE ERROR OF THE MEAN OF TWO FRAMES AND THE B HEADER
    for extraction in constructions:
        sky = extraction["subtractedFrame"]
        assert sky.uncertainty.array[16, 8] == pytest.approx(FRAME_SIGMA / np.sqrt(2), rel=1e-5)
        assert sky.header["ESO SEQ CUMOFF Y"] == -3.0


def test_a_standard_star_is_extracted_twice_once_flattened_and_once_not(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    _, constructions, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=True)

    # ASSERT
    assert [construction["notFlattened"] for construction in constructions] == [False, True]
    assert len(captured["stack_extractions"]) == 2
    assert [len(call["args"][0]) for call in captured["stack_extractions"]] == [1, 1]


def test_the_ab_and_ba_differences_are_still_written_and_the_combined_frame_only_in_debug(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT
    assert (tmp_path / "OBJECT_VIS_AB_1.fits").exists()
    assert (tmp_path / "OBJECT_VIS_BA_1.fits").exists()
    assert not list(tmp_path.glob("*COMBINED*"))


def test_debug_writes_the_combined_frame_as_a_product(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    recipe, _, _ = _stacked_run(log, tmp_path, monkeypatch, isStandard=False, debug=True)

    # ASSERT
    assert (tmp_path / "OBJECT_VIS_AB_COMBINED_1.fits").exists()
    labels = recipe.products["product_label"].tolist()
    assert labels.count("AB_COMBINED") == 1


def _qc_values(recipe: Any) -> dict[str, object]:
    return dict(zip(recipe.qc["qc_name"], recipe.qc["qc_value"], strict=True))


def test_the_qcs_record_a_combined_pass_and_no_median_filled_shift(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    recipe, _, _ = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT
    values = _qc_values(recipe)
    assert int(values["NOD COMBINED"]) == 1
    assert int(values["NOD COMBINE FALLBACK NUM"]) == 0


def test_a_throw_shorter_than_the_aperture_extracts_the_differences_separately_and_warns(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT: A 20 PIXEL APERTURE IS 4 ARCSEC AND THE THROW IS 3 ARCSEC
    recipe, constructions, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=False, slitLength=20)

    # ASSERT
    assert len(constructions) == 2
    assert [len(call["args"][0]) for call in captured["stack_extractions"]] == [2]
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("shorter than the extraction aperture" in message for message in warnings)
    values = _qc_values(recipe)
    assert int(values["NOD COMBINED"]) == 0
    assert int(values["NOD COMBINE FALLBACK NUM"]) == 0


def test_without_a_measured_position_the_qcs_record_no_combination(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    recipe, _, _ = _stacked_run(log, tmp_path, monkeypatch, isStandard=False, positions=None)

    # ASSERT
    values = _qc_values(recipe)
    assert int(values["NOD COMBINED"]) == 0
    assert int(values["NOD COMBINE FALLBACK NUM"]) == 0


@pytest.mark.parametrize(("isStandard", "expectedExtractions"), [(False, 2), (True, 4)])
def test_without_a_measured_position_each_pass_extracts_the_ab_and_ba_differences_and_warns(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    isStandard: bool,
    expectedExtractions: int,
) -> None:
    # ARRANGE AND ACT: THE ALIGNMENT STUB OF THE SHARED COLLABORATORS REPORTS NO POSITION
    _, constructions, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=isStandard, positions=None)

    # ASSERT
    assert len(constructions) == expectedExtractions
    assert {len(call["args"][0]) for call in captured["stack_extractions"]} == {2}
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("extracted separately" in message for message in warnings)


def test_the_alignment_measures_groups_of_one_frame(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    _set_aperture(recipe)
    _patch(recipe, monkeypatch, tmp_path)
    received = _install_alignment(monkeypatch, 1.5, -1.5)
    _real_cycle(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert received[0]["measureSingleFrames"] is True


def test_both_passes_of_a_standard_receive_the_same_plan(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    recipe.generateReponseCurve = True
    _set_aperture(recipe)
    _, captured = _patch(recipe, monkeypatch, tmp_path)
    _install_alignment(monkeypatch, 1.5, -1.5)
    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", _StubResponse)

    # ACT
    recipe.produce_product()

    # ASSERT: THE PLAN IS THE B POSITION MINUS THE A POSITION, FOR THE ONE ORDER OF THE STUB MAP
    firstCycle, secondCycle = captured["extract_cycle"]
    assert firstCycle["shiftAndAdd"] is secondCycle["shiftAndAdd"]
    assert firstCycle["shiftAndAdd"].orderShifts == pytest.approx({1: -3.0})


def test_the_per_cycle_path_still_extracts_the_ab_and_ba_differences_of_each_cycle(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    _extract_each_cycle(recipe)
    _, captured = _patch(recipe, monkeypatch, tmp_path)
    constructions = _real_cycle(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT: TWO CYCLES, EACH WITH AN A-B AND A B-A EXTRACTION, STACKED AS TWO SPECTRA
    assert len(constructions) == 4
    assert [len(call["args"][0]) for call in captured["stack_extractions"]] == [2]
    assert not captured["align"]
    assert not list(Path(tmp_path).glob("*COMBINED*"))


def _flat_scaling_plan_run(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, notFlattened: bool
) -> dict[str, Any]:
    """Run one cycle with a plan and a master flat that doubles every frame, and return the extractor arguments."""
    recipe = _cycle_recipe(log, tmp_path)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    monkeypatch.setattr(recipe, "detrend", lambda **kwargs: kwargs["inputFrame"].multiply(2.0))
    plan = shift_and_add_plan({}, SimpleNamespace(order=np.zeros((3, 3))))
    aFrame, bFrame = _ab_frames()

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=1,
        orderTablePath="orders.fits",
        notFlattened=notFlattened,
        masterFlat=object(),
        shiftAndAdd=plan,
    )

    assert len(extractCalls) == 1
    return extractCalls[0]


def test_a_flat_fielded_pass_extracts_the_combination_of_the_flat_fielded_differences(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    arguments = _flat_scaling_plan_run(log, tmp_path, monkeypatch, notFlattened=False)

    # ASSERT: THE FLAT DOUBLES THE ERROR OF THE COMBINATION THAT THE TRACE IS FOUND ON
    assert arguments["unflattenedFrame"].uncertainty.array[1, 1] == pytest.approx(FRAME_SIGMA)
    assert arguments["skySubtractedFrame"].uncertainty.array[1, 1] == pytest.approx(2 * FRAME_SIGMA)
    # THE SKY FRAME IS BUILT FROM THE FLAT-FIELDED A AND B FRAMES, WHICH THE STUB FLAT DOUBLES
    assert arguments["subtractedFrame"].uncertainty.array[1, 1] == pytest.approx(2 * FRAME_SIGMA / np.sqrt(2))


def test_a_not_flattened_pass_extracts_the_one_unflattened_combination(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    arguments = _flat_scaling_plan_run(log, tmp_path, monkeypatch, notFlattened=True)

    # ASSERT
    assert arguments["skySubtractedFrame"] is arguments["unflattenedFrame"]
    assert arguments["skySubtractedFrame"].uncertainty.array[1, 1] == pytest.approx(FRAME_SIGMA)
    assert arguments["subtractedFrame"].uncertainty.array[1, 1] == pytest.approx(FRAME_SIGMA / np.sqrt(2))
    assert arguments["notFlattened"] is True


def test_a_cycle_without_a_plan_still_makes_two_extractions(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _cycle_recipe(log, tmp_path)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()

    # ACT
    mergedA, mergedB, _ = recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits"
    )

    # ASSERT
    assert len(extractCalls) == 2
    assert mergedA is not False
    assert mergedB is not False


@pytest.mark.parametrize(
    ("header", "slitLength", "dispersionAxis", "expected"),
    [
        ({}, 15, "y", 15.0),
        ({"WIN_BINX": 2, "WIN_BINY": 4}, 15, "y", 3.75),
        ({"WIN_BINX": "x", "WIN_BINY": 4}, 15, "y", 15.0),
    ],
)
def test_the_aperture_uses_the_binning_in_the_header_and_one_by_one_without_it(
    log: Any, tmp_path: Path, header: dict[str, object], slitLength: float, dispersionAxis: str, expected: float
) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    recipe.recipeSettings = {**recipe.recipeSettings, "horne-extraction-slit-length": slitLength}
    recipe.detectorParams = {"dispersion-axis": dispersionAxis}
    frame = SimpleNamespace(header=header)

    # ACT
    aperture = recipe._horne_aperture_pixels(frame)

    # ASSERT
    assert aperture == pytest.approx(expected)
