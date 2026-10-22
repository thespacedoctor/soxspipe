"""Characterization of the `soxs_nod.produce_product` branches the orchestration tests miss.

`tests/integration/test_observing_recipe_orchestration.py` already covers one
balanced AB pair through the stacked path. These tests pin the rest of
`produce_product` before it is split (400 lines): the per-cycle path taken by
several distinct nodding offsets, the flux-standard input that renames the
recipe, the master flat, the raw-frame names recorded per cycle, the
rejection messages, the response curve and flux calibration. They reuse that
module's routing collection, recipe configuration and collaborator stubs, so
every file that drives a nodding recipe builds it the same way.

Where the pinned behaviour looks like a defect, the test name says so and the
behaviour is pinned as it is today rather than corrected here.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.table import Table

from soxspipe.recipes import soxs_nod
from tests.factories import response_product_fits, synthetic_ccd
from tests.integration.test_observing_recipe_orchestration import (
    RouteCollection,
    _configure_nodding_recipe,
    _empty_products,
    _patch_nodding_collaborators,
    _route,
    _write_prepared_frame,
)

pytestmark = pytest.mark.integration

# THE TECHNIQUE EVERY NODDING SCIENCE FRAME IS TAKEN WITH.
NOD_TECHNIQUE = "ECHELLE,SLIT,NODDING"

# THE CALIBRATION ROUTES EVERY REDUCTION IN THIS MODULE NEEDS.
CALIBRATION_NAMES = ("ORDER_TAB_VIS", "DISP_TAB_VIS", "DISP_IMAGE_VIS")

# THE HEADER CARD THE A/B SPLIT READS. POSITIVE IS A, EVERYTHING ELSE IS B.
CUMULATIVE_OFFSET = "ESO SEQ CUMOFF Y"


def _nod_header(offset: float, mjd: float, **extra: object) -> dict[str, object]:
    """Return the header overrides for one nodding frame."""
    return {CUMULATIVE_OFFSET: offset, "MJD-OBS": mjd, **extra}


def _calibration_routes(tmpPath: Path) -> dict[frozenset[tuple[str, object]], list[str]]:
    """Return the order-table and dispersion-map routes."""
    return {_route(PRO_CATG=name): [str(tmpPath / f"{name}.fits")] for name in CALIBRATION_NAMES}


def _nod_recipe(
    log: Any,
    tmpPath: Path,
    headers: list[dict[str, object]],
    *,
    dprType: str = "OBJECT",
    extraRoutes: dict[frozenset[tuple[str, object]], list[str]] | None = None,
    stem: str = "nod",
) -> tuple[soxs_nod, list[Path]]:
    """Return an unconstructed nod recipe routed to one frame per header."""
    paths = [
        _write_prepared_frame(tmpPath / f"{stem}-{index}_pre.fits", seed=200 + index, headerOverrides=header)
        for index, header in enumerate(headers)
    ]
    recipe = soxs_nod.__new__(soxs_nod)
    _configure_nodding_recipe(
        recipe,
        log=log,
        tmp_path=tmpPath,
        objectPaths=paths,
        technique=NOD_TECHNIQUE,
    )
    routes = {
        _route(DPR_TYPE=dprType, DPR_TECH=NOD_TECHNIQUE): [str(path) for path in paths],
        **_calibration_routes(tmpPath),
        **(extraRoutes or {}),
    }
    recipe.inputFrames = RouteCollection(routes)
    return recipe, paths


def _patch(
    recipe: soxs_nod,
    monkeypatch: pytest.MonkeyPatch,
    tmpPath: Path,
    alignFallbackCount: int = 0,
) -> tuple[list[str], dict[str, Any]]:
    """Stub the shared collaborators and return the call log and captured arguments."""
    calls, captured, _ = _patch_nodding_collaborators(
        recipe,
        monkeypatch=monkeypatch,
        plotPath=tmpPath / "OBJECT_VIS_MERGED_QC.pdf",
        extractionPath=tmpPath / "OBJECT_VIS_EXTRACTED.fits",
        alignFallbackCount=alignFallbackCount,
    )
    return calls, captured


def _extract_each_cycle(recipe: soxs_nod) -> soxs_nod:
    """Turn `stack_before_extraction` off, so the recipe takes the per-cycle path."""
    recipe.recipeSettings = {**recipe.recipeSettings, "stack_before_extraction": False}
    return recipe


# ONE BALANCED AB PAIR AT ONE NODDING LOCATION, WHICH TAKES THE STACKED PATH.
SINGLE_PAIR = [_nod_header(3.0, 60000.0), _nod_header(-3.0, 60000.1)]

# TWO AB PAIRS AT TWO DISTINCT A OFFSETS, WHICH TAKES THE PER-CYCLE PATH.
TWO_LOCATIONS = [
    _nod_header(3.0, 60000.0),
    _nod_header(5.0, 60000.2),
    _nod_header(-3.0, 60000.1),
    _nod_header(-5.0, 60000.3),
]


def test_with_the_setting_off_each_cycle_is_extracted_before_one_final_stack(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With `stack_before_extraction` off, each cycle is extracted on its own, with no frame stacking at all."""
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    _extract_each_cycle(recipe)
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    productPath, _ = recipe.produce_product()

    # ASSERT
    assert productPath is None
    assert calls == [
        "keywords",
        "keywords",
        "extract_cycle",
        "keywords",
        "keywords",
        "extract_cycle",
        "stack_extractions",
        "plot",
        "report",
        "clean_up",
    ]
    assert [entry["locationSetIndex"] for entry in captured["extract_cycle"]] == [1, 2]
    assert (
        "print",
        "# PROCESSING 2 AB NODDING CYCLES WITH 2 UNIQUE PAIRS OF OFFSET LOCATIONS",
    ) in log.messages
    assert ("print", "Processing AB Nodding Sequence 1") in log.messages
    assert ("print", "Processing AB Nodding Sequence 2") in log.messages


def test_each_cycle_pairs_the_a_and_b_frames_in_observation_time_order(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both nodding sequences are sorted by `MJD-OBS` before they are paired."""
    # ARRANGE: THE A FRAMES ARRIVE IN REVERSE TIME ORDER, THE B FRAMES IN TIME ORDER.
    headers = [
        _nod_header(5.0, 60000.4),
        _nod_header(3.0, 60000.0),
        _nod_header(-3.0, 60000.1),
        _nod_header(-5.0, 60000.5),
    ]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    _extract_each_cycle(recipe)
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    cycles = captured["extract_cycle"]
    assert [entry["aFrame"].header["MJD-OBS"] for entry in cycles] == [60000.0, 60000.4]
    assert [entry["bFrame"].header["MJD-OBS"] for entry in cycles] == [60000.1, 60000.5]
    assert [entry["aFrame"].header[CUMULATIVE_OFFSET] for entry in cycles] == [3.0, 5.0]


def test_the_frame_name_lists_are_not_sorted_with_their_frames(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A frame whose file order differs from its time order is stamped with another frame's name.

    The frames are sorted by `MJD-OBS`; the parallel name lists are not, so
    the name recorded for a cycle belongs to whichever frame was read in that
    position. Pinned as found, not corrected.
    """
    # ARRANGE: THE FIRST A FRAME READ IS THE LAST IN TIME, AND NO FRAME CARRIES
    # `ARCFILE` OR `ORIGFILE`, SO THE RECORDED NAMES COME FROM THE FILE NAMES.
    headers = [
        _nod_header(5.0, 60000.4),
        _nod_header(3.0, 60000.0),
        _nod_header(-3.0, 60000.1),
        _nod_header(-5.0, 60000.5),
    ]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    _extract_each_cycle(recipe)
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT: CYCLE 1 PROCESSES THE MJD 60000.0 FRAME BUT RECORDS `nod-0.fits`,
    # WHICH IS THE MJD 60000.4 FRAME'S FILE.
    assert [entry["rawFrames"] for entry in captured["keywords"]] == [
        ["nod-0.fits", "nod-2.fits"],
        ["nod-0.fits", "nod-2.fits"],
        ["nod-1.fits", "nod-3.fits"],
        ["nod-1.fits", "nod-3.fits"],
    ]
    assert captured["keywords"][0]["frame"].header["MJD-OBS"] == 60000.0


@pytest.mark.parametrize(
    ("extraHeader", "expectedRawFrames"),
    [
        (
            lambda index: {"ARCFILE": f"ARC-{index}.fits", "ORIGFILE": f"ORIG-{index}.fits"},
            [["ARC-0.fits", "ARC-2.fits"], ["ARC-1.fits", "ARC-3.fits"]],
        ),
        (
            lambda index: {"ORIGFILE": f"ORIG-{index}.fits"},
            [["ORIG-0.fits", "ORIG-2.fits"], ["ORIG-1.fits", "ORIG-3.fits"]],
        ),
        (
            lambda index: {},
            [["nod-0.fits", "nod-2.fits"], ["nod-1.fits", "nod-3.fits"]],
        ),
    ],
)
def test_each_cycle_records_its_raw_frame_names_in_priority_order(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    extraHeader: Any,
    expectedRawFrames: list[list[str]],
) -> None:
    """Per cycle, both frames are stamped with the A and B `ARCFILE`, else `ORIGFILE`, else file name."""
    # ARRANGE
    headers = [{**header, **extraHeader(index)} for index, header in enumerate(TWO_LOCATIONS)]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    _extract_each_cycle(recipe)
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    stamped = [entry["rawFrames"] for entry in captured["keywords"]]
    assert stamped == [
        expectedRawFrames[0],
        expectedRawFrames[0],
        expectedRawFrames[1],
        expectedRawFrames[1],
    ]
    assert all(stamped[index] is stamped[index + 1] for index in (0, 2))


def test_the_per_cycle_spectra_are_concatenated_in_cycle_order_for_the_final_stack(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Each cycle's A and B spectra are appended to their own frame, and the last cycle's joins win."""
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    _extract_each_cycle(recipe)
    calls, captured = _patch(recipe, monkeypatch, tmp_path)
    cycleIndex = 0

    def fake_extract_cycle(**kwargs: object) -> tuple[pd.DataFrame, pd.DataFrame, dict[int, int]]:
        nonlocal cycleIndex
        cycleIndex += 1
        calls.append("extract_cycle")
        captured["extract_cycle"].append(kwargs)
        spectrumA = pd.DataFrame({"WAVE": [500.0 + cycleIndex]})
        spectrumB = pd.DataFrame({"WAVE": [600.0 + cycleIndex]})
        return spectrumA, spectrumB, {"10-11": cycleIndex}

    monkeypatch.setattr(recipe, "process_single_ab_nodding_cycle", fake_extract_cycle)

    # ACT
    recipe.produce_product()

    # ASSERT
    stackArguments = captured["stack_extractions"][0]
    spectrumA, spectrumB = stackArguments["args"][0]
    assert spectrumA["WAVE"].tolist() == [501.0, 502.0]
    assert spectrumB["WAVE"].tolist() == [601.0, 602.0]
    assert stackArguments["orderJoins"] == {"10-11": 2}


def test_one_unique_offset_takes_the_stacked_path_with_a_singular_banner(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two AB pairs at one offset location are stacked, and the banner drops its plural."""
    # ARRANGE
    headers = [
        _nod_header(3.0, 60000.0),
        _nod_header(3.0, 60000.2),
        _nod_header(-3.0, 60000.1),
        _nod_header(-3.0, 60000.3),
    ]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert (
        "print",
        "# PROCESSING 2 AB NODDING CYCLES WITH 1 UNIQUE PAIR OF OFFSET LOCATIONS",
    ) in log.messages
    assert calls.count("stack") == 2
    assert calls.count("extract_cycle") == 1
    assert [len(entry["frames"]) for entry in captured["stack"]] == [2, 2]
    assert [entry["recipe"] for entry in captured["stack"]] == ["soxs_nod", "soxs_nod"]
    assert captured["stack"][0]["ignore_input_masks"] is False
    assert captured["stack"][0]["post_stack_clipping"] is False


def test_a_zero_offset_frame_is_counted_as_a_b_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An offset of exactly zero falls through the positive test into the B sequence."""
    # ARRANGE
    headers = [_nod_header(3.0, 60000.0), _nod_header(0.0, 60000.1)]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["stack"][1]["frames"][0].header[CUMULATIVE_OFFSET] == 0.0
    assert len(captured["stack"][1]["frames"]) == 1


def test_an_unbalanced_nodding_sequence_is_rejected_by_count(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unequal A and B counts raise before any frame is stacked, and log the same sentence."""
    # ARRANGE
    headers = [_nod_header(3.0, 60000.0), _nod_header(5.0, 60000.1), _nod_header(-3.0, 60000.2)]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    calls, _ = _patch(recipe, monkeypatch, tmp_path)
    expected = "Found 2 A frames and 1 B frames. The number of A and B frames must be the same for nodding reductions."

    # ACT / ASSERT
    with pytest.raises(Exception) as raised:
        recipe.produce_product()

    assert str(raised.value) == expected
    assert ("error", expected) in log.messages
    assert calls == []


def test_a_sequence_with_no_positive_offsets_is_rejected_as_unbalanced_first(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The "no positive offset" error cannot be reached: zero A frames fails the balance test first."""
    # ARRANGE
    headers = [_nod_header(-3.0, 60000.0), _nod_header(-5.0, 60000.1)]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    _patch(recipe, monkeypatch, tmp_path)

    # ACT / ASSERT
    with pytest.raises(Exception) as raised:
        recipe.produce_product()

    assert str(raised.value) == (
        "Found 0 A frames and 2 B frames. The number of A and B frames must be the same for nodding reductions."
    )


def test_an_empty_inventory_fails_at_the_first_quicklook(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With no science frames of any type, the recipe raises `IndexError` before the A/B split.

    This is the reachable side of the "no positive offset" branch: a balanced
    sequence with no A frames has no frames at all, and an empty inventory
    fails earlier, at `allObjectFrames[0]`.
    """
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, [])
    calls, _ = _patch(recipe, monkeypatch, tmp_path)

    # ACT / ASSERT
    with pytest.raises(IndexError):
        recipe.produce_product()

    assert calls == []


def test_flux_standard_frames_rename_the_recipe_once_and_skip_later_types(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With no object frames, flux standards are used, the recipe becomes `-std` once, and tellurics are ignored."""
    # ARRANGE
    telluricPaths = [
        _write_prepared_frame(
            tmp_path / f"telluric-{index}.fits",
            seed=300 + index,
            headerOverrides=_nod_header(offset, 60001.0 + index),
        )
        for index, offset in enumerate((9.0, -9.0))
    ]
    recipe, _ = _nod_recipe(
        log,
        tmp_path,
        SINGLE_PAIR,
        dprType="STD,FLUX",
        extraRoutes={_route(DPR_TYPE="STD,TELLURIC", DPR_TECH=NOD_TECHNIQUE): [str(p) for p in telluricPaths]},
    )
    recipe.productDir = str(tmp_path / "products" / "soxs-nod")
    settingsReads: list[str] = []
    stdSettings = {"use_flat": False, "source": "std"}
    monkeypatch.setattr(recipe, "get_recipe_settings", lambda: settingsReads.append(recipe.recipeName) or stdSettings)
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeName == "soxs-nod-std"
    assert settingsReads == ["soxs-nod-std"]
    assert recipe.recipeSettings is stdSettings
    assert recipe.productDir == str(tmp_path / "products" / "soxs-nod-std")
    assert recipe.masterHeaderFrame.header[CUMULATIVE_OFFSET] == 3.0
    assert [entry["frames"][0].header[CUMULATIVE_OFFSET] for entry in captured["stack"]] == [3.0, -3.0]


def test_object_frames_win_over_standard_frames_and_leave_the_recipe_name(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The type loop stops at the first type that matched, so standards are never read."""
    # ARRANGE
    standardPaths = [
        _write_prepared_frame(
            tmp_path / f"standard-{index}.fits",
            seed=310 + index,
            headerOverrides=_nod_header(offset, 60002.0 + index),
        )
        for index, offset in enumerate((7.0, -7.0))
    ]
    recipe, _ = _nod_recipe(
        log,
        tmp_path,
        SINGLE_PAIR,
        extraRoutes={_route(DPR_TYPE="STD,FLUX", DPR_TECH=NOD_TECHNIQUE): [str(p) for p in standardPaths]},
    )
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeName == "soxs-nod"
    assert [entry["frames"][0].header[CUMULATIVE_OFFSET] for entry in captured["stack"]] == [3.0, -3.0]


@pytest.mark.parametrize(
    ("headers", "useFlat", "withFlat", "expectFlat", "expectedCycles", "stackFirst"),
    [
        (SINGLE_PAIR, True, True, True, 1, True),
        (SINGLE_PAIR, True, False, False, 1, True),
        (SINGLE_PAIR, False, True, False, 1, True),
        (TWO_LOCATIONS, True, True, True, 1, True),
        (TWO_LOCATIONS, True, True, True, 2, False),
        (TWO_LOCATIONS, True, False, False, 2, False),
    ],
)
def test_the_master_flat_reaches_extraction_only_when_present_and_enabled(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    headers: list[dict[str, object]],
    useFlat: bool,
    withFlat: bool,
    expectFlat: bool,
    expectedCycles: int,
    stackFirst: bool,
) -> None:
    """Each extraction receives the read master flat only when `use_flat` is set and a flat was supplied."""
    # ARRANGE
    flatRoutes = {}
    flatPath = tmp_path / "MASTER_FLAT_VIS.fits"
    if withFlat:
        _write_prepared_frame(flatPath, seed=400, headerOverrides={})
        flatRoutes = {_route(PRO_CATG="MASTER_FLAT_VIS"): [str(flatPath)]}
    recipe, _ = _nod_recipe(log, tmp_path, headers, extraRoutes=flatRoutes)
    recipe.recipeSettings = {"use_flat": useFlat, "stack_before_extraction": stackFirst}
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    flats = [entry["masterFlat"] for entry in captured["extract_cycle"]]
    assert len(flats) == expectedCycles
    for flat in flats:
        if expectFlat:
            np.testing.assert_array_equal(flat.data, synthetic_ccd(seed=400, prepared=True).data)
            assert str(flat.unit) == "electron"
            assert flat.uncertainty is not None
        else:
            assert flat is False
    if expectFlat:
        assert all(flat is flats[0] for flat in flats)


def test_a_multi_location_standard_is_stacked_so_it_builds_one_response_curve(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Distinct offsets no longer send a standard down the per-cycle path, which cannot build a response."""
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    recipe.generateReponseCurve = True
    calls, captured = _patch(recipe, monkeypatch, tmp_path)
    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", _StubResponse)

    # ACT
    recipe.produce_product()

    # ASSERT: ONE EXTRACTION OF THE FLATTENED STACK AND ONE OF THE UNFLATTENED STACK
    assert calls.count("extract_cycle") == 2
    assert calls.count("stack_extractions") == 2
    assert captured["clean_up"] == [{"forceFail": False}]


def test_a_single_location_standard_builds_a_response_from_an_unflattened_extraction(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The stacked path re-extracts unflattened, stacks it, and hands both extractions to the response."""
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    recipe.generateReponseCurve = True
    calls, captured = _patch(recipe, monkeypatch, tmp_path)
    responseArguments: list[dict[str, object]] = []
    responseQc = pd.DataFrame({"qc_name": ["RESPONSE"]})
    responseProducts = pd.DataFrame({"product_label": ["RESPONSE_TABLE"]})

    class FakeResponse:
        def __init__(self, **kwargs: object) -> None:
            responseArguments.append(kwargs)
            calls.append("response")

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, str]:
            calls.append("response_get")
            return responseQc, responseProducts, "synthetic response failure"

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", FakeResponse)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls == [
        "stack",
        "stack",
        "keywords",
        "keywords",
        "extract_cycle",
        "stack_extractions",
        "extract_cycle",
        "stack_extractions",
        "response",
        "response_get",
        "plot",
        "report",
        "clean_up",
    ]
    firstCycle, secondCycle = captured["extract_cycle"]
    assert "notFlattened" not in firstCycle
    assert secondCycle["notFlattened"] is True
    # THE ALIGNMENT STUB REPORTS NO OBJECT POSITION, SO NEITHER PASS GETS A SHIFT-AND-ADD PLAN (DY-1674)
    assert firstCycle["shiftAndAdd"] is None
    assert secondCycle["shiftAndAdd"] is None
    assert secondCycle["aFrame"] is firstCycle["aFrame"]
    assert secondCycle["bFrame"] is firstCycle["bFrame"]
    assert captured["stack_extractions"][1]["notFlattened"] is True
    arguments = responseArguments[0]
    assert arguments["recipeName"] == "soxs-nod"
    assert arguments["sofName"] == "OBJECT_VIS"
    assert arguments["stdExtractionPath"] == str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")
    assert arguments["stdNotFlatExtractionPath"] == str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")
    assert arguments["orderJoins"] == {"10-11": 11}
    assert arguments["startNightDate"] == "2024-01-02"
    assert recipe.qc is responseQc
    assert ("print", "# CALCULATING RESPONSE FUNCTION\n") in log.messages
    assert captured["clean_up"] == [{"forceFail": "synthetic response failure"}]


def test_a_response_table_flux_calibrates_the_stack_and_draws_a_second_plot(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Flux calibration uses the mean airmass of every A and B frame, adds its products, and draws a second plot."""
    # ARRANGE
    responsePath = response_product_fits(tmp_path / "RESP_TAB_VIS.fits")
    headers = [
        {**SINGLE_PAIR[0], "HIERARCH ESO TEL AIRM START": 1.1, "HIERARCH ESO TEL AIRM END": 1.3, "EXPTIME": 300.0},
        {**SINGLE_PAIR[1], "HIERARCH ESO TEL AIRM START": 2.1, "HIERARCH ESO TEL AIRM END": 2.3, "EXPTIME": 60.0},
    ]
    recipe, _ = _nod_recipe(
        log,
        tmp_path,
        headers,
        extraRoutes={_route(PRO_CATG="RESP_TAB_VIS"): [str(responsePath)]},
    )
    recipe.detectorParams = {"extinction": "extinction.dat"}
    calls, captured = _patch(recipe, monkeypatch, tmp_path)
    fluxPath = tmp_path / "OBJECT_VIS_FLUXCAL.fits"
    Table({"WAVE": [500.0], "FLUX_CALIBRATED": [2.5]}).write(fluxPath)
    fluxProducts = _empty_products()
    fluxProducts.loc[0] = ["soxs-nod", "FLUXCAL", "f.fits", "FITS", "o", "r", "d", "p", "PROD"]
    calibratorArguments: list[dict[str, object]] = []

    class FakeFluxCalibrator:
        def __init__(self, **kwargs: object) -> None:
            calibratorArguments.append(kwargs)
            calls.append("flux_calibrator")

        def calibrate(self) -> tuple[str, pd.DataFrame]:
            calls.append("calibrate")
            return str(fluxPath), fluxProducts

    monkeypatch.setattr(import_module("soxspipe.commonutils.flux_calibration"), "flux_calibration", FakeFluxCalibrator)
    monkeypatch.setattr(import_module("soxspipe.recipes.soxs_nod"), "get_calibrations_path", lambda **_: "/cal")

    # ACT
    productPath, _ = recipe.produce_product()

    # ASSERT: THE NOD RECIPE RETURNS `None` WHATEVER IT CALIBRATED.
    assert productPath is None
    assert calls[-6:] == ["flux_calibrator", "calibrate", "plot", "plot", "report", "clean_up"]
    arguments = calibratorArguments[0]
    assert arguments["responseFunction"] == str(responsePath)
    # BOTH NOD POSITIONS HOLD THE OBJECT: (300 S x 1.2 + 60 S x 2.2) / 360 S (DY-1695)
    assert arguments["airmass"] == pytest.approx(492.0 / 360.0)
    assert arguments["exptime"] == 300.0
    assert arguments["extinctionPath"] == "/cal/extinction.dat"
    assert arguments["arm"] == "VIS"
    assert arguments["recipeName"] == "soxs-nod"
    assert arguments["header"][CUMULATIVE_OFFSET] == 3.0
    assert arguments["extractedSpectrum"] is captured["plot"][0]["merged_orders"]
    assert ("print", "# PERFORMING FLUX CALIBRATION\n") in log.messages
    assert ("print", "# FLUX CALIBRATION COMPLETED\n") in log.messages
    assert "FLUXCAL" in recipe.products["product_label"].tolist()
    firstPlot, secondPlot = captured["plot"]
    assert firstPlot["fluxCalibrated"] is False
    assert secondPlot["fluxCalibrated"] is True
    calibrated = secondPlot["merged_orders"]
    assert calibrated["WAVE"].unit == u.nm
    assert list(calibrated["FLUX_COUNTS"]) == [2.5]
    assert list(calibrated["SNR"]) == [20.0]


def test_the_merged_spectrum_plot_reads_the_template_and_date_stack_extractions_set(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`stack_extractions` writes `filenameTemplate` and `dateObs`, and the plot that follows reads them.

    The collaborator stub every other test in this module uses does not write either attribute, so those
    tests would pass even if the flow between the two blocks broke. This one makes the stub behave like
    the real method and pins that the plot receives what the stack wrote, in that order.
    """
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    calls, captured = _patch(recipe, monkeypatch, tmp_path)
    recipe.filenameTemplate = "STALE.fits"
    recipe.dateObs = "1999-01-01T00:00:00.000"
    stacked = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})

    def stacking_that_stamps_the_recipe(*args: object, **kwargs: object) -> tuple[pd.DataFrame, str]:
        calls.append("stack_extractions")
        captured["stack_extractions"].append({"args": args, **kwargs})
        recipe.filenameTemplate = "2024-01-02_nod.fits"
        recipe.dateObs = "2024-01-02T03:04:05.678"
        return stacked.copy(), str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")

    monkeypatch.setattr(recipe, "stack_extractions", stacking_that_stamps_the_recipe)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls.index("stack_extractions") < calls.index("plot")
    plotArguments = captured["plot"][0]
    assert plotArguments["filenameTemplate"] == "2024-01-02_nod.fits"
    assert plotArguments["dateObs"] == "2024-01-02T03:04:05.678"


@pytest.mark.parametrize("withBias", [True, False], ids=["sof-names-a-bias", "sof-names-no-bias"])
def test_the_master_bias_the_sof_names_is_held_for_the_difference_frames(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    withBias: bool,
) -> None:
    """A VIS nod reduction reads the SOF's master bias for the A-B errors, and warns when there is none (DY-1390)."""
    # ARRANGE
    biasRoutes = {}
    biasPath = tmp_path / "MASTER_BIAS_VIS.fits"
    if withBias:
        _write_prepared_frame(biasPath, seed=500, headerOverrides={})
        biasRoutes = {_route(PRO_CATG="MASTER_BIAS_VIS"): [str(biasPath)]}
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR, extraRoutes=biasRoutes)
    _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    warnings = [message for level, message in log.messages if "MASTER_BIAS_VIS" in message and level == "warning"]
    if withBias:
        np.testing.assert_array_equal(recipe.masterBias.data, synthetic_ccd(seed=500, prepared=True).data)
        assert warnings == []
    else:
        assert recipe.masterBias is False
        assert len(warnings) == 1


def test_a_jittered_set_is_aligned_group_by_group_then_stacked_and_extracted_once(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With `stack_before_extraction` on, distinct offsets no longer choose the per-cycle path (DY-1612)."""
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT: A ALIGNS TO A WITH B AS PARTNER, THEN B ALIGNS TO B WITH A AS PARTNER
    assert [[frame.header["MJD-OBS"] for frame in entry["frames"]] for entry in captured["align"]] == [
        [60000.0, 60000.2],
        [60000.1, 60000.3],
    ]
    assert [[frame.header["MJD-OBS"] for frame in entry["partnerFrames"]] for entry in captured["align"]] == [
        [60000.1, 60000.3],
        [60000.0, 60000.2],
    ]
    assert [entry["frameNames"] for entry in captured["align"]] == [
        ["nod-0.fits", "nod-1.fits"],
        ["nod-2.fits", "nod-3.fits"],
    ]
    assert captured["align"][0]["twoDMapPath"] == str(tmp_path / "DISP_IMAGE_VIS.fits")
    assert captured["align"][0]["dispersionMap"] == str(tmp_path / "DISP_TAB_VIS.fits")
    assert captured["align"][0]["recipeName"] == "soxs-nod"
    assert calls == [
        "stack",
        "stack",
        "keywords",
        "keywords",
        "extract_cycle",
        "stack_extractions",
        "plot",
        "report",
        "clean_up",
    ]
    assert all(frame.header["ALIGNED"] for entry in captured["stack"] for frame in entry["frames"])
    assert [len(entry["frames"]) for entry in captured["stack"]] == [2, 2]


def test_turning_stack_before_extraction_off_extracts_each_cycle_even_without_jitter(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The setting alone chooses the path: two pairs at one offset are extracted cycle by cycle."""
    # ARRANGE
    headers = [
        _nod_header(3.0, 60000.0),
        _nod_header(3.0, 60000.2),
        _nod_header(-3.0, 60000.1),
        _nod_header(-3.0, 60000.3),
    ]
    recipe, _ = _nod_recipe(log, tmp_path, headers)
    recipe.recipeSettings = {"use_flat": False, "stack_before_extraction": False}
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["align"] == []
    assert captured["align_plot"] == []
    assert calls.count("stack") == 0
    assert [entry["locationSetIndex"] for entry in captured["extract_cycle"]] == [1, 2]


def test_a_standard_is_aligned_and_stacked_even_when_the_setting_is_off_and_the_log_says_so(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A standard must still give one response curve, so it stacks whatever the setting says."""
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    recipe.recipeSettings = {"use_flat": False, "stack_before_extraction": False}
    recipe.generateReponseCurve = True
    calls, captured = _patch(recipe, monkeypatch, tmp_path)
    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", _StubResponse)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["align"]) == 2
    assert calls.count("stack") == 2
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("stack_before_extraction" in message and "response" in message for message in warnings)


class _StubResponse:
    """A response-curve builder that reports success and leaves the tables as they are."""

    def __init__(self, **kwargs: Any) -> None:
        self.qcTable = kwargs["qcTable"]
        self.productsTable = kwargs["productsTable"]

    def get(self) -> tuple[Any, Any, bool]:
        return self.qcTable, self.productsTable, False


def test_the_fallback_count_is_a_header_qc_and_the_shifts_get_one_pdf_product(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both groups' fallbacks add up to `SLIT ALIGN FALLBACK NUM`, and one QC plot covers both groups."""
    # ARRANGE: THE FAKE ALIGNMENT REPORTS 2 FALLBACKS FOR EACH GROUP
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    recipe.dateObs = False
    _, captured = _patch(recipe, monkeypatch, tmp_path, alignFallbackCount=2)

    # ACT
    recipe.produce_product()

    # ASSERT
    qcRow = recipe.qc.loc[recipe.qc["qc_name"] == "SLIT ALIGN FALLBACK NUM"].iloc[0]
    assert qcRow["qc_value"] == 4
    assert qcRow["to_header"]
    assert qcRow["obs_date_utc"] == "2024-01-02T03:04:05.678"
    assert len(captured["align_plot"]) == 1
    plotArguments = captured["align_plot"][0]
    assert plotArguments["shifts"]["frame"].tolist() == ["nod-0.fits", "nod-1.fits", "nod-2.fits", "nod-3.fits"]
    assert plotArguments["filePath"] == str(tmp_path / "qc" / "OBJECT_VIS_SLIT_ALIGN_QC_PLOT.pdf")
    plotRow = recipe.products.loc[recipe.products["product_label"] == "SLIT_ALIGNMENT_QC_PLOT"].iloc[0]
    assert plotRow["file_type"] == "PDF"
    assert plotRow["file_path"] == plotArguments["filePath"]


@pytest.mark.parametrize(
    ("instrument", "arm", "perOrder"),
    [
        ("SOXS", "VIS", True),
        ("SOXS", "NIR", False),
        ("XSH", "UVB", False),
        ("XSH", "VIS", False),
        ("XSH", "NIR", False),
    ],
)
def test_only_soxs_vis_gets_a_shift_for_every_order(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    instrument: str,
    arm: str,
    perOrder: bool,
) -> None:
    """SOXS NIR and every X-Shooter arm get one shift per frame (DY-1612 ruling)."""
    # ARRANGE
    armRoutes = {
        _route(PRO_CATG=f"{name}_{arm}"): [str(tmp_path / f"{name}_{arm}.fits")]
        for name in ("ORDER_TAB", "DISP_TAB", "DISP_IMAGE")
    }
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS, extraRoutes=armRoutes)
    recipe.inst = instrument
    recipe.arm = arm
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert [entry["perOrder"] for entry in captured["align"]] == [perOrder, perOrder]
