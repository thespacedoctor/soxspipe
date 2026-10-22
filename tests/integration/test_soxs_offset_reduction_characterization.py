"""Characterization of the `soxs_offset.produce_product` branches the orchestration tests miss.

`tests/integration/test_observing_recipe_orchestration.py` already covers the
balanced single-location reduction, the multi-location (jitter) reduction and
the unbalanced inventory. These tests pin the rest of `produce_product`
before it is split: the flux-standard input that renames the recipe, the
master flat, the raw-frame names recorded per cycle, both response-curve
paths, flux calibration, and the empty inventory. They reuse that module's
routing collection, recipe configuration and collaborator stubs, so both files
build the recipe the same way.
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

from soxspipe.recipes import soxs_offset
from soxspipe.recipes.soxs_offset import STACKED_LOCATION_SET
from tests.factories import synthetic_ccd
from tests.integration.test_observing_recipe_orchestration import (
    RouteCollection,
    _configure_nodding_recipe,
    _empty_products,
    _patch_nodding_collaborators,
    _route,
    _write_prepared_frame,
)

pytestmark = pytest.mark.integration

# THE TECHNIQUE EVERY OFFSET SCIENCE FRAME IS TAKEN WITH.
OFFSET_TECHNIQUE = "ECHELLE,SLIT,OFFSET"

# THE CALIBRATION ROUTES EVERY REDUCTION IN THIS MODULE NEEDS.
CALIBRATION_NAMES = ("ORDER_TAB_VIS", "DISP_TAB_VIS", "DISP_IMAGE_VIS")


def _offset_header(ra: float, dec: float, mjd: float, **extra: object) -> dict[str, object]:
    """Return the header overrides for one offset frame."""
    return {
        "HIERARCH ESO SEQ FIXOFF RA": ra,
        "HIERARCH ESO SEQ FIXOFF DEC": dec,
        "MJD-OBS": mjd,
        **extra,
    }


def _calibration_routes(tmpPath: Path) -> dict[frozenset[tuple[str, object]], list[str]]:
    """Return the order-table and dispersion-map routes."""
    return {_route(PRO_CATG=name): [str(tmpPath / f"{name}.fits")] for name in CALIBRATION_NAMES}


def _offset_recipe(
    log: Any,
    tmpPath: Path,
    headers: list[dict[str, object]],
    *,
    dprType: str = "OBJECT",
    extraRoutes: dict[frozenset[tuple[str, object]], list[str]] | None = None,
    stem: str = "offset",
) -> tuple[soxs_offset, list[Path]]:
    """Return an unconstructed offset recipe routed to one frame per header."""
    paths = [
        _write_prepared_frame(tmpPath / f"{stem}-{index}_pre.fits", seed=200 + index, headerOverrides=header)
        for index, header in enumerate(headers)
    ]
    recipe = soxs_offset.__new__(soxs_offset)
    _configure_nodding_recipe(
        recipe,
        log=log,
        tmp_path=tmpPath,
        objectPaths=paths,
        technique=OFFSET_TECHNIQUE,
    )
    routes = {
        _route(DPR_TYPE=dprType, DPR_TECH=OFFSET_TECHNIQUE): [str(path) for path in paths],
        **_calibration_routes(tmpPath),
        **(extraRoutes or {}),
    }
    recipe.inputFrames = RouteCollection(routes)
    return recipe, paths


def _patch(
    recipe: soxs_offset,
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


def _extract_each_cycle(recipe: soxs_offset) -> soxs_offset:
    """Turn `stack_before_extraction` off, so the recipe takes the per-cycle path."""
    recipe.recipeSettings = {**recipe.recipeSettings, "stack_before_extraction": False}
    return recipe


# A BALANCED SINGLE-LOCATION PAIR: ONE ON FRAME AND ONE OFF FRAME.
SINGLE_PAIR = [_offset_header(-2.0, 0.0, 60000.0), _offset_header(0.0, 0.0, 60000.1)]

# TWO ON/OFF PAIRS AT TWO DISTINCT OFFSET LOCATIONS, WHICH TAKES THE PER-CYCLE PATH.
TWO_LOCATIONS = [
    _offset_header(-3.0, 1.0, 60000.0),
    _offset_header(-4.0, 2.0, 60000.1),
    _offset_header(0.0, 1.0, 60000.2),
    _offset_header(0.0, 2.0, 60000.3),
]


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
            headerOverrides=_offset_header(ra, 0.0, 60001.0 + index),
        )
        for index, ra in enumerate((-9.0, 0.0))
    ]
    recipe, standardPaths = _offset_recipe(
        log,
        tmp_path,
        SINGLE_PAIR,
        dprType="STD,FLUX",
        extraRoutes={_route(DPR_TYPE="STD,TELLURIC", DPR_TECH=OFFSET_TECHNIQUE): [str(p) for p in telluricPaths]},
    )
    recipe.productDir = str(tmp_path / "products" / "soxs-offset")
    settingsReads: list[str] = []
    stdSettings = {"use_flat": False, "source": "std"}
    monkeypatch.setattr(recipe, "get_recipe_settings", lambda: settingsReads.append(recipe.recipeName) or stdSettings)
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeName == "soxs-offset-std"
    assert settingsReads == ["soxs-offset-std"]
    assert recipe.recipeSettings is stdSettings
    assert recipe.productDir == str(tmp_path / "products" / "soxs-offset-std")
    assert recipe.masterHeaderFrame.header["HIERARCH ESO SEQ FIXOFF RA"] == -2.0
    stackedRa = [entry["frames"][0].header["HIERARCH ESO SEQ FIXOFF RA"] for entry in captured["stack"]]
    assert stackedRa == [-2.0, 0.0]
    assert [len(entry["frames"]) for entry in captured["stack"]] == [1, 1]
    assert len(standardPaths) == 2


@pytest.mark.parametrize(
    ("headers", "useFlat", "withFlat", "expectFlat", "stackFirst"),
    [
        (SINGLE_PAIR, True, True, True, True),
        (SINGLE_PAIR, True, False, False, True),
        (SINGLE_PAIR, False, True, False, True),
        (TWO_LOCATIONS, True, True, True, True),
        (TWO_LOCATIONS, True, True, True, False),
        (TWO_LOCATIONS, True, False, False, False),
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
    stackFirst: bool,
) -> None:
    """Each extraction receives the read master flat only when `use_flat` is set and a flat was supplied."""
    # ARRANGE
    flatRoutes = {}
    flatPath = tmp_path / "MASTER_FLAT_VIS.fits"
    if withFlat:
        _write_prepared_frame(flatPath, seed=400, headerOverrides={})
        flatRoutes = {_route(PRO_CATG="MASTER_FLAT_VIS"): [str(flatPath)]}
    recipe, _ = _offset_recipe(log, tmp_path, headers, extraRoutes=flatRoutes)
    recipe.recipeSettings = {"use_flat": useFlat, "stack_before_extraction": stackFirst}
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    flats = [entry["masterFlat"] for entry in captured["extract_cycle"]]
    assert len(flats) == (2 if headers is TWO_LOCATIONS and not stackFirst else 1)
    for flat in flats:
        if expectFlat:
            np.testing.assert_array_equal(flat.data, synthetic_ccd(seed=400, prepared=True).data)
            assert str(flat.unit) == "electron"
            assert flat.uncertainty is not None
        else:
            assert flat is False
    if expectFlat:
        assert all(flat is flats[0] for flat in flats)


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
            [["offset-0.fits", "offset-2.fits"], ["offset-1.fits", "offset-3.fits"]],
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
    """Per cycle, both frames are stamped with the ON and OFF `ARCFILE`, else `ORIGFILE`, else file name."""
    # ARRANGE
    headers = [{**header, **extraHeader(index)} for index, header in enumerate(TWO_LOCATIONS)]
    recipe, _ = _offset_recipe(log, tmp_path, headers)
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


def test_a_single_location_standard_builds_a_response_from_an_unflattened_extraction(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The stacked path re-extracts unflattened, stacks it, and hands both extractions to the response."""
    # ARRANGE
    recipe, _ = _offset_recipe(log, tmp_path, SINGLE_PAIR)
    recipe.generateReponseCurve = True
    calls, captured = _patch(recipe, monkeypatch, tmp_path)
    responseArguments: list[dict[str, object]] = []

    class FakeResponse:
        def __init__(self, **kwargs: object) -> None:
            responseArguments.append(kwargs)
            calls.append("response")

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, str]:
            calls.append("response_get")
            return recipe.qc, recipe.products, "synthetic response failure"

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", FakeResponse)

    # ACT
    productPath, _ = recipe.produce_product()

    # ASSERT
    assert productPath == str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")
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
    assert secondCycle["aFrame"] is firstCycle["aFrame"]
    assert secondCycle["bFrame"] is firstCycle["bFrame"]
    assert captured["stack_extractions"][1]["notFlattened"] is True
    arguments = responseArguments[0]
    assert arguments["stdExtractionPath"] == str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")
    assert arguments["stdNotFlatExtractionPath"] == str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")
    assert arguments["orderJoins"] == {"10-11": 11}
    assert arguments["recipeName"] == "soxs-offset"
    assert captured["clean_up"] == [{"forceFail": "synthetic response failure"}]


def test_the_response_of_a_standard_uses_the_mean_airmass_of_the_on_frames_only(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The OFF frames are sky only, so their airmass does not enter the response curve (DY-1713)."""
    # ARRANGE
    headers = [
        {**SINGLE_PAIR[0], "HIERARCH ESO TEL AIRM START": 1.1, "HIERARCH ESO TEL AIRM END": 1.3, "EXPTIME": 300.0},
        {**SINGLE_PAIR[1], "HIERARCH ESO TEL AIRM START": 1.7, "HIERARCH ESO TEL AIRM END": 1.9, "EXPTIME": 60.0},
    ]
    recipe, _ = _offset_recipe(log, tmp_path, headers)
    recipe.generateReponseCurve = True
    _patch(recipe, monkeypatch, tmp_path)
    responseArguments: list[dict[str, object]] = []

    class FakeResponse:
        def __init__(self, **kwargs: object) -> None:
            responseArguments.append(kwargs)

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
            return recipe.qc, recipe.products, False

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", FakeResponse)

    # ACT
    recipe.produce_product()

    # ASSERT: THE ON FRAME'S START AND END MIDPOINT; WITH THE OFF FRAME IT WOULD BE (300 x 1.2 + 60 x 1.8) / 360
    assert responseArguments[0]["airmass"] == pytest.approx(1.2)


def test_a_response_table_flux_calibrates_the_stack_and_returns_the_calibrated_path(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Flux calibration uses the mean airmass of the ON frames only, adds its products, and draws a second plot."""
    # ARRANGE
    responsePath = tmp_path / "RESP_TAB_VIS.fits"
    headers = [
        {**SINGLE_PAIR[0], "HIERARCH ESO TEL AIRM START": 1.1, "HIERARCH ESO TEL AIRM END": 1.3, "EXPTIME": 300.0},
        {**SINGLE_PAIR[1], "HIERARCH ESO TEL AIRM START": 1.7, "HIERARCH ESO TEL AIRM END": 1.9, "EXPTIME": 60.0},
    ]
    recipe, _ = _offset_recipe(
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
    fluxProducts.loc[0] = ["soxs-offset", "FLUXCAL", "f.fits", "FITS", "o", "r", "d", "p", "PROD"]
    calibratorArguments: list[dict[str, object]] = []

    class FakeFluxCalibrator:
        def __init__(self, **kwargs: object) -> None:
            calibratorArguments.append(kwargs)
            calls.append("flux_calibrator")

        def calibrate(self) -> tuple[str, pd.DataFrame]:
            calls.append("calibrate")
            return str(fluxPath), fluxProducts

    monkeypatch.setattr(import_module("soxspipe.commonutils.flux_calibration"), "flux_calibration", FakeFluxCalibrator)
    monkeypatch.setattr(import_module("soxspipe.recipes.soxs_offset"), "get_calibrations_path", lambda **_: "/cal")

    # ACT
    productPath, _ = recipe.produce_product()

    # ASSERT
    assert productPath == str(fluxPath)
    assert calls[-7:] == ["stack_extractions", "flux_calibrator", "calibrate", "plot", "plot", "report", "clean_up"]
    arguments = calibratorArguments[0]
    assert arguments["responseFunction"] == str(responsePath)
    # THE OFF FRAME IS SKY ONLY, SO ONLY THE ON FRAME'S START AND END MIDPOINT COUNTS (DY-1695)
    assert arguments["airmass"] == pytest.approx(1.2)
    assert arguments["exptime"] == 300.0
    assert arguments["extinctionPath"] == "/cal/extinction.dat"
    assert arguments["header"]["HIERARCH ESO SEQ FIXOFF RA"] == -2.0
    assert arguments["extractedSpectrum"] is captured["plot"][0]["merged_orders"]
    assert "FLUXCAL" in recipe.products["product_label"].tolist()
    firstPlot, secondPlot = captured["plot"]
    assert firstPlot["fluxCalibrated"] is False
    assert secondPlot["fluxCalibrated"] is True
    calibrated = secondPlot["merged_orders"]
    assert calibrated["WAVE"].unit == u.nm
    assert list(calibrated["FLUX_COUNTS"]) == [2.5]
    assert list(calibrated["SNR"]) == [20.0]


def test_an_empty_inventory_fails_at_the_first_quicklook(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With no science frames of any type, the recipe raises `IndexError` before the ON/OFF split."""
    # ARRANGE
    recipe, _ = _offset_recipe(log, tmp_path, [])
    calls, _ = _patch(recipe, monkeypatch, tmp_path)

    # ACT / ASSERT
    with pytest.raises(IndexError):
        recipe.produce_product()

    assert calls == []


def test_each_cycle_and_the_stacked_pair_extract_under_their_own_location_set(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Cycles are numbered from one, and the stacked pair takes the stack's own index so its products never collide."""
    # ARRANGE
    cyclesRecipe, _ = _offset_recipe(log, tmp_path, TWO_LOCATIONS, stem="cycles")
    _extract_each_cycle(cyclesRecipe)
    _, cyclesCaptured = _patch(cyclesRecipe, monkeypatch, tmp_path)
    stackRecipe, _ = _offset_recipe(log, tmp_path, SINGLE_PAIR, stem="stack")
    stackRecipe.generateReponseCurve = False
    _, stackCaptured = _patch(stackRecipe, monkeypatch, tmp_path)

    # ACT
    cyclesRecipe.produce_product()
    stackRecipe.produce_product()

    # ASSERT
    assert [call["locationSetIndex"] for call in cyclesCaptured["extract_cycle"]] == [1, 2]
    assert [call["locationSetIndex"] for call in stackCaptured["extract_cycle"]] == [STACKED_LOCATION_SET]


@pytest.mark.parametrize("withBias", [True, False], ids=["sof-names-a-bias", "sof-names-no-bias"])
def test_the_master_bias_the_sof_names_is_held_for_the_difference_frames(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    withBias: bool,
) -> None:
    """A VIS offset reduction reads the SOF's master bias for the A-B errors, and warns when there is none (DY-1390)."""
    # ARRANGE
    biasRoutes = {}
    biasPath = tmp_path / "MASTER_BIAS_VIS.fits"
    if withBias:
        _write_prepared_frame(biasPath, seed=500, headerOverrides={})
        biasRoutes = {_route(PRO_CATG="MASTER_BIAS_VIS"): [str(biasPath)]}
    recipe, _ = _offset_recipe(log, tmp_path, SINGLE_PAIR, extraRoutes=biasRoutes)
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


def test_a_jittered_set_aligns_the_on_frames_and_stacks_the_off_frames_unaligned(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With `stack_before_extraction` on, only the ON frames are aligned, with the OFF frames as their partners."""
    # ARRANGE
    recipe, _ = _offset_recipe(log, tmp_path, TWO_LOCATIONS)
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["align"]) == 1
    alignment = captured["align"][0]
    assert [frame.header["MJD-OBS"] for frame in alignment["frames"]] == [60000.0, 60000.1]
    assert [frame.header["MJD-OBS"] for frame in alignment["partnerFrames"]] == [60000.2, 60000.3]
    assert alignment["frameNames"] == ["offset-0.fits", "offset-1.fits"]
    assert alignment["recipeName"] == "soxs-offset"
    assert alignment["perOrder"] is True
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
    onStack, offStack = captured["stack"]
    assert all(frame.header["ALIGNED"] for frame in onStack["frames"])
    assert [len(onStack["frames"]), len(offStack["frames"])] == [2, 2]
    assert not any("ALIGNED" in frame.header for frame in offStack["frames"])
    assert [entry["locationSetIndex"] for entry in captured["extract_cycle"]] == [STACKED_LOCATION_SET]


def test_turning_stack_before_extraction_off_extracts_each_offset_cycle_even_without_jitter(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The setting alone chooses the path: two ON/OFF pairs at one offset are extracted cycle by cycle."""
    # ARRANGE
    headers = [
        _offset_header(-3.0, 1.0, 60000.0),
        _offset_header(-3.0, 1.0, 60000.1),
        _offset_header(0.0, 1.0, 60000.2),
        _offset_header(0.0, 1.0, 60000.3),
    ]
    recipe, _ = _offset_recipe(log, tmp_path, headers)
    _extract_each_cycle(recipe)
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["align"] == []
    assert calls.count("stack") == 0
    assert [entry["locationSetIndex"] for entry in captured["extract_cycle"]] == [1, 2]


def test_a_multi_location_standard_is_aligned_stacked_and_builds_a_response_even_with_the_setting_off(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A standard always stacks, so the per-cycle path's forced failure no longer ends a jittered standard."""
    # ARRANGE
    recipe, _ = _offset_recipe(log, tmp_path, TWO_LOCATIONS)
    _extract_each_cycle(recipe)
    recipe.generateReponseCurve = True
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    class FakeResponse:
        def __init__(self, **kwargs: object) -> None:
            calls.append("response")

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
            return recipe.qc, recipe.products, False

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", FakeResponse)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["align"]) == 1
    assert calls.count("extract_cycle") == 2
    assert "response" in calls
    assert captured["clean_up"] == [{"forceFail": False}]
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("stack_before_extraction" in message and "response" in message for message in warnings)


def test_the_offset_fallback_count_is_a_header_qc_and_the_on_shifts_get_one_pdf_product(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The ON group's fallbacks become `SLIT ALIGN FALLBACK NUM`, and its shifts one QC plot product."""
    # ARRANGE
    recipe, _ = _offset_recipe(log, tmp_path, TWO_LOCATIONS)
    recipe.dateObs = False
    _, captured = _patch(recipe, monkeypatch, tmp_path, alignFallbackCount=3)

    # ACT
    recipe.produce_product()

    # ASSERT
    qcRow = recipe.qc.loc[recipe.qc["qc_name"] == "SLIT ALIGN FALLBACK NUM"].iloc[0]
    assert qcRow["qc_value"] == 3
    assert qcRow["to_header"]
    assert qcRow["soxspipe_recipe"] == "soxs-offset"
    assert qcRow["obs_date_utc"] == "2024-01-02T03:04:05.678"
    assert captured["align_plot"][0]["shifts"]["frame"].tolist() == ["offset-0.fits", "offset-1.fits"]
    plotRow = recipe.products.loc[recipe.products["product_label"] == "SLIT_ALIGNMENT_QC_PLOT"].iloc[0]
    assert plotRow["file_path"] == str(tmp_path / "qc" / "OBJECT_VIS_SLIT_ALIGN_QC_PLOT.pdf")
