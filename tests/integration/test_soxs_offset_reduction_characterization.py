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


def _patch(recipe: soxs_offset, monkeypatch: pytest.MonkeyPatch, tmpPath: Path) -> tuple[list[str], dict[str, Any]]:
    """Stub the shared collaborators and return the call log and captured arguments."""
    calls, captured, _ = _patch_nodding_collaborators(
        recipe,
        monkeypatch=monkeypatch,
        plotPath=tmpPath / "OBJECT_VIS_MERGED_QC.pdf",
        extractionPath=tmpPath / "OBJECT_VIS_EXTRACTED.fits",
    )
    return calls, captured


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
    ("headers", "useFlat", "withFlat", "expectFlat"),
    [
        (SINGLE_PAIR, True, True, True),
        (SINGLE_PAIR, True, False, False),
        (SINGLE_PAIR, False, True, False),
        (TWO_LOCATIONS, True, True, True),
        (TWO_LOCATIONS, True, False, False),
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
) -> None:
    """Each extraction receives the read master flat only when `use_flat` is set and a flat was supplied."""
    # ARRANGE
    flatRoutes = {}
    flatPath = tmp_path / "MASTER_FLAT_VIS.fits"
    if withFlat:
        _write_prepared_frame(flatPath, seed=400, headerOverrides={})
        flatRoutes = {_route(PRO_CATG="MASTER_FLAT_VIS"): [str(flatPath)]}
    recipe, _ = _offset_recipe(log, tmp_path, headers, extraRoutes=flatRoutes)
    recipe.recipeSettings = {"use_flat": useFlat}
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    flats = [entry["masterFlat"] for entry in captured["extract_cycle"]]
    assert len(flats) == (1 if headers is SINGLE_PAIR else 2)
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


def test_a_multi_location_standard_forces_failure_instead_of_building_a_response(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A response curve is not built on the per-cycle path; the run is marked failed instead."""
    # ARRANGE
    recipe, _ = _offset_recipe(log, tmp_path, TWO_LOCATIONS)
    recipe.generateReponseCurve = True
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    productPath, _ = recipe.produce_product()

    # ASSERT
    assert productPath == str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")
    assert calls.count("extract_cycle") == 2
    assert calls.count("stack_extractions") == 1
    assert captured["clean_up"] == [{"forceFail": True}]


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
    assert arguments["orderJoins"] == {10: 11}
    assert arguments["recipeName"] == "soxs-offset"
    assert captured["clean_up"] == [{"forceFail": "synthetic response failure"}]


def test_a_response_table_flux_calibrates_the_stack_and_returns_the_calibrated_path(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Flux calibration reads the first ON frame's header, adds its products, and draws a second QC plot."""
    # ARRANGE
    responsePath = tmp_path / "RESP_TAB_VIS.fits"
    headers = [
        {**SINGLE_PAIR[0], "HIERARCH ESO TEL AIRM END": 1.3, "EXPTIME": 300.0},
        {**SINGLE_PAIR[1], "HIERARCH ESO TEL AIRM END": 1.9, "EXPTIME": 60.0},
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
    assert arguments["airmass"] == 1.3
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
