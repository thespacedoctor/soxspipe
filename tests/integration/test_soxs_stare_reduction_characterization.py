"""Characterization of the `soxs_stare.produce_product` branches the orchestration tests miss.

`tests/integration/test_observing_recipe_orchestration.py` already covers the
object-frame reduction with sky subtraction, the non-object fallbacks, the
response-curve path and flux calibration. These tests pin the rest of
`produce_product` before it is split: the master bias, master dark, NIR dark
and PAE flat-lamp-as-dark reads, the `use_flat` setting and its `KeyError`
fallback, a sky subtraction that returns no model, the unflattening by the
master flat, and both `quicklook_image` calls. They reuse that module's
routing collection and frame writer, so both files build the recipe the same
way.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import StdDevUncertainty

from soxspipe.recipes import soxs_stare
from tests.factories import pipeline_settings, qc_table, synthetic_ccd
from tests.integration.test_observing_recipe_orchestration import (
    RouteCollection,
    _empty_products,
    _keyword,
    _route,
    _write_prepared_frame,
)

pytestmark = pytest.mark.integration

# THE TECHNIQUE EVERY STARE SCIENCE FRAME IS TAKEN WITH.
STARE_TECHNIQUE = "ECHELLE,SLIT,STARE"

# SEEDS OF THE CALIBRATION FRAMES, SO A TEST CAN TELL WHICH FILE WAS READ.
BIAS_SEED = 501
MASTER_DARK_SEED = 502
NIR_DARK_SEED = 503
PAE_DARK_SEED = 504
FLAT_SEED = 505


def _calibration_path(tmpPath: Path, name: str, seed: int) -> Path:
    """Write one prepared calibration frame and return its path."""
    return _write_prepared_frame(tmpPath / f"{name}.fits", seed=seed, headerOverrides={})


def _stare_recipe(
    log: Any,
    tmpPath: Path,
    *,
    arm: str = "VIS",
    extraRoutes: dict[frozenset[tuple[str, object]], list[str]] | None = None,
    recipeSettings: dict[str, object] | None = None,
    settingsOverrides: dict[str, object] | None = None,
) -> soxs_stare:
    """Return an unconstructed stare recipe routed to one object frame and the calibration tables."""
    objectPath = _write_prepared_frame(
        tmpPath / "stare-object.fits",
        seed=500,
        headerOverrides={"DPR_TYPE": "OBJECT", "DPR_TECH": STARE_TECHNIQUE},
    )
    routes = {
        _route(DPR_TYPE="OBJECT", DPR_TECH=STARE_TECHNIQUE): [str(objectPath)],
        _route(PRO_CATG=f"ORDER_TAB_{arm}"): [str(tmpPath / f"ORDER_TAB_{arm}.fits")],
        _route(PRO_CATG=f"DISP_TAB_{arm}"): [str(tmpPath / f"DISP_TAB_{arm}.fits")],
        _route(PRO_CATG=f"DISP_IMAGE_{arm}"): [str(tmpPath / f"DISP_IMAGE_{arm}.fits")],
        **(extraRoutes or {}),
    }
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = arm
    recipe.kw = _keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RouteCollection(routes)
    recipe.recipeName = "soxs-stare"
    recipe.recipeSettings = (
        recipeSettings if recipeSettings is not None else {"use_flat": True, "sky-subtraction": {"subtract_sky": False}}
    )
    recipe.settings = {**pipeline_settings(tmpPath), **(settingsOverrides or {})}
    recipe.generateReponseCurve = False
    recipe.sofName = f"OBJECT_{arm}"
    recipe.startNightDate = "2024-01-02"
    recipe.workspaceRootPath = str(tmpPath)
    recipe.qcDir = str(tmpPath / "qc")
    recipe.filenameTemplate = f"OBJECT_{arm}.fits"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.qc = qc_table()
    recipe.products = _empty_products()
    return recipe


def _patch_stare(
    recipe: soxs_stare,
    monkeypatch: pytest.MonkeyPatch,
    tmpPath: Path,
    *,
    detrended: Any = None,
    skyResult: tuple[Any, Any, Any] | None = None,
) -> dict[str, list[dict[str, object]]]:
    """Stub every collaborator of `produce_product` and return the captured arguments."""
    captured: dict[str, list[dict[str, object]]] = {
        "stack": [],
        "detrend": [],
        "sky": [],
        "write": [],
        "extractor": [],
        "quicklook": [],
    }

    def fake_stack(**kwargs: object) -> Any:
        result = synthetic_ccd(seed=510, prepared=True)
        captured["stack"].append({**kwargs, "result": result})
        return result

    def fake_detrend(**kwargs: object) -> Any:
        result = detrended if detrended is not None else synthetic_ccd(seed=511, prepared=True)
        captured["detrend"].append({**kwargs, "result": result})
        return result

    def fake_write(**kwargs: object) -> str:
        path = str(tmpPath / str(kwargs["filename"]))
        captured["write"].append({**kwargs, "path": path})
        return path

    monkeypatch.setattr(recipe, "clip_and_stack", fake_stack)
    monkeypatch.setattr(recipe, "detrend", fake_detrend)
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **_: None)
    monkeypatch.setattr(recipe, "_write", fake_write)
    monkeypatch.setattr(recipe, "report_output", lambda: recipe.qc)
    monkeypatch.setattr(recipe, "clean_up", lambda **_: None)

    class FakeSkySubtractor:
        def __init__(self, **kwargs: object) -> None:
            captured["sky"].append(kwargs)

        def subtract(self) -> tuple[Any, Any, Any, pd.DataFrame, pd.DataFrame]:
            assert skyResult is not None
            return (*skyResult, recipe.qc, recipe.products)

    class FakeExtractor:
        def __init__(self, **kwargs: object) -> None:
            captured["extractor"].append(kwargs)

        def extract(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, int], str]:
            spectrum = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})
            return recipe.qc, recipe.products, spectrum, {10: 11}, str(tmpPath / "extracted.fits")

    stareModule = import_module("soxspipe.recipes.soxs_stare")
    monkeypatch.setattr(stareModule, "subtract_sky", FakeSkySubtractor)
    monkeypatch.setattr(stareModule, "generic_quality_checks", lambda **kwargs: kwargs["qcTable"])
    monkeypatch.setattr(stareModule, "spectroscopic_image_quality_checks", lambda **kwargs: kwargs["qcTable"])
    commonutils = import_module("soxspipe.commonutils")
    monkeypatch.setattr(commonutils, "horne_extraction", FakeExtractor)
    toolkit = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: captured["quicklook"].append(kwargs))
    monkeypatch.setattr(
        toolkit,
        "plot_merged_spectrum_qc",
        lambda **kwargs: (kwargs["products"], str(tmpPath / "merged.pdf")),
    )
    return captured


def _assert_read_of(frame: Any, seed: int) -> None:
    """Assert that `frame` is the prepared calibration frame written with `seed`, read in electrons."""
    expected = synthetic_ccd(seed=seed, prepared=True)
    np.testing.assert_array_equal(frame.data, expected.data)
    np.testing.assert_array_equal(frame.uncertainty.array, expected.uncertainty.array)
    np.testing.assert_array_equal(frame.mask, expected.mask)
    assert frame.unit == u.electron


def test_master_bias_dark_and_flat_are_read_and_passed_to_detrend(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A UVB/VIS reduction reads the master bias, master dark and master flat and detrends with all three."""
    # ARRANGE
    routes = {
        _route(PRO_CATG="MASTER_BIAS_VIS"): [str(_calibration_path(tmp_path, "MASTER_BIAS_VIS", BIAS_SEED))],
        _route(PRO_CATG="MASTER_DARK_VIS"): [str(_calibration_path(tmp_path, "MASTER_DARK_VIS", MASTER_DARK_SEED))],
        _route(PRO_CATG="MASTER_FLAT_VIS"): [str(_calibration_path(tmp_path, "MASTER_FLAT_VIS", FLAT_SEED))],
    }
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    detrendArgs = captured["detrend"][0]
    _assert_read_of(detrendArgs["master_bias"], BIAS_SEED)
    _assert_read_of(detrendArgs["dark"], MASTER_DARK_SEED)
    _assert_read_of(detrendArgs["master_flat"], FLAT_SEED)
    assert detrendArgs["inputFrame"] is captured["stack"][0]["result"]
    assert detrendArgs["order_table"] == str(tmp_path / "ORDER_TAB_VIS.fits")


def test_a_master_dark_suppresses_the_nir_dark_read(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When a master dark was read, the NIR off-frame route is not consulted and the master dark is used."""
    # ARRANGE
    routes = {
        _route(PRO_CATG="MASTER_DARK_NIR"): [str(_calibration_path(tmp_path, "MASTER_DARK_NIR", MASTER_DARK_SEED))],
        _route(DPR_TYPE="OBJECT", DPR_TECH="IMAGE"): [str(_calibration_path(tmp_path, "NIR_DARK", NIR_DARK_SEED))],
    }
    recipe = _stare_recipe(log, tmp_path, arm="NIR", extraRoutes=routes)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    _assert_read_of(captured["detrend"][0]["dark"], MASTER_DARK_SEED)
    assert captured["detrend"][0]["master_bias"] is False
    consulted = [filters for kind, filters in recipe.inputFrames.calls if kind == "files_filtered"]
    assert {"DPR_TYPE": "OBJECT", "DPR_TECH": "IMAGE"} not in consulted


def test_without_a_master_dark_the_nir_off_frame_is_the_dark(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With no master dark, the OBJECT/IMAGE off-frame is read as the dark."""
    # ARRANGE
    routes = {
        _route(DPR_TYPE="OBJECT", DPR_TECH="IMAGE"): [str(_calibration_path(tmp_path, "NIR_DARK", NIR_DARK_SEED))],
    }
    recipe = _stare_recipe(log, tmp_path, arm="NIR", extraRoutes=routes)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    _assert_read_of(captured["detrend"][0]["dark"], NIR_DARK_SEED)


def test_pae_reads_the_flat_lamp_image_as_the_dark_over_a_master_dark(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """In PAE mode a FLAT,LAMP/IMAGE frame replaces whatever dark was read before it."""
    # ARRANGE
    routes = {
        _route(PRO_CATG="MASTER_DARK_VIS"): [str(_calibration_path(tmp_path, "MASTER_DARK_VIS", MASTER_DARK_SEED))],
        _route(DPR_TYPE="FLAT,LAMP", DPR_TECH="IMAGE"): [str(_calibration_path(tmp_path, "PAE_DARK", PAE_DARK_SEED))],
    }
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes, settingsOverrides={"PAE": True})
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    _assert_read_of(captured["detrend"][0]["dark"], PAE_DARK_SEED)
    assert recipe.subtractSky is False


@pytest.mark.parametrize(
    ("recipeSettings", "expectFlat"),
    [
        ({"use_flat": True, "sky-subtraction": {"subtract_sky": False}}, True),
        ({"use_flat": False, "sky-subtraction": {"subtract_sky": False}}, False),
        ({"sky-subtraction": {"subtract_sky": False}}, False),
    ],
    ids=["use_flat-true", "use_flat-false", "use_flat-missing"],
)
def test_the_master_flat_reaches_detrend_only_when_use_flat_is_set(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    recipeSettings: dict[str, object],
    expectFlat: bool,
) -> None:
    """`use_flat` false, or missing from the recipe settings, drops a master flat that was read."""
    # ARRANGE
    routes = {
        _route(PRO_CATG="MASTER_FLAT_VIS"): [str(_calibration_path(tmp_path, "MASTER_FLAT_VIS", FLAT_SEED))],
    }
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes, recipeSettings=recipeSettings)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    masterFlat = captured["detrend"][0]["master_flat"]
    if expectFlat:
        _assert_read_of(masterFlat, FLAT_SEED)
    else:
        assert masterFlat is False


def test_a_sky_subtraction_with_no_model_falls_back_to_the_detrended_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A `None` sky model turns sky subtraction off, widens the Horne slit and extracts the detrended frame."""
    # ARRANGE
    recipeSettings = {"use_flat": False, "sky-subtraction": {"subtract_sky": True}}
    recipe = _stare_recipe(log, tmp_path, recipeSettings=recipeSettings)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=(None, None, None))

    # ACT
    productPath, _ = recipe.produce_product()

    # ASSERT
    assert productPath is None
    assert captured["write"] == []
    assert recipe.subtractSky is False
    assert recipe.recipeSettings["sky-subtraction"]["subtract_sky"] is False
    assert recipe.recipeSettings["horne-extraction-slit-length"] == 30.0
    detrended = captured["detrend"][0]["result"]
    assert captured["sky"][0]["objectFrame"] is detrended
    extractorArgs = captured["extractor"][0]
    assert extractorArgs["skySubtractedFrame"] is detrended
    assert extractorArgs["subtractedFrame"] is False
    # WITH SKY SUBTRACTION NOW OFF, THE UNFLATTENED FRAME IS THE RAW STACK
    assert extractorArgs["unflattenedFrame"] is captured["stack"][0]["result"]


def test_sky_subtracted_frame_is_unflattened_by_the_master_flat(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With sky subtraction and a master flat, the unflattened frame is the product with the flat."""
    # ARRANGE
    routes = {
        _route(PRO_CATG="MASTER_FLAT_VIS"): [str(_calibration_path(tmp_path, "MASTER_FLAT_VIS", FLAT_SEED))],
    }
    recipeSettings = {"use_flat": True, "sky-subtraction": {"subtract_sky": True}}
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes, recipeSettings=recipeSettings)
    detrended = synthetic_ccd(seed=520, prepared=True)
    detrended.uncertainty = StdDevUncertainty(np.full(detrended.shape, 2.5), unit=u.electron)
    skyModel = synthetic_ccd(seed=521, prepared=True)
    skySubtracted = synthetic_ccd(seed=522, prepared=True)
    residuals = synthetic_ccd(seed=523, prepared=True)
    captured = _patch_stare(
        recipe,
        monkeypatch,
        tmp_path,
        detrended=detrended,
        skyResult=(skyModel, skySubtracted, residuals),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    extractorArgs = captured["extractor"][0]
    unflattened = extractorArgs["unflattenedFrame"]
    flatData = synthetic_ccd(seed=FLAT_SEED, prepared=True).data
    np.testing.assert_allclose(unflattened.data, skySubtracted.data * flatData, rtol=1e-12, atol=0)
    np.testing.assert_array_equal(unflattened.uncertainty.array, np.full(detrended.shape, 2.5))
    assert unflattened.header is skySubtracted.header
    assert extractorArgs["skySubtractedFrame"] is skySubtracted
    assert extractorArgs["subtractedFrame"] is skyModel


def test_sky_subtracted_frame_without_a_flat_is_its_own_unflattened_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With sky subtraction and no master flat, the sky-subtracted frame is passed as the unflattened frame."""
    # ARRANGE
    recipeSettings = {"use_flat": True, "sky-subtraction": {"subtract_sky": True}}
    recipe = _stare_recipe(log, tmp_path, recipeSettings=recipeSettings)
    skySubtracted = synthetic_ccd(seed=531, prepared=True)
    captured = _patch_stare(
        recipe,
        monkeypatch,
        tmp_path,
        skyResult=(synthetic_ccd(seed=530, prepared=True), skySubtracted, synthetic_ccd(seed=532, prepared=True)),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["extractor"][0]["unflattenedFrame"] is skySubtracted


def test_both_quicklook_images_are_rendered_with_their_arguments(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`produce_product` renders the detrended frame, then the leftover debug surface plot of the raw stack (DY-112)."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    first, second = captured["quicklook"]
    assert first == {
        "log": recipe.log,
        "CCDObject": captured["detrend"][0]["result"],
        "show": False,
        "ext": False,
        "stdWindow": 3,
        "title": False,
        "surfacePlot": False,
        "dispMap": str(tmp_path / "DISP_TAB_VIS.fits"),
        "dispMapImage": str(tmp_path / "DISP_IMAGE_VIS.fits"),
        "settings": recipe.settings,
        "skylines": False,
    }
    assert second == {
        "log": recipe.log,
        "CCDObject": captured["stack"][0]["result"],
        "show": True,
        "ext": "data",
        "stdWindow": 3,
        "title": False,
        "surfacePlot": True,
    }
