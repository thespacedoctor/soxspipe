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

from soxspipe.commonutils.slit_illumination import (
    FLAT_EDGE_COLUMNS,
    OFFSET_COLUMNS,
    PROFILE_COLUMNS,
    STATUS_FAILED,
    STATUS_FITTED,
    _read_disp_map,
    slit_illumination_result,
)
from soxspipe.commonutils.slit_illumination_mask import overlap_exclude_mask
from soxspipe.recipes import soxs_stare
from tests.factories import pipeline_settings, qc_table, response_product_fits, synthetic_ccd
from tests.factories.slit_illumination import write_map_images
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
    traceFallback: bool = False,
) -> dict[str, list[dict[str, object]]]:
    """Stub every collaborator of `produce_product` and return the captured arguments."""
    captured: dict[str, list[dict[str, object]]] = {
        "stack": [],
        "detrend": [],
        "sky": [],
        "write": [],
        "extractor": [],
        "quicklook": [],
        "cleanUp": [],
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
    monkeypatch.setattr(recipe, "clean_up", lambda **kwargs: captured["cleanUp"].append(kwargs))

    class FakeSkySubtractor:
        def __init__(self, **kwargs: object) -> None:
            captured["sky"].append(kwargs)

        def subtract(self) -> tuple[Any, Any, Any, pd.DataFrame, pd.DataFrame]:
            assert skyResult is not None
            return (*skyResult, recipe.qc, recipe.products)

    class FakeExtractor:
        # NO SLIT-EDGE SPILL, MATCHING horne_extraction.slitEdgeOrders
        slitEdgeOrders: list = []
        extractionNotes: list = []
        # NO SKY SHIFTS TO HAND ON, MATCHING horne_extraction.flatPassHandoff BEFORE extract() TUNES (DY-1680)
        flatPassHandoff = None

        def __init__(self, **kwargs: object) -> None:
            captured["extractor"].append(kwargs)
            self.traceFallback = traceFallback

        def extract(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, int], str]:
            spectrum = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})
            return recipe.qc, recipe.products, spectrum, {"10-11": 11}, str(tmpPath / "extracted.fits")

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
    skySubtracted.uncertainty = StdDevUncertainty(np.full(skySubtracted.shape, 1.5), unit=u.electron)
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
    # THE SKY-SUBTRACTED ERROR IS SCALED BACK BY THE FLAT, NOT TAKEN FROM THE DETRENDED FRAME (DY-1391)
    np.testing.assert_allclose(unflattened.uncertainty.array, 1.5 * np.abs(flatData), rtol=1e-6, atol=0)
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
    """`produce_product` renders the detrended frame, then the debug surface plot of the raw stack (DY-112)."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path)
    recipe.debug = False
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
        "show": False,
        "ext": "data",
        "stdWindow": 3,
        "title": False,
        "surfacePlot": True,
    }


def test_raw_stack_surface_plot_is_shown_only_in_debug_mode(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The raw-stack surface plot blocks on a window, so `produce_product` shows it only when debugging."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path)
    recipe.debug = True
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["quicklook"][1]["show"] is True


# THE SLIT-ILLUMINATION BLOCK THE SOXS DEFAULT SETTINGS CARRY WITH THE MEASUREMENT SWITCHED ON
ILLUMINATION_BLOCK = {
    "measure-offset": True,
    "overlap-margin-arcsec": 0.45,
    "profile-bin-arcsec": 0.1,
    "max-offset-arcsec": 3.0,
    "max-order-scatter-arcsec": 0.3,
}
ILLUMINATION_QC_NAMES = ("SLIT ILLUM OFFSET O11", "SLIT ILLUM OFFSET O12", "SLIT ILLUM APPLIED", "SLIT ILLUM N ORDERS")
MAP_SHAPE = (32, 32)
OFFSET_ARCSEC = {11: 0.8, 12: -0.4}


def _illumination_recipe_settings(block: dict[str, object] | None, *, useFlat: bool = True) -> dict[str, object]:
    recipeSettings: dict[str, object] = {"use_flat": useFlat, "sky-subtraction": {"subtract_sky": True}}
    if block is not None:
        recipeSettings["slit-illumination"] = block
    return recipeSettings


def _write_illumination_map(tmpPath: Path, arm: str = "VIS") -> Path:
    """Write a 2D map with order 11 in the top half and order 12 in the bottom half, SLIT running -6 to 6 arcsec."""
    slit = np.tile(np.linspace(-6.0, 6.0, MAP_SHAPE[1]), (MAP_SHAPE[0], 1))
    orderMap = np.where(np.arange(MAP_SHAPE[0])[:, None] < MAP_SHAPE[0] // 2, 11.0, 12.0) * np.ones(MAP_SHAPE)
    wavelength = np.full(MAP_SHAPE, 600.0)
    return write_map_images(tmpPath / f"DISP_IMAGE_{arm}.fits", wavelength, slit, orderMap)


def _illumination_result(*, isApplied: bool = True, hasProfiles: bool = True) -> slit_illumination_result:
    status = STATUS_FITTED if isApplied else STATUS_FAILED
    offsets = pd.DataFrame(
        [
            {
                "order": order,
                "offset": OFFSET_ARCSEC[order] if isApplied else np.nan,
                "offsetErr": 0.02,
                "scienceWidth": 0.15,
                "nBins": 100,
                "status": status,
                "reason": "" if isApplied else "too few skyline bins",
                "borrowedFrom": np.nan,
            }
            for order in OFFSET_ARCSEC
        ],
        columns=OFFSET_COLUMNS,
    )
    flatEdges = pd.DataFrame(
        [{"order": order, "lower": -5.0, "upper": 5.0, "width": 0.15, "reason": ""} for order in OFFSET_ARCSEC],
        columns=FLAT_EDGE_COLUMNS,
    )
    slit = np.linspace(-6.0, 6.0, 40)
    profiles = pd.DataFrame(
        {"order": 11, "slit": slit, "profile": 1.0, "profileErr": 0.1, "nLines": 5}, columns=PROFILE_COLUMNS
    )
    return slit_illumination_result(
        offsets=offsets,
        flatEdges=flatEdges,
        profiles=profiles if hasProfiles else pd.DataFrame(columns=PROFILE_COLUMNS),
        isApplied=isApplied,
        reason="" if isApplied else "no order has a usable offset",
        overlapMarginArcsec=0.45,
    )


def _stub_measurement(
    monkeypatch: pytest.MonkeyPatch, result: slit_illumination_result | None
) -> list[dict[str, object]]:
    """Replace `measure_slit_illumination` with a double that records its arguments; None makes any call a failure."""
    calls: list[dict[str, object]] = []

    def fake_measure(
        log: Any,
        frame: Any,
        masterFlat: Any,
        twoDMapPath: Any,
        settings: Any,
        arm: Any,
        binx: Any,
        biny: Any,
        illuminationSettings: Any = None,
        partners: Any = None,
    ) -> slit_illumination_result:
        calls.append(
            {
                "frame": frame,
                "masterFlat": masterFlat,
                "twoDMapPath": twoDMapPath,
                "arm": arm,
                "binx": binx,
                "biny": biny,
                "illuminationSettings": illuminationSettings,
                "partners": partners,
            }
        )
        assert result is not None, "the slit-illumination measurement must not run"
        return result

    monkeypatch.setattr(
        import_module("soxspipe.commonutils.slit_illumination"), "measure_slit_illumination", fake_measure
    )
    return calls


def _illumination_qc(recipe: soxs_stare) -> pd.DataFrame:
    return recipe.qc.loc[recipe.qc["qc_name"].str.startswith("SLIT ILLUM")].set_index("qc_name")


def test_without_a_slit_illumination_block_nothing_is_measured_and_the_sky_gets_no_mask(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Develop's behaviour: no block means no measurement, no QC, no plot and `excludeMask=None`."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path, recipeSettings=_illumination_recipe_settings(None))
    calls = _stub_measurement(monkeypatch, None)
    captured = _patch_stare(
        recipe, monkeypatch, tmp_path, skyResult=tuple(synthetic_ccd(seed=540 + n, prepared=True) for n in range(3))
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls == []
    assert captured["sky"][0]["excludeMask"] is None
    assert len(_illumination_qc(recipe)) == 0
    assert not recipe.products["product_label"].eq("SLIT_ILLUMINATION_QC_PLOT").any()


@pytest.mark.parametrize(
    ("measureOffset", "isWarned"),
    [(False, False), ("yes", True), (1, True), (None, True)],
    ids=["false", "string", "int", "none"],
)
def test_measure_offset_is_off_unless_it_is_the_boolean_true(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    measureOffset: object,
    isWarned: bool,
) -> None:
    """Only `measure-offset: True` switches the measurement on; any other value is off, and a non-bool warns."""
    # ARRANGE
    block = {**ILLUMINATION_BLOCK, "measure-offset": measureOffset}
    recipe = _stare_recipe(log, tmp_path, recipeSettings=_illumination_recipe_settings(block))
    calls = _stub_measurement(monkeypatch, None)
    captured = _patch_stare(
        recipe, monkeypatch, tmp_path, skyResult=tuple(synthetic_ccd(seed=540 + n, prepared=True) for n in range(3))
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls == []
    assert captured["sky"][0]["excludeMask"] is None
    assert len(_illumination_qc(recipe)) == 0
    warnings = [text for level, text in log.messages if level == "warning" and "measure-offset" in text]
    assert bool(warnings) is isWarned


def test_a_measured_offset_masks_the_sky_model_and_records_the_qcs_and_plot(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With the block on, the measurement gets the right inputs, its mask reaches the sky, and QCs and plot are kept."""
    # ARRANGE
    routes = {_route(PRO_CATG="MASTER_FLAT_VIS"): [str(_calibration_path(tmp_path, "MASTER_FLAT_VIS", FLAT_SEED))]}
    recipe = _stare_recipe(
        log, tmp_path, extraRoutes=routes, recipeSettings=_illumination_recipe_settings(ILLUMINATION_BLOCK)
    )
    mapPath = _write_illumination_map(tmp_path)
    Path(recipe.qcDir).mkdir()
    result = _illumination_result()
    calls = _stub_measurement(monkeypatch, result)
    detrended = synthetic_ccd(seed=550, prepared=True, headerOverrides={"WIN_BINX": 1, "WIN_BINY": 1})
    captured = _patch_stare(
        recipe,
        monkeypatch,
        tmp_path,
        detrended=detrended,
        skyResult=tuple(synthetic_ccd(seed=560 + n, prepared=True) for n in range(3)),
    )
    _, slit, orderMap = _read_disp_map(mapPath, 1, 1, recipe.kw)
    expectedMask = overlap_exclude_mask(log, slit, orderMap, result.flatEdges, result.offsets, 0.3)

    # ACT
    recipe.produce_product()

    # ASSERT
    [call] = calls
    assert call["frame"] is detrended
    _assert_read_of(call["masterFlat"], FLAT_SEED)
    assert (call["twoDMapPath"], call["arm"], call["binx"], call["biny"]) == (str(mapPath), "VIS", 1, 1)
    assert call["illuminationSettings"] == ILLUMINATION_BLOCK
    assert call["partners"] is None
    assert expectedMask.any()
    np.testing.assert_array_equal(captured["sky"][0]["excludeMask"], expectedMask)
    qcs = _illumination_qc(recipe)
    assert set(qcs.index) == set(ILLUMINATION_QC_NAMES)
    assert qcs.loc["SLIT ILLUM OFFSET O11", "qc_value"] == 0.8
    assert qcs.loc["SLIT ILLUM OFFSET O12", "qc_value"] == -0.4
    assert qcs.loc["SLIT ILLUM OFFSET O11", "qc_unit"] == "arcsec"
    assert qcs.loc["SLIT ILLUM APPLIED", "qc_value"] == 1
    assert qcs.loc["SLIT ILLUM N ORDERS", "qc_value"] == 2
    assert qcs["to_header"].all()
    plotPath = Path(recipe.qcDir) / "OBJECT_VIS_SLIT_ILLUM_QC_PLOT.pdf"
    assert plotPath.read_bytes().startswith(b"%PDF")
    [row] = recipe.products.loc[recipe.products["product_label"] == "SLIT_ILLUMINATION_QC_PLOT"].to_dict("records")
    assert (row["file_name"], row["file_path"], row["file_type"], row["label"]) == (
        plotPath.name,
        str(plotPath),
        "PDF",
        "QC",
    )


def test_a_master_flat_that_was_switched_off_is_measured_as_no_flat(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`use_flat` false leaves `master_flat` False, which must reach the measurement as None, never False."""
    # ARRANGE
    routes = {_route(PRO_CATG="MASTER_FLAT_VIS"): [str(_calibration_path(tmp_path, "MASTER_FLAT_VIS", FLAT_SEED))]}
    recipeSettings = _illumination_recipe_settings(ILLUMINATION_BLOCK, useFlat=False)
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes, recipeSettings=recipeSettings)
    _write_illumination_map(tmp_path)
    calls = _stub_measurement(monkeypatch, _illumination_result(isApplied=False, hasProfiles=False))
    _patch_stare(
        recipe, monkeypatch, tmp_path, skyResult=tuple(synthetic_ccd(seed=570 + n, prepared=True) for n in range(3))
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls[0]["masterFlat"] is None


@pytest.mark.parametrize(
    ("arm", "headerBinning", "expected"),
    [
        ("VIS", {"WIN_BINX": 2, "WIN_BINY": 4}, (2, 4)),
        ("VIS", {}, (1, 1)),
        ("NIR", {"WIN_BINX": 2, "WIN_BINY": 2}, (1, 1)),
    ],
    ids=["vis-binned", "vis-no-keyword", "nir-ignores-keyword"],
)
def test_the_measurement_uses_the_header_binning_except_on_the_nir_arm_or_without_the_keyword(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    arm: str,
    headerBinning: dict[str, object],
    expected: tuple[int, int],
) -> None:
    """Binning comes from `WIN_BINX`/`WIN_BINY` of the detrended frame, as `soxs_order_centres` reads it."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path, arm=arm, recipeSettings=_illumination_recipe_settings(ILLUMINATION_BLOCK))
    calls = _stub_measurement(monkeypatch, _illumination_result(isApplied=False, hasProfiles=False))
    _patch_stare(
        recipe,
        monkeypatch,
        tmp_path,
        detrended=synthetic_ccd(seed=580, arm=arm, prepared=True, headerOverrides=headerBinning),
        skyResult=tuple(synthetic_ccd(seed=581 + n, prepared=True) for n in range(3)),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert (calls[0]["binx"], calls[0]["biny"]) == expected


def test_a_measurement_that_is_not_applied_records_zero_and_gives_the_sky_no_mask(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A not-applied result writes APPLIED 0 and N ORDERS 0, no offset QCs and `excludeMask=None`; the plot stays."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path, recipeSettings=_illumination_recipe_settings(ILLUMINATION_BLOCK))
    _write_illumination_map(tmp_path)
    Path(recipe.qcDir).mkdir()
    _stub_measurement(monkeypatch, _illumination_result(isApplied=False))
    captured = _patch_stare(
        recipe, monkeypatch, tmp_path, skyResult=tuple(synthetic_ccd(seed=590 + n, prepared=True) for n in range(3))
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["sky"][0]["excludeMask"] is None
    qcs = _illumination_qc(recipe)
    assert set(qcs.index) == {"SLIT ILLUM APPLIED", "SLIT ILLUM N ORDERS"}
    assert qcs.loc["SLIT ILLUM APPLIED", "qc_value"] == 0
    assert qcs.loc["SLIT ILLUM N ORDERS", "qc_value"] == 0
    assert (Path(recipe.qcDir) / "OBJECT_VIS_SLIT_ILLUM_QC_PLOT.pdf").exists()


def test_no_plot_is_written_when_the_measurement_has_no_profiles(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A measurement that never got as far as profiles (no flat, unreadable map) has nothing to plot."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path, recipeSettings=_illumination_recipe_settings(ILLUMINATION_BLOCK))
    _stub_measurement(monkeypatch, _illumination_result(isApplied=False, hasProfiles=False))
    _patch_stare(
        recipe, monkeypatch, tmp_path, skyResult=tuple(synthetic_ccd(seed=595 + n, prepared=True) for n in range(3))
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert not recipe.products["product_label"].eq("SLIT_ILLUMINATION_QC_PLOT").any()
    assert not (Path(recipe.qcDir) / "OBJECT_VIS_SLIT_ILLUM_QC_PLOT.pdf").exists()


def test_a_wide_rotator_spread_between_input_frames_is_reported(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Frames whose ABSROT differs by more than the limit get an info note; the stacked frame's ABSROT is logged."""
    # ARRANGE
    objectRoute = _route(DPR_TYPE="OBJECT", DPR_TECH=STARE_TECHNIQUE)
    paths = [
        str(
            _write_prepared_frame(
                tmp_path / f"stare-object-{n}.fits",
                seed=600 + n,
                headerOverrides={"DPR_TYPE": "OBJECT", "DPR_TECH": STARE_TECHNIQUE, "ABSROT": angle},
            )
        )
        for n, angle in enumerate((175.0, -175.0, 150.0))
    ]
    recipe = _stare_recipe(
        log,
        tmp_path,
        extraRoutes={objectRoute: paths},
        recipeSettings=_illumination_recipe_settings(ILLUMINATION_BLOCK),
    )
    _stub_measurement(monkeypatch, _illumination_result(isApplied=False, hasProfiles=False))
    _patch_stare(
        recipe,
        monkeypatch,
        tmp_path,
        detrended=synthetic_ccd(seed=610, prepared=True, headerOverrides={"ABSROT": 172.5}),
        skyResult=tuple(synthetic_ccd(seed=611 + n, prepared=True) for n in range(3)),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    spreads = [
        text for level, text in log.messages if level == "info" and "ABSROT of the stacked input frames spans" in text
    ]
    assert len(spreads) == 1
    assert "35.0" in spreads[0]
    assert any(level == "info" and "172.5" in text and "ABSROT" in text for level, text in log.messages)


def test_a_narrow_rotator_spread_or_missing_rotator_keywords_raise_no_warning(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Frames within the limit, and frames with no ABSROT at all, neither warn nor fail the recipe."""
    # ARRANGE
    objectRoute = _route(DPR_TYPE="OBJECT", DPR_TECH=STARE_TECHNIQUE)
    overrides = [{"ABSROT": 179.0}, {"ABSROT": -179.0}, {}]
    paths = [
        str(
            _write_prepared_frame(
                tmp_path / f"stare-object-{n}.fits",
                seed=620 + n,
                headerOverrides={"DPR_TYPE": "OBJECT", "DPR_TECH": STARE_TECHNIQUE, **override},
            )
        )
        for n, override in enumerate(overrides)
    ]
    recipe = _stare_recipe(
        log,
        tmp_path,
        extraRoutes={objectRoute: paths},
        recipeSettings=_illumination_recipe_settings(ILLUMINATION_BLOCK),
    )
    _stub_measurement(monkeypatch, _illumination_result(isApplied=False, hasProfiles=False))
    _patch_stare(
        recipe, monkeypatch, tmp_path, skyResult=tuple(synthetic_ccd(seed=630 + n, prepared=True) for n in range(3))
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert not [text for _, text in log.messages if "ABSROT of the stacked input frames spans" in text]


def test_the_slit_illumination_qcs_are_never_checked_against_an_acceptable_range(log: Any) -> None:
    """Flexure must not fail a science recipe: a QC with no `qc-acceptable-ranges` entry always passes."""
    # ARRANGE
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.inst = "XSH"
    recipe.recipeName = "soxs-stare"
    recipe.dateObs = "2024-01-02T03:04:05.678"
    recipe.recipeSettings = {"qc-acceptable-ranges": {"RON": [0, 1]}}
    recipe.qc = qc_table()
    for name, value in (("SLIT ILLUM OFFSET O11", 2.9), ("SLIT ILLUM APPLIED", 0), ("SLIT ILLUM N ORDERS", 0)):
        recipe.add_qc(qcName=name, qcValue=value, qcComment="x", reductionDateUtc="2024-01-02T04:05:06")

    # ACT
    recipe.flag_poor_data()

    # ASSERT
    flags = recipe.qc.set_index("qc_name")["qc_flag"]
    assert flags["RON"] == "fail"
    assert flags[[name for name in flags.index if name.startswith("SLIT ILLUM")]].eq("pass").all()


@pytest.mark.parametrize("subtractSky", [False, True], ids=["sky-off", "no-model"])
def test_with_no_sky_model_both_extraction_frames_carry_the_overlap_mask(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    subtractSky: bool,
) -> None:
    """Sky off, or a sky model that comes back None, still masks the skyline-free slit ends in both Horne frames."""
    # ARRANGE
    recipeSettings = {
        "use_flat": False,
        "sky-subtraction": {"subtract_sky": subtractSky},
        "slit-illumination": ILLUMINATION_BLOCK,
    }
    recipe = _stare_recipe(log, tmp_path, recipeSettings=recipeSettings)
    mapPath = _write_illumination_map(tmp_path)
    Path(recipe.qcDir).mkdir()
    result = _illumination_result()
    _stub_measurement(monkeypatch, result)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=(None, None, None))
    _, slit, orderMap = _read_disp_map(mapPath, 1, 1, recipe.kw)
    excludeMask = overlap_exclude_mask(log, slit, orderMap, result.flatEdges, result.offsets, 0.3)
    stackBefore = synthetic_ccd(seed=510, prepared=True)

    # ACT
    recipe.produce_product()

    # ASSERT
    stack = captured["stack"][0]["result"]
    detrended = captured["detrend"][0]["result"]
    extractorArgs = captured["extractor"][0]
    expectedUnflattened = np.asarray(stackBefore.mask, dtype=bool) | excludeMask
    expectedSkySubtracted = np.asarray(detrended.mask, dtype=bool) | excludeMask
    assert excludeMask.any()
    assert extractorArgs["subtractedFrame"] is False
    np.testing.assert_array_equal(extractorArgs["skySubtractedFrame"].mask, expectedSkySubtracted)
    np.testing.assert_array_equal(extractorArgs["unflattenedFrame"].mask, expectedUnflattened)
    # NEITHER INPUT FRAME IS MUTATED
    assert extractorArgs["unflattenedFrame"] is not stack
    assert extractorArgs["skySubtractedFrame"] is not detrended
    np.testing.assert_array_equal(stack.mask, stackBefore.mask)
    np.testing.assert_array_equal(detrended.mask, synthetic_ccd(seed=511, prepared=True).mask)


def test_with_no_sky_model_and_no_applied_offset_the_extraction_frames_are_the_input_frames(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A measurement that is not applied gives no mask, so the no-model frames stay the very objects develop used."""
    # ARRANGE
    recipeSettings = {
        "use_flat": False,
        "sky-subtraction": {"subtract_sky": True},
        "slit-illumination": ILLUMINATION_BLOCK,
    }
    recipe = _stare_recipe(log, tmp_path, recipeSettings=recipeSettings)
    _stub_measurement(monkeypatch, _illumination_result(isApplied=False, hasProfiles=False))
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=(None, None, None))

    # ACT
    recipe.produce_product()

    # ASSERT
    extractorArgs = captured["extractor"][0]
    assert extractorArgs["unflattenedFrame"] is captured["stack"][0]["result"]
    assert extractorArgs["skySubtractedFrame"] is captured["detrend"][0]["result"]


NO_TRACE_STANDARD_MESSAGE = (
    "No object trace was found in the standard-star frame, so no response curve was made. "
    "The spectrum is a robust boxcar extraction; check the pointing and conditions."
)


def _spy_response_and_plot(recipe: soxs_stare, monkeypatch: pytest.MonkeyPatch) -> dict[str, list[tuple[Any, ...]]]:
    """Replace the response-curve, flux-calibration and merged-plot steps with call recorders."""
    calls: dict[str, list[tuple[Any, ...]]] = {"response": [], "fluxCal": [], "plot": []}
    monkeypatch.setattr(recipe, "_generate_response_curve", lambda *args, **_: calls["response"].append(args) or False)
    monkeypatch.setattr(
        recipe, "_flux_calibrate_spectrum", lambda *args: calls["fluxCal"].append(args) or str(Path("fluxcal.fits"))
    )
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: calls["plot"].append(args))
    return calls


def test_the_stare_extractor_is_asked_for_the_no_trace_boxcar_fallback(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    recipe.produce_product()

    assert captured["extractor"][0]["noTraceFallback"] is True


def test_a_standard_star_without_an_object_trace_makes_no_response_curve_and_fails_the_recipe(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path)
    recipe.generateReponseCurve = True
    captured = _patch_stare(recipe, monkeypatch, tmp_path, traceFallback=True)
    calls = _spy_response_and_plot(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls["response"] == []
    assert captured["cleanUp"] == [{"forceFail": NO_TRACE_STANDARD_MESSAGE}]
    assert ("error", NO_TRACE_STANDARD_MESSAGE) in log.messages


def test_a_standard_star_without_an_object_trace_still_plots_the_merged_spectrum_with_the_failure_note(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _stare_recipe(log, tmp_path)
    recipe.generateReponseCurve = True
    _patch_stare(recipe, monkeypatch, tmp_path, traceFallback=True)
    calls = _spy_response_and_plot(recipe, monkeypatch)

    recipe.produce_product()

    assert len(calls["plot"]) == 1
    mergedSpectrum, orderJoins, filePathFluxcal, notes, isBoxcar = calls["plot"][0]
    assert mergedSpectrum["WAVE"].tolist() == [500.0]
    assert filePathFluxcal is None
    assert ("ERROR", NO_TRACE_STANDARD_MESSAGE) in notes
    assert isBoxcar is True


def test_a_standard_star_with_an_object_trace_still_builds_the_response_curve(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _stare_recipe(log, tmp_path)
    recipe.generateReponseCurve = True
    captured = _patch_stare(recipe, monkeypatch, tmp_path)
    calls = _spy_response_and_plot(recipe, monkeypatch)

    recipe.produce_product()

    assert len(calls["response"]) == 1
    assert captured["cleanUp"] == [{"forceFail": False}]
    assert calls["plot"][0][3] == []
    assert calls["plot"][0][4] is False


def test_a_science_frame_without_an_object_trace_is_still_flux_calibrated_when_a_response_exists(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    routes = {_route(PRO_CATG="RESP_TAB_VIS"): [str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))]}
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, traceFallback=True)
    calls = _spy_response_and_plot(recipe, monkeypatch)

    recipe.produce_product()

    assert len(calls["fluxCal"]) == 1
    assert calls["fluxCal"][0][0] == str(tmp_path / "RESP_TAB_VIS.fits")
    assert calls["response"] == []
    assert captured["cleanUp"] == [{"forceFail": False}]
    assert calls["plot"][0][2] == "fluxcal.fits"


def test_a_standard_star_without_an_object_trace_is_flux_calibrated_and_not_failed_when_a_response_is_supplied(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A supplied response table is applied to the boxcar like it is to a science frame (no response curve is built)."""
    routes = {_route(PRO_CATG="RESP_TAB_VIS"): [str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))]}
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes)
    recipe.generateReponseCurve = True
    captured = _patch_stare(recipe, monkeypatch, tmp_path, traceFallback=True)
    calls = _spy_response_and_plot(recipe, monkeypatch)

    recipe.produce_product()

    assert len(calls["fluxCal"]) == 1
    assert calls["response"] == []
    assert captured["cleanUp"] == [{"forceFail": False}]
    assert calls["plot"][0][3] == []


def test_a_science_frame_without_an_object_trace_or_a_response_does_not_fail_the_recipe(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, traceFallback=True)
    calls = _spy_response_and_plot(recipe, monkeypatch)

    recipe.produce_product()

    assert captured["cleanUp"] == [{"forceFail": False}]
    assert calls["plot"][0][3] == []


@pytest.mark.parametrize("isBoxcar", [False, True])
def test_both_merged_spectrum_plots_are_told_whether_the_spectrum_is_a_boxcar(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    isBoxcar: bool,
) -> None:
    from astropy.table import Table

    recipe = _stare_recipe(log, tmp_path)
    recipe.dateObs = "2024-01-02T03:04:05"
    toolkit = import_module("soxspipe.commonutils.toolkit")
    plotKwargs: list[dict[str, object]] = []
    monkeypatch.setattr(
        toolkit,
        "plot_merged_spectrum_qc",
        lambda **kwargs: plotKwargs.append(kwargs) or (kwargs["products"], str(tmp_path / "merged.pdf")),
    )
    fluxCalPath = tmp_path / "fluxcal.fits"
    Table({"WAVE": [500.0], "FLUX_CALIBRATED": [1.0]}).write(fluxCalPath)
    merged = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0], "SKY_COUNTS": [1.0]})

    recipe._plot_merged_spectrum_qcs(merged, {"10-11": 11}, str(fluxCalPath), [], isBoxcar)

    assert [kwargs["boxcarExtraction"] for kwargs in plotKwargs] == [isBoxcar, isBoxcar]


def test_the_measured_slit_offsets_reach_the_extractor_for_the_faint_order_trace(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The usable offsets of the slit-illumination measurement go to the Horne extraction as `slitOffsetsArcsec`."""
    # ARRANGE
    routes = {_route(PRO_CATG="MASTER_FLAT_VIS"): [str(_calibration_path(tmp_path, "MASTER_FLAT_VIS", FLAT_SEED))]}
    recipe = _stare_recipe(
        log, tmp_path, extraRoutes=routes, recipeSettings=_illumination_recipe_settings(ILLUMINATION_BLOCK)
    )
    _write_illumination_map(tmp_path)
    Path(recipe.qcDir).mkdir()
    _stub_measurement(monkeypatch, _illumination_result())
    captured = _patch_stare(
        recipe,
        monkeypatch,
        tmp_path,
        detrended=synthetic_ccd(seed=550, prepared=True, headerOverrides={"WIN_BINX": 1, "WIN_BINY": 1}),
        skyResult=tuple(synthetic_ccd(seed=560 + n, prepared=True) for n in range(3)),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["extractor"][0]["slitOffsetsArcsec"] == OFFSET_ARCSEC


def test_a_frame_whose_offsets_are_not_measured_hands_the_extractor_no_offsets_from_an_earlier_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Per-frame mode reuses the recipe object, so a stale offset dictionary must not outlive its frame."""
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path, recipeSettings=_illumination_recipe_settings(None))
    recipe.slitIlluminationOffsets = {11: 9.9}
    captured = _patch_stare(
        recipe, monkeypatch, tmp_path, skyResult=tuple(synthetic_ccd(seed=540 + n, prepared=True) for n in range(3))
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["extractor"][0]["slitOffsetsArcsec"] is None
    assert recipe.slitIlluminationOffsets is None
