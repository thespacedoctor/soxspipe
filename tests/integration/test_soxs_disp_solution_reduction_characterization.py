"""Characterization of the `soxs_disp_solution.produce_product` branches the orchestration test misses.

`tests/integration/test_mapping_recipe_orchestration.py` already covers one
path: a SOXS VIS reduction with a master bias, no dark, `tune-pipeline` off and
no `polyOrders`. That leaves the whole back half of `produce_product`
uncovered -- lines 251, 265, 272, 288, 307-315, 318-371 and 375-378 at the
branch point, the largest single block being the 54-line `tune-pipeline`
branch.

These tests pin that region before anything restructures it: the master-dark
read, the NIR lamp-off dark, the X-Shooter filter variants, the last-frame-wins
precedence of the frame loops, the intermediate-product write, the parameter
tuning branch and the `polyOrders` override. They reuse the orchestration
module's routing collection so both files build the recipe the same way.
"""

from __future__ import annotations

import re
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_disp_solution
from soxspipe.recipes.base_recipe import base_recipe
from tests.factories import pipeline_settings, prepared_fits, qc_row, qc_table
from tests.integration.test_mapping_recipe_orchestration import (
    RoutedFrameCollection,
    _empty_products,
    _route,
)

pytestmark = pytest.mark.integration

# THE FIVE DIGITS THE TUNING GRID PERMUTES, FOUR AT A TIME.
TUNING_DIGITS = [2, 3, 4, 5, 6]
TUNING_PERMUTATIONS = len(TUNING_DIGITS) ** 4

# THE RENDERED SHAPE OF THE REDUCTION TIMESTAMP, WHOLE SECONDS AND NO OFFSET.
TIMESTAMP_PATTERN = r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}"

DISPERSION_MODULE = import_module("soxspipe.recipes.soxs_disp_solution")


def _declared_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the empty QC and product tables a real recipe starts from.

    The integration factories declare a different column order from
    `base_recipe`, so a test that pins column order or cell dtype has to start
    from the recipe's own declaration.
    """
    return base_recipe._empty_qc_and_product_tables(base_recipe, pd)


def _configure_recipe(
    recipe: soxs_disp_solution,
    log: Any,
    tmpPath: Path,
    *,
    arm: str,
    instrument: str,
    routes: dict[frozenset[tuple[str, object]], tuple[str, ...]],
    settingOverrides: dict[str, Any] | None = None,
    polyOrders: Any = False,
) -> None:
    """Populate an unconstructed recipe with everything `produce_product` reads."""
    recipe.log = log
    recipe.arm = arm
    recipe.inst = instrument
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RoutedFrameCollection(filePathsByFilters=routes)
    recipe.recipeName = "soxs-disp-solution"
    recipe.settings = pipeline_settings(
        tmpPath,
        overrides={"tune-pipeline": False, **(settingOverrides or {})},
    )
    recipe.recipeSettings = {}
    recipe.workspaceRootPath = str(tmpPath)
    recipe.sofName = "synthetic-dispersion"
    recipe.startNightDate = "2024-01-02"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.polyOrders = polyOrders
    recipe.qc = qc_table()
    recipe.products = _empty_products()


def _patch_reduction(
    recipe: soxs_disp_solution,
    monkeypatch: pytest.MonkeyPatch,
    *,
    productPath: Path,
    fittedQc: pd.DataFrame,
) -> tuple[list[str], dict[str, Any]]:
    """Stub the collaborators outside the recipe and return the call log and captures."""
    calls: list[str] = []
    captured: dict[str, Any] = {}

    def fake_detrend(**kwargs: Any) -> Any:
        captured["detrend"] = kwargs
        calls.append("detrend")
        return kwargs["inputFrame"].copy()

    monkeypatch.setattr(recipe, "detrend", fake_detrend)
    monkeypatch.setattr(
        recipe,
        "update_fits_keywords",
        lambda **kwargs: calls.append("keywords"),
    )

    class FakeDispersionMap:
        def __init__(self, **kwargs: object) -> None:
            captured.setdefault("mapArguments", []).append(kwargs)
            calls.append("create_dispersion_map")

        def get(self) -> tuple[object, ...]:
            calls.append("map_get")
            return (
                str(productPath),
                None,
                (),
                fittedQc,
                _empty_products(),
                "LINE-DETECTION-TABLE",
            )

    monkeypatch.setattr(DISPERSION_MODULE, "create_dispersion_map", FakeDispersionMap)
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))
    return calls, captured


def _fitted_qc() -> pd.DataFrame:
    """Return the QC row the stubbed dispersion map hands back."""
    return qc_row(
        recipeName="soxs-disp-solution",
        name="DISPERSION_RESIDUAL",
        value=0.012,
        unit="pixel",
        comment="Synthetic dispersion-fit residual",
    )


# ---------------------------------------------------------------------------
# THE CALIBRATION AND PINHOLE FRAME READS
# ---------------------------------------------------------------------------


def test_a_master_dark_is_read_and_reaches_the_detrend(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A `MASTER_DARK` for the arm is read and handed to `detrend`."""
    # ARRANGE
    biasPath = prepared_fits(tmp_path / "MASTER_BIAS_VIS.fits", seed=11)
    darkPath = prepared_fits(tmp_path / "MASTER_DARK_VIS.fits", seed=13)
    pinholePath = prepared_fits(tmp_path / "PINHOLE_VIS.fits", seed=12)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(PRO_CATG="MASTER_BIAS_VIS"): (str(biasPath),),
            _route(PRO_CATG="MASTER_DARK_VIS"): (str(darkPath),),
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
    )
    calls, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["dark"] is not False
    assert captured["detrend"]["master_bias"] is not False
    assert calls[0] == "detrend"


def test_the_nir_lamp_off_frame_overrides_a_master_dark(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The lamp-off `IMAGE` read runs after the master-dark read, so it wins."""
    # ARRANGE
    darkPath = prepared_fits(tmp_path / "MASTER_DARK_NIR.fits", seed=13)
    lampOffPath = prepared_fits(tmp_path / "LAMP_OFF_NIR.fits", seed=14)
    pinholePath = prepared_fits(tmp_path / "PINHOLE_NIR.fits", seed=12)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="NIR",
        instrument="SOXS",
        routes={
            _route(PRO_CATG="MASTER_DARK_NIR"): (str(darkPath),),
            _route(DPR_TYPE="WAVE,LAMP", DPR_TECH="IMAGE"): (str(lampOffPath),),
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
    )
    calls, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_NIR.fits",
        fittedQc=_fitted_qc(),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["master_bias"] is False
    assert captured["detrend"]["dark"].data[0][0] == pytest.approx(
        _read_first_pixel(lampOffPath), rel=1e-12
    )
    assert "create_dispersion_map" in calls


def test_the_last_matching_pinhole_frame_is_the_one_reduced(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The SOXS filter list is walked in order, so a `WAVE,LAMP` frame wins over `LAMP,WAVE`."""
    # ARRANGE
    firstPath = prepared_fits(tmp_path / "PINHOLE_A_VIS.fits", seed=21)
    secondPath = prepared_fits(tmp_path / "PINHOLE_B_VIS.fits", seed=22)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(firstPath),),
            _route(DPR_TYPE="WAVE,LAMP", DPR_TECH="ECHELLE,PINHOLE"): (str(secondPath),),
        },
    )
    _, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["inputFrame"].data[0][0] == pytest.approx(
        _read_first_pixel(secondPath), rel=1e-12
    )


def test_an_xshooter_reduction_routes_the_fmtchk_filters(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A non-SOXS instrument reads `LAMP,FMTCHK` for both the dark and the pinhole frame."""
    # ARRANGE
    lampOffPath = prepared_fits(tmp_path / "XSH_LAMP_OFF.fits", seed=31)
    pinholePath = prepared_fits(tmp_path / "XSH_PINHOLE.fits", seed=32)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="NIR",
        instrument="XSH",
        routes={
            _route(DPR_TYPE="LAMP,FMTCHK", DPR_TECH="IMAGE"): (str(lampOffPath),),
            _route(DPR_TYPE="LAMP,FMTCHK", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
    )
    _, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_NIR.fits",
        fittedQc=_fitted_qc(),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["inputFrame"].data[0][0] == pytest.approx(
        _read_first_pixel(pinholePath), rel=1e-12
    )
    assert captured["detrend"]["dark"].data[0][0] == pytest.approx(
        _read_first_pixel(lampOffPath), rel=1e-12
    )
    routedFilters = [filters for name, filters in recipe.inputFrames.calls if name == "files_filtered"]
    assert {"DPR_TYPE": "LAMP,FMTCHK", "DPR_TECH": "IMAGE"} in routedFilters
    assert {"DPR_TYPE": "LAMP,FMTCHK", "DPR_TECH": "ECHELLE,PINHOLE"} in routedFilters
    assert {"DPR_TYPE": "WAVE,LAMP", "DPR_TECH": "IMAGE"} not in routedFilters


def _read_first_pixel(path: Path) -> float:
    """Return the first flux pixel of a prepared frame, for identity assertions."""
    from astropy import units as u
    from astropy.nddata import CCDData

    frame = CCDData.read(
        str(path),
        hdu=0,
        unit=u.electron,
        hdu_uncertainty="ERRS",
        hdu_mask="QUAL",
        hdu_flags="FLAGS",
        key_uncertainty_type="UTYPE",
    )
    return float(frame.data[0][0])


# ---------------------------------------------------------------------------
# THE INTERMEDIATE PRODUCT WRITE
# ---------------------------------------------------------------------------


def test_saving_intermediate_products_writes_the_calibrated_pinhole_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With the setting on, the calibrated frame is written and its path printed."""
    # ARRANGE
    pinholePath = prepared_fits(tmp_path / "PINHOLE_VIS.fits", seed=12)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
        settingOverrides={"save-intermediate-products": True},
    )
    calls, _ = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
    )
    writeArguments: dict[str, Any] = {}

    def fake_write(**kwargs: object) -> str:
        writeArguments.update(kwargs)
        calls.append("write")
        return str(tmp_path / "PINHOLE_VIS_pre.fits")

    monkeypatch.setattr(recipe, "_write", fake_write)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert writeArguments["frame"] is recipe.pinholeFrame
    assert writeArguments["filedir"] == str(tmp_path)
    assert writeArguments["filename"] is False
    assert writeArguments["overwrite"] is True
    assert writeArguments["product"] is False
    assert calls.index("write") < calls.index("create_dispersion_map")
    printed = [message for level, message in log.messages if level == "print"]
    assert f"\nCalibrated single pinhole frame: {tmp_path / 'PINHOLE_VIS_pre.fits'}\n" in printed


# ---------------------------------------------------------------------------
# THE PARAMETER TUNING BRANCH
# ---------------------------------------------------------------------------


def _tuning_recipe(
    log: Any,
    tmpPath: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[soxs_disp_solution, list[str], dict[str, Any], dict[str, Any]]:
    """Return a recipe configured for the `tune-pipeline` branch, with tuning stubbed."""
    pinholePath = prepared_fits(tmpPath / "PINHOLE_VIS.fits", seed=12)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmpPath,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
        settingOverrides={"tune-pipeline": True},
    )
    calls, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmpPath / "DISP_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
    )

    tuningArguments: dict[str, Any] = {}

    def fake_fmultiprocess(**kwargs: object) -> list[object]:
        tuningArguments.update(kwargs)
        calls.append("fmultiprocess")
        return []

    fundamentals = import_module("fundamentals")
    monkeypatch.setattr(fundamentals, "fmultiprocess", fake_fmultiprocess)
    return recipe, calls, captured, tuningArguments


def test_tuning_fits_the_line_list_once_then_sweeps_the_polynomial_grid(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The tuning branch fits once for the line list, then tunes over 625 permutations."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-run"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    (tuningDirectory / "residuals.txt").write_text("stale residuals\n")
    recipe, calls, captured, tuningArguments = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    returnedPath, returnedQc = recipe.produce_product()

    # ASSERT
    assert returnedPath is None
    pd.testing.assert_frame_equal(returnedQc, _fitted_qc())
    assert calls == ["detrend", "keywords", "create_dispersion_map", "map_get", "fmultiprocess"]
    assert not (tuningDirectory / "residuals.txt").exists()

    assert tuningArguments["function"] is DISPERSION_MODULE.parameterTuning
    assert len(tuningArguments["inputArray"]) == TUNING_PERMUTATIONS
    assert tuningArguments["inputArray"][0] == (2, 2, 2, 2)
    assert tuningArguments["inputArray"][-1] == (6, 6, 6, 6)
    assert tuningArguments["poolSize"] == 100
    assert tuningArguments["timeout"] == 3600
    assert tuningArguments["turnOffMP"] is recipe.debug
    assert tuningArguments["mute"] is True
    assert tuningArguments["progressBar"] is True
    assert tuningArguments["lineDetectionTable"] == "LINE-DETECTION-TABLE"
    assert tuningArguments["pinholeFrame"] is recipe.pinholeFrame


def test_tuning_leaves_the_products_and_quality_control_tables_untouched(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No product row is recorded and neither `report_output` nor `clean_up` runs."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-untouched"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    recipe, calls, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)
    originalProducts = recipe.products.copy(deep=True)
    originalQc = recipe.qc.copy(deep=True)

    # ACT
    recipe.produce_product()

    # ASSERT
    pd.testing.assert_frame_equal(recipe.products, originalProducts)
    pd.testing.assert_frame_equal(recipe.qc, originalQc)
    assert "report" not in calls
    assert "clean_up" not in calls


def test_tuning_without_a_stale_residuals_file_records_the_failed_removal(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A missing `residuals.txt` is swallowed and logged rather than raised."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-clean"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    recipe, _, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    debugged = [message for level, message in log.messages if level == "debug"]
    assert any("`os.remove('residuals.txt')` failed" in message for message in debugged)


def test_tuning_announces_itself_on_standard_output(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The tuning banner is a bare `print`, not a logger call."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-banner"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    recipe, _, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert "TUNING SOXSPIPE\n" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# THE `polyOrders` OVERRIDE AND THE PRODUCT ROW
# ---------------------------------------------------------------------------


def test_poly_orders_overrides_the_recipe_settings_degrees(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A four-digit `polyOrders` splits into the order and wavelength degree pairs."""
    # ARRANGE
    pinholePath = prepared_fits(tmp_path / "PINHOLE_VIS.fits", seed=12)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
        polyOrders=3456,
    )
    _, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeSettings["order-deg"] == [3, 4]
    assert recipe.recipeSettings["wavelength-deg"] == [5, 6]
    assert recipe.polyOrders == [3, 4, 5, 6]
    assert captured["mapArguments"][0]["recipeSettings"] is recipe.recipeSettings


def test_the_product_row_is_appended_to_the_tables_the_recipe_declares(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The row's columns and their order are the ones `base_recipe` declares."""
    # ARRANGE
    pinholePath = prepared_fits(tmp_path / "PINHOLE_VIS.fits", seed=12)
    productPath = tmp_path / "DISP_TAB_VIS.fits"
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
    )
    declaredQc, declaredProducts = _declared_tables()
    recipe.qc = declaredQc
    recipe.products = declaredProducts
    _patch_reduction(recipe, monkeypatch, productPath=productPath, fittedQc=_fitted_qc())

    # ACT
    recipe.produce_product()

    # ASSERT
    # `base_recipe` DOES NOT DECLARE `product_desc`, SO THE ROW'S OWN KEY ORDER
    # DECIDES WHERE THAT ONE COLUMN LANDS: AT THE END, AFTER `label`.
    row = recipe.products.iloc[-1]
    assert list(recipe.products.columns) == [*declaredProducts.columns, "product_desc"]
    assert row["soxspipe_recipe"] == "soxs-disp-solution"
    assert row["product_label"] == "DISP_MAP"
    assert row["file_name"] == productPath.name
    assert row["file_type"] == "FITS Table"
    assert row["obs_date_utc"] == "2024-01-02T03:04:05.678"
    assert row["product_desc"] == "VIS first pass dispersion solution"
    assert row["file_path"] == str(productPath)
    assert row["label"] == "PROD"
    assert recipe.dateObs == "2024-01-02T03:04:05.678"


def test_the_reduction_timestamp_renders_to_whole_seconds(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The product row's reduction timestamp holds the format both revisions render.

    This module writes exactly one product row, so there is no second row to
    share a value with. The contract that holds against the helper conversion
    is therefore the rendered format; the companion test added with the
    conversion pins that the clock is read exactly once.
    """
    # ARRANGE
    pinholePath = prepared_fits(tmp_path / "PINHOLE_VIS.fits", seed=12)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
    )
    _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert re.fullmatch(TIMESTAMP_PATTERN, recipe.products.iloc[-1]["reduction_date_utc"])


def test_the_reduction_reads_the_clock_exactly_once(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One product row means one clock read, which a format assertion cannot catch.

    This test holds only after the helper adoption, because it stubs the clock
    at the name the module imported. It wraps the real `utcnow_string` rather
    than replacing it, so the rendered format is still the production one.
    """
    # ARRANGE
    pinholePath = prepared_fits(tmp_path / "PINHOLE_VIS.fits", seed=12)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="VIS",
        instrument="SOXS",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        },
    )
    _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "DISP_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
    )
    realUtcnowString = DISPERSION_MODULE.utcnow_string
    clockReads: list[str] = []

    def counting_utcnow_string(**kwargs: object) -> str:
        rendered = realUtcnowString(**kwargs)
        clockReads.append(rendered)
        return rendered

    monkeypatch.setattr(DISPERSION_MODULE, "utcnow_string", counting_utcnow_string)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(clockReads) == 1
    assert recipe.products.iloc[-1]["reduction_date_utc"] == clockReads[0]
    assert re.fullmatch(TIMESTAMP_PATTERN, clockReads[0])
