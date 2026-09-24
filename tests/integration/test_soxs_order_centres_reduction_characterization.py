"""Characterization of the `soxs_order_centres.produce_product` branches the orchestration test misses.

`tests/integration/test_mapping_recipe_orchestration.py` already covers one
path: a SOXS VIS reduction with a master bias, no dark, no binning keywords,
`tune-pipeline` off and no `polyOrders`. That leaves lines 253, 267-269, 291,
318-320, 324-325, 331-384 and 387-390 uncovered at the branch point, the
largest single block being the 54-line `tune-pipeline` branch.

These tests pin that region before anything restructures it: the master-dark
read, the NIR lamp-off dark, the X-Shooter filter variants, the last-frame-wins
precedence of the frame loops, the intermediate-product write, the binning
read, the `polyOrders` override and the product row. The tuning branch and the
failure points are in `test_soxs_order_centres_tuning_characterization.py`,
which imports the builders defined here. Both files reuse the orchestration
module's routing collection so they build the recipe the same way.
"""

from __future__ import annotations

import re
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_order_centres
from soxspipe.recipes.base_recipe import base_recipe
from tests.factories import pipeline_settings, prepared_fits, qc_row, qc_table
from tests.integration.test_mapping_recipe_orchestration import (
    RoutedFrameCollection,
    _empty_products,
    _route,
)

pytestmark = pytest.mark.integration


# THE RENDERED SHAPE OF THE REDUCTION TIMESTAMP, WHOLE SECONDS AND NO OFFSET.
TIMESTAMP_PATTERN = r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}"

ORDER_MODULE = import_module("soxspipe.recipes.soxs_order_centres")

# THE KEYWORD ARGUMENTS EVERY CONTINUUM DETECTOR IS BUILT WITH, BESIDES THE
# FRAME, THE TABLE AND THE BINNING.
DETECTOR_KEYWORDS = {
    "recipeName",
    "settings",
    "recipeSettings",
    "qcTable",
    "productsTable",
    "sofName",
    "startNightDate",
    "log",
    "traceFrame",
    "dispersion_map",
    "binx",
    "biny",
}


def _declared_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the empty QC and product tables a real recipe starts from.

    The integration factories declare a different column order from
    `base_recipe`, so a test that pins column order has to start from the
    recipe's own declaration.
    """
    return base_recipe._empty_qc_and_product_tables(base_recipe, pd)


def _configure_recipe(
    recipe: soxs_order_centres,
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
    recipe.recipeName = "soxs-order-centres"
    recipe.settings = pipeline_settings(
        tmpPath,
        overrides={"tune-pipeline": False, **(settingOverrides or {})},
    )
    recipe.recipeSettings = {"detect-continuum": {}}
    recipe.workspaceRootPath = str(tmpPath)
    recipe.sofName = "synthetic-order-centres"
    recipe.startNightDate = "2024-01-02"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.polyOrders = polyOrders
    recipe.qc = qc_table()
    recipe.products = _empty_products()


def _vis_recipe(
    log: Any,
    tmpPath: Path,
    **options: Any,
) -> soxs_order_centres:
    """Return a SOXS VIS recipe routed to one pinhole frame and one dispersion table."""
    pinholePath = prepared_fits(tmpPath / "PINHOLE_VIS.fits", seed=12)
    dispersionPath = tmpPath / "DISP_TAB_VIS.fits"
    dispersionPath.touch()
    routes = {
        _route(DPR_TYPE="FLAT,LAMP", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
        _route(PRO_CATG="DISP_TAB_VIS"): (str(dispersionPath),),
        **options.pop("extraRoutes", {}),
    }
    recipe = soxs_order_centres.__new__(soxs_order_centres)
    _configure_recipe(recipe, log, tmpPath, arm="VIS", instrument="SOXS", routes=routes, **options)
    return recipe


def _patch_reduction(
    recipe: soxs_order_centres,
    monkeypatch: pytest.MonkeyPatch,
    *,
    productPath: Path | None,
    fittedQc: pd.DataFrame,
    fittedProducts: pd.DataFrame | None = None,
    header: dict[str, object] | None = None,
) -> tuple[list[str], dict[str, Any]]:
    """Stub the collaborators outside the recipe and return the call log and captures."""
    calls: list[str] = []
    captured: dict[str, Any] = {"detectors": []}

    def fake_detrend(**kwargs: Any) -> Any:
        captured["detrend"] = kwargs
        calls.append("detrend")
        calibrated = kwargs["inputFrame"].copy()
        for key, value in (header or {}).items():
            calibrated.header[key] = value
        return calibrated

    monkeypatch.setattr(recipe, "detrend", fake_detrend)
    monkeypatch.setattr(
        recipe,
        "update_fits_keywords",
        lambda **kwargs: calls.append("keywords"),
    )

    class FakeContinuumDetector:
        def __init__(self, **kwargs: object) -> None:
            captured["detectors"].append(kwargs)
            calls.append("detect_continuum")

        def sample_trace(self) -> tuple[object, float]:
            calls.append("sample_trace")
            return "ORDER-PIXEL-TABLE", 87.5

        def get(self) -> tuple[object, ...]:
            calls.append("detector_get")
            return (
                None if productPath is None else str(productPath),
                fittedQc,
                _empty_products() if fittedProducts is None else fittedProducts,
                "POLY-TABLE",
                "PIXEL-TABLE",
                "META-TABLE",
            )

    monkeypatch.setattr(ORDER_MODULE, "detect_continuum", FakeContinuumDetector)
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))
    return calls, captured


def _fitted_qc() -> pd.DataFrame:
    """Return the QC row the stubbed continuum detector hands back."""
    return qc_row(
        recipeName="soxs-order-centres",
        name="ORDER_CENTRE_RESIDUAL",
        value=0.021,
        unit="pixel",
        comment="Synthetic order-centre residual",
    )


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
# THE CALIBRATION AND ORDER-DEFINITION FRAME READS
# ---------------------------------------------------------------------------


def test_a_master_dark_is_read_and_reaches_the_detrend(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A `MASTER_DARK` for the arm is read and handed to `detrend` with the bias."""
    # ARRANGE
    biasPath = prepared_fits(tmp_path / "MASTER_BIAS_VIS.fits", seed=11)
    darkPath = prepared_fits(tmp_path / "MASTER_DARK_VIS.fits", seed=13)
    recipe = _vis_recipe(
        log,
        tmp_path,
        extraRoutes={
            _route(PRO_CATG="MASTER_BIAS_VIS"): (str(biasPath),),
            _route(PRO_CATG="MASTER_DARK_VIS"): (str(darkPath),),
        },
    )
    calls, captured = _patch_reduction(
        recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc()
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["master_bias"].data[0][0] == pytest.approx(_read_first_pixel(biasPath), rel=1e-12)
    assert captured["detrend"]["dark"].data[0][0] == pytest.approx(_read_first_pixel(darkPath), rel=1e-12)
    assert calls[0] == "detrend"


def test_the_soxs_lamp_off_frame_overrides_a_master_dark(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The lamp-off `IMAGE` read runs after the master-dark read, so it wins."""
    # ARRANGE
    darkPath = prepared_fits(tmp_path / "MASTER_DARK_NIR.fits", seed=13)
    lampOffPath = prepared_fits(tmp_path / "LAMP_OFF_NIR.fits", seed=14)
    pinholePath = prepared_fits(tmp_path / "PINHOLE_NIR.fits", seed=12)
    recipe = soxs_order_centres.__new__(soxs_order_centres)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="NIR",
        instrument="SOXS",
        routes={
            _route(PRO_CATG="MASTER_DARK_NIR"): (str(darkPath),),
            _route(DPR_TYPE="FLAT,LAMP", DPR_TECH="IMAGE"): (str(lampOffPath),),
            _route(DPR_TYPE="FLAT,LAMP", DPR_TECH="ECHELLE,PINHOLE"): (str(pinholePath),),
            _route(PRO_CATG="DISP_TAB_NIR"): (str(tmp_path / "DISP_TAB_NIR.fits"),),
        },
    )
    _, captured = _patch_reduction(
        recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_NIR.fits", fittedQc=_fitted_qc()
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["master_bias"] is False
    assert captured["detrend"]["dark"].data[0][0] == pytest.approx(_read_first_pixel(lampOffPath), rel=1e-12)
    assert captured["detrend"]["inputFrame"].data[0][0] == pytest.approx(_read_first_pixel(pinholePath), rel=1e-12)


def test_the_last_matching_soxs_order_frame_is_the_one_reduced(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The SOXS filter list is walked in order, so its last match wins."""
    # ARRANGE
    firstPath = prepared_fits(tmp_path / "PINHOLE_A_VIS.fits", seed=21)
    lastPath = prepared_fits(tmp_path / "PINHOLE_B_VIS.fits", seed=22)
    recipe = _vis_recipe(
        log,
        tmp_path,
        extraRoutes={
            _route(DPR_TYPE="LAMP,DFLAT", DPR_TECH="ECHELLE,PINHOLE"): (str(firstPath),),
            _route(DPR_TYPE="LAMP,DFLAT", DPR_TECH="ECHELLE,SLIT"): (str(lastPath),),
        },
    )
    _, captured = _patch_reduction(
        recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc()
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["inputFrame"].data[0][0] == pytest.approx(_read_first_pixel(lastPath), rel=1e-12)
    routedFilters = [filters for name, filters in recipe.inputFrames.calls if name == "files_filtered"]
    orderFilters = [filters for filters in routedFilters if "DPR_TECH" in filters and filters["DPR_TECH"] != "IMAGE"]
    assert orderFilters == [
        {"DPR_TYPE": "FLAT,LAMP", "DPR_TECH": "ECHELLE,PINHOLE"},
        {"DPR_TYPE": "LAMP,FLAT", "DPR_TECH": "ECHELLE,PINHOLE"},
        {"DPR_TYPE": "LAMP,DFLAT", "DPR_TECH": "ECHELLE,PINHOLE"},
        {"DPR_TYPE": "FLAT,LAMP", "DPR_TECH": "ECHELLE,SLIT"},
        {"DPR_TYPE": "LAMP,FLAT", "DPR_TECH": "ECHELLE,SLIT"},
        {"DPR_TYPE": "LAMP,DFLAT", "DPR_TECH": "ECHELLE,SLIT"},
    ]


def test_an_xshooter_reduction_routes_the_orderdef_filters(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A non-SOXS instrument reads `LAMP,ORDERDEF` for the dark and walks the three ORDERDEF lamps."""
    # ARRANGE
    lampOffPath = prepared_fits(tmp_path / "XSH_LAMP_OFF.fits", seed=31)
    qthPath = prepared_fits(tmp_path / "XSH_QTH.fits", seed=32)
    d2Path = prepared_fits(tmp_path / "XSH_D2.fits", seed=33)
    recipe = soxs_order_centres.__new__(soxs_order_centres)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="UVB",
        instrument="XSH",
        routes={
            _route(DPR_TYPE="LAMP,ORDERDEF", DPR_TECH="IMAGE"): (str(lampOffPath),),
            _route(DPR_TYPE="LAMP,QORDERDEF", DPR_TECH="ECHELLE,PINHOLE"): (str(qthPath),),
            _route(DPR_TYPE="LAMP,DORDERDEF", DPR_TECH="ECHELLE,PINHOLE"): (str(d2Path),),
            _route(PRO_CATG="DISP_TAB_UVB"): (str(tmp_path / "DISP_TAB_UVB.fits"),),
        },
    )
    _, captured = _patch_reduction(
        recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_UVB.fits", fittedQc=_fitted_qc()
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detrend"]["inputFrame"].data[0][0] == pytest.approx(_read_first_pixel(d2Path), rel=1e-12)
    assert captured["detrend"]["dark"].data[0][0] == pytest.approx(_read_first_pixel(lampOffPath), rel=1e-12)
    routedFilters = [filters for name, filters in recipe.inputFrames.calls if name == "files_filtered"]
    assert routedFilters == [
        {"PRO_CATG": "MASTER_BIAS_UVB"},
        {"PRO_CATG": "MASTER_DARK_UVB"},
        {"DPR_TYPE": "LAMP,ORDERDEF", "DPR_TECH": "IMAGE"},
        {"DPR_TYPE": "LAMP,ORDERDEF", "DPR_TECH": "ECHELLE,PINHOLE"},
        {"DPR_TYPE": "LAMP,QORDERDEF", "DPR_TECH": "ECHELLE,PINHOLE"},
        {"DPR_TYPE": "LAMP,DORDERDEF", "DPR_TECH": "ECHELLE,PINHOLE"},
        {"PRO_CATG": "DISP_TAB_UVB"},
    ]


def test_the_last_dispersion_table_listed_is_the_one_traced(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The dispersion-table loop keeps its last match, as a path, not a read table."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="DISP_TAB_VIS")] = ("first.fits", "last.fits")
    _, captured = _patch_reduction(
        recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc()
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["detectors"][0]["dispersion_map"] == "last.fits"


# ---------------------------------------------------------------------------
# THE INTERMEDIATE PRODUCT WRITE AND THE BINNING READ
# ---------------------------------------------------------------------------


def test_saving_intermediate_products_writes_the_calibrated_order_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With the setting on, the calibrated frame is written and its path printed."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, settingOverrides={"save-intermediate-products": True})
    calls, _ = _patch_reduction(recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc())
    writeArguments: dict[str, Any] = {}

    def fake_write(*args: object, **kwargs: object) -> str:
        writeArguments["args"] = args
        writeArguments.update(kwargs)
        calls.append("write")
        return str(tmp_path / "PINHOLE_VIS_pre.fits")

    monkeypatch.setattr(recipe, "_write", fake_write)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert writeArguments["args"] == (recipe.orderFrame, str(tmp_path))
    assert writeArguments["filename"] is False
    assert writeArguments["overwrite"] is True
    assert writeArguments["product"] is False
    assert calls[:4] == ["detrend", "keywords", "write", "detect_continuum"]
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == [f"\nCalibrated single pinhole frame frame saved to {tmp_path / 'PINHOLE_VIS_pre.fits'}\n"]


@pytest.mark.parametrize(
    ("arm", "header", "expected"),
    [
        ("VIS", {"WIN_BINX": "2", "WIN_BINY": 3}, (2, 3)),
        ("UVB", {"WIN_BINX": 1.0, "WIN_BINY": 2.0}, (1, 2)),
        ("VIS", {}, (1, 1)),
        ("NIR", {"WIN_BINX": 2, "WIN_BINY": 2}, (1, 1)),
    ],
)
def test_the_binning_comes_from_the_calibrated_frame_header_off_the_nir_arm(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    arm: str,
    header: dict[str, object],
    expected: tuple[int, int],
) -> None:
    """`WIN_BINX` and `WIN_BINY` are read as integers, except on NIR or when absent."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    recipe.arm = arm
    recipe.inputFrames.filePathsByFilters[_route(PRO_CATG=f"DISP_TAB_{arm}")] = ("DISP_TAB.fits",)
    _, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "ORDER_TAB.fits",
        fittedQc=_fitted_qc(),
        header=header,
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    detector = captured["detectors"][0]
    assert (detector["binx"], detector["biny"]) == expected
    assert type(detector["binx"]) is int
    assert type(detector["biny"]) is int


# ---------------------------------------------------------------------------
# THE `polyOrders` OVERRIDE, THE TABLES AND THE PRODUCT ROW
# ---------------------------------------------------------------------------


def test_poly_orders_overrides_the_continuum_degrees(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A two-digit `polyOrders` splits into the order and dispersion-axis degrees."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, polyOrders=35)
    _, captured = _patch_reduction(
        recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc()
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeSettings == {"detect-continuum": {"order-deg": 3, "disp-axis-deg": 5}}
    assert recipe.polyOrders == [3, 5]
    assert captured["detectors"][0]["recipeSettings"] is recipe.recipeSettings


def test_a_long_poly_orders_uses_only_its_first_two_digits(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Extra digits are split into the list and then ignored."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, polyOrders=3456)
    _patch_reduction(recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc())

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeSettings == {"detect-continuum": {"order-deg": 3, "disp-axis-deg": 4}}
    assert recipe.polyOrders == [3, 4, 5, 6]


def test_the_detector_tables_are_appended_before_the_product_row(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The detector's QC and product rows land first, then the recipe's own row."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    detectorProducts = pd.DataFrame([{"product_label": "DETECTOR_ROW", "file_name": "detector.fits"}])
    originalQc = recipe.qc.copy(deep=True)
    calls, _ = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmp_path / "ORDER_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
        fittedProducts=detectorProducts,
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert list(recipe.products["product_label"]) == ["DETECTOR_ROW", "ORDER_CENTRES"]
    assert list(recipe.products.index) == [0, 1]
    pd.testing.assert_frame_equal(recipe.qc, pd.concat([originalQc, _fitted_qc()]))
    assert calls == ["detrend", "keywords", "detect_continuum", "detector_get", "report", "clean_up"]


def test_the_product_row_is_appended_to_the_tables_the_recipe_declares(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The row's columns and their order are the ones `base_recipe` declares."""
    # ARRANGE
    productPath = tmp_path / "ORDER_TAB_VIS.fits"
    recipe = _vis_recipe(log, tmp_path)
    declaredQc, declaredProducts = _declared_tables()
    recipe.qc = declaredQc
    recipe.products = declaredProducts
    _patch_reduction(recipe, monkeypatch, productPath=productPath, fittedQc=_fitted_qc())

    # ACT
    returnedPath, returnedQc = recipe.produce_product()

    # ASSERT
    # `base_recipe` DOES NOT DECLARE `product_desc`, SO THE ROW'S OWN KEY ORDER
    # DECIDES WHERE THAT ONE COLUMN LANDS: AT THE END, AFTER `label`.
    row = recipe.products.iloc[-1]
    assert list(recipe.products.columns) == [*declaredProducts.columns, "product_desc"]
    assert row["soxspipe_recipe"] == "soxs-order-centres"
    assert row["product_label"] == "ORDER_CENTRES"
    assert row["file_name"] == productPath.name
    assert row["file_type"] == "FITS"
    assert row["obs_date_utc"] == "2024-01-02T03:04:05.678"
    assert row["product_desc"] == "VIS order centre traces"
    assert row["file_path"] == str(productPath)
    assert row["label"] == "PROD"
    assert recipe.dateObs == "2024-01-02T03:04:05.678"
    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc


def test_a_full_reduction_logs_only_its_own_entry_and_exit(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With every collaborator stubbed, the only debug events are the method's own pair."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    _patch_reduction(recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc())

    # ACT
    recipe.produce_product()

    # ASSERT
    debugged = [message for level, message in log.messages if level == "debug"]
    assert debugged == [
        "starting the ``produce_product`` method",
        "completed the ``produce_product`` method",
    ]


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
    recipe = _vis_recipe(log, tmp_path)
    _patch_reduction(recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc())

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
    recipe = _vis_recipe(log, tmp_path)
    _patch_reduction(recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc())
    realUtcnowString = ORDER_MODULE.utcnow_string
    clockReads: list[str] = []

    def counting_utcnow_string(**kwargs: object) -> str:
        rendered = realUtcnowString(**kwargs)
        clockReads.append(rendered)
        return rendered

    monkeypatch.setattr(ORDER_MODULE, "utcnow_string", counting_utcnow_string)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(clockReads) == 1
    assert recipe.products.iloc[-1]["reduction_date_utc"] == clockReads[0]
    assert re.fullmatch(TIMESTAMP_PATTERN, clockReads[0])
