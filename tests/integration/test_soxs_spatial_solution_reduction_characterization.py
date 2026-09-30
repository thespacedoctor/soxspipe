"""Characterization of the `soxs_spatial_solution.produce_product` branches the orchestration test misses.

`tests/integration/test_mapping_recipe_orchestration.py` already covers one
path: a SOXS VIS reduction with a master bias and a master flat, no dark, no
slit arc, `use_flat` on, `debug` off, `tune-pipeline` off and no `polyOrders`.

These tests pin the rest of `produce_product` before anything restructures
it: the dark reads, the `use_flat` switch, the multi-pinhole and slit-arc
filters for both instruments, last-match and first-match precedence, the
intermediate-product write, the `debug` switch, the `polyOrders` override, the
two product rows and the shared timestamp, and the quick-look call. The tuning
branch and the failure points are in
`test_soxs_spatial_solution_tuning_characterization.py`, which imports the
builders defined here. Both files reuse the orchestration module's routing
collection so they build the recipe the same way.
"""

from __future__ import annotations

import re
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_spatial_solution
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

# THE RECIPE IMPORTS BOTH COLLABORATORS INSIDE `produce_product`, SO THEY ARE
# REPLACED ON THE MODULES THEY ARE IMPORTED FROM.
COMMONUTILS_MODULE = import_module("soxspipe.commonutils")
TOOLKIT_MODULE = import_module("soxspipe.commonutils.toolkit")

# THE DATE-OBS EVERY PREPARED FACTORY FRAME CARRIES.
FACTORY_DATE_OBS = "2024-01-02T03:04:05.678"


def _declared_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the empty QC and product tables a real recipe starts from."""
    return base_recipe._empty_qc_and_product_tables(base_recipe, pd)


def _configure_recipe(
    recipe: soxs_spatial_solution,
    log: Any,
    tmpPath: Path,
    *,
    arm: str,
    instrument: str,
    routes: dict[frozenset[tuple[str, object]], tuple[str, ...]],
    settingOverrides: dict[str, Any] | None = None,
    polyOrders: Any = False,
    useFlat: bool = True,
    debug: bool = False,
    create2DMap: bool = True,
) -> None:
    """Populate an unconstructed recipe with everything `produce_product` reads."""
    recipe.log = log
    recipe.arm = arm
    recipe.inst = instrument
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RoutedFrameCollection(filePathsByFilters=routes)
    recipe.recipeName = "soxs-spat-solution"
    recipe.settings = pipeline_settings(
        tmpPath,
        overrides={"tune-pipeline": False, **(settingOverrides or {})},
    )
    recipe.recipeSettings = {"use_flat": useFlat}
    recipe.workspaceRootPath = str(tmpPath)
    recipe.create2DMap = create2DMap
    recipe.polyOrders = polyOrders
    recipe.sofName = "synthetic-spatial-solution"
    recipe.startNightDate = "2024-01-02"
    recipe.debug = debug
    recipe.turnOffMP = True
    recipe.qc = qc_table()
    recipe.products = _empty_products()


def _vis_recipe(log: Any, tmpPath: Path, **options: Any) -> soxs_spatial_solution:
    """Return a SOXS VIS recipe routed to one multi-pinhole frame and both tables."""
    pinholePath = prepared_fits(tmpPath / "MULTI_PINHOLE_VIS.fits", seed=43)
    routes = {
        _route(DPR_TYPE="WAVE,LAMP", DPR_TECH="ECHELLE,MULTI-PINHOLE"): (str(pinholePath),),
        _route(PRO_CATG="ORDER_TAB_VIS"): (str(tmpPath / "ORDER_TAB_VIS.fits"),),
        _route(PRO_CATG="DISP_TAB_VIS"): (str(tmpPath / "DISP_TAB_VIS.fits"),),
        **options.pop("extraRoutes", {}),
    }
    recipe = soxs_spatial_solution.__new__(soxs_spatial_solution)
    _configure_recipe(recipe, log, tmpPath, arm="VIS", instrument="SOXS", routes=routes, **options)
    return recipe


def _patch_reduction(
    recipe: soxs_spatial_solution,
    monkeypatch: pytest.MonkeyPatch,
    *,
    mapPath: Path,
    mapImagePath: Path | None,
    header: dict[str, object] | None = None,
) -> tuple[list[str], dict[str, Any]]:
    """Stub the collaborators outside the recipe and return the call log and captures."""
    calls: list[str] = []
    captured: dict[str, Any] = {"detrend": [], "maps": []}

    def fake_detrend(**kwargs: Any) -> Any:
        captured["detrend"].append(kwargs)
        calls.append("detrend")
        calibrated = kwargs["inputFrame"].copy()
        for key, value in (header or {}).items():
            calibrated.header[key] = value
        return calibrated

    monkeypatch.setattr(recipe, "detrend", fake_detrend)
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **kwargs: calls.append("keywords"))

    class FakeDispersionMap:
        def __init__(self, **kwargs: object) -> None:
            captured["maps"].append(kwargs)
            calls.append("create_dispersion_map")

        def get(self) -> tuple[object, ...]:
            calls.append("map_get")
            return (
                str(mapPath),
                None if mapImagePath is None else str(mapImagePath),
                "RESIDUAL-PLOTS",
                _fitted_qc(),
                pd.DataFrame([{"product_label": "MAP_ROW", "file_name": "map-row.fits"}]),
                "LINE-DETECTION-TABLE",
            )

    def fake_quicklook(**kwargs: object) -> None:
        captured["quicklook"] = kwargs
        calls.append("quicklook")

    monkeypatch.setattr(COMMONUTILS_MODULE, "create_dispersion_map", FakeDispersionMap)
    monkeypatch.setattr(TOOLKIT_MODULE, "quicklook_image", fake_quicklook)
    monkeypatch.setattr(recipe, "report_output", lambda: calls.append("report") or recipe.qc)
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))
    return calls, captured


def _fitted_qc() -> pd.DataFrame:
    """Return the QC row the stubbed dispersion map hands back."""
    return qc_row(
        recipeName="soxs-spat-solution",
        name="SPATIAL_RESIDUAL",
        value=0.034,
        unit="pixel",
        comment="Synthetic spatial-fit residual",
    )


def _first_pixel(path: Path) -> float:
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


def _routed_filters(recipe: soxs_spatial_solution) -> list[dict[str, object]]:
    """Return every filter set `files_filtered` was asked for, in order."""
    return [filters for name, filters in recipe.inputFrames.calls if name == "files_filtered"]


# ---------------------------------------------------------------------------
# THE CALIBRATION FRAME READS
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("useFlat", [True, False])
def test_the_calibrations_reach_both_detrends_and_use_flat_gates_the_flat(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    useFlat: bool,
) -> None:
    """The bias, dark and flat are read once and handed to both the pinhole and the arc detrend."""
    # ARRANGE
    biasPath = prepared_fits(tmp_path / "MASTER_BIAS_VIS.fits", seed=11)
    darkPath = prepared_fits(tmp_path / "MASTER_DARK_VIS.fits", seed=12)
    flatPath = prepared_fits(tmp_path / "MASTER_FLAT_VIS.fits", seed=13)
    arcPath = prepared_fits(tmp_path / "SLIT_ARC_VIS.fits", seed=14)
    recipe = _vis_recipe(
        log,
        tmp_path,
        useFlat=useFlat,
        extraRoutes={
            _route(PRO_CATG="MASTER_BIAS_VIS"): (str(biasPath),),
            _route(PRO_CATG="MASTER_DARK_VIS"): (str(darkPath),),
            _route(PRO_CATG="MASTER_FLAT_VIS"): (str(flatPath),),
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,SLIT"): (str(arcPath),),
        },
    )
    calls, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls[:3] == ["detrend", "detrend", "keywords"]
    pinhole, arc = captured["detrend"]
    assert pinhole["inputFrame"].data[0][0] == pytest.approx(_first_pixel(tmp_path / "MULTI_PINHOLE_VIS.fits"))
    assert arc["inputFrame"].data[0][0] == pytest.approx(_first_pixel(arcPath))
    for detrended in (pinhole, arc):
        assert detrended["master_bias"].data[0][0] == pytest.approx(_first_pixel(biasPath))
        assert detrended["dark"].data[0][0] == pytest.approx(_first_pixel(darkPath))
        assert detrended["order_table"] == str(tmp_path / "ORDER_TAB_VIS.fits")
        if useFlat:
            assert detrended["master_flat"].data[0][0] == pytest.approx(_first_pixel(flatPath))
        else:
            assert detrended["master_flat"] is False


@pytest.mark.parametrize(("instrument", "lampOffType"), [("SOXS", "WAVE,LAMP"), ("XSH", "LAMP,WAVE")])
def test_the_lamp_off_frame_overrides_a_master_dark_per_instrument(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    instrument: str,
    lampOffType: str,
) -> None:
    """The lamp-off `IMAGE` read runs after the master-dark read, so it wins."""
    # ARRANGE
    darkPath = prepared_fits(tmp_path / "MASTER_DARK_NIR.fits", seed=21)
    lampOffPath = prepared_fits(tmp_path / "LAMP_OFF_NIR.fits", seed=22)
    pinholePath = prepared_fits(tmp_path / "MULTI_PINHOLE_NIR.fits", seed=23)
    recipe = soxs_spatial_solution.__new__(soxs_spatial_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="NIR",
        instrument=instrument,
        routes={
            _route(PRO_CATG="MASTER_DARK_NIR"): (str(darkPath),),
            _route(DPR_TYPE=lampOffType, DPR_TECH="IMAGE"): (str(lampOffPath),),
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,MULTI-PINHOLE"): (str(pinholePath),),
            _route(PRO_CATG="ORDER_TAB_NIR"): ("ORDER_TAB_NIR.fits",),
            _route(PRO_CATG="DISP_TAB_NIR"): ("DISP_TAB_NIR.fits",),
        },
    )
    _, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    (pinhole,) = captured["detrend"]
    assert pinhole["master_bias"] is False
    assert pinhole["dark"].data[0][0] == pytest.approx(_first_pixel(lampOffPath))
    assert pinhole["inputFrame"].data[0][0] == pytest.approx(_first_pixel(pinholePath))


def test_a_soxs_reduction_walks_both_arc_spellings_and_keeps_the_last_match(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """SOXS tries `WAVE,LAMP` then `LAMP,WAVE` for the pinhole and the slit arc; the last read wins."""
    # ARRANGE
    lastPinholePath = prepared_fits(tmp_path / "PINHOLE_LAST.fits", seed=31)
    firstArcPath = prepared_fits(tmp_path / "ARC_FIRST.fits", seed=32)
    lastArcPath = prepared_fits(tmp_path / "ARC_LAST.fits", seed=33)
    recipe = _vis_recipe(
        log,
        tmp_path,
        extraRoutes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,MULTI-PINHOLE"): (str(lastPinholePath),),
            _route(DPR_TYPE="WAVE,LAMP", DPR_TECH="ECHELLE,SLIT"): (str(firstArcPath),),
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,SLIT"): (str(lastArcPath),),
        },
    )
    _, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    pinhole, arc = captured["detrend"]
    assert pinhole["inputFrame"].data[0][0] == pytest.approx(_first_pixel(lastPinholePath))
    assert arc["inputFrame"].data[0][0] == pytest.approx(_first_pixel(lastArcPath))
    assert _routed_filters(recipe) == [
        {"PRO_CATG": "MASTER_BIAS_VIS"},
        {"PRO_CATG": "MASTER_DARK_VIS"},
        {"DPR_TYPE": "WAVE,LAMP", "DPR_TECH": "IMAGE"},
        {"PRO_CATG": "MASTER_FLAT_VIS"},
        {"DPR_TYPE": "WAVE,LAMP", "DPR_TECH": "ECHELLE,MULTI-PINHOLE"},
        {"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "ECHELLE,MULTI-PINHOLE"},
        {"DPR_TYPE": "WAVE,LAMP", "DPR_TECH": "ECHELLE,SLIT"},
        {"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "ECHELLE,SLIT"},
        {"PRO_CATG": "ORDER_TAB_VIS"},
        {"PRO_CATG": "DISP_TAB_VIS"},
    ]


def test_an_xshooter_reduction_reads_only_the_lamp_wave_spelling(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """X-Shooter uses one filter each for the lamp-off, the pinhole and the slit arc."""
    # ARRANGE
    pinholePath = prepared_fits(tmp_path / "XSH_PINHOLE.fits", seed=41)
    recipe = soxs_spatial_solution.__new__(soxs_spatial_solution)
    _configure_recipe(
        recipe,
        log,
        tmp_path,
        arm="UVB",
        instrument="XSH",
        routes={
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,MULTI-PINHOLE"): (str(pinholePath),),
            _route(PRO_CATG="ORDER_TAB_UVB"): ("ORDER_TAB_UVB.fits",),
            _route(PRO_CATG="DISP_TAB_UVB"): ("DISP_TAB_UVB.fits",),
        },
    )
    _, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _routed_filters(recipe) == [
        {"PRO_CATG": "MASTER_BIAS_UVB"},
        {"PRO_CATG": "MASTER_DARK_UVB"},
        {"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "IMAGE"},
        {"PRO_CATG": "MASTER_FLAT_UVB"},
        {"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "ECHELLE,MULTI-PINHOLE"},
        {"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "ECHELLE,SLIT"},
        {"PRO_CATG": "ORDER_TAB_UVB"},
        {"PRO_CATG": "DISP_TAB_UVB"},
    ]
    assert captured["maps"][0]["arcFrame"] is None
    assert recipe.slit_arc is None


def test_the_first_order_table_and_the_last_dispersion_table_are_used(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The order table is indexed `[0]` through `filter`; the dispersion table loop keeps its last match."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="ORDER_TAB_VIS")] = ("order-first.fits", "order-last.fits")
    recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="DISP_TAB_VIS")] = ("disp-first.fits", "disp-last.fits")
    _, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["maps"][0]["orderTable"] == "order-first.fits"
    assert captured["maps"][0]["firstGuessMap"] == "disp-last.fits"
    assert captured["detrend"][0]["order_table"] == "order-first.fits"
    assert ("filter", {"PRO_CATG": "ORDER_TAB_VIS"}) in recipe.inputFrames.calls


# ---------------------------------------------------------------------------
# THE SLIT ARC, THE DATE AND THE INTERMEDIATE WRITE
# ---------------------------------------------------------------------------


def test_a_slit_arc_is_detrended_and_handed_to_the_map_as_the_arc_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The detrended arc replaces the raw one on `self.slit_arc` and reaches the dispersion map."""
    # ARRANGE
    arcPath = prepared_fits(tmp_path / "SLIT_ARC_VIS.fits", seed=51)
    recipe = _vis_recipe(
        log,
        tmp_path,
        extraRoutes={_route(DPR_TYPE="WAVE,LAMP", DPR_TECH="ECHELLE,SLIT"): (str(arcPath),)},
    )
    _, captured = _patch_reduction(
        recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None, header={"CALIBRATED": True}
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.slit_arc.header["CALIBRATED"] is True
    assert captured["maps"][0]["arcFrame"] is recipe.slit_arc
    assert captured["maps"][0]["pinholeFrame"] is recipe.multiPinholeFrame


def test_the_observation_date_comes_from_the_raw_pinhole_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`self.dateObs` is read before the detrend, so a calibrated header does not change it."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    _patch_reduction(
        recipe,
        monkeypatch,
        mapPath=tmp_path / "MAP.fits",
        mapImagePath=None,
        header={"DATE-OBS": "2030-01-01T00:00:00.000"},
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.dateObs == FACTORY_DATE_OBS
    assert recipe.products.iloc[-1]["obs_date_utc"] == FACTORY_DATE_OBS


def test_saving_intermediate_products_writes_the_calibrated_pinhole_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With the setting on, the calibrated frame is written after the keywords and its path printed."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, settingOverrides={"save-intermediate-products": True})
    calls, _ = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)
    writeArguments: dict[str, Any] = {}

    def fake_write(*args: object, **kwargs: object) -> str:
        writeArguments["args"] = args
        writeArguments.update(kwargs)
        calls.append("write")
        return str(tmp_path / "PINHOLE_pre.fits")

    monkeypatch.setattr(recipe, "_write", fake_write)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert writeArguments["args"] == (recipe.multiPinholeFrame, str(tmp_path))
    assert (writeArguments["filename"], writeArguments["overwrite"], writeArguments["product"]) == (False, True, False)
    assert calls[:4] == ["detrend", "keywords", "write", "create_dispersion_map"]
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == [f"\nCalibrated multi pinhole frame frame saved to {tmp_path / 'PINHOLE_pre.fits'}\n"]


# ---------------------------------------------------------------------------
# THE MAP CALL: `debug`, `polyOrders` AND THE ARGUMENTS
# ---------------------------------------------------------------------------


def test_the_dispersion_map_receives_the_reduction_state(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every keyword the full map is built with is pinned."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    _, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    (arguments,) = captured["maps"]
    assert set(arguments) == {
        "log",
        "settings",
        "recipeSettings",
        "pinholeFrame",
        "firstGuessMap",
        "orderTable",
        "qcTable",
        "productsTable",
        "sofName",
        "create2DMap",
        "startNightDate",
        "arcFrame",
        "debug",
        "turnOffMP",
    }
    assert arguments["create2DMap"] is True
    assert arguments["debug"] is False
    assert arguments["turnOffMP"] is True
    assert arguments["startNightDate"] == "2024-01-02"
    assert arguments["sofName"] == "synthetic-spatial-solution"
    assert arguments["recipeSettings"] is recipe.recipeSettings


def test_debug_turns_off_the_2d_map_and_drops_even_a_detrended_slit_arc(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """In debug mode the arc is still detrended, then replaced by False before the map is built."""
    # ARRANGE
    arcPath = prepared_fits(tmp_path / "SLIT_ARC_VIS.fits", seed=61)
    recipe = _vis_recipe(
        log,
        tmp_path,
        debug=True,
        extraRoutes={_route(DPR_TYPE="WAVE,LAMP", DPR_TECH="ECHELLE,SLIT"): (str(arcPath),)},
    )
    calls, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert calls.count("detrend") == 2
    assert recipe.create2DMap is False
    assert recipe.slit_arc is False
    assert captured["maps"][0]["create2DMap"] is False
    assert captured["maps"][0]["arcFrame"] is False
    assert captured["maps"][0]["debug"] is True


def test_poly_orders_overrides_all_three_degree_pairs(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A six-digit `polyOrders` splits into order, wavelength and slit degree pairs."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, polyOrders=345435)
    _, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeSettings == {
        "use_flat": True,
        "order-deg": [3, 4],
        "wavelength-deg": [5, 4],
        "slit-deg": [3, 5],
    }
    assert recipe.polyOrders == [3, 4, 5, 4, 3, 5]
    assert captured["maps"][0]["recipeSettings"] is recipe.recipeSettings


@pytest.mark.parametrize(
    ("polyOrders", "slitDegrees"),
    [(3454, []), (34543, [3]), (3454356, [3, 5, 6])],
)
def test_the_slit_pair_takes_every_digit_after_the_fourth(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    polyOrders: int,
    slitDegrees: list[int],
) -> None:
    """The constructor never checks the digit count, so the slit "pair" is whatever is left."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, polyOrders=polyOrders)
    _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.recipeSettings["slit-deg"] == slitDegrees
    assert recipe.recipeSettings["wavelength-deg"] == [5, 4]


# ---------------------------------------------------------------------------
# THE PRODUCT ROWS, THE QUICK-LOOK AND THE RETURN VALUE
# ---------------------------------------------------------------------------


def test_both_product_rows_land_in_the_declared_tables_with_their_own_file_path(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The map's own rows land first, then SPAT_SOL and 2D_MAP, each recording the file it wrote."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    declaredQc, declaredProducts = _declared_tables()
    recipe.qc = declaredQc
    recipe.products = declaredProducts
    mapPath = tmp_path / "SPAT_SOL_VIS.fits"
    mapImagePath = tmp_path / "2D_MAP_VIS.fits"
    _patch_reduction(recipe, monkeypatch, mapPath=mapPath, mapImagePath=mapImagePath)

    # ACT
    recipe.produce_product()

    # ASSERT
    # `base_recipe` DOES NOT DECLARE `product_desc`, SO THE ROWS' OWN KEY ORDER
    # DECIDES WHERE THAT ONE COLUMN LANDS: AT THE END, AFTER `label`.
    assert list(recipe.products.columns) == [*declaredProducts.columns, "product_desc"]
    assert list(recipe.products["product_label"]) == ["MAP_ROW", "SPAT_SOL", "2D_MAP"]
    assert list(recipe.products.index) == [0, 1, 2]
    spatial, image = recipe.products.iloc[1], recipe.products.iloc[2]
    assert (spatial["file_name"], image["file_name"]) == (mapPath.name, mapImagePath.name)
    assert spatial["product_desc"] == "VIS full dispersion-spatial solution"
    assert image["product_desc"] == "VIS 2D detector map of wavelength, slit position and order"
    assert (spatial["file_path"], image["file_path"]) == (str(mapPath), str(mapImagePath))
    for row in (spatial, image):
        assert row["soxspipe_recipe"] == "soxs-spat-solution"
        assert row["file_type"] == "FITS"
        assert row["label"] == "PROD"
        assert row["obs_date_utc"] == FACTORY_DATE_OBS


def test_both_product_rows_share_one_reduction_timestamp(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The two rows are stamped from one value, in whole seconds."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=tmp_path / "IMAGE.fits")

    # ACT
    recipe.produce_product()

    # ASSERT
    stamps = list(recipe.products.iloc[-2:]["reduction_date_utc"])
    assert stamps[0] == stamps[1]
    assert re.fullmatch(TIMESTAMP_PATTERN, stamps[0])


def test_without_a_2d_map_image_only_the_spatial_row_is_recorded_and_none_is_returned(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The method returns the map *image* path, so a run with no image returns None."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    calls, _ = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    returnedPath, returnedQc = recipe.produce_product()

    # ASSERT
    assert returnedPath is None
    assert returnedQc is recipe.qc
    assert list(recipe.products["product_label"]) == ["MAP_ROW", "SPAT_SOL"]
    assert calls[-5:] == ["create_dispersion_map", "map_get", "quicklook", "report", "clean_up"]


def test_the_quicklook_overlays_both_maps_on_the_calibrated_pinhole_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The quick-look call is pinned argument by argument, and the image path is returned."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    mapPath, mapImagePath = tmp_path / "MAP.fits", tmp_path / "IMAGE.fits"
    _, captured = _patch_reduction(recipe, monkeypatch, mapPath=mapPath, mapImagePath=mapImagePath)

    # ACT
    returnedPath, _ = recipe.produce_product()

    # ASSERT
    assert returnedPath == str(mapImagePath)
    assert captured["quicklook"] == {
        "log": recipe.log,
        "CCDObject": recipe.multiPinholeFrame,
        "show": False,
        "ext": False,
        "stdWindow": 1,
        "title": "Multi-pinhole Frame Overlaid with Dispersion Solution",
        "surfacePlot": True,
        "dispMap": str(mapPath),
        "dispMapImage": str(mapImagePath),
        "settings": recipe.settings,
        "skylines": False,
    }


def test_the_map_tables_are_merged_into_the_quality_control_table(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The map's QC rows are appended to the recipe's QC table before the report."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    originalQc = recipe.qc.copy(deep=True)
    _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    pd.testing.assert_frame_equal(recipe.qc, pd.concat([originalQc, _fitted_qc()]))


def test_a_full_reduction_logs_only_its_own_entry_and_exit(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With every collaborator stubbed, the only debug events are the method's own pair."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT
    recipe.produce_product()

    # ASSERT
    debugged = [message for level, message in log.messages if level == "debug"]
    assert debugged == [
        "starting the ``produce_product`` method",
        "completed the ``produce_product`` method",
    ]


def test_the_reduction_reads_the_clock_exactly_once_for_both_rows(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two product rows, one clock read: the helper must receive the value, not mint its own.

    This test holds only after the helper adoption, because it stubs the clock
    at the name the module imported. It wraps the real `utcnow_string` rather
    than replacing it, so the rendered format is still the production one.
    """
    # ARRANGE
    spatialModule = import_module("soxspipe.recipes.soxs_spatial_solution")
    recipe = _vis_recipe(log, tmp_path)
    _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=tmp_path / "IMAGE.fits")
    realUtcnowString = spatialModule.utcnow_string
    clockReads: list[str] = []

    def counting_utcnow_string(**kwargs: object) -> str:
        rendered = realUtcnowString(**kwargs)
        clockReads.append(rendered)
        return rendered

    monkeypatch.setattr(spatialModule, "utcnow_string", counting_utcnow_string)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(clockReads) == 1
    assert list(recipe.products.iloc[-2:]["reduction_date_utc"]) == [clockReads[0], clockReads[0]]
    assert re.fullmatch(TIMESTAMP_PATTERN, clockReads[0])
