"""Mapping-recipe orchestration contracts with expensive fitting isolated."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import (
    soxs_disp_solution,
    soxs_order_centres,
    soxs_spatial_solution,
    soxs_straighten,
)
from tests.factories import (
    pipeline_settings,
    prepared_fits,
    product_table,
    qc_row,
    qc_table,
)

pytestmark = pytest.mark.integration


@dataclass
class RoutedFrameCollection:
    """Resolve image-collection queries against exact synthetic inventories."""

    filePathsByFilters: dict[frozenset[tuple[str, object]], tuple[str, ...]]
    activeFilters: dict[str, object] = field(default_factory=dict)
    calls: list[tuple[str, dict[str, object]]] = field(default_factory=list)

    def files_filtered(
        self,
        *,
        include_path: bool,
        **filters: object,
    ) -> list[str]:
        """Return paths matching the active and immediate filters."""
        combinedFilters = {**self.activeFilters, **filters}
        self.calls.append(("files_filtered", combinedFilters))
        return list(self.filePathsByFilters.get(frozenset(combinedFilters.items()), ()))

    def filter(self, **filters: object) -> RoutedFrameCollection:
        """Return a collection view carrying the requested filters."""
        self.calls.append(("filter", filters))
        return replace(self, activeFilters={**self.activeFilters, **filters})


def _empty_products() -> Any:
    return product_table().iloc[0:0].copy()


def _route(**filters: object) -> frozenset[tuple[str, object]]:
    return frozenset(filters.items())


def _assert_product_row(
    row: pd.Series,
    *,
    recipeName: str,
    productLabel: str,
    fileName: str,
    fileType: str,
    description: str,
    filePath: str | None,
) -> None:
    assert set(row.index) == {
        "soxspipe_recipe",
        "product_label",
        "file_name",
        "file_type",
        "obs_date_utc",
        "reduction_date_utc",
        "product_desc",
        "file_path",
        "label",
    }
    assert row["soxspipe_recipe"] == recipeName
    assert row["product_label"] == productLabel
    assert row["file_name"] == fileName
    assert row["file_type"] == fileType
    assert row["obs_date_utc"] == "2024-01-02T03:04:05.678"
    assert isinstance(row["reduction_date_utc"], str)
    assert row["reduction_date_utc"]
    assert row["product_desc"] == description
    assert row["file_path"] == filePath
    assert row["label"] == "PROD"


def test_dispersion_solution_produce_product_calibrates_maps_and_reports(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Dispersion orchestration detrends a real FITS frame before fitting."""
    biasPath = prepared_fits(tmp_path / "MASTER_BIAS_VIS.fits", seed=11)
    pinholePath = prepared_fits(tmp_path / "PINHOLE_VIS.fits", seed=12)
    productPath = tmp_path / "DISP_TAB_VIS.fits"
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RoutedFrameCollection(
        filePathsByFilters={
            _route(PRO_CATG="MASTER_BIAS_VIS"): (str(biasPath),),
            _route(DPR_TYPE="LAMP,WAVE", DPR_TECH="ECHELLE,PINHOLE"): (
                str(pinholePath),
            ),
        }
    )
    recipe.recipeName = "soxs-disp-solution"
    recipe.settings = pipeline_settings(tmp_path, overrides={"tune-pipeline": False})
    recipe.recipeSettings = {}
    recipe.workspaceRootPath = str(tmp_path)
    recipe.sofName = "synthetic-dispersion"
    recipe.startNightDate = "2024-01-02"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.polyOrders = False
    recipe.qc = qc_table()
    originalQc = recipe.qc.copy(deep=True)
    fittedQc = qc_row(
        recipeName="soxs-disp-solution",
        name="DISPERSION_RESIDUAL",
        value=0.012,
        unit="pixel",
        comment="Synthetic dispersion-fit residual",
    )
    expectedQc = pd.concat([originalQc, fittedQc])
    recipe.products = _empty_products()
    calls: list[str] = []
    calibratedFrames: list[object] = []

    def fake_detrend(**kwargs: object) -> object:
        assert kwargs["master_bias"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        assert kwargs["dark"] is False
        calibratedFrame = kwargs["inputFrame"].copy()
        calibratedFrames.append(calibratedFrame)
        calls.append("detrend")
        return calibratedFrame

    def fake_update_fits_keywords(**kwargs: object) -> None:
        assert kwargs == {"frame": calibratedFrames[0]}
        assert kwargs["frame"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        calls.append("keywords")

    monkeypatch.setattr(recipe, "detrend", fake_detrend)
    monkeypatch.setattr(recipe, "update_fits_keywords", fake_update_fits_keywords)

    class FakeDispersionMap:
        def __init__(self, **kwargs: object) -> None:
            assert kwargs["pinholeFrame"] is calibratedFrames[0]
            assert kwargs["settings"] is recipe.settings
            assert kwargs["qcTable"] is recipe.qc
            assert kwargs["productsTable"] is recipe.products
            calls.append("create_dispersion_map")

        def get(self) -> tuple[object, ...]:
            calls.append("map_get")
            return str(productPath), None, (), fittedQc, _empty_products(), None

    dispersionModule = import_module("soxspipe.recipes.soxs_disp_solution")
    monkeypatch.setattr(dispersionModule, "create_dispersion_map", FakeDispersionMap)
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc
    pd.testing.assert_frame_equal(returnedQc, expectedQc)
    _assert_product_row(
        recipe.products.iloc[-1],
        recipeName="soxs-disp-solution",
        productLabel="DISP_MAP",
        fileName=productPath.name,
        fileType="FITS Table",
        description="VIS first pass dispersion solution",
        filePath=str(productPath),
    )
    assert calls == [
        "detrend",
        "keywords",
        "create_dispersion_map",
        "map_get",
        "report",
        "clean_up",
    ]


def test_order_centres_produce_product_preserves_qc_and_records_order_table(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Order-centre fitting returns its table and records the product contract."""
    masterBiasPath = prepared_fits(tmp_path / "MASTER_BIAS_VIS.fits", seed=31)
    orderFramePath = prepared_fits(tmp_path / "ORDER_DEFINITION_VIS.fits", seed=32)
    dispersionPath = tmp_path / "DISP_TAB_VIS.fits"
    dispersionPath.touch()
    productPath = tmp_path / "ORDER_TAB_VIS.fits"
    inputFrames = RoutedFrameCollection(
        filePathsByFilters={
            _route(PRO_CATG="MASTER_BIAS_VIS"): (str(masterBiasPath),),
            _route(DPR_TYPE="FLAT,LAMP", DPR_TECH="ECHELLE,PINHOLE"): (
                str(orderFramePath),
            ),
            _route(PRO_CATG="DISP_TAB_VIS"): (str(dispersionPath),),
        }
    )
    recipe = soxs_order_centres.__new__(soxs_order_centres)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.detectorParams = {}
    recipe.inputFrames = inputFrames
    recipe.recipeName = "soxs-order-centres"
    recipe.settings = pipeline_settings(
        tmp_path,
        overrides={"tune-pipeline": False},
    )
    recipe.workspaceRootPath = str(tmp_path)
    recipe.recipeSettings = {"detect-continuum": {}}
    recipe.polyOrders = False
    recipe.sofName = "synthetic-order-centres"
    recipe.startNightDate = "2024-01-02"
    recipe.qc = qc_table()
    originalQc = recipe.qc.copy(deep=True)
    recipe.products = _empty_products()
    calls: list[str] = []
    calibratedFrames: list[object] = []
    expectedFittedQc = qc_row(
        recipeName="soxs-order-centres",
        name="ORDER_CENTRE_RESIDUAL",
        value=0.021,
        unit="pixel",
        comment="Synthetic order-centre residual",
    )
    expectedQc = pd.concat([originalQc, expectedFittedQc])
    fittedQc = expectedFittedQc
    fittedProducts = _empty_products()

    def fake_detrend(**kwargs: object) -> object:
        assert kwargs["master_bias"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        assert kwargs["dark"] is False
        calibratedFrame = kwargs["inputFrame"].copy()
        calibratedFrames.append(calibratedFrame)
        calls.append("detrend")
        return calibratedFrame

    def fake_update_fits_keywords(**kwargs: object) -> None:
        assert kwargs == {"frame": calibratedFrames[0]}
        assert kwargs["frame"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        calls.append("keywords")

    monkeypatch.setattr(recipe, "detrend", fake_detrend)
    monkeypatch.setattr(recipe, "update_fits_keywords", fake_update_fits_keywords)

    class FakeContinuumDetector:
        def __init__(self, **kwargs: object) -> None:
            assert kwargs["traceFrame"] is calibratedFrames[0]
            assert kwargs["dispersion_map"] == str(dispersionPath)
            assert kwargs["settings"] is recipe.settings
            assert kwargs["qcTable"] is recipe.qc
            assert kwargs["productsTable"] is recipe.products
            calls.append("detect_continuum")

        def get(self) -> tuple[object, ...]:
            return (
                str(productPath),
                fittedQc,
                fittedProducts,
                object(),
                object(),
                object(),
            )

    orderModule = import_module("soxspipe.recipes.soxs_order_centres")
    monkeypatch.setattr(orderModule, "detect_continuum", FakeContinuumDetector)
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc
    pd.testing.assert_frame_equal(returnedQc, expectedQc)
    _assert_product_row(
        recipe.products.iloc[-1],
        recipeName="soxs-order-centres",
        productLabel="ORDER_CENTRES",
        fileName=productPath.name,
        fileType="FITS",
        description="VIS order centre traces",
        filePath=str(productPath),
    )
    assert calls == ["detrend", "keywords", "detect_continuum", "report", "clean_up"]


def test_spatial_solution_produce_product_preserves_qc_and_records_maps(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Spatial fitting returns the 2D map and records both mapping products."""
    masterBiasPath = prepared_fits(tmp_path / "MASTER_BIAS_VIS.fits", seed=41)
    masterFlatPath = prepared_fits(tmp_path / "MASTER_FLAT_VIS.fits", seed=42)
    pinholePath = prepared_fits(tmp_path / "MULTI_PINHOLE_VIS.fits", seed=43)
    orderTablePath = tmp_path / "ORDER_TAB_VIS.fits"
    dispersionPath = tmp_path / "DISP_TAB_VIS.fits"
    orderTablePath.touch()
    dispersionPath.touch()
    mapPath = tmp_path / "DISP_IMAGE_VIS.fits"
    mapImagePath = tmp_path / "2D_MAP_VIS.fits"
    inputFrames = RoutedFrameCollection(
        filePathsByFilters={
            _route(PRO_CATG="MASTER_BIAS_VIS"): (str(masterBiasPath),),
            _route(PRO_CATG="MASTER_FLAT_VIS"): (str(masterFlatPath),),
            _route(DPR_TYPE="WAVE,LAMP", DPR_TECH="ECHELLE,MULTI-PINHOLE"): (
                str(pinholePath),
            ),
            _route(PRO_CATG="ORDER_TAB_VIS"): (str(orderTablePath),),
            _route(PRO_CATG="DISP_TAB_VIS"): (str(dispersionPath),),
        }
    )
    recipe = soxs_spatial_solution.__new__(soxs_spatial_solution)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.detectorParams = {}
    recipe.inputFrames = inputFrames
    recipe.recipeName = "soxs-spat-solution"
    recipe.settings = pipeline_settings(
        tmp_path,
        overrides={"tune-pipeline": False},
    )
    recipe.workspaceRootPath = str(tmp_path)
    recipe.recipeSettings = {"use_flat": True}
    recipe.create2DMap = True
    recipe.polyOrders = False
    recipe.sofName = "synthetic-spatial-solution"
    recipe.startNightDate = "2024-01-02"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.qc = qc_table()
    originalQc = recipe.qc.copy(deep=True)
    recipe.products = _empty_products()
    calls: list[str] = []
    calibratedFrames: list[object] = []
    expectedFittedQc = qc_row(
        recipeName="soxs-spat-solution",
        name="SPATIAL_RESIDUAL",
        value=0.034,
        unit="pixel",
        comment="Synthetic spatial-fit residual",
    )
    expectedQc = pd.concat([originalQc, expectedFittedQc])
    fittedQc = expectedFittedQc
    fittedProducts = _empty_products()

    def fake_detrend(**kwargs: object) -> object:
        assert kwargs["master_bias"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        assert kwargs["dark"] is False
        assert kwargs["master_flat"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        assert kwargs["order_table"] == str(orderTablePath)
        calibratedFrame = kwargs["inputFrame"].copy()
        calibratedFrames.append(calibratedFrame)
        calls.append("detrend")
        return calibratedFrame

    def fake_update_fits_keywords(**kwargs: object) -> None:
        assert kwargs == {"frame": calibratedFrames[0]}
        assert kwargs["frame"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        calls.append("keywords")

    monkeypatch.setattr(recipe, "detrend", fake_detrend)
    monkeypatch.setattr(recipe, "update_fits_keywords", fake_update_fits_keywords)

    class FakeDispersionMap:
        def __init__(self, **kwargs: object) -> None:
            assert kwargs["pinholeFrame"] is calibratedFrames[0]
            assert kwargs["firstGuessMap"] == str(dispersionPath)
            assert kwargs["orderTable"] == str(orderTablePath)
            assert kwargs["settings"] is recipe.settings
            assert kwargs["qcTable"] is recipe.qc
            assert kwargs["productsTable"] is recipe.products
            assert kwargs["create2DMap"] is True
            assert kwargs["arcFrame"] is None
            calls.append("create_dispersion_map")

        def get(self) -> tuple[object, ...]:
            return (
                str(mapPath),
                str(mapImagePath),
                (),
                fittedQc,
                fittedProducts,
                object(),
            )

    commonutilsModule = import_module("soxspipe.commonutils")
    toolkitModule = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(commonutilsModule, "create_dispersion_map", FakeDispersionMap)
    monkeypatch.setattr(
        toolkitModule,
        "quicklook_image",
        lambda **kwargs: calls.append("quicklook"),
    )
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(mapImagePath)
    assert returnedQc is recipe.qc
    pd.testing.assert_frame_equal(returnedQc, expectedQc)
    outputProducts = recipe.products.tail(2)
    _assert_product_row(
        outputProducts.iloc[0],
        recipeName="soxs-spat-solution",
        productLabel="SPAT_SOL",
        fileName=mapPath.name,
        fileType="FITS",
        description="VIS full dispersion-spatial solution",
        filePath=None,
    )
    _assert_product_row(
        outputProducts.iloc[1],
        recipeName="soxs-spat-solution",
        productLabel="2D_MAP",
        fileName=mapImagePath.name,
        fileType="FITS",
        description="VIS 2D detector map of wavelength, slit position and order",
        filePath=None,
    )
    assert calls == [
        "detrend",
        "keywords",
        "create_dispersion_map",
        "quicklook",
        "report",
        "clean_up",
    ]


def test_straighten_produce_product_reports_then_cleans_up(log: Any) -> None:
    """Straighten retains its public placeholder result and cleanup order."""
    recipe = soxs_straighten.__new__(soxs_straighten)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = lambda keyword: keyword
    recipe.detectorParams = {}
    calls: list[str] = []
    expectedQc = qc_table()
    recipe.report_output = lambda: calls.append("report") or expectedQc
    recipe.clean_up = lambda: calls.append("clean_up")

    productPath, returnedQc = recipe.produce_product()

    assert productPath is None
    assert returnedQc is expectedQc
    pd.testing.assert_frame_equal(returnedQc, expectedQc)
    assert calls == ["report", "clean_up"]
