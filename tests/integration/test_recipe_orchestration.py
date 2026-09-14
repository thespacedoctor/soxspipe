"""Recipe orchestration contracts with scientific collaborators isolated."""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_mbias, soxs_mdark, soxs_mflat
from tests.factories import (
    pipeline_settings,
    prepared_fits,
    product_table,
    qc_row,
    qc_table,
    synthetic_ccd,
)

pytestmark = pytest.mark.integration


@dataclass
class FakeFrameCollection:
    """Small image collection that records recipe queries."""

    ccdFrames: list[object] = field(default_factory=list)
    filePathsByToken: dict[str, list[str]] = field(default_factory=dict)
    calls: list[tuple[str, object]] = field(default_factory=list)
    activePaths: list[str] = field(default_factory=list)

    def ccds(self, **kwargs: object) -> list[object]:
        self.calls.append(("ccds", kwargs))
        return list(self.ccdFrames)

    def files_filtered(
        self,
        *,
        include_path: bool,
        **filters: object,
    ) -> list[str]:
        self.calls.append(("files_filtered", filters))
        if self.activePaths:
            return list(self.activePaths)
        for token, paths in self.filePathsByToken.items():
            if token in filters.values():
                return list(paths)
        return []

    def filter(self, **filters: object) -> FakeFrameCollection:
        self.calls.append(("filter", filters))
        selectedPaths = next(
            (
                paths
                for token, paths in self.filePathsByToken.items()
                if token in filters.values()
            ),
            [],
        )
        return FakeFrameCollection(
            ccdFrames=self.ccdFrames,
            filePathsByToken=self.filePathsByToken,
            calls=self.calls,
            activePaths=list(selectedPaths),
        )


def _empty_products() -> Any:
    return product_table().iloc[0:0].copy()


def _empty_qc() -> Any:
    return qc_table().iloc[0:0].copy()


def _assert_product_row(
    row: pd.Series,
    *,
    recipeName: str,
    productLabel: str,
    fileType: str,
    productPath: Path,
    description: str,
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
    assert row["file_name"] == productPath.name
    assert row["file_type"] == fileType
    assert row["obs_date_utc"] == "2024-01-02T03:04:05.678"
    assert isinstance(row["reduction_date_utc"], str)
    assert row["reduction_date_utc"]
    assert row["product_desc"] == description
    assert row["file_path"] == str(productPath)
    assert row["label"] == "PROD"


def test_master_bias_produce_product_runs_collaborators_and_records_product(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Master-bias orchestration returns its path and records its FITS product."""
    rawFrames = [synthetic_ccd(seed=1), synthetic_ccd(seed=2)]
    noiseFrame = synthetic_ccd(seed=3, prepared=True)
    stackedFrame = noiseFrame.copy()
    expectedMasterRon = float(stackedFrame.data.std())
    recipe = soxs_mbias.__new__(soxs_mbias)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.detectorParams = {}
    recipe.inputFrames = FakeFrameCollection(ccdFrames=rawFrames)
    recipe.recipeName = "soxs-mbias"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.workspaceRootPath = str(tmp_path)
    recipe.qc = qc_table()
    originalQc = recipe.qc.copy(deep=True)
    recipe.products = _empty_products()
    productPath = tmp_path / "MASTER_BIAS_VIS.fits"
    calls: list[str] = []
    qcArguments: dict[str, object] = {}
    expectedAddedQc = {
        "qc_periodic_pattern_noise": qc_row(
            recipeName="soxs-mbias",
            name="BIAS PATTERN",
            value=0.4,
            unit="electron",
            comment="Synthetic periodic pattern amplitude",
        ),
        "qc_ron": qc_row(
            recipeName="soxs-mbias",
            name="MASTER RON",
            value=1.2,
            unit="electron",
            comment="Synthetic master read noise",
        ),
        "qc_bias_structure": qc_row(
            recipeName="soxs-mbias",
            name="BIAS STRUCTURE",
            value=0.8,
            unit="electron",
            comment="Synthetic bias structure",
        ),
        "generic_quality_checks": qc_row(
            recipeName="soxs-mbias",
            name="BADPIX NUM",
            value=5,
            unit="",
            comment="Synthetic bad-pixel count",
        ),
        "qc_median_flux_level": qc_row(
            recipeName="soxs-mbias",
            name="MEDIAN FLUX",
            value=100.0,
            unit="electron",
            comment="Synthetic median bias level",
        ),
    }

    def fake_subtract(frame: object) -> tuple[float, float, object]:
        assert any(frame is rawFrame for rawFrame in rawFrames)
        return 100.0, 3.0, noiseFrame.copy()

    def fake_stack(**kwargs: object) -> object:
        assert len(kwargs["frames"]) == 2
        assert kwargs["recipe"] == "soxs_mbias"
        assert kwargs["post_stack_clipping"] is True
        calls.append("stack")
        return stackedFrame

    monkeypatch.setattr(recipe, "subtract_mean_flux_level", fake_subtract)
    monkeypatch.setattr(recipe, "clip_and_stack", fake_stack)

    def append_qc(methodName: str, *args: object, **kwargs: object) -> None:
        calls.append(methodName)
        qcArguments[methodName] = {"args": args, **kwargs}
        recipe.qc = pd.concat(
            [recipe.qc, expectedAddedQc[methodName]], ignore_index=True
        )

    monkeypatch.setattr(
        recipe,
        "qc_periodic_pattern_noise",
        lambda *args, **kwargs: append_qc("qc_periodic_pattern_noise", *args, **kwargs),
    )
    monkeypatch.setattr(
        recipe,
        "qc_ron",
        lambda *args, **kwargs: append_qc("qc_ron", *args, **kwargs),
    )
    monkeypatch.setattr(
        recipe,
        "qc_bias_structure",
        lambda *args, **kwargs: append_qc("qc_bias_structure", *args, **kwargs),
    )

    def fake_median_qc(*args: object, **kwargs: object) -> float:
        append_qc("qc_median_flux_level", *args, **kwargs)
        return 100.0

    monkeypatch.setattr(recipe, "qc_median_flux_level", fake_median_qc)

    def fake_update_fits_keywords(**kwargs: object) -> None:
        assert kwargs["frame"] is stackedFrame
        assert kwargs["frame"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        assert set(kwargs) == {"frame"}
        calls.append("update_fits_keywords")

    def fake_write(**kwargs: object) -> str:
        assert kwargs["frame"] is stackedFrame
        assert kwargs["filedir"] == recipe.workspaceRootPath
        assert kwargs["filename"] is False
        assert kwargs["overwrite"] is True
        calls.append("write")
        return str(productPath)

    monkeypatch.setattr(recipe, "update_fits_keywords", fake_update_fits_keywords)
    biasModule = import_module("soxspipe.recipes.soxs_mbias")

    def fake_generic_quality_checks(**kwargs: object) -> pd.DataFrame:
        assert kwargs == {
            "log": log,
            "frame": stackedFrame,
            "settings": recipe.settings,
            "recipeName": "soxs-mbias",
            "qcTable": recipe.qc,
        }
        calls.append("generic_quality_checks")
        return pd.concat(
            [kwargs["qcTable"], expectedAddedQc["generic_quality_checks"]],
            ignore_index=True,
        )

    monkeypatch.setattr(
        biasModule,
        "generic_quality_checks",
        fake_generic_quality_checks,
    )
    monkeypatch.setattr(
        recipe,
        "_write",
        fake_write,
    )
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc
    pd.testing.assert_frame_equal(
        returnedQc.reset_index(drop=True),
        pd.concat(
            [
                originalQc,
                expectedAddedQc["qc_periodic_pattern_noise"],
                expectedAddedQc["qc_ron"],
                expectedAddedQc["qc_bias_structure"],
                expectedAddedQc["generic_quality_checks"],
                expectedAddedQc["qc_median_flux_level"],
            ],
            ignore_index=True,
        ),
    )
    assert qcArguments["qc_periodic_pattern_noise"] == {
        "args": (),
        "frames": recipe.inputFrames,
    }
    ronArgs = qcArguments["qc_ron"]
    assert ronArgs["args"] == ()
    assert ronArgs["frameType"] == "MBIAS"
    assert ronArgs["frameName"] == "master bias"
    assert ronArgs["masterFrame"] is stackedFrame
    assert ronArgs["rawRon"] == pytest.approx(3.0)
    assert ronArgs["masterRon"] == pytest.approx(expectedMasterRon)
    assert qcArguments["qc_bias_structure"] == {
        "args": (stackedFrame,),
    }
    medianArgs = qcArguments["qc_median_flux_level"]
    assert medianArgs == {
        "args": (),
        "frame": stackedFrame,
        "frameType": "MBIAS",
        "frameName": "master bias",
        "medianFlux": pytest.approx(100.0),
    }
    _assert_product_row(
        recipe.products.iloc[-1],
        recipeName="soxs-mbias",
        productLabel="MBIAS",
        fileType="FITS",
        productPath=productPath,
        description="VIS Master bias frame",
    )
    assert calls == [
        "stack",
        "qc_periodic_pattern_noise",
        "qc_ron",
        "qc_bias_structure",
        "generic_quality_checks",
        "qc_median_flux_level",
        "update_fits_keywords",
        "write",
        "report",
        "clean_up",
    ]


def test_master_dark_produce_product_preserves_qc_and_records_product(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Master-dark orchestration returns the written path and product metadata."""
    rawFrames = [synthetic_ccd(seed=21), synthetic_ccd(seed=22)]
    noiseFrame = synthetic_ccd(seed=23, prepared=True)
    stackedFrame = noiseFrame.copy()
    recipe = soxs_mdark.__new__(soxs_mdark)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.detectorParams = {}
    recipe.inputFrames = FakeFrameCollection(ccdFrames=rawFrames)
    recipe.recipeName = "soxs-mdark"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.workspaceRootPath = str(tmp_path)
    recipe.qc = qc_table()
    originalQc = recipe.qc.copy(deep=True)
    recipe.products = _empty_products()
    productPath = tmp_path / "MASTER_DARK_VIS.fits"
    calls: list[str] = []
    genericQc = qc_row(
        recipeName="soxs-mdark",
        name="HOTPIX NUM",
        value=4,
        unit="",
        comment="Synthetic hot-pixel count",
    )
    medianQc = qc_row(
        recipeName="soxs-mdark",
        name="MEDIAN FLUX",
        value=95.0,
        unit="electron",
        comment="Synthetic median dark level",
    )
    ronQc = qc_row(
        recipeName="soxs-mdark",
        name="MASTER RON",
        value=1.0,
        unit="electron",
        comment="Synthetic master dark read noise",
    )

    def fake_subtract(frame: object) -> tuple[float, float, object]:
        assert any(frame is rawFrame for rawFrame in rawFrames)
        return 95.0, 3.0, noiseFrame.copy()

    def fake_stack(**kwargs: object) -> object:
        assert len(kwargs["frames"]) == 2
        assert kwargs["recipe"] == "soxs_mdark"
        assert kwargs["post_stack_clipping"] is True
        calls.append("stack")
        return stackedFrame

    monkeypatch.setattr(recipe, "subtract_mean_flux_level", fake_subtract)
    monkeypatch.setattr(recipe, "clip_and_stack", fake_stack)
    monkeypatch.setattr(
        recipe,
        "qc_median_flux_level",
        lambda **kwargs: (
            calls.append("median_qc"),
            setattr(
                recipe,
                "qc",
                pd.concat([recipe.qc, medianQc], ignore_index=True),
            ),
            95.0,
        )[-1],
    )
    monkeypatch.setattr(
        recipe,
        "qc_ron",
        lambda **kwargs: (
            calls.append("ron_qc"),
            setattr(
                recipe,
                "qc",
                pd.concat([recipe.qc, ronQc], ignore_index=True),
            ),
            (3.0, 1.0),
        )[-1],
    )

    def fake_update_fits_keywords(**kwargs: object) -> None:
        assert kwargs["frame"] is stackedFrame
        assert kwargs["frame"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        assert set(kwargs) == {"frame"}
        calls.append("keywords")

    def fake_write(**kwargs: object) -> str:
        assert kwargs["frame"] is stackedFrame
        assert kwargs["filedir"] == recipe.workspaceRootPath
        assert kwargs["filename"] is False
        assert kwargs["overwrite"] is True
        calls.append("write")
        return str(productPath)

    monkeypatch.setattr(recipe, "update_fits_keywords", fake_update_fits_keywords)
    darkModule = import_module("soxspipe.recipes.soxs_mdark")
    monkeypatch.setattr(
        darkModule,
        "generic_quality_checks",
        lambda **kwargs: (
            calls.append("generic_qc"),
            pd.concat([kwargs["qcTable"], genericQc], ignore_index=True),
        )[-1],
    )
    monkeypatch.setattr(
        recipe,
        "_write",
        fake_write,
    )
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc
    pd.testing.assert_frame_equal(
        returnedQc.reset_index(drop=True),
        pd.concat([originalQc, genericQc, medianQc, ronQc], ignore_index=True),
    )
    _assert_product_row(
        recipe.products.iloc[-1],
        recipeName="soxs-mdark",
        productLabel="MDARK",
        fileType="FITS",
        productPath=productPath,
        description="VIS Master dark frame",
    )
    assert calls == [
        "stack",
        "generic_qc",
        "median_qc",
        "ron_qc",
        "keywords",
        "write",
        "report",
        "clean_up",
    ]


def test_single_lamp_master_flat_records_stable_product_and_preserves_qc(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A single-lamp flat follows the public production workflow."""
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=31)
    flatFrame = synthetic_ccd(seed=32, prepared=True)
    firstStack = synthetic_ccd(seed=33, prepared=True)
    secondStack = synthetic_ccd(seed=34, prepared=True)
    maskedFrame = synthetic_ccd(seed=35, prepared=True)
    productPath = tmp_path / "MASTER_FLAT_VIS.fits"
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.inputFrames = FakeFrameCollection(
        filePathsByToken={"ORDER_TAB_VIS": [str(orderPath)]}
    )
    recipe.recipeName = "soxs-mflat"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.recipeSettings = {"subtract_background": False}
    recipe.sofName = "MASTER_FLAT_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.binRatioX = 1
    recipe.binRatioY = 1
    recipe.calibratedFlatFiles = ["flat-1.fits"]
    recipe.dFlatFiles = []
    recipe.qFlatFiles = []
    recipe.domeFlatFiles = []
    recipe.qc = qc_table()
    originalQc = recipe.qc.copy(deep=True)
    recipe.products = _empty_products()
    calls: list[str] = []
    normaliseCount = 0
    stackCount = 0
    keywordCount = 0
    writeCount = 0
    qcArguments: dict[str, dict[str, object]] = {}
    maskQc = qc_row(
        recipeName="soxs-mflat",
        name="LOW SENS PIX",
        value=7,
        unit="pixel",
        comment="Synthetic low-sensitivity pixel count",
    )
    genericQc = qc_row(
        recipeName="soxs-mflat",
        name="COLDPIX NUM",
        value=6,
        unit="",
        comment="Synthetic cold-pixel count",
    )
    spectroscopicQc = qc_row(
        recipeName="soxs-mflat",
        name="ORDER MEDIAN",
        value=99.5,
        unit="electron",
        comment="Synthetic order median flux",
    )

    monkeypatch.setattr(
        recipe,
        "calibrate_frame_set",
        lambda: calls.append("calibrate") or ([flatFrame], [], [], []),
    )

    def fake_normalise(*args: object, **kwargs: object) -> list[object]:
        nonlocal normaliseCount
        assert args[0] == [flatFrame]
        assert kwargs["orderTablePath"] == str(orderPath)
        assert kwargs["lamp"] == ""
        if normaliseCount == 0:
            assert "firstPassMasterFlat" not in kwargs
        else:
            assert kwargs["firstPassMasterFlat"] is firstStack
        normaliseCount += 1
        calls.append("normalise")
        return [flatFrame.copy()]

    def fake_stack(**kwargs: object) -> object:
        nonlocal stackCount
        assert len(kwargs["frames"]) == 1
        assert kwargs["recipe"] == "soxs_mflat"
        assert kwargs["ignore_input_masks"] is False
        assert kwargs["post_stack_clipping"] is False
        result = firstStack if stackCount == 0 else secondStack
        stackCount += 1
        calls.append("stack")
        return result

    def fake_update_fits_keywords(**kwargs: object) -> None:
        nonlocal keywordCount
        expectedFrame = secondStack if keywordCount == 0 else maskedFrame
        assert kwargs["frame"] is expectedFrame
        assert kwargs["frame"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        if keywordCount == 0:
            assert kwargs["rawFrames"] == ["flat-1.fits"]
        else:
            assert set(kwargs) == {"frame"}
        keywordCount += 1
        calls.append("keywords")

    monkeypatch.setattr(recipe, "normalise_flats", fake_normalise)
    monkeypatch.setattr(recipe, "clip_and_stack", fake_stack)
    monkeypatch.setattr(recipe, "update_fits_keywords", fake_update_fits_keywords)

    def fake_mask(**kwargs: object) -> tuple[object, pd.DataFrame]:
        assert kwargs["frame"] is secondStack
        assert kwargs["orderTablePath"] == str(orderPath)
        assert kwargs["returnMedianOrderFlux"] is True
        assert kwargs["writeQC"] is True
        calls.append("mask")
        recipe.qc = pd.concat([recipe.qc, maskQc], ignore_index=True)
        return maskedFrame, pd.DataFrame({"order": [10], "medianFlux": [100.0]})

    monkeypatch.setattr(recipe, "mask_low_sens_pixels", fake_mask)

    def fake_write(*args: object, **kwargs: object) -> str:
        nonlocal writeCount
        assert args[0].header["DATE-OBS"] == maskedFrame.header["DATE-OBS"]
        assert args[0].shape == maskedFrame.shape
        assert args[1] == recipe.settings["workspace-root-dir"]
        assert kwargs["overwrite"] is True
        if writeCount == 0:
            assert kwargs["filename"] == "MASTER_FLAT_VIS.fits"
        else:
            assert args[0] is maskedFrame
            assert "filename" not in kwargs
        writeCount += 1
        calls.append("write")
        return str(productPath)

    monkeypatch.setattr(recipe, "_write", fake_write)
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(recipe, "clean_up", lambda: calls.append("clean_up"))

    class FakeEdges:
        def __init__(self, **kwargs: object) -> None:
            assert kwargs["flatFrame"] is secondStack
            assert kwargs["orderCentreTable"] == str(orderPath)
            assert kwargs["settings"] is recipe.settings
            assert kwargs["qcTable"] is recipe.qc
            assert kwargs["productsTable"] is recipe.products
            calls.append("edges")

        def get(self) -> tuple[Any, Any, pd.DataFrame]:
            calls.append("edge_get")
            products = product_table().iloc[0:0].copy()
            products.loc[0] = {
                "soxspipe_recipe": "soxs-mflat",
                "product_label": "ORDER_LOC",
                "file_name": orderPath.name,
                "file_type": "FITS",
                "obs_date_utc": "2024-01-02T03:04:05.678",
                "reduction_date_utc": "2024-01-02T04:05:06.789",
                "product_desc": "Synthetic order locations",
                "file_path": str(orderPath),
                "label": "PROD",
            }
            return (
                products,
                recipe.qc.copy(deep=True),
                pd.DataFrame({"order": [10], "count": [32]}),
            )

    flatModule = import_module("soxspipe.recipes.soxs_mflat")
    monkeypatch.setattr(flatModule, "quicklook_image", lambda **kwargs: None)
    monkeypatch.setattr(flatModule, "detect_order_edges", FakeEdges)
    monkeypatch.setattr(
        flatModule,
        "unpack_order_table",
        lambda **kwargs: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )

    def fake_generic_qc(**kwargs: object) -> pd.DataFrame:
        assert kwargs["qcTable"] is recipe.qc
        calls.append("generic_qc")
        qcArguments["generic_qc"] = kwargs
        return pd.concat([kwargs["qcTable"], genericQc], ignore_index=True)

    def fake_spectroscopic_qc(**kwargs: object) -> pd.DataFrame:
        assert kwargs["qcTable"] is recipe.qc
        calls.append("spectroscopic_qc")
        qcArguments["spectroscopic_qc"] = kwargs
        return pd.concat([kwargs["qcTable"], spectroscopicQc], ignore_index=True)

    monkeypatch.setattr(flatModule, "generic_quality_checks", fake_generic_qc)
    monkeypatch.setattr(
        flatModule,
        "spectroscopic_image_quality_checks",
        fake_spectroscopic_qc,
    )

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc
    pd.testing.assert_frame_equal(
        recipe.qc.reset_index(drop=True),
        pd.concat(
            [originalQc, maskQc, originalQc, genericQc, spectroscopicQc],
            ignore_index=True,
        ),
    )
    genericArgs = qcArguments["generic_qc"]
    assert genericArgs["frame"] is maskedFrame
    assert genericArgs["settings"] is recipe.settings
    assert genericArgs["recipeName"] == "soxs-mflat"
    assert genericArgs["qcTable"] is not recipe.qc
    pd.testing.assert_frame_equal(
        genericArgs["qcTable"].reset_index(drop=True),
        pd.concat([originalQc, maskQc, originalQc], ignore_index=True),
    )
    spectroscopicArgs = qcArguments["spectroscopic_qc"]
    assert spectroscopicArgs["frame"] is maskedFrame
    assert spectroscopicArgs["settings"] is recipe.settings
    assert spectroscopicArgs["recipeName"] == "soxs-mflat"
    assert spectroscopicArgs["orderTablePath"] == str(orderPath)
    pd.testing.assert_frame_equal(
        spectroscopicArgs["qcTable"].reset_index(drop=True),
        pd.concat([originalQc, maskQc, originalQc, genericQc], ignore_index=True),
    )
    _assert_product_row(
        recipe.products.loc[recipe.products["product_label"] == "ORDER_LOC"].iloc[0],
        recipeName="soxs-mflat",
        productLabel="ORDER_LOC",
        fileType="FITS",
        productPath=orderPath,
        description="Synthetic order locations",
    )
    _assert_product_row(
        recipe.products.iloc[-1],
        recipeName="soxs-mflat",
        productLabel="MFLAT",
        fileType="FITS",
        productPath=productPath,
        description="VIS master spectroscopic flat frame",
    )
    assert calls == [
        "calibrate",
        "normalise",
        "stack",
        "normalise",
        "stack",
        "keywords",
        "edges",
        "edge_get",
        "mask",
        "write",
        "keywords",
        "write",
        "generic_qc",
        "spectroscopic_qc",
        "report",
        "clean_up",
    ]


def test_multi_lamp_master_flat_stitches_independent_lamp_products(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Three lamp sets produce tagged flats before the public UV stitch step."""
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=81)
    productPath = tmp_path / "MASTER_FLAT_VIS.fits"
    frames = [synthetic_ccd(seed=seed, prepared=True) for seed in range(82, 85)]
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.inputFrames = FakeFrameCollection(
        filePathsByToken={"ORDER_TAB_VIS": [str(orderPath)]}
    )
    recipe.recipeName = "soxs-mflat"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.recipeSettings = {"subtract_background": False}
    recipe.sofName = "MASTER_FLAT_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.binRatioX = 1
    recipe.binRatioY = 1
    recipe.calibratedFlatFiles = ["orderdef.fits"]
    recipe.dFlatFiles = ["dorderdef.fits"]
    recipe.qFlatFiles = ["qorderdef.fits"]
    recipe.domeFlatFiles = []
    recipe.qc = qc_table()
    recipe.products = _empty_products()
    calls: list[tuple[str, str]] = []
    normaliseCalls = 0

    monkeypatch.setattr(
        recipe,
        "calibrate_frame_set",
        lambda: ([frames[0]], [frames[1]], [frames[2]], []),
    )

    def fake_normalise(*args: object, **kwargs: object) -> list[object]:
        nonlocal normaliseCalls
        lamp = str(kwargs["lamp"])
        calls.append(("normalise", lamp))
        frameIndex = {"": 0, "_DLAMP": 1, "_QLAMP": 2}[lamp]
        assert args[0] == [frames[frameIndex]]
        assert kwargs["orderTablePath"] == str(orderPath)
        if normaliseCalls % 2:
            assert "firstPassMasterFlat" in kwargs
        else:
            assert "firstPassMasterFlat" not in kwargs
        normaliseCalls += 1
        return [frames[frameIndex].copy()]

    monkeypatch.setattr(recipe, "normalise_flats", fake_normalise)
    monkeypatch.setattr(
        recipe,
        "clip_and_stack",
        lambda **kwargs: kwargs["frames"][0].copy(),
    )
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **kwargs: None)
    monkeypatch.setattr(
        recipe,
        "mask_low_sens_pixels",
        lambda **kwargs: (
            kwargs["frame"].copy(),
            pd.DataFrame({"order": [10], "medianFlux": [100.0]}),
        ),
    )
    monkeypatch.setattr(recipe, "_write", lambda *args, **kwargs: str(productPath))
    monkeypatch.setattr(recipe, "report_output", lambda: recipe.qc)
    monkeypatch.setattr(recipe, "clean_up", lambda: None)

    def fake_stitch(medianFlux: pd.DataFrame, *, orderTablePath: str) -> object:
        assert orderTablePath == str(orderPath)
        assert set(medianFlux.columns) == {"order", "_DLAMP", "_QLAMP"}
        calls.append(("stitch", ""))
        return frames[0].copy()

    monkeypatch.setattr(recipe, "stitch_uv_mflats", fake_stitch)
    flatModule = import_module("soxspipe.recipes.soxs_mflat")
    monkeypatch.setattr(flatModule, "quicklook_image", lambda **kwargs: None)
    monkeypatch.setattr(
        flatModule,
        "unpack_order_table",
        lambda **kwargs: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(
        flatModule,
        "generic_quality_checks",
        lambda **kwargs: kwargs["qcTable"],
    )
    monkeypatch.setattr(
        flatModule,
        "spectroscopic_image_quality_checks",
        lambda **kwargs: kwargs["qcTable"],
    )

    class FakeEdges:
        def __init__(self, **kwargs: object) -> None:
            self.tag = str(kwargs["tag"])
            calls.append(("edges", self.tag))

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            label = f"ORDER_LOC{self.tag}"
            products = recipe.products.copy()
            products.loc[len(products)] = {
                "soxspipe_recipe": recipe.recipeName,
                "product_label": label,
                "file_name": orderPath.name,
                "file_type": "FITS",
                "obs_date_utc": "2024-01-02T03:04:05.678",
                "reduction_date_utc": "2024-01-02T04:05:06.789",
                "product_desc": "Synthetic order locations",
                "file_path": str(orderPath),
                "label": "PROD",
            }
            return products, recipe.qc.copy(), pd.DataFrame({"order": [10], "count": [32]})

    monkeypatch.setattr(flatModule, "detect_order_edges", FakeEdges)

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc
    assert calls == [
        ("normalise", ""),
        ("normalise", ""),
        ("edges", ""),
        ("normalise", "_DLAMP"),
        ("normalise", "_DLAMP"),
        ("edges", "_DLAMP"),
        ("normalise", "_QLAMP"),
        ("normalise", "_QLAMP"),
        ("edges", "_QLAMP"),
        ("stitch", ""),
    ]
    assert set(recipe.products["product_label"]) == {
        "ORDER_LOC",
        "ORDER_LOC_DLAMP",
        "ORDER_LOC_QLAMP",
        "MFLAT",
        "MFLAT_DLAMP",
        "MFLAT_QLAMP",
    }
