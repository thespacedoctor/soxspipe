"""Observing-recipe orchestration contracts with science isolated."""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_nod, soxs_offset, soxs_stare
from tests.factories import (
    pipeline_settings,
    product_table,
    qc_row,
    qc_table,
    synthetic_ccd,
)

pytestmark = pytest.mark.integration


@dataclass
class RouteCollection:
    """Route image-collection filters to deterministic temporary files."""

    routes: dict[frozenset[tuple[str, object]], list[str]]
    activePaths: list[str] | None = None
    calls: list[tuple[str, dict[str, object]]] = field(default_factory=list)

    def files_filtered(
        self,
        *,
        include_path: bool,
        **filters: object,
    ) -> list[str]:
        self.calls.append(("files_filtered", filters))
        if self.activePaths is not None:
            return list(self.activePaths)
        return list(self.routes.get(frozenset(filters.items()), []))

    def filter(self, **filters: object) -> RouteCollection:
        self.calls.append(("filter", filters))
        return RouteCollection(
            routes=self.routes,
            activePaths=list(self.routes.get(frozenset(filters.items()), [])),
            calls=self.calls,
        )


def _route(**filters: object) -> frozenset[tuple[str, object]]:
    return frozenset(filters.items())


def _write_prepared_frame(
    destination: Path,
    *,
    seed: int,
    headerOverrides: dict[str, object],
) -> Path:
    frame = synthetic_ccd(
        seed=seed,
        prepared=True,
        headerOverrides=headerOverrides,
    )
    hdus = frame.to_hdu(hdu_mask="QUAL", hdu_uncertainty="ERRS")
    hdus[0].name = "FLUX"
    hdus.writeto(destination)
    return destination


def _empty_products() -> pd.DataFrame:
    return product_table().iloc[0:0].copy()


def _keyword(keyword: str) -> str:
    return {
        "DATE_OBS": "DATE-OBS",
        "NOD_CUMULATIVE_OFFSETY": "ESO SEQ CUMOFF Y",
        "OFFSET_RA": "HIERARCH ESO SEQ FIXOFF RA",
        "OFFSET_DEC": "HIERARCH ESO SEQ FIXOFF DEC",
    }.get(keyword, keyword)


def _spectrum_product(products: pd.DataFrame, path: Path) -> pd.DataFrame:
    row = {
        "soxspipe_recipe": "observing-recipe",
        "product_label": "MERGED_SPECTRUM",
        "file_name": path.name,
        "file_type": "FITS",
        "obs_date_utc": "2024-01-02T03:04:05.678",
        "reduction_date_utc": "2024-01-02T04:05:06.789",
        "product_desc": "Synthetic merged spectrum",
        "file_path": str(path),
        "label": "PROD",
    }
    return pd.concat([products, pd.DataFrame([row])], ignore_index=True)


def _configure_nodding_recipe(
    recipe: Any,
    *,
    log: Any,
    tmp_path: Path,
    objectPaths: list[Path],
    technique: str,
) -> None:
    routes = {
        _route(DPR_TYPE="OBJECT", DPR_TECH=technique): [
            str(path) for path in objectPaths
        ],
        _route(PRO_CATG="ORDER_TAB_VIS"): [str(tmp_path / "ORDER_TAB_VIS.fits")],
        _route(PRO_CATG="DISP_TAB_VIS"): [str(tmp_path / "DISP_TAB_VIS.fits")],
        _route(PRO_CATG="DISP_IMAGE_VIS"): [str(tmp_path / "DISP_IMAGE_VIS.fits")],
    }
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = _keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RouteCollection(routes)
    recipe.recipeName = "soxs-offset" if "OFFSET" in technique else "soxs-nod"
    recipe.recipeSettings = {"use_flat": False}
    recipe.settings = pipeline_settings(tmp_path)
    recipe.generateReponseCurve = False
    recipe.sofName = "OBJECT_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.productDir = str(tmp_path)
    recipe.qcDir = str(tmp_path / "qc")
    recipe.filenameTemplate = "OBJECT_VIS.fits"
    recipe.dateObs = "2024-01-02T03:04:05.678"
    recipe.debug = False
    recipe.qc = qc_table()
    recipe.products = _empty_products()


def _patch_nodding_collaborators(
    recipe: Any,
    *,
    monkeypatch: pytest.MonkeyPatch,
    plotPath: Path,
    extractionPath: Path,
) -> tuple[list[str], dict[str, list[dict[str, object]]], pd.DataFrame]:
    calls: list[str] = []
    captured: dict[str, list[dict[str, object]]] = {
        "stack": [],
        "keywords": [],
        "extract_cycle": [],
        "stack_extractions": [],
        "plot": [],
        "clean_up": [],
    }
    spectrum = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})
    extractionQc = qc_row(
        recipeName=recipe.recipeName,
        name="EXTRACTED ORDER COUNT",
        value=2,
        unit="count",
        comment="Synthetic extracted-order count",
    )

    def fake_stack(**kwargs: object) -> Any:
        calls.append("stack")
        captured["stack"].append(kwargs)
        return kwargs["frames"][0].copy()

    monkeypatch.setattr(
        recipe,
        "clip_and_stack",
        fake_stack,
    )

    def fake_keywords(**kwargs: object) -> None:
        calls.append("keywords")
        captured["keywords"].append(kwargs)

    monkeypatch.setattr(
        recipe,
        "update_fits_keywords",
        fake_keywords,
    )

    def fake_extract_cycle(
        **kwargs: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame, dict[int, int]]:
        calls.append("extract_cycle")
        captured["extract_cycle"].append(kwargs)
        recipe.qc = pd.concat([recipe.qc, extractionQc], ignore_index=True)
        return spectrum.copy(), spectrum.copy(), {10: 11}

    monkeypatch.setattr(
        recipe,
        "process_single_ab_nodding_cycle",
        fake_extract_cycle,
    )

    def fake_stack_extractions(
        *args: object, **kwargs: object
    ) -> tuple[pd.DataFrame, str]:
        calls.append("stack_extractions")
        captured["stack_extractions"].append({"args": args, **kwargs})
        return spectrum.copy(), str(extractionPath)

    monkeypatch.setattr(
        recipe,
        "stack_extractions",
        fake_stack_extractions,
    )
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(
        recipe,
        "clean_up",
        lambda **kwargs: (
            calls.append("clean_up"),
            captured["clean_up"].append(kwargs),
        )[-1],
    )

    toolkit = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)

    def fake_plot(**kwargs: object) -> tuple[pd.DataFrame, str]:
        calls.append("plot")
        captured["plot"].append(kwargs)
        return _spectrum_product(kwargs["products"], plotPath), str(plotPath)

    monkeypatch.setattr(
        toolkit,
        "plot_merged_spectrum_qc",
        fake_plot,
    )
    return calls, captured, extractionQc


def _assert_merged_product(products: pd.DataFrame, plotPath: Path) -> None:
    """Assert the complete stable merged-spectrum product contract."""
    assert products.iloc[-1].to_dict() == {
        "soxspipe_recipe": "observing-recipe",
        "product_label": "MERGED_SPECTRUM",
        "file_name": plotPath.name,
        "file_type": "FITS",
        "obs_date_utc": "2024-01-02T03:04:05.678",
        "reduction_date_utc": "2024-01-02T04:05:06.789",
        "product_desc": "Synthetic merged spectrum",
        "file_path": str(plotPath),
        "label": "PROD",
    }


def test_nod_cycle_extracts_both_difference_frames_without_flattening(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """An AB cycle preserves its two directional difference-spectrum contracts."""
    import soxspipe.commonutils as commonutils
    nod_module = import_module("soxspipe.recipes.soxs_nod")

    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.recipeName = "soxs-nod"
    recipe.sofName = "SYNTHETIC"
    recipe.debug = False
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.products = _empty_products()
    recipe.settings = {"instrument": "soxs"}
    recipe.recipeSettings = {"save_single_frame_extractions": False}
    recipe.twoDMap = "two-d-map.fits"
    recipe.dispMap = "dispersion-map.fits"
    recipe.startNightDate = "2024-01-02"
    recipe.turnOffMP = True
    recipe.arm = "VIS"
    recipe.kw = _keyword
    recipe.productDir = str(tmp_path)
    calls: list[dict[str, object]] = []

    class Extractor:
        def __init__(self, **kwargs: object) -> None:
            calls.append(kwargs)

        def extract(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, int], str]:
            spectrum = pd.DataFrame({"WAVE": [500.0], "FLUX": [1.0]})
            return recipe.qc, recipe.products, spectrum, {10: 1}, "extraction.fits"

    monkeypatch.setattr(commonutils, "horne_extraction", Extractor)
    monkeypatch.setattr(
        nod_module,
        "generic_quality_checks",
        lambda **kwargs: kwargs["qcTable"],
    )
    monkeypatch.setattr(
        nod_module,
        "spectroscopic_image_quality_checks",
        lambda **kwargs: kwargs["qcTable"],
    )
    aFrame = synthetic_ccd(shape=(3, 3), seed=11, prepared=True)
    bFrame = synthetic_ccd(shape=(3, 3), seed=12, prepared=True)

    spectrumA, spectrumB, joins = recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=1,
        orderTablePath="orders.fits",
        notFlattened=True,
    )

    assert len(calls) == 2
    assert calls[0]["skySubtractedFrame"].data[0, 0] == pytest.approx(
        aFrame.data[0, 0] - bFrame.data[0, 0]
    )
    assert calls[1]["skySubtractedFrame"].data[0, 0] == pytest.approx(
        bFrame.data[0, 0] - aFrame.data[0, 0]
    )
    assert spectrumA["WAVE"].tolist() == [500.0]
    assert spectrumB["WAVE"].tolist() == [500.0]
    assert joins == {10: 1}


def test_nod_stack_extractions_writes_merged_fits_and_ascii_products(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Merged AB spectra retain their median flux, SNR, and product records."""
    from astropy import units as u
    from astropy.io import fits
    from soxspipe.commonutils import toolkit

    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.recipeName = "soxs-nod"
    recipe.sofName = "SYNTHETIC"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.productDir = str(tmp_path)
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.products = _empty_products()
    recipe.masterHeaderFrame = synthetic_ccd(shape=(3, 3), prepared=True)
    recipe.update_fits_keywords = lambda **_: None
    monkeypatch.setattr(
        toolkit,
        "add_snr_efficiency_qcs",
        lambda **kwargs: kwargs["qcTable"],
    )
    first = pd.DataFrame(
        {
            "WAVE": [500.12344 * u.nm, 501.0 * u.nm],
            "FLUX_COUNTS": [9.0, 16.0],
            "VARIANCE": [9.0, 4.0],
            "FLUX_DENSITY_COUNTS": [1.0, 2.0],
        }
    )
    second = pd.DataFrame(
        {
            "WAVE": [500.12345 * u.nm, 501.0 * u.nm],
            "FLUX_COUNTS": [15.0, 20.0],
            "VARIANCE": [9.0, 4.0],
            "FLUX_DENSITY_COUNTS": [3.0, 4.0],
        }
    )

    stacked, fitsPath = recipe.stack_extractions([first, second], orderJoins={10: 1})

    assert Path(fitsPath).is_file()
    assert Path(fitsPath.replace(".fits", ".txt")).is_file()
    assert stacked["FLUX_COUNTS"].tolist() == [12.0, 18.0]
    assert stacked["SNR"].tolist() == [4.0, 9.0]
    assert recipe.products["product_label"].tolist() == [
        "EXTRACTED_MERGED_TABLE",
        "EXTRACTED_MERGED_ASCII",
    ]
    with fits.open(fitsPath) as hdus:
        assert hdus[0].header["ESO PRO CATG"] == "SCI_SLIT_FLUX_VIS"


def test_nod_success_records_qc_and_characterizes_none_product_path(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A balanced AB pair records QC while retaining the legacy None path."""
    frameA = _write_prepared_frame(
        tmp_path / "nod-a.fits",
        seed=41,
        headerOverrides={"ESO SEQ CUMOFF Y": 3.0},
    )
    frameB = _write_prepared_frame(
        tmp_path / "nod-b.fits",
        seed=42,
        headerOverrides={"ESO SEQ CUMOFF Y": -3.0},
    )
    recipe = soxs_nod.__new__(soxs_nod)
    _configure_nodding_recipe(
        recipe,
        log=log,
        tmp_path=tmp_path,
        objectPaths=[frameA, frameB],
        technique="ECHELLE,SLIT,NODDING",
    )
    originalQc = recipe.qc.copy(deep=True)
    plotPath = tmp_path / "OBJECT_VIS_MERGED_QC.pdf"
    calls, captured, extractionQc = _patch_nodding_collaborators(
        recipe,
        monkeypatch=monkeypatch,
        plotPath=plotPath,
        extractionPath=tmp_path / "OBJECT_VIS_EXTRACTED.fits",
    )

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath is None
    expectedQc = pd.concat([originalQc, extractionQc], ignore_index=True)
    pd.testing.assert_frame_equal(returnedQc, expectedQc)
    _assert_merged_product(recipe.products, plotPath)
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
    assert [entry["recipe"] for entry in captured["stack"]] == [
        "soxs_nod",
        "soxs_nod",
    ]
    assert captured["stack"][0]["frames"][0].header["ESO SEQ CUMOFF Y"] == 3.0
    assert captured["stack"][1]["frames"][0].header["ESO SEQ CUMOFF Y"] == -3.0
    extractArgs = captured["extract_cycle"][0]
    assert extractArgs["aFrame"] is captured["keywords"][0]["frame"]
    assert extractArgs["bFrame"] is captured["keywords"][1]["frame"]
    assert extractArgs["locationSetIndex"] == 1
    assert extractArgs["orderTablePath"] == str(tmp_path / "ORDER_TAB_VIS.fits")
    assert extractArgs["masterFlat"] is False
    plotArgs = captured["plot"][0]
    assert plotArgs["filenameTemplate"] == "OBJECT_VIS.fits"
    assert plotArgs["settings"] is recipe.settings
    assert plotArgs["qcTable"] is recipe.qc
    assert captured["clean_up"] == [{"forceFail": False}]


def test_offset_success_returns_extraction_and_records_qc(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A balanced ON/OFF pair returns its path and complete extraction QC."""
    frameOn = _write_prepared_frame(
        tmp_path / "offset-on.fits",
        seed=43,
        headerOverrides={
            "HIERARCH ESO SEQ FIXOFF RA": -2.0,
            "HIERARCH ESO SEQ FIXOFF DEC": 0.0,
        },
    )
    frameOff = _write_prepared_frame(
        tmp_path / "offset-off.fits",
        seed=44,
        headerOverrides={
            "HIERARCH ESO SEQ FIXOFF RA": 0.0,
            "HIERARCH ESO SEQ FIXOFF DEC": 0.0,
        },
    )
    recipe = soxs_offset.__new__(soxs_offset)
    _configure_nodding_recipe(
        recipe,
        log=log,
        tmp_path=tmp_path,
        objectPaths=[frameOn, frameOff],
        technique="ECHELLE,SLIT,OFFSET",
    )
    originalQc = recipe.qc.copy(deep=True)
    extractionPath = tmp_path / "OBJECT_VIS_EXTRACTED.fits"
    plotPath = tmp_path / "OBJECT_VIS_MERGED_QC.pdf"
    calls, captured, extractionQc = _patch_nodding_collaborators(
        recipe,
        monkeypatch=monkeypatch,
        plotPath=plotPath,
        extractionPath=extractionPath,
    )

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == str(extractionPath)
    expectedQc = pd.concat([originalQc, extractionQc], ignore_index=True)
    pd.testing.assert_frame_equal(returnedQc, expectedQc)
    _assert_merged_product(recipe.products, plotPath)
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
    assert [entry["recipe"] for entry in captured["stack"]] == [
        "soxs_offset",
        "soxs_offset",
    ]
    assert (
        captured["stack"][0]["frames"][0].header["HIERARCH ESO SEQ FIXOFF RA"] == -2.0
    )
    assert captured["stack"][1]["frames"][0].header["HIERARCH ESO SEQ FIXOFF RA"] == 0.0
    extractArgs = captured["extract_cycle"][0]
    assert extractArgs["aFrame"] is captured["keywords"][0]["frame"]
    assert extractArgs["bFrame"] is captured["keywords"][1]["frame"]
    assert extractArgs["locationSetIndex"] == 1
    assert extractArgs["orderTablePath"] == str(tmp_path / "ORDER_TAB_VIS.fits")
    assert extractArgs["masterFlat"] is False
    stackArgs = captured["stack_extractions"][0]
    assert len(stackArgs["args"][0]) == 1
    assert stackArgs["orderJoins"] == {10: 11}
    plotArgs = captured["plot"][0]
    assert plotArgs["merged_orders"].equals(
        pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})
    )
    assert plotArgs["filenameTemplate"] == "OBJECT_VIS.fits"
    assert plotArgs["settings"] is recipe.settings
    assert plotArgs["qcTable"] is recipe.qc
    assert captured["clean_up"] == [{"forceFail": False}]


def test_offset_multiple_locations_extracts_each_pair_before_final_stack(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Distinct offset locations retain one directional extraction per pair."""
    headers = [
        {"HIERARCH ESO SEQ FIXOFF RA": -3.0, "HIERARCH ESO SEQ FIXOFF DEC": 1.0},
        {"HIERARCH ESO SEQ FIXOFF RA": -4.0, "HIERARCH ESO SEQ FIXOFF DEC": 2.0},
        {"HIERARCH ESO SEQ FIXOFF RA": 0.0, "HIERARCH ESO SEQ FIXOFF DEC": 1.0},
        {"HIERARCH ESO SEQ FIXOFF RA": 0.0, "HIERARCH ESO SEQ FIXOFF DEC": 2.0},
    ]
    paths = [
        _write_prepared_frame(
            tmp_path / f"offset-{index}.fits",
            seed=90 + index,
            headerOverrides={**header, "MJD-OBS": 60_000.0 + index},
        )
        for index, header in enumerate(headers)
    ]
    recipe = soxs_offset.__new__(soxs_offset)
    _configure_nodding_recipe(
        recipe,
        log=log,
        tmp_path=tmp_path,
        objectPaths=paths,
        technique="ECHELLE,SLIT,OFFSET",
    )
    calls, captured, _ = _patch_nodding_collaborators(
        recipe,
        monkeypatch=monkeypatch,
        plotPath=tmp_path / "OBJECT_VIS_MERGED_QC.pdf",
        extractionPath=tmp_path / "OBJECT_VIS_EXTRACTED.fits",
    )

    productPath, _ = recipe.produce_product()

    assert productPath == str(tmp_path / "OBJECT_VIS_EXTRACTED.fits")
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
    assert len(captured["stack_extractions"][0]["args"]) == 1
    assert len(captured["stack_extractions"][0]["args"][0][0]) == 2
    assert captured["clean_up"] == [{"forceFail": False}]


@pytest.mark.parametrize(
    ("recipeClass", "technique", "offsets", "message"),
    [
        (
            soxs_nod,
            "ECHELLE,SLIT,NODDING",
            [{"ESO SEQ CUMOFF Y": 3.0}, {"ESO SEQ CUMOFF Y": 2.0}],
            "Found 2 A frames and 0 B frames. The number of A and B frames must be the same for nodding reductions.",
        ),
        (
            soxs_offset,
            "ECHELLE,SLIT,OFFSET",
            [
                {
                    "HIERARCH ESO SEQ FIXOFF RA": -2.0,
                    "HIERARCH ESO SEQ FIXOFF DEC": 0.0,
                },
                {
                    "HIERARCH ESO SEQ FIXOFF RA": -1.0,
                    "HIERARCH ESO SEQ FIXOFF DEC": 0.0,
                },
            ],
            "Found 2 ON frames and 0 OFF frames. The number of ON and OFF frames must be the same for offset reductions.",
        ),
    ],
)
def test_unbalanced_observing_inventory_fails_with_counts(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    recipeClass: type[Any],
    technique: str,
    offsets: list[dict[str, object]],
    message: str,
) -> None:
    """Unbalanced observing inventories fail before scientific processing."""
    objectPaths = [
        _write_prepared_frame(
            tmp_path / f"unbalanced-{index}.fits",
            seed=50 + index,
            headerOverrides=header,
        )
        for index, header in enumerate(offsets)
    ]
    recipe = recipeClass.__new__(recipeClass)
    _configure_nodding_recipe(
        recipe,
        log=log,
        tmp_path=tmp_path,
        objectPaths=objectPaths,
        technique=technique,
    )
    toolkit = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)

    with pytest.raises(Exception) as error:
        recipe.produce_product()

    assert type(error.value) is Exception
    assert str(error.value) == message


def test_stare_success_returns_last_sky_path_and_records_products(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Stare returns the residual path and records all sky products."""
    objectPath = _write_prepared_frame(
        tmp_path / "stare-object.fits",
        seed=61,
        headerOverrides={
            "DPR_TYPE": "OBJECT",
            "DPR_TECH": "ECHELLE,SLIT,STARE",
        },
    )
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = _keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RouteCollection(
        {
            _route(DPR_TYPE="OBJECT", DPR_TECH="ECHELLE,SLIT,STARE"): [str(objectPath)],
            _route(PRO_CATG="ORDER_TAB_VIS"): [str(tmp_path / "ORDER_TAB_VIS.fits")],
            _route(PRO_CATG="DISP_TAB_VIS"): [str(tmp_path / "DISP_TAB_VIS.fits")],
            _route(PRO_CATG="DISP_IMAGE_VIS"): [str(tmp_path / "DISP_IMAGE_VIS.fits")],
        }
    )
    recipe.recipeName = "soxs-stare"
    recipe.recipeSettings = {
        "use_flat": False,
        "sky-subtraction": {"subtract_sky": True},
    }
    recipe.settings = pipeline_settings(tmp_path)
    recipe.generateReponseCurve = False
    recipe.sofName = "OBJECT_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.workspaceRootPath = str(tmp_path)
    recipe.qcDir = str(tmp_path / "qc")
    recipe.filenameTemplate = "OBJECT_VIS.fits"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.qc = qc_table()
    recipe.products = _empty_products()
    originalQc = recipe.qc.copy(deep=True)
    genericQc = qc_row(
        recipeName="soxs-stare",
        name="BADPIX NUM",
        value=4,
        unit="pixel",
        comment="Synthetic bad-pixel count",
    )
    spectroscopicQc = qc_row(
        recipeName="soxs-stare",
        name="ORDER CENTRE RESIDUAL",
        value=0.25,
        unit="pixel",
        comment="Synthetic order-centre residual",
    )
    combinedFrame = synthetic_ccd(seed=62, prepared=True)
    skyModel = synthetic_ccd(seed=63, prepared=True)
    skySubtracted = synthetic_ccd(seed=64, prepared=True)
    residuals = synthetic_ccd(seed=65, prepared=True)
    calls: list[str] = []
    captured: dict[str, list[dict[str, object]]] = {
        "stack": [],
        "detrend": [],
        "keywords": [],
        "sky": [],
        "write": [],
        "generic_qc": [],
        "spectroscopic_qc": [],
        "extractor": [],
        "plot": [],
        "clean_up": [],
    }

    def fake_stack(**kwargs: object) -> Any:
        calls.append("stack")
        result = combinedFrame.copy()
        captured["stack"].append({**kwargs, "result": result})
        return result

    monkeypatch.setattr(
        recipe,
        "clip_and_stack",
        fake_stack,
    )

    def fake_detrend(**kwargs: object) -> Any:
        calls.append("detrend")
        result = combinedFrame.copy()
        captured["detrend"].append({**kwargs, "result": result})
        return result

    monkeypatch.setattr(
        recipe,
        "detrend",
        fake_detrend,
    )

    def fake_keywords(**kwargs: object) -> None:
        calls.append("keywords")
        captured["keywords"].append(kwargs)

    monkeypatch.setattr(
        recipe,
        "update_fits_keywords",
        fake_keywords,
    )
    writtenPaths: list[str] = []

    def fake_write(*args: object, **kwargs: object) -> str:
        calls.append("write")
        path = str(tmp_path / str(kwargs["filename"]))
        writtenPaths.append(path)
        captured["write"].append({"args": args, **kwargs, "path": path})
        return path

    monkeypatch.setattr(recipe, "_write", fake_write)
    monkeypatch.setattr(
        recipe,
        "report_output",
        lambda: calls.append("report") or recipe.qc,
    )
    monkeypatch.setattr(
        recipe,
        "clean_up",
        lambda **kwargs: (
            calls.append("clean_up"),
            captured["clean_up"].append(kwargs),
        )[-1],
    )

    stareModule = import_module("soxspipe.recipes.soxs_stare")

    class FakeSkySubtractor:
        def __init__(self, **kwargs: object) -> None:
            calls.append("sky")
            captured["sky"].append(kwargs)

        def subtract(self) -> tuple[Any, Any, Any, pd.DataFrame, pd.DataFrame]:
            calls.append("sky_subtract")
            return skyModel, skySubtracted, residuals, recipe.qc, recipe.products

    class FakeExtractor:
        def __init__(self, **kwargs: object) -> None:
            calls.append("extractor")
            captured["extractor"].append(kwargs)

        def extract(self) -> tuple[Any, Any, pd.DataFrame, dict[int, int], str]:
            calls.append("extract")
            spectrum = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})
            return (
                recipe.qc,
                recipe.products,
                spectrum,
                {10: 11},
                str(tmp_path / "OBJECT_VIS_EXTRACTED.fits"),
            )

    monkeypatch.setattr(stareModule, "subtract_sky", FakeSkySubtractor)

    def fake_generic_qc(**kwargs: object) -> pd.DataFrame:
        calls.append("generic_qc")
        captured["generic_qc"].append(kwargs)
        assert kwargs["qcTable"] is recipe.qc
        return pd.concat([kwargs["qcTable"], genericQc], ignore_index=True)

    def fake_spectroscopic_qc(**kwargs: object) -> pd.DataFrame:
        calls.append("spectroscopic_qc")
        captured["spectroscopic_qc"].append(kwargs)
        assert kwargs["qcTable"] is recipe.qc
        return pd.concat([kwargs["qcTable"], spectroscopicQc], ignore_index=True)

    monkeypatch.setattr(stareModule, "generic_quality_checks", fake_generic_qc)
    monkeypatch.setattr(
        stareModule,
        "spectroscopic_image_quality_checks",
        fake_spectroscopic_qc,
    )
    commonutils = import_module("soxspipe.commonutils")
    monkeypatch.setattr(commonutils, "horne_extraction", FakeExtractor)
    toolkit = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)
    plotPath = tmp_path / "OBJECT_VIS_MERGED_QC.pdf"

    def fake_plot(**kwargs: object) -> tuple[pd.DataFrame, str]:
        calls.append("plot")
        captured["plot"].append(kwargs)
        return _spectrum_product(kwargs["products"], plotPath), str(plotPath)

    monkeypatch.setattr(
        toolkit,
        "plot_merged_spectrum_qc",
        fake_plot,
    )

    returnedPath, returnedQc = recipe.produce_product()

    assert returnedPath == writtenPaths[-1]
    assert Path(returnedPath).name == "OBJECT_VIS_SKYSUB_RESIDUALS.fits"
    expectedQc = pd.concat([originalQc, genericQc, spectroscopicQc], ignore_index=True)
    pd.testing.assert_frame_equal(returnedQc, expectedQc)
    assert recipe.products["product_label"].tolist() == [
        "SKY_SUBTRACTED_OBJECT",
        "SKY_MODEL",
        "SKY_SUB_RESIDUALS",
        "MERGED_SPECTRUM",
    ]
    expectedSkyProducts = [
        {
            "soxspipe_recipe": "soxs-stare",
            "product_label": "SKY_SUBTRACTED_OBJECT",
            "file_name": "OBJECT_VIS_SKYSUB.fits",
            "file_type": "FITS",
            "obs_date_utc": "2024-01-02T03:04:05.678",
            "product_desc": "The sky-subtracted object",
            "file_path": str(tmp_path / "OBJECT_VIS_SKYSUB.fits"),
            "label": "PROD",
        },
        {
            "soxspipe_recipe": "soxs-stare",
            "product_label": "SKY_MODEL",
            "file_name": "OBJECT_VIS_SKYMODEL.fits",
            "file_type": "FITS",
            "obs_date_utc": "2024-01-02T03:04:05.678",
            "product_desc": "The sky background model",
            "file_path": str(tmp_path / "OBJECT_VIS_SKYMODEL.fits"),
            "label": "PROD",
        },
        {
            "soxspipe_recipe": "soxs-stare",
            "product_label": "SKY_SUB_RESIDUALS",
            "file_name": "OBJECT_VIS_SKYSUB_RESIDUALS.fits",
            "file_type": "FITS",
            "obs_date_utc": "2024-01-02T03:04:05.678",
            "product_desc": "The sky subtraction residuals",
            "file_path": str(tmp_path / "OBJECT_VIS_SKYSUB_RESIDUALS.fits"),
            "label": "PROD",
        },
    ]
    for row, expected in zip(
        recipe.products.iloc[:3].to_dict("records"), expectedSkyProducts
    ):
        reductionDate = row.pop("reduction_date_utc")
        assert row == expected
        assert len(reductionDate) == 19
        assert reductionDate[4] == "-" and reductionDate[10] == "T"
    _assert_merged_product(recipe.products, plotPath)
    assert calls == [
        "stack",
        "detrend",
        "keywords",
        "sky",
        "sky_subtract",
        "write",
        "write",
        "write",
        "generic_qc",
        "spectroscopic_qc",
        "extractor",
        "extract",
        "plot",
        "report",
        "clean_up",
    ]
    stackArgs = captured["stack"][0]
    assert stackArgs["recipe"] == "soxs_stare"
    assert stackArgs["frames"][0].header["DATE-OBS"] == recipe.dateObs
    detrendArgs = captured["detrend"][0]
    assert detrendArgs["inputFrame"] is stackArgs["result"]
    assert captured["keywords"][0]["frame"] is detrendArgs["result"]
    assert detrendArgs["master_bias"] is False
    assert detrendArgs["dark"] is False
    assert detrendArgs["master_flat"] is False
    assert detrendArgs["order_table"] == str(tmp_path / "ORDER_TAB_VIS.fits")
    skyArgs = captured["sky"][0]
    assert skyArgs["objectFrame"] is captured["keywords"][0]["frame"]
    assert skyArgs["twoDMap"] == str(tmp_path / "DISP_IMAGE_VIS.fits")
    assert skyArgs["dispMap"] == str(tmp_path / "DISP_TAB_VIS.fits")
    assert skyArgs["settings"] is recipe.settings
    pd.testing.assert_frame_equal(skyArgs["qcTable"], originalQc)
    expectedFrames = [skySubtracted, skyModel, residuals]
    expectedFilenames = [
        "OBJECT_VIS_SKYSUB.fits",
        "OBJECT_VIS_SKYMODEL.fits",
        "OBJECT_VIS_SKYSUB_RESIDUALS.fits",
    ]
    assert all(
        entry["frame"] is expected
        for entry, expected in zip(captured["write"], expectedFrames)
    )
    assert [entry["filename"] for entry in captured["write"]] == expectedFilenames
    assert all(entry["filedir"] == str(tmp_path) for entry in captured["write"])
    assert all(entry["overwrite"] is True for entry in captured["write"])
    assert captured["write"][0]["maskToZero"] is True
    for entry in captured["write"]:
        assert entry["frame"].header["DATE-OBS"] == "2024-01-02T03:04:05.678"
        assert str(entry["frame"].unit) == "electron"
        assert entry["frame"].uncertainty is not None
    extractorArgs = captured["extractor"][0]
    assert extractorArgs["skySubtractedFrame"] is skySubtracted
    assert extractorArgs["subtractedFrame"] is skyModel
    assert extractorArgs["twoDMapPath"] == str(tmp_path / "DISP_IMAGE_VIS.fits")
    assert extractorArgs["dispersionMap"] == str(tmp_path / "DISP_TAB_VIS.fits")
    assert extractorArgs["recipeSettings"] is recipe.recipeSettings
    assert extractorArgs["productsTable"] is captured["plot"][0]["products"]
    plotArgs = captured["plot"][0]
    assert plotArgs["filenameTemplate"] == "OBJECT_VIS.fits"
    assert plotArgs["qcDir"] == str(tmp_path / "qc")
    assert plotArgs["settings"] is recipe.settings
    assert captured["clean_up"] == [{"forceFail": False}]


def test_stare_uses_flux_standard_frames_when_object_frames_are_absent(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A flux standard is the established fallback input for response generation."""
    standardPath = _write_prepared_frame(
        tmp_path / "flux-standard.fits",
        seed=71,
        headerOverrides={
            "DPR_TYPE": "STD,FLUX",
            "DPR_TECH": "ECHELLE,SLIT,STARE",
        },
    )
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = _keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RouteCollection(
        {
            _route(DPR_TYPE="OBJECT", DPR_TECH="ECHELLE,SLIT,STARE"): [],
            _route(DPR_TYPE="OBJECT,ASYNC", DPR_TECH="ECHELLE,SLIT,STARE"): [],
            _route(DPR_TYPE="STD,FLUX", DPR_TECH="ECHELLE,SLIT,STARE"): [str(standardPath)],
        }
    )
    recipe.recipeSettings = {"sky-subtraction": {"subtract_sky": False}}
    recipe.settings = pipeline_settings(tmp_path)
    recipe.generateReponseCurve = False

    class StopAfterInputSelection(Exception):
        """Terminate the test at the first downstream orchestration boundary."""

    def stop_after_stack(**kwargs: object) -> None:
        assert len(kwargs["frames"]) == 1
        raise StopAfterInputSelection()

    monkeypatch.setattr(recipe, "clip_and_stack", stop_after_stack)

    with pytest.raises(StopAfterInputSelection):
        recipe.produce_product()

    assert recipe.generateReponseCurve is True


@pytest.mark.parametrize(
    ("frameType", "technique", "settingsOverrides", "expectedSkySubtraction"),
    [
        ("STD,TELLURIC", "ECHELLE,SLIT,STARE", {}, True),
        ("LAMP,FLAT", "ECHELLE,PINHOLE", {"PAE": True}, False),
        ("STD,FLUX", "ECHELLE,SLIT,NODDING", {}, True),
    ],
)
def test_stare_selects_supported_non_object_fallbacks(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    frameType: str,
    technique: str,
    settingsOverrides: dict[str, object],
    expectedSkySubtraction: bool,
) -> None:
    """Each documented non-object frame type reaches the shared stack boundary."""
    fallbackPath = _write_prepared_frame(
        tmp_path / f"{frameType.replace(',', '_')}.fits",
        seed=74,
        headerOverrides={"DPR_TYPE": frameType, "DPR_TECH": technique},
    )
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = _keyword
    recipe.detectorParams = {}
    routes = {
        _route(DPR_TYPE="OBJECT", DPR_TECH="ECHELLE,SLIT,STARE"): [],
        _route(DPR_TYPE="OBJECT,ASYNC", DPR_TECH="ECHELLE,SLIT,STARE"): [],
        _route(DPR_TYPE="STD,FLUX", DPR_TECH="ECHELLE,SLIT,STARE"): [],
        _route(DPR_TYPE="STD,TELLURIC", DPR_TECH="ECHELLE,SLIT,STARE"): [],
        _route(DPR_TYPE="LAMP,FLAT", DPR_TECH="ECHELLE,PINHOLE"): [],
        _route(DPR_TYPE="STD,FLUX", DPR_TECH="ECHELLE,SLIT,NODDING"): [],
        _route(DPR_TYPE=frameType, DPR_TECH=technique): [str(fallbackPath)],
    }
    recipe.inputFrames = RouteCollection(routes)
    recipe.recipeSettings = {"sky-subtraction": {"subtract_sky": True}}
    recipe.settings = {**pipeline_settings(tmp_path), **settingsOverrides}
    recipe.generateReponseCurve = False

    class StopAfterInputSelection(Exception):
        """Terminate at the first shared downstream collaborator."""

    def stop_after_stack(**kwargs: object) -> None:
        assert len(kwargs["frames"]) == 1
        raise StopAfterInputSelection()

    monkeypatch.setattr(recipe, "clip_and_stack", stop_after_stack)

    with pytest.raises(StopAfterInputSelection):
        recipe.produce_product()

    assert recipe.subtractSky is expectedSkySubtraction


def test_stare_generates_a_response_from_a_flux_standard(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A flux standard takes the unflattened extraction through response creation."""
    standardPath = _write_prepared_frame(
        tmp_path / "response-standard.fits",
        seed=72,
        headerOverrides={
            "DPR_TYPE": "STD,FLUX",
            "DPR_TECH": "ECHELLE,SLIT,STARE",
        },
    )
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = _keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RouteCollection(
        {
            _route(DPR_TYPE="OBJECT", DPR_TECH="ECHELLE,SLIT,STARE"): [],
            _route(DPR_TYPE="OBJECT,ASYNC", DPR_TECH="ECHELLE,SLIT,STARE"): [],
            _route(DPR_TYPE="STD,FLUX", DPR_TECH="ECHELLE,SLIT,STARE"): [str(standardPath)],
            _route(PRO_CATG="ORDER_TAB_VIS"): [str(tmp_path / "ORDER_TAB_VIS.fits")],
            _route(PRO_CATG="DISP_TAB_VIS"): [str(tmp_path / "DISP_TAB_VIS.fits")],
            _route(PRO_CATG="DISP_IMAGE_VIS"): [str(tmp_path / "DISP_IMAGE_VIS.fits")],
        }
    )
    recipe.recipeName = "soxs-stare"
    recipe.recipeSettings = {"use_flat": False, "sky-subtraction": {"subtract_sky": False}}
    recipe.settings = pipeline_settings(tmp_path)
    recipe.generateReponseCurve = False
    recipe.sofName = "STANDARD_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.workspaceRootPath = str(tmp_path)
    recipe.qcDir = str(tmp_path / "qc")
    recipe.filenameTemplate = "STANDARD_VIS.fits"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.qc = qc_table()
    recipe.products = _empty_products()
    combinedFrame = synthetic_ccd(seed=73, prepared=True)
    calls: list[str] = []
    extractionArguments: list[dict[str, object]] = []
    responseArguments: list[dict[str, object]] = []

    monkeypatch.setattr(recipe, "clip_and_stack", lambda **_: combinedFrame.copy())
    monkeypatch.setattr(recipe, "detrend", lambda **_: combinedFrame.copy())
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **_: calls.append("keywords"))
    monkeypatch.setattr(recipe, "report_output", lambda: calls.append("report") or recipe.qc)
    monkeypatch.setattr(
        recipe,
        "clean_up",
        lambda **kwargs: calls.append(f"clean_up:{kwargs['forceFail']}"),
    )

    class FakeExtractor:
        def __init__(self, **kwargs: object) -> None:
            extractionArguments.append(kwargs)
            calls.append("extractor")

        def extract(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, int], str]:
            calls.append("extract")
            extractionNumber = len(extractionArguments)
            spectrum = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0], "SKY_COUNTS": [3.0]})
            return (
                recipe.qc,
                recipe.products,
                spectrum,
                {10: 11},
                str(tmp_path / f"extracted-{extractionNumber}.fits"),
            )

    class FakeResponse:
        def __init__(self, **kwargs: object) -> None:
            responseArguments.append(kwargs)
            calls.append("response")

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, str]:
            calls.append("response_get")
            return recipe.qc, recipe.products, "synthetic response mismatch"

    commonutils = import_module("soxspipe.commonutils")
    toolkit = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(commonutils, "horne_extraction", FakeExtractor)
    monkeypatch.setattr(commonutils, "response_function", FakeResponse)
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **_: None)
    monkeypatch.setattr(
        toolkit,
        "plot_merged_spectrum_qc",
        lambda **kwargs: (kwargs["products"], str(tmp_path / "merged.pdf")),
    )

    productPath, returnedQc = recipe.produce_product()

    assert productPath is None
    assert returnedQc is recipe.qc
    assert recipe.generateReponseCurve is True
    assert calls == [
        "keywords",
        "extractor",
        "extract",
        "extractor",
        "extract",
        "response",
        "response_get",
        "report",
        "clean_up:synthetic response mismatch",
    ]
    assert extractionArguments[0].get("notFlattened", False) is False
    assert extractionArguments[1]["notFlattened"] is True
    assert responseArguments[0]["stdExtractionPath"] == str(tmp_path / "extracted-1.fits")
    assert responseArguments[0]["stdNotFlatExtractionPath"] == str(tmp_path / "extracted-2.fits")
    assert responseArguments[0]["orderJoins"] == {10: 11}


def test_stare_flux_calibrates_an_extracted_spectrum_when_response_exists(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Stare applies a supplied response and creates calibrated QC output."""
    from astropy.table import Table

    objectPath = _write_prepared_frame(
        tmp_path / "stare-standard.fits",
        seed=66,
        headerOverrides={
            "DPR_TYPE": "OBJECT",
            "DPR_TECH": "ECHELLE,SLIT,STARE",
            "HIERARCH ESO TEL AIRM END": 1.2,
            "EXPTIME": 120.0,
        },
    )
    responsePath = tmp_path / "RESP_TAB_VIS.fits"
    responsePath.touch()
    fluxPath = tmp_path / "OBJECT_VIS_FLUXCAL.fits"
    Table({"WAVE": [500.0], "FLUX_CALIBRATED": [2.5]}).write(fluxPath)
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = _keyword
    recipe.detectorParams = {}
    recipe.inputFrames = RouteCollection(
        {
            _route(DPR_TYPE="OBJECT", DPR_TECH="ECHELLE,SLIT,STARE"): [str(objectPath)],
            _route(PRO_CATG="ORDER_TAB_VIS"): [str(tmp_path / "ORDER_TAB_VIS.fits")],
            _route(PRO_CATG="DISP_TAB_VIS"): [str(tmp_path / "DISP_TAB_VIS.fits")],
            _route(PRO_CATG="DISP_IMAGE_VIS"): [str(tmp_path / "DISP_IMAGE_VIS.fits")],
            _route(PRO_CATG="RESP_TAB_VIS"): [str(responsePath)],
        }
    )
    recipe.recipeName = "soxs-stare"
    recipe.recipeSettings = {"use_flat": False, "sky-subtraction": {"subtract_sky": False}}
    recipe.settings = pipeline_settings(tmp_path)
    recipe.generateReponseCurve = False
    recipe.sofName = "OBJECT_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.workspaceRootPath = str(tmp_path)
    recipe.qcDir = str(tmp_path / "qc")
    recipe.filenameTemplate = "OBJECT_VIS.fits"
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.qc = qc_table()
    recipe.products = _empty_products()
    combinedFrame = synthetic_ccd(seed=67, prepared=True)
    combinedFrame.header["HIERARCH ESO TEL AIRM END"] = 1.2
    combinedFrame.header["EXPTIME"] = 120.0
    calls: list[str] = []

    monkeypatch.setattr(recipe, "clip_and_stack", lambda **_: combinedFrame.copy())
    monkeypatch.setattr(recipe, "detrend", lambda **_: combinedFrame.copy())
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **_: calls.append("keywords"))
    monkeypatch.setattr(recipe, "report_output", lambda: calls.append("report") or recipe.qc)
    monkeypatch.setattr(recipe, "clean_up", lambda **kwargs: calls.append(f"clean_up:{kwargs['forceFail']}"))

    class FakeExtractor:
        def __init__(self, **kwargs: object) -> None:
            calls.append("extractor")

        def extract(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, int], str]:
            calls.append("extract")
            spectrum = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0], "SKY_COUNTS": [3.0]})
            return recipe.qc, recipe.products, spectrum, {10: 11}, str(tmp_path / "extracted.fits")

    class FakeFluxCalibrator:
        def __init__(self, **kwargs: object) -> None:
            assert kwargs["responseFunction"] == str(responsePath)
            assert kwargs["airmass"] == 1.2
            assert kwargs["exptime"] == 120.0
            calls.append("flux_calibrator")

        def calibrate(self) -> tuple[str, pd.DataFrame]:
            calls.append("calibrate")
            return str(fluxPath), _empty_products()

    commonutils = import_module("soxspipe.commonutils")
    stareModule = import_module("soxspipe.recipes.soxs_stare")
    toolkit = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(commonutils, "horne_extraction", FakeExtractor)
    monkeypatch.setattr(commonutils, "flux_calibration", FakeFluxCalibrator)
    monkeypatch.setattr(stareModule, "detector_lookup", lambda **_: type("Lookup", (), {"get": lambda _, __: {"extinction": "extinction.dat"}})())
    monkeypatch.setattr(stareModule, "get_calibrations_path", lambda **_: str(tmp_path))
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **_: None)

    def fake_plot(**kwargs: object) -> tuple[pd.DataFrame, str]:
        calls.append(f"plot:{kwargs['fluxCalibrated']}")
        return kwargs["products"], str(tmp_path / "merged.pdf")

    monkeypatch.setattr(toolkit, "plot_merged_spectrum_qc", fake_plot)

    productPath, returnedQc = recipe.produce_product()

    assert productPath is None
    assert returnedQc is recipe.qc
    assert calls == [
        "keywords",
        "extractor",
        "extract",
        "flux_calibrator",
        "calibrate",
        "plot:False",
        "plot:True",
        "report",
        "clean_up:False",
    ]
