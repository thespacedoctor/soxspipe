"""File-name contracts for the continuum-trace residual PDF and OBJECT_TAB table."""

from __future__ import annotations

from pathlib import Path

import importlib

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.recipes.soxs_offset import STACKED_LOCATION_SET
from tests.factories import instrument_header, pipeline_settings, qc_table

pytestmark = pytest.mark.unit

SOF_NAME = "SOXS_VIS_FLAT_SOF"
POLY_ORDERS = "01"
QC_DIR_NAME = "qc"
PRODUCT_DIR_NAME = "products"


def _detector(log: object, tmp_path: Path, *, recipeName: str, noddingSequence: str = "") -> detect_continuum:
    """Return a continuum detector holding the state the file-naming code reads."""
    settings = pipeline_settings(tmp_path)
    settings["tune-pipeline"] = False
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.settings = settings
    detector.kw = keyword_lookup(log=log, settings=settings).get
    detector.recipeName = recipeName
    detector.noddingSequence = noddingSequence
    detector.sofName = SOF_NAME
    detector.lampTag = False
    detector.inst = "SOXS"
    detector.arm = "VIS"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.orderDeg = 0
    detector.axisBDeg = 1
    detector.sliceLength = 10
    detector.exptime = 60.0
    detector.lamp = ""
    detector.debug = False
    detector.recipeSettings = {}
    detector.qc = qc_table()
    detector.dateObs = "2024-01-02T03:04:05.678"
    detector.products = pd.DataFrame()
    detector.qcDir = str(tmp_path / QC_DIR_NAME)
    detector.productDir = str(tmp_path / PRODUCT_DIR_NAME)
    Path(detector.qcDir).mkdir(exist_ok=True)
    detector.detectorParams = {"rotate-qc-plot": 0, "flip-qc-plot": 0, "dispersion-axis": "x"}
    detector.traceFrame = CCDData(
        np.full((32, 32), 20.0),
        unit=u.electron,
        meta=instrument_header(arm="VIS", overrides={"EXPTIME": 60.0}),
    )
    return detector


def _plot_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return a one-order pixel table, polynomial table and clipped-sample table for ``plot_results``."""
    yValues = np.arange(3.0, 30.0, 3.0)
    pixels = pd.DataFrame(
        {
            "order": np.full(len(yValues), 10.0),
            "cont_x": np.full(len(yValues), 16.0),
            "cont_y": yValues,
            "cont_x_fit_res": np.zeros(len(yValues)),
            "wavelength": np.linspace(500.0, 510.0, len(yValues)),
            "gauss_stddev_fit": np.full(len(yValues), 2.0),
        }
    )
    polynomials = pd.DataFrame([{"cent_00": 16.0, "cent_01": 0.0, "std_00": 2.0, "std_01": 0.0}])
    clipped = pd.DataFrame({"cont_x": [17.0], "cont_y": [6.0], "fit_x": [16.0], "fit_y": [6.0]})
    return pixels, polynomials, clipped


def _residual_pdf_name(detector: detect_continuum, monkeypatch: pytest.MonkeyPatch) -> str:
    """Run the real ``plot_results`` and return the name of the PDF it writes."""
    monkeypatch.setattr("soxspipe.commonutils.toolkit.qc_settings_plot_tables", lambda **_: None)
    pixels, polynomials, clipped = _plot_inputs()
    filePath, _ = detector.plot_results(orderPixelTable=pixels, orderPolyTable=polynomials, clippedData=clipped)
    assert Path(filePath).read_bytes().startswith(b"%PDF")
    return Path(filePath).name


def _order_table_name(detector: detect_continuum, monkeypatch: pytest.MonkeyPatch) -> str:
    """Run the real ``write_order_table_to_file`` and return the name of the table it targets."""
    written: list[str] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: written.append(kwargs["filePath"]),
    )
    detector.write_order_table_to_file(
        frame=detector.traceFrame,
        orderPolyTable=pd.DataFrame(),
        orderMetaTable=pd.DataFrame(),
    )
    return Path(written[0]).name


@pytest.mark.parametrize(
    ("recipeName", "noddingSequence", "expectedName"),
    [
        ("soxs-stare", "", f"{SOF_NAME}_OBJECT_TRACE_residuals_{POLY_ORDERS}.pdf"),
        ("soxs-nod", "_A1", f"{SOF_NAME}_OBJECT_TRACE_residuals_A1_{POLY_ORDERS}.pdf"),
        ("soxs-order-centres", "", f"{SOF_NAME}_residuals_{POLY_ORDERS}.pdf"),
        ("soxs-mflat", "", f"{SOF_NAME}_OBJECT_TRACE_residuals_{POLY_ORDERS}.pdf"),
        ("soxs-nod-std", "_A1", f"{SOF_NAME}_OBJECT_TRACE_residuals_A1_{POLY_ORDERS}.pdf"),
        ("soxs-offset-std", "_B2", f"{SOF_NAME}_OBJECT_TRACE_residuals_B2_{POLY_ORDERS}.pdf"),
        ("soxs-offset", "_A1", f"{SOF_NAME}_OBJECT_TRACE_residuals_A1_{POLY_ORDERS}.pdf"),
        ("soxs-offset", "_B2", f"{SOF_NAME}_OBJECT_TRACE_residuals_B2_{POLY_ORDERS}.pdf"),
    ],
)
def test_residual_pdf_name_follows_recipe_and_nodding_sequence(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    recipeName: str,
    noddingSequence: str,
    expectedName: str,
) -> None:
    detector = _detector(log, tmp_path, recipeName=recipeName, noddingSequence=noddingSequence)

    assert _residual_pdf_name(detector, monkeypatch) == expectedName


@pytest.mark.parametrize(
    ("recipeName", "noddingSequence", "expectedName"),
    [
        ("soxs-stare", "", f"{SOF_NAME}_OBJTRACE.fits"),
        ("soxs-nod", "_A1", f"{SOF_NAME}_OBJTRACE_A1.fits"),
        ("soxs-mflat", "", "SOXS_VIS_OLOC_SOF.fits"),
        ("soxs-order-centres", "", f"{SOF_NAME}.fits"),
        ("soxs-nod-std", "_A1", f"{SOF_NAME}_OBJTRACE_A1.fits"),
        ("soxs-offset-std", "_B2", f"{SOF_NAME}_OBJTRACE_B2.fits"),
        ("soxs-offset", "_A1", f"{SOF_NAME}_OBJTRACE_A1.fits"),
        ("soxs-offset", "_B2", f"{SOF_NAME}_OBJTRACE_B2.fits"),
    ],
)
def test_order_table_name_follows_recipe_and_nodding_sequence(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    recipeName: str,
    noddingSequence: str,
    expectedName: str,
) -> None:
    detector = _detector(log, tmp_path, recipeName=recipeName, noddingSequence=noddingSequence)

    assert _order_table_name(detector, monkeypatch) == expectedName


def test_offset_sequences_write_distinct_residual_pdfs_and_order_tables(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    detectors = [
        _detector(log, tmp_path, recipeName="soxs-offset", noddingSequence=sequence) for sequence in ("_A1", "_B2")
    ]

    pdfNames = [_residual_pdf_name(detector, monkeypatch) for detector in detectors]
    tableNames = [_order_table_name(detector, monkeypatch) for detector in detectors]

    assert len(set(pdfNames)) == 2
    assert len(set(tableNames)) == 2
    assert f"{SOF_NAME}.fits" not in tableNames
    assert all("_A1" in name for name in (pdfNames[0], tableNames[0]))
    assert all("_B2" in name for name in (pdfNames[1], tableNames[1]))


def _fit_trace_stub(clipped: pd.DataFrame) -> object:
    """Return a stand-in for ``fit_global_polynomial`` that reports zero residuals."""

    def fit_trace(**kwargs: object) -> tuple[list[float], pd.DataFrame, pd.DataFrame]:
        fitted = kwargs["pixelList"].copy()
        fitted[f"{kwargs['axisACol']}_fit_res"] = 0.0
        return [16.0, 0.0], fitted, clipped

    return fit_trace


def test_offset_trace_residual_products_point_to_different_files(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    productsPerSequence = []
    monkeypatch.setattr("soxspipe.commonutils.toolkit.qc_settings_plot_tables", lambda **_: None)
    monkeypatch.setattr("soxspipe.commonutils.phase3.write_fits_table_to_disk", lambda **_: None)
    pixels, _, clipped = _plot_inputs()
    for sequence in ("_A1", "_B2"):
        detector = _detector(log, tmp_path, recipeName="soxs-offset", noddingSequence=sequence)
        detector.orderPixelTable = pixels.drop(columns=["cont_x_fit_res"]).assign(mask=False)
        detector.coeff_dict = {}
        monkeypatch.setattr(detector, "fit_global_polynomial", _fit_trace_stub(clipped))
        _, _, products, *_ = detector.get()
        productsPerSequence.append(products)

    residualRows = pd.concat(productsPerSequence)
    residualRows = residualRows[residualRows["product_label"].str.startswith("OBJECT_TRACE_RES")]
    assert residualRows["product_label"].tolist() == ["OBJECT_TRACE_RES_A1", "OBJECT_TRACE_RES_B2"]
    assert residualRows["file_path"].nunique() == 2


OFFSET_HEADER = {"HIERARCH ESO SEQ FIXOFF RA": -2.0, "HIERARCH ESO SEQ FIXOFF DEC": 0.0}


def _constructed_detector(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    recipeName: str,
    headerOverrides: dict[str, object],
    locationSetIndex: object = False,
) -> detect_continuum:
    """Build a continuum detector through its real constructor."""
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.utility_setup",
        lambda **_: (str(tmp_path / QC_DIR_NAME), str(tmp_path / PRODUCT_DIR_NAME)),
    )
    monkeypatch.setattr(
        importlib.import_module("soxspipe.commonutils.detect_continuum"), "get_calibration_lamp", lambda **_: ""
    )
    (tmp_path / QC_DIR_NAME).mkdir(exist_ok=True)
    frame = CCDData(
        np.full((32, 32), 20.0),
        unit=u.electron,
        meta=instrument_header(arm="VIS", overrides={"EXPTIME": 60.0, **headerOverrides}),
    )
    detector = detect_continuum(
        log=log,
        traceFrame=frame,
        dispersion_map="dispersion.fits",
        settings={**pipeline_settings(tmp_path), "tune-pipeline": False},
        recipeSettings={"detect-continuum": {"disp-axis-deg": 1, "order-deg": 0}},
        recipeName=recipeName,
        qcTable=qc_table(),
        productsTable=pd.DataFrame(),
        sofName=SOF_NAME,
        locationSetIndex=locationSetIndex,
        startNightDate="2024-01-02",
    )
    detector.detectorParams = {"rotate-qc-plot": 0, "flip-qc-plot": 0, "dispersion-axis": "x"}
    detector.axisBDeg = 1
    detector.orderDeg = 0
    return detector


@pytest.mark.parametrize(
    ("recipeName", "headerOverrides", "locationSetIndex", "expectedSequence"),
    [
        ("soxs-nod", {"HIERARCH ESO SEQ CUMOFF Y": 2.0}, 2, "_A2"),
        ("soxs-nod", {"HIERARCH ESO SEQ CUMOFF Y": -2.0}, 1, "_B1"),
        ("soxs-nod", {}, 1, ""),
        ("soxs-stare", {"HIERARCH ESO SEQ CUMOFF Y": 2.0}, False, "_A"),
        ("soxs-stare", {}, False, ""),
    ],
)
def test_nodding_sequence_for_nod_and_stare_is_unchanged(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    recipeName: str,
    headerOverrides: dict[str, object],
    locationSetIndex: object,
    expectedSequence: str,
) -> None:
    detector = _constructed_detector(
        log,
        monkeypatch,
        tmp_path,
        recipeName=recipeName,
        headerOverrides=headerOverrides,
        locationSetIndex=locationSetIndex,
    )

    assert detector.noddingSequence == expectedSequence


@pytest.mark.parametrize(
    ("locationSetIndex", "expectedSequence"),
    [(1, "_1"), (2, "_2"), (STACKED_LOCATION_SET, "_STACK")],
)
def test_offset_frames_without_cumoff_y_take_their_sequence_from_the_location_set(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    locationSetIndex: object,
    expectedSequence: str,
) -> None:
    detector = _constructed_detector(
        log,
        monkeypatch,
        tmp_path,
        recipeName="soxs-offset",
        headerOverrides=OFFSET_HEADER,
        locationSetIndex=locationSetIndex,
    )

    assert detector.noddingSequence == expectedSequence


def test_two_offset_cycles_and_the_stacked_pair_write_three_sets_of_trace_files(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    detectors = [
        _constructed_detector(
            log,
            monkeypatch,
            tmp_path,
            recipeName="soxs-offset",
            headerOverrides=OFFSET_HEADER,
            locationSetIndex=locationSetIndex,
        )
        for locationSetIndex in (1, 2, STACKED_LOCATION_SET)
    ]

    pdfNames = [_residual_pdf_name(detector, monkeypatch) for detector in detectors]
    tableNames = [_order_table_name(detector, monkeypatch) for detector in detectors]

    assert len(set(pdfNames)) == 3
    assert len(set(tableNames)) == 3
    assert f"{SOF_NAME}.fits".upper() not in tableNames
