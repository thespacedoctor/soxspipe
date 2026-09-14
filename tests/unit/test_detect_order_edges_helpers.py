"""Analytic contracts for order-edge flux and position measurements."""

from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe.commonutils.detect_order_edges import detect_order_edges
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from tests.factories import instrument_header, pipeline_settings, qc_table

pytestmark = pytest.mark.unit
orderEdgesModule = importlib.import_module("soxspipe.commonutils.detect_order_edges")


def _edge_detector(log: object) -> detect_order_edges:
    """Return the minimal detector state needed by the measurement methods."""
    detector = detect_order_edges.__new__(detect_order_edges)
    detector.log = log
    detector.axisA = "x"
    detector.axisB = "y"
    detector.arm = "VIS"
    detector.flatFrame = object()
    detector.sliceWidth = 3
    detector.sliceLength = 40
    detector.minThresholdPercenage = 0.25
    detector.maxThresholdPercenage = 0.75
    detector.debug = False
    return detector


def test_flux_threshold_uses_central_order_cross_section(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _edge_detector(log)
    captured: dict[str, float] = {}
    crossSection = np.array([0.0] * 5 + [10.0] * 21 + [0.0] * 5)
    monkeypatch.setattr(
        orderEdgesModule,
        "cut_image_slice",
        lambda **kwargs: (captured.update(kwargs) or crossSection, 0, 0),
    )
    pixels = pd.DataFrame(
        {
            "order": [10, 10, 10],
            "xcoord_centre": [100.0, 101.0, 102.0],
            "ycoord": [50.0, 51.0, 52.0],
        }
    )

    result = detector.determine_order_flux_threshold(
        pd.Series({"order": 10}),
        pixels,
    )

    assert (captured["x"], captured["y"], captured["sliceAxis"]) == (
        101.0,
        51.0,
        "x",
    )
    assert result["minThreshold"] == pytest.approx(2.5)
    assert result["maxThreshold"] == pytest.approx(7.5)
    assert result["maxvalue"] == pytest.approx(10.0)


def test_flux_threshold_leaves_order_unchanged_when_slice_is_out_of_bounds(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _edge_detector(log)
    monkeypatch.setattr(
        orderEdgesModule,
        "cut_image_slice",
        lambda **kwargs: (None, None, None),
    )
    order = pd.Series({"order": 10})
    pixels = pd.DataFrame(
        {"order": [10], "xcoord_centre": [1.0], "ycoord": [2.0]}
    )

    result = detector.determine_order_flux_threshold(order, pixels)

    pd.testing.assert_series_equal(result, order)


def test_edge_positions_interpolate_threshold_crossings(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _edge_detector(log)
    crossSection = np.array([0.0] * 10 + [10.0] * 20 + [0.0] * 10)
    monkeypatch.setattr(
        orderEdgesModule,
        "cut_image_slice",
        lambda **kwargs: (crossSection, 0, 0),
    )
    order = pd.Series(
        {"order": 10, "xcoord_centre": 50.0, "ycoord": 30.0}
    )

    result = detector.determine_lower_upper_edge_pixel_positions(order)

    assert result["xcoord_lower"] == pytest.approx(40.25)
    assert result["xcoord_upper"] == pytest.approx(58.75)


def test_get_records_edge_fits_qc_and_product_contracts(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    detector = _edge_detector(log)
    detector.recipeSettings = {
        "slice-length-for-edge-detection": 20,
        "slice-width-for-edge-detection": 3,
        "min-percentage-threshold-for-edge-detection": 25,
        "max-percentage-threshold-for-edge-detection": 75,
    }
    detector.orderCentreTable = "centres.fits"
    detector.binx = 1
    detector.biny = 1
    detector.pixelDelta = 1
    detector.slit = ""
    detector.axisAbin = 1
    detector.axisBbin = 1
    detector.axisBDeg = 0
    detector.orderDeg = 0
    detector.extendToEdges = True
    detector.flatFrame = SimpleNamespace(data=np.zeros((5, 10)))
    detector.qc = pd.DataFrame()
    detector.products = pd.DataFrame()
    detector.recipeName = "soxs-mflat"
    detector.dateObs = "2024-01-01T00:00:00"
    detector.tag = ""
    orderPixels = pd.DataFrame(
        {
            "order": [10, 10],
            "xcoord_centre": [5.0, 5.0],
            "ycoord": [1.0, 3.0],
        }
    )
    orderMeta = pd.DataFrame({"order": [10]})
    monkeypatch.setattr(
        orderEdgesModule,
        "unpack_order_table",
        lambda **kwargs: (pd.DataFrame({"existing": [1]}), orderPixels, orderMeta),
    )

    def with_thresholds(row: pd.Series, **kwargs: object) -> pd.Series:
        return pd.Series({**row.to_dict(), "minThreshold": 2.5, "maxThreshold": 7.5})

    def with_edge_positions(row: pd.Series) -> pd.Series:
        return pd.Series({**row.to_dict(), "xcoord_upper": 8.0, "xcoord_lower": 2.0})

    detector.determine_order_flux_threshold = with_thresholds
    detector.determine_lower_upper_edge_pixel_positions = with_edge_positions

    def fit_edge(
        *,
        pixelList: pd.DataFrame,
        axisACol: str,
        **kwargs: object,
    ) -> tuple[list[float], pd.DataFrame, pd.DataFrame]:
        fitted = pixelList.copy()
        fitted["x_fit"] = fitted[axisACol]
        fitted["x_fit_res"] = [-0.25, 0.25]
        return [3.0], fitted, fitted.iloc[:0].copy()

    detector.fit_global_polynomial = fit_edge
    detector.write_order_table_to_file = lambda **kwargs: str(tmp_path / "edges.fits")
    detector.plot_results = lambda **kwargs: str(tmp_path / "edges.pdf")

    products, qc, detectionCounts = detector.get()

    assert detectionCounts.loc[10, "count"] == 2
    assert qc["qc_name"].tolist() == ["X RES MIN", "X RES MAX", "X RES SD"]
    assert products["product_label"].tolist() == ["ORDER_LOC", "ORDER_LOC_RES"]
    assert products["file_name"].tolist() == ["edges.fits", "edges.pdf"]


def test_plot_results_writes_order_edge_diagnostic(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Render the public edge-fit diagnostic from compact analytic data."""
    detector = _edge_detector(log)
    detector.axisAbin = 1
    detector.axisBbin = 1
    detector.axisBDeg = 0
    detector.orderDeg = 0
    detector.exptime = 1.0
    detector.slit = ""
    detector.tag = ""
    detector.sofName = "synthetic"
    detector.qcDir = str(tmp_path)
    detector.recipeSettings = {}
    detector.qc = qc_table()
    detector.detectorParams = {"rotate-qc-plot": False, "flip-qc-plot": False}
    detector.flatFrame = CCDData(np.full((32, 32), 20.0), unit=u.electron)
    yValues = np.arange(3.0, 30.0, 3.0)
    lower = pd.DataFrame(
        {
            "order": np.full(len(yValues), 10.0),
            "xcoord_lower": np.full(len(yValues), 8.0),
            "xcoord_lower_fit_res": np.zeros(len(yValues)),
            "ycoord": yValues,
        }
    )
    upper = pd.DataFrame(
        {
            "order": np.full(len(yValues), 10.0),
            "xcoord_upper": np.full(len(yValues), 20.0),
            "xcoord_upper_fit_res": np.zeros(len(yValues)),
            "ycoord": yValues,
        }
    )
    coefficients = pd.DataFrame([{"edgeup_0": 20.0, "edgelow_0": 8.0}])
    metadata = pd.DataFrame({"order": [10], "ymin": [0.0], "ymax": [31.0]})
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.qc_settings_plot_tables",
        lambda **_: None,
    )

    outputPath = detector.plot_results(
        upper,
        lower,
        coefficients,
        metadata,
        upper.iloc[:1],
        lower.iloc[:1],
    )

    assert Path(outputPath).read_bytes().startswith(b"%PDF")
    assert Path(outputPath).name == "synthetic_ORD_LOC.pdf"


def test_constructor_resolves_vis_metadata_and_qc_directories(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Construct the public detector with its FITS metadata contract intact."""
    qcPath = tmp_path / "qc"
    productPath = tmp_path / "products"
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.utility_setup",
        lambda **_: (str(qcPath), str(productPath)),
    )
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.quicklook_image",
        lambda **_: None,
    )
    header = instrument_header(
        arm="VIS",
        overrides={
            "EXPTIME": 60.0,
            "ESO INS VISE NAME": "1.0x11",
        },
    )
    frame = CCDData(np.ones((32, 32)), unit=u.electron, meta=header)

    detector = detect_order_edges(
        log=log,
        flatFrame=frame,
        orderCentreTable="centres.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"disp-axis-deg": 2, "order-deg": 1},
        startNightDate="2024-01-02",
    )

    assert detector.arm == "VIS"
    assert detector.slit == "1.0x11"
    assert detector.axisA == "x"
    assert detector.axisB == "y"
    assert detector.pixelDelta == 25
    assert (detector.qcDir, detector.productDir) == (str(qcPath), str(productPath))


def test_write_order_table_preserves_mflat_fits_contract(
    log: object,
    tmp_path: Path,
) -> None:
    """Write compact order data through the public FITS product boundary."""
    detector = _edge_detector(log)
    detector.kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get
    detector.settings = pipeline_settings(tmp_path)
    detector.productDir = str(tmp_path)
    detector.sofName = "synthetic_flat"
    detector.recipeName = "soxs-mflat"
    detector.lampTag = False
    detector.inst = "SOXS"
    detector.qc = qc_table()
    header = instrument_header(arm="VIS")
    frame = CCDData(np.ones((8, 8)), unit=u.electron, meta=header)

    outputPath = detector.write_order_table_to_file(
        frame,
        pd.DataFrame({"edgeup_0": [20.0], "edgelow_0": [8.0]}),
        pd.DataFrame({"order": [10], "ymin": [0.0], "ymax": [7.0]}),
    )

    with fits.open(outputPath) as product:
        assert Path(outputPath).name == "SYNTHETIC_OLOC.fits"
        assert product[0].header["ESO SEQ ARM"] == "VIS"
        assert product[0].header["HIERARCH ESO PRO CATG"] == "ORDER_TAB_VIS"
        assert product[1].data["edgeup_0"][0] == pytest.approx(20.0)
        assert product[2].data["order"][0] == 10
