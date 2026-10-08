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
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

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
    pixels = pd.DataFrame({"order": [10], "xcoord_centre": [1.0], "ycoord": [2.0]})

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
    order = pd.Series({"order": 10, "xcoord_centre": 50.0, "ycoord": 30.0})

    result = detector.determine_lower_upper_edge_pixel_positions(order)

    assert result["xcoord_lower"] == pytest.approx(40.25)
    assert result["xcoord_upper"] == pytest.approx(58.75)


def _get_detector(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    orderPixels: pd.DataFrame,
) -> detect_order_edges:
    """Return a detector with everything `get()` needs except the polynomial fit."""
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
    orderMeta = pd.DataFrame({"order": orderPixels["order"].unique()})
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
    detector.write_order_table_to_file = lambda **kwargs: str(tmp_path / "edges.fits")
    detector.plot_results = lambda **kwargs: str(tmp_path / "edges.pdf")
    return detector


def test_get_records_edge_fits_qc_and_product_contracts(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    orderPixels = pd.DataFrame(
        {
            "order": [10, 10],
            "xcoord_centre": [5.0, 5.0],
            "ycoord": [1.0, 3.0],
        }
    )
    detector = _get_detector(log, monkeypatch, tmp_path, orderPixels)

    def fit_edge(
        *,
        pixelList: pd.DataFrame,
        axisACol: str,
        **kwargs: object,
    ) -> tuple[list[float], pd.DataFrame, pd.DataFrame]:
        fitted = pixelList.copy()
        fitted["x_fit"] = fitted[axisACol]
        fitted["x_fit_res"] = [-0.123456789, 0.123456789]
        return [3.0], fitted, fitted.iloc[:0].copy()

    detector.fit_global_polynomial = fit_edge

    products, qc, detectionCounts = detector.get()

    assert detectionCounts.loc[10, "count"] == 2
    assert qc["qc_name"].tolist() == ["X RES MIN", "X RES MAX", "X RES SD"]
    assert qc["qc_value"].tolist() == [0.12346, 0.12346, 0.0]
    assert pd.api.types.is_numeric_dtype(qc["qc_value"])
    assert products["product_label"].tolist() == ["ORDER_LOC", "ORDER_LOC_RES"]
    assert products["file_name"].tolist() == ["edges.fits", "edges.pdf"]


def test_get_fits_lower_edge_using_only_rows_with_finite_lower_positions(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    # ARRANGE: ORDER 11 HAS A MISSING CENTRE, SO ITS EDGE POSITIONS COME OUT AS NAN
    orderPixels = pd.DataFrame(
        {
            "order": [10, 10, 11, 11],
            "xcoord_centre": [5.0, 5.0, 5.0, np.nan],
            "ycoord": [1.0, 3.0, 1.0, 3.0],
        }
    )
    detector = _get_detector(log, monkeypatch, tmp_path, orderPixels)
    fitInputs: dict[str, pd.DataFrame] = {}

    def fit_edge(
        *,
        pixelList: pd.DataFrame,
        axisACol: str,
        **kwargs: object,
    ) -> tuple[list[float], pd.DataFrame, pd.DataFrame]:
        if pixelList[axisACol].isna().any():
            raise ValueError("non-finite values in fit input")
        fitInputs[axisACol] = pixelList
        fitted = pixelList.copy()
        fitted["x_fit"] = fitted[axisACol]
        fitted["x_fit_res"] = 0.0
        return [3.0], fitted, fitted.iloc[:0].copy()

    detector.fit_global_polynomial = fit_edge

    # ACT
    detector.get()

    # ASSERT
    assert fitInputs["xcoord_lower"]["ycoord"].tolist() == [1.0, 3.0]
    assert fitInputs["xcoord_upper"]["ycoord"].tolist() == [1.0, 3.0]


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


def _plot_inputs(
    *,
    axisA: str,
    lowerEdge: float,
    upperEdge: float,
    axisBLength: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return constant-edge pixel tables, coefficients and order metadata for ``plot_results``."""
    axisB = "y" if axisA == "x" else "x"
    axisBValues = np.arange(3.0, axisBLength - 2.0, 3.0)
    lower = pd.DataFrame(
        {
            "order": np.full(len(axisBValues), 10.0),
            f"{axisA}coord_lower": np.full(len(axisBValues), lowerEdge),
            f"{axisA}coord_lower_fit_res": np.zeros(len(axisBValues)),
            f"{axisB}coord": axisBValues,
        }
    )
    upper = pd.DataFrame(
        {
            "order": np.full(len(axisBValues), 10.0),
            f"{axisA}coord_upper": np.full(len(axisBValues), upperEdge),
            f"{axisA}coord_upper_fit_res": np.zeros(len(axisBValues)),
            f"{axisB}coord": axisBValues,
        }
    )
    coefficients = pd.DataFrame([{"edgeup_0": upperEdge, "edgelow_0": lowerEdge}])
    metadata = pd.DataFrame({"order": [10], f"{axisB}min": [0.0], f"{axisB}max": [axisBLength - 1.0]})
    return upper, lower, coefficients, metadata


def _capture_saved_figure(monkeypatch: pytest.MonkeyPatch) -> dict[str, Figure]:
    """Record the figure that ``plot_results`` saves, while still writing the PDF."""
    import matplotlib.pyplot as plt

    captured: dict[str, Figure] = {}
    realSavefig = plt.savefig

    def capture(*args: object, **kwargs: object) -> None:
        captured["figure"] = plt.gcf()
        realSavefig(*args, **kwargs)

    monkeypatch.setattr(plt, "savefig", capture)
    return captured


def _plot_detector(
    log: object, tmp_path: Path, *, shape: tuple[int, int], rotate: int, flip: int
) -> detect_order_edges:
    """Return a detector ready for ``plot_results`` on a flat of ``shape``."""
    detector = _edge_detector(log)
    detector.axisAbin = 1
    detector.axisBbin = 1
    detector.axisBDeg = 0
    detector.orderDeg = 0
    detector.exptime = 2.0
    detector.slit = "5.0x11"
    detector.tag = "QLAMP"
    detector.sofName = "synthetic"
    detector.qcDir = str(tmp_path)
    detector.recipeSettings = {}
    detector.qc = qc_table()
    detector.detectorParams = {"rotate-qc-plot": rotate, "flip-qc-plot": flip}
    detector.flatFrame = CCDData(np.full(shape, 20.0), unit=u.electron)
    return detector


def test_plot_results_embeds_image_without_resampling(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Embed the flat at native resolution so vector edge markers line up with its pixels when zoomed."""
    detector = _plot_detector(log, tmp_path, shape=(32, 32), rotate=0, flip=0)
    upper, lower, coefficients, metadata = _plot_inputs(axisA="x", lowerEdge=8.0, upperEdge=20.0, axisBLength=32)
    captured = _capture_saved_figure(monkeypatch)

    detector.plot_results(upper, lower, coefficients, metadata, upper.iloc[:1], lower.iloc[:1])

    images = [image for ax in captured["figure"].axes for image in ax.get_images()]
    assert len(images) == 2
    assert [image.get_interpolation() for image in images] == ["none", "none"]


@pytest.mark.parametrize("axisAbin", [1, 2])
def test_flip_only_edges_map_to_flipped_pixel_rows(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    axisAbin: int,
) -> None:
    """A flip without rotation sends binned row ``r`` to row ``N - 1 - r``, for markers and fit lines alike."""
    rows = 64
    detector = _plot_detector(log, tmp_path, shape=(rows, 40), rotate=0, flip=1)
    detector.axisA = "y"
    detector.axisB = "x"
    detector.arm = "NIR"
    detector.axisAbin = axisAbin
    upper, lower, coefficients, metadata = _plot_inputs(
        axisA="y", lowerEdge=20.0 * axisAbin, upperEdge=40.0 * axisAbin, axisBLength=40
    )
    captured = _capture_saved_figure(monkeypatch)

    detector.plot_results(upper, lower, coefficients, metadata, upper.iloc[:1], lower.iloc[:1])

    topAxis, fitAxis = captured["figure"].axes[:2]
    markerRows = np.unique(topAxis.collections[0].get_offsets()[:, 1])
    assert markerRows.tolist() == [rows - 1 - 40.0, rows - 1 - 20.0]
    fitRows = sorted({float(np.unique(line.get_ydata())[0]) for line in fitAxis.get_lines()})
    assert fitRows == [rows - 1 - 40.0, rows - 1 - 20.0]


def test_plot_results_tables_do_not_overlap(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Keep the QC table and the mflat settings table apart on a wide VIS-shaped flat."""
    detector = _plot_detector(log, tmp_path, shape=(410, 82), rotate=90, flip=1)
    detector.recipeSettings = {
        "subtract_background": True,
        "stacked-clipping-sigma": 15,
        "stacked-clipping-iterations": 3,
        "centre-order-window": 10,
        "slice-length-for-edge-detection": 70,
        "slice-width-for-edge-detection": 5,
        "min-percentage-threshold-for-edge-detection": 20,
        "max-percentage-threshold-for-edge-detection": 50,
        "disp-axis-deg": 3,
        "order-deg": 3,
        "poly-fitting-residual-clipping-sigma": 5,
        "poly-clipping-iteration-limit": 2,
        "low-sensitivity-clipping-sigma": 2,
        "order-edge-min-flat-fraction": 0.5,
    }
    detector.qc = pd.concat(
        [
            qc_table().assign(qc_name=name, qc_value=value, qc_unit=unit, qc_comment=comment)
            for name, value, unit, comment in [
                ("ORDEXP10", 234.3, "electrons", "[e-] 10th percentile inter-order flux"),
                ("ORDEXP50", 5457.21, "electrons", "[e-] 50th percentile inter-order flux"),
                ("ORDEXP90", 40640.352, "electrons", "[e-] 90th percentile inter-order flux"),
                ("X RES MIN", 0.00206, "pixels", "[px] Minimum residual in order edge fit along x-axis"),
                ("X RES MAX", 2.24219, "pixels", "[px] Maximum residual in order edge fit along x-axis"),
                ("X RES SD", 0.3453, "pixels", "[px] Std-dev of residual order edge fit along x-axis"),
            ]
        ],
        ignore_index=True,
    )
    upper, lower, coefficients, metadata = _plot_inputs(axisA="x", lowerEdge=30.0, upperEdge=50.0, axisBLength=410)
    captured = _capture_saved_figure(monkeypatch)

    detector.plot_results(upper, lower, coefficients, metadata, upper.iloc[:1], lower.iloc[:1])

    figure = captured["figure"]
    # PLOT_RESULTS CLOSES THE FIGURE, SO ATTACH AN AGG CANVAS RATHER THAN RELY ON THE ACTIVE BACKEND
    renderer = FigureCanvasAgg(figure).get_renderer()
    figure.draw(renderer)
    tables = [table for ax in figure.axes for table in ax.tables]
    assert len(tables) == 2
    qcExtent, settingsExtent = (table.get_window_extent(renderer) for table in tables)
    assert not qcExtent.overlaps(settingsExtent)


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
