"""Analytic continuum-trace fitting contracts."""

from __future__ import annotations

import importlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe.commonutils.detect_continuum import detect_continuum
from tests.factories import instrument_header, pipeline_settings, qc_table
from tests.unit._plot_spies import image_interpolations, probe_on_savefig, table_extents

pytestmark = pytest.mark.unit
continuumModule = importlib.import_module("soxspipe.commonutils.detect_continuum")


def _detector(log: object, *, orderDeg: int = 0, axisBDeg: int = 1) -> detect_continuum:
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.arm = "VIS"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.orderDeg = orderDeg
    detector.axisBDeg = axisBDeg
    detector.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": 2.0,
        "poly-clipping-iteration-limit": 4,
    }
    return detector


def test_calculate_residuals_matches_exact_global_trace(log: object) -> None:
    detector = _detector(log, orderDeg=1, axisBDeg=1)
    pixels = pd.DataFrame(
        {
            "order": [10.0, 10.0, 11.0],
            "cont_y": [0.0, 1.0, 2.0],
        }
    )
    expectedFit = 1.0 + 2.0 * pixels["cont_y"] + 3.0 * pixels["order"] + 4.0 * pixels["order"] * pixels["cont_y"]
    pixels["cont_x"] = expectedFit - np.array([1.0, -2.0, 0.0])

    residuals, mean, std, median, fitted = detector.calculate_residuals(
        orderPixelTable=pixels,
        coeff=[1.0, 2.0, 3.0, 4.0],
        axisACol="cont_x",
        axisBCol="cont_y",
        orderCol="order",
    )

    np.testing.assert_allclose(fitted, expectedFit, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(residuals, [1.0, -2.0, 0.0], rtol=1e-12, atol=1e-12)
    assert mean == pytest.approx(1.0)
    assert std == pytest.approx(np.std([1.0, -2.0, 0.0]))
    assert median == pytest.approx(1.0)


def test_calculate_residuals_records_numeric_rounded_qcs(log: object) -> None:
    detector = _detector(log, orderDeg=0, axisBDeg=0)
    detector.recipeName = "soxs-order-centres"
    detector.dateObs = "2024-01-01T00:00:00"
    detector.qc = pd.DataFrame()
    pixels = pd.DataFrame({"order": [10.0] * 3, "cont_y": [0.0, 1.0, 2.0]})
    offsets = np.array([0.1234567, -0.2, 0.0])
    pixels["cont_x"] = 5.0 - offsets

    detector.calculate_residuals(
        orderPixelTable=pixels,
        coeff=[5.0],
        axisACol="cont_x",
        axisBCol="cont_y",
        orderCol="order",
        writeQCs=True,
    )

    resQc = detector.qc[detector.qc["qc_name"].str.contains("RES")]
    assert resQc["qc_name"].tolist() == ["X RES MIN", "X RES MAX", "X RES SD", "X RES MEDIAN"]
    assert resQc["qc_value"].tolist() == [-0.2, 0.123, round(float(np.std(offsets)), 3), 0.108]
    assert pd.api.types.is_numeric_dtype(detector.qc["qc_value"])
    assert resQc["to_header"].all()


def test_global_fit_removes_nan_and_outlier_rows(log: object) -> None:
    detector = _detector(log)
    yValues = np.arange(12, dtype=float)
    pixels = pd.DataFrame(
        {
            "source": np.arange(12),
            "order": np.full(12, 10.0),
            "cont_y": yValues,
            "cont_x": 1.0 + 2.0 * yValues,
        }
    )
    pixels.loc[10, "cont_x"] = 100.0
    pixels.loc[11, "cont_y"] = np.nan

    coefficients, fitted, clipped = detector.fit_global_polynomial(pixels.copy())

    np.testing.assert_allclose(coefficients, [1.0, 2.0], rtol=1e-6, atol=1e-6)
    assert {10, 11}.issubset(set(clipped["source"]))
    assert {10, 11}.isdisjoint(set(fitted["source"]))
    np.testing.assert_allclose(fitted["cont_x_fit_res"], 0.0, rtol=0, atol=1e-6)


def test_global_fit_preserves_empty_input_failure(log: object) -> None:
    detector = _detector(log)
    empty = pd.DataFrame(columns=["order", "cont_y", "cont_x"])

    with pytest.raises(ValueError, match="cannot set a frame with no defined index"):
        detector.fit_global_polynomial(empty)


def test_get_stops_before_fitting_when_trace_sampling_is_insufficient(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The public continuum seam preserves its no-solution result for sparse samples."""
    detector = _detector(log)
    detector.orderPixelTable = False
    detector.coeff_dict = {}
    detector.qc = pd.DataFrame({"qc_name": []})
    detector.products = pd.DataFrame({"product_label": []})
    monkeypatch.setattr(
        detector,
        "sample_trace",
        lambda: (pd.DataFrame(), 9.0),
    )

    result = detector.get()

    assert result == (None, detector.qc, detector.products, None, None, None)


def test_get_returns_no_solution_after_exhausting_polynomial_retries(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A repeatedly failing global fit keeps the established no-solution contract."""
    detector = _detector(log, orderDeg=2, axisBDeg=2)
    detector.orderPixelTable = pd.DataFrame(
        {
            "order": [10.0, 10.0],
            "cont_x": [1.0, 2.0],
            "cont_y": [3.0, 4.0],
            "mask": [False, False],
        }
    )
    detector.coeff_dict = {}
    detector.recipeName = "soxs-order-centres"
    detector.inst = "SOXS"
    detector.dateObs = "2021-01-01T00:00:00"
    detector.settings = {"tune-pipeline": False}
    detector.qc = pd.DataFrame({"qc_name": []})
    detector.products = pd.DataFrame({"product_label": []})
    monkeypatch.setattr(
        detector,
        "fit_global_polynomial",
        lambda **kwargs: (_ for _ in ()).throw(AttributeError("no fit")),
    )

    result = detector.get()

    assert result == (None, detector.qc, detector.products, None, None, None)
    assert detector.axisBDeg == -1
    assert detector.orderDeg == 0
    assert detector.recipeSettings["disp-axis-deg"] == -1
    assert detector.recipeSettings["order-deg"] == 0


def test_get_records_fitted_trace_qc_product_and_order_table(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A converged public trace fit preserves its QC and product contracts."""
    detector = _detector(log, orderDeg=1, axisBDeg=1)
    detector.orderPixelTable = pd.DataFrame(
        {
            "order": [10.0, 11.0],
            "cont_x": [2.0, 3.0],
            "cont_y": [4.0, 5.0],
            "gauss_stddev": [1.0, 1.5],
            "mask": [False, False],
        }
    )
    detector.coeff_dict = {}
    detector.recipeName = "soxs-order-centres"
    detector.inst = "SOXS"
    detector.dateObs = "2024-01-02T00:00:00"
    detector.settings = {"tune-pipeline": False}
    detector.qc = pd.DataFrame({"qc_name": []})
    detector.products = pd.DataFrame({"product_label": []})
    detector.noddingSequence = ""
    detector.traceFrame = object()
    fittedCalls: list[str] = []

    def fit_trace(**kwargs: object) -> tuple[list[float], pd.DataFrame, pd.DataFrame]:
        axisACol = str(kwargs["axisACol"])
        fittedCalls.append(axisACol)
        fitted = kwargs["pixelList"].copy()
        fitted[f"{axisACol}_fit_res"] = [0.0, 0.5]
        clipped = pd.DataFrame({"order": [10.0]}) if axisACol == "cont_x" else pd.DataFrame()
        return [1.0, 2.0, 3.0, 4.0], fitted, clipped

    monkeypatch.setattr(detector, "fit_global_polynomial", fit_trace)
    monkeypatch.setattr(
        detector,
        "plot_results",
        lambda **kwargs: "/qc/continuum.pdf",
    )
    monkeypatch.setattr(detector, "build_order_metadata", lambda **kwargs: pd.DataFrame({"axis": ["x"]}))
    writeCalls: list[dict] = []
    monkeypatch.setattr(
        detector,
        "write_order_table_to_file",
        lambda **kwargs: writeCalls.append(kwargs) or "/products/orders.fits",
    )

    result = detector.get()

    orderPath, qc, products, orderPoly, orderPixels, orderMeta = result
    assert orderPath == "/products/orders.fits"
    assert fittedCalls == ["cont_x", "gauss_stddev"]
    assert qc["qc_name"].tolist() == ["SAMPLES CLIP NUM", "SAMPLES CLIP FRAC"]
    assert qc["qc_value"].tolist() == [1, 0.5]
    assert pd.api.types.is_numeric_dtype(qc["qc_value"])
    assert products["product_label"].tolist() == ["ORDER_CENTRES_RES"]
    assert products["file_name"].tolist() == ["continuum.pdf"]
    assert orderPoly.loc[0, "cent_11"] == 4.0
    assert orderPoly.loc[0, "std_11"] == 4.0
    assert "mask" not in orderPixels
    assert orderMeta.to_dict("list") == {"axis": ["x"]}
    assert writeCalls[0]["orderMetaTable"] is orderMeta


def test_gaussian_slice_fit_recovers_trace_and_leaves_low_signal_unfitted(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _detector(log)
    detector.traceFrame = object()
    detector.sliceWidth = 3
    detector.detectorParams = {"dispersion-axis": "x"}
    detector.debug = False
    detector.peakSigmaLimit = 3.0
    detector.sliceAxis = "x"
    detector.sliceAntiAxis = "y"
    samplePositions = np.arange(21, dtype=float)
    gaussian = 20.0 * np.exp(-0.5 * ((samplePositions - 10.0) / 2.0) ** 2)
    slices = iter([(gaussian, 5, 7.0), (np.zeros(21), 5, 8.0)])
    monkeypatch.setattr(
        continuumModule,
        "cut_image_slice",
        lambda **kwargs: next(slices),
    )
    pixels = pd.DataFrame({"fit_x": [15.0, 16.0], "fit_y": [7.0, 8.0]})

    result = detector.fit_1d_gaussian_to_slices(pixels.copy(), sliceLength=21)

    np.testing.assert_allclose(result.loc[0, "cont_x"], 15.0, rtol=0, atol=1e-6)
    np.testing.assert_allclose(result.loc[0, "cont_y"], 7.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.loc[0, "gauss_mean"], 10.0, rtol=0, atol=1e-6)
    np.testing.assert_allclose(result.loc[0, "gauss_stddev"], 2.0, rtol=0, atol=1e-6)
    assert np.isnan(result.loc[1, "cont_x"])
    assert np.isnan(result.loc[1, "gauss_mean"])


def test_gaussian_slice_fit_reports_the_anti_axis_centre_of_the_collapsed_rows(log: object) -> None:
    # ARRANGE: A HORIZONTAL TRACE ON ROW 30, SO THE ONE-ROW SLICE COLLAPSES ONLY ROW 30
    detector = _detector(log)
    detector.sliceWidth = 1
    detector.detectorParams = {"dispersion-axis": "x"}
    detector.debug = False
    detector.peakSigmaLimit = 3.0
    detector.sliceAxis = "x"
    detector.sliceAntiAxis = "y"
    columns = np.arange(80, dtype=float)
    frame = np.tile(20.0 * np.exp(-0.5 * ((columns - 40.0) / 2.0) ** 2), (60, 1))
    detector.traceFrame = np.ma.masked_array(frame, mask=np.zeros(frame.shape, dtype=bool))
    pixels = pd.DataFrame({"fit_x": [40.0], "fit_y": [30.0]})

    # ACT
    result = detector.fit_1d_gaussian_to_slices(pixels, sliceLength=21)

    # ASSERT: PIXEL CENTRES ARE AT INTEGER INDICES, SO ROW 30 IS REPORTED AS 30.0 AND NOT 30.5
    assert result.loc[0, "cont_y"] == pytest.approx(30.0, abs=1e-12)


def test_constructor_resolves_trace_orientation_and_nodding_metadata(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Construct a continuum detector from a real compact FITS header."""
    qcPath = tmp_path / "qc"
    productPath = tmp_path / "products"
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.utility_setup",
        lambda **_: (str(qcPath), str(productPath)),
    )
    monkeypatch.setattr(
        continuumModule,
        "get_calibration_lamp",
        lambda **_: "QTH",
    )
    header = instrument_header(
        arm="VIS",
        overrides={"EXPTIME": 60.0, "HIERARCH ESO SEQ CUMOFF Y": 2.0},
    )
    frame = CCDData(np.ones((32, 32)), unit=u.electron, meta=header)

    detector = detect_continuum(
        log=log,
        traceFrame=frame,
        dispersion_map="dispersion.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"detect-continuum": {"disp-axis-deg": 2, "order-deg": 1}},
        recipeName="soxs-order-centres",
        startNightDate="2024-01-02",
    )

    assert detector.noddingSequence == "_A"
    assert detector.lamp == "QTH"
    assert (detector.sliceAxis, detector.sliceAntiAxis) == ("x", "y")
    assert (detector.axisA, detector.axisB) == ("x", "y")
    assert detector.coeff_dict == {"degorder_cent": 1, "degy_cent": 2}
    assert (detector.qcDir, detector.productDir) == (str(qcPath), str(productPath))


def test_create_pixel_arrays_samples_each_spectral_order(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Sample spectral-format limits and retain dispersion-map binning metadata."""
    detector = _detector(log)
    detector.settings = {"instrument": "soxs"}
    detector.dispersion_map = str(tmp_path / "dispersion.fits")
    detector.recipeSettings["order-sample-count"] = 5
    detector.kw = lambda key: {"WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    fits.PrimaryHDU(header=fits.Header({"BINX": 2, "BINY": 3})).writeto(detector.dispersion_map)
    monkeypatch.setattr(
        continuumModule,
        "read_spectral_format",
        lambda **_: ([10, 11], [500.0, 600.0], [510.0, 620.0], [0.0, 10.0], [10.0, 30.0]),
    )
    monkeypatch.setattr(
        continuumModule,
        "dispersion_map_to_pixel_arrays",
        lambda *, orderPixelTable, **_: orderPixelTable.assign(
            fit_x=orderPixelTable["order"], fit_y=orderPixelTable["wavelength"]
        ),
    )

    pixels, binx, biny = detector.create_pixel_arrays()

    assert pixels.groupby("order").size().to_dict() == {10.0: 5, 11.0: 10}
    assert pixels["slit_position"].eq(0.0).all()
    assert (binx, biny) == (2, 3)


def test_sample_trace_records_detection_qc_for_analytic_trace(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The public sampler records a fully detected analytic object trace."""
    detector = _detector(log)
    detector.inst = "XSH"
    detector.recipeName = "soxs-stare"
    detector.traceFrame = CCDData(np.ones((32, 32)), unit=u.electron)
    detector.kw = lambda key: {"DPR_TYPE": "DPR TYPE", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    detector.recipeSettings.update({"slice-length": 12, "peak-sigma-limit": 3.0, "slice-width": 5})
    detector.debug = False
    detector.qc = pd.DataFrame()
    detector.dateObs = "2024-01-02T03:04:05"
    samples = pd.DataFrame(
        {
            "order": np.full(50, 10.0),
            "wavelength": np.linspace(500.0, 510.0, 50),
            "fit_x": np.full(50, 16.0),
            "fit_y": np.arange(4.0, 54.0),
        }
    )
    monkeypatch.setattr(detector, "create_pixel_arrays", lambda: (samples.copy(), 1, 1))
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.quicklook_image",
        lambda **_: None,
    )

    def fit_slices(
        orderPixelTable: pd.DataFrame,
        **_: object,
    ) -> pd.DataFrame:
        return orderPixelTable.assign(
            cont_x=orderPixelTable["fit_x"] - np.linspace(1.9, 2.1, len(orderPixelTable)),
            cont_y=orderPixelTable["fit_y"],
            gauss_stddev=1.5,
        )

    monkeypatch.setattr(detector, "fit_1d_gaussian_to_slices", fit_slices)

    result, detectionPercentage = detector.sample_trace()

    assert detectionPercentage == pytest.approx(100.0)
    assert result["pre-clipped"].eq(False).all()
    assert detector.qc["qc_name"].tolist() == [
        "SAMPLES TOT NUM",
        "SAMPLES DET NUM",
        "SAMPLES DET FRAC",
    ]
    assert detector.qc["qc_value"].tolist() == [50, 50, 1.0]
    assert pd.api.types.is_numeric_dtype(detector.qc["qc_value"])


def test_sample_trace_keeps_soxs_vis_order_groups_separate_until_fitted(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The SOXS VIS sampler fits its two order pairs before recombining them."""
    detector = _detector(log)
    detector.inst = "SOXS"
    detector.arm = "VIS"
    detector.recipeName = "soxs-mflat"
    detector.traceFrame = CCDData(np.ones((32, 32)), unit=u.electron)
    detector.kw = lambda key: {"DPR_TYPE": "DPR TYPE", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    detector.recipeSettings.update({"slice-length": 12, "peak-sigma-limit": 3.0, "slice-width": 5})
    detector.debug = False
    detector.qc = pd.DataFrame()
    detector.dateObs = "2024-01-02T03:04:05"
    orders = np.repeat([1.0, 2.0, 3.0, 4.0], 50)
    samples = pd.DataFrame(
        {
            "order": orders,
            "wavelength": np.linspace(500.0, 700.0, len(orders)),
            "fit_x": np.full(len(orders), 16.0),
            "fit_y": np.arange(len(orders), dtype=float),
        }
    )
    monkeypatch.setattr(detector, "create_pixel_arrays", lambda: (samples.copy(), 1, 1))
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None)
    monkeypatch.setattr(
        detector,
        "fit_1d_gaussian_to_slices",
        lambda orderPixelTable, **_: orderPixelTable.assign(
            cont_x=orderPixelTable["fit_x"] - np.linspace(1.9, 2.1, len(orderPixelTable)),
            cont_y=orderPixelTable["fit_y"],
            gauss_stddev=1.5,
        ),
    )

    result, detectionPercentage = detector.sample_trace()

    assert detectionPercentage == pytest.approx(100.0)
    assert result.groupby("order").size().to_dict() == {1.0: 50, 2.0: 50, 3.0: 50, 4.0: 50}


def _plot_ready_detector(tmp_path: Path, log: object, *, orderDeg: int = 0, axisBDeg: int = 1) -> detect_continuum:
    """Return a detector carrying just the state `plot_results` reads, on a 32x32 trace frame."""
    detector = _detector(log, orderDeg=orderDeg, axisBDeg=axisBDeg)
    detector.recipeName = "soxs-order-centres"
    detector.sofName = "synthetic"
    detector.qcDir = str(tmp_path)
    detector.settings = {"tune-pipeline": False}
    detector.exptime = 1.0
    detector.lamp = False
    detector.sliceLength = 4
    detector.debug = False
    detector.qc = pd.DataFrame({"qc_name": []})
    detector.detectorParams = {
        "dispersion-axis": "x",
        "flip-qc-plot": False,
        "rotate-qc-plot": False,
    }
    detector.traceFrame = CCDData(np.full((32, 32), 20.0), unit=u.electron)
    return detector


def _trace_pixels(orders: list[float], slope: float, intercept: float) -> pd.DataFrame:
    """Return a residual-free pixel table with one straight trace per order."""
    yValues = np.arange(4.0, 28.0)
    traceY = np.tile(yValues, len(orders))
    return pd.DataFrame(
        {
            "order": np.repeat(orders, len(yValues)),
            "cont_x": intercept + slope * traceY,
            "cont_y": traceY,
            "cont_x_fit_res": np.zeros(len(traceY)),
            "gauss_stddev_fit": np.full(len(traceY), 1.5),
            "wavelength": 500.0 + traceY,
        }
    )


def _no_continuum_clipped() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cont_x": [np.nan],
            "cont_y": [10.0],
            "fit_x": [10.0],
            "fit_y": [10.0],
        }
    )


@pytest.fixture
def savefig_text_colours(monkeypatch: pytest.MonkeyPatch):
    """Record `{label: colour}` for every axis of the current figure each time `plt.savefig` is called.

    `plot_results` clears every axis after saving, so the labels must be captured at save time.
    """
    snapshots = []
    realSavefig = plt.savefig

    def recorder(*args, **kwargs):
        snapshots.append([{t.get_text(): t.get_color() for t in ax.texts} for ax in plt.gcf().axes])
        return realSavefig(*args, **kwargs)

    monkeypatch.setattr(plt, "savefig", recorder)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.qc_settings_plot_tables", lambda **_: None)
    yield snapshots
    plt.close("all")


def test_plot_results_writes_continuum_diagnostic_and_order_limits(
    tmp_path, log: object, savefig_text_colours: list
) -> None:
    """Render the fitted trace diagnostic from a compact analytic order table."""
    # ARRANGE
    detector = _plot_ready_detector(tmp_path, log)
    pixels = _trace_pixels([10.0, 11.0], slope=0.2, intercept=8.0)
    coefficients = pd.DataFrame([{"cent_0": 8.0, "cent_1": 0.2, "std_0": 1.5, "std_1": 0.0}])

    # ACT
    outputPath = detector.plot_results(pixels, coefficients, _no_continuum_clipped())
    orderMetadata = detector.build_order_metadata(pixels, coefficients)

    # ASSERT
    assert Path(outputPath).read_bytes().startswith(b"%PDF")
    assert orderMetadata["order"].tolist() == [10.0, 11.0]
    assert (orderMetadata["xmax"] > orderMetadata["xmin"]).all()


def test_plot_results_pairs_each_order_with_its_own_colour_when_an_earlier_fit_is_off_detector(
    tmp_path, log: object, savefig_text_colours: list
) -> None:
    """Order 10 is fitted at x=-2 (off-detector) so the middle panel skips it; orders 11 and 12 keep their colours."""
    # ARRANGE
    detector = _plot_ready_detector(tmp_path, log, orderDeg=1, axisBDeg=1)
    pixels = _trace_pixels([10.0, 11.0, 12.0], slope=0.0, intercept=16.0)
    coefficients = pd.DataFrame(
        [
            {
                "cent_00": -82.0,
                "cent_01": 0.0,
                "cent_10": 8.0,
                "cent_11": 0.0,
                "std_00": 1.5,
                "std_01": 0.0,
                "std_10": 0.0,
                "std_11": 0.0,
            }
        ]
    )

    # ACT
    detector.plot_results(pixels, coefficients, _no_continuum_clipped())
    orderMetadata = detector.build_order_metadata(pixels, coefficients)

    # ASSERT
    midrow, bottomleft, bottomright, fwhmaxis = savefig_text_colours[0][1:5]
    assert orderMetadata["order"].tolist() == [11.0, 12.0]
    assert set(midrow) == {"11", "12"}
    for panel in (bottomleft, bottomright, fwhmaxis):
        assert {"10", "11", "12"} <= set(panel)
        assert panel["11"] == midrow["11"]
        assert panel["12"] == midrow["12"]
    assert len({midrow["11"], midrow["12"]}) == 2


# TWO ORDERS WHOSE CENTRES ARE -82 + 8 * ORDER, SO ORDERS 11 AND 12 SIT AT 6 AND 14 PIXELS, AND 12 AND 13 AT 14 AND 22
FLIP_COEFFICIENTS = {
    "cent_00": -82.0,
    "cent_01": 0.0,
    "cent_10": 8.0,
    "cent_11": 0.0,
    "std_00": 1.5,
    "std_01": 0.0,
    "std_10": 0.0,
    "std_11": 0.0,
}
FLIP_ROWS = 64
FLIP_COLUMNS = 40


def _nir_flip_detector(tmp_path: Path, log: object) -> detect_continuum:
    """Return a NIR-like detector (dispersion along y, flipped but not rotated) on a 64x40 trace frame."""
    detector = _plot_ready_detector(tmp_path, log, orderDeg=1, axisBDeg=1)
    detector.arm = "NIR"
    detector.axisA = "y"
    detector.axisB = "x"
    detector.detectorParams = {"dispersion-axis": "y", "flip-qc-plot": 1, "rotate-qc-plot": 0}
    detector.traceFrame = CCDData(np.full((FLIP_ROWS, FLIP_COLUMNS), 20.0), unit=u.electron)
    detector.qc = qc_table()
    return detector


def _nir_trace_pixels(orders: list[float], centres: list[float]) -> pd.DataFrame:
    """Return a residual-free pixel table with one constant-row trace per order, dispersed along x."""
    xValues = np.arange(4.0, 36.0)
    return pd.DataFrame(
        {
            "order": np.repeat(orders, len(xValues)),
            "cont_y": np.repeat(centres, len(xValues)),
            "cont_x": np.tile(xValues, len(orders)),
            "cont_y_fit_res": np.zeros(len(xValues) * len(orders)),
            "gauss_stddev_fit": np.full(len(xValues) * len(orders), 1.5),
            "wavelength": 500.0 + np.tile(xValues, len(orders)),
        }
    )


def _nir_clipped_data() -> pd.DataFrame:
    """Return one sample with no continuum found and one clipped peak."""
    return pd.DataFrame(
        {
            "cont_x": [np.nan, 20.0],
            "cont_y": [10.0, 30.0],
            "fit_x": [10.0, 20.0],
            "fit_y": [10.0, 25.0],
        }
    )


def _finite_rows(collection: object) -> list[float]:
    """Return the sorted distinct finite y-values of a scatter collection."""
    rows = np.asarray(collection.get_offsets())[:, 1]
    return sorted({float(row) for row in rows[np.isfinite(rows)]})


def _probe_figure(fig: plt.Figure) -> dict[str, object]:
    """Report the figure properties the continuum QC plot tests assert on."""
    extents = table_extents(fig)
    toprow, midrow = fig.axes[:2]
    return {
        "interpolations": image_interpolations(fig),
        "tableExtents": extents,
        "peakRows": [_finite_rows(collection) for collection in toprow.collections],
        "notFoundRows": [list(line.get_ydata()) for line in toprow.get_lines()],
        "fitRows": sorted({float(np.unique(line.get_ydata())[0]) for line in midrow.get_lines()}),
        "bandRows": sorted(
            (float(path.vertices[:, 1].min()), float(path.vertices[:, 1].max()))
            for collection in midrow.collections
            for path in collection.get_paths()
        ),
    }


@pytest.fixture
def savefig_probes(monkeypatch: pytest.MonkeyPatch):
    """Probe the current figure each time `plt.savefig` is called, before `plot_results` clears its axes."""
    probes = probe_on_savefig(monkeypatch, _probe_figure)
    yield probes
    plt.close("all")


def test_plot_results_embeds_images_without_resampling(tmp_path, log: object, savefig_probes: list) -> None:
    """Embed the trace frame at native resolution so vector markers line up with its pixels when zoomed."""
    # ARRANGE
    detector = _nir_flip_detector(tmp_path, log)
    pixels = _nir_trace_pixels([12.0, 13.0], [14.0, 22.0])

    # ACT
    detector.plot_results(pixels, pd.DataFrame([FLIP_COEFFICIENTS]), _nir_clipped_data())

    # ASSERT
    assert savefig_probes[0]["interpolations"] == ["none", "none"]


def test_plot_results_flip_only_maps_markers_and_fit_lines_to_flipped_pixel_rows(
    tmp_path, log: object, savefig_probes: list
) -> None:
    """A flip without rotation sends row ``r`` to ``N - 1 - r`` for peaks, missing-continuum bars, fits and bands."""
    # ARRANGE
    detector = _nir_flip_detector(tmp_path, log)
    pixels = _nir_trace_pixels([12.0, 13.0], [14.0, 22.0])

    # ACT
    detector.plot_results(pixels, pd.DataFrame([FLIP_COEFFICIENTS]), _nir_clipped_data())

    # ASSERT
    probe = savefig_probes[0]
    last = FLIP_ROWS - 1
    assert probe["peakRows"][0] == [last - 22.0, last - 14.0]
    assert probe["peakRows"][1] == [last - 30.0, last - 10.0]
    assert probe["notFoundRows"] == [[last - 10.0 - 2, last - 10.0 + 2]]
    assert probe["fitRows"] == [last - 22.0, last - 14.0]
    assert probe["bandRows"] == [
        (last - 22.0 - 4.5, last - 22.0 + 4.5),
        (last - 14.0 - 4.5, last - 14.0 + 4.5),
    ]


def test_plot_results_leaves_the_callers_pixel_and_clipped_tables_unflipped(tmp_path, log: object) -> None:
    """Flip the plotted rows without rewriting the detector coordinates the caller passed in."""
    # ARRANGE
    detector = _nir_flip_detector(tmp_path, log)
    pixels = _nir_trace_pixels([12.0, 13.0], [14.0, 22.0])
    clipped = _nir_clipped_data()
    pixelsBefore, clippedBefore = pixels.copy(), clipped.copy()

    # ACT
    detector.plot_results(pixels, pd.DataFrame([FLIP_COEFFICIENTS]), clipped)
    plt.close("all")

    # ASSERT
    pd.testing.assert_frame_equal(pixels, pixelsBefore)
    pd.testing.assert_frame_equal(clipped, clippedBefore)


def test_plot_results_tables_do_not_overlap(tmp_path, log: object, savefig_probes: list) -> None:
    """Keep the QC table and the settings table apart on a tall VIS-shaped trace frame."""
    # ARRANGE
    detector = _plot_ready_detector(tmp_path, log)
    detector.detectorParams = {"dispersion-axis": "x", "flip-qc-plot": 1, "rotate-qc-plot": 90}
    detector.traceFrame = CCDData(np.full((410, 82), 20.0), unit=u.electron)
    detector.recipeSettings.update({f"setting-number-{i}": i for i in range(14)})
    detector.qc = pd.concat(
        [qc_table().assign(qc_name=f"QC {i}", qc_comment=f"[px] synthetic QC number {i}") for i in range(8)],
        ignore_index=True,
    )
    pixels = _trace_pixels([10.0, 11.0], slope=0.2, intercept=8.0)
    coefficients = pd.DataFrame([{"cent_0": 8.0, "cent_1": 0.2, "std_0": 1.5, "std_1": 0.0}])

    # ACT
    detector.plot_results(pixels, coefficients, _no_continuum_clipped())

    # ASSERT
    qcExtent, settingsExtent = savefig_probes[0]["tableExtents"]
    assert not qcExtent.overlaps(settingsExtent)


def test_plot_results_order_metadata_for_unflipped_frame_is_pinned(tmp_path, log: object) -> None:
    """Pin the order limits written into the order-table product: axis A holds the fitted detector pixel (DY-1219)."""
    # ARRANGE
    detector = _plot_ready_detector(tmp_path, log, orderDeg=1, axisBDeg=1)
    detector.qc = qc_table()
    pixels = _trace_pixels([11.0, 12.0], slope=0.0, intercept=16.0)

    # ACT
    orderMetadata = detector.build_order_metadata(pixels, pd.DataFrame([FLIP_COEFFICIENTS]))

    # ASSERT
    assert orderMetadata.to_dict("list") == {
        "order": [11.0, 12.0],
        "ymin": [0, 0],
        "ymax": [30, 30],
        "xmin": [6.0, 14.0],
        "xmax": [6.0, 14.0],
    }


def test_plot_results_order_metadata_for_flip_only_frame_is_pinned(tmp_path, log: object) -> None:
    """Pin the order limits written into the order-table product for a flipped, unrotated frame."""
    # ARRANGE
    detector = _nir_flip_detector(tmp_path, log)
    pixels = _nir_trace_pixels([12.0, 13.0], [14.0, 22.0])

    # ACT
    orderMetadata = detector.build_order_metadata(pixels, pd.DataFrame([FLIP_COEFFICIENTS]))

    # ASSERT
    assert orderMetadata.to_dict("list") == {
        "order": [12.0, 13.0],
        "xmin": [0, 0],
        "xmax": [39, 39],
        "ymin": [14.0, 22.0],
        "ymax": [14.0, 22.0],
    }


# THREE ORIENTATIONS WITH A NON-SQUARE TRACE FRAME, SO A MIRRORED AXIS-A VALUE NEVER EQUALS THE DETECTOR VALUE
SLOPED_COEFFICIENTS = {**FLIP_COEFFICIENTS, "cent_01": 0.25}
ORIENTATIONS = {
    "flip-only-nir": {"arm": "NIR", "axisA": "y", "axisB": "x", "flip": 1, "rotate": 0, "shape": (64, 40)},
    "rotate-and-flip-vis": {"arm": "VIS", "axisA": "x", "axisB": "y", "flip": 1, "rotate": 90, "shape": (96, 80)},
    "no-rotate-no-flip": {"arm": "VIS", "axisA": "x", "axisB": "y", "flip": 0, "rotate": 0, "shape": (96, 80)},
}
# ORDERS 11 AND 12 SIT AT -82 + 8 * ORDER = 6 AND 14 PIXELS AT AXIS B = 0, RISING 0.25 PIXELS PER AXIS-B PIXEL
ORIENTATION_ORDERS = [11.0, 12.0]


def _oriented_detector(
    tmp_path: Path, log: object, orientation: dict, *, flip: int | None = None, rotate: int | None = None
) -> detect_continuum:
    """Return a detector for one named orientation, optionally overriding its QC-plot flip and rotation."""
    detector = _plot_ready_detector(tmp_path, log, orderDeg=1, axisBDeg=1)
    detector.arm = orientation["arm"]
    detector.axisA = orientation["axisA"]
    detector.axisB = orientation["axisB"]
    detector.detectorParams = {
        "dispersion-axis": orientation["axisB"],
        "flip-qc-plot": orientation["flip"] if flip is None else flip,
        "rotate-qc-plot": orientation["rotate"] if rotate is None else rotate,
    }
    detector.traceFrame = CCDData(np.full(orientation["shape"], 20.0), unit=u.electron)
    detector.qc = qc_table()
    return detector


def _oriented_pixels(orientation: dict) -> pd.DataFrame:
    """Return a residual-free pixel table for the orientation's two orders, in detector coordinates."""
    axisA, axisB = orientation["axisA"], orientation["axisB"]
    pixels = _trace_pixels(ORIENTATION_ORDERS, slope=0.0, intercept=0.0)
    return pixels.rename(columns={"cont_x": f"cont_{axisA}", "cont_y": f"cont_{axisB}"}).rename(
        columns={"cont_x_fit_res": f"cont_{axisA}_fit_res"}
    )


def _expected_axis_a_limits(orientation: dict) -> tuple[list[float], list[float]]:
    """Return the min and max detector axis-A positions of each order's fit over the sampled axis-B grid."""
    axisBLength = orientation["shape"][1] if orientation["axisB"] == "x" else orientation["shape"][0]
    lastSample = 3 * ((axisBLength - 1) // 3)
    lows = [-82.0 + 8.0 * order for order in ORIENTATION_ORDERS]
    return lows, [low + 0.25 * lastSample for low in lows]


def _order_metadata(detector: detect_continuum, orientation: dict) -> pd.DataFrame:
    return detector.build_order_metadata(_oriented_pixels(orientation), pd.DataFrame([SLOPED_COEFFICIENTS]))


@pytest.mark.parametrize("name", list(ORIENTATIONS))
def test_order_metadata_axis_a_limits_are_the_fitted_detector_positions(tmp_path, log: object, name: str) -> None:
    """Publish ``{axisA}min/max`` as the min and max fitted axis-A pixel, whatever the QC plot orientation."""
    # ARRANGE
    orientation = ORIENTATIONS[name]
    detector = _oriented_detector(tmp_path, log, orientation)
    expectedMin, expectedMax = _expected_axis_a_limits(orientation)

    # ACT
    orderMetadata = _order_metadata(detector, orientation)

    # ASSERT
    assert orderMetadata[f"{orientation['axisA']}min"].tolist() == pytest.approx(expectedMin)
    assert orderMetadata[f"{orientation['axisA']}max"].tolist() == pytest.approx(expectedMax)


# THE AXIS-B LIMITS ARE THE FIRST AND LAST SAMPLE OF np.arange(0, AXIS-B LENGTH, 3)
EXPECTED_AXIS_B_LIMITS = {"flip-only-nir": 39, "rotate-and-flip-vis": 93, "no-rotate-no-flip": 93}


@pytest.mark.parametrize("name", list(ORIENTATIONS))
def test_order_metadata_keeps_the_order_list_and_axis_b_limits(tmp_path, log: object, name: str) -> None:
    """Pin the order list and ``{axisB}min/max``, which this fix must leave alone."""
    # ARRANGE
    orientation = ORIENTATIONS[name]
    detector = _oriented_detector(tmp_path, log, orientation)
    axisB = orientation["axisB"]

    # ACT
    orderMetadata = _order_metadata(detector, orientation)

    # ASSERT
    assert orderMetadata["order"].tolist() == ORIENTATION_ORDERS
    assert orderMetadata[f"{axisB}min"].tolist() == [0, 0]
    assert orderMetadata[f"{axisB}max"].tolist() == [EXPECTED_AXIS_B_LIMITS[name]] * 2
    assert list(orderMetadata.columns) == [
        "order",
        f"{axisB}min",
        f"{axisB}max",
        f"{orientation['axisA']}min",
        f"{orientation['axisA']}max",
    ]


@pytest.mark.parametrize("name", list(ORIENTATIONS))
def test_order_metadata_is_computed_without_plotting(
    tmp_path, log: object, monkeypatch: pytest.MonkeyPatch, name: str
) -> None:
    """Build the order metadata with the plotting code stubbed out to fail on any call."""
    # ARRANGE
    orientation = ORIENTATIONS[name]
    detector = _oriented_detector(tmp_path, log, orientation)
    expectedMin, expectedMax = _expected_axis_a_limits(orientation)

    def plotting_forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError("the order metadata must not call plot_results")

    monkeypatch.setattr(detector, "plot_results", plotting_forbidden)
    monkeypatch.setattr(plt, "figure", plotting_forbidden)

    # ACT
    orderMetadata = _order_metadata(detector, orientation)

    # ASSERT
    assert orderMetadata[f"{orientation['axisA']}min"].tolist() == pytest.approx(expectedMin)
    assert orderMetadata[f"{orientation['axisA']}max"].tolist() == pytest.approx(expectedMax)


@pytest.mark.parametrize("flip, rotate", [(0, 0), (1, 0), (1, 90), (0, 90), (1, 270)])
def test_order_metadata_does_not_change_with_the_qc_plot_orientation(
    tmp_path, log: object, flip: int, rotate: int
) -> None:
    """Change ``flip-qc-plot`` and ``rotate-qc-plot`` on one frame and get identical order metadata."""
    # ARRANGE
    orientation = ORIENTATIONS["no-rotate-no-flip"]
    reference = _order_metadata(_oriented_detector(tmp_path, log, orientation), orientation)
    detector = _oriented_detector(tmp_path, log, orientation, flip=flip, rotate=rotate)

    # ACT
    orderMetadata = _order_metadata(detector, orientation)

    # ASSERT
    pd.testing.assert_frame_equal(orderMetadata, reference)


def test_plot_results_returns_only_the_plot_path(tmp_path, log: object, savefig_text_colours: list) -> None:
    """Keep ``plot_results`` to drawing: it returns the PDF path and no metadata."""
    # ARRANGE
    orientation = ORIENTATIONS["no-rotate-no-flip"]
    detector = _oriented_detector(tmp_path, log, orientation)

    # ACT
    outputPath = detector.plot_results(
        _oriented_pixels(orientation), pd.DataFrame([SLOPED_COEFFICIENTS]), _no_continuum_clipped()
    )

    # ASSERT
    assert isinstance(outputPath, str)
    assert Path(outputPath).read_bytes().startswith(b"%PDF")


def _develop_fit_curves(detector: detect_continuum, orderPixelTable: pd.DataFrame, orderPolyTable: pd.DataFrame):
    """Evaluate the order fits with the zip-and-filter code that `plot_results` carried before DY-1219."""
    axisALength = detector.traceFrame.data.shape[1] if detector.axisB == "y" else detector.traceFrame.data.shape[0]
    axisBLength = detector.traceFrame.data.shape[0] if detector.axisB == "y" else detector.traceFrame.data.shape[1]
    axisBLineList = np.arange(0, axisBLength, 3)
    poly = continuumModule.chebyshev_order_xy_polynomials(
        log=detector.log,
        axisBCol=detector.axisB,
        orderCol="order",
        orderDeg=detector.orderDeg,
        axisBDeg=detector.axisBDeg,
    ).poly
    row = orderPolyTable.iloc[-1]
    centCoeff = [float(v) for k, v in row.items() if "cent_" in k]
    stdCoeff = [float(v) for k, v in row.items() if "std_" in k]
    df = pd.DataFrame({detector.axisB: axisBLineList})
    curves = []
    for o in orderPixelTable["order"].unique():
        df["order"] = o
        xfit = poly(df, *centCoeff)
        stdfit = poly(df, *stdCoeff)
        try:
            xfit, yfit, stdfit = zip(
                *[
                    (x, y, std)
                    for x, y, std in zip(xfit, axisBLineList, stdfit, strict=True)
                    if x > 0 and x < axisALength - 10
                ],
                strict=True,
            )
        except ValueError:
            continue
        curves.append((o, np.array(xfit), np.array(yfit), np.array(stdfit)))
    return curves


def test_fitted_order_curves_match_the_pre_refactor_evaluation_exactly(tmp_path, log: object) -> None:
    """Keep the plotted curves bit-identical, including partly and wholly filtered-out orders."""
    # ARRANGE
    orientation = ORIENTATIONS["no-rotate-no-flip"]
    detector = _oriented_detector(tmp_path, log, orientation)
    # ORDER 5 IS ALWAYS OFF THE LOW EDGE, 10 CROSSES IT, 11 AND 12 ARE ON THE DETECTOR, 20 IS ALWAYS OFF THE HIGH EDGE
    orders = [5.0, 10.0, 11.0, 12.0, 20.0]
    pixels = _trace_pixels(orders, slope=0.0, intercept=0.0)
    coefficients = pd.DataFrame([SLOPED_COEFFICIENTS])
    expected = _develop_fit_curves(detector, pixels, coefficients)

    # ACT
    curves = detector._fitted_order_curves(pixels, coefficients)

    # ASSERT
    assert [curve[0] for curve in curves] == [10.0, 11.0, 12.0]
    assert len(curves) == len(expected)
    for curve, oldCurve in zip(curves, expected, strict=True):
        assert curve[0] == oldCurve[0]
        for new, old in zip(curve[1:], oldCurve[1:], strict=True):
            assert np.array_equal(new, old)
    assert len(curves[0][1]) < len(curves[1][1])


def test_order_colours_assigns_a_stable_distinct_colour_per_order_and_wraps_after_the_cycle() -> None:
    # ARRANGE
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    orders = [float(o) for o in range(10, 10 + len(cycle) + 2)]

    # ACT
    colours = continuumModule._order_colours(orders)

    # ASSERT
    assert list(colours) == orders
    assert [colours[o] for o in orders[: len(cycle)]] == cycle
    assert len({colours[o] for o in orders[: len(cycle)]}) == len(cycle)
    assert colours[orders[len(cycle)]] == cycle[0]
    assert colours[orders[len(cycle) + 1]] == cycle[1]
    assert colours == continuumModule._order_colours(orders)
