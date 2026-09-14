"""Analytic continuum-trace fitting contracts."""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe.commonutils.detect_continuum import detect_continuum
from tests.factories import instrument_header, pipeline_settings

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
    expectedFit = (
        1.0
        + 2.0 * pixels["cont_y"]
        + 3.0 * pixels["order"]
        + 4.0 * pixels["order"] * pixels["cont_y"]
    )
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
    assert set(fitted["source"]) == set(range(10))
    assert set(clipped["source"]) == {10, 11}
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
        lambda **kwargs: ("/qc/continuum.pdf", pd.DataFrame({"axis": ["x"]})),
    )
    monkeypatch.setattr(detector, "write_order_table_to_file", lambda **kwargs: "/products/orders.fits")

    result = detector.get()

    orderPath, qc, products, orderPoly, orderPixels, orderMeta = result
    assert orderPath == "/products/orders.fits"
    assert fittedCalls == ["cont_x", "gauss_stddev"]
    assert qc["qc_name"].tolist() == ["SAMPLES CLIP NUM", "SAMPLES CLIP FRAC"]
    assert qc["qc_value"].tolist() == [1, "0.500"]
    assert products["product_label"].tolist() == ["ORDER_CENTRES_RES"]
    assert products["file_name"].tolist() == ["continuum.pdf"]
    assert orderPoly.loc[0, "cent_11"] == 4.0
    assert orderPoly.loc[0, "std_11"] == 4.0
    assert "mask" not in orderPixels
    assert orderMeta.to_dict("list") == {"axis": ["x"]}


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
    fits.PrimaryHDU(header=fits.Header({"BINX": 2, "BINY": 3})).writeto(
        detector.dispersion_map
    )
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
    detector.recipeSettings.update(
        {"slice-length": 12, "peak-sigma-limit": 3.0, "slice-width": 5}
    )
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
    assert detector.qc["qc_value"].tolist() == [50, 50, "1.000"]


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
    detector.recipeSettings.update(
        {"slice-length": 12, "peak-sigma-limit": 3.0, "slice-width": 5}
    )
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


def test_plot_results_writes_continuum_diagnostic_and_order_limits(
    tmp_path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Render the fitted trace diagnostic from a compact analytic order table."""
    detector = _detector(log)
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
    yValues = np.arange(4.0, 28.0)
    orders = np.repeat([10.0, 11.0], len(yValues))
    traceY = np.tile(yValues, 2)
    traceX = 8.0 + 0.2 * traceY
    pixels = pd.DataFrame(
        {
            "order": orders,
            "cont_x": traceX,
            "cont_y": traceY,
            "cont_x_fit_res": np.zeros(len(orders)),
            "gauss_stddev_fit": np.full(len(orders), 1.5),
            "wavelength": 500.0 + traceY,
        }
    )
    clipped = pd.DataFrame(
        {
            "cont_x": [np.nan],
            "cont_y": [10.0],
            "fit_x": [10.0],
            "fit_y": [10.0],
        }
    )
    coefficients = pd.DataFrame(
        [{"cent_0": 8.0, "cent_1": 0.2, "std_0": 1.5, "std_1": 0.0}]
    )
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.qc_settings_plot_tables",
        lambda **_: None,
    )

    outputPath, orderMetadata = detector.plot_results(pixels, coefficients, clipped)

    assert Path(outputPath).read_bytes().startswith(b"%PDF")
    assert orderMetadata["order"].tolist() == [10.0, 11.0]
    assert (orderMetadata["xmax"] > orderMetadata["xmin"]).all()
