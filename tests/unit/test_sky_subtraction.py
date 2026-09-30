"""Analytic residual contracts for sky subtraction."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from scipy.interpolate import splrep

from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header

pytestmark = pytest.mark.unit


def _subtractor(log: object) -> subtract_sky:
    subtractor = object.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "residual_floor_percentile": 50,
            "noise_sigma": 0,
        }
    }
    return subtractor


def test_constructor_prepares_vis_sky_subtraction_metadata_and_workspace(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Initialize the public sky-subtraction workflow from a prepared VIS frame."""
    import importlib

    import soxspipe.commonutils.toolkit as toolkit

    module = importlib.import_module("soxspipe.commonutils.subtract_sky")
    header = instrument_header()
    header["ESO INS VISE NAME"] = "SLIT1.0"
    header["ESO DET BINX"] = 2
    header["ESO DET BINY"] = 1
    frame = CCDData(np.ones((3, 3)), unit=u.electron, meta=header)
    mapTable = pd.DataFrame(
        {"order": [10, 10], "wavelength": [500.0, 501.0], "slit_position": [0.0, 0.1]}
    )

    class DetectorLookup:
        def __init__(self, **_: object) -> None:
            return None

        def get(self, arm: str) -> dict[str, object]:
            assert arm == "VIS"
            return {"dispersion-axis": "x", "slit_length": 11.0}

    monkeypatch.setattr(module, "detector_lookup", DetectorLookup)
    monkeypatch.setattr(
        module,
        "twoD_disp_map_image_to_dataframe",
        lambda **_: (mapTable.copy(), np.zeros((3, 3), dtype=bool)),
    )
    monkeypatch.setattr(module, "quicklook_image", lambda **_: None)
    monkeypatch.setattr(
        toolkit,
        "utility_setup",
        lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")),
    )

    subtractor = subtract_sky(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={"sky-subtraction": {"bspline_order": 3}},
        objectFrame=frame,
        twoDMap="two-d-map.fits",
        qcTable=pd.DataFrame(),
        productsTable=pd.DataFrame(),
        dispMap="dispersion-map.fits",
        sofName="SYNTHETIC",
        startNightDate="2024-01-02",
    )

    assert subtractor.filenameTemplate == "SYNTHETIC.fits"
    assert subtractor.arm == "VIS"
    assert subtractor.axisA == "x"
    assert subtractor.axisB == "y"
    assert subtractor.binx == 2
    assert subtractor.biny == 1
    assert subtractor.mapDF.equals(mapTable)


@pytest.mark.parametrize("dispersionAxis", ["x", "y"])
def test_placeholder_images_receive_model_flux_and_normalised_residuals(
    log: object, dispersionAxis: str
) -> None:
    subtractor = _subtractor(log)
    subtractor.objectFrame = CCDData(
        np.arange(9, dtype=float).reshape(3, 3),
        unit="electron",
        mask=np.zeros((3, 3), dtype=bool),
    )
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.detectorParams = {"dispersion-axis": dispersionAxis}
    pixels = pd.DataFrame(
        {
            "x": [1],
            "y": [2],
            "sky_model": [4.0],
            "sky_subtracted_flux": [6.0],
            "error": [2.0],
        }
    )

    model, subtracted, residuals = subtractor.create_placeholder_images()
    model, subtracted, residuals = subtractor.add_data_to_placeholder_images(
        pixels, model, subtracted, residuals
    )

    indices = (2, 1) if dispersionAxis == "x" else (1, 2)
    assert model.data[indices] == 4.0
    assert subtracted.data[indices] == 6.0
    assert residuals.data[indices] == 3.0
    assert np.count_nonzero(np.isfinite(model.data)) == 1


def test_subtract_runs_one_order_workflow_and_returns_finite_product_frames(
    log: object,
) -> None:
    """The public subtraction workflow fills aligned model and residual outputs."""
    subtractor = _subtractor(log)
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.debug = False
    subtractor.stopSubtraction = False
    subtractor.detectorParams = {"dispersion-axis": "x"}
    subtractor.dateObs = "2024-01-02T03:04:05"
    subtractor.objectFrame = CCDData(
        np.full((2, 2), 10.0),
        unit=u.electron,
        mask=np.zeros((2, 2), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((2, 2))),
    )
    subtractor.mapDF = pd.DataFrame(
        {
            "order": [10, 10],
            "x": [0, 1],
            "y": [0, 1],
            "sky_model": [4.0, 5.0],
            "sky_subtracted_flux": [6.0, 5.0],
            "error": [2.0, 1.0],
        }
    )
    subtractor.qc = pd.DataFrame()
    subtractor.products = pd.DataFrame()
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "bspline_order": 3,
            "clip-slit-edge-fraction": 0.1,
            "aggressive_object_masking": False,
            "sky_model_qc_plot": False,
        }
    }
    subtractor.get_over_sampled_sky_from_order = lambda frame, **_: frame.copy()
    subtractor.clip_object_slit_positions = lambda orders, **_: orders

    def fit(frame: pd.DataFrame) -> tuple[pd.DataFrame, tuple[object, object, int], np.ndarray, np.ndarray, float]:
        return frame.copy(), (None, None, 3), np.array([501.0]), np.array([1.0, 2.0]), 0.5

    subtractor.fit_bspline_curve_to_sky = fit

    model, subtracted, residuals, qc, products = subtractor.subtract()

    np.testing.assert_allclose(model.data, [[4.0, 0.0], [0.0, 5.0]])
    np.testing.assert_allclose(subtracted.data, [[6.0, 0.0], [0.0, 5.0]])
    np.testing.assert_allclose(residuals.data, [[3.0, np.nan], [np.nan, 5.0]], equal_nan=True)
    assert qc.empty
    assert products.empty


def test_subtract_records_order_and_frame_qc_products_when_requested(
    log: object,
    tmp_path: Path,
) -> None:
    """Sky-model diagnostics are registered with their stable product labels."""
    subtractor = _subtractor(log)
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.debug = True
    subtractor.stopSubtraction = False
    subtractor.detectorParams = {"dispersion-axis": "x"}
    subtractor.dateObs = "2024-01-02T03:04:05"
    subtractor.objectFrame = CCDData(
        np.full((2, 2), 10.0),
        unit=u.electron,
        mask=np.zeros((2, 2), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((2, 2))),
    )
    subtractor.mapDF = pd.DataFrame(
        {
            "order": [10, 10],
            "x": [0, 1],
            "y": [0, 1],
            "sky_model": [4.0, 5.0],
            "sky_subtracted_flux": [6.0, 5.0],
            "error": [2.0, 1.0],
        }
    )
    subtractor.qc = pd.DataFrame()
    subtractor.products = pd.DataFrame()
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "bspline_order": 3,
            "clip-slit-edge-fraction": 0.1,
            "aggressive_object_masking": False,
            "sky_model_qc_plot": True,
        }
    }
    orderPdf = tmp_path / "sky-order.pdf"
    comparisonPdf = tmp_path / "sky-comparison.pdf"
    orderPdf.touch()
    comparisonPdf.touch()
    subtractor.get_over_sampled_sky_from_order = lambda frame, **_: frame.copy()
    subtractor.clip_object_slit_positions = lambda orders, **_: orders
    subtractor.fit_bspline_curve_to_sky = lambda frame: (
        frame.copy(),
        (None, None, 3),
        np.array([501.0]),
        np.array([1.0, 2.0]),
        0.5,
    )
    subtractor.plot_sky_sampling = lambda **_: str(orderPdf)
    subtractor.plot_image_comparison = lambda *_: str(comparisonPdf)

    _, _, _, _, products = subtractor.subtract()

    assert products["product_label"].tolist() == [
        "SKY_MODEL_QC_PLOTS",
        "SKY SUBTRACTION QUICKLOOK",
    ]
    assert products["file_path"].tolist() == [str(orderPdf), str(comparisonPdf)]
    assert products["label"].tolist() == ["QC", "QC"]


def test_oversampled_sky_clips_requested_slit_edges_and_bad_pixels(
    log: object,
) -> None:
    """Sampling marks detector defects and requested slit-edge exclusions."""
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_clipping_sigma": 4,
            "percentile_clipping_iterations": 2,
            "percentile_rolling_window_size": 7,
        }
    )
    pixels = pd.DataFrame(
        {
            "slit_position": [-1.0, 0.0, 1.0],
            "mask": [False, True, False],
            "flux": [10.0, 11.0, 12.0],
        }
    )
    calls: list[dict[str, object]] = []

    def record_clipping(**kwargs: object) -> pd.DataFrame:
        calls.append(kwargs)
        return cast(pd.DataFrame, kwargs["imageMapOrderDF"]).copy()

    subtractor.rolling_window_clipping = record_clipping

    clipped = subtractor.get_over_sampled_sky_from_order(
        pixels.copy(),
        clipBPs=True,
        clipSlitEdge=0.25,
    )

    assert clipped["flagged_edge_clipped"].tolist() == [True, False, True]
    assert clipped["flagged_bad_pixel_clipped"].tolist() == [False, True, False]
    assert clipped["flagged_all_clipped"].tolist() == [True, True, True]
    assert calls[0]["windowSize"] == 7
    assert calls[0]["sigma_clip_limit"] == 4
    assert calls[0]["max_iterations"] == 2


def test_rolling_window_clipping_rejects_an_isolated_object_peak(log: object) -> None:
    """An isolated high-flux sample is excluded from subsequent sky fitting."""
    subtractor = _subtractor(log)
    subtractor.arm = "UVB"
    subtractor.debug = False
    wavelengths = np.arange(101, dtype=float)
    pixels = pd.DataFrame(
        {
            "order": np.full(wavelengths.size, 10),
            "wavelength": wavelengths,
            "flux": np.where(wavelengths == 50, 100.0, 10.0),
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
        }
    )

    clipped = subtractor.rolling_window_clipping(
        pixels,
        windowSize=11,
        sigma_clip_limit=3,
        max_iterations=3,
    )

    assert clipped.loc[50, "flagged_object_clipped"]
    assert clipped.loc[50, "flagged_all_clipped"]
    assert not clipped.loc[49, "flagged_object_clipped"]


def test_clip_object_slit_positions_updates_masks_and_local_noise(log: object) -> None:
    """Object candidates are excluded before local sky-noise measurements."""
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_rolling_window_size": 5,
            "noise_rolling_window_size": 30,
        }
    )
    slitPositions = np.linspace(-1.0, 1.0, 40)
    pixels = pd.DataFrame(
        {
            "slit_position": slitPositions,
            "flagged_object_clipped": slitPositions > 0.7,
            "flagged_all_clipped": False,
            "flux_minus_smoothed_residual": np.sin(slitPositions),
            "flux_percentile_smoothed": 100.0 + slitPositions,
        }
    )

    result = subtractor.clip_object_slit_positions([pixels])

    assert result[0] is pixels
    assert pixels.loc[pixels["slit_position"] > 0.7, "flagged_all_clipped"].all()
    assert pixels.loc[pixels["slit_position"] <= 0.7, "flagged_all_clipped"].eq(False).all()
    assert pixels.loc[10, "residual_windowed_std"] > 0
    assert pixels.loc[20, "flux_windowed_long_median"] == pytest.approx(99.948718, abs=1e-6)


def test_cross_dispersion_normaliser_keeps_the_current_unity_flux_contract(log: object) -> None:
    """The fitted illumination profile currently preserves a unity correction."""
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"]["slit_illumination_order"] = 1
    subtractor.binx = 1
    subtractor.biny = 1
    subtractor.debug = False
    slitPositions = np.linspace(-1.0, 1.0, 20)
    pixels = pd.DataFrame(
        {
            "order": np.full(slitPositions.size, 10),
            "slit_position": slitPositions,
            "flagged_all_clipped": False,
            "sky_subtracted_flux": 20.0 + 3.0 * slitPositions,
            "residual_windowed_std": np.full(slitPositions.size, 0.5),
        }
    )

    corrected = subtractor.cross_dispersion_flux_normaliser(pixels)

    assert corrected is pixels
    np.testing.assert_allclose(corrected["slit_normalisation_ratio"], 1.0)


@pytest.mark.parametrize("arm", ["UVB", "VIS"])
def test_image_comparison_writes_recipe_named_pdf(
    tmp_path: Path, log: object, arm: str
) -> None:
    subtractor = _subtractor(log)
    subtractor.arm = arm
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.mapDF = pd.DataFrame({"x": [0, 1], "y": [0, 1]})
    subtractor.filenameTemplate = "science.fits"
    subtractor.qcDir = str(tmp_path)
    frame = CCDData(
        np.arange(16, dtype=float).reshape(4, 4),
        unit="electron",
        mask=np.zeros((4, 4), dtype=bool),
    )

    outputPath = subtractor.plot_image_comparison(frame, frame.copy(), frame.copy())

    assert Path(outputPath).read_bytes().startswith(b"%PDF")


def test_sky_sampling_plot_writes_a_complete_order_diagnostic(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Render the complete order diagnostic with a deterministic fitted sky."""
    subtractor = _subtractor(log)
    subtractor.arm = "VIS"
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.filenameTemplate = "science.fits"
    subtractor.qcDir = str(tmp_path)
    subtractor.detectorParams = {
        "dispersion-axis": "x",
        "flip-qc-plot": False,
        "rotate-qc-plot": False,
    }
    subtractor.objectFrame = CCDData(
        np.full((16, 16), 12.0),
        unit=u.electron,
        mask=np.zeros((16, 16), dtype=bool),
    )
    coordinates = np.arange(256)
    wavelength = np.linspace(500.0, 510.0, 256)
    flags = coordinates % 53 == 0
    pixels = pd.DataFrame(
        {
            "x": coordinates % 16,
            "y": coordinates // 16,
            "wavelength": wavelength,
            "flux": np.full(256, 12.0),
            "slit_position": np.linspace(-1.0, 1.0, 256),
            "residual_windowed_std": np.ones(256),
            "residual_global_sigma": np.where(flags, 2.0, 0.0),
            "flux_percentile_smoothed": np.full(256, 12.0),
            "sky_model": np.full(256, 12.0),
            "sky_subtracted_flux": np.zeros(256),
            "sky_subtracted_flux_weighted": np.zeros(256),
            "flagged_all_clipped": flags,
            "flagged_object_clipped": flags,
            "flagged_bad_pixel_clipped": flags,
            "flagged_edge_clipped": flags,
            "flagged_bspline_clipped": flags,
        }
    )
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None
    )

    outputPath = subtractor.plot_sky_sampling(
        order=10,
        imageMapOrderDF=pixels,
        knotLocations=np.array([502.0, 508.0]),
    )

    assert Path(outputPath).read_bytes().startswith(b"%PDF")


def test_calculate_residuals_evaluates_polynomial_and_propagates_nan(
    log: object,
) -> None:
    subtractor = _subtractor(log)
    pixels = pd.DataFrame(
        {
            "order_pow_0": [1.0, 1.0, 1.0],
            "wavelength_pow_0": [1.0, 1.0, 1.0],
            "wavelength_pow_1": [2.0, 3.0, 4.0],
            "slit_position_pow_0": [1.0, 1.0, 1.0],
            "sky_subtracted_flux": [4.0, 8.0, np.nan],
        }
    )

    mean, std, median, result = subtractor.calculate_residuals(
        pixels.copy(),
        fluxcoeff=[1.0, 2.0],
        orderDeg=0,
        wavelengthDeg=1,
        slitDeg=0,
    )

    np.testing.assert_allclose(
        result["fit_sky_subtracted_flux"],
        [5.0, 7.0, 9.0],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        result["residuals_sky_subtracted_flux"],
        [1.0, -1.0, np.nan],
        rtol=1e-12,
        atol=1e-12,
        equal_nan=True,
    )
    assert mean == pytest.approx(0.0)
    assert std == pytest.approx(1.0)
    assert np.isnan(median)


def test_object_clipping_marks_flagged_pixels_and_refreshes_local_noise(
    log: object,
) -> None:
    """The sky sampler excludes object pixels while retaining sampled sky rows."""
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_rolling_window_size": 5,
            "noise_rolling_window_size": 30,
        }
    )
    slitPositions = np.linspace(-1.0, 1.0, 64)
    pixels = pd.DataFrame(
        {
            "slit_position": slitPositions,
            "flagged_object_clipped": slitPositions > 0.5,
            "flagged_all_clipped": False,
            "flux_minus_smoothed_residual": np.sin(slitPositions),
            "flux_percentile_smoothed": 10.0 + np.cos(slitPositions),
            "residual_windowed_std": np.nan,
            "residual_windowed_long_median": np.nan,
            "flux_windowed_long_median": np.nan,
        }
    )

    [clipped] = subtractor.clip_object_slit_positions([pixels.copy()])

    assert clipped.loc[clipped["flagged_object_clipped"], "flagged_all_clipped"].all()
    skyRows = ~clipped["flagged_all_clipped"]
    assert clipped.loc[skyRows, "residual_windowed_std"].notna().sum() > 40
    assert clipped.loc[skyRows, "residual_windowed_long_median"].notna().sum() >= 20
    assert clipped.loc[skyRows, "flux_windowed_long_median"].notna().sum() >= 20


def test_residual_floor_normalises_errors_and_flags_local_excess(
    log: object,
) -> None:
    subtractor = _subtractor(log)
    wavelength = np.linspace(500.0, 501.0, 32)
    model = np.full(32, 8.0)
    pixels = pd.DataFrame(
        {
            "wavelength": wavelength,
            "flux": np.full(32, 10.0),
            "error": np.ones(32),
            "flagged_all_clipped": False,
            "residual_windowed_long_median": np.zeros(32),
            "flux_windowed_long_median": np.full(32, 10.0),
            "sky_residual_rolling_average": np.full(32, 2.0),
        }
    )
    pixels.loc[5, "error"] = 0.0
    pixels.loc[10, "sky_residual_rolling_average"] = 3.0
    spline = splrep(wavelength, model, k=1)

    result, residualFloor = subtractor.determine_residual_floor(
        pixels.copy(),
        spline,
        iteration=1,
    )

    assert residualFloor == 5
    assert result.loc[5, "sky_residuals"] == 1.0
    assert np.isfinite(result["sky_residuals"]).all()
    assert result.loc[10, "flagged_sky_line"] == "line"
    assert result.loc[0, "flagged_sky_line"] == False
    assert not result["flagged_noisy_region"].any()


def test_rectify_order_builds_a_grid_and_masks_clipped_detector_pixels(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Rectification preserves sampled flux while removing explicitly clipped pixels."""
    import importlib

    import soxspipe.commonutils.toolkit as toolkit

    module = importlib.import_module("soxspipe.commonutils.subtract_sky")
    subtractor = _subtractor(log)
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.arm = "VIS"
    subtractor.dispMap = "synthetic-dispersion-map.fits"
    subtractor.settings = {"instrument": "soxs"}
    subtractor.kw = lambda key: key
    subtractor.detectorParams = {"slit_length": 4.0}
    pixels = pd.DataFrame(
        {
            "x": [0, 1, 0],
            "y": [0, 0, 1],
            "flux": [15.0, 25.0, 0.0],
            "flagged_all_clipped": [False, True, False],
        }
    )

    monkeypatch.setattr(
        module,
        "read_spectral_format",
        lambda **_: (np.array([10]), np.array([500.0]), np.array([510.0])),
    )

    def map_to_detector(**kwargs: object) -> pd.DataFrame:
        table = cast(pd.DataFrame, kwargs["orderPixelTable"]).copy()
        table["fit_x"] = np.arange(len(table), dtype=float)
        table["fit_y"] = 0.0
        return table

    monkeypatch.setattr(module, "dispersion_map_to_pixel_arrays", map_to_detector)
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **_: None)

    image = subtractor._rectify_order(
        order=10,
        imageMapOrder=pixels,
        remove_clipped=True,
        conserve_flux=True,
    )

    assert image.shape == (2, 1)
    assert image[0, 0] == pytest.approx(15.0)
    assert np.isnan(image[1, 0])


def test_over_sampled_sky_marks_bad_pixels_and_slit_edges_before_clipping(
    log: object,
) -> None:
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_clipping_sigma": 4,
            "percentile_clipping_iterations": 2,
            "percentile_rolling_window_size": 9,
        }
    )
    pixels = pd.DataFrame(
        {
            "slit_position": [-1.0, 0.0, 1.0],
            "mask": [False, True, False],
        }
    )
    received: dict[str, object] = {}

    def record_clipping(**kwargs: object) -> pd.DataFrame:
        received.update(kwargs)
        return cast(pd.DataFrame, kwargs["imageMapOrderDF"])

    subtractor.rolling_window_clipping = record_clipping

    result = subtractor.get_over_sampled_sky_from_order(
        pixels,
        clipBPs=True,
        clipSlitEdge=0.25,
    )

    assert result["flagged_edge_clipped"].tolist() == [True, False, True]
    assert result["flagged_bad_pixel_clipped"].tolist() == [False, True, False]
    assert result["flagged_all_clipped"].tolist() == [True, True, True]
    assert received["windowSize"] == 9
    assert received["sigma_clip_limit"] == 4
    assert received["max_iterations"] == 2


def test_rolling_window_clipping_marks_an_isolated_bright_object_pixel(
    log: object,
) -> None:
    """Flag a bright synthetic object while retaining the surrounding sky pixels."""
    subtractor = _subtractor(log)
    subtractor.debug = False
    subtractor.stopSubtraction = False
    pixels = pd.DataFrame(
        {
            "order": [10] * 11,
            "flux": [10.0, 10.0, 10.0, 10.0, 10.0, 100.0, 10.0, 10.0, 10.0, 10.0, 10.0],
            "flagged_all_clipped": [False] * 11,
            "flagged_object_clipped": [False] * 11,
        }
    )

    result = subtractor.rolling_window_clipping(
        pixels,
        windowSize=5,
        sigma_clip_limit=2.0,
        max_iterations=3,
    )

    assert result.loc[5, "flagged_object_clipped"]
    assert result.loc[5, "flagged_all_clipped"]
    assert not result.loc[4, "flagged_object_clipped"]
    assert "residual_global_sigma" in result


def test_subtract_stops_without_products_when_order_sampling_requests_abort(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The public subtraction seam preserves its established early-abort result."""
    subtractor = _subtractor(log)
    subtractor.mapDF = pd.DataFrame(
        {
            "order": [10, 10, 11, 11],
            "mask": [False, False, False, False],
        }
    )
    subtractor.qc = pd.DataFrame({"qc_name": []})
    subtractor.products = pd.DataFrame({"product_label": []})
    subtractor.stopSubtraction = False
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "bspline_order": 3,
            "clip-slit-edge-fraction": 0.1,
        }
    )
    placeholder = object()
    calls: list[tuple[int, bool, float]] = []
    monkeypatch.setattr(
        subtractor,
        "create_placeholder_images",
        lambda: (placeholder, placeholder, placeholder),
    )

    def abort_after_first_order(
        imageMapOrder: pd.DataFrame,
        *,
        clipBPs: bool,
        clipSlitEdge: float,
    ) -> pd.DataFrame:
        calls.append((int(imageMapOrder["order"].iloc[0]), clipBPs, clipSlitEdge))
        subtractor.stopSubtraction = True
        return imageMapOrder

    monkeypatch.setattr(subtractor, "get_over_sampled_sky_from_order", abort_after_first_order)

    result = subtractor.subtract()

    assert result == (None, None, None, subtractor.qc, subtractor.products)
    assert calls == [(10, True, 0.1)]


@pytest.mark.parametrize("writeQCPlot", [False, True])
def test_subtract_assembles_modelled_order_and_optionally_registers_qc_plot(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    writeQCPlot: bool,
) -> None:
    """A fitted synthetic order is placed into output images and optional QC products."""
    subtractor = _subtractor(log)
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.debug = False
    subtractor.dateObs = "2024-01-02T00:00:00"
    subtractor.detectorParams = {"dispersion-axis": "x"}
    subtractor.objectFrame = CCDData(
        np.zeros((3, 3)),
        unit=u.electron,
        mask=np.zeros((3, 3), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((3, 3)), unit=u.electron),
    )
    subtractor.mapDF = pd.DataFrame(
        {"order": [10], "x": [1], "y": [2], "mask": [False]}
    )
    subtractor.qc = pd.DataFrame({"qc_name": []})
    subtractor.products = pd.DataFrame({"product_label": []})
    subtractor.stopSubtraction = False
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "aggressive_object_masking": False,
            "bspline_order": 3,
            "clip-slit-edge-fraction": 0.1,
            "sky_model_qc_plot": writeQCPlot,
        }
    )
    fittedOrder = pd.DataFrame(
        {
            "x": [1],
            "y": [2],
            "sky_model": [4.0],
            "sky_subtracted_flux": [6.0],
            "error": [2.0],
        }
    )
    monkeypatch.setattr(
        subtractor,
        "get_over_sampled_sky_from_order",
        lambda imageMapOrder, **kwargs: imageMapOrder,
    )
    monkeypatch.setattr(
        subtractor,
        "clip_object_slit_positions",
        lambda orders, **kwargs: orders,
    )
    monkeypatch.setattr(
        subtractor,
        "fit_bspline_curve_to_sky",
        lambda imageMapOrder: (fittedOrder, object(), [1, 2], np.array([2.0, 4.0]), 1.5),
    )
    plotCalls: list[tuple[CCDData, CCDData, CCDData]] = []
    monkeypatch.setattr(
        subtractor,
        "plot_image_comparison",
        lambda original, model, subtracted: (
            plotCalls.append((original, model, subtracted)) or "/qc/sky-model.pdf"
        ),
    )

    model, subtracted, residuals, qc, products = subtractor.subtract()

    assert model.data[2, 1] == 4.0
    assert subtracted.data[2, 1] == 6.0
    assert residuals.data[2, 1] == 3.0
    assert np.count_nonzero(model.data) == 1
    assert qc.empty
    assert products["product_label"].tolist() == (
        ["SKY SUBTRACTION QUICKLOOK"] if writeQCPlot else []
    )
    assert len(plotCalls) == int(writeQCPlot)


def test_object_clipping_marks_object_pixels_and_updates_local_noise(
    log: object,
) -> None:
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_rolling_window_size": 3,
            "noise_rolling_window_size": 30,
        }
    )
    pixels = pd.DataFrame(
        {
            "slit_position": [-1.0, -0.5, 0.0, 0.5, 1.0],
            "flagged_all_clipped": [False, False, False, False, False],
            "flagged_object_clipped": [False, False, True, False, False],
            "flux_minus_smoothed_residual": [1.0, 2.0, 3.0, 4.0, 5.0],
            "flux_percentile_smoothed": [2.0, 3.0, 4.0, 5.0, 6.0],
        }
    )

    result = subtractor.clip_object_slit_positions([pixels])

    assert result[0]["flagged_all_clipped"].tolist() == [
        False,
        False,
        True,
        False,
        False,
    ]
    assert np.isnan(result[0].loc[2, "residual_windowed_std"])
    assert result[0].loc[1, "residual_windowed_std"] == pytest.approx(
        np.sqrt(7 / 3)
    )
