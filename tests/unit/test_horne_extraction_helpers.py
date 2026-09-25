"""Analytic contracts for Horne-extraction helper calculations."""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.horne_extraction import (
    compute_extractions,
    extract_single_order,
    fit_object_profile,
    generate_masks,
    horne_extraction,
)
from tests.factories import instrument_header, product_table

pytestmark = pytest.mark.unit


def _extractor(log: object) -> horne_extraction:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    return extractor


def _configure_synthetic_orchestration(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    products: bool | pd.DataFrame,
) -> tuple[horne_extraction, pd.DataFrame, dict[str, object]]:
    """Build an extractor with deterministic collaborators for public orchestration tests."""
    import fundamentals

    import soxspipe.commonutils.toolkit as toolkit

    transformerModule = importlib.import_module(
        "soxspipe.commonutils.image_transformer"
    )
    frame = CCDData(
        np.full((2, 3), 10.0),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros((2, 3), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((2, 3)), unit=u.electron),
    )
    orderSlice = pd.DataFrame({"order": [10, 10]})
    extraction = pd.DataFrame(
        {
            "wavelengthMean": [500.0, 501.0, 502.0],
            "pixelScaleNm": [0.1, 4.0, 0.1],
            "extractedFluxOptimal": [1.0, 2.0, 3.0],
        }
    )
    merged = pd.DataFrame(
        {
            "WAVE": [500.0, 502.0],
            "FLUX_COUNTS": [1.0, 3.0],
            "SNR": [1.5, 2.5],
            "FLUX_DENSITY_COUNTS": [0.1, 0.3],
        }
    )
    captured: dict[str, object] = {}

    class FakeTransformer:
        """Minimal rectification collaborator with deterministic order slices."""

        def __init__(self, **kwargs: object) -> None:
            captured["transformer"] = kwargs

        def cache_image(self, name: str, image: np.ndarray, **kwargs: object) -> None:
            captured.setdefault("cached", []).append(name)

        def get_order_slices(self) -> list[pd.DataFrame]:
            return [orderSlice]

        def get_order_wavelength_ranges(self) -> list[tuple[float, float]]:
            return [(499.0, 503.0)]

        def get_order_rectified(self) -> list[dict[str, np.ndarray]]:
            return [{}]

    extractor = _extractor(log)
    extractor.orderPixelTable = pd.DataFrame({"order": [10]})
    extractor.qc = pd.DataFrame()
    extractor.products = products
    extractor.kw = lambda key: key
    extractor.arm = "VIS"
    extractor.twoDMap = {
        "WAVELENGTH": CCDData(np.ones((2, 3)), unit=u.nm),
        "SLIT": CCDData(np.zeros((2, 3)), unit=u.one),
    }
    extractor.skySubtractedFrame = frame
    extractor.subtractedFrame = None
    extractor.settings = {"instrument": "soxs"}
    extractor.twoDMapPath = "two-d-map.fits"
    extractor.dispersionMap = "dispersion-map.fits"
    extractor.slitHalfLength = 2
    extractor.debug = False
    extractor.turnOffMP = True
    extractor.clippingSigma = 3.0
    extractor.clippingIterationLimit = 2
    extractor.globalClippingSigma = 3.0
    extractor.axisA = "x"
    extractor.axisB = "y"
    extractor.recipeSettings = {"horne-extraction-profile-poly-order": 2}
    extractor.recipeName = "soxs-stare"
    extractor.dateObs = "2024-01-02T03:04:05"
    monkeypatch.setattr(transformerModule, "image_transformer", FakeTransformer)
    monkeypatch.setattr(fundamentals, "fmultiprocess", lambda **kwargs: [extraction])
    monkeypatch.setattr(
        toolkit,
        "add_snr_efficiency_qcs",
        lambda **kwargs: kwargs["qcTable"],
    )
    extractor.plot_extracted_spectrum_qc = lambda **kwargs: captured.update(kwargs)
    extractor.tune_wavelength_calibration_to_skylines = lambda orders, arm: orders
    extractor.merge_extracted_orders = lambda orders: (merged.copy(), {"1011": 501.0})
    return extractor, extraction, captured


def test_order_array_extraction_coerces_numeric_values_and_missing_object_flux(
    log: object,
) -> None:
    extractor = _extractor(log)
    order = pd.DataFrame(
        {
            "wavelength_shifted": [500.0, "invalid"],
            "skyFlux": [4.0, "invalid"],
        }
    )

    wavelength, sky, objectFlux = extractor._extract_order_arrays(order)

    np.testing.assert_allclose(wavelength, [500.0, np.nan], equal_nan=True)
    np.testing.assert_allclose(sky, [4.0, np.nan], equal_nan=True)
    assert objectFlux.isna().all()


def test_extract_orchestrates_synthetic_orders_without_writing_products(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Extract synthetic order slices while leaving product writing disabled."""
    extractor, extraction, captured = _configure_synthetic_orchestration(
        log, monkeypatch, products=False
    )

    qc, products, mergedSpectrum, joins, filePath = extractor.extract()

    assert qc.empty
    assert products is False
    assert filePath is False
    assert joins == {"1011": 501.0}
    assert mergedSpectrum["WAVE"].tolist() == [500.0, 502.0]
    assert captured["transformer"]
    assert captured["cached"] == ["fluxRaw", "variance"]
    assert captured["extractions"][0].equals(extraction)


def test_extract_returns_empty_results_when_no_order_trace_exists(log: object) -> None:
    """Return the compatibility tuple without collaborators when no trace was detected."""
    extractor = _extractor(log)
    extractor.orderPixelTable = None
    extractor.qc = pd.DataFrame({"qc_name": ["existing"]})
    extractor.products = product_table().iloc[0:0]

    result = extractor.extract()

    assert result[0].equals(extractor.qc)
    assert result[1].empty
    assert result[2:] == (None, None, None)


def test_constructor_prepares_vis_extraction_after_trace_detection(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Initialize the public extractor with the trace and order-table contracts intact."""
    import soxspipe.commonutils as commonutils
    import soxspipe.commonutils.toolkit as toolkit
    from soxspipe.commonutils.base_util import base_util

    header = instrument_header()
    header["INSTRUME"] = "SOXS"
    header["MJDOBS"] = 60311.0
    header["HIERARCH ESO SEQ CUMOFF Y"] = 1.0
    frame = CCDData(
        np.ones((3, 3)),
        unit=u.electron,
        meta=header,
        mask=np.zeros((3, 3), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((3, 3)), unit=u.electron),
    )
    orderPixels = pd.DataFrame(
        {
            "order": [10, 10],
            "xcoord_centre": [4.0, 5.0],
        }
    )
    captured: dict[str, object] = {}

    def initialize_base(
        self: horne_extraction,
        receivedLog: object,
        settings: dict[str, object],
        **_: object,
    ) -> None:
        self.log = receivedLog
        self.settings = settings
        self.dispersionMap = "dispersion-map.fits"
        self.binx = 1
        self.biny = 1
        self.detectorParams = {"dispersion-axis": "x"}
        self.imageMap = pd.DataFrame(
            {"wavelength": [500.0, 0.0], "slit_position": [0.0, 0.0]}
        )
        self.twoDMap = {
            "WAVELENGTH": CCDData(
                np.ones((3, 3)), unit=u.nm, meta={"MJDOBS": 60310.5}
            )
        }
        self.kw = lambda name: name
        self.arm = "VIS"
        self.axisA = "x"

    class TraceDetector:
        """Synthetic trace detector exposing the public result tuple."""

        def __init__(self, **kwargs: object) -> None:
            captured.update(kwargs)

        def get(self) -> tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            return (
                "products/TRACE.fits",
                pd.DataFrame({"qc_name": ["TRACE"]}),
                product_table().iloc[0:0],
                pd.DataFrame(),
                orderPixels.copy(),
                pd.DataFrame(),
            )

    monkeypatch.setattr(base_util, "__init__", initialize_base)
    monkeypatch.setattr(commonutils, "detect_continuum", TraceDetector)
    monkeypatch.setattr(
        toolkit,
        "utility_setup",
        lambda **_: ("products/qc", "products"),
    )
    monkeypatch.setattr(
        toolkit,
        "unpack_order_table",
        lambda **_: (pd.DataFrame(), orderPixels.copy(), pd.DataFrame()),
    )

    extractor = horne_extraction(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={
            "horne-extraction-slit-length": 8,
            "horne-extraction-profile-clipping-sigma": 3.0,
            "horne-extraction-profile-clipping-iteration-count": 2,
            "horne-extraction-profile-global-clipping-sigma": 4.0,
        },
        skySubtractedFrame=frame,
        unflattenedFrame=frame,
        twoDMapPath="two-d-map.fits",
        recipeName="soxs-stare",
        qcTable=pd.DataFrame(),
        productsTable=product_table().iloc[0:0],
        dispersionMap="dispersion-map.fits",
        sofName="SYNTHETIC",
        locationSetIndex=2,
        startNightDate="2024-01-02",
        turnOffMP=True,
    )

    assert extractor.noddingSequence == "_A2"
    assert extractor.filenameTemplate == "SYNTHETIC.fits"
    assert extractor.slitHalfLength == 4
    assert extractor.orderPixelTable["order"].tolist() == [10, 10]
    assert extractor.imageMap["wavelength"].tolist() == [500.0]
    assert captured["recipeName"] == "soxs-stare"
    assert captured["locationSetIndex"] == 2


def test_extract_writes_order_and_merged_product_contracts(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Write stable extracted-spectrum product metadata and an ASCII companion file."""
    import soxspipe.commonutils.phase3 as phase3

    extractor, _, _ = _configure_synthetic_orchestration(
        log, monkeypatch, products=product_table().iloc[0:0]
    )
    extractor.filenameTemplate = "SYNTHETIC.fits"
    extractor.productDir = str(tmp_path)
    extractor.noddingSequence = ""
    extractor.notFlattened = ""
    calls: list[dict[str, object]] = []
    monkeypatch.setattr(phase3, "write_fits_table_to_disk", lambda **kwargs: calls.append(kwargs))

    _, products, mergedSpectrum, joins, filePath = extractor.extract()

    assert len(calls) == 2
    assert calls[0]["filePath"].endswith("SYNTHETIC_EXTRACTED_ORDERS.fits")
    assert calls[1]["filePath"].endswith("SYNTHETIC_EXTRACTED_MERGED.fits")
    assert calls[1]["header"]["PRO_TYPE"] == "REDUCED"
    assert calls[1]["header"]["PRO_CATG"] == "SCI_SLIT_FLUX_VIS"
    assert products["product_label"].tolist() == [
        "EXTRACTED_ORDERS_TABLE",
        "EXTRACTED_MERGED_ASCII",
        "EXTRACTED_MERGED_TABLE",
    ]
    assert (tmp_path / "SYNTHETIC_EXTRACTED_MERGED.txt").is_file()
    assert filePath == str(tmp_path / "SYNTHETIC_EXTRACTED_MERGED.fits")
    assert mergedSpectrum["WAVE"].tolist() == [500.0, 502.0]
    assert joins == {"1011": 501.0}


def test_plot_extracted_spectrum_qc_writes_pdf_and_product_record(
    log: object,
    tmp_path: Path,
) -> None:
    """Write the public extracted-spectrum diagnostic with its product contract."""
    extractor = _extractor(log)
    wavelengths = np.linspace(500.0, 523.0, 24)
    extractor.arm = "VIS"
    extractor.settings = {"PAE": False}
    extractor.skylinesDF = pd.DataFrame({"WAVELENGTH": [510.0], "FLUX": [3.0]})
    extractor.filenameTemplate = "SYNTHETIC.fits"
    extractor.noddingSequence = "_AB"
    extractor.notFlattened = ""
    extractor.qcDir = str(tmp_path)
    extractor.dateObs = "2024-01-02T03:04:05"
    extractor.products = product_table().iloc[0:0]
    extractions = [
        pd.DataFrame(
            {
                "order": [10] * len(wavelengths),
                "wavelengthMean": wavelengths,
                "extractedFluxOptimal": np.linspace(2.0, 4.0, len(wavelengths)),
                "extractedFluxBoxcarRobust": np.linspace(1.0, 3.0, len(wavelengths)),
                "varianceSpectrum": np.ones(len(wavelengths)),
                "skyFlux": np.linspace(0.2, 0.5, len(wavelengths)),
                "snr": np.linspace(3.0, 6.0, len(wavelengths)),
            }
        )
    ]

    result = extractor.plot_extracted_spectrum_qc(extractions)

    expectedPath = tmp_path / "SYNTHETIC_EXTRACTED_ORDERS_QC_PLOT_AB.pdf"
    assert result is None
    assert expectedPath.is_file()
    assert expectedPath.read_bytes().startswith(b"%PDF")
    assert extractor.products["product_label"].tolist() == [
        "EXTRACTED_ORDERS_QC_PLOT_AB"
    ]
    assert extractor.products.loc[0, "file_path"] == str(expectedPath)


def test_skyline_matching_and_clipped_shift_reject_outliers(log: object) -> None:
    extractor = _extractor(log)
    matchedPixels, matchedShifts = extractor._match_peaks_to_skylines(
        np.array([100.0, 101.0, 200.0]),
        np.array([0, 2]),
        np.array([100.2, 199.8]),
        tolerance=0.5,
    )
    clippedWave, clippedShifts, medianShift = extractor._compute_clipped_shift(
        list(range(21)),
        [0.1] * 20 + [10.0],
    )

    assert matchedPixels == [100.2, 199.8]
    assert matchedShifts == pytest.approx([0.2, -0.2])
    assert clippedWave.tolist() == [20]
    assert clippedShifts.tolist() == [10.0]
    assert medianShift == pytest.approx(0.1)


def test_local_skylines_filters_wavelengths_and_projects_single_order(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    extractor = _extractor(log)
    extractor.dispersionMap = "synthetic-dispersion-map.fits"
    extractor.skylinesDF = pd.DataFrame(
        {
            "WAVELENGTH": [499.0, 500.0, 501.0, 502.0],
            "ISOLATED": [True, False, True, True],
        }
    )
    received: dict[str, object] = {}

    def project_skylines(**kwargs: object) -> pd.DataFrame:
        received.update(kwargs)
        projected = kwargs["orderPixelTable"].copy()
        projected["fit_x"] = projected["wavelength"] - 400.0
        return projected

    import soxspipe.commonutils as commonutils

    monkeypatch.setattr(commonutils, "dispersion_map_to_pixel_arrays", project_skylines)

    allSkylines, allCalibration = extractor._get_local_skylines_for_order(
        500.0, 501.0, "all", "wavelength"
    )
    orderSkylines, orderCalibration = extractor._get_local_skylines_for_order(
        500.0, 501.0, 12, "fit_x"
    )

    assert allSkylines.tolist() == [500.0, 501.0]
    assert allCalibration.tolist() == [501.0]
    assert orderSkylines.tolist() == [100.0, 101.0]
    assert orderCalibration.tolist() == [101.0]
    assert received["dispersionMapPath"] == "synthetic-dispersion-map.fits"
    assert received["orderPixelTable"]["order"].tolist() == [12, 12]


def test_order_shift_qc_records_rounded_pixel_correction(log: object) -> None:
    extractor = _extractor(log)
    extractor.recipeName = "soxs-stare"
    extractor.dateObs = "2026-09-10T12:00:00"
    extractor.qc = pd.DataFrame()

    extractor._record_order_shift_qc(12, 0.12349)

    assert extractor.qc.loc[0, "qc_name"] == "SKY SHIFT O12"
    assert extractor.qc.loc[0, "qc_value"] == pytest.approx(0.123)
    assert extractor.qc.loc[0, "qc_unit"] == "pixels"
    # NUMPY BOOL, NOT PYTHON BOOL, NOW THE COLUMN CARRIES A REAL BOOL DTYPE
    assert bool(extractor.qc.loc[0, "to_header"]) is True


def test_sky_peak_detection_retains_original_flux_and_peak_coordinates(
    log: object,
) -> None:
    extractor = _extractor(log)
    wavelength = pd.Series(np.arange(500.0, 521.0))
    sky = pd.Series(np.ones(21))
    sky.iloc[10] = 20.0
    objectFlux = pd.Series(np.arange(21.0))

    originalSky, normalisedSky, peaks, peakWavelengths, peakObjectFlux = (
        extractor._detect_sky_peaks(wavelength, sky, objectFlux, sky.notna())
    )

    np.testing.assert_allclose(originalSky, sky)
    assert peaks.tolist() == [10]
    assert normalisedSky[10] > normalisedSky[0]
    assert peakWavelengths[10] == 510.0
    assert peakObjectFlux[10] == 10.0


def test_order_merge_resamples_nir_flux_and_preserves_variance(log: object) -> None:
    extractor = _extractor(log)
    extractor.arm = "NIR"
    extractor.kw = lambda key: key
    extractedOrders = pd.DataFrame(
        {
            "order": [10, 10, 10],
            "wavelengthMean": [1000.00, 1000.06, 1000.12],
            "pixelScaleNm": [0.06, 0.06, 0.06],
            "extractedFluxOptimal": [10.0, 16.0, 22.0],
            "extractedFluxBoxcar": [20.0, 30.0, 40.0],
            "extractedFluxBoxcarRobust": [20.0, 30.0, 40.0],
            "varianceSpectrum": [4.0, 4.0, 4.0],
            "skyFlux": [1.0, 2.0, 3.0],
        }
    )

    merged, joins = extractor.merge_extracted_orders(extractedOrders)

    np.testing.assert_allclose(merged["WAVE"].values.value, [1000.02, 1000.08])
    np.testing.assert_allclose(merged["FLUX_COUNTS"].values.value, [12.0, 18.0])
    np.testing.assert_allclose(merged["VARIANCE"].values.value, [4.0, 4.0])
    np.testing.assert_allclose(merged["SKY_COUNTS"].values.value, [4.0 / 3, 7.0 / 3])
    np.testing.assert_allclose(merged["SNR"].values.value, [6.0, 9.0])
    assert joins == {}


def test_order_merge_helpers_weight_signal_and_choose_lower_residual(
    log: object,
) -> None:
    extractor = _extractor(log)
    overlap = pd.DataFrame(
        {
            "flux_resampled": [10.0, -4.0],
            "snr": [2.0, 1.0],
        }
    )

    weighted = extractor.weighted_average(overlap.copy())
    selected = extractor.residual_merge(overlap.copy())

    assert weighted.tolist() == pytest.approx([16.0 / 3.0, 16.0 / 3.0])
    assert selected.tolist() == pytest.approx([-4.0, -4.0])


def test_wavelength_tuning_uses_order_four_shift_for_unmatched_vis_order(
    log: object,
) -> None:
    extractor = _extractor(log)
    extractor.debug = False
    extraction = pd.DataFrame(
        {
            "order": [4, 11],
            "wavelengthMean": [500.0, 600.0],
            "pixelScaleNm": [0.1, 0.1],
        }
    )
    measuredShifts = iter([0.1, 0.1, 0.1, 0.0, 0.0, 0.0])
    recordedOrders: list[tuple[int, float]] = []

    extractor._extract_order_arrays = lambda order: (
        pd.Series([1.0]),
        pd.Series([1.0]),
        pd.Series([1.0]),
    )
    extractor._detect_sky_peaks = lambda *args: (
        np.array([1.0]),
        np.array([1.0]),
        np.array([0]),
        np.array([1.0]),
        np.array([1.0]),
    )
    extractor._get_local_skylines_for_order = lambda *args: (
        np.array([1.0]),
        np.array([1.0]),
    )
    extractor._match_peaks_to_skylines = lambda *args, **kwargs: ([1.0], [0.0])
    extractor._compute_clipped_shift = lambda *args: (
        np.array([]),
        np.array([]),
        next(measuredShifts),
    )
    extractor._record_order_shift_qc = lambda order, shift: recordedOrders.append(
        (order, shift)
    )

    result = extractor.tune_wavelength_calibration_to_skylines(extraction, "VIS")

    assert result["wavelengthMean"].tolist() == pytest.approx([500.3, 600.3])
    assert recordedOrders == [(4, pytest.approx(0.3)), (11, pytest.approx(0.0))]


def test_compute_extractions_returns_sorted_optimal_and_boxcar_spectra() -> None:
    crossDispersionSlices = pd.DataFrame(
        {
            "pixelScaleNm": [1.0, 1.0, 1.0],
        }
    )
    rectifiedImages = {
        "fluxRaw": np.array([[2.0, 4.0, 6.0], [4.0, 6.0, 8.0]]),
        "objectProfile": np.full((2, 3), 0.5),
        "variance": np.ones((2, 3)),
        "mask": np.zeros((2, 3), dtype=bool),
        "wavelength": np.array([[502.0, 501.0, 500.0], [502.0, 501.0, 500.0]]),
        "fluxSky": np.full((2, 3), 2.0),
    }

    result = compute_extractions(crossDispersionSlices, rectifiedImages, order=10)

    assert result["wavelengthMean"].tolist() == [500.0, 501.0, 502.0]
    assert result["extractedFluxOptimal"].tolist() == [14.0, 10.0, 6.0]
    assert result["extractedFluxBoxcar"].tolist() == [14.0, 10.0, 6.0]
    assert result["skyFlux"].tolist() == [2.0, 2.0, 2.0]
    np.testing.assert_allclose(result["varianceSpectrum"], 2.0)


def test_mask_generation_calculates_pixel_scale_and_propagates_bad_pixels() -> None:
    slices = pd.DataFrame(index=range(3))
    images = {
        "bpMask": np.array([[False, True, False], [False, False, False]]),
        "fluxRaw": np.ones((2, 3)),
        "wavelength": np.array([[500.0, 501.0, 502.0], [500.0, 501.0, 502.0]]),
    }

    resultSlices, resultImages = generate_masks(slices, images)

    np.testing.assert_allclose(
        resultSlices["pixelScaleNm"], [np.nan, 1.0, np.nan], equal_nan=True
    )
    assert resultImages["mask"].tolist() == [
        [False, True, False],
        [False, False, False],
    ]
    assert resultSlices.loc[1, "mask"].tolist() == [True, False]


def test_profile_fitting_normalises_a_symmetric_slit_profile() -> None:
    slices = pd.DataFrame(index=range(5))
    images = {
        "fluxRaw": np.array(
            [[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0]]
        ),
        "mask": np.zeros((2, 5), dtype=bool),
    }

    resultSlices, resultImages = fit_object_profile(
        slices,
        images,
        slitHalfLength=1,
        clippingSigma=3.0,
        clippingIterationLimit=2,
        hornePolyOrder=0,
        axisB="x",
        order=10,
        debug=False,
        plt=None,
    )

    np.testing.assert_allclose(resultImages["objectProfile"], 0.5)
    np.testing.assert_allclose(np.stack(resultSlices["objectProfile"]), 0.5)


def test_single_order_extraction_returns_sorted_science_columns(log: object) -> None:
    slices = pd.DataFrame({"order": [10] * 5})
    images = {
        "bpMask": np.zeros((2, 5), dtype=bool),
        "fluxRaw": np.array(
            [[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0]]
        ),
        "variance": np.ones((2, 5)),
        "wavelength": np.array(
            [[504.0, 503.0, 502.0, 501.0, 500.0], [504.0, 503.0, 502.0, 501.0, 500.0]]
        ),
    }

    result = extract_single_order(
        (slices, images),
        log,
        slitHalfLength=1,
        clippingSigma=3.0,
        clippingIterationLimit=2,
        globalClippingSigma=3.0,
        axisA="y",
        axisB="x",
        hornePolyOrder=0,
    )

    assert result is not None
    assert result.columns.tolist() == [
        "order",
        "wavelengthMean",
        "pixelScaleNm",
        "varianceSpectrum",
        "snr",
        "extractedFluxOptimal",
        "extractedFluxBoxcar",
        "extractedFluxBoxcarRobust",
        "skyFlux",
    ]
    assert result["wavelengthMean"].tolist() == [501.0, 502.0, 503.0]
    assert result["extractedFluxOptimal"].tolist() == [8.0, 6.0, 4.0]
    np.testing.assert_allclose(result["varianceSpectrum"], 2.0)
    assert result["skyFlux"].isna().all()
