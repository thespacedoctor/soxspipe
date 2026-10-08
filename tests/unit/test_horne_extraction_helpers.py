"""Analytic contracts for Horne-extraction helper calculations."""

from __future__ import annotations

import importlib
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.horne_extraction import (
    _slit_edge_truncated_columns,
    _unsupported_pixel_mask,
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

    transformerModule = importlib.import_module("soxspipe.commonutils.image_transformer")
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
            "order": [10, 10, 10],
            "wavelengthMean": [500.0, 501.0, 502.0],
            "pixelScaleNm": [0.1, 4.0, 0.1],
            "extractedFluxOptimal": [1.0, 2.0, 3.0],
            "slitEdgeTruncated": [False, False, False],
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
            captured["transformerInstance"] = self

        def cache_image(self, name: str, image: np.ndarray, **kwargs: object) -> None:
            captured.setdefault("cached", []).append(name)
            captured.setdefault("cacheKwargs", {})[name] = kwargs

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
    extractor.notFlattened = ""
    monkeypatch.setattr(transformerModule, "image_transformer", FakeTransformer)
    monkeypatch.setattr(fundamentals, "fmultiprocess", lambda **kwargs: [extraction])
    monkeypatch.setattr(
        toolkit,
        "add_snr_efficiency_qcs",
        lambda **kwargs: kwargs["qcTable"],
    )
    extractor.plot_extracted_spectrum_qc = lambda **kwargs: captured.update(kwargs)
    extractor.plot_slit_drift_qc = lambda transformer: captured.update(slitDriftTransformer=transformer)
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
    extractor, extraction, captured = _configure_synthetic_orchestration(log, monkeypatch, products=False)

    qc, products, mergedSpectrum, joins, filePath = extractor.extract()

    assert qc["qc_name"].tolist() == ["N ORDERS SLIT EDGE"]
    assert qc.loc[0, "qc_value"] == 0
    assert products is False
    assert filePath is False
    assert joins == {"1011": 501.0}
    assert mergedSpectrum["WAVE"].tolist() == [500.0, 502.0]
    assert captured["transformer"]
    assert captured["cached"] == ["fluxRaw", "variance"]
    assert captured["cacheKwargs"]["variance"]["associatedMask"] is extractor.skySubtractedFrame.mask
    assert captured["slitDriftTransformer"] is captured["transformerInstance"]
    assert captured["extractions"][0].equals(extraction.drop(columns=["slitEdgeTruncated"]))
    assert extractor.slitEdgeOrders == []


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


def _construct_extractor(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    *,
    headerOverrides: dict[str, object],
    recipeName: str,
    locationSetIndex: object,
    slitLength: float = 8,
    binning: int = 1,
) -> tuple[horne_extraction, dict[str, object]]:
    """Build an extractor through its real constructor, with synthetic trace detection."""
    import soxspipe.commonutils as commonutils
    import soxspipe.commonutils.toolkit as toolkit
    from soxspipe.commonutils.base_util import base_util

    header = instrument_header()
    header["INSTRUME"] = "SOXS"
    header["MJDOBS"] = 60311.0
    for keyword, value in headerOverrides.items():
        header[keyword] = value
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
        self.binx = binning
        self.biny = binning
        self.detectorParams = {"dispersion-axis": "x"}
        self.imageMap = pd.DataFrame({"wavelength": [500.0, 0.0], "slit_position": [0.0, 0.0]})
        self.twoDMap = {"WAVELENGTH": CCDData(np.ones((3, 3)), unit=u.nm, meta={"MJDOBS": 60310.5})}
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
            "horne-extraction-slit-length": slitLength,
            "horne-extraction-profile-clipping-sigma": 3.0,
            "horne-extraction-profile-clipping-iteration-count": 2,
            "horne-extraction-profile-global-clipping-sigma": 4.0,
        },
        skySubtractedFrame=frame,
        unflattenedFrame=frame,
        twoDMapPath="two-d-map.fits",
        recipeName=recipeName,
        qcTable=pd.DataFrame(),
        productsTable=product_table().iloc[0:0],
        dispersionMap="dispersion-map.fits",
        sofName="SYNTHETIC",
        locationSetIndex=locationSetIndex,
        startNightDate="2024-01-02",
        turnOffMP=True,
    )

    return extractor, captured


def test_constructor_prepares_vis_extraction_after_trace_detection(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Initialize the public extractor with the trace and order-table contracts intact."""
    extractor, captured = _construct_extractor(
        log,
        monkeypatch,
        headerOverrides={"HIERARCH ESO SEQ CUMOFF Y": 1.0},
        recipeName="soxs-stare",
        locationSetIndex=2,
    )

    assert extractor.noddingSequence == "_A2"
    assert extractor.filenameTemplate == "SYNTHETIC.fits"
    assert extractor.slitHalfLength == 4
    assert extractor.orderPixelTable["order"].tolist() == [10, 10]
    assert extractor.imageMap["wavelength"].tolist() == [500.0]
    assert captured["recipeName"] == "soxs-stare"
    assert captured["locationSetIndex"] == 2


def test_odd_slit_length_setting_keeps_its_half_pixel_in_the_half_length(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep half of an odd slit-length setting as x.5 so the extracted slit is as long as the setting."""
    extractor, _ = _construct_extractor(
        log,
        monkeypatch,
        headerOverrides={},
        recipeName="soxs-stare",
        locationSetIndex=1,
        slitLength=15,
    )

    assert extractor.slitHalfLength == 7.5


def test_binned_slit_half_length_is_scaled_without_rounding(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Scale the slit half length by the binning and leave whole-row rounding to the image transformer."""
    extractor, _ = _construct_extractor(
        log,
        monkeypatch,
        headerOverrides={},
        recipeName="soxs-stare",
        locationSetIndex=1,
        slitLength=15,
        binning=2,
    )

    assert extractor.slitHalfLength == 3.75


@pytest.mark.parametrize(
    ("locationSetIndex", "expectedSequence"),
    [(1, "_1"), (2, "_2"), ("STACK", "_STACK")],
)
def test_offset_extraction_without_cumoff_y_is_suffixed_by_its_location_set(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    locationSetIndex: object,
    expectedSequence: str,
) -> None:
    """Real offset frames carry no CUMOFF Y, so each cycle and the stack take their own product suffix."""
    extractor, captured = _construct_extractor(
        log,
        monkeypatch,
        headerOverrides={"HIERARCH ESO SEQ FIXOFF RA": -2.0, "HIERARCH ESO SEQ FIXOFF DEC": 0.0},
        recipeName="soxs-offset",
        locationSetIndex=locationSetIndex,
    )

    assert extractor.noddingSequence == expectedSequence
    assert captured["locationSetIndex"] == locationSetIndex


def test_constructor_names_products_without_a_sof_name(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Use the input frame to name products when no SOF name is available."""
    import soxspipe.commonutils as commonutils
    import soxspipe.commonutils.toolkit as toolkit
    from soxspipe.commonutils.base_util import base_util

    extraction_module = importlib.import_module("soxspipe.commonutils.horne_extraction")

    frame = CCDData(
        np.ones((3, 3)),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros((3, 3), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((3, 3)), unit=u.electron),
    )
    named_frames: list[CCDData] = []

    def initialize_base(
        self: horne_extraction,
        receivedLog: object,
        settings: dict[str, object],
        **_: object,
    ) -> None:
        self.log = receivedLog
        self.settings = settings
        self.dispersionMap = None
        self.binx = 1
        self.biny = 1
        self.detectorParams = {"dispersion-axis": "x"}
        self.imageMap = pd.DataFrame({"wavelength": [500.0], "slit_position": [0.0]})

    def name_frame(**kwargs: object) -> str:
        named_frames.append(kwargs["frame"])
        return "OBJECT_STARE.fits"

    class MissingTrace:
        def __init__(self, **_: object) -> None:
            pass

        def get(self) -> tuple[None, None, None, None, None, None]:
            return (None, None, None, None, None, None)

    monkeypatch.setattr(base_util, "__init__", initialize_base)
    monkeypatch.setattr(extraction_module, "filenamer", name_frame)
    monkeypatch.setattr(commonutils, "detect_continuum", MissingTrace)
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: ("products/qc", "products"))

    extractor = horne_extraction(
        log=log,
        settings={},
        recipeSettings={
            "horne-extraction-slit-length": 20,
            "horne-extraction-profile-clipping-sigma": 3,
            "horne-extraction-profile-clipping-iteration-count": 5,
            "horne-extraction-profile-global-clipping-sigma": 25,
        },
        skySubtractedFrame=frame,
        unflattenedFrame=frame,
        twoDMapPath=None,
        sofName=False,
    )

    assert extractor.filenameTemplate == "OBJECT_STARE.fits"
    assert named_frames == [frame]


def test_extract_writes_order_and_merged_product_contracts(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Write stable extracted-spectrum product metadata and an ASCII companion file."""
    import soxspipe.commonutils.phase3 as phase3

    extractor, _, _ = _configure_synthetic_orchestration(log, monkeypatch, products=product_table().iloc[0:0])
    extractor.filenameTemplate = "SYNTHETIC.fits"
    extractor.productDir = str(tmp_path)
    extractor.noddingSequence = ""
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


@pytest.mark.parametrize("recipeName", ["soxs-stare", "soxs-nod", "soxs-offset"])
def test_extract_labels_every_product_row_with_the_running_recipe(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    recipeName: str,
) -> None:
    """Record the recipe the extraction was built with on each product row it adds."""
    import soxspipe.commonutils.phase3 as phase3

    extractor, _, _ = _configure_synthetic_orchestration(log, monkeypatch, products=product_table().iloc[0:0])
    extractor.recipeName = recipeName
    extractor.filenameTemplate = "SYNTHETIC.fits"
    extractor.productDir = str(tmp_path)
    extractor.noddingSequence = ""
    monkeypatch.setattr(phase3, "write_fits_table_to_disk", lambda **kwargs: None)

    _, products, _, _, _ = extractor.extract()

    assert products["product_label"].tolist() == [
        "EXTRACTED_ORDERS_TABLE",
        "EXTRACTED_MERGED_ASCII",
        "EXTRACTED_MERGED_TABLE",
    ]
    assert products["soxspipe_recipe"].tolist() == [recipeName] * 3


@pytest.mark.parametrize("recipeName", ["soxs-stare", "soxs-nod", "soxs-offset"])
def test_plot_extracted_spectrum_qc_labels_product_row_with_the_running_recipe(
    log: object,
    tmp_path: Path,
    recipeName: str,
) -> None:
    """Record the recipe the extraction was built with on the QC plot row."""
    extractor = _extractor(log)
    wavelengths = np.linspace(500.0, 523.0, 24)
    extractor.recipeName = recipeName
    extractor.arm = "VIS"
    extractor.settings = {"PAE": False}
    extractor.skylinesDF = pd.DataFrame({"WAVELENGTH": [510.0], "FLUX": [3.0]})
    extractor.filenameTemplate = "SYNTHETIC.fits"
    extractor.noddingSequence = ""
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

    extractor.plot_extracted_spectrum_qc(extractions)

    assert extractor.products["product_label"].tolist() == ["EXTRACTED_ORDERS_QC_PLOT"]
    assert extractor.products["soxspipe_recipe"].tolist() == [recipeName]


def test_plot_extracted_spectrum_qc_writes_pdf_and_product_record(
    log: object,
    tmp_path: Path,
) -> None:
    """Write the public extracted-spectrum diagnostic with its product contract."""
    extractor = _extractor(log)
    wavelengths = np.linspace(500.0, 523.0, 24)
    extractor.recipeName = "soxs-nod"
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
    assert extractor.products["product_label"].tolist() == ["EXTRACTED_ORDERS_QC_PLOT_AB"]
    assert extractor.products.loc[0, "file_path"] == str(expectedPath)


def _slit_drift_extractor(
    log: object,
    tmp_path: Path,
    products: bool | pd.DataFrame,
    notFlattened: str = "",
) -> tuple[horne_extraction, SimpleNamespace]:
    """Build an extractor and a transformer stand-in carrying the slit-drift attributes."""
    extractor = _extractor(log)
    extractor.filenameTemplate = "SYNTHETIC.fits"
    extractor.noddingSequence = "_A1"
    extractor.notFlattened = notFlattened
    extractor.qcDir = str(tmp_path)
    extractor.dateObs = "2024-01-02T03:04:05"
    extractor.recipeName = "soxs-nod"
    extractor.products = products
    extractor.arm = "VIS"
    orderPixelTable = pd.DataFrame(
        {
            "order": [10, 10, 10, 11],
            "wavelength": [600.0, 601.0, 602.0, 500.0],
            "slit_position": [1.1, 1.9, 3.2, 0.5],
        }
    )
    transformer = SimpleNamespace(
        orderPixelTable=orderPixelTable,
        uniqueOrders=[10, 11],
        orderSlitCentreCoeffs=[np.array([1.0, -599.0]), np.array([0.25])],
        wlMinMax=[(600.0, 602.0), (500.0, 501.0)],
        orderSlitCentreFallback=[False, True],
        globalSlitCentreArcsec=0.25,
    )
    return extractor, transformer


def test_plot_slit_drift_qc_writes_pdf_and_product_record(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor, transformer = _slit_drift_extractor(log, tmp_path, product_table().iloc[0:0])

    # ACT
    result = extractor.plot_slit_drift_qc(transformer)

    # ASSERT
    expectedPath = tmp_path / "SYNTHETIC_SLIT_DRIFT_QC_PLOT_A1.pdf"
    assert result is None
    assert expectedPath.read_bytes().startswith(b"%PDF")
    assert extractor.products["product_label"].tolist() == ["SLIT_DRIFT_QC_PLOT_A1"]
    assert extractor.products.loc[0, "file_path"] == str(expectedPath)
    assert extractor.products.loc[0, "file_name"] == expectedPath.name
    assert extractor.products.loc[0, "file_type"] == "PDF"
    assert extractor.products.loc[0, "label"] == "QC"
    assert extractor.products.loc[0, "soxspipe_recipe"] == "soxs-nod"


def test_plot_slit_drift_qc_writes_file_but_no_row_when_products_disabled(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor, transformer = _slit_drift_extractor(log, tmp_path, products=False)

    # ACT
    extractor.plot_slit_drift_qc(transformer)

    # ASSERT
    assert (tmp_path / "SYNTHETIC_SLIT_DRIFT_QC_PLOT_A1.pdf").is_file()
    assert extractor.products is False


def test_plot_slit_drift_qc_writes_nothing_for_the_not_flattened_re_extraction(log: object, tmp_path: Path) -> None:
    # ARRANGE
    products = product_table().iloc[0:0]
    extractor, transformer = _slit_drift_extractor(log, tmp_path, products, notFlattened="_NOTFLAT")

    # ACT
    extractor.plot_slit_drift_qc(transformer)

    # ASSERT
    assert list(tmp_path.glob("*.pdf")) == []
    assert extractor.products.empty


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

    allSkylines, allCalibration = extractor._get_local_skylines_for_order(500.0, 501.0, "all", "wavelength")
    orderSkylines, orderCalibration = extractor._get_local_skylines_for_order(500.0, 501.0, 12, "fit_x")

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

    originalSky, normalisedSky, peaks, peakWavelengths, peakObjectFlux = extractor._detect_sky_peaks(
        wavelength, sky, objectFlux, sky.notna()
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


def test_order_merge_renumbers_soxs_vis_orders_in_place_and_keeps_their_dtype(log: object) -> None:
    extractor = _extractor(log)
    extractor.arm = "VIS"
    extractor.kw = lambda key: key
    traceOrders = [1, 2, 3, 4]
    starts = {4: 700.0, 1: 600.0, 3: 500.0, 2: 400.0}
    rows = []
    for order in traceOrders:
        for step in range(3):
            rows.append(
                {
                    "order": order,
                    "wavelengthMean": starts[order] + 0.02 * step,
                    "pixelScaleNm": 0.02,
                    "extractedFluxOptimal": 10.0,
                    "extractedFluxBoxcar": 10.0,
                    "extractedFluxBoxcarRobust": 10.0,
                    "varianceSpectrum": 1.0,
                    "skyFlux": 1.0,
                }
            )
    extractedOrders = pd.DataFrame(rows).astype({"order": np.int16})

    extractor.merge_extracted_orders(extractedOrders)

    assert extractedOrders["order"].dtype == np.int16
    renumbered = extractedOrders.groupby("order")["wavelengthMean"].min().to_dict()
    assert renumbered == {1: 700.0, 2: 600.0, 3: 500.0, 4: 400.0}


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
    extractor._record_order_shift_qc = lambda order, shift: recordedOrders.append((order, shift))

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

    np.testing.assert_allclose(resultSlices["pixelScaleNm"], [np.nan, 1.0, np.nan], equal_nan=True)
    assert resultImages["mask"].tolist() == [
        [False, True, False],
        [False, False, False],
    ]
    assert resultSlices.loc[1, "mask"].tolist() == [True, False]


def test_profile_fitting_normalises_a_symmetric_slit_profile() -> None:
    slices = pd.DataFrame(index=range(5))
    images = {
        "fluxRaw": np.array([[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0]]),
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


def _fit_profile(images: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Fit a constant-in-dispersion object profile to rectified order images."""
    slices = pd.DataFrame(index=range(images["fluxRaw"].shape[1]))
    _, fitted = fit_object_profile(
        slices,
        images,
        slitHalfLength=images["fluxRaw"].shape[0] // 2,
        clippingSigma=3.0,
        clippingIterationLimit=5,
        hornePolyOrder=0,
        axisB="x",
        order=10,
        debug=False,
        plt=None,
    )
    return fitted


def _edge_order_images(
    onOrderColumns: int, columns: int = 40, rows: int = 13, offOrderRows: int = 7
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Build an order whose Gaussian object leaves the order after `onOrderColumns`.

    Returns the rectified images and the true on-slit flux per column.
    """
    rowIndex = np.arange(rows)[:, np.newaxis]
    shape = np.exp(-0.5 * ((rowIndex - 3.0) / 1.5) ** 2)
    columnFlux = 1000.0 + 10.0 * np.arange(columns)
    flux = shape * columnFlux[np.newaxis, :]
    offOrder = np.zeros((rows, columns), dtype=bool)
    offOrder[:offOrderRows, onOrderColumns:] = True
    flux[offOrder] = np.nan
    visibleFlux = shape[offOrderRows:, 0].sum() * columnFlux
    images = {
        "fluxRaw": flux,
        "mask": offOrder.copy(),
        "variance": np.where(offOrder, 1e12, 1.0),
        "wavelength": np.tile(np.arange(columns, dtype=float) + 500.0, (rows, 1)),
    }
    return images, visibleFlux


def test_profile_fitting_excludes_off_order_pixels_from_the_extracted_flux() -> None:
    images, visibleFlux = _edge_order_images(onOrderColumns=2)
    slices = pd.DataFrame({"pixelScaleNm": np.ones(40)})

    fitted = _fit_profile(images)
    result = compute_extractions(slices, fitted, order=10)

    edgeColumns = result["wavelengthMean"] >= 502.0
    assert edgeColumns.sum() == 38
    assert np.isfinite(result["extractedFluxOptimal"]).all()
    np.testing.assert_allclose(
        result.loc[edgeColumns, "extractedFluxOptimal"],
        visibleFlux[2:],
        rtol=1e-2,
    )


def test_profile_fitting_gives_zero_weight_to_off_order_pixels() -> None:
    images, _ = _edge_order_images(onOrderColumns=2)

    fitted = _fit_profile(images)

    assert (fitted["objectProfile"][images["mask"]] == 0).all()
    np.testing.assert_allclose(fitted["objectProfile"].sum(axis=0), 1.0)


def test_profile_fitting_drops_a_column_with_no_on_order_pixels() -> None:
    images, _ = _edge_order_images(onOrderColumns=2, offOrderRows=13)
    slices = pd.DataFrame({"pixelScaleNm": np.ones(40)})

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        fitted = _fit_profile(images)
    result = compute_extractions(slices, fitted, order=10)

    assert np.isnan(fitted["objectProfile"][:, 2:]).all()
    assert result["wavelengthMean"].tolist() == [500.0, 501.0]


def test_profile_fitting_keeps_weight_on_a_masked_pixel_with_finite_flux() -> None:
    rowIndex = np.arange(7)[:, np.newaxis]
    shape = np.exp(-0.5 * ((rowIndex - 3.0) / 1.0) ** 2)
    flux = shape * (100.0 + np.arange(10.0))[np.newaxis, :]
    mask = np.zeros(flux.shape, dtype=bool)
    mask[3, 4] = True
    images = {"fluxRaw": flux, "mask": mask}

    fitted = _fit_profile(images)

    profile = fitted["objectProfile"]
    assert profile[3, 4] > 0
    assert profile[3, 4] == pytest.approx(profile[3, 5])
    np.testing.assert_allclose(profile.sum(axis=0), 1.0)


def _centred_window_images(
    rows: int = 13, columns: int = 80, crhPixels: tuple[tuple[int, int], ...] = ()
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Build a well-sampled, centred Gaussian object window with masked cosmic-ray hits.

    Returns the rectified images and the true normalised cross-slit profile.
    """
    rowIndex = np.arange(rows)[:, np.newaxis]
    shape = np.exp(-0.5 * ((rowIndex - rows // 2) / 1.5) ** 2)
    flux = shape * (1000.0 + 5.0 * np.arange(columns))[np.newaxis, :]
    mask = np.zeros((rows, columns), dtype=bool)
    for row, column in crhPixels:
        mask[row, column] = True
        flux[row, column] *= 50.0
    return {"fluxRaw": flux, "mask": mask}, shape[:, 0] / shape[:, 0].sum()


def test_profile_fitting_keeps_weight_on_isolated_cosmic_ray_pixels_in_a_centred_window() -> None:
    crhPixels = ((6, 10), (5, 30), (7, 31), (6, 55), (8, 70))
    images, truth = _centred_window_images(crhPixels=crhPixels)

    fitted = _fit_profile(images)

    profile = fitted["objectProfile"]
    for row, column in crhPixels:
        assert profile[row, column] > 0
        assert profile[row, column] == pytest.approx(truth[row], rel=5e-2)
    np.testing.assert_allclose(profile.sum(axis=0), 1.0)
    np.testing.assert_allclose(profile, np.tile(truth[:, np.newaxis], (1, 80)), rtol=5e-2, atol=1e-3)


def _slit_edge_images(
    columns: int = 80, rows: int = 13, offOrderRows: int = 4, sparseRows: tuple[int, ...] = (4, 5, 6)
) -> dict[str, np.ndarray]:
    """Build a window whose object sits on the slit edge.

    The top rows are off the order (NaN). The next rows are finite but masked, apart from
    one usable pixel in ten (none at either end), as order-edge flat flags leave them.
    """
    rowIndex = np.arange(rows)[:, np.newaxis]
    shape = np.exp(-0.5 * ((rowIndex - 5.0) / 1.5) ** 2)
    flux = shape * (1000.0 + 5.0 * np.arange(columns))[np.newaxis, :]
    mask = np.zeros((rows, columns), dtype=bool)
    mask[:offOrderRows, :] = True
    flux[:offOrderRows, :] = np.nan
    for row in sparseRows:
        mask[row, :] = np.arange(columns) % 10 != 5
    return {"fluxRaw": flux, "mask": mask}


def test_profile_fitting_gives_no_weight_to_masked_pixels_in_rows_with_sparse_support() -> None:
    images = _slit_edge_images()

    fitted = _fit_profile(images)

    profile = fitted["objectProfile"]
    assert (profile[:4, :] == 0).all()
    assert (profile[images["mask"]] == 0).all()
    assert (profile[~images["mask"]] > 0).all()
    np.testing.assert_allclose(profile.sum(axis=0), 1.0)


def test_unsupported_pixel_mask_flags_off_order_and_masked_pixels_in_sparse_rows_only() -> None:
    images = _slit_edge_images()
    images["mask"][9, 20] = True

    unsupported = _unsupported_pixel_mask(images["fluxRaw"], images["mask"])

    assert unsupported[:4, :].all()
    assert unsupported[4:7, :][images["mask"][4:7, :]].all()
    assert not unsupported[9, 20]
    assert not unsupported[~images["mask"]].any()


def test_slit_edge_truncated_columns_flags_an_object_cut_off_by_the_slit_edge() -> None:
    images = _slit_edge_images()
    unsupported = _unsupported_pixel_mask(images["fluxRaw"], images["mask"])

    truncated = _slit_edge_truncated_columns(images["fluxRaw"], images["mask"], unsupported)

    assert truncated.shape == (80,)
    assert truncated.all()


def test_slit_edge_truncated_columns_ignores_a_centred_object_with_off_order_edge_rows() -> None:
    images, _ = _centred_window_images()
    images["fluxRaw"][:2, :] = np.nan
    images["mask"][:2, :] = True
    unsupported = _unsupported_pixel_mask(images["fluxRaw"], images["mask"])

    truncated = _slit_edge_truncated_columns(images["fluxRaw"], images["mask"], unsupported)

    assert unsupported[:2, :].all()
    assert not truncated.any()


def test_slit_edge_truncated_columns_ignores_isolated_cosmic_rays() -> None:
    images, _ = _centred_window_images(crhPixels=((6, 10), (2, 30), (11, 31)))
    unsupported = _unsupported_pixel_mask(images["fluxRaw"], images["mask"])

    truncated = _slit_edge_truncated_columns(images["fluxRaw"], images["mask"], unsupported)

    assert not unsupported.any()
    assert not truncated.any()


def test_slit_edge_truncated_columns_flags_the_bottom_edge_and_skips_fully_unsupported_columns() -> None:
    images = _slit_edge_images()
    images["fluxRaw"] = images["fluxRaw"][::-1].copy()
    images["mask"] = images["mask"][::-1].copy()
    images["fluxRaw"][:, 7] = np.nan
    images["mask"][:, 7] = True
    unsupported = _unsupported_pixel_mask(images["fluxRaw"], images["mask"])

    truncated = _slit_edge_truncated_columns(images["fluxRaw"], images["mask"], unsupported)

    assert not truncated[7]
    assert np.delete(truncated, 7).all()


@pytest.mark.parametrize("flipRows", [False, True], ids=["top-edge", "bottom-edge"])
def test_slit_edge_truncated_columns_looks_past_a_cosmic_ray_on_the_row_next_to_the_slit_edge(flipRows: bool) -> None:
    images = _slit_edge_images(sparseRows=())
    images["mask"][4, 30] = True
    if flipRows:
        images = {key: value[::-1].copy() for key, value in images.items()}
    unsupported = _unsupported_pixel_mask(images["fluxRaw"], images["mask"])

    truncated = _slit_edge_truncated_columns(images["fluxRaw"], images["mask"], unsupported)

    assert not unsupported[4 if not flipRows else 8, 30]
    assert truncated[30]
    assert truncated.all()


def test_profile_fitting_is_unchanged_when_every_pixel_is_on_order() -> None:
    rowFractions = np.array([[0.25], [0.5], [0.25]])
    images = {
        "fluxRaw": rowFractions * (10.0 + np.arange(6.0))[np.newaxis, :],
        "mask": np.zeros((3, 6), dtype=bool),
    }

    fitted = _fit_profile(images)

    np.testing.assert_allclose(fitted["objectProfile"], np.tile(rowFractions, (1, 6)))


def test_single_order_extraction_returns_sorted_science_columns(log: object) -> None:
    slices = pd.DataFrame({"order": [10] * 5})
    images = {
        "bpMask": np.zeros((2, 5), dtype=bool),
        "fluxRaw": np.array([[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0]]),
        "variance": np.ones((2, 5)),
        "wavelength": np.array([[504.0, 503.0, 502.0, 501.0, 500.0], [504.0, 503.0, 502.0, 501.0, 500.0]]),
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
        "slitEdgeTruncated",
    ]
    assert result["slitEdgeTruncated"].tolist() == [False, False, False]
    assert result["wavelengthMean"].tolist() == [501.0, 502.0, 503.0]
    assert result["extractedFluxOptimal"].tolist() == [8.0, 6.0, 4.0]
    np.testing.assert_allclose(result["varianceSpectrum"], 2.0)
    assert result["skyFlux"].isna().all()


def test_single_order_extraction_marks_slices_where_the_object_spills_over_the_slit_edge(
    log: object,
) -> None:
    images = _slit_edge_images(columns=80)
    bpMask = images["mask"].astype(int)
    bpMask[:4, :] = 0
    images["fluxRaw"][:4, :] = np.nan
    images.update(
        {
            "bpMask": bpMask,
            "variance": np.ones((13, 80)),
            "wavelength": np.tile(600.0 - 0.1 * np.arange(80), (13, 1)),
        }
    )
    slices = pd.DataFrame({"order": [10] * 80})

    result = extract_single_order(
        (slices, images),
        log,
        slitHalfLength=6,
        clippingSigma=3.0,
        clippingIterationLimit=5,
        globalClippingSigma=3.0,
        axisA="y",
        axisB="x",
        hornePolyOrder=0,
    )

    assert result is not None
    assert result["slitEdgeTruncated"].dtype == bool
    assert result["slitEdgeTruncated"].mean() > 0.9


def _order_extraction(order: int, truncatedFraction: float, columns: int = 10, wlStart: float = 500.0) -> pd.DataFrame:
    """Build an order's extraction frame with the given share of slit-edge-truncated slices."""
    truncatedCount = round(truncatedFraction * columns)
    return pd.DataFrame(
        {
            "order": [order] * columns,
            "wavelengthMean": wlStart + np.arange(columns, dtype=float),
            "slitEdgeTruncated": [True] * truncatedCount + [False] * (columns - truncatedCount),
        }
    )


def _spill_extractor(log: object, notFlattened: str = "", qc: object = None) -> horne_extraction:
    extractor = _extractor(log)
    extractor.notFlattened = notFlattened
    extractor.qc = pd.DataFrame() if qc is None else qc
    extractor.recipeName = "soxs-stare"
    extractor.dateObs = "2024-01-02T03:04:05"
    extractor.arm = "NIR"
    return extractor


def test_slit_edge_summary_records_spilled_orders_logs_one_error_and_adds_the_qc_row(log: object) -> None:
    extractor = _spill_extractor(log)
    extractions = [
        _order_extraction(11, 0.9, wlStart=538.04),
        _order_extraction(12, 0.1, wlStart=600.0),
        _order_extraction(13, 1.0, wlStart=671.0),
    ]

    cleaned = extractor._summarise_slit_edge_spill(extractions)

    assert [o[0] for o in extractor.slitEdgeOrders] == [11, 13]
    assert extractor.slitEdgeOrders[0] == (11, 538.04, 547.04, pytest.approx(0.9))
    assert extractor.slitEdgeOrders[1][3] == pytest.approx(1.0)
    assert all("slitEdgeTruncated" not in e.columns for e in cleaned)
    errors = [message for level, message in log.messages if level == "error"]
    assert len(errors) == 1
    assert "slit edge" in errors[0]
    assert "538.0-547.0 nm" in errors[0]
    assert "671.0-680.0 nm" in errors[0]
    assert "90%" in errors[0]
    assert "unreliable" in errors[0]
    row = extractor.qc.iloc[0]
    assert (row["qc_name"], row["qc_value"]) == ("N ORDERS SLIT EDGE", 2)
    assert row["qc_comment"] == "Number of orders where the object spills over the slit edge"


def test_slit_edge_summary_names_soxs_vis_orders_as_the_products_number_them(log: object) -> None:
    extractor = _spill_extractor(log)
    extractor.arm = "VIS"
    extractions = [_order_extraction(1, 1.0, wlStart=538.0), _order_extraction(4, 1.0, wlStart=671.0)]

    extractor._summarise_slit_edge_spill(extractions)

    assert [o[0] for o in extractor.slitEdgeOrders] == [2, 1]
    errors = [message for level, message in log.messages if level == "error"]
    assert "order(s) 2 (538.0-547.0 nm" in errors[0]
    assert "1 (671.0-680.0 nm" in errors[0]


def test_slit_edge_summary_reports_a_clean_extraction_with_a_zero_count_and_no_error(log: object) -> None:
    extractor = _spill_extractor(log)
    extractions = [_order_extraction(11, 0.2), _order_extraction(12, 0.0)]

    extractor._summarise_slit_edge_spill(extractions)

    assert extractor.slitEdgeOrders == []
    assert [level for level, _ in log.messages if level == "error"] == []
    assert extractor.qc["qc_value"].tolist() == [0]


def test_slit_edge_summary_stays_silent_for_the_unflattened_re_extraction(log: object) -> None:
    extractor = _spill_extractor(log, notFlattened="_NOTFLAT")

    extractor._summarise_slit_edge_spill([_order_extraction(11, 1.0)])

    assert [o[0] for o in extractor.slitEdgeOrders] == [11]
    assert [level for level, _ in log.messages if level == "error"] == []
    assert extractor.qc.empty


def test_slit_edge_summary_leaves_a_missing_qc_table_alone(log: object) -> None:
    extractor = _spill_extractor(log, qc=False)

    extractor._summarise_slit_edge_spill([_order_extraction(11, 1.0)])

    assert extractor.qc is False
    assert len([level for level, _ in log.messages if level == "error"]) == 1


def test_extract_reports_slit_edge_spill_from_the_per_order_extractions(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    extractor, extraction, _ = _configure_synthetic_orchestration(log, monkeypatch, products=False)
    extraction["slitEdgeTruncated"] = True

    qc, _, _, _, _ = extractor.extract()

    assert extractor.slitEdgeOrders == [(10, 500.0, 502.0, 1.0)]
    assert qc["qc_value"].tolist() == [1]
    assert len([level for level, _ in log.messages if level == "error"]) == 1
