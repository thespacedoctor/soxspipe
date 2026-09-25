"""Public dispersion-map behavior with synthetic projection collaborators."""

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
from astropy.table import Table

from soxspipe.commonutils.create_dispersion_map import (
    create_dispersion_map,
    measure_line_position,
)

pytestmark = pytest.mark.unit


def test_create_new_static_line_list_expands_ranked_atlas_lines_to_detector_positions(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A first-pass map produces a stable detector-positioned line list."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.settings = {"instrument": "soxs"}
    mapper.detectorParams = {"line-atlas": "atlas.fits"}
    mapper.firstGuessMap = False
    calibrationRoot = tmp_path / "calibrations"
    calibrationRoot.mkdir()
    Table(
        {
            "ion": [b"ArI", b"NeI", b"ArI"],
            "source": [b"NIST", b"NIST", b"NIST"],
            "wave": [500.123456, 501.0, 650.0],
            "amplitude": [20.0, 10.0, 100.0],
        }
    ).write(calibrationRoot / "atlas.fits")
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    toolkitModule = importlib.import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(module, "get_calibrations_path", lambda **_: str(calibrationRoot))
    monkeypatch.setattr(
        toolkitModule,
        "read_spectral_format",
        lambda **_: (np.array([10]), np.array([499.0]), np.array([502.0])),
    )

    def project_lines(*, orderPixelTable: pd.DataFrame, **_: object) -> pd.DataFrame:
        return orderPixelTable.assign(
            fit_x=orderPixelTable["wavelength"] - 490.0,
            fit_y=orderPixelTable["order"] + orderPixelTable["slit_index"],
        )

    monkeypatch.setattr(commonutilsModule, "dispersion_map_to_pixel_arrays", project_lines)

    result = mapper.create_new_static_line_list("first-pass-map.fits")

    assert result["ion"].tolist() == ["ArI", "NeI"]
    assert result["slit_index"].tolist() == [4, 4]
    assert result["slit_position"].tolist() == [0.0, 0.0]
    np.testing.assert_allclose(result["detector_x"], [10.12346, 11.0])
    np.testing.assert_allclose(result["detector_y"], [14.0, 14.0])
    assert (tmp_path / "Hg.fits").exists()


def test_map_to_image_combines_orders_and_writes_three_extension_fits_product(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Serial map conversion sums order images into the established FITS layout."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.settings = {"instrument": "soxs"}
    mapper.kw = lambda key: key
    mapper.detectorParams = {}
    mapper.recipeSettings = {"map_to_image_displacement_threshold": 0.3}
    mapper.debug = False
    mapper.turnOffMP = True
    mapper.productDir = str(tmp_path)
    mapper.dispMapHeader = fits.Header({"ORIGIN": "synthetic"})
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    toolkitModule = importlib.import_module("soxspipe.commonutils.toolkit")
    fundamentalsModule = importlib.import_module("fundamentals")
    monkeypatch.setattr(
        module,
        "read_spectral_format",
        lambda **_: (np.array([10, 11]), np.array([500.0, 600.0]), np.array([501.0, 601.0])),
    )
    monkeypatch.setattr(toolkitModule, "frame_to_32", lambda frame: frame)

    def ccd(value: float) -> CCDData:
        return CCDData(np.full((2, 2), value), unit=u.electron)

    mapper.create_placeholder_images = lambda reverse: (ccd(0.0), ccd(0.0), ccd(7.0))
    received: dict[str, object] = {}

    def run_orders(**kwargs: object) -> list[tuple[CCDData, CCDData]]:
        received.update(kwargs)
        return [(ccd(1.0), ccd(10.0)), (ccd(2.0), ccd(20.0))]

    monkeypatch.setattr(fundamentalsModule, "fmultiprocess", run_orders)

    result = mapper.map_to_image("synthetic-map.fits", orders=[11, 10])

    assert received["turnOffMP"] is True
    assert received["inputArray"] == [(10, 500.0, 501.0), (11, 600.0, 601.0)]
    with fits.open(result) as hdus:
        assert [hdu.name for hdu in hdus] == ["WAVELENGTH", "SLIT", "ORDER"]
        np.testing.assert_allclose(hdus[0].data, 30.0)
        np.testing.assert_allclose(hdus[1].data, 3.0)
        np.testing.assert_allclose(hdus[2].data, 7.0)
        assert hdus[0].header["PRO_CATG"] == "DISP_IMAGE_VIS"


def test_convert_and_fit_places_nearest_valid_projection_in_detector_maps(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Projected grid points fill valid pixels and retain unmatched candidates."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.dispersionMapPath = "synthetic.fits"
    mapper.map_to_image_displacement_threshold = 0.3
    slitMap = SimpleNamespace(data=np.full((3, 3), np.nan, dtype=float))
    wavelengthMap = SimpleNamespace(data=np.full((3, 3), np.nan, dtype=float))
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")

    def project(*, orderPixelTable: pd.DataFrame, **_: object) -> pd.DataFrame:
        return orderPixelTable.assign(fit_x=[1.1, 3.6], fit_y=[1.1, 0.2])

    monkeypatch.setattr(module, "dispersion_map_to_pixel_arrays", project)

    remaining, remainingCount = mapper.convert_and_fit(
        order=10,
        bigWlArray=np.array([500.0, 501.0]),
        bigSlitArray=np.array([0.0, 1.0]),
        slitMap=slitMap,
        wlMap=wavelengthMap,
        iteration=1,
    )

    assert remainingCount == 1
    assert remaining["wavelength"].tolist() == [501.0]
    assert wavelengthMap.data[1, 1] == pytest.approx(500.0)
    assert slitMap.data[1, 1] == pytest.approx(0.0)


def test_order_to_image_stops_after_a_terminal_projection_pass(log: object) -> None:
    """One terminal projection pass returns the mapper's two detector images."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.detectorParams = {"slit_length": 2.0}
    mapper.create_placeholder_images = lambda order: (
        SimpleNamespace(data=np.full((2, 2), np.nan)),
        SimpleNamespace(data=np.ones((2, 2))),
        SimpleNamespace(data=np.zeros((2, 2))),
    )
    calls: list[dict[str, object]] = []

    def convert_and_fit(**kwargs: object) -> tuple[pd.DataFrame, int]:
        calls.append(kwargs)
        return pd.DataFrame(), 0

    mapper.convert_and_fit = convert_and_fit

    slitMap, wavelengthMap = mapper.order_to_image((10, 500.0, 502.0))

    assert len(calls) == 1
    assert calls[0]["order"] == 10
    assert calls[0]["iteration"] == 1
    assert slitMap.data.shape == (2, 2)
    assert wavelengthMap.data.shape == (2, 2)


def test_update_static_line_list_projects_each_unique_line_at_each_slit_position(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Expand each unique predicted line through every configured slit position."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.uniqueSlitPos = [-0.5, 0.5]
    inputLines = pd.DataFrame(
        {
            "order": [10, 10, 11],
            "wavelength": [500.0, 500.0, 600.0],
            "ignored": ["first", "duplicate", "second"],
        }
    )
    commonutilsModule = importlib.import_module("soxspipe.commonutils")

    def project_lines(
        *, orderPixelTable: pd.DataFrame, **_: object
    ) -> pd.DataFrame:
        return orderPixelTable.assign(
            fit_x=orderPixelTable["order"] * 10 + orderPixelTable["slit_index"],
            fit_y=orderPixelTable["wavelength"] + orderPixelTable["slit_position"],
        )

    monkeypatch.setattr(commonutilsModule, "dispersion_map_to_pixel_arrays", project_lines)

    result = mapper.update_static_line_list_detector_positions(
        inputLines,
        "synthetic-dispersion-map.fits",
    )

    assert result.to_dict("records") == [
        {
            "wavelength": 500.0,
            "order": 10,
            "slit_index": 0,
            "slit_position": -0.5,
            "detector_x": 100,
            "detector_y": 499.5,
        },
        {
            "wavelength": 600.0,
            "order": 11,
            "slit_index": 0,
            "slit_position": -0.5,
            "detector_x": 110,
            "detector_y": 599.5,
        },
        {
            "wavelength": 500.0,
            "order": 10,
            "slit_index": 1,
            "slit_position": 0.5,
            "detector_x": 101,
            "detector_y": 500.5,
        },
        {
            "wavelength": 600.0,
            "order": 11,
            "slit_index": 1,
            "slit_position": 0.5,
            "detector_x": 111,
            "detector_y": 600.5,
        },
    ]


@pytest.mark.parametrize(
    ("brightest", "expectedX"),
    [(False, 1.0), (True, 3.0)],
)
def test_measure_line_position_selects_nearest_or_brightest_source_and_handles_iraf_failure(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    brightest: bool,
    expectedX: float,
) -> None:
    """Single-source mode chooses by policy and degrades cleanly without IRAF FWHM."""
    photutilsModule = importlib.import_module("photutils")
    sources = Table(
        {
            "xcentroid": [1.0, 3.0],
            "ycentroid": [1.0, 3.0],
            "sharpness": [0.1, 0.2],
            "roundness1": [0.1, 0.2],
            "roundness2": [0.1, 0.2],
            "npix": [4, 5],
            "sky": [0.0, 0.0],
            "peak": [10.0, 100.0],
            "flux": [10.0, 100.0],
        }
    )
    monkeypatch.setattr(photutilsModule, "DAOStarFinder", lambda **_: lambda *args, **kwargs: sources)
    monkeypatch.setattr(
        photutilsModule,
        "IRAFStarFinder",
        lambda **_: (_ for _ in ()).throw(RuntimeError("unavailable")),
    )
    stamp = CCDData(np.arange(36.0).reshape(6, 6), unit=u.adu)

    result = measure_line_position(
        (stamp, 4, 10, 8, 14),
        log,
        windowHalf=2,
        iraf=True,
        sigmaLimit=3,
        iteration=1,
        brightest=brightest,
        returnAll=False,
    )

    assert len(result) == 1
    assert result[0]["observed_x"] == pytest.approx(expectedX + 4.0)
    assert np.isnan(result[0]["fwhm_pin_px"])


def test_write_map_to_file_serializes_per_axis_coefficients_and_metadata(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Serialize asymmetric axis degrees without writing a real FITS product."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.kw = lambda key: key
    mapper.sofName = "SYNTHETIC"
    mapper.pinholeFrame = type(
        "Frame",
        (),
        {"header": {"DPR_TECH": "ARC", "DPR_CATG": "CALIB", "DPR_TYPE": "LAMP"}},
    )()
    mapper.productDir = "/products"
    mapper.settings = {}
    mapper.qc = pd.DataFrame(
        columns=["qc_name", "qc_value", "qc_unit", "qc_comment"]
    )
    phase3Module = importlib.import_module("soxspipe.commonutils.phase3")
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        phase3Module,
        "write_fits_table_to_disk",
        lambda **kwargs: captured.update(kwargs),
    )

    result = mapper.write_map_to_file(
        xcoeff=[1.0, 2.0],
        ycoeff=[3.0, 4.0, 5.0, 6.0],
        orderDeg=[1, 0],
        wavelengthDeg=[0, 1],
        slitDeg=[0, 1],
    )

    table = captured["tables"][0].to_pandas()
    assert result == "/products/SYNTHETIC.fits"
    assert table["axis"].tolist() == ["x", "y"]
    assert table["order_deg"].tolist() == [1, 0]
    assert table["wavelength_deg"].tolist() == [0, 1]
    assert table["slit_deg"].tolist() == [0, 1]
    assert table.loc[0, "c000"] == pytest.approx(1.0)
    assert table.loc[0, "c100"] == pytest.approx(2.0)
    assert table.loc[1, "c000"] == pytest.approx(3.0)
    assert table.loc[1, "c001"] == pytest.approx(4.0)
    assert table.loc[1, "c010"] == pytest.approx(5.0)
    assert table.loc[1, "c011"] == pytest.approx(6.0)
    assert captured["header"]["PRO_TECH"] == "ECHELLE,MULTI-PINHOLE"
    assert captured["header"]["SEQ_ARM"] == "VIS"
    assert captured["header"]["PRO_TYPE"] == "REDUCED"
    assert captured["header"]["PRO_CATG"] == "DISP_TAB_VIS"


def _fitting_mapper(log: object) -> create_dispersion_map:
    """Return a mapper configured for one linear synthetic fitting pass."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.recipeName = "soxs-disp-solution"
    mapper.settings = {}
    mapper.detectorParams = {}
    mapper.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": 4.0,
        "poly-clipping-iteration-limit": 2,
    }
    return mapper


def _synthetic_fitted_lines() -> pd.DataFrame:
    """Return two retained lines and one explicitly dropped line."""
    return pd.DataFrame(
        {
            "order": [10.0, 10.0, 11.0],
            "wavelength": [500.0, 501.0, 600.0],
            "slit_position": [0.0, 0.0, 0.0],
            "observed_x": [100.0, 101.0, 200.0],
            "observed_y": [20.0, 21.0, 30.0],
            "dropped": [False, False, True],
        }
    )


def test_fit_polynomials_returns_axis_error_when_first_fit_fails(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Report an x-axis fitting failure before residual processing begins."""
    mapper = _fitting_mapper(log)
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    optimiseModule = importlib.import_module("scipy.optimize")
    monkeypatch.setattr(commonutilsModule, "get_cached_coeffs", lambda **_: ([0.0], [0.0]))
    monkeypatch.setattr(optimiseModule, "curve_fit", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("fit failed")))

    result = mapper.fit_polynomials(
        _synthetic_fitted_lines(),
        wavelengthDeg=0,
        orderDeg=0,
        slitDeg=0,
    )

    assert result == ("xerror", None, None, None)


def test_fit_polynomials_returns_y_axis_error_after_a_successful_x_fit(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Preserve the distinct y-fit failure result used by degree fallback."""
    mapper = _fitting_mapper(log)
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    optimiseModule = importlib.import_module("scipy.optimize")
    monkeypatch.setattr(commonutilsModule, "get_cached_coeffs", lambda **_: ([0.0], [0.0]))
    fitResults = iter([([1.0], None), RuntimeError("y fit failed")])

    def fit(*args: object, **kwargs: object) -> tuple[list[float], None]:
        result = next(fitResults)
        if isinstance(result, Exception):
            raise result
        return result

    monkeypatch.setattr(optimiseModule, "curve_fit", fit)

    result = mapper.fit_polynomials(
        _synthetic_fitted_lines(),
        wavelengthDeg=0,
        orderDeg=0,
        slitDeg=0,
    )

    assert result == (None, "yerror", None, None)


def test_fit_polynomials_retains_dropped_lines_and_runs_final_qc_fit(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Fit retained lines, retain dropped rows, and request final QC residuals."""
    mapper = _fitting_mapper(log)
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    optimiseModule = importlib.import_module("scipy.optimize")
    monkeypatch.setattr(commonutilsModule, "get_cached_coeffs", lambda **_: ([0.0], [0.0]))
    fitResults = iter([([1.5], None), ([2.5], None)])
    monkeypatch.setattr(optimiseModule, "curve_fit", lambda *args, **kwargs: next(fitResults))
    qcCalls: list[bool] = []

    def calculate_residuals(**kwargs: object) -> tuple[float, float, float, pd.DataFrame]:
        qcCalls.append(bool(kwargs["writeQCs"]))
        table = kwargs["orderPixelTable"].copy()
        return 0.0, 0.0, 0.0, table.assign(
            residuals_x=0.0,
            residuals_y=0.0,
            residuals_xy=0.0,
        )

    mapper.calculate_residuals = calculate_residuals

    xcoeff, ycoeff, retainedLines, clippedLines = mapper.fit_polynomials(
        _synthetic_fitted_lines(),
        wavelengthDeg=0,
        orderDeg=0,
        slitDeg=0,
    )

    assert xcoeff == [1.5]
    assert ycoeff == [2.5]
    assert retainedLines["wavelength"].tolist() == [500.0, 501.0]
    assert clippedLines["wavelength"].tolist() == [600.0]
    assert qcCalls == [False, True]


def test_detect_pinhole_arc_lines_uses_shifted_positions_and_collects_measurements(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Build bounded stamps and retain each candidate measurement set per line."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.pinholeFrameMasked = np.arange(48.0).reshape(6, 8)
    mapper.windowHalf = 2
    mapper.debug = False
    mapper.turnOffMP = False
    inputLines = pd.DataFrame(
        {
            "detector_x": [7.0, 0.0],
            "detector_y": [5.0, 0.0],
            "detector_x_shifted": [1.1, 5.8],
            "detector_y_shifted": [0.5, 4.4],
        }
    )
    fundamentalsModule = importlib.import_module("fundamentals")
    captured: dict[str, object] = {}
    metricNames = [
        "sharpness",
        "roundness1",
        "roundness2",
        "npix",
        "sky",
        "peak",
        "flux",
        "observed_x",
        "observed_y",
        "fwhm_pin_px",
    ]
    measurements = [
        [dict.fromkeys(metricNames, 10.0)],
        [dict.fromkeys(metricNames, 20.0)],
    ]

    def measure_stamps(**kwargs: object) -> list[list[dict[str, float]]]:
        captured.update(kwargs)
        return measurements

    monkeypatch.setattr(fundamentalsModule, "fmultiprocess", measure_stamps)

    result = mapper.detect_pinhole_arc_lines(
        inputLines,
        iraf=False,
        sigmaLimit=5,
        iteration=2,
        brightest=True,
        exclude_border=True,
        returnAll=True,
    )

    stamps = captured["inputArray"]
    assert [(stamp[1], stamp[2], stamp[3], stamp[4]) for stamp in stamps] == [
        (0, 3, 0, 2),
        (4, 8, 2, 6),
    ]
    np.testing.assert_array_equal(stamps[0][0], np.arange(48.0).reshape(6, 8)[0:2, 0:3])
    assert captured["turnOffMP"] is True
    assert captured["iraf"] is False
    assert captured["sigmaLimit"] == 5
    assert captured["iteration"] == 2
    assert captured["brightest"] is True
    assert captured["exclude_border"] is True
    assert captured["returnAll"] is True
    assert result["observed_x"].tolist() == [[10.0], [20.0]]
    assert result["fwhm_pin_px"].tolist() == [[10.0], [20.0]]


def test_get_writes_the_fitted_map_and_returns_workflow_products(log: object) -> None:
    """Run the public map workflow with deterministic detection and output seams."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.lineDetectionTable = pd.DataFrame(
        {
            "wavelength": [500.0, 501.0],
            "order": [10, 10],
            "slit_position": [0.0, 0.0],
            "observed_x": [100.0, 101.0],
        }
    )
    mapper.recipeName = "soxs-disp-solution"
    mapper.recipeSettings = {"mph_line_set_min": 1}
    mapper.settings = {"tune-pipeline": False}
    mapper.firstGuessMap = False
    mapper.orderTable = False
    mapper.create2DMap = False
    mapper.qc = pd.DataFrame()
    mapper.products = pd.DataFrame()
    mapper.dateObs = "2025-01-02T03:04:05"

    predictedLines = mapper.lineDetectionTable.copy()
    clippedLines = predictedLines.iloc[0:0].copy()
    calls: list[str] = []
    mapper._initialize_recipe_settings = lambda: (False, False, 1, 1, 0)
    mapper.get_predicted_line_list = lambda: predictedLines.copy()
    mapper._prepare_pinhole_frame = lambda: "prepared-frame"
    mapper._generate_quicklook_and_debug_plots = lambda *args: calls.append("quicklook")
    mapper._clip_on_measured_line_metrics = lambda table: table
    mapper.fit_polynomials = lambda **kwargs: (
        np.array([1.0]),
        np.array([2.0]),
        predictedLines.copy(),
        clippedLines.copy(),
    )
    mapper._get_output_filenames = lambda: ("good.fits", "missing.fits")
    mapper._write_qc_metrics = lambda *args: calls.append("metrics")
    mapper._write_line_list_qc = lambda *args: calls.append("line-qc")
    mapper._prepare_line_list_columns = lambda *args: (predictedLines.copy(), predictedLines.copy())
    mapper._write_fitted_lines_file = lambda *args: calls.append("fitted")
    mapper._write_missing_lines_file = lambda *args: calls.append("missing")
    mapper.write_map_to_file = lambda *args: "/tmp/synthetic-dispersion-map.fits"
    mapper._create_dispersion_map_qc_plot = lambda **kwargs: ["residuals.pdf"]

    result = mapper.get()

    assert result[0:3] == ("/tmp/synthetic-dispersion-map.fits", None, ["residuals.pdf"])
    assert result[3] is mapper.qc
    assert result[4] is mapper.products
    assert result[5].equals(mapper.lineDetectionTable)
    assert calls == ["quicklook", "metrics", "line-qc", "fitted", "missing"]


def test_get_detects_and_shifts_lines_before_fitting(log: object) -> None:
    """Run the three-pass detector branch before the deterministic fitting seam."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.lineDetectionTable = False
    mapper.recipeName = "soxs-disp-solution"
    mapper.recipeSettings = {"mph_line_set_min": 1}
    mapper.settings = {"tune-pipeline": False}
    mapper.firstGuessMap = False
    mapper.orderTable = False
    mapper.create2DMap = False
    mapper.qc = pd.DataFrame()
    mapper.products = pd.DataFrame()
    mapper.dateObs = "2025-01-02T03:04:05"
    mapper.arm = "VIS"
    mapper.inst = "SOXS"
    mapper.debug = False

    lines = pd.DataFrame(
        {
            "wavelength": [500.0, 501.0],
            "order": [10, 10],
            "slit_position": [0.0, 0.0],
            "detector_x": [100.0, 101.0],
            "detector_y": [20.0, 21.0],
            "observed_x": [100.2, 101.2],
            "observed_y": [20.1, 21.1],
        }
    )
    detectedIterations: list[tuple[int, bool]] = []
    mapper._initialize_recipe_settings = lambda: (False, False, 1, 1, 0)
    mapper.get_predicted_line_list = lambda: lines.copy()
    mapper._prepare_pinhole_frame = lambda: "prepared-frame"
    mapper._generate_quicklook_and_debug_plots = lambda *args: None
    mapper._get_detection_parameters = lambda iteration, tight: (3, 5, False, False, False)

    def detect(*args: object, **kwargs: object) -> pd.DataFrame:
        detectedIterations.append((int(kwargs["iteration"]), bool(kwargs["iraf"])))
        return lines.copy()

    mapper.detect_pinhole_arc_lines = detect
    mapper._explode_multiple_detections = lambda table: table
    mapper._calculate_position_differences = lambda table: table
    mapper._find_and_apply_cluster_shift = lambda table, mask: (table, 0.2, 0.1)
    mapper._calculate_order_shift_statistics = lambda table, mask: (0.2, 0.1, 0.01, 0.01, 0.3, 0.02)
    mapper._clip_on_measured_line_metrics = lambda table: table
    mapper.fit_polynomials = lambda **kwargs: (
        np.array([1.0]),
        np.array([2.0]),
        lines.copy(),
        lines.iloc[0:0].copy(),
    )
    mapper._get_output_filenames = lambda: ("good.fits", "missing.fits")
    mapper._write_qc_metrics = lambda *args: None
    mapper._write_line_list_qc = lambda *args: None
    mapper._prepare_line_list_columns = lambda *args: (lines.copy(), lines.copy())
    mapper._write_fitted_lines_file = lambda *args: None
    mapper._write_missing_lines_file = lambda *args: None
    mapper.write_map_to_file = lambda *args: "/tmp/synthetic-dispersion-map.fits"
    mapper._create_dispersion_map_qc_plot = lambda **kwargs: []

    result = mapper.get()

    assert detectedIterations == [(0, False), (1, True), (2, True)]
    assert result[0] == "/tmp/synthetic-dispersion-map.fits"
    assert result[5]["shift_group"].tolist() == [10, 10]
    assert result[5]["detector_x_shifted"].tolist() == [100.0, 101.0]


def test_get_bootstraps_a_spatial_map_then_writes_its_image_product(log: object) -> None:
    """Bootstrap once, retain complete pinhole sets, and return the map image."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    lines = pd.DataFrame(
        {
            "wavelength": [500.0] * 9,
            "order": [13] * 9,
            "slit_position": np.linspace(-1.0, 1.0, 9),
            "observed_x": np.arange(100.0, 109.0),
        }
    )
    mapper.lineDetectionTable = lines.copy()
    mapper.recipeName = "soxs-spat-solution"
    mapper.recipeSettings = {"mph_line_set_min": 1}
    mapper.settings = {"tune-pipeline": False}
    mapper.firstGuessMap = "first-guess.fits"
    mapper.orderTable = "orders.fits"
    mapper.create2DMap = True
    mapper.qc = pd.DataFrame()
    mapper.products = pd.DataFrame()
    mapper.dateObs = "2025-01-02T03:04:05"
    mapper.inst = "SOXS"

    mapWrites: list[str] = []
    mapper._initialize_recipe_settings = lambda: (True, False, 1, 1, 1)
    mapper.get_predicted_line_list = lambda: lines.copy()
    mapper._prepare_pinhole_frame = lambda: "prepared-frame"
    mapper._generate_quicklook_and_debug_plots = lambda *args: None
    mapper._clip_on_measured_line_metrics = lambda table: table
    mapper.fit_polynomials = lambda **kwargs: (
        np.array([1.0]),
        np.array([2.0]),
        lines.copy(),
        lines.iloc[0:0].copy(),
    )
    mapper.create_new_static_line_list = lambda **kwargs: lines.copy()
    mapper._get_output_filenames = lambda: ("good.fits", "missing.fits")
    mapper._write_qc_metrics = lambda *args: None
    mapper._write_line_list_qc = lambda *args: None
    mapper._prepare_line_list_columns = lambda *args: (lines.copy(), lines.copy())
    mapper._write_fitted_lines_file = lambda *args: None
    mapper._write_missing_lines_file = lambda *args: None

    def write_map(*args: object) -> str:
        mapWrites.append("map")
        return "/tmp/synthetic-dispersion-map.fits"

    mapper.write_map_to_file = write_map
    mapper.map_to_image = lambda **kwargs: "/tmp/synthetic-dispersion-map_IMAGE.fits"
    mapper._create_dispersion_map_qc_plot = lambda **kwargs: ["residuals.pdf"]

    result = mapper.get()

    assert mapWrites == ["map", "map"]
    assert result[0:3] == (
        "/tmp/synthetic-dispersion-map.fits",
        "/tmp/synthetic-dispersion-map_IMAGE.fits",
        ["residuals.pdf"],
    )
    assert mapper.minpin == 9
    assert mapper.qc["qc_name"].tolist() == ["PINHOLE COUNT MIN"]
    assert mapper.qc["qc_value"].tolist() == [9]


def test_map_to_image_writes_wavelength_slit_and_order_extensions(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: object,
) -> None:
    """Combine per-order synthetic maps into the public three-extension FITS product."""
    from astropy import units as u
    from astropy.io import fits
    from astropy.nddata import CCDData

    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.kw = lambda name: name
    mapper.detectorParams = {}
    mapper.arm = "VIS"
    mapper.recipeSettings = {"map_to_image_displacement_threshold": 0.25}
    mapper.settings = {}
    mapper.debug = False
    mapper.turnOffMP = True
    mapper.productDir = str(tmp_path)
    mapper.dispMapHeader = fits.Header({"ORIGIN": "synthetic"})
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    fundamentalsModule = importlib.import_module("fundamentals")
    monkeypatch.setattr(module, "read_spectral_format", lambda **_: ([10, 11], [500.0, 600.0], [501.0, 601.0]))
    monkeypatch.setattr(
        fundamentalsModule,
        "fmultiprocess",
        lambda **kwargs: [kwargs["function"](item) for item in kwargs["inputArray"]],
    )

    def placeholders(*, reverse: bool = False, order: int | bool = False) -> tuple[CCDData, CCDData, CCDData]:
        if reverse:
            return (
                CCDData(np.full((3, 3), 1000.0), unit=u.adu),
                CCDData(np.full((3, 3), 2000.0), unit=u.adu),
                CCDData(np.full((3, 3), 3000.0), unit=u.adu),
            )
        value = float(order)
        return tuple(CCDData(np.full((3, 3), value), unit=u.adu) for _ in range(3))

    mapper.create_placeholder_images = placeholders
    mapper.order_to_image = lambda info: (
        CCDData(np.full((3, 3), float(info[0])), unit=u.adu),
        CCDData(np.full((3, 3), float(info[1])), unit=u.adu),
    )

    result = mapper.map_to_image("synthetic-map.fits", orders=[11])

    assert result == str(tmp_path / "synthetic-map_IMAGE.fits")
    with fits.open(result) as hdul:
        assert [hdu.header["EXTNAME"] for hdu in hdul] == ["WAVELENGTH", "SLIT", "ORDER"]
        np.testing.assert_allclose(hdul[0].data, 2600.0)
        np.testing.assert_allclose(hdul[1].data, 1011.0)
        np.testing.assert_allclose(hdul[2].data, 3000.0)
        assert hdul[0].header["PRO_CATG"] == "DISP_IMAGE_VIS"


def test_convert_and_fit_keeps_nearest_model_value_per_detector_pixel(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Map analytic samples into the closest detector pixels and report remaining rows."""
    from astropy import units as u
    from astropy.nddata import CCDData

    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.dispersionMapPath = "synthetic-map.fits"
    mapper.map_to_image_displacement_threshold = 0.5
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")

    def project(**kwargs: object) -> pd.DataFrame:
        table = kwargs["orderPixelTable"].copy()
        return table.assign(fit_x=[1.1, 1.2, 8.6], fit_y=[1.1, 1.2, 8.6])

    monkeypatch.setattr(module, "dispersion_map_to_pixel_arrays", project)
    slitMap = CCDData(np.full((3, 3), np.nan), unit=u.adu)
    wavelengthMap = CCDData(np.full((3, 3), np.nan), unit=u.adu)

    remaining, remainingCount = mapper.convert_and_fit(
        order=10,
        bigWlArray=np.array([500.0, 501.0, 502.0]),
        bigSlitArray=np.array([0.0, 0.1, 0.2]),
        slitMap=slitMap,
        wlMap=wavelengthMap,
        iteration=1,
    )

    assert remainingCount == 1
    assert remaining[["pixel_x", "pixel_y"]].values.tolist() == [[9, 9]]
    assert wavelengthMap.data[1, 1] == pytest.approx(500.0)
    assert slitMap.data[1, 1] == pytest.approx(0.0)


@pytest.mark.parametrize(
    "reverse, expectedInside, expectedOutside, expectedOrder",
    [(False, np.nan, 0.0, np.nan), (True, 0.0, np.nan, 10.0)],
)
def test_create_placeholder_images_marks_only_order_footprints(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    reverse: bool,
    expectedInside: float,
    expectedOutside: float,
    expectedOrder: float,
) -> None:
    """Build detector-sized placeholder maps with stable order-footprint semantics."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.kw = lambda name: name
    mapper.orderTable = "synthetic-orders.fits"
    mapper.detectorParams = {
        "science-pixels": {"columns": {"start": 0, "end": 5}, "rows": {"start": 0, "end": 4}},
        "dispersion-axis": "x",
    }
    mapper.axisA = "x"
    mapper.axisB = "y"
    mapper.debug = False
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    pixelTable = pd.DataFrame(
        {
            "order": [10, 10, 11],
            "ycoord": [1, 2, 1],
            "xcoord_edgeup": [4.0, 3.0, 2.0],
            "xcoord_edgelow": [1.0, 0.0, 0.0],
        }
    )
    monkeypatch.setattr(module, "unpack_order_table", lambda **_: (pd.DataFrame(), pixelTable, pd.DataFrame()))

    slitMap, wavelengthMap, orderMap = mapper.create_placeholder_images(reverse=reverse, order=10)

    assert wavelengthMap.shape == (4, 5)
    assert slitMap.data[1, 1] == pytest.approx(expectedInside, nan_ok=True)
    assert wavelengthMap.data[0, 0] == pytest.approx(expectedOutside, nan_ok=True)
    assert orderMap.data[1, 1] == pytest.approx(expectedOrder, nan_ok=True)


def test_predicted_line_list_loads_cleans_and_selects_the_mid_slit(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: object,
) -> None:
    """Load the calibration table through the public prediction workflow."""
    from types import SimpleNamespace

    from astropy.table import Table

    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.inst = "SOXS"
    mapper.firstGuessMap = False
    mapper.settings = {"instrument": "soxs"}
    mapper.kw = lambda name: name
    mapper.pinholeFrame = SimpleNamespace(
        header={"DPR_TECH": "ECHELLE,PINHOLE", "WIN_BINX": 1, "WIN_BINY": 1}
    )
    mapper.detectorParams = {
        "mid_slit_index": 1,
        "science-pixels": {
            "columns": {"start": 0, "end": 8},
            "rows": {"start": 0, "end": 6},
        },
        "predicted pinhole lines": {"single": {"1x1": "predicted.fits"}},
    }
    calibrationPath = tmp_path / "calibrations"
    calibrationPath.mkdir()
    Table(
        {
            "ORDER": [10, 10, 10, 10],
            "WAVELENGTH": [500.0, 500.0, 501.0, 502.0],
            "slit_index": [1, 1, 0, 1],
            "slit_position": [0.0, 0.0, -1.0, 0.0],
            "detector_x": [10.0, 10.0, 11.0, 12.0],
            "detector_y": [20.0, 20.0, 21.0, 22.0],
            "ion": ["Ar", "Ar", "Ne", "Kr"],
            "delete": [0, 0, 0, 1],
        }
    ).write(calibrationPath / "predicted.fits")
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    monkeypatch.setattr(module, "get_calibrations_path", lambda **_: str(calibrationPath))

    result = mapper.get_predicted_line_list()

    assert result.to_dict("records") == [
        {
            "order": 10.0,
            "wavelength": 500.0,
            "slit_index": 1,
            "slit_position": 0.0,
            "detector_x": 10.0,
            "detector_y": 20.0,
            "ion": "Ar",
            "delete": 0,
        }
    ]


def test_first_guess_corrections_shift_complete_multi_pinhole_sets(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Apply fitted first-guess offsets while removing incomplete line groups."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.firstGuessMap = "first-guess.fits"
    mapper.detectorParams = {"mid_slit_index": 1}
    source = pd.DataFrame(
        {
            "wavelength": [500.0, 500.0, 501.0, 501.0, 600.0],
            "order": [10, 10, 10, 10, 11],
            "slit_index": [0, 1, 0, 1, 1],
            "detector_x": [10.0, 10.0, 20.0, 20.0, 30.0],
            "detector_y": [15.0, 15.0, 25.0, 25.0, 35.0],
        }
    )
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    monkeypatch.setattr(
        module,
        "dispersion_map_to_pixel_arrays",
        lambda **kwargs: kwargs["orderPixelTable"].assign(
            fit_x=kwargs["orderPixelTable"]["detector_x"] - 0.5,
            fit_y=kwargs["orderPixelTable"]["detector_y"] + 0.25,
        ),
    )

    result = mapper._apply_first_guess_corrections(source)

    assert set(result["wavelength"]) == {500.0, 501.0}
    assert result["detector_x_shifted"].tolist() == pytest.approx([9.5, 9.5, 19.5, 19.5])
    assert result["detector_y_shifted"].tolist() == pytest.approx([15.25, 15.25, 25.25, 25.25])
    assert "fit_x" not in result
    assert "shift_x" not in result


@pytest.mark.parametrize(("debug", "flip"), [(False, False), (True, False), (False, True)])
def test_qc_plot_writes_single_pinhole_residual_product(
    debug: bool,
    flip: bool,
    log: object,
    tmp_path: Path,
) -> None:
    """A single-pinhole solution creates the established residual QC product."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.kw = lambda key: key
    mapper.detectorParams = {"rotate-qc-plot": 0, "flip-qc-plot": flip}
    mapper.pinholeFrameMasked = CCDData(
        np.arange(72, dtype=float).reshape(6, 12),
        unit=u.adu,
        mask=np.zeros((6, 12), dtype=bool),
    )
    mapper.arcFrame = False
    mapper.firstGuessMap = False
    mapper.debug = debug
    mapper.axisA = "x"
    mapper.axisB = "y"
    mapper.meanFrameFlux = 20.0
    mapper.stdFrameFlux = 5.0
    mapper.sofName = "SYNTHETIC"
    mapper.qcDir = str(tmp_path)
    mapper.products = pd.DataFrame()
    mapper.recipeName = "soxs-disp-solution"
    mapper.dateObs = "2025-01-01T00:00:00"
    mapper.qc = pd.DataFrame(
        columns=["qc_name", "qc_value", "qc_unit", "qc_comment"]
    )
    mapper.recipeSettings = {"sample_setting": 1}
    mapper.exptime = 1.0
    mapper.settings = {"tune-pipeline": False}
    base = {
        "order": [10, 10, 11, 11],
        "wavelength": [500.0, 501.0, 600.0, 601.0],
        "observed_x": [1.0, 2.0, 3.0, 4.0],
        "observed_y": [1.0, 2.0, 2.0, 3.0],
        "detector_x": [1.1, 2.1, 3.1, 4.1],
        "detector_y": [1.1, 2.1, 2.1, 3.1],
        "detector_x_shifted": [1.0, 2.0, 3.0, 4.0],
        "detector_y_shifted": [1.0, 2.0, 2.0, 3.0],
        "fit_x": [1.05, 2.05, 3.05, 4.05],
        "fit_y": [1.05, 2.05, 2.05, 3.05],
        "residuals_x": [0.05, -0.05, 0.1, -0.1],
        "residuals_y": [0.05, -0.05, -0.1, 0.1],
        "fwhm_pin_px": [2.0, 2.1, 1.9, 2.0],
        "R_pin": [4000.0, 4100.0, 3900.0, 4050.0],
    }
    fitted = pd.DataFrame(base)
    clipped = pd.DataFrame(base).assign(dropped=False)
    missing = pd.DataFrame(base).iloc[:1].copy()

    result = mapper._create_dispersion_map_qc_plot(
        xcoeff=[1.0],
        ycoeff=[1.0],
        orderDeg=1,
        wavelengthDeg=1,
        slitDeg=0,
        orderPixelTable=fitted,
        missingLines=missing,
        allClippedLines=clipped,
    )

    assert result == "SYNTHETIC_RESIDUALS_110.pdf"
    assert (tmp_path / result).is_file()
    assert mapper.products.iloc[0]["product_label"] == "DISP_MAP_RES"


def test_fit_polynomials_clips_outliers_and_retains_previously_dropped_lines(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The iterative fit retains both input and residual-clipped line records."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.settings = {}
    mapper.recipeName = "soxs-disp-solution"
    mapper.detectorParams = {}
    mapper.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": 1.0,
        "poly-clipping-iteration-limit": 3,
    }
    source = pd.DataFrame(
        {
            "order": [10, 10, 10, 10, 10],
            "wavelength": [500.0, 501.0, 502.0, 503.0, 504.0],
            "slit_position": [0.0] * 5,
            "observed_x": [1.0, 2.0, 3.0, 4.0, 5.0],
            "observed_y": [5.0, 4.0, 3.0, 2.0, 1.0],
            "dropped": [False, False, False, False, True],
        }
    )
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    optimizeModule = importlib.import_module("scipy.optimize")
    monkeypatch.setattr(
        commonutilsModule,
        "get_cached_coeffs",
        lambda **_: (np.array([0.0]), np.array([0.0])),
    )
    monkeypatch.setattr(
        optimizeModule,
        "curve_fit",
        lambda _, xdata, ydata, p0, **__: (np.asarray(p0), np.eye(len(p0))),
    )

    def residuals(**kwargs: object) -> tuple[float, float, float, pd.DataFrame]:
        table = kwargs["orderPixelTable"].copy()
        table["residuals_x"] = 0.1
        table["residuals_y"] = 0.1
        table["residuals_xy"] = np.hypot(0.1, 0.1)
        if len(table) > 3:
            table.loc[table.index[-1], "residuals_x"] = 10.0
            table.loc[table.index[-1], "residuals_xy"] = np.hypot(10.0, 0.1)
        return 0.1, 0.0, 0.1, table

    mapper.calculate_residuals = residuals

    xcoeff, ycoeff, fitted, clipped = mapper.fit_polynomials(
        source,
        wavelengthDeg=0,
        orderDeg=0,
        slitDeg=0,
    )

    np.testing.assert_allclose(xcoeff, [0.0])
    np.testing.assert_allclose(ycoeff, [0.0])
    assert len(fitted) == 3
    assert sorted(clipped["wavelength"].tolist()) == [503.0, 504.0]


def test_fit_polynomials_clips_an_entire_multi_pinhole_arc_line_set(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Reject every pinhole measurement for an anomalous arc-line group."""
    mapper = _fitting_mapper(log)
    mapper.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": 1.0,
        "poly-clipping-iteration-limit": 3,
        "poly-clipping-pinhole-sets": True,
    }
    source = pd.DataFrame(
        {
            "order": [10] * 6,
            "wavelength": [500.0, 500.0, 501.0, 501.0, 502.0, 502.0],
            "slit_position": [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0],
            "observed_x": [1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
            "observed_y": [1.0, 1.0, 2.0, 2.0, 3.0, 3.0],
            "dropped": [False] * 6,
            "R_pin": [4000.0] * 6,
        }
    )
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    optimizeModule = importlib.import_module("scipy.optimize")
    monkeypatch.setattr(commonutilsModule, "get_cached_coeffs", lambda **_: ([0.0], [0.0]))
    monkeypatch.setattr(
        optimizeModule,
        "curve_fit",
        lambda _, xdata, ydata, p0, **__: (np.asarray(p0), np.eye(len(p0))),
    )

    def residuals(**kwargs: object) -> tuple[float, float, float, pd.DataFrame]:
        table = kwargs["orderPixelTable"].copy()
        table["residuals_x"] = np.where(table["wavelength"] == 502.0, 20.0, 0.1)
        table["residuals_y"] = 0.1
        table["residuals_xy"] = np.hypot(table["residuals_x"], table["residuals_y"])
        return 0.1, 0.0, 0.1, table

    mapper.calculate_residuals = residuals

    _, _, retained, clipped = mapper.fit_polynomials(
        source,
        wavelengthDeg=0,
        orderDeg=0,
        slitDeg=0,
    )

    assert retained["wavelength"].tolist() == [500.0, 500.0, 501.0, 501.0]
    assert clipped["wavelength"].tolist() == [502.0, 502.0]


def test_qc_plot_with_an_arc_frame_measures_slit_geometry(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Arc-frame QC includes the order-edge slit and inter-order diagnostics."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.kw = lambda key: key
    mapper.detectorParams = {"rotate-qc-plot": 0, "flip-qc-plot": False}
    mapper.pinholeFrameMasked = CCDData(
        np.arange(72, dtype=float).reshape(6, 12),
        unit=u.adu,
        mask=np.zeros((6, 12), dtype=bool),
    )
    mapper.arcFrame = CCDData(
        np.ones((6, 12)), unit=u.adu, meta={"SLIT_VIS": "SLIT0.5x11"}
    )
    mapper.firstGuessMap = False
    mapper.debug = False
    mapper.axisA = "x"
    mapper.axisB = "y"
    mapper.meanFrameFlux = 20.0
    mapper.stdFrameFlux = 5.0
    mapper.sofName = "ARC"
    mapper.qcDir = str(tmp_path)
    mapper.orderTable = "synthetic-order-table.fits"
    mapper.products = pd.DataFrame()
    mapper.recipeName = "soxs-disp-solution"
    mapper.dateObs = "2025-01-01T00:00:00"
    mapper.qc = pd.DataFrame(
        columns=["qc_name", "qc_value", "qc_unit", "qc_comment"]
    )
    mapper.recipeSettings = {"sample_setting": 1}
    mapper.exptime = 1.0
    mapper.settings = {"tune-pipeline": False}
    lines = pd.DataFrame(
        {
            "order": [10, 10, 11, 11],
            "wavelength": [500.0, 501.0, 600.0, 601.0],
            "observed_x": [1.0, 2.0, 3.0, 4.0],
            "observed_y": [1.0, 2.0, 2.0, 3.0],
            "detector_x": [1.1, 2.1, 3.1, 4.1],
            "detector_y": [1.1, 2.1, 2.1, 3.1],
            "detector_x_shifted": [1.0, 2.0, 3.0, 4.0],
            "detector_y_shifted": [1.0, 2.0, 2.0, 3.0],
            "fit_x": [1.05, 2.05, 3.05, 4.05],
            "fit_y": [1.05, 2.05, 2.05, 3.05],
            "residuals_x": [0.05, -0.05, 0.1, -0.1],
            "residuals_y": [0.05, -0.05, -0.1, 0.1],
            "fwhm_slit_px": [2.0, 2.1, 1.9, 2.0],
            "R_slit": [4000.0, 4100.0, 3900.0, 4050.0],
        }
    )
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    geometry = pd.DataFrame(
        {
            "order": [10, 11],
            "xcoord_edgeup": [4.0, 5.0],
            "xcoord_edgelow": [0.0, 1.0],
            "xcoord_centre": [2.0, 3.0],
            "ycoord": [1.0, 1.0],
        }
    )
    mapData = pd.DataFrame(
        {
            "order": [10, 10, 11, 11],
            "x": [0.0, 4.0, 1.0, 5.0],
            "y": [1.0, 1.0, 1.0, 1.0],
            "wavelength": [500.0, 501.0, 600.0, 601.0],
            "slit_position": [-1.0, 1.0, -1.0, 1.0],
            "flux": [1.0] * 4,
            "pixelScale": [1.0] * 4,
        }
    )
    monkeypatch.setattr(
        module,
        "unpack_order_table",
        lambda **_: (pd.DataFrame(), geometry.copy(), pd.DataFrame()),
    )
    monkeypatch.setattr(
        module,
        "twoD_disp_map_image_to_dataframe",
        lambda **_: (mapData.copy(), np.zeros((1, 1), dtype=bool)),
    )

    result = mapper._create_dispersion_map_qc_plot(
        xcoeff=[1.0],
        ycoeff=[1.0],
        orderDeg=1,
        wavelengthDeg=1,
        slitDeg=0,
        orderPixelTable=lines,
        missingLines=lines.iloc[:1].copy(),
        allClippedLines=lines.assign(dropped=False),
        dispMapImage=False,
    )

    assert result == "ARC_RESIDUALS_110.pdf"
    assert (tmp_path / result).is_file()


def test_metric_clipping_marks_line_sets_with_each_pre_fit_diagnostic(
    log: object,
) -> None:
    """Multi-pinhole line metrics add the established clipping diagnostic flags."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.firstGuessMap = True
    mapper.debug = False
    rows = []
    for wavelength in range(500, 508):
        for slitIndex in (0, 1):
            rows.append(
                {
                    "wavelength": float(wavelength),
                    "order": 10,
                    "ion": "Ar",
                    "dropped": False,
                    "fwhm_pin_px": 2.0 + slitIndex * 0.1,
                    "flux": 100.0 + wavelength,
                    "peak": 200.0 + wavelength,
                    "x_diff": 0.1 * slitIndex,
                    "y_diff": 0.2 * slitIndex,
                    "xy_diff": 0.3 * slitIndex,
                }
            )
    rows[-1]["fwhm_pin_px"] = 30.0

    result = mapper._clip_on_measured_line_metrics(pd.DataFrame(rows))

    assert {
        "droppedOnScatter",
        "droppedOnFWHM",
        "droppedOnFlux",
        "droppedOnPeak",
    }.issubset(result.columns)
    assert result["dropped"].dtype == bool


@pytest.mark.parametrize("returnAll", [True, False])
def test_measure_line_position_detects_a_synthetic_arc_feature(
    log: object,
    returnAll: bool,
) -> None:
    """Photutils line detection returns detector coordinates and source metrics."""
    y, x = np.mgrid[:11, :11]
    data = 1000.0 * np.exp(-((x - 5.0) ** 2 + (y - 5.0) ** 2) / 2.0)
    stamp = CCDData(data, unit=u.adu, mask=np.zeros_like(data, dtype=bool))

    detected = measure_line_position(
        (stamp, 10, 21, 20, 31),
        log=log,
        windowHalf=5,
        iraf=False,
        sigmaLimit=3,
        iteration=0,
        returnAll=returnAll,
    )

    assert len(detected) == 1
    assert detected[0]["observed_x"] == pytest.approx(15.0, abs=0.2)
    assert detected[0]["observed_y"] == pytest.approx(25.0, abs=0.2)
    assert detected[0]["flux"] > 0
