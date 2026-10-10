"""Order tables, extraction and the dispersion-map fit never write into a pandas slice (DY-1462).

Each test runs under ``mode.chained_assignment = "raise"`` and keeps any parent table alive, because pandas only flags
a write to a slice while the frame it was cut from still exists. The output assertions pin today's results, so they
also pass on the code before DY-1462 with the option off.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from soxspipe.commonutils.create_dispersion_map import create_dispersion_map
from soxspipe.commonutils.horne_extraction import compute_extractions
from soxspipe.commonutils.toolkit import unpack_order_table
from tests.factories import order_table_fits, workspace_organiser

pytestmark = pytest.mark.unit

PANDAS_DEFAULT_CHAINED_ASSIGNMENT = "warn"
RAISE = "raise"


def _order_table(tmp_path: Path) -> Path:
    """Two orders whose centre is x = 1 + 2y + 3 order + 4 order y, with a constant trace width of 2."""
    polynomials = pd.DataFrame(
        [
            {
                "degorder_cent": 1,
                "degy_cent": 1,
                "cent_00": 1.0,
                "cent_01": 2.0,
                "cent_10": 3.0,
                "cent_11": 4.0,
                "std_00": 2.0,
                "std_01": 0.0,
                "std_10": 0.0,
                "std_11": 0.0,
            }
        ]
    )
    metadata = pd.DataFrame({"order": [10, 11], "ymin": [0.0, 2.0], "ymax": [8.0, 10.0]})
    return order_table_fits(tmp_path / "orders.fits", polynomials=polynomials, metadata=metadata)


def _bare_mapper(log: object) -> create_dispersion_map:
    """A dispersion mapper with only a logger; each test sets the state its method reads."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    return mapper


@pytest.mark.parametrize("order", [None, 11])
def test_unpacking_a_binned_order_table_never_writes_into_a_slice(tmp_path: Path, log: object, order: int) -> None:
    # ARRANGE
    orderPath = _order_table(tmp_path)

    # ACT
    with pd.option_context("mode.chained_assignment", RAISE):
        _, pixelTable, metadataTable = unpack_order_table(
            log=log, orderTablePath=str(orderPath), pixelDelta=1, binx=2, biny=2, order=order
        )

    # ASSERT: ONLY THE INTEGER BINNED ROWS SURVIVE, AND BOTH AXES ARE SCALED BY THE BINNING
    order11 = pixelTable.loc[pixelTable["order"] == 11]
    assert list(order11["ycoord"]) == [1, 2, 3, 4]
    unbinnedY = order11["ycoord"] * 2
    np.testing.assert_allclose(
        order11["xcoord_centre"], (1 + 2 * unbinnedY + 3 * 11 + 4 * 11 * unbinnedY) / 2, rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(pixelTable["std"], 1.0, rtol=1e-12, atol=1e-12)
    meta11 = metadataTable.loc[metadataTable["order"] == 11].iloc[0]
    assert (meta11["ymin"], meta11["ymax"]) == (1.0, 5.0)
    assert set(pixelTable["order"]) == ({10, 11} if order is None else {11})


def test_extraction_drops_zero_wavelength_and_nan_rows_without_writing_into_a_slice() -> None:
    # ARRANGE: FOUR SLICES, ONE WITH ZERO WAVELENGTH AND ONE WITH NO PIXEL SCALE
    wavelengths = np.array([503.0, 501.0, -1.0, 502.0])
    profile = np.tile(np.array([[0.25], [0.5], [0.25]]), (1, wavelengths.size))
    orderImages = {
        "fluxRaw": profile * np.array([40.0, 80.0, 20.0, 60.0]),
        "objectProfile": profile,
        "variance": np.full_like(profile, 4.0),
        "wavelength": np.tile(wavelengths, (3, 1)),
        "mask": np.zeros_like(profile, dtype=bool),
    }
    slices = pd.DataFrame({"order": np.full(4, 12), "pixelScaleNm": [0.1, 0.1, 0.1, np.nan]})

    # ACT
    with pd.option_context("mode.chained_assignment", RAISE):
        result = compute_extractions(slices, orderImages, order=12)

    # ASSERT
    np.testing.assert_array_equal(result["wavelengthMean"], [501.0, 503.0])
    np.testing.assert_allclose(result["extractedFluxOptimal"], [80.0, 40.0], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(result["varianceSpectrum"], [32.0 / 3.0, 32.0 / 3.0], rtol=1e-12, atol=1e-12)


def test_first_guess_corrections_drop_incomplete_pinhole_sets_without_writing_into_a_slice(
    log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: TWO COMPLETE TWO-PINHOLE LINES AND ONE LINE SEEN THROUGH ONE PINHOLE ONLY
    mapper = _bare_mapper(log)
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

    # ACT
    with pd.option_context("mode.chained_assignment", RAISE):
        result = mapper._apply_first_guess_corrections(source)

    # ASSERT
    assert result["wavelength"].tolist() == [500.0, 500.0, 501.0, 501.0]
    assert result["detector_x_shifted"].tolist() == pytest.approx([9.5, 9.5, 19.5, 19.5])
    assert "droppedOnMissing" not in result


def test_dispersion_fit_of_a_table_with_dropped_lines_never_writes_into_a_slice(
    log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: ONE LINE IS ALREADY DROPPED AND ONE LINE HAS A LARGE RESIDUAL, AS IN A REAL ARC FIT
    mapper = _bare_mapper(log)
    mapper.arm = "VIS"
    mapper.settings = {}
    mapper.recipeName = "soxs-disp-solution"
    mapper.detectorParams = {}
    mapper.recipeSettings = {"poly-fitting-residual-clipping-sigma": 1.0, "poly-clipping-iteration-limit": 3}
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
    monkeypatch.setattr(
        importlib.import_module("soxspipe.commonutils"),
        "get_cached_coeffs",
        lambda **_: (np.array([0.0]), np.array([0.0])),
    )
    monkeypatch.setattr(
        importlib.import_module("scipy.optimize"),
        "curve_fit",
        lambda _, xdata, ydata, p0, **__: (np.asarray(p0), np.eye(len(p0))),
    )

    def residuals_in_place(**kwargs: object) -> tuple[float, float, float, pd.DataFrame]:
        # THE REAL calculate_residuals WRITES ITS COLUMNS ONTO THE TABLE IT IS GIVEN AND RETURNS THAT SAME TABLE
        table = kwargs["orderPixelTable"]
        table["residuals_x"] = 0.1
        table["residuals_y"] = 0.1
        table["residuals_xy"] = np.hypot(0.1, 0.1)
        if len(table) > 3:
            table.loc[table.index[-1], "residuals_x"] = 10.0
            table.loc[table.index[-1], "residuals_xy"] = np.hypot(10.0, 0.1)
        return 0.1, 0.0, 0.1, table

    mapper.calculate_residuals = residuals_in_place

    # ACT
    with pd.option_context("mode.chained_assignment", RAISE):
        _, _, fitted, clipped = mapper.fit_polynomials(source, wavelengthDeg=0, orderDeg=0, slitDeg=0)

    # ASSERT
    assert fitted["wavelength"].tolist() == [500.0, 501.0, 502.0]
    assert sorted(clipped["wavelength"].tolist()) == [503.0, 504.0]
    assert "sigma_clipped" not in source


def test_grouping_raw_frames_leaves_the_global_chained_assignment_option_at_the_pandas_default(
    tmp_path: Path, log: object
) -> None:
    # ARRANGE: THE OPTION STARTS AT THE PANDAS DEFAULT, AS IN A FRESH soxspipe reduce PROCESS
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.filterKeywordsExtras = ["mjd-obs"]
    rawFrames = pd.DataFrame(
        {
            "eso seq arm": ["VIS", "VIS"],
            "eso dpr type": ["BIAS", "BIAS"],
            "file": ["bias-2.fits", "bias-1.fits"],
            "filepath": ["./raw/bias-2.fits", "./raw/bias-1.fits"],
            "date-obs": ["2024-01-02T03:05:06.900", "2024-01-02T03:04:05.100"],
            "mjd-obs": [60311.12855, 60311.12784],
        }
    )

    # ACT
    with pd.option_context("mode.chained_assignment", PANDAS_DEFAULT_CHAINED_ASSIGNMENT):
        grouped = organiser._group_raw_frames(rawFrames, ["eso seq arm", "eso dpr type"])
        optionAfterGrouping = pd.get_option("mode.chained_assignment")

    # ASSERT
    assert optionAfterGrouping == PANDAS_DEFAULT_CHAINED_ASSIGNMENT
    assert grouped.loc[0, "counts"] == 2
    assert grouped.loc[0, "date-obs"] == "20240102T030405"


LINE_LIST_COLUMNS = [
    "wavelength", "order", "slit_index", "slit_position", "detector_x", "detector_y", "observed_x", "observed_y",
    "x_diff", "y_diff", "fit_x", "fit_y", "residuals_x", "residuals_y", "residuals_xy", "sigma_clipped", "sharpness",
    "roundness1", "roundness2", "npix", "sky", "peak", "flux", "fwhm_pin_px", "R_pin", "pixelScaleNm",
    "detector_x_shifted", "detector_y_shifted", "R_slit", "fwhm_slit_px",
]  # fmt: skip


def _line_list(wavelengths: list[float]) -> pd.DataFrame:
    """A fitted-line table holding every column the line-list writer keeps, plus one it drops."""
    table = pd.DataFrame({column: np.zeros(len(wavelengths)) for column in LINE_LIST_COLUMNS})
    table["wavelength"] = wavelengths
    table["order"] = 10
    table["sigma_clipped"] = False
    table["scratch"] = 1.0
    return table


@pytest.mark.parametrize("clippedWavelengths", [[], [501.5]])
def test_line_list_columns_are_trimmed_and_sorted_without_writing_into_the_fitted_table(
    log: object, clippedWavelengths: list[float]
) -> None:
    # ARRANGE: THE CALLER KEEPS ITS FITTED TABLE ALIVE, AS create_dispersion_map.get DOES
    mapper = _bare_mapper(log)
    fitted = _line_list([502.0, 500.0, 501.0])
    clipped = _line_list(clippedWavelengths)

    # ACT
    with pd.option_context("mode.chained_assignment", RAISE):
        goodAndClipped, good = mapper._prepare_line_list_columns(fitted, clipped)

    # ASSERT
    assert good["wavelength"].tolist() == [500.0, 501.0, 502.0]
    assert list(good.columns) == LINE_LIST_COLUMNS
    assert goodAndClipped["wavelength"].tolist() == sorted([500.0, 501.0, 502.0, *clippedWavelengths])
    assert goodAndClipped["sigma_clipped"].tolist() == [w in clippedWavelengths for w in goodAndClipped["wavelength"]]
    assert fitted["wavelength"].tolist() == [502.0, 500.0, 501.0]


def test_dispersion_map_writes_its_missing_lines_without_writing_into_the_detection_table(
    tmp_path: Path, log: object
) -> None:
    # ARRANGE: THREE PREDICTED LINES, ONE NOT DETECTED, WITH THE FIT AND THE OTHER WRITERS STUBBED OUT
    mapper = _bare_mapper(log)
    mapper.lineDetectionTable = pd.DataFrame(
        {
            "wavelength": [502.0, 500.0, 501.0],
            "order": [10, 10, 10],
            "slit_index": [0, 0, 0],
            "slit_position": [0.0, 0.0, 0.0],
            "detector_x": [12.0, 10.0, 11.0],
            "detector_y": [22.0, 20.0, 21.0],
            "observed_x": [np.nan, 100.0, np.nan],
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
    mapper.arm = "VIS"
    mapper.qcDir = str(tmp_path)
    fittedLines = mapper.lineDetectionTable.iloc[1:2].copy()
    mapper._initialize_recipe_settings = lambda: (False, False, 1, 1, 0)
    mapper.get_predicted_line_list = lambda: mapper.lineDetectionTable.copy()
    mapper._prepare_pinhole_frame = lambda: "prepared-frame"
    mapper._generate_quicklook_and_debug_plots = lambda *args: None
    mapper._clip_on_measured_line_metrics = lambda table: table
    mapper.fit_polynomials = lambda **kwargs: (np.array([1.0]), np.array([2.0]), fittedLines, fittedLines.iloc[0:0])
    mapper._get_output_filenames = lambda: ("good.fits", "missing.fits")
    mapper._write_qc_metrics = lambda *args: None
    mapper._write_line_list_qc = lambda *args: None
    mapper._prepare_line_list_columns = lambda *args: (fittedLines, fittedLines)
    mapper._write_fitted_lines_file = lambda *args: None
    mapper.write_map_to_file = lambda *args: str(tmp_path / "map.fits")
    mapper._create_dispersion_map_qc_plot = lambda **kwargs: ["residuals.pdf"]

    # ACT
    with pd.option_context("mode.chained_assignment", RAISE):
        mapper.get()

    # ASSERT: THE TWO UNDETECTED LINES ARE WRITTEN, SORTED BY WAVELENGTH
    missing = Table.read(tmp_path / "missing.fits").to_pandas()
    assert missing["wavelength"].tolist() == [501.0, 502.0]
    assert missing["detector_x"].tolist() == [11.0, 12.0]
