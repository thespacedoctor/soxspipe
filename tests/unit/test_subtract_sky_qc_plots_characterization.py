"""Characterization of the sky-sampling and image-comparison QC plots in subtract_sky (DY-40)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy.nddata import CCDData
from scipy.interpolate import splrep

import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.unit._plot_spies import quiet_show, spy_figures

pytestmark = pytest.mark.unit

FLAG_COLUMNS = (
    "flagged_all_clipped",
    "flagged_object_clipped",
    "flagged_bad_pixel_clipped",
    "flagged_edge_clipped",
    "flagged_bspline_clipped",
)


def _subtractor(
    log: Any, outputPath: Path, *, dispersionAxis: str, rotate: int | bool, shape: tuple[int, int] = (16, 16)
) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    # THE CONSTRUCTOR SWAPS THE AXES FOR A Y-DISPERSION ARM
    subtractor.axisA, subtractor.axisB = ("x", "y") if dispersionAxis == "x" else ("y", "x")
    subtractor.filenameTemplate = "science.fits"
    subtractor.qcDir = str(outputPath)
    subtractor.detectorParams = {
        "dispersion-axis": dispersionAxis,
        "flip-qc-plot": False,
        "rotate-qc-plot": rotate,
    }
    subtractor.objectFrame = CCDData(
        np.arange(np.prod(shape), dtype=float).reshape(shape),
        unit="electron",
        mask=np.zeros(shape, dtype=bool),
    )
    return subtractor


def _order_strip(
    *, pixelCount: int = 64, columns: int = 4, flagged: bool = True, largeScatter: bool = False
) -> pd.DataFrame:
    """An order `columns` pixels wide running down the rows, with every ninth pixel flagged."""
    index = np.arange(pixelCount)
    pixels = pd.DataFrame(
        {
            "x": index % columns,
            "y": index // columns,
            "wavelength": np.linspace(500.0, 510.0, pixelCount),
            "flux": 12.0 + np.sin(index),
            "slit_position": np.linspace(-1.0, 1.0, pixelCount),
            "residual_windowed_std": np.full(pixelCount, 3000.0 if largeScatter else 1.0),
            "residual_global_sigma": np.cos(index),
            "flux_percentile_smoothed": np.full(pixelCount, 12.0),
            "sky_model": 100.0 + index,
            "sky_subtracted_flux": np.sin(index / 3.0),
            "sky_subtracted_flux_weighted": np.sin(index / 3.0),
        }
    )
    for column in FLAG_COLUMNS:
        pixels[column] = (index % 9 == 0) if flagged else False
    return pixels


@pytest.fixture
def figures(monkeypatch: pytest.MonkeyPatch) -> list:
    # THE SURFACE QUICKLOOK OF THE OBJECT FRAME IS NOT PART OF THIS PLOT'S OUTPUT
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)
    quiet_show(monkeypatch)
    return spy_figures(monkeypatch)


def _unmasked_rows_and_columns(image: np.ma.MaskedArray) -> tuple[list[int], list[int]]:
    rows, columns = np.nonzero(~np.ma.getmaskarray(image))
    return sorted(set(rows.tolist())), sorted(set(columns.tolist()))


def _positions(pixels: np.ndarray) -> set[tuple[int, int]]:
    """The (row, column) positions of the True entries of a boolean array."""
    return {(int(row), int(column)) for row, column in np.argwhere(pixels)}


def _strip_positions(strip: pd.DataFrame) -> set[tuple[int, int]]:
    """The (row, column) detector positions of a strip's pixels."""
    return {(int(row), int(column)) for row, column in strip[["y", "x"]].to_numpy()}


def _assert_every_panel_follows_the_order(figure: Any, strip: pd.DataFrame) -> None:
    """Each image panel holds its data on the order's own (row, column) pixels."""
    rows, columns = sorted(set(strip["y"])), sorted(set(strip["x"]))
    assert _unmasked_rows_and_columns(figure.axes[0].images[0].get_array()) == (rows, columns)
    # IMAGES 1-4 OF THE CLIPPED-PIXEL PANEL ARE THE OBJECT, BAD-PIXEL, EDGE AND B-SPLINE OVERLAYS
    objectOverlay = figure.axes[3].images[1].get_array()
    assert _positions(~np.ma.getmaskarray(objectOverlay)) == _strip_positions(
        strip.loc[strip["flagged_object_clipped"]]
    )
    # THE LAST CLIPPED-PIXEL IMAGE MASKS THE ORDER ITSELF, LEAVING EVERYTHING ELSE VISIBLE
    orderOutline = figure.axes[3].images[-1].get_array()
    assert _positions(np.ma.getmaskarray(orderOutline)) == _strip_positions(strip)
    for panel in (figure.axes[5].images[0], figure.axes[6].images[0]):
        assert _unmasked_rows_and_columns(panel.get_array()) == (rows, columns)
    skyModelPanel = figure.axes[5].images[0].get_array()
    dataRows, dataColumns = np.nonzero(skyModelPanel.data)
    assert sorted(set(dataRows.tolist())) == rows
    assert sorted(set(dataColumns.tolist())) == columns
    assert skyModelPanel.data[0, :4].tolist() == [100.0, 101.0, 102.0, 103.0]
    assert skyModelPanel.data[:4, 0].tolist() == [100.0, 104.0, 108.0, 112.0]
    assert float(skyModelPanel.data.sum()) == 8416.0
    skySubPanel = figure.axes[6].images[0].get_array()
    assert skySubPanel.data[3, 2] == pytest.approx(np.sin(14 / 3.0))


@pytest.mark.parametrize("dispersionAxis", ["x", "y"])
def test_every_sky_qc_panel_places_its_pixels_on_the_true_detector_rows_and_columns(
    log: Any, tmp_path: Path, figures: list, capsys: pytest.CaptureFixture[str], dispersionAxis: str
) -> None:
    """An order in columns 0-3 fills columns 0-3 of every panel, whichever way the arm disperses."""
    outputPath = tmp_path / f"{dispersionAxis}_dispersion_plots"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis=dispersionAxis, rotate=False)
    strip = _order_strip()

    filePath = subtractor.plot_sky_sampling(order=11, imageMapOrderDF=strip, knotLocations=np.array([502.0, 508.0]))

    assert filePath == f"{outputPath}/science_SKYMODEL_QC_PLOTS_ORDER_11.pdf"
    assert Path(filePath).read_bytes().startswith(b"%PDF")
    assert capsys.readouterr().out == "DEBUG: median 0.19867914959516253, std 0.7723658536192702\n"
    [figure] = figures
    assert len(figure.axes) == 9
    _assert_every_panel_follows_the_order(figure, strip)


def test_a_y_dispersion_order_on_a_non_square_frame_renders_its_panels(log: Any, tmp_path: Path, figures: list) -> None:
    """A 20-pixel-wide order on a 16-row, 24-column y-dispersion frame plots without error."""
    outputPath = tmp_path / "non_square_plots"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis="y", rotate=False, shape=(16, 24))
    strip = _order_strip(pixelCount=80, columns=20)

    filePath = subtractor.plot_sky_sampling(order=11, imageMapOrderDF=strip, knotLocations=np.array([502.0]))

    assert Path(filePath).read_bytes().startswith(b"%PDF")
    [figure] = figures
    skyModelPanel = figure.axes[5].images[0].get_array()
    assert skyModelPanel.shape == (16, 24)
    assert _unmasked_rows_and_columns(skyModelPanel) == (list(range(4)), list(range(20)))
    assert skyModelPanel.data[3, 19] == 179.0


def test_a_rotated_clipped_pixel_panel_draws_every_layer_on_the_rotated_order(
    log: Any, tmp_path: Path, figures: list
) -> None:
    """With `rotate-qc-plot` set, every clipped-pixel panel layer is transposed like the raw frame.

    The rotation `flipud(rot90(image, 1))` is a transpose, so detector pixel
    (row, column) lands at (column, row). The order-outline layer must follow,
    or it whites out the order wherever the untransposed outline misses it.
    """
    outputPath = tmp_path / "rotated_clipped_pixel_panel"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis="x", rotate=90, shape=(16, 24))
    strip = _order_strip()
    transposedStrip = {(column, row) for row, column in _strip_positions(strip)}

    subtractor.plot_sky_sampling(order=11, imageMapOrderDF=strip, knotLocations=np.array([502.0, 508.0]))

    [figure] = figures
    clippedPixelPanel = figure.axes[3]
    # RAW FRAME, FOUR FLAG OVERLAYS, THEN THE ORDER OUTLINE LAST
    assert [image.get_array().shape for image in clippedPixelPanel.images] == [(24, 16)] * 6
    orderOutline = clippedPixelPanel.images[-1].get_array()
    assert _positions(np.ma.getmaskarray(orderOutline)) == transposedStrip


def test_a_rotated_spline_plot_without_clipped_rows_falls_back_on_every_limit(
    log: Any, tmp_path: Path, figures: list, capsys: pytest.CaptureFixture[str]
) -> None:
    """With no clipped rows, huge scatter and under 200 rows, each limit takes its fallback.

    The residual panel's std and median are NaN (replaced by 500 and 0), the
    flux panel floor drops below -3000 (replaced by -300) and the weighted
    residual limits are NaN, so `set_ylim` raises and is logged.
    """
    outputPath = tmp_path / "rotated_plots"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis="x", rotate=90)
    pixels = _order_strip(flagged=False, largeScatter=True)
    spline = splrep(pixels["wavelength"], pixels["flux"], k=3, s=len(pixels) * 4)

    filePath = subtractor.plot_sky_sampling(
        order=12, imageMapOrderDF=pixels, tck=spline, knotLocations=np.array([502.0, 508.0])
    )

    assert filePath == f"{outputPath}/science_SKYMODEL_QC_PLOTS_ORDER_12.pdf"
    assert capsys.readouterr().out == "DEBUG: median nan, std nan\n"
    [figure] = figures
    assert figure.axes[1].get_ylim() == (
        pytest.approx(-300.0, rel=1e-12, abs=0),
        pytest.approx(12.999911860107268, rel=1e-12, abs=0),
    )
    assert figure.axes[2].get_ylim() == (
        pytest.approx(-1500.0, rel=1e-12, abs=0),
        pytest.approx(3500.0, rel=1e-12, abs=0),
    )
    assert figure.axes[4].get_ylim() == (
        pytest.approx(-18.0, rel=1e-12, abs=0),
        pytest.approx(14.399999999999999, rel=1e-12, abs=0),
    )
    modelPanel = figure.axes[4]
    assert [line.get_label() for line in modelPanel.lines] == ["sky model"]
    assert modelPanel.lines[0].get_xdata().shape == (1000000,)
    assert [collection.get_label() for collection in modelPanel.collections] == ["unclipped pixels", "knots"]
    assert (
        "debug",
        "plot_sky_sampling: `ninerow.set_ylim(mean - 10 * std, mean +...` failed, continuing: "
        "Axis limits cannot be NaN or Inf",
    ) in log.messages
    assert figure._suptitle.get_text() == "VIS sky model: order 12"


def test_image_comparison_falls_back_to_plain_statistics_when_nan_statistics_fail(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failing `np.nanstd` on the object and sky-subtracted panels is logged and replaced."""
    outputPath = tmp_path / "comparison_plots"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis="x", rotate=False, shape=(4, 4))
    subtractor.mapDF = pd.DataFrame({"x": [0, 1], "y": [0, 1]})
    frame = subtractor.objectFrame
    realNanstd = np.nanstd
    calls: list[int] = []

    def fail_on_object_and_subtracted_panels(values: Any, *args: Any, **kwargs: Any) -> Any:
        calls.append(len(calls) + 1)
        # CALLS 1 AND 3 ARE INSIDE THE GUARDED BLOCKS; CALL 2 (SKY MODEL) IS NOT
        if len(calls) in (1, 3):
            raise TypeError("statistics unavailable")
        return realNanstd(values, *args, **kwargs)

    monkeypatch.setattr(np, "nanstd", fail_on_object_and_subtracted_panels)

    filePath = subtractor.plot_image_comparison(frame, frame.copy(), frame.copy())

    assert filePath == f"{outputPath}/science_skysub_quicklook.pdf"
    assert Path(filePath).read_bytes().startswith(b"%PDF")
    assert calls == [1, 2, 3]
    fallbackMessage = (
        "debug",
        "plot_image_comparison: `std = np.nanstd(maskedDataValues)` failed, continuing: statistics unavailable",
    )
    assert log.messages.count(fallbackMessage) == 2


@pytest.mark.parametrize("dispersionAxis", ["x", "y"])
def test_image_comparison_order_mask_covers_the_true_detector_pixels_of_the_order(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dispersionAxis: str
) -> None:
    """The statistics of the object panel are taken over the order's own (row, column) pixels."""
    outputPath = tmp_path / f"{dispersionAxis}_comparison_plots"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis=dispersionAxis, rotate=False, shape=(6, 8))
    subtractor.mapDF = pd.DataFrame({"x": [7, 6], "y": [1, 4]})
    frame = subtractor.objectFrame
    realNanstd = np.nanstd
    objectPanelValues: list[np.ndarray] = []

    def record_object_panel_values(values: Any, *args: Any, **kwargs: Any) -> Any:
        objectPanelValues.append(np.array(values, dtype=float))
        return realNanstd(values, *args, **kwargs)

    monkeypatch.setattr(np, "nanstd", record_object_panel_values)

    subtractor.plot_image_comparison(frame, frame.copy(), frame.copy())

    assert _positions(~np.isnan(objectPanelValues[0])) == {(1, 7), (4, 6)}
