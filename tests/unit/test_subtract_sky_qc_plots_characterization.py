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


def test_a_y_dispersion_order_draws_its_sky_model_panel_transposed(
    log: Any, tmp_path: Path, figures: list, capsys: pytest.CaptureFixture[str]
) -> None:
    """The order mask follows columns 0-3, but the sky-model values are written into rows 0-3."""
    # SUSPICIOUS: SKY-MODEL PANEL IS TRANSPOSED FOR Y-DISPERSION ARMS, FILED AS DY-597
    outputPath = tmp_path / "y_dispersion_plots"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis="y", rotate=False)

    filePath = subtractor.plot_sky_sampling(
        order=11, imageMapOrderDF=_order_strip(), knotLocations=np.array([502.0, 508.0])
    )

    assert filePath == f"{outputPath}/science_SKYMODEL_QC_PLOTS_ORDER_11.pdf"
    assert Path(filePath).read_bytes().startswith(b"%PDF")
    assert capsys.readouterr().out == "DEBUG: median 0.19867914959516253, std 0.7723658536192702\n"
    [figure] = figures
    assert len(figure.axes) == 9
    assert _unmasked_rows_and_columns(figure.axes[0].images[0].get_array()) == (list(range(16)), [0, 1, 2, 3])
    skyModelPanel = figure.axes[5].images[0].get_array()
    assert _unmasked_rows_and_columns(skyModelPanel) == (list(range(16)), [0, 1, 2, 3])
    dataRows, dataColumns = np.nonzero(skyModelPanel.data)
    assert sorted(set(dataRows.tolist())) == [0, 1, 2, 3]
    assert sorted(set(dataColumns.tolist())) == list(range(16))
    assert skyModelPanel.data[0, :4].tolist() == [100.0, 104.0, 108.0, 112.0]
    assert skyModelPanel.data[:4, 0].tolist() == [100.0, 101.0, 102.0, 103.0]
    assert float(skyModelPanel.data.sum()) == 8416.0


def test_a_y_dispersion_order_on_a_non_square_frame_raises_index_error(
    log: Any, tmp_path: Path, figures: list
) -> None:
    """A 20-pixel-wide order on a 16-row frame indexes past the last row once transposed."""
    # SUSPICIOUS: NON-SQUARE Y-DISPERSION FRAMES CRASH THE QC PLOT, FILED AS DY-597
    outputPath = tmp_path / "non_square_plots"
    outputPath.mkdir()
    subtractor = _subtractor(log, outputPath, dispersionAxis="y", rotate=False, shape=(16, 24))

    with pytest.raises(IndexError, match="index 18 is out of bounds for axis 0 with size 16"):
        subtractor.plot_sky_sampling(
            order=11, imageMapOrderDF=_order_strip(pixelCount=80, columns=20), knotLocations=np.array([502.0])
        )

    assert list(outputPath.iterdir()) == []


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
