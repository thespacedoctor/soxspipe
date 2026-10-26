"""Behavioral contracts for toolkit plotting paths."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils import toolkit
from tests.factories import product_table

pytestmark = pytest.mark.unit


def _frame(instrument: str = "SOXS") -> CCDData:
    data = np.arange(64, dtype=float).reshape(8, 8)
    return CCDData(
        data,
        unit=u.electron,
        mask=data % 7 == 0,
        uncertainty=StdDevUncertainty(np.ones((8, 8))),
        meta={"INSTRUME": instrument},
    )


def test_quicklook_returns_before_importing_plotting_when_disabled(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    monkeypatch.setattr(plt, "figure", lambda *args, **kwargs: pytest.fail("plotted"))

    assert toolkit.quicklook_image(log, _frame(), show=False) is None


@pytest.mark.parametrize("extension", ["data", "mask", "uncertainty"])
def test_quicklook_renders_each_ccd_extension_without_showing(tmp_path: Path, log: object, extension: str) -> None:
    outputPath = tmp_path / f"{extension}.pdf"

    result = toolkit.quicklook_image(
        log,
        _frame("XSHOOTER"),
        show=False,
        ext=extension,
        title=extension,
        saveToPath=outputPath,
    )

    assert result is None
    assert outputPath.read_bytes().startswith(b"%PDF")


@pytest.mark.parametrize("instrument", ["SOXS", "XSHOOTER", "OTHER"])
def test_quicklook_surface_plot_supports_each_orientation(tmp_path: Path, log: object, instrument: str) -> None:
    outputPath = tmp_path / f"surface-{instrument}.pdf"

    toolkit.quicklook_image(
        log,
        np.arange(64, dtype=float).reshape(8, 8) + 20,
        show=False,
        inst=instrument,
        surfacePlot=True,
        saveToPath=outputPath,
    )

    assert outputPath.is_file()


def test_quicklook_uses_default_instrument_for_plain_array(tmp_path: Path, log: object) -> None:
    outputPath = tmp_path / "array.pdf"

    toolkit.quicklook_image(
        log,
        np.zeros((8, 8)),
        show=False,
        saveToPath=outputPath,
    )

    assert outputPath.is_file()


def test_merged_spectrum_plot_skips_work_when_no_product_table_is_supplied(
    log: object,
) -> None:
    """QC plotting retains the no-product shortcut used by lightweight callers."""
    products, path = toolkit.plot_merged_spectrum_qc(
        pd.DataFrame(),
        False,
        log,
        "unused",
        "unused.fits",
        False,
        "2024-01-02T03:04:05",
        "VIS",
        "soxs-stare",
    )

    assert products is False
    assert path is None


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_merged_spectrum_plot_writes_qc_pdf_and_product_record(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    fluxCalibrated: bool,
) -> None:
    """Merged spectra generate a local plot with sky lines, joins, and SNR labels."""
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda *_: pd.DataFrame({"WAVELENGTH": [501.0, 503.0], "ISOLATED": [True, False]}),
    )
    merged = pd.DataFrame(
        {
            "WAVE": [500.0, 501.0, 502.0, 503.0, 504.0],
            "FLUX_COUNTS": [10.0, 12.0, 11.0, 13.0, 12.0],
            "SNR": [5.0, 6.0, 5.5, 7.0, 6.5],
            "SKY_COUNTS": [11.0, 12.0, 13.0, 12.0, 11.0],
        }
    )
    qc = pd.DataFrame(
        {
            "qc_name": ["SNR MEDIAN", "SNR MEDIAN O12", "SNR MEDIAN O11"],
            "qc_order": ["GLOBAL", 12, 11],
            "qc_value": [6.0, 5.0, 7.0],
        }
    )

    products, outputPath = toolkit.plot_merged_spectrum_qc(
        merged,
        product_table().iloc[0:0].copy(),
        log,
        str(tmp_path),
        "synthetic.fits",
        "_AB",
        "2024-01-02T03:04:05",
        "VIS",
        "soxs-nod",
        orderJoins={"11-12": 502.0},
        fluxCalibrated=fluxCalibrated,
        qcTable=qc,
        settings={},
    )

    assert Path(outputPath).read_bytes().startswith(b"%PDF")
    product = products.iloc[-1]
    assert product["soxspipe_recipe"] == "soxs-nod"
    assert product["label"] == "QC"
    assert product["file_path"] == outputPath
    assert ("FLUXCALIBRATED" in product["product_label"]) is fluxCalibrated


@pytest.mark.parametrize("useSkylines", [False, True])
def test_dispersion_grid_builds_order_edges_and_cross_lines(
    monkeypatch: pytest.MonkeyPatch, log: object, useSkylines: bool
) -> None:
    mapFrame = pd.DataFrame(
        {
            "order": [10, 10, 10, 11, 11, 11],
            "wavelength": [500.0, 501.0, 502.0, 600.0, 601.0, 602.0],
            "slit_position": [-1.0, 0.0, 1.0, -1.0, 0.0, 1.0],
        }
    )
    interOrderMask = np.array([[False, True]])
    monkeypatch.setattr(
        toolkit,
        "twoD_disp_map_image_to_dataframe",
        lambda **kwargs: (mapFrame, interOrderMask),
    )
    monkeypatch.setattr(
        toolkit,
        "dispersion_map_to_pixel_arrays",
        lambda **kwargs: kwargs["orderPixelTable"].assign(
            fit_x=lambda table: table["wavelength"],
            fit_y=lambda table: table["slit_position"],
        ),
    )
    skylines = pd.DataFrame({"WAVELENGTH": [501.0, 601.0]}) if useSkylines else False

    result, mask = toolkit.create_dispersion_solution_grid_lines_for_plot(
        log,
        "solution.fits",
        "image.fits",
        _frame(),
        lambda name: name,
        skylines=skylines,
        slitPositions=(-0.5, 0.5),
    )

    assert mask is interOrderMask
    assert result["order"].unique().tolist() == [10.0, 11.0]
    assert {"line", "fit_x", "fit_y"}.issubset(result.columns)
    if useSkylines:
        assert {501.0, 601.0}.issubset(result["wavelength"])


@pytest.mark.parametrize(
    ("shape", "expectedRatio"),
    [((50, 100), 0.5), ((100, 100), 1.0), ((300, 100), 1.5), ((1000, 100), 1.5)],
)
def test_qc_image_panel_height_scales_with_aspect_and_caps_tall_images(
    shape: tuple[int, int], expectedRatio: float
) -> None:
    height = toolkit.qc_image_panel_height(np.zeros(shape))

    assert height == pytest.approx(toolkit.QC_PLOT_PANEL_WIDTH_INCHES * expectedRatio)


@pytest.mark.parametrize(("qcRows", "expectedQcRows"), [(0, 1), (3, 4), (10, 11), (25, 11)])
def test_qc_table_heights_cap_qc_rows_and_pad_settings(qcRows: int, expectedQcRows: int) -> None:
    qc = pd.DataFrame({"qc_name": [f"QC{i}" for i in range(qcRows)]})
    settings = {"a": 1, "b": 2, "c": 3}

    qcHeight, settingsHeight = toolkit.qc_table_heights(qc, settings)

    assert qcHeight == pytest.approx(toolkit.QC_TABLE_ROW_HEIGHT_INCHES * expectedQcRows)
    assert settingsHeight == pytest.approx(toolkit.QC_TABLE_ROW_HEIGHT_INCHES * 5)


def test_legend_below_axis_anchors_legend_under_the_axis_with_given_font_size() -> None:
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="edge")

    legend = toolkit.legend_below_axis(ax, fig, fontsize=6)

    renderer = fig.canvas.get_renderer()
    fig.draw(renderer)
    assert legend.get_window_extent(renderer).y1 < ax.get_window_extent(renderer).y0
    assert legend.get_texts()[0].get_fontsize() == 6
    plt.close(fig)


def test_legend_below_axis_defaults_to_fontsize_four() -> None:
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="edge")

    legend = toolkit.legend_below_axis(ax, fig)

    assert legend.get_texts()[0].get_fontsize() == 4
    plt.close(fig)


def _figure_with_two_axes() -> plt.Figure:
    fig, axes = plt.subplots(2, 1, constrained_layout=True)
    axes[0].plot([0, 1], [0, 1])
    return fig


def test_add_notes_below_figure_adds_no_text_when_there_are_no_notes() -> None:
    fig = _figure_with_two_axes()

    toolkit.add_notes_below_figure(fig, [])
    toolkit.add_notes_below_figure(fig, None)

    assert fig.texts == []
    plt.close(fig)


def test_add_notes_below_figure_draws_one_text_per_note_in_order() -> None:
    fig = _figure_with_two_axes()

    toolkit.add_notes_below_figure(fig, ["WARNING: first", "ERROR: second"])

    assert [text.get_text() for text in fig.texts] == ["WARNING: first", "ERROR: second"]
    plt.close(fig)


def test_add_notes_below_figure_colours_errors_red_and_warnings_orange() -> None:
    from matplotlib.colors import to_rgb

    fig = _figure_with_two_axes()

    toolkit.add_notes_below_figure(fig, ["WARNING: w", "ERROR: e", "[A1] ERROR: tagged"])

    warningRgb, errorRgb, taggedErrorRgb = (to_rgb(text.get_color()) for text in fig.texts)
    assert errorRgb[0] > 0.7 and errorRgb[1] < 0.3 and errorRgb[2] < 0.3
    assert taggedErrorRgb == errorRgb
    assert warningRgb[0] > 0.7 and 0.4 < warningRgb[1] < 0.7 and warningRgb[2] < 0.2
    plt.close(fig)


def test_add_notes_below_figure_stacks_each_note_below_the_previous_one() -> None:
    fig = _figure_with_two_axes()
    toolkit.add_notes_below_figure(fig, ["WARNING: first", "WARNING: second", "WARNING: third"])
    renderer = fig.canvas.get_renderer()
    fig.draw(renderer)

    extents = [text.get_window_extent(renderer) for text in fig.texts]

    assert extents[0].y0 > extents[1].y0 > extents[2].y0
    plt.close(fig)


def test_add_notes_below_figure_wraps_a_long_note_onto_several_lines() -> None:
    fig = _figure_with_two_axes()
    longNote = "WARNING: " + "word " * 100

    toolkit.add_notes_below_figure(fig, [longNote])

    lines = fig.texts[0].get_text().split("\n")
    assert len(lines) > 1
    assert all(len(line) <= toolkit.NOTE_WRAP_WIDTH for line in lines)
    plt.close(fig)


def test_add_notes_below_figure_does_not_overlap_a_wrapped_note_with_the_next_note() -> None:
    fig = _figure_with_two_axes()
    toolkit.add_notes_below_figure(fig, ["WARNING: " + "word " * 100, "WARNING: next"])
    renderer = fig.canvas.get_renderer()
    fig.draw(renderer)

    first, second = (text.get_window_extent(renderer) for text in fig.texts)

    assert first.y0 >= second.y1 - 1e-6
    plt.close(fig)


def test_add_notes_below_figure_truncates_past_the_maximum_and_counts_the_rest() -> None:
    fig = _figure_with_two_axes()
    notes = [f"WARNING: note {i}" for i in range(toolkit.MAX_PLOT_NOTES + 3)]

    toolkit.add_notes_below_figure(fig, notes)

    texts = [text.get_text() for text in fig.texts]
    assert len(texts) == toolkit.MAX_PLOT_NOTES + 1
    assert texts[:-1] == notes[: toolkit.MAX_PLOT_NOTES]
    assert texts[-1] == "… and 3 more (see the recipe log)"
    plt.close(fig)


def test_add_notes_below_figure_keeps_every_note_when_exactly_at_the_maximum() -> None:
    fig = _figure_with_two_axes()
    notes = [f"WARNING: note {i}" for i in range(toolkit.MAX_PLOT_NOTES)]

    toolkit.add_notes_below_figure(fig, notes)

    assert [text.get_text() for text in fig.texts] == notes
    plt.close(fig)


def test_add_notes_below_figure_leaves_the_axes_count_unchanged() -> None:
    fig = _figure_with_two_axes()

    toolkit.add_notes_below_figure(fig, ["WARNING: first"])

    assert len(fig.axes) == 2
    plt.close(fig)


def test_add_notes_below_figure_extends_the_tight_bounding_box_to_include_the_notes() -> None:
    fig = _figure_with_two_axes()
    renderer = fig.canvas.get_renderer()
    heightBefore = fig.get_tightbbox(renderer).height

    toolkit.add_notes_below_figure(fig, ["WARNING: first", "ERROR: second"])

    assert fig.get_tightbbox(renderer).height > heightBefore
    plt.close(fig)


def _first_note_top_and_lowest_tight_bottom(fig: plt.Figure) -> tuple[float, float]:
    """Return the first note's top edge and the lowest axes tight-bbox bottom, in display pixels."""
    fig.draw_without_rendering()
    lowestAxes = min(ax.get_tightbbox().y0 for ax in fig.axes)
    noteTop = fig.texts[0].get_window_extent().y1
    return noteTop, lowestAxes


def _merged_style_figure() -> plt.Figure:
    fig = plt.figure(figsize=(14, 10), constrained_layout=True, dpi=180)
    gs = fig.add_gridspec(5, 1, height_ratios=[3, 1, 1, 1, 0])
    for row in range(4):
        ax = fig.add_subplot(gs[row, :])
        ax.plot([0, 1], [0, 1])
        ax.set_xlabel("wavelength (nm)")
    return fig


@pytest.mark.parametrize("figureBuilder", [_figure_with_two_axes, _merged_style_figure])
def test_add_notes_below_figure_puts_the_first_note_just_below_the_lowest_axes_tight_bbox(figureBuilder) -> None:
    fig = figureBuilder()

    toolkit.add_notes_below_figure(fig, ["ERROR: " + "word " * 80])

    noteTop, lowestAxes = _first_note_top_and_lowest_tight_bottom(fig)
    gapPoints = (lowestAxes - noteTop) * 72 / fig.dpi
    assert gapPoints >= -0.5
    assert gapPoints <= 20
    plt.close(fig)


def test_add_notes_below_figure_works_on_a_figure_with_a_non_agg_canvas() -> None:
    from matplotlib.backends.backend_pdf import FigureCanvasPdf

    fig = _merged_style_figure()
    FigureCanvasPdf(fig)

    toolkit.add_notes_below_figure(fig, ["ERROR: " + "word " * 80])

    assert len(fig.texts) == 1
    noteTop, lowestAxes = _first_note_top_and_lowest_tight_bottom(fig)
    assert noteTop <= lowestAxes + 0.5
    plt.close(fig)
