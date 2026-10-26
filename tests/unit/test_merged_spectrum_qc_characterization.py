"""Characterization tests freezing the current rendering behaviour of
`plot_merged_spectrum_qc` in `soxspipe/commonutils/toolkit.py`, ahead of the DY-79 split of
that function into smaller helpers.

These tests pin current behaviour, including behaviour that looks like a
defect (noted inline as BUG-LIKE). They pin the matplotlib artists produced
-- figure sizes, axes, line data, colours, text and axis limits -- by spying
on `matplotlib.pyplot.figure` and inspecting the closed-but-inspectable
`Figure` objects it returns. The PDF output and the savefig/clf/close call
sequence are pinned in `tests/unit/test_toolkit_characterization.py`.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.stats import sigma_clip, sigma_clipped_stats

from soxspipe.commonutils import toolkit
from tests.unit._plot_spies import probe_on_savefig, spy_figures
from tests.unit._qc_plot_builders import call_plot_merged_spectrum_qc, stub_merged_skylines, synthetic_merged_orders

pytestmark = pytest.mark.unit


def _expected_flux_ylims(flux: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    """Independently recompute the top-panel and middle-panel ylims from the
    same sigma-clip/sigma-clipped-stats calls the production code uses."""
    arrayMask = sigma_clip(flux, sigma_lower=3, sigma_upper=15.0, maxiters=1, cenfunc="mean", stdfunc="std")
    _, _, std = sigma_clipped_stats(flux, sigma=5.0, stdfunc="std", cenfunc="mean", maxiters=3)
    topYlim = (arrayMask.min() - 3 * std, arrayMask.max() + 3 * std)
    middleYlim = (max(arrayMask.min() * 0.5, 0), arrayMask.max() * 2)
    return topYlim, middleYlim


def _expected_snr_ylim(snr: np.ndarray) -> tuple[float, float]:
    mean, _, std = sigma_clipped_stats(snr, sigma=5.0, stdfunc="std", cenfunc="mean", maxiters=3)
    return (0.0, mean + 4 * std)


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_figure_layout_size_dpi_and_axis_order(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    """The figure is 14x10in at 180dpi with exactly four axes, created in
    top/middle/bottom(SNR)/sky order."""
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    figure = figures[-1]
    assert tuple(figure.get_size_inches()) == (14.0, 10.0)
    assert figure.dpi == 180
    assert len(figure.axes) == 4
    top, middle, bottom, sky = figure.axes
    assert middle.get_yscale() == "log"
    assert sky.get_yscale() == "log"
    assert bottom.get_ylabel() == "SNR"
    assert sky.get_ylabel() == "sky flux ($e^{-}$)"


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_xlim_matches_wave_range_on_every_panel(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    for axis in figures[-1].axes:
        assert axis.get_xlim() == pytest.approx((500.0, 600.0), rel=1e-12)


def test_plot_merged_spectrum_qc_wave_as_quantity_column_raises_type_error(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """NOT PRACTICAL AS A `.value`-BRANCH CHARACTERIZATION TEST: a pandas
    column of per-row `astropy.units.Quantity` scalars (the shape the
    `try: ...min().value` fallback appears to defend against) cannot reach
    that fallback at all -- `ax.plot()` itself raises `TypeError` while
    building the line path, before any `set_xlim()` call runs. Pinned as a
    crash, not a `.value` codepath."""
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()
    merged["WAVE"] = pd.Series([value * u.nm for value in merged["WAVE"]])

    # THE TypeError FIRES AFTER THE FIGURE IS CREATED BUT BEFORE THE FUNCTION'S
    # OWN plt.close("all"), SO CLOSE IT HERE TO KEEP LATER TESTS ISOLATED
    try:
        with pytest.raises(TypeError):
            call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged)
    finally:
        plt.close("all")


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_ylim_and_yscale_pin_sigma_clipped_statistics(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    top, middle, bottom, sky = figures[-1].axes
    expectedTopYlim, expectedMiddleYlim = _expected_flux_ylims(merged["FLUX_COUNTS"].to_numpy())
    assert top.get_ylim() == pytest.approx(expectedTopYlim, rel=1e-12)
    assert middle.get_ylim() == pytest.approx(expectedMiddleYlim, rel=1e-12)
    assert middle.get_yscale() == "log"

    expectedBottomYlim = _expected_snr_ylim(merged["SNR"].to_numpy())
    assert bottom.get_ylim() == pytest.approx(expectedBottomYlim, rel=1e-12)

    expectedSkyYlim = (10.0, merged["SKY_COUNTS"].max() * 1.1)
    assert sky.get_ylim() == pytest.approx(expectedSkyYlim, rel=1e-12)
    assert sky.get_yscale() == "log"


def test_plot_merged_spectrum_qc_sky_panel_has_no_line_and_pinned_default_ylim_when_sky_counts_absent(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """When `SKY_COUNTS` is absent, the sky flux line is never drawn and the
    axis keeps matplotlib's own autoscaled default for an empty log-scale
    axis (recorded here as a literal, since it is not derived from the
    input data at all)."""
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders(includeSky=False)

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

    sky = figures[-1].axes[3]
    # ONLY THE TWO SKYLINE AXVLINES ARE PRESENT; NO SKY-FLUX LINE IS DRAWN
    assert len(sky.lines) == 2
    assert sky.get_ylim() == pytest.approx((0.8912509381337456, 11.220184543019636), rel=1e-9)


@pytest.mark.parametrize(
    ("fluxCalibrated", "expectedLabelSuffix", "expectedTopColour", "expectedSkyColour"),
    [
        (False, "flux ($e^{-}$)", "#dc322f", "#859900"),
        (True, "flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)", "#2aa198", "#2aa198"),
    ],
)
def test_plot_merged_spectrum_qc_labels_title_and_line_colours(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    fluxCalibrated: bool,
    expectedLabelSuffix: str,
    expectedTopColour: str,
    expectedSkyColour: str,
) -> None:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    top, middle, bottom, sky = figures[-1].axes
    assert top.get_ylabel() == expectedLabelSuffix
    assert top.get_title() == "Optimally Extracted Order-Merged Object Spectrum (VIS)\nsynthetic"
    assert bottom.get_xlabel() == "wavelength (nm)"

    topFluxLine = top.lines[0]
    assert topFluxLine.get_color() == expectedTopColour
    skyFluxLine = sky.lines[0]
    assert skyFluxLine.get_color() == expectedSkyColour
    snrLine = bottom.lines[0]
    assert snrLine.get_color() == "black"


def _snr_text_for(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object, qcTable: pd.DataFrame, arm: str) -> str:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)

    call_plot_merged_spectrum_qc(
        monkeypatch, tmp_path, log, synthetic_merged_orders(), qcTable=qcTable, orderJoins=False, arm=arm
    )

    bottom = figures[-1].axes[2]
    snrTexts = [text for text in bottom.texts if ":" in text.get_text()]
    assert len(snrTexts) == 1
    return snrTexts[0].get_text()


def test_plot_merged_spectrum_qc_snr_text_lists_soxs_vis_orders_blue_to_red_with_their_bands(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """GLOBAL comes first, then the VIS orders by wavelength (u, g, r, i), each labelled with its trace number
    and band; duplicate `(qc_name, qc_order)` rows keep the *last* occurrence's value."""
    qcTable = pd.DataFrame(
        {
            "qc_name": [
                "SNR MEDIAN",
                "SNR MEDIAN O4",
                "SNR MEDIAN O1",
                "SNR MEDIAN O3",
                "SNR MEDIAN O2",
                "SNR MEDIAN O3",
            ],
            "qc_order": [np.nan, 4, 1, 3, 2, 3],
            # THE SECOND O3 ROW (VALUE 9.0) MUST WIN OVER THE FIRST (5.0)
            "qc_value": [6.0, 4.0, 3.0, 5.0, 1.0, 9.0],
        }
    )

    text = _snr_text_for(monkeypatch, tmp_path, log, qcTable, arm="VIS")

    assert text == "GLOBAL: 6\n2 (u): 1\n3 (g): 9\n1 (r): 3\n4 (i): 4"


def test_plot_merged_spectrum_qc_snr_text_lists_other_arm_orders_by_number(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    qcTable = pd.DataFrame(
        {
            "qc_name": ["SNR MEDIAN", "SNR MEDIAN O12", "SNR MEDIAN O11"],
            "qc_order": [np.nan, 12, 11],
            "qc_value": [6.0, 5.0, 7.0],
        }
    )

    text = _snr_text_for(monkeypatch, tmp_path, log, qcTable, arm="NIR")

    assert text == "GLOBAL: 6\n11: 7\n12: 5"


def test_plot_merged_spectrum_qc_snr_text_absent_when_qc_table_is_false(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=False, orderJoins=False)

    bottom = figures[-1].axes[2]
    assert len(bottom.texts) == 0


def test_plot_merged_spectrum_qc_snr_text_absent_when_no_snr_median_rows(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()
    qcTable = pd.DataFrame({"qc_name": ["OTHER"], "qc_order": ["GLOBAL"], "qc_value": [1.0]})

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=qcTable, orderJoins=False)

    bottom = figures[-1].axes[2]
    assert len(bottom.texts) == 0


def test_plot_merged_spectrum_qc_order_joins_draw_axvline_and_text_on_three_panels_not_sky(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Each order join draws a dashed black `axvline` and an "ORDER JOIN"
    text on the top/middle/bottom panels -- but *not* the sky panel -- with
    the text placed at `(v + 5, 0.9 * ylim[1])` using each panel's own ylim
    at that point."""
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins={"11-12": 550.0})

    top, middle, bottom, sky = figures[-1].axes
    expectedTopYlim, expectedMiddleYlim = _expected_flux_ylims(merged["FLUX_COUNTS"].to_numpy())
    expectedBottomYlim = _expected_snr_ylim(merged["SNR"].to_numpy())

    for panel, expectedYlimTop in (
        (top, expectedTopYlim[1]),
        (middle, expectedMiddleYlim[1]),
        (bottom, expectedBottomYlim[1]),
    ):
        joinLines = [line for line in panel.lines if line.get_linestyle() == "--"]
        assert len(joinLines) == 1
        assert joinLines[0].get_color() == "black"
        assert joinLines[0].get_xdata()[0] == pytest.approx(550.0)

        joinTexts = [text for text in panel.texts if text.get_text() == "ORDER JOIN"]
        assert len(joinTexts) == 1
        x, y = joinTexts[0].get_position()
        assert x == pytest.approx(555.0)
        assert y == pytest.approx(0.9 * expectedYlimTop, rel=1e-12)

    assert not any(line.get_linestyle() == "--" for line in sky.lines)
    assert not any(text.get_text() == "ORDER JOIN" for text in sky.texts)


def test_plot_merged_spectrum_qc_skylines_isolated_vs_other_colour_and_label(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Isolated skylines draw blue lines labelled "calibration skyline" (only
    on the top panel); non-isolated ones draw grey lines labelled
    "skyline" (also only on the top panel)."""
    figures = spy_figures(monkeypatch)
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame({"WAVELENGTH": [520.0, 540.0], "ISOLATED": [True, False]}),
    )
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

    top = figures[-1].axes[0]
    isolatedLine = next(line for line in top.lines if line.get_xdata()[0] == pytest.approx(520.0))
    otherLine = next(line for line in top.lines if line.get_xdata()[0] == pytest.approx(540.0))
    assert isolatedLine.get_color() == "blue"
    assert isolatedLine.get_label() == "calibration skyline"
    assert otherLine.get_color() == "grey"
    assert otherLine.get_label() == "skyline"

    middle = figures[-1].axes[1]
    middleIsolatedLine = next(line for line in middle.lines if line.get_xdata()[0] == pytest.approx(520.0))
    # AN EMPTY-STRING LABEL IS TREATED BY MATPLOTLIB AS "NO LABEL", SO IT COMES
    # BACK AS AN AUTO-GENERATED "_childN" NAME RATHER THAN "" -- PIN THAT
    assert middleIsolatedLine.get_label().startswith("_child")


def test_plot_merged_spectrum_qc_skylines_drops_non_numeric_wavelengths(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`WAVELENGTH` values that fail `pd.to_numeric` are coerced to `NaN` and
    dropped, rather than raising or plotting at a bogus position."""
    figures = spy_figures(monkeypatch)
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame({"WAVELENGTH": [520.0, "not-a-number"], "ISOLATED": [True, True]}),
    )
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

    top = figures[-1].axes[0]
    skylineLines = [line for line in top.lines if line.get_color() == "blue"]
    assert len(skylineLines) == 1
    assert skylineLines[0].get_xdata()[0] == pytest.approx(520.0)


def test_plot_merged_spectrum_qc_debug_calls_show_before_savefig(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    order: list[str] = []
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: order.append("show"))
    realSavefig = plt.savefig

    def _savefig_spy(*args: object, **kwargs: object) -> object:
        order.append("savefig")
        return realSavefig(*args, **kwargs)

    monkeypatch.setattr(plt, "savefig", _savefig_spy)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, debug=True)

    assert order == ["show", "savefig"]


@pytest.mark.parametrize(
    ("noddingSequence", "expectedSuffix"),
    [(False, ""), ("_AB", "_AB")],
)
def test_plot_merged_spectrum_qc_filename_suffix_reflects_nodding_sequence(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, noddingSequence, expectedSuffix: str
) -> None:
    spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()

    _, outputPath = call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, noddingSequence=noddingSequence)

    assert outputPath == str(tmp_path / f"synthetic_EXTRACTED_MERGED_QC_PLOT{expectedSuffix}.pdf")


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_prints_extraction_notes_under_the_figure(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    probes = probe_on_savefig(monkeypatch, lambda figure: [text.get_text() for text in figure.texts])
    stub_merged_skylines(monkeypatch)
    merged = synthetic_merged_orders()
    notes = ["[A1] WARNING: Order skipped: 12", "[B1] ERROR: spill"]

    call_plot_merged_spectrum_qc(
        monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated, extractionNotes=notes
    )

    assert probes == [notes]


def test_plot_merged_spectrum_qc_with_notes_keeps_the_four_axes(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)

    call_plot_merged_spectrum_qc(
        monkeypatch, tmp_path, log, synthetic_merged_orders(), extractionNotes=["WARNING: note"]
    )

    assert len(figures[-1].axes) == 4


@pytest.mark.parametrize("notes", [None, []])
def test_plot_merged_spectrum_qc_without_notes_draws_no_figure_text(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, notes: object
) -> None:
    probes = probe_on_savefig(monkeypatch, lambda figure: list(figure.texts))
    stub_merged_skylines(monkeypatch)

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, synthetic_merged_orders(), extractionNotes=notes)

    assert probes == [[]]


@pytest.mark.parametrize(
    ("overrides", "expectedTitle"),
    [
        ({}, "Optimally Extracted Order-Merged Object Spectrum (VIS)\nsynthetic"),
        ({"boxcarExtraction": True}, "Robust Boxcar Extracted Order-Merged Object Spectrum (VIS)\nsynthetic"),
    ],
    ids=["default-optimal", "boxcar-fallback"],
)
def test_plot_merged_spectrum_qc_title_names_the_extraction_method(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    overrides: dict[str, object],
    expectedTitle: str,
) -> None:
    figures = spy_figures(monkeypatch)
    stub_merged_skylines(monkeypatch)

    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, synthetic_merged_orders(), **overrides)

    assert figures[-1].axes[0].get_title() == expectedTitle
