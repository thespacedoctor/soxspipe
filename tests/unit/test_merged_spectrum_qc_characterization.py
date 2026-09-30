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
from tests.factories import product_table
from tests.unit._plot_spies import spy_figures

pytestmark = pytest.mark.unit


def _synthetic_merged_orders(*, includeSky: bool = True) -> pd.DataFrame:
    """~50-row synthetic merged spectrum with a couple of outliers in FLUX and
    SNR, deterministic (no RNG) so sigma-clipped statistics can be
    recomputed independently in each test."""
    wave = np.linspace(500.0, 600.0, 50)
    flux = 10.0 + 2.0 * np.sin(np.linspace(0.0, 10.0, 50))
    flux[5] = 500.0
    flux[40] = -500.0
    snr = 5.0 + 0.1 * np.arange(50)
    snr[10] = 200.0
    columns = {"WAVE": wave, "FLUX_COUNTS": flux, "SNR": snr}
    if includeSky:
        columns["SKY_COUNTS"] = 10.0 + np.linspace(0.0, 5.0, 50)
    return pd.DataFrame(columns)


def _stub_skylines(monkeypatch: pytest.MonkeyPatch, extra: tuple = ()) -> None:
    wavelengths = [520.0, 540.0, *extra]
    isolated = [True, False] + [False] * len(extra)
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame({"WAVELENGTH": wavelengths, "ISOLATED": isolated}),
    )


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


def _call_plot_merged_spectrum_qc(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: object,
    merged: pd.DataFrame,
    **overrides: object,
) -> tuple[pd.DataFrame, str]:
    kwargs: dict[str, object] = {
        "products": product_table().iloc[0:0].copy(),
        "log": log,
        "qcDir": str(tmp_path),
        "filenameTemplate": "synthetic.fits",
        "noddingSequence": "_AB",
        "dateObs": "2024-01-02T03:04:05",
        "arm": "VIS",
        "recipeName": "soxs-nod",
        "orderJoins": {"11-12": 550.0},
        "fluxCalibrated": False,
        "qcTable": False,
        "settings": {},
    }
    kwargs.update(overrides)
    return toolkit.plot_merged_spectrum_qc(merged, **kwargs)


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_figure_layout_size_dpi_and_axis_order(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    """The figure is 14x10in at 180dpi with exactly four axes, created in
    top/middle/bottom(SNR)/sky order."""
    figures = spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

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
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

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
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()
    merged["WAVE"] = pd.Series([value * u.nm for value in merged["WAVE"]])

    # THE TypeError FIRES AFTER THE FIGURE IS CREATED BUT BEFORE THE FUNCTION'S
    # OWN plt.close("all"), SO CLOSE IT HERE TO KEEP LATER TESTS ISOLATED
    try:
        with pytest.raises(TypeError):
            _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged)
    finally:
        plt.close("all")


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_ylim_and_yscale_pin_sigma_clipped_statistics(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    figures = spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

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
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders(includeSky=False)

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

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
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

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


def test_plot_merged_spectrum_qc_snr_text_sorts_global_vis_letters_then_numeric_orders(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The SNR text box sorts `GLOBAL` first, then VIS arm letters in
    `_vis_rank` order, then numeric orders ascending, then any order name
    that is neither a VIS letter nor numeric (the `_order_key` fallback
    `except` branch); duplicate `(qc_name, qc_order)` rows keep the *last*
    occurrence's value."""
    figures = spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()
    qcTable = pd.DataFrame(
        {
            "qc_name": ["SNR MEDIAN"] * 8,
            "qc_order": [np.nan, "12", "11", "12", "u", "g", "r", "extra"],
            # THE SECOND "12" ROW (VALUE 9.0) MUST WIN OVER THE FIRST (5.0)
            "qc_value": [6.0, 5.0, 7.0, 9.0, 1.0, 2.0, 3.0, 4.0],
        }
    )

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=qcTable, orderJoins=False)

    bottom = figures[-1].axes[2]
    snrTexts = [text for text in bottom.texts if ":" in text.get_text()]
    assert len(snrTexts) == 1
    assert snrTexts[0].get_text() == "GLOBAL: 6\nu: 1\ng: 2\nr: 3\n11: 7\n12: 9\nextra: 4"


def test_plot_merged_spectrum_qc_snr_text_absent_when_qc_table_is_false(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    figures = spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=False, orderJoins=False)

    bottom = figures[-1].axes[2]
    assert len(bottom.texts) == 0


def test_plot_merged_spectrum_qc_snr_text_absent_when_no_snr_median_rows(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    figures = spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()
    qcTable = pd.DataFrame({"qc_name": ["OTHER"], "qc_order": ["GLOBAL"], "qc_value": [1.0]})

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=qcTable, orderJoins=False)

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
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins={"11-12": 550.0})

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
        lambda log, settings, arm: pd.DataFrame(
            {"WAVELENGTH": [520.0, 540.0], "ISOLATED": [True, False]}
        ),
    )
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

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
        lambda log, settings, arm: pd.DataFrame(
            {"WAVELENGTH": [520.0, "not-a-number"], "ISOLATED": [True, True]}
        ),
    )
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

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
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, debug=True)

    assert order == ["show", "savefig"]


@pytest.mark.parametrize(
    ("noddingSequence", "expectedSuffix"),
    [(False, ""), ("_AB", "_AB")],
)
def test_plot_merged_spectrum_qc_filename_suffix_reflects_nodding_sequence(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, noddingSequence, expectedSuffix: str
) -> None:
    spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _, outputPath = _call_plot_merged_spectrum_qc(
        monkeypatch, tmp_path, log, merged, noddingSequence=noddingSequence
    )

    assert outputPath == str(tmp_path / f"synthetic_EXTRACTED_MERGED_QC_PLOT{expectedSuffix}.pdf")
