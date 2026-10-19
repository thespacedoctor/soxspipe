"""Behaviour contracts for the response QC figure."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from soxspipe.commonutils.response_fit import fit_response_curve
from soxspipe.commonutils.response_plot import draw_response_qc_figure

pytestmark = pytest.mark.unit

FIT_SETTINGS = {
    "telluric_regions": [[600.0, 720.0]],
    "stellar_lines": [],
    "line_half_width_nm": 2.5,
    "bin_width_nm": 5.0,
    "min_bin_samples": 5,
    "knot_spacing_nm": 60.0,
    "clip_sigma_low": 4.0,
    "clip_sigma_high": 2.5,
    "max_iteration": 10,
    "telluric_wing_max_nm": 60.0,
    "telluric_wing_ref_nm": 40.0,
    "telluric_wing_threshold": 0.03,
}


def _fit_with_wings_and_a_clipped_bin() -> object:
    """Fit a response with raised telluric wings and one outlying bin."""
    wave = np.linspace(350.0, 850.0, 10000)
    ratio = np.exp(1.0 - 0.002 * (wave - 350.0))
    distance = np.where(wave < 600.0, 600.0 - wave, np.where(wave > 720.0, wave - 720.0, np.inf))
    ratio = ratio * (1.0 + 0.25 * np.exp(-distance / 10.0))
    ratio[(wave > 450.0) & (wave < 455.0)] *= 3.0
    return fit_response_curve(wave, ratio, FIT_SETTINGS)


def _fit_panel(fit: object) -> object:
    wave = np.linspace(350.0, 850.0, 500)
    flux = np.full(wave.size, 1.0)
    fig = draw_response_qc_figure(fit, FIT_SETTINGS, wave, flux, flux, wave, None, "STAR", [])
    return fig.axes[2]


def test_wing_bins_get_their_own_legend_entry_and_marker_series() -> None:
    # ARRANGE
    fit = _fit_with_wings_and_a_clipped_bin()

    # ACT
    panel = _fit_panel(fit)
    series = {line.get_label(): line for line in panel.get_lines()}
    plt.close("all")

    # ASSERT
    assert "telluric wing bins" in series
    assert len(series["telluric wing bins"].get_xdata()) == int(fit.wingMask.sum())
    assert series["telluric wing bins"].get_color() != series["clipped bins"].get_color()


def test_the_clipped_bins_series_excludes_wing_bins() -> None:
    # ARRANGE
    fit = _fit_with_wings_and_a_clipped_bin()
    clippedCount = int(np.count_nonzero(~fit.keptMask & ~fit.wingMask))

    # ACT
    panel = _fit_panel(fit)
    series = {line.get_label(): line for line in panel.get_lines()}
    plt.close("all")

    # ASSERT
    assert fit.wingMask.any()
    assert clippedCount > 0
    assert len(series["clipped bins"].get_xdata()) == clippedCount
