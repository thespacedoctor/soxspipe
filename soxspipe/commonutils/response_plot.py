#!/usr/bin/env python
"""
*QC figure for the spline fit of the instrument response*

:Author:
    David Young
"""

import os
from typing import Any

os.environ["TERM"] = "vt100"

PLOT_RESIDUAL_LIMIT = 0.3
# HEADROOM ABOVE THE HIGHEST BINNED RESPONSE ON THE RATIO PANELS
PLOT_RATIO_HEADROOM = 1.5


def _shade_masked_regions(axis: Any, fitSettings: Any) -> None:
    """Shade the telluric regions and stellar-line windows that were left out of the fit."""
    halfWidth = fitSettings["line_half_width_nm"]
    for lower, upper in fitSettings["telluric_regions"]:
        axis.axvspan(lower, upper, color="orange", alpha=0.2, linewidth=0)
    for line in fitSettings["stellar_lines"]:
        axis.axvspan(line - halfWidth, line + halfWidth, color="purple", alpha=0.15, linewidth=0)


def _style_axis(axis: Any, title: str, xlabel: str, ylabel: str) -> None:
    axis.set_title(title, fontsize=12)
    axis.set_xlabel(xlabel, fontsize=9)
    axis.set_ylabel(ylabel, fontsize=9)
    axis.tick_params(axis="both", which="major", labelsize=9)


def _plot_fit_panels(raw: Any, fitPanel: Any, fit: Any, fitSettings: Any, stdExtWave: Any) -> None:
    """Draw the raw ratio and the binned fit with its model curve."""
    import numpy as np

    finite = np.isfinite(fit.ratio) & np.isfinite(fit.wave)
    model = fit.evaluate(stdExtWave)
    # SCALE TO THE BINNED RESPONSE: TELLURIC BANDS AND NEAR-ZERO COUNTS SEND THE RAW RATIO OFF-SCALE
    binnedResponse = np.exp(fit.binLog) if fit.binLog.size else np.array([1.0])
    topLimit = float(np.nanmax(binnedResponse)) * PLOT_RATIO_HEADROOM

    raw.scatter(fit.wave[finite], fit.ratio[finite], marker=".", s=1, c="grey", alpha=0.5)
    _shade_masked_regions(raw, fitSettings)
    raw.set_ylim(0, topLimit)
    _style_axis(raw, "Raw ratio", "wavelength (nm)", "Ratio $\\frac{F_{\\lambda}}{F_c}$")

    kept = fit.keptMask
    fitPanel.scatter(fit.wave[finite], fit.ratio[finite], marker=".", s=1, c="grey", alpha=0.3)
    fitPanel.plot(fit.binWave[kept], np.exp(fit.binLog[kept]), ".", c="blue", ms=4, label="kept bins")
    fitPanel.plot(fit.binWave[~kept], np.exp(fit.binLog[~kept]), "x", c="red", ms=5, label="clipped bins")
    fitPanel.plot(stdExtWave, model, c="red", linewidth=1.5, alpha=0.8, label="response curve")
    _shade_masked_regions(fitPanel, fitSettings)
    fitPanel.set_ylim(0, topLimit)
    fitPanel.legend(fontsize=7)
    _style_axis(fitPanel, "Fitted response function", "wavelength (nm)", "absolute-extracted flux ratio")


def _plot_selfcal_panels(selfCal: Any, residual: Any, fit: Any, stdExtWave: Any, stdExtFlux: Any, absFlux: Any) -> None:
    """Draw the self-calibrated standard and its relative residuals against the catalogue."""
    fluxCalib = stdExtFlux * fit.evaluate(stdExtWave) * 10**-17
    selfCal.plot(stdExtWave, fluxCalib, linewidth=0.2)
    _style_axis(selfCal, "Self calibration of std star", "wavelength (nm)", "flux ($\\mathrm{erg/cm^{2}/s/angstom}$)")

    residual.plot(stdExtWave, (fluxCalib * 10**17 - absFlux) / absFlux, linewidth=0.2)
    residual.axhline(0.0, color="black", linewidth=0.5)
    residual.set_ylim(-PLOT_RESIDUAL_LIMIT, PLOT_RESIDUAL_LIMIT)
    _style_axis(residual, "Relative residuals", "wavelength (nm)", "residual")


def _plot_efficiency_panel(axis: Any, stdExtWaveNotFlat: Any, stdEfficiencyEstimate: Any) -> None:
    """Draw the end-to-end efficiency estimate."""
    import numpy as np

    axis.plot(stdExtWaveNotFlat, stdEfficiencyEstimate, linewidth=0.2)
    axis.set_ylim(0, np.min([np.nanmax(stdEfficiencyEstimate), 1.0]))
    _style_axis(axis, "Efficiency (end-to-end)", "wavelength (nm)", "Efficiency")


def _title_qc_warnings(fig: Any, qcWarnings: Any) -> None:
    """Write the out-of-tolerance QC messages in red at the top of the figure."""
    if qcWarnings:
        fig.suptitle("QC OUT OF TOLERANCE\n" + "\n".join(qcWarnings), color="red", fontsize=8)


def draw_response_qc_figure(
    fit: Any,
    fitSettings: Any,
    stdExtWave: Any,
    stdExtFlux: Any,
    absFlux: Any,
    stdExtWaveNotFlat: Any,
    stdEfficiencyEstimate: Any,
    objName: str,
    qcWarnings: Any,
) -> Any:
    """*draw the response curve QC figure*

    **Key Arguments:**

    - ``fit`` -- the ``ResponseFit``
    - ``fitSettings`` -- the arm settings, which hold the masked regions
    - ``stdExtWave`` -- extracted standard star wavelengths (nm)
    - ``stdExtFlux`` -- extracted standard star flux (counts per second, extinction corrected)
    - ``absFlux`` -- catalogue flux of the standard at ``stdExtWave`` (1e-17 erg/cm2/s/A)
    - ``stdExtWaveNotFlat`` -- wavelengths (nm) of the efficiency estimate
    - ``stdEfficiencyEstimate`` -- the efficiency estimate, or None
    - ``objName`` -- name of the standard star
    - ``qcWarnings`` -- out-of-tolerance QC messages, drawn in red in the figure title

    **Return:**

    - ``fig`` -- the matplotlib figure
    """
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 22), constrained_layout=True, dpi=120)
    gs = fig.add_gridspec(30 if stdEfficiencyEstimate is not None else 25, 4)
    axes = [fig.add_subplot(gs[5 * i : 5 * i + 4, :]) for i in range(5)]
    if stdEfficiencyEstimate is not None:
        axes.append(fig.add_subplot(gs[25:29, :]))

    axes[0].plot(stdExtWave, absFlux * 10**-17, linewidth=0.2)
    _style_axis(
        axes[0],
        f"{objName} absolute flux spectrum",
        "wavelength (nm)",
        "flux ($\\mathrm{erg/cm^{2}/s/angstom}$)",
    )
    _plot_fit_panels(axes[1], axes[2], fit, fitSettings, stdExtWave)
    _plot_selfcal_panels(axes[3], axes[4], fit, stdExtWave, stdExtFlux, absFlux)

    if stdEfficiencyEstimate is not None:
        _plot_efficiency_panel(axes[5], stdExtWaveNotFlat, stdEfficiencyEstimate)

    _title_qc_warnings(fig, qcWarnings)
    return fig
