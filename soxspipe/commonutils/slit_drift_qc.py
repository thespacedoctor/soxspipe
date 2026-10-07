"""Slit-position drift QC: the object trace's slit position against wavelength, order by order.

Before rectification each order is centred on a polynomial fit of the object trace's slit position as a function of
wavelength. This module builds the plotting series for those fits and draws the QC figure.
"""

import numpy as np

SLIT_DRIFT_FIT_SAMPLES = 100
FALLBACK_SPAN_ALPHA = 0.15


def slit_drift_series(
    orderPixelTable,
    orders,
    centreCoeffs,
    wlMinMax,
    fallbackFlags,
    fitSamples=SLIT_DRIFT_FIT_SAMPLES,
):
    """*build the per-order trace, fit and residual arrays used by the slit-drift QC plot*

    **Key Arguments:**

    - ``orderPixelTable`` -- dataframe with ``order``, ``slit_position`` and ``wavelength`` columns
    - ``orders`` -- the order numbers, parallel to ``centreCoeffs``, ``wlMinMax`` and ``fallbackFlags``
    - ``centreCoeffs`` -- per-order slit-centre polynomial coefficients (highest power first, as for ``np.polyval``)
    - ``wlMinMax`` -- per-order ``(wlmin, wlmax)`` wavelength range in nm
    - ``fallbackFlags`` -- per-order flag, *True* where the order fell back to the mean slit position
    - ``fitSamples`` -- number of wavelengths at which each fit curve is sampled. Default *SLIT_DRIFT_FIT_SAMPLES*

    **Return:**

    - ``series`` -- list of dictionaries sorted blue to red, each with ``order``, ``isFallback``, ``traceWavelength``,
      ``traceSlit``, ``residual`` (trace point minus fit), ``fitWavelength`` and ``fitSlit``

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_drift_qc import slit_drift_series
    series = slit_drift_series(orderPixelTable, orders, centreCoeffs, wlMinMax, fallbackFlags)
    ```
    """
    series = []
    for order, coeffs, (wlmin, wlmax), isFallback in zip(orders, centreCoeffs, wlMinMax, fallbackFlags, strict=True):
        # SAME VALIDITY RULE AS THE TRANSFORMER: FINITE ON BOTH COLUMNS
        orderTrace = orderPixelTable.loc[orderPixelTable["order"] == order, ["slit_position", "wavelength"]]
        wavelength = orderTrace["wavelength"].to_numpy(dtype=float)
        slit = orderTrace["slit_position"].to_numpy(dtype=float)
        isFinite = np.isfinite(wavelength) & np.isfinite(slit)
        wavelength, slit = wavelength[isFinite], slit[isFinite]
        fitWavelength = np.linspace(wlmin, wlmax, fitSamples)
        series.append(
            {
                "order": order,
                "isFallback": bool(isFallback),
                "traceWavelength": wavelength,
                "traceSlit": slit,
                "residual": slit - np.polyval(coeffs, wavelength),
                "fitWavelength": fitWavelength,
                "fitSlit": np.polyval(coeffs, fitWavelength),
            }
        )
    return sorted(series, key=lambda entry: entry["fitWavelength"][0])


def _draw_order(mainAxis, residualAxis, entry, colour):
    """*draw one order's trace points, fit curve and residuals*"""
    mainAxis.scatter(entry["traceWavelength"], entry["traceSlit"], s=1, color=colour, rasterized=True)
    mainAxis.plot(entry["fitWavelength"], entry["fitSlit"], color=colour, lw=1)
    residualAxis.scatter(entry["traceWavelength"], entry["residual"], s=1, color=colour, rasterized=True)
    if not entry["isFallback"]:
        return
    wlmin, wlmax = entry["fitWavelength"][0], entry["fitWavelength"][-1]
    mainAxis.axvspan(wlmin, wlmax, color="red", alpha=FALLBACK_SPAN_ALPHA)
    mainAxis.text(
        (wlmin + wlmax) / 2,
        0.98,
        f"order {entry['order']}: no usable trace fit, fell back to mean",
        transform=mainAxis.get_xaxis_transform(),
        rotation=90,
        ha="center",
        va="top",
        fontsize=7,
        color="red",
    )


def plot_slit_drift_qc(series, meanSlitArcsec, title, filePath):
    """*plot the slit position of the object trace against wavelength, with the per-order fits and residuals*

    **Key Arguments:**

    - ``series`` -- the list returned by ``slit_drift_series``
    - ``meanSlitArcsec`` -- the mean slit position over all orders (the single centre used before per-order traces)
    - ``title`` -- the figure title
    - ``filePath`` -- the PDF path to write

    **Return:**

    - ``filePath`` -- the path the plot was written to

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_drift_qc import plot_slit_drift_qc
    plot_slit_drift_qc(series, meanSlitArcsec, "title", "/path/to/plot.pdf")
    ```
    """
    import matplotlib.pyplot as plt

    fig, (mainAxis, residualAxis) = plt.subplots(
        2, 1, sharex=True, figsize=(16, 8), gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05}
    )
    try:
        colours = plt.cm.rainbow(np.linspace(0, 1, max(len(series), 1)))
        for entry, colour in zip(series, colours, strict=False):
            _draw_order(mainAxis, residualAxis, entry, colour)
        mainAxis.axhline(
            meanSlitArcsec,
            ls="--",
            color="black",
            lw=1,
            label=f"mean slit position over all orders (former single centre) = {meanSlitArcsec:.3f} arcsec",
        )
        residualAxis.axhline(0, color="black", lw=0.5)
        mainAxis.set_ylabel("slit position (arcsec)", fontsize=10)
        residualAxis.set_ylabel("residual (arcsec)", fontsize=10)
        residualAxis.set_xlabel("wavelength (nm)", fontsize=10)
        mainAxis.set_title(title, fontsize=10)
        mainAxis.legend(fontsize=8, loc="best")
        fig.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    return filePath
