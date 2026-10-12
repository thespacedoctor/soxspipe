"""Slit-position drift QC: the object trace's slit position against wavelength, order by order.

Before rectification each order is centred on a polynomial fit of the object trace's slit position as a function of
wavelength. This module builds the plotting series for those fits and draws the QC figure.
"""

import numpy as np

SLIT_DRIFT_FIT_SAMPLES = 100
FALLBACK_SPAN_ALPHA = 0.15
SLIT_DRIFT_EDGE_TRIM_FRACTION = 0.05


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


def _valid_samples(entry):
    """*the fit samples of an order whose wavelength and slit position are both finite*"""
    wavelength, slit = entry["fitWavelength"], entry["fitSlit"]
    isValid = np.isfinite(wavelength) & np.isfinite(slit)
    return wavelength[isValid], slit[isValid]


def _trimmed_fit_value(wavelength, slit, edge):
    """*the fitted slit position at ``SLIT_DRIFT_EDGE_TRIM_FRACTION`` of the range in from the blue or red edge*"""
    wlMin, wlMax = wavelength.min(), wavelength.max()
    offset = SLIT_DRIFT_EDGE_TRIM_FRACTION * (wlMax - wlMin)
    atWavelength = wlMin + offset if edge == "blue" else wlMax - offset
    return float(np.interp(atWavelength, wavelength, slit))


def slit_drift_arcsec(series):
    """*the slit-position drift of the object trace, reddest minus bluest wavelength*

    The fit is read slightly inside the two ends of the wavelength range, ``SLIT_DRIFT_EDGE_TRIM_FRACTION`` of each
    end order's range in from its edge, as the polynomial fit is least reliable at the extreme edges. Orders that
    fell back to the mean slit position carry no trace information and are skipped, as are fit samples that are
    not finite. The end orders are the ones with the smallest and largest finite wavelength.

    **Key Arguments:**

    - ``series`` -- the list returned by ``slit_drift_series``

    **Return:**

    - ``driftArcsec`` -- fitted slit position near the red end of the reddest order minus that near the blue end of the
      bluest order, in arcsec (signed). *None* if every order fell back or the result is not finite

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_drift_qc import slit_drift_arcsec
    driftArcsec = slit_drift_arcsec(series)
    ```
    """
    samples = [_valid_samples(entry) for entry in series if not entry["isFallback"]]
    samples = [(wavelength, slit) for wavelength, slit in samples if wavelength.size]
    if not samples:
        return None
    # THE BLUE AND RED ENDS ARE CHOSEN BY WAVELENGTH, NOT BY POSITION IN THE LIST, AS ORDERS CAN BE NESTED
    blue = min(samples, key=lambda sample: sample[0].min())
    red = max(samples, key=lambda sample: sample[0].max())
    drift = _trimmed_fit_value(*red, "red") - _trimmed_fit_value(*blue, "blue")
    return float(drift) if np.isfinite(drift) else None


def _curve_zero_crossings(wavelength, slit):
    """*wavelengths where a sampled curve crosses zero, with linear interpolation between samples*

    A sample is valid only if its wavelength and slit position are both finite. Crossings are found only at and
    between consecutive valid samples, and a run of exact-zero valid samples gives one crossing.
    """
    isValid = np.isfinite(wavelength) & np.isfinite(slit)
    crossings = []
    for index in range(len(slit)):
        if not isValid[index]:
            continue
        wl, value = wavelength[index], slit[index]
        isPreviousValid = index > 0 and isValid[index - 1]
        # A RUN OF EXACT-ZERO SAMPLES IS ONE CROSSING, AT ITS BLUEST SAMPLE
        if value == 0 and not (isPreviousValid and slit[index - 1] == 0):
            crossings.append(float(wl))
        if isPreviousValid and value * slit[index - 1] < 0:
            previousWl, previousValue = wavelength[index - 1], slit[index - 1]
            fraction = previousValue / (previousValue - value)
            crossings.append(float(previousWl + fraction * (wl - previousWl)))
    return crossings


def slit_zero_crossings(series):
    """*the wavelengths where the fitted object trace crosses slit position 0*

    Crossings are found within each order only, by linear interpolation between fit samples. A sign change
    between one order and the next gives no crossing, as there is no fit between the orders to cross zero in, and
    orders are never concatenated, as overlapping orders would give false sign flips. Orders that fell back to the
    mean slit position are skipped.

    **Key Arguments:**

    - ``series`` -- the list returned by ``slit_drift_series``

    **Return:**

    - ``crossings`` -- crossing wavelengths in nm, sorted blue to red. Empty if the trace never crosses

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_drift_qc import slit_zero_crossings
    crossings = slit_zero_crossings(series)
    ```
    """
    crossings = []
    for entry in series:
        if not entry["isFallback"]:
            crossings.extend(_curve_zero_crossings(entry["fitWavelength"], entry["fitSlit"]))
    return sorted(crossings)


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


def _annotate_zero_crossing(mainAxis, residualAxis, zeroWavelength):
    """*mark the wavelength where the trace crosses slit position 0 on both panels*"""
    mainAxis.axhline(0, color="grey", lw=0.5)
    for axis in (mainAxis, residualAxis):
        axis.axvline(zeroWavelength, ls=":", color="grey", lw=1)
    mainAxis.text(
        zeroWavelength,
        0.02,
        f"trace crosses slit 0 at {zeroWavelength:.1f} nm",
        transform=mainAxis.get_xaxis_transform(),
        rotation=90,
        ha="right",
        va="bottom",
        fontsize=7,
        color="grey",
    )


def _annotate_metrics(mainAxis, residualAxis, driftArcsec, zeroWavelength):
    """*add the drift legend entry and the zero-crossing marks, for whichever metrics are given*"""
    if driftArcsec is not None:
        mainAxis.plot([], [], linestyle="none", label=f"slit drift (red - blue) = {driftArcsec:+.3f} arcsec")
    if zeroWavelength is not None:
        _annotate_zero_crossing(mainAxis, residualAxis, zeroWavelength)


def plot_slit_drift_qc(series, meanSlitArcsec, title, filePath, driftArcsec=None, zeroWavelength=None):
    """*plot the slit position of the object trace against wavelength, with the per-order fits and residuals*

    **Key Arguments:**

    - ``series`` -- the list returned by ``slit_drift_series``
    - ``meanSlitArcsec`` -- the mean slit position over all orders (the single centre used before per-order traces)
    - ``title`` -- the figure title
    - ``filePath`` -- the PDF path to write
    - ``driftArcsec`` -- the slit-position drift (red minus blue) to quote in the legend. Default *None* (not shown)
    - ``zeroWavelength`` -- the wavelength (nm) at which to mark the trace crossing slit position 0. Default *None*
      (not shown)

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
        _annotate_metrics(mainAxis, residualAxis, driftArcsec, zeroWavelength)
        mainAxis.set_ylabel("slit position (arcsec)", fontsize=10)
        residualAxis.set_ylabel("residual (arcsec)", fontsize=10)
        residualAxis.set_xlabel("wavelength (nm)", fontsize=10)
        mainAxis.set_title(title, fontsize=10)
        mainAxis.legend(fontsize=8, loc="best")
        fig.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    return filePath
