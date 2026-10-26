"""Slit-position drift QC: the object trace's slit position against wavelength, order by order.

Before rectification each order is centred on a polynomial fit of the object trace's slit position as a function of
wavelength. This module builds the plotting series for those fits and draws the QC figure.
"""

import math

import numpy as np

from soxspipe.commonutils.slit_illumination_qc import CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES, read_rotator_angle

SLIT_DRIFT_FIT_SAMPLES = 100
FALLBACK_SPAN_ALPHA = 0.15
SLIT_DRIFT_EDGE_TRIM_FRACTION = 0.05
SECONDS_PER_MINUTE = 60.0
SECONDS_PER_HOUR = 3600.0
SECONDS_PER_DAY = 86400.0
MINUTES_PER_HOUR = 60.0
HOURS_PER_DAY = 24.0
OFFSET_DECIMALS = 1
FULL_TURN_DEGREES = 360.0
HALF_TURN_DEGREES = 180.0


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


def _signed(value):
    """*a value as signed text with one decimal, where a value that rounds to zero reads +0.0 and never -0.0*"""
    return f"{round(value, OFFSET_DECIMALS) + 0.0:+.{OFFSET_DECIMALS}f}"


def format_time_offset(seconds):
    """*a signed time offset as text, in minutes, hours or days according to its size*

    **Key Arguments:**

    - ``seconds`` -- the time offset in seconds

    **Return:**

    - ``text`` -- the signed offset with one decimal and a unit: ``min`` if the rounded offset is under an hour, ``h``
      if under a day, otherwise ``d``

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_drift_qc import format_time_offset
    text = format_time_offset(2520.0)  # "+42.0 min"
    ```
    """
    # THE UNIT IS CHOSEN AFTER ROUNDING, SO 3599.99 s READS +1.0 h AND NOT +60.0 min
    minutes = round(abs(seconds) / SECONDS_PER_MINUTE, OFFSET_DECIMALS)
    if minutes < MINUTES_PER_HOUR:
        return f"{_signed(seconds / SECONDS_PER_MINUTE)} min"
    hours = round(abs(seconds) / SECONDS_PER_HOUR, OFFSET_DECIMALS)
    if hours < HOURS_PER_DAY:
        return f"{_signed(seconds / SECONDS_PER_HOUR)} h"
    return f"{_signed(seconds / SECONDS_PER_DAY)} d"


def rotator_angle_offset(scienceAngle, pinholeAngle):
    """*the rotator angle of the science frame minus that of the multi-pinhole frame, wrapped to (-180, 180] degrees*

    **Key Arguments:**

    - ``scienceAngle`` -- the science frame rotator angle in degrees
    - ``pinholeAngle`` -- the multi-pinhole frame rotator angle in degrees

    **Return:**

    - ``offset`` -- the signed angle difference in degrees, taking the shorter way round the circle

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_drift_qc import rotator_angle_offset
    offset = rotator_angle_offset(179.0, -179.0)  # -2.0
    ```
    """
    # PYTHON'S MODULO TAKES THE SIGN OF THE DIVISOR, GIVING [-180, 180), SO +180 DEGREES IS MAPPED BACK FROM -180
    offset = (scienceAngle - pinholeAngle + HALF_TURN_DEGREES) % FULL_TURN_DEGREES - HALF_TURN_DEGREES
    if offset == -HALF_TURN_DEGREES:
        return HALF_TURN_DEGREES
    return offset


def _observation_mjd(header, kw):
    """*the observation time of a header as an MJD, from MJD-OBS or else DATE-OBS, or None if neither is usable*"""
    from astropy.time import Time

    try:
        value = header[kw("MJDOBS")]
        # A BOOLEAN IS NOT A DATE, EVEN THOUGH float() ACCEPTS IT, SO IT FALLS THROUGH TO DATE-OBS
        if not isinstance(value, bool | np.bool_) and math.isfinite(float(value)):
            return float(value)
    except (KeyError, TypeError, ValueError, OverflowError):
        pass
    try:
        dateObs = str(header[kw("DATE_OBS")])
        # A DATE WITH NO TIME OF DAY CANNOT GIVE A TIME OFFSET
        if "T" not in dateObs:
            return None
        return float(Time(dateObs, format="isot", scale="utc").mjd)
    except (KeyError, TypeError, ValueError):
        return None


def pinhole_offset_labels(scienceHeader, pinholeHeader, kw):
    """*the legend labels for the science frame's time and rotator-angle offsets from the multi-pinhole frame*

    An offset that cannot be worked out because a value is missing, unreadable or not finite is left out, and the
    other is still given. The exception is the multi-pinhole rotator angle: daytime calibration frames carry no
    ABSROT and are taken with the rotator parked at ``CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES``, so that angle is used
    when the multi-pinhole header has no usable one, and the label says it was assumed. A missing science angle still
    drops the rotator entry. This function never raises.

    **Key Arguments:**

    - ``scienceHeader`` -- the science frame header
    - ``pinholeHeader`` -- the multi-pinhole frame header
    - ``kw`` -- the keyword lookup, mapping an alias such as ``MJDOBS`` to the header keyword

    **Return:**

    - ``labels`` -- a list of up to two legend labels, the obs-time offset first and then the rotator-angle offset,
      both science minus multi-pinhole

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_drift_qc import pinhole_offset_labels
    labels = pinhole_offset_labels(scienceHeader, pinholeHeader, self.kw)
    ```
    """
    labels = []
    scienceMjd = _observation_mjd(scienceHeader, kw)
    pinholeMjd = _observation_mjd(pinholeHeader, kw)
    if scienceMjd is not None and pinholeMjd is not None:
        seconds = (scienceMjd - pinholeMjd) * SECONDS_PER_DAY
        if math.isfinite(seconds):
            labels.append(f"science - multi-pinhole obs time = {format_time_offset(seconds)}")
    rotatorLabel = _rotator_label(scienceHeader, pinholeHeader, kw("ABSROT"))
    if rotatorLabel is not None:
        labels.append(rotatorLabel)
    return labels


def _rotator_label(scienceHeader, pinholeHeader, key):
    """*the rotator-angle offset legend label, or None if the science angle or the offset is not usable*

    A multi-pinhole header with no usable rotator angle is a daytime calibration frame, so the parked angle
    ``CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES`` is used and the label says so.
    """
    scienceAngle = read_rotator_angle(scienceHeader, key)
    if scienceAngle is None:
        return None
    pinholeAngle = read_rotator_angle(pinholeHeader, key)
    isAssumed = pinholeAngle is None
    if isAssumed:
        pinholeAngle = CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES
    offset = rotator_angle_offset(scienceAngle, pinholeAngle)
    if not math.isfinite(offset):
        return None
    label = f"science - multi-pinhole rotator angle (ABSROT) = {_signed(offset)} deg"
    if isAssumed:
        label += f" (multi-pinhole ABSROT assumed {CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES:.1f} deg)"
    return label


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


def _annotate_metrics(mainAxis, residualAxis, driftArcsec, zeroWavelength, extraLabels):
    """*add the drift legend entry, any extra legend entries and the zero-crossing marks, for whichever are given*"""
    if driftArcsec is not None:
        mainAxis.plot([], [], linestyle="none", label=f"slit drift (red - blue) = {driftArcsec:+.3f} arcsec")
    for label in extraLabels:
        mainAxis.plot([], [], linestyle="none", label=label)
    if zeroWavelength is not None:
        _annotate_zero_crossing(mainAxis, residualAxis, zeroWavelength)


def plot_slit_drift_qc(
    series, meanSlitArcsec, title, filePath, driftArcsec=None, zeroWavelength=None, extraLabels=None
):
    """*plot the slit position of the object trace against wavelength, with the per-order fits and residuals*

    **Key Arguments:**

    - ``series`` -- the list returned by ``slit_drift_series``
    - ``meanSlitArcsec`` -- the mean slit position over all orders (the single centre used before per-order traces)
    - ``title`` -- the figure title
    - ``filePath`` -- the PDF path to write
    - ``driftArcsec`` -- the slit-position drift (red minus blue) to quote in the legend. Default *None* (not shown)
    - ``zeroWavelength`` -- the wavelength (nm) at which to mark the trace crossing slit position 0. Default *None*
      (not shown)
    - ``extraLabels`` -- further legend entries, each drawn as text with no marker. Default *None* (none added)

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
        _annotate_metrics(mainAxis, residualAxis, driftArcsec, zeroWavelength, extraLabels or ())
        mainAxis.set_ylabel("slit position (arcsec)", fontsize=10)
        residualAxis.set_ylabel("residual (arcsec)", fontsize=10)
        residualAxis.set_xlabel("wavelength (nm)", fontsize=10)
        mainAxis.set_title(title, fontsize=10)
        mainAxis.legend(fontsize=8, loc="best")
        fig.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    return filePath
