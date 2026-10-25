"""The anchor geometry of the faint-order trace: slit scale, the clipped curve, the anchor edge and the binned fit.

``faint_order_trace`` decides which orders are faint and how their rows are replaced. The pieces it shares with the
anchor measurement live here: the signed slit scale of an order, the MAD-clipped curve of an order's centre shifts, the
slit position and width of an anchor at the edge it shares with a faint order, the search window, and the fit of one
wavelength bin. The module imports numpy only and takes the dispersion map as a ``pixelsAt`` function.
"""

from dataclasses import dataclass

import numpy as np

from soxspipe.commonutils.faint_order_profile import FIT_MARGIN_SIGMAS, collapse_bin_profile, fit_fixed_width_gaussian

# THE SMALLEST NUMBER OF ACCEPTED SAMPLES THAT FIXES A DEGREE-2 CURVE WITH ONE SAMPLE TO SPARE
MIN_FIT_SAMPLES = 4
# THE DEGREE OF THE CENTRE-SHIFT CURVE THAT SEPARATES SCATTER FROM A SMOOTH TRACE
SCATTER_FIT_DEGREE = 2
# MAD-CLIPPING OF THE CENTRE-SHIFT CURVE: SIGMA LIMIT AND ITERATION COUNT
CLIP_SIGMA = 3.0
CLIP_ITERATIONS = 5
# THE WIDTH OF A BRIGHT TRACE IS TAKEN FROM THIS FRACTION OF THE ORDER NEAREST THE SHARED EDGE, AND FROM THE WHOLE
# ORDER WHEN THAT HOLDS FEWER THAN `MIN_EDGE_SAMPLES` ACCEPTED SAMPLES
EDGE_FRACTION = 0.25
MIN_EDGE_SAMPLES = 3
# THE CLIPPING WIDTH NEVER FALLS BELOW THIS (PIXELS), SO NOISE-FREE DATA ARE NOT CLIPPED AT ROUNDING ERROR
MIN_CLIP_SIGMA = 1e-9
# THE SEARCH WINDOW AROUND THE ANCHOR IS THREE WIDTHS EITHER SIDE, NEVER UNDER 1.5 ARCSEC AND NEVER OVER 1.8 ARCSEC, SO
# IT CANNOT REACH A SOURCE 3.7 ARCSEC AWAY
WINDOW_WIDTHS = 3.0
WINDOW_MIN_ARCSEC = 1.5
WINDOW_MAX_ARCSEC = 1.8
# SEEING WIDENS AS WAVELENGTH TO THIS POWER, SO A WIDTH IS CARRIED ACROSS THE ORDER BOUNDARY WITH IT
WIDTH_EXPONENT = -0.2
# A BIN IS SAMPLED EVERY DISPERSION PIXEL, WITH AT LEAST THIS MANY SAMPLES
MIN_BIN_TRACK_POINTS = 2
# CONVERTS A MEDIAN ABSOLUTE DEVIATION TO A STANDARD DEVIATION FOR GAUSSIAN NOISE
MAD_TO_SIGMA = 1.4826


def clipped_polynomial_residuals(x, y, degree):
    """*the residuals of a MAD-clipped polynomial fit of ``y`` against ``x``*

    The fit is a linear least-squares solve on a rescaled ``x``, and it is repeated on the samples inside
    ``CLIP_SIGMA`` robust standard deviations of the previous curve. The residuals are for all the samples.

    **Key Arguments:**

    - ``x`` -- the abscissae, such as wavelengths
    - ``y`` -- the values to fit, such as centre shifts
    - ``degree`` -- the polynomial degree

    **Return:**

    - ``curve`` -- the fitted curve, a ``numpy.polynomial.Polynomial``
    - ``residuals`` -- ``y`` minus the curve, one per sample

    A ``ValueError`` is raised when there are fewer than ``MIN_FIT_SAMPLES`` samples.
    """
    if len(x) < MIN_FIT_SAMPLES:
        raise ValueError(f"A clipped polynomial fit needs at least {MIN_FIT_SAMPLES} samples, got {len(x)}.")
    keep = np.ones(len(x), dtype=bool)
    for _ in range(CLIP_ITERATIONS):
        curve = np.polynomial.Polynomial.fit(x[keep], y[keep], degree)
        residuals = y - curve(x)
        centre = np.median(residuals[keep])
        sigma = MAD_TO_SIGMA * np.median(np.abs(residuals[keep] - centre))
        newKeep = np.abs(residuals - centre) <= CLIP_SIGMA * max(sigma, MIN_CLIP_SIGMA)
        if newKeep.sum() < degree + 1 or np.array_equal(newKeep, keep):
            break
        keep = newKeep
    return curve, residuals


def pixels_per_arcsec(pixelsAt, order, wavelengths):
    """*the signed cross-dispersion pixels per slit arcsec of an order, from the map at slit +1 and -1 arcsec*

    The sign is negative where the slit runs against the pixel axis, as it does on the VIS detector.

    **Key Arguments:**

    - ``pixelsAt`` -- the function ``pixelsAt(order, wavelengths, slitArcsec)``, giving the (axisA, axisB) pixels
    - ``order`` -- the order number
    - ``wavelengths`` -- the wavelengths (nm) at which to measure the scale

    **Return:**

    - ``scale`` -- one signed scale per wavelength
    """
    upper, _ = pixelsAt(order, wavelengths, 1.0)
    lower, _ = pixelsAt(order, wavelengths, -1.0)
    return (np.asarray(upper, dtype=float) - np.asarray(lower, dtype=float)) / 2.0


def bright_edge_anchor(rows, edgeNm, pixelsAt, *, axisA):
    """*the slit position and width of an anchor order's trace at the wavelength edge it shares with a faint order*

    The slit position of each accepted sample is its measured centre minus the map's slit-zero pixel at the sample's
    wavelength, over the pixel scale. The ``fit_<axisA>`` and ``centre_shift`` columns cannot give it: on the second
    sampling pass ``detect_continuum`` moves ``fit_<axisA>`` of each order onto the order's median trace, so
    ``centre_shift`` is measured from that trace and not from the map's slit zero. The offset of the samples from the
    slit-zero pixel is fitted against wavelength with a MAD-clipped degree-2 curve and evaluated at ``edgeNm``. The
    width is the median Gaussian standard deviation of the accepted samples in the quarter of the order's wavelength
    range nearest the edge.

    **Key Arguments:**

    - ``rows`` -- the sampled rows of one anchor order
    - ``edgeNm`` -- the wavelength (nm) of the edge it shares with the faint order
    - ``pixelsAt`` -- the function ``pixelsAt(order, wavelengths, slitArcsec)``, giving the (axisA, axisB) pixels
    - ``axisA`` -- the cross-dispersion axis, ``"x"`` or ``"y"``

    **Return:**

    - ``slitArcsec`` -- the trace's slit position at the edge, in arcsec from the map's slit zero
    - ``sigmaArcsec`` -- the trace's Gaussian standard deviation near the edge, in arcsec

    A ``ValueError`` is raised when the order has fewer than ``MIN_FIT_SAMPLES`` accepted samples. An order that can
    anchor has at least that many.
    """
    order = rows["order"].iloc[0]
    accepted = rows.loc[rows[f"cont_{axisA}"].notna() & ~rows["pre-clipped"].astype(bool)]
    if len(accepted.index) < MIN_FIT_SAMPLES:
        raise ValueError(
            f"The anchor order {order} needs at least {MIN_FIT_SAMPLES} accepted samples, got {len(accepted.index)}."
        )
    wavelength = accepted["wavelength"].to_numpy(float)
    zeroPixels, _ = pixelsAt(order, wavelength, 0.0)
    offsetPx = accepted[f"cont_{axisA}"].to_numpy(float) - np.asarray(zeroPixels, dtype=float)
    curve, _ = clipped_polynomial_residuals(wavelength, offsetPx, SCATTER_FIT_DEGREE)
    scale = float(pixels_per_arcsec(pixelsAt, order, np.array([edgeNm]))[0])

    span = rows["wavelength"].max() - rows["wavelength"].min()
    isNearEdge = np.abs(wavelength - edgeNm) <= EDGE_FRACTION * span
    nearEdge = accepted.loc[isNearEdge] if isNearEdge.sum() >= MIN_EDGE_SAMPLES else accepted
    sigmaPx = float(np.median(nearEdge["gauss_stddev"].to_numpy(float)))
    return float(curve(edgeNm) / scale), abs(sigmaPx / scale)


def search_window_arcsec(sigmaArcsec):
    """*the half-width of the search window around the anchor, in arcsec*

    The window is three widths either side of the anchor, but never under 1.5 arcsec and never over 1.8 arcsec, so it
    cannot reach a source 3.7 arcsec away.

    **Key Arguments:**

    - ``sigmaArcsec`` -- the Gaussian standard deviation of the anchor trace, in arcsec

    **Return:**

    - ``windowArcsec`` -- the half-width of the window, in arcsec
    """
    return min(max(WINDOW_WIDTHS * sigmaArcsec, WINDOW_MIN_ARCSEC), WINDOW_MAX_ARCSEC)


@dataclass(frozen=True)
class bin_geometry:
    """Where to look for a trace: order, slit position, width, search window, shared edge and SNR limit."""

    order: float
    slitArcsec: float
    sigmaArcsec: float
    windowArcsec: float
    edgeNm: float
    minSnr: float


def fit_bin(image, pixelsAt, geometry, lowNm, highNm, axisA):
    """*the centre, in pixels from the anchor track, of the trace in one wavelength bin, or None*"""
    order, slitArcsec, sigmaArcsec = geometry.order, geometry.slitArcsec, geometry.sigmaArcsec
    windowArcsec, edgeNm, minSnr = geometry.windowArcsec, geometry.edgeNm, geometry.minSnr
    _, ends = pixelsAt(order, np.array([lowNm, highNm]), slitArcsec)
    pointCount = max(int(abs(ends[1] - ends[0])) + 1, MIN_BIN_TRACK_POINTS)
    wavelength = np.linspace(lowNm, highNm, pointCount)
    middleNm = float((lowNm + highNm) / 2.0)
    scale = abs(float(pixels_per_arcsec(pixelsAt, order, np.array([middleNm]))[0]))
    sigmaPx = sigmaArcsec * scale * (middleNm / edgeNm) ** WIDTH_EXPONENT
    windowPx = windowArcsec * scale
    trackA, trackB = pixelsAt(order, wavelength, slitArcsec)
    halfLength = int(np.ceil(windowPx + FIT_MARGIN_SIGMAS * sigmaPx)) + 1
    offsets, profile = collapse_bin_profile(image, trackA, trackB, halfLength, axisA)
    centre = fit_fixed_width_gaussian(profile, offsets, sigmaPx, windowPx, minSnr)
    return None if centre is None else (float(wavelength.mean()), centre)
