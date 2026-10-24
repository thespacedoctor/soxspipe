"""Anchor the trace of a faint VIS order to its bright neighbour and refit it on binned profiles.

The slice-by-slice Gaussian fit of ``detect_continuum`` follows a bright order well and a faint order badly: at a
signal-to-noise ratio near 0.5 it fits noise, and a degree-5 polynomial through that noise wanders by arcseconds. The
functions here find the faint orders, take the slit position of the adjacent bright order as an anchor, collapse the
faint order into a few wavelength bins and fit one fixed-width Gaussian per bin inside a narrow window around the
anchor. The module needs numpy and pandas only, and runs no least-squares optimiser.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays

# CODE FALLBACKS FOR SETTINGS FILES WRITTEN BEFORE THE FAINT-ORDER KEYS EXISTED. THE TRACE IS OFF BY DEFAULT, SO OLD
# SETTINGS FILES AND THE -STD RECIPES DO NOT CHANGE
FAINT_TRACE_DEFAULTS = {
    "faint-order-trace": False,
    "faint-order-min-detection-fraction": 0.5,
    "faint-order-max-scatter-px": 1.5,
    "faint-order-bin-count": 15,
    "faint-order-min-bin-snr": 3.0,
    "same-source-max-join-arcsec": 2.0,
}


@dataclass(frozen=True)
class faint_trace_settings:
    """The ``detect-continuum`` settings of the faint-order trace.

    ``isEnabled`` switches the whole feature. ``minDetectionFraction`` and ``maxScatterPx`` flag a faint order,
    ``binCount`` and ``minBinSnr`` shape the binned fit, and ``maxJoinArcsec`` is the same-source warning limit.
    """

    isEnabled: bool
    minDetectionFraction: float
    maxScatterPx: float
    binCount: int
    minBinSnr: float
    maxJoinArcsec: float

    @classmethod
    def from_block(cls, block):
        """*read the faint-order settings from a ``detect-continuum`` block, with code fallbacks*

        **Key Arguments:**

        - ``block`` -- the ``detect-continuum`` settings dictionary, or *None*

        **Return:**

        - ``settings`` -- the frozen settings, off when the block has no ``faint-order-trace`` key
        """
        values = {**FAINT_TRACE_DEFAULTS, **(block or {})}
        return cls(
            isEnabled=bool(values["faint-order-trace"]),
            minDetectionFraction=float(values["faint-order-min-detection-fraction"]),
            maxScatterPx=float(values["faint-order-max-scatter-px"]),
            binCount=int(values["faint-order-bin-count"]),
            minBinSnr=float(values["faint-order-min-bin-snr"]),
            maxJoinArcsec=float(values["same-source-max-join-arcsec"]),
        )


# THE SOXS VIS BAND GROUPS, EACH FROM THE BLUER ORDER TO THE REDDER: (u, g) AND (r, i). THIS MATCHES THE LIST IN
# `toolkit.SOXS_VIS_ORDERS_BLUE_TO_RED`, WHICH THIS MODULE DOES NOT IMPORT TO STAY LIGHT
VIS_BAND_GROUPS = ((2, 3), (1, 4))
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
# CONVERTS A MEDIAN ABSOLUTE DEVIATION TO A STANDARD DEVIATION FOR GAUSSIAN NOISE
MAD_TO_SIGMA = 1.4826


@dataclass(frozen=True)
class order_trace_stats:
    """How well the slice-by-slice fit found one order.

    ``detectionFraction`` is the fraction of the order's samples with a detected continuum. ``scatterPx`` is the MAD
    scatter of the centre shift of the accepted samples about a degree-2 curve, or NaN when too few were accepted.
    """

    detectionFraction: float
    scatterPx: float


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
    """
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


def _scatter_px(accepted):
    """*the MAD scatter of the centre shift of the accepted samples about a clipped degree-2 curve*

    **Key Arguments:**

    - ``accepted`` -- the accepted rows of one order

    **Return:**

    - ``scatterPx`` -- the scatter in pixels, or NaN when fewer than ``MIN_FIT_SAMPLES`` rows were accepted
    """
    if len(accepted.index) < MIN_FIT_SAMPLES:
        return np.nan
    _, residuals = clipped_polynomial_residuals(
        accepted["wavelength"].to_numpy(float), accepted["centre_shift"].to_numpy(float), SCATTER_FIT_DEGREE
    )
    return float(MAD_TO_SIGMA * np.median(np.abs(residuals - np.median(residuals))))


def find_faint_orders(table, axisA, settings):
    """*find the orders whose slice-by-slice trace cannot be trusted*

    An order is faint when the continuum was detected in too small a fraction of its samples, or when the centre
    shifts of its accepted samples scatter about a smooth curve by more than the limit, or when too few samples were
    accepted to tell.

    **Key Arguments:**

    - ``table`` -- the sampled table: ``order``, ``wavelength``, ``cont_<axisA>``, ``centre_shift``, ``pre-clipped``
    - ``axisA`` -- the cross-dispersion axis of the detector, ``"x"`` or ``"y"``
    - ``settings`` -- the ``faint_trace_settings``

    **Return:**

    - ``faintOrders`` -- the ``order_trace_stats`` of each faint order, keyed by order number
    """
    faintOrders = {}
    for order, rows in table.groupby("order", sort=True):
        isDetected = rows[f"cont_{axisA}"].notna()
        accepted = rows.loc[isDetected & ~rows["pre-clipped"].astype(bool)]
        stats = order_trace_stats(float(isDetected.mean()), _scatter_px(accepted))
        isScattered = np.isnan(stats.scatterPx) or stats.scatterPx > settings.maxScatterPx
        if stats.detectionFraction < settings.minDetectionFraction or isScattered:
            faintOrders[order] = stats
    return faintOrders


def anchor_neighbours(table, faintOrders):
    """*pair each faint order with the bright order next to it in its band group*

    The SOXS VIS orders come in two band groups, u and g, then r and i, and a trace is never carried across the
    boundary between the groups. A faint order whose partner is also faint, or is not in the table, gets no anchor.

    **Key Arguments:**

    - ``table`` -- the sampled order table
    - ``faintOrders`` -- the faint orders, keyed by order number

    **Return:**

    - ``anchors`` -- the bright partner of each faint order that has one, keyed by the faint order
    """
    present = set(table["order"].unique())
    anchors = {}
    for group in VIS_BAND_GROUPS:
        for faint, partner in (group, group[::-1]):
            isUsable = faint in faintOrders and partner in present and partner not in faintOrders
            if isUsable:
                anchors[float(faint)] = float(partner)
    return anchors


def make_pixels_at(log, dispersionMap, mapBinning, axisA, axisB):
    """*build the function that turns an order, wavelengths and a slit position into frame pixels*

    The dispersion map gives pixels at the map's binning, so they are divided by ``mapBinning`` to land on the frame.
    Samples off the detector are kept, because a search window can reach past the edge of an order.

    **Key Arguments:**

    - ``log`` -- logger
    - ``dispersionMap`` -- path to the dispersion map
    - ``mapBinning`` -- the frame's binning over the map's binning, ``{"x": ..., "y": ...}``
    - ``axisA`` -- the cross-dispersion axis, ``"x"`` or ``"y"``
    - ``axisB`` -- the dispersion axis, ``"x"`` or ``"y"``

    **Return:**

    - ``pixelsAt`` -- the function ``pixelsAt(order, wavelengths, slitArcsec)``, giving the (axisA, axisB) pixels
    """

    def pixels_at(order, wavelengths, slitArcsec):
        wavelengths = np.asarray(wavelengths, dtype=float)
        table = pd.DataFrame(
            {
                "order": float(order),
                "wavelength": wavelengths,
                "slit_position": np.broadcast_to(np.asarray(slitArcsec, dtype=float), wavelengths.shape),
            }
        )
        table = dispersion_map_to_pixel_arrays(
            log=log, dispersionMapPath=dispersionMap, orderPixelTable=table, removeOffDetectorLocation=False
        )
        return (
            table[f"fit_{axisA}"].to_numpy(float) / mapBinning[axisA],
            table[f"fit_{axisB}"].to_numpy(float) / mapBinning[axisB],
        )

    return pixels_at


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


def bright_edge_anchor(rows, edgeNm, pixelsAt):
    """*the slit position and width of a bright order's trace at the wavelength edge it shares with a faint order*

    The centre shift of the accepted samples is fitted against wavelength with a MAD-clipped degree-2 curve and
    evaluated at ``edgeNm``. The shift is the distance of the map's slit-zero track from the trace, so the slit position
    is minus the shift over the pixel scale. The width is the median Gaussian standard deviation of the accepted
    samples in the quarter of the order's wavelength range nearest the edge.

    **Key Arguments:**

    - ``rows`` -- the sampled rows of one bright order
    - ``edgeNm`` -- the wavelength (nm) of the edge it shares with the faint order
    - ``pixelsAt`` -- the function ``pixelsAt(order, wavelengths, slitArcsec)``, giving the (axisA, axisB) pixels

    **Return:**

    - ``slitArcsec`` -- the trace's slit position at the edge, in arcsec
    - ``sigmaArcsec`` -- the trace's Gaussian standard deviation near the edge, in arcsec
    """
    accepted = rows.loc[rows["centre_shift"].notna() & ~rows["pre-clipped"].astype(bool)]
    wavelength = accepted["wavelength"].to_numpy(float)
    curve, _ = clipped_polynomial_residuals(wavelength, accepted["centre_shift"].to_numpy(float), SCATTER_FIT_DEGREE)
    scale = float(pixels_per_arcsec(pixelsAt, rows["order"].iloc[0], np.array([edgeNm]))[0])

    span = rows["wavelength"].max() - rows["wavelength"].min()
    isNearEdge = np.abs(wavelength - edgeNm) <= EDGE_FRACTION * span
    nearEdge = accepted.loc[isNearEdge] if isNearEdge.sum() >= MIN_EDGE_SAMPLES else accepted
    sigmaPx = float(np.median(nearEdge["gauss_stddev"].to_numpy(float)))
    return float(-curve(edgeNm) / scale), abs(sigmaPx / scale)
