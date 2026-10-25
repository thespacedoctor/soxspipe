"""Anchor the trace of a faint VIS order to its anchor order and refit it on binned profiles.

The slice-by-slice Gaussian fit of ``detect_continuum`` follows a bright order well and a faint order badly: at a
signal-to-noise ratio near 0.5 it fits noise, and a degree-5 polynomial through that noise wanders by arcseconds. The
functions here find the faint orders, take the slit position of the anchor order (the adjacent bright order, or the
better of two faint ones) as an anchor, collapse the faint order into a few wavelength bins and fit one fixed-width
Gaussian per bin inside a narrow window around the anchor. The module imports numpy only, takes its tables as pandas
data-frames, reads the dispersion map through a ``pixelsAt`` function from ``faint_order_map``, and runs no
least-squares optimiser.
"""

import math
from dataclasses import dataclass
from numbers import Real

import numpy as np

from soxspipe.commonutils.faint_order_anchor import (
    MAD_TO_SIGMA,
    MIN_FIT_SAMPLES,
    SCATTER_FIT_DEGREE,
    WIDTH_EXPONENT,
    bin_geometry,
    bright_edge_anchor,
    clipped_polynomial_residuals,
    fit_bin,
    pixels_per_arcsec,
    profile_edge_slit,
    search_window_arcsec,
)

# CODE FALLBACKS FOR SETTINGS FILES WRITTEN BEFORE THE FAINT-ORDER KEYS EXISTED. THE TRACE IS OFF BY DEFAULT, SO OLD
# SETTINGS FILES AND THE -STD RECIPES DO NOT CHANGE
FAINT_TRACE_DEFAULTS = {
    "faint-order-trace": False,
    "faint-order-min-detection-fraction": 0.5,
    "faint-order-max-scatter-px": 1.5,
    "faint-order-max-anchor-scatter-px": 3.0,
    "faint-order-bin-count": 15,
    "faint-order-min-bin-snr": 3.0,
    "same-source-max-join-arcsec": 2.0,
}


def _checked_bool(values, key):
    """*the setting ``key`` as a bool, raising a ``ValueError`` that names it when it is not a real bool*"""
    value = values[key]
    if not isinstance(value, bool):
        raise ValueError(f"The setting `{key}` must be true or false, not {value!r}.")
    return value


def _checked_positive(values, key):
    """*the setting ``key`` as a float, raising a ``ValueError`` that names it unless it is finite and above zero*"""
    value = values[key]
    isNumber = isinstance(value, Real) and not isinstance(value, bool)
    if not isNumber or not math.isfinite(value) or value <= 0:
        raise ValueError(f"The setting `{key}` must be a finite number above zero, not {value!r}.")
    return float(value)


def _checked_fraction(values, key):
    """*the setting ``key`` as a float, raising a ``ValueError`` that names it unless it is in (0, 1]*"""
    value = _checked_positive(values, key)
    if value > 1:
        raise ValueError(f"The setting `{key}` must be a fraction above zero and at most one, not {values[key]!r}.")
    return value


def _checked_positive_integer(values, key):
    """*the setting ``key`` as an int, raising a ``ValueError`` that names it unless it is an integer above zero*"""
    value = values[key]
    if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
        raise ValueError(f"The setting `{key}` must be a whole number above zero, not {value!r}.")
    return value


@dataclass(frozen=True)
class faint_trace_settings:
    """The ``detect-continuum`` settings of the faint-order trace.

    ``isEnabled`` switches the whole feature. ``minDetectionFraction`` and ``maxScatterPx`` flag a faint order.
    ``maxAnchorScatterPx`` is the largest scatter at which a faint order can still anchor its neighbour, provided it
    meets ``minDetectionFraction``. ``binCount`` and ``minBinSnr`` shape the binned fit, and ``maxJoinArcsec`` is the
    same-source warning limit.
    """

    isEnabled: bool
    minDetectionFraction: float
    maxScatterPx: float
    maxAnchorScatterPx: float
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

        A key that is missing or *None* takes its default. A malformed value raises a ``ValueError`` that names the key:
        the switch must be a bool, the bin count a whole number above zero, the detection fraction above zero and at
        most one, and every other number finite and above zero.
        """
        block = block or {}
        values = {
            key: default if block.get(key) is None else block[key] for key, default in FAINT_TRACE_DEFAULTS.items()
        }
        return cls(
            isEnabled=_checked_bool(values, "faint-order-trace"),
            minDetectionFraction=_checked_fraction(values, "faint-order-min-detection-fraction"),
            maxScatterPx=_checked_positive(values, "faint-order-max-scatter-px"),
            maxAnchorScatterPx=_checked_positive(values, "faint-order-max-anchor-scatter-px"),
            binCount=_checked_positive_integer(values, "faint-order-bin-count"),
            minBinSnr=_checked_positive(values, "faint-order-min-bin-snr"),
            maxJoinArcsec=_checked_positive(values, "same-source-max-join-arcsec"),
        )


# THE SOXS VIS BAND GROUPS, EACH FROM THE BLUER ORDER TO THE REDDER: (U, G) AND (R, I). THIS MATCHES THE LIST IN
# `toolkit.SOXS_VIS_ORDERS_BLUE_TO_RED`, WHICH THIS MODULE DOES NOT IMPORT TO STAY LIGHT
VIS_BAND_GROUPS = ((2, 3), (1, 4))
VIS_ORDERS_BLUE_TO_RED = tuple(order for group in VIS_BAND_GROUPS for order in group)
# WHERE AN ANCHOR'S SLIT POSITION CAME FROM: ITS CLIPPED CURVE, ITS EDGE PROFILE, OR THE CURVE AFTER THE PROFILE FIT
# FAILED
ANCHOR_CURVE = "curve"
ANCHOR_PROFILE = "profile"
ANCHOR_FALLBACK = "curve after profile failed"
# THE EDGE PROFILE OF A FAINT ANCHOR COVERS THIS MANY FAINT-ORDER BIN WIDTHS NEXT TO THE SHARED EDGE. ON VV4580SGR THE
# ANCHOR'S FIRST 5.5 NM ARE MASKED AND ONE BIN WIDTH IS 6.2 NM, SO A SINGLE WIDTH HOLDS NO DATA AND TWO FAIL THE
# SIGNAL-TO-NOISE LIMIT. THREE (18.6 NM) PASS, AND THE CENTRE MOVES BY ABOUT 0.05 ARCSEC BETWEEN 15 AND 60 NM
EDGE_PROFILE_BIN_WIDTHS = 3.0
# THE POLYNOMIAL THROUGH THE BIN CENTRES: NO CURVE BELOW THIS MANY KEPT BINS, A STRAIGHT LINE BELOW THE SECOND LIMIT
# AND A PARABOLA FROM IT
MIN_CURVE_BINS = 4
MIN_PARABOLA_BINS = 7


@dataclass(frozen=True)
class order_trace_stats:
    """How well the slice-by-slice fit found one order.

    ``detectionFraction`` is the fraction of the order's samples with a detected continuum. ``scatterPx`` is the MAD
    scatter of the centre shift of the accepted samples about a degree-2 curve, or NaN when too few were accepted.
    ``canAnchor`` is *True* when the order is detected often enough and scatters little enough to anchor a faint
    neighbour.
    """

    detectionFraction: float
    scatterPx: float
    canAnchor: bool = False


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
        detectionFraction, scatterPx = float(isDetected.mean()), _scatter_px(accepted)
        isScattered = np.isnan(scatterPx) or scatterPx > settings.maxScatterPx
        if detectionFraction < settings.minDetectionFraction or isScattered:
            canAnchor = detectionFraction >= settings.minDetectionFraction and scatterPx <= settings.maxAnchorScatterPx
            faintOrders[order] = order_trace_stats(detectionFraction, scatterPx, bool(canAnchor))
    return faintOrders


def _scatter_rank(stats):
    """*the sort key of an order's scatter, with a missing estimate ranked worst*"""
    return np.inf if stats is None or np.isnan(stats.scatterPx) else stats.scatterPx


def _anchor_for_group(group, faintOrders, present):
    """*the (faint order, anchor order) pair of a band group, or None when the group has no transfer*

    When one order is faint, its partner is the anchor. When both are faint, the anchor is the one that ``canAnchor``,
    and the one with the lower scatter when both can. When neither can, the group has no transfer.
    """
    bluer, redder = group
    if bluer not in present or redder not in present:
        return None
    isBluerFaint, isRedderFaint = bluer in faintOrders, redder in faintOrders
    if isBluerFaint != isRedderFaint:
        return (bluer, redder) if isBluerFaint else (redder, bluer)
    if not isBluerFaint:
        return None
    candidates = [order for order in group if faintOrders[order] is not None and faintOrders[order].canAnchor]
    if not candidates:
        return None
    better = min(candidates, key=lambda order: _scatter_rank(faintOrders[order]))
    return (redder, bluer) if better == bluer else (bluer, redder)


def anchor_neighbours(table, faintOrders):
    """*pair each faint order with the order in its band group that anchors it*

    The SOXS VIS orders come in two band groups, u and g, then r and i, and a trace is never carried across the
    boundary between the groups. A faint order whose partner is bright is anchored by that partner. When both orders
    of a group are faint, the best of them that ``canAnchor`` (the one with the lower scatter when both can) anchors
    the other and keeps its own sampled rows. When neither can, or when a partner is not in the table, the group gets
    no transfer.

    **Key Arguments:**

    - ``table`` -- the sampled order table
    - ``faintOrders`` -- the ``order_trace_stats`` of each faint order, keyed by order number

    **Return:**

    - ``anchors`` -- the anchor of each order to transfer, keyed by that order
    """
    present = set(table["order"].unique())
    anchors = {}
    for group in VIS_BAND_GROUPS:
        pair = _anchor_for_group(group, faintOrders, present)
        if pair is not None:
            anchors[float(pair[0])] = float(pair[1])
    return anchors


ANCHOR_NOTES = {
    ANCHOR_CURVE: "",
    ANCHOR_PROFILE: ", anchor from the edge profile",
    ANCHOR_FALLBACK: ", edge profile fit failed so the anchor is the curve estimate",
}


@dataclass(frozen=True)
class transfer_report:
    """What the transfer did for one faint order.

    ``binsKept`` of ``binCount`` bins gave a centre, and ``degree`` is the degree of the curve through them. An order
    that kept no bin is left as it was measured. ``slitArcsec`` is the anchor's slit position and ``anchorSource``
    says how it was found: from the anchor's clipped curve (``"curve"``), from a profile collapsed at the shared edge
    because the anchor is itself faint (``"profile"``), or from the curve after that profile fit failed.
    """

    order: float
    anchorOrder: float
    binsKept: int
    binCount: int
    degree: int
    slitArcsec: float
    anchorSource: str = ANCHOR_CURVE

    @property
    def is_transferred(self):
        """*True when at least one bin was kept, so the order's rows were replaced*"""
        return self.binsKept > 0

    def describe(self, orderLabel, anchorLabel):
        """*the one-line account of the transfer for the log and the extraction notes*

        **Key Arguments:**

        - ``orderLabel`` -- the faint order's label, such as ``"2 (u)"``
        - ``anchorLabel`` -- the anchor order's label, such as ``"3 (g)"``

        **Return:**

        - ``line`` -- the sentence
        """
        suffix = ANCHOR_NOTES[self.anchorSource]
        if not self.is_transferred:
            return f"Order {orderLabel} from anchor {anchorLabel}: no bin kept, trace left as measured{suffix}"
        head = f"Order {orderLabel} trace transferred from {anchorLabel}"
        return f"{head}: {self.binsKept}/{self.binCount} bins, deg {self.degree}{suffix}"


def _shared_edge_nm(brightRows, anchorOrder, order):
    """*the wavelength at the end of the anchor order that faces the faint order*"""
    group = next(group for group in VIS_BAND_GROUPS if order in group)
    isBrightRedder = group.index(anchorOrder) > group.index(order)
    wavelength = brightRows["wavelength"]
    return float(wavelength.min() if isBrightRedder else wavelength.max())


def _bin_edges(wavelength, binCount):
    """*the wavelength edges of ``binCount`` bins that hold equal numbers of samples, covering the whole order*"""
    ordered = np.sort(wavelength)
    groups = np.array_split(ordered, min(binCount, len(ordered)))
    inner = [(left[-1] + right[0]) / 2.0 for left, right in zip(groups[:-1], groups[1:], strict=True)]
    return np.array([ordered[0], *inner, ordered[-1]])


def _offset_curve(binWavelength, binCentre, marginNm):
    """*the curve of the trace's offset from the anchor track against wavelength, and its degree*

    The polynomial is used up to ``marginNm`` (one bin width) beyond the first and last kept bins, which covers the
    ends of the order when its end bins are kept. Beyond that the curve is flat, so it cannot swing away where
    nothing was measured.
    """
    count = len(binWavelength)
    degree = 0 if count < MIN_CURVE_BINS else (1 if count < MIN_PARABOLA_BINS else 2)
    low, high = binWavelength.min() - marginNm, binWavelength.max() + marginNm
    curve = np.polynomial.Polynomial.fit(binWavelength, binCentre, degree)
    return (lambda wavelength: curve(np.clip(wavelength, low, high))), degree


def _edge_profile_slit(table, order, anchorOrder, curveEstimate, *, image, pixelsAt, settings, axisA):
    """*the slit position of a faint anchor from its edge profile, and where it came from*

    The profile is collapsed from the anchor's rows within one faint-order bin width of the shared edge, and a fixed
    width Gaussian is fitted around the curve estimate. When that fit fails the curve estimate stands.
    """
    slitArcsec, sigmaArcsec, edgeNm = curveEstimate
    widthNm = (
        EDGE_PROFILE_BIN_WIDTHS * float(np.ptp(table.loc[table["order"] == order, "wavelength"])) / settings.binCount
    )
    geometry = bin_geometry(
        anchorOrder, slitArcsec, sigmaArcsec, search_window_arcsec(sigmaArcsec), edgeNm, settings.minBinSnr
    )
    measured = profile_edge_slit(
        image,
        pixelsAt,
        geometry,
        anchorWavelengths=table.loc[table["order"] == anchorOrder, "wavelength"].to_numpy(float),
        widthNm=widthNm,
        axisA=axisA,
    )
    return (slitArcsec, ANCHOR_FALLBACK) if measured is None else (measured, ANCHOR_PROFILE)


def _anchor_slit_arcsec(table, order, anchorOrder, *, isAnchorFaint, image, pixelsAt, slitOffsets, settings, axisA):
    """*the anchor of a faint order: slit position, width, shared edge and where the slit position came from*

    A bright anchor takes its slit position from its clipped curve at the shared edge. An anchor that is faint itself
    takes it from a profile collapsed next to the edge, because its curve through scattered samples is not reliable
    there. The slit-illumination offsets of the two orders are then applied.
    """
    brightRows = table.loc[table["order"] == anchorOrder]
    edgeNm = _shared_edge_nm(brightRows, anchorOrder, order)
    slitArcsec, sigmaArcsec = bright_edge_anchor(brightRows, edgeNm, pixelsAt, axisA=axisA)
    source = ANCHOR_CURVE
    if isAnchorFaint:
        slitArcsec, source = _edge_profile_slit(
            table,
            order,
            anchorOrder,
            (slitArcsec, sigmaArcsec, edgeNm),
            image=image,
            pixelsAt=pixelsAt,
            settings=settings,
            axisA=axisA,
        )
    offsets = slitOffsets or {}
    if np.isfinite(offsets.get(order, np.nan)) and np.isfinite(offsets.get(anchorOrder, np.nan)):
        slitArcsec = slitArcsec - offsets[anchorOrder] + offsets[order]
    return slitArcsec, sigmaArcsec, edgeNm, source


def _synthetic_rows(rows, geometry, offsetCurve, pixelsAt, axisA, axisB):
    """*the replacement values of a faint order's rows: the anchor track plus the fitted offset curve*"""
    wavelength = rows["wavelength"].to_numpy(float)
    scale = pixels_per_arcsec(pixelsAt, geometry.order, wavelength)
    slitArcsec = geometry.slitArcsec + offsetCurve(wavelength) / scale
    contA, contB = pixelsAt(geometry.order, wavelength, slitArcsec)
    widthPx = geometry.sigmaArcsec * np.abs(scale) * (wavelength / geometry.edgeNm) ** WIDTH_EXPONENT
    return {
        f"cont_{axisA}": contA,
        f"cont_{axisB}": contB,
        "gauss_stddev": widthPx,
        "gauss_amplitude": np.nan,
        "gauss_mean": np.nan,
        "pre-clipped": False,
        "centre_shift": rows[f"fit_{axisA}"].to_numpy(float) - contA,
        "faint_transfer": True,
    }


def _transfer_one_order(
    table, order, anchorOrder, *, isAnchorFaint, image, pixelsAt, slitOffsets, settings, axisA, axisB
):
    """*the rows, the replacement values (None when no bin was kept) and the report for one faint order*"""
    rows = table.loc[table["order"] == order]
    slitArcsec, sigmaArcsec, edgeNm, source = _anchor_slit_arcsec(
        table,
        order,
        anchorOrder,
        isAnchorFaint=isAnchorFaint,
        image=image,
        pixelsAt=pixelsAt,
        slitOffsets=slitOffsets,
        settings=settings,
        axisA=axisA,
    )
    geometry = bin_geometry(
        order, slitArcsec, sigmaArcsec, search_window_arcsec(sigmaArcsec), edgeNm, settings.minBinSnr
    )
    edges = _bin_edges(rows["wavelength"].to_numpy(float), settings.binCount)
    fits = [
        fit_bin(image, pixelsAt, geometry, low, high, axisA) for low, high in zip(edges[:-1], edges[1:], strict=True)
    ]
    kept = np.array([fit for fit in fits if fit is not None], dtype=float).reshape(-1, 2)
    slit, binCount = float(slitArcsec), len(fits)
    if len(kept) == 0:
        return rows.index, None, transfer_report(order, anchorOrder, 0, binCount, 0, slit, source)
    curve, degree = _offset_curve(kept[:, 0], kept[:, 1], (edges[-1] - edges[0]) / binCount)
    values = _synthetic_rows(rows, geometry, curve, pixelsAt, axisA, axisB)
    return rows.index, values, transfer_report(order, anchorOrder, len(kept), binCount, degree, slit, source)


def transfer_faint_orders(table, image, pixelsAt, *, slitOffsets, settings, axisA, axisB):
    """*replace the rows of each faint order by a trace anchored to its anchor order and fitted on binned profiles*

    Each order that has an anchor in its band group (see ``anchor_neighbours``) is split into wavelength bins. Every
    bin is collapsed along the anchor track into one profile, and a Gaussian of the anchor's width is fitted inside a
    narrow window around the anchor. A low-degree curve through the kept bins gives the trace, and the order's rows
    are replaced one for one with exact samples of that curve, marked in a new ``faint_transfer`` column. An order
    that keeps no bin has no measurement to transfer, so its rows stay as they were measured. When no order is faint,
    none has an anchor, or none keeps a bin, the table is returned as it is.

    **Key Arguments:**

    - ``table`` -- the sampled order table
    - ``image`` -- the 2D frame, NaN where masked
    - ``pixelsAt`` -- the function ``pixelsAt(order, wavelengths, slitArcsec)``, giving the (axisA, axisB) pixels
    - ``slitOffsets`` -- the slit-illumination offset of each order in arcsec, or *None*
    - ``settings`` -- the ``faint_trace_settings``
    - ``axisA`` -- the cross-dispersion axis, ``"x"`` or ``"y"``
    - ``axisB`` -- the dispersion axis, ``"x"`` or ``"y"``

    **Return:**

    - ``table`` -- the same table object when nothing was transferred, otherwise an edited copy
    - ``reports`` -- one ``transfer_report`` per order that has an anchor, including those that kept no bin
    """
    faintOrders = find_faint_orders(table, axisA, settings)
    anchors = anchor_neighbours(table, faintOrders)
    outcomes = [
        _transfer_one_order(
            table,
            order,
            anchorOrder,
            isAnchorFaint=anchorOrder in faintOrders,
            image=image,
            pixelsAt=pixelsAt,
            slitOffsets=slitOffsets,
            settings=settings,
            axisA=axisA,
            axisB=axisB,
        )
        for order, anchorOrder in anchors.items()
    ]
    reports = [report for _, _, report in outcomes]
    if not any(report.is_transferred for report in reports):
        return table, reports
    result = table.copy()
    result["faint_transfer"] = False
    for index, values, _ in outcomes:
        for column, value in (values or {}).items():
            result.loc[index, column] = value
    return result, reports


@dataclass(frozen=True)
class join_mismatch:
    """Two adjacent orders that put the trace at different slit positions where they join.

    ``differenceArcsec`` is the slit position in the redder order minus the one in the bluer order, after the
    slit-illumination offsets are removed.
    """

    blueOrder: float
    redOrder: float
    differenceArcsec: float


def _trace_slit_arcsec(table, order, edge, trace, pixelsAt):
    """*the slit position of the fitted trace at one end of an order, in arcsec from the map's slit-zero track*"""
    wavelength = table.loc[table["order"] == order, "wavelength"]
    edgeNm = np.array([wavelength.min() if edge == "blue" else wavelength.max()])
    mapA, mapB = pixelsAt(order, edgeNm, 0.0)
    scale = pixels_per_arcsec(pixelsAt, order, edgeNm)
    return float(((np.asarray(trace(order, mapB)) - mapA) / scale)[0])


def join_slit_mismatches(table, trace, pixelsAt, slitOffsets, limitArcsec):
    """*find adjacent VIS orders whose fitted traces sit at different slit positions where the orders join*

    The traces of the four orders must follow one source. At each join, from the bluer order's red end to the redder
    order's blue end, the fitted trace is converted to a slit position through the dispersion map. The
    slit-illumination offsets are removed when both orders have one. A difference above ``limitArcsec`` is reported.

    **Key Arguments:**

    - ``table`` -- the fitted order table, with ``order`` and ``wavelength``
    - ``trace`` -- the function ``trace(order, axisBPixels)`` giving the fitted cross-dispersion pixel
    - ``pixelsAt`` -- the function ``pixelsAt(order, wavelengths, slitArcsec)``, giving the (axisA, axisB) pixels
    - ``slitOffsets`` -- the slit-illumination offset of each order in arcsec, or *None*
    - ``limitArcsec`` -- the largest difference that is not reported

    **Return:**

    - ``mismatches`` -- one ``join_mismatch`` for each join over the limit, from blue to red
    """
    present = set(table["order"].unique())
    offsets = slitOffsets or {}
    mismatches = []
    for blue, red in zip(VIS_ORDERS_BLUE_TO_RED[:-1], VIS_ORDERS_BLUE_TO_RED[1:], strict=True):
        if blue not in present or red not in present:
            continue
        difference = _trace_slit_arcsec(table, red, "blue", trace, pixelsAt) - _trace_slit_arcsec(
            table, blue, "red", trace, pixelsAt
        )
        if blue in offsets and red in offsets:
            difference = difference - (offsets[red] - offsets[blue])
        if abs(difference) > limitArcsec:
            mismatches.append(join_mismatch(float(blue), float(red), float(difference)))
    return mismatches
