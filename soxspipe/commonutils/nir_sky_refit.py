"""One global sky-line fit across the SOXS NIR orders, in place of a separate median shift for each order (DY-1968).

The NIR orders move together on the detector, as one rigid shift that changes smoothly with the order. Measured one
order at a time, an order with few lines, or a line aliased onto its neighbour, gets a shift that its neighbours
contradict. This module fits every matched line at once:

    dy(order, t) = sum_i a_i m^i + sum_j b_j t^j + c_o + s_o t

``dy`` is the detector motion along the dispersion axis in px. ``m`` is the order normalised to [-1, 1] over the
orders of the frame, and ``t`` the position along the order, from -1 at its bluest wavelength to +1 at its reddest. The
``a`` and ``b`` terms are global. ``c_o`` is a per-order offset: the arc solution of an order can sit several px
from the detector shift of its neighbours (orders 10 and 11), so each order with enough lines keeps the part of its
offset that is statistically significant. ``s_o`` is a per-order slope, fitted only to an order with enough lines that
span enough of it, and shrunk towards zero by its line count, so that an order with a shape error in its arc solution
(order 10) is corrected where the data support it.

Every function here is pure: it reads arrays and returns new objects, and holds no recipe state. The wiring that
detects and matches the sky lines lives in ``horne_extraction``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import NamedTuple

import numpy as np
import pandas as pd

from soxspipe.commonutils.polynomials import solve_linear_poly_coefficients
from soxspipe.commonutils.ug_wavelength_transfer import detector_px_to_nm_shift, nm_shift_to_detector_px

# THE SKY-TUNING SETTINGS KEYS AND THEIR DEFAULTS, USED WHEN A (OLDER) WORKSPACE SETTINGS FILE HAS NO BLOCK
DEFAULT_SETTINGS = {
    "global-fit": True,
    "order-deg": 1,
    "along-order-deg": 1,
    "slope-min-lines": 8,
    "clip-sigma": 3.0,
}
MAX_POLY_DEGREE = 3
# A STRAIGHT LINE THROUGH TWO POINTS FITS EXACTLY, SO A PER-ORDER SLOPE NEEDS AT LEAST THREE LINES
MIN_SLOPE_LINES = 3

# FEWER MATCHED LINES THAN THIS CANNOT PIN DOWN THE GLOBAL MODEL, SO TUNING FALLS BACK TO ONE SHIFT PER ORDER
MIN_GLOBAL_LINES = 10
# THE LINES KEPT AFTER CLIPPING MUST FIT THE MODEL TO THIS RMS (px). A WELL-MATCHED FRAME FITS TO 0.1-0.3 px, WHILE
# RANDOM MATCHES WITHIN THE TIGHT TOLERANCE (UP TO 1.7 px) SCATTER BY ABOUT 1 px EVEN AFTER CLIPPING
MAX_GLOBAL_RMS_PX = 0.6
# AN ORDER'S LINES MUST SPAN THIS FRACTION OF THE ORDER FOR ITS OWN SLOPE TO BE TRUSTED
SLOPE_MIN_SPAN = 0.5
# THE PER-ORDER SLOPE IS MULTIPLIED BY n / (n + SLOPE_SHRINK_LINES) FOR n LINES: HALF OF IT AT 4 LINES, 0.83 AT 20
SLOPE_SHRINK_LINES = 4
# THE MATCH TOLERANCE IS THIS MANY TIMES THE ARC-SOLUTION RMS, BUT NEVER MORE THAN MAX_TIGHT_TOLERANCE_NM. A WIDER
# TOLERANCE ALIASED ORDER 15 ONTO A NEIGHBOURING LINE
TIGHT_TOLERANCE_ARC_RMS = 2.5
MAX_TIGHT_TOLERANCE_NM = 0.1

MAX_CLIP_ITERATIONS = 10
# THE GLOBAL MODEL AND THE PER-ORDER OFFSETS ARE FITTED IN TURN THIS MANY TIMES
OFFSET_ITERATIONS = 3
# AN ORDER NEEDS AT LEAST THIS MANY LINES FOR AN OFFSET OF ITS OWN. FEWER TAKE THE GLOBAL VALUE
MIN_OFFSET_LINES = 4
# AN ORDER'S OWN FIRST-PASS SHIFT IS A USABLE SEED ONLY IF IT KEPT AT LEAST THIS MANY LINES
MIN_OWN_SEED_LINES = 2
# AN ORDER KEEPS z^2 / (z^2 + OFFSET_SHRINK_Z^2) OF ITS MEDIAN OFFSET, z BEING THE OFFSET IN STANDARD ERRORS OF THE
# MEDIAN: HALF AT 3 SIGMA, 0.92 AT 10, AND A NOISE-LEVEL OFFSET (z ~ 1) KEEPS ABOUT A TENTH
OFFSET_SHRINK_Z = 3.0
# AN ORDER IS REPORTED AS TAKING ITS OWN OFFSET (``source``) WHEN THAT OFFSET IS AT LEAST THIS MANY px
SOURCE_OFFSET_PX = 0.1
# STANDARD ERROR OF A MEDIAN, IN UNITS OF sigma / sqrt(n), FOR NORMAL SCATTER
MEDIAN_SE_FACTOR = 1.2533
# THE FIRST MATCH USES THE RIGID SEED, THE LATER ONES THE PREVIOUS FIT
MATCH_ITERATIONS = 3
# THE FLOOR OF THE ROBUST SCATTER (px), SO A NOISE-FREE FIT DOES NOT CLIP ITS OWN ROUNDING ERROR
MIN_CLIP_SCALE_PX = 0.02
# MEDIAN ABSOLUTE DEVIATION TO STANDARD DEVIATION, FOR NORMAL SCATTER
MAD_TO_SIGMA = 1.4826
# SINGULAR VALUES BELOW THIS FRACTION OF THE LARGEST (OF THE COLUMN-SCALED DESIGN MATRIX) MAKE IT RANK-DEFICIENT
RANK_RCOND = 1.0e-9


@dataclass(frozen=True)
class sky_tuning_settings:
    """*the ``sky-tuning`` settings of a SOXS NIR extraction recipe*

    **Key Arguments:**

    - ``globalFit`` -- fit all orders together. *False* restores one median shift per order
    - ``orderDeg`` -- degree of the polynomial in the normalised order
    - ``alongOrderDeg`` -- degree of the polynomial in the position along the order
    - ``slopeMinLines`` -- the fewest lines an order needs for its own slope
    - ``clipSigma`` -- the robust sigma clipping limit
    """

    globalFit: bool = True
    orderDeg: int = 1
    alongOrderDeg: int = 1
    slopeMinLines: int = 8
    clipSigma: float = 3.0


@dataclass(frozen=True)
class order_geometry:
    """*where an order sits in wavelength and how its pixels map to wavelength*

    **Key Arguments:**

    - ``centreNm`` -- the middle of the order's wavelength range (nm), where ``t`` is 0
    - ``halfRangeNm`` -- half the order's wavelength range (nm), so ``t`` is -1 and +1 at its ends
    - ``dispersionNmPerPx`` -- the signed median dispersion d(wavelength)/d(pixel) of the order (nm/px)
    - ``medianPixelScaleNm`` -- the median ``pixelScaleNm`` of the order's rows, which the local scale is set against
    """

    centreNm: float
    halfRangeNm: float
    dispersionNmPerPx: float
    medianPixelScaleNm: float


@dataclass(frozen=True)
class order_sky_correction:
    """*the wavelength correction of one order: a shift at its centre and a polynomial along it*

    **Key Arguments:**

    - ``order`` -- the order
    - ``shiftPx`` -- the detector motion (px) at the order's centre
    - ``alongPx`` -- the coefficients (px) of ``t``, ``t^2``, ... along the order, the first with ``perOrderSlopePx``
      added
    - ``perOrderSlopePx`` -- the shrunk slope (px per unit ``t``) fitted to this order alone, 0 when unsupported
    - ``nLines`` -- the lines of this order kept by the fit
    - ``geometry`` -- the order's ``order_geometry``
    - ``perOrderOffsetPx`` -- the part of ``shiftPx`` that is this order's own offset from the global model, 0 when
      unsupported. Default 0
    - ``source`` -- where ``shiftPx`` came from: ``global`` (the model), ``own-offset`` (the model plus an offset of at
      least ``SOURCE_OFFSET_PX`` from the order's own lines) or ``own-seed`` (the order's own measured shift, kept
      because it has too few lines for an offset, or because no line matched about either seed). Default ``global``
    - ``globalSeedLines`` -- the lines the order matched about the global seed within ``MAX_TIGHT_TOLERANCE_NM`` in the
      first round, or *None* if it supplied none. Default *None*
    - ``ownSeedLines`` -- the same about the order's own seed, or *None* if it had no usable own seed. Default *None*
    """

    order: int
    shiftPx: float
    alongPx: tuple[float, ...]
    perOrderSlopePx: float
    nLines: int
    geometry: order_geometry
    perOrderOffsetPx: float = 0.0
    source: str = "global"
    globalSeedLines: int | None = None
    ownSeedLines: int | None = None

    @property
    def shift_nm(self) -> float:
        """The wavelength shift (nm) at the centre of the order, the value recorded as ``SKY SHIFT O<n>``."""
        return detector_px_to_nm_shift(self.shiftPx, self.geometry.dispersionNmPerPx)

    @property
    def slope_px(self) -> float:
        """The linear coefficient (px per unit ``t``) of the correction along the order."""
        return self.alongPx[0] if self.alongPx else 0.0

    def correction_nm(self, wavelengthNm, pixelScaleNm=None) -> np.ndarray:
        """*the wavelength shift (nm) to add to each row of the order*

          The detector motion at the row's position along the order, times the row's local signed dispersion. The local
          dispersion is the order's median dispersion scaled by the row's ``pixelScaleNm`` over the median.

          **Key Arguments:**

          - ``wavelengthNm`` -- the rows' wavelengths (nm), before the correction
          - ``pixelScaleNm`` -- the rows' ``pixelScaleNm``. Default *None* uses the order's median for every row, and so
        does a row whose scale is not finite and positive

          **Return:**

          - ``shiftNm`` -- the shift (nm) of each row
        """
        geometry = self.geometry
        wavelengths = np.asarray(wavelengthNm, dtype=float)
        t = (wavelengths - geometry.centreNm) / geometry.halfRangeNm
        shiftPx = self.shiftPx + sum(c * t ** (j + 1) for j, c in enumerate(self.alongPx))
        if pixelScaleNm is None:
            dispersion = geometry.dispersionNmPerPx
        else:
            scales = np.asarray(pixelScaleNm, dtype=float)
            # A ROW WITH NO USABLE PIXEL SCALE TAKES THE ORDER'S MEDIAN, SO IT KEEPS A FINITE WAVELENGTH
            scale = np.where(np.isfinite(scales) & (scales > 0), scales, geometry.medianPixelScaleNm)
            scale = scale / geometry.medianPixelScaleNm
            dispersion = geometry.dispersionNmPerPx * scale
        return -dispersion * shiftPx


@dataclass(frozen=True)
class sky_refit_result:
    """*the outcome of ``fit_nir_sky_model``*

    **Key Arguments:**

    - ``corrections`` -- ``{order: order_sky_correction}`` for every order with a geometry, empty on a fallback
    - ``coefficients`` -- the global coefficients (px): ``m^0 .. m^orderDeg`` then ``t^1 .. t^alongOrderDeg``
    - ``nLines`` -- the lines kept by the fit (the lines offered, on a fallback)
    - ``rmsPx`` -- the rms (px) of the kept lines about the full model, NaN on a fallback before the fit
    - ``clipFraction`` -- the fraction of lines clipped, NaN when there was no fit
    - ``lines`` -- one row per line offered: ``order``, ``catalogueNm``, ``observedNm``, and once fitted ``residualPx``
      (the offset before the model), ``modelPx`` (the model at the line) and ``clipped``
    - ``fallbackReason`` -- why the fit cannot be used, or *None*
    """

    corrections: dict[int, order_sky_correction]
    coefficients: tuple[float, ...]
    nLines: int
    rmsPx: float
    clipFraction: float
    lines: pd.DataFrame
    fallbackReason: str | None


def _integer_setting(block, key, lowest, highest=None):
    """*the integer setting ``key`` of the block, or raise ``ValueError`` if it is not an integer in range*"""
    value = block.get(key, DEFAULT_SETTINGS[key])
    isInteger = isinstance(value, int) and not isinstance(value, bool)
    if not isInteger or value < lowest or (highest is not None and value > highest):
        bounds = f">= {lowest}" if highest is None else f"between {lowest} and {highest}"
        raise ValueError(f"the sky-tuning setting {key} must be a whole number {bounds}, got {value!r}")
    return value


def read_sky_tuning_settings(block) -> sky_tuning_settings:
    """*read the ``sky-tuning`` settings block of a recipe, with defaults for a missing block or key*

    **Key Arguments:**

    - ``block`` -- the ``sky-tuning`` mapping, or *None* when the settings file has none (an older workspace)

    **Return:**

    - ``settings`` -- the validated ``sky_tuning_settings``

    **Raises:**

    - ``ValueError`` -- on an unknown key, or a value that cannot drive the fit
    """
    block = {} if block is None else block
    if not isinstance(block, Mapping):
        raise ValueError(f"the sky-tuning settings must be a mapping of keys to values, got {block!r}")
    unknown = sorted(set(block) - set(DEFAULT_SETTINGS))
    if unknown:
        raise ValueError(f"the sky-tuning settings have no key {unknown[0]!r} (known keys: {sorted(DEFAULT_SETTINGS)})")
    globalFit = block.get("global-fit", DEFAULT_SETTINGS["global-fit"])
    if not isinstance(globalFit, bool):
        raise ValueError(f"the sky-tuning setting global-fit must be True or False, got {globalFit!r}")
    clipSigma = block.get("clip-sigma", DEFAULT_SETTINGS["clip-sigma"])
    isNumber = isinstance(clipSigma, (int, float)) and not isinstance(clipSigma, bool)
    if not (isNumber and math.isfinite(clipSigma) and clipSigma > 0):
        raise ValueError(f"the sky-tuning setting clip-sigma must be a finite number > 0, got {clipSigma!r}")
    return sky_tuning_settings(
        globalFit=globalFit,
        orderDeg=_integer_setting(block, "order-deg", 0, MAX_POLY_DEGREE),
        alongOrderDeg=_integer_setting(block, "along-order-deg", 0, MAX_POLY_DEGREE),
        slopeMinLines=_integer_setting(block, "slope-min-lines", MIN_SLOPE_LINES),
        clipSigma=float(clipSigma),
    )


def geometry_from_rows(wavelengths, pixelScales, dispersionNmPerPx) -> order_geometry:
    """*the ``order_geometry`` of an order from its extracted rows*

    **Key Arguments:**

    - ``wavelengths`` -- the order's ``wavelengthMean`` values (nm)
    - ``pixelScales`` -- the order's ``pixelScaleNm`` values (nm)
    - ``dispersionNmPerPx`` -- the signed median dispersion of the order (nm/px)

    **Return:**

    - ``geometry`` -- the order's ``order_geometry``

    **Raises:**

    - ``ValueError`` -- if the order has no finite wavelength range, no positive pixel scale or no usable dispersion
    """
    values = pd.to_numeric(pd.Series(wavelengths), errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or values.max() == values.min():
        raise ValueError("the order has no wavelength range to place sky lines in")
    scales = pd.to_numeric(pd.Series(pixelScales), errors="coerce").to_numpy(dtype=float)
    scales = scales[np.isfinite(scales)]
    geometry = order_geometry(
        centreNm=float((values.max() + values.min()) / 2),
        halfRangeNm=float((values.max() - values.min()) / 2),
        dispersionNmPerPx=float(dispersionNmPerPx),
        medianPixelScaleNm=float(np.median(scales)) if scales.size else float("nan"),
    )
    if not _is_usable_geometry(geometry):
        raise ValueError("the order has no positive pixel scale or no usable dispersion to place sky lines with")
    return geometry


def _is_usable_geometry(geometry) -> bool:
    """*whether every number of an ``order_geometry`` is finite, with a non-zero dispersion and positive scales*"""
    values = (geometry.centreNm, geometry.halfRangeNm, geometry.dispersionNmPerPx, geometry.medianPixelScaleNm)
    return bool(
        all(math.isfinite(v) for v in values)
        and geometry.halfRangeNm > 0
        and geometry.dispersionNmPerPx != 0
        and geometry.medianPixelScaleNm > 0
    )


def match_tolerance_nm(arcRmsPx, dispersionNmPerPx) -> float:
    """*the tight catalogue match tolerance (nm) of an order*

    **Key Arguments:**

    - ``arcRmsPx`` -- the rms (px) of the arc-solution fit along the dispersion, or *None* / NaN when the arc has none
    - ``dispersionNmPerPx`` -- the order's signed dispersion (nm/px)

    **Return:**

    - ``toleranceNm`` -- ``TIGHT_TOLERANCE_ARC_RMS`` times the arc rms in nm, at most ``MAX_TIGHT_TOLERANCE_NM``
    """
    if arcRmsPx is None or not math.isfinite(arcRmsPx):
        return MAX_TIGHT_TOLERANCE_NM
    return min(TIGHT_TOLERANCE_ARC_RMS * arcRmsPx * abs(dispersionNmPerPx), MAX_TIGHT_TOLERANCE_NM)


def match_peaks_to_catalogue(predictedPeakNm, catalogueNm, toleranceNm):
    """*pair each observed sky-line centre with the catalogue line that is its nearest, and it theirs*

    A pair must lie within the tolerance and each be the other's nearest neighbour, so one peak is never matched to two
    catalogue lines (or the reverse) and a line that is not the nearest one is not matched at all.

    **Key Arguments:**

    - ``predictedPeakNm`` -- the observed line centres (nm), moved by the current model
    - ``catalogueNm`` -- the catalogue wavelengths (nm)
    - ``toleranceNm`` -- the largest separation (nm) of a pair

    **Return:**

    - ``peakIndex`` -- the matched peaks, ascending
    - ``catalogueIndex`` -- the catalogue line of each matched peak
    """
    peaks = np.asarray(predictedPeakNm, dtype=float)
    catalogue = np.asarray(catalogueNm, dtype=float)
    if peaks.size == 0 or catalogue.size == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    separation = np.abs(peaks[:, None] - catalogue[None, :])
    nearestCatalogue = separation.argmin(axis=1)
    nearestPeak = separation.argmin(axis=0)
    peakIndex = np.arange(peaks.size)
    isMutual = nearestPeak[nearestCatalogue] == peakIndex
    isClose = separation[peakIndex, nearestCatalogue] < toleranceNm
    keep = isMutual & isClose
    return peakIndex[keep], nearestCatalogue[keep]


def _weighted_median(values, weights):
    """*the smallest value whose cumulative weight reaches half the total*"""
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(np.asarray(weights, dtype=float)[order])
    return float(np.asarray(values, dtype=float)[order][np.searchsorted(cumulative, cumulative[-1] / 2)])


def global_seed_px(shiftsNm, lineCounts, dispersions):
    """*the rigid detector shift (px) that seeds the global line matching*

    The weighted median, by line count, of the first-pass shifts of the orders that hold at least the median number of
    lines. An order with few lines, or one aliased onto a neighbouring line, cannot move it.

    **Key Arguments:**

    - ``shiftsNm`` -- ``{order: first-pass wavelength shift (nm)}``
    - ``lineCounts`` -- ``{order: sky lines kept in that first pass}``
    - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``

    **Return:**

    - ``seedPx`` -- the seed (px), or *None* if no order with a dispersion kept a line
    """
    seeds = own_seeds_px(shiftsNm, lineCounts, dispersions)
    if not seeds:
        return None
    threshold = float(np.median([lineCounts[o] for o in seeds]))
    rich = [o for o in seeds if lineCounts[o] >= threshold]
    return _weighted_median([seeds[o] for o in rich], [lineCounts[o] for o in rich])


def own_seeds_px(shiftsNm, lineCounts, dispersions):
    """*the first-pass shift (px) of each order that measured one, the seed that order found for itself*

    An order has a usable seed only if its first pass kept at least ``MIN_OWN_SEED_LINES`` lines, it has a dispersion,
    and its shift is finite and not 0.0, which is what the first pass leaves for an order it did not measure.

    **Key Arguments:**

    - ``shiftsNm`` -- ``{order: first-pass wavelength shift (nm)}``
    - ``lineCounts`` -- ``{order: sky lines kept in that first pass}``
    - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``

    **Return:**

    - ``seeds`` -- ``{order: shift (px)}`` for each order with a usable seed
    """
    return {
        int(o): nm_shift_to_detector_px(shiftsNm[o], dispersions[o])
        for o in shiftsNm
        if lineCounts.get(o, 0) >= MIN_OWN_SEED_LINES
        and o in dispersions
        and math.isfinite(shiftsNm[o])
        and shiftsNm[o] != 0.0
    }


def _is_rank_deficient(design) -> bool:
    """*whether the column-scaled design matrix has a column that the others explain*"""
    scale = np.linalg.norm(design, axis=0)
    scale[scale == 0] = 1.0
    singular = np.linalg.svd(design / scale, compute_uv=False)
    return bool(design.shape[0] < design.shape[1] or singular.min() < RANK_RCOND * singular.max())


def _solve(design, y):
    """*the least-squares coefficients of the design matrix, or None if it is rank-deficient*"""
    if _is_rank_deficient(design):
        return None
    return solve_linear_poly_coefficients(
        lambda matrix, *coefficients: matrix @ np.asarray(coefficients), design, y, design.shape[1]
    )


def _is_within_clip(design, y, coefficients, kept, clipSigma):
    """*the lines within ``clipSigma`` robust sigmas of the model, judged against the scatter of the kept lines*"""
    residual = y - design @ coefficients
    keptResidual = residual[kept]
    mad = np.median(np.abs(keptResidual - np.median(keptResidual)))
    return np.abs(residual) <= clipSigma * max(MAD_TO_SIGMA * mad, MIN_CLIP_SCALE_PX)


def _robust_fit(design, y, clipSigma):
    """*fit by least squares with iterative clipping; return the coefficients and the kept lines, or None*"""
    kept = np.ones(len(y), dtype=bool)
    coefficients = _solve(design, y)
    if coefficients is None:
        return None
    for _ in range(MAX_CLIP_ITERATIONS):
        newKept = _is_within_clip(design, y, coefficients, kept, clipSigma)
        if np.array_equal(newKept, kept) or newKept.sum() <= design.shape[1]:
            break
        newCoefficients = _solve(design[newKept], y[newKept])
        if newCoefficients is None:
            break
        kept, coefficients = newKept, newCoefficients
    return coefficients, kept


def _global_design(m, t, settings):
    """*the design matrix of the global model: ``m^0 .. m^orderDeg`` then ``t^1 .. t^alongOrderDeg``*"""
    columns = [m**i for i in range(settings.orderDeg + 1)] + [t**j for j in range(1, settings.alongOrderDeg + 1)]
    return np.column_stack(columns)


def _line_residuals_px(lines, geometry):
    """*the offset of each line from its catalogue wavelength as detector motion (px), on the order's local scale*"""
    residual = np.empty(len(lines))
    orders = lines["order"].to_numpy()
    shiftNm = (lines["catalogueNm"] - lines["observedNm"]).to_numpy(dtype=float)
    for o in np.unique(orders):
        rows = orders == o
        residual[rows] = nm_shift_to_detector_px(shiftNm[rows], geometry[o].dispersionNmPerPx)
        if "pixelScaleNm" in lines.columns:
            localScale = lines.loc[rows, "pixelScaleNm"].to_numpy(dtype=float)
            residual[rows] *= geometry[o].medianPixelScaleNm / localScale
    return residual


def _normalised_orders(orders):
    """*the function mapping an order number to [-1, 1] over the orders of the frame*"""
    low, high = min(orders), max(orders)
    half = (high - low) / 2 or 1.0
    return lambda order: (np.asarray(order, dtype=float) - (high + low) / 2) / half


def _order_slope(t, residual, settings):
    """*the shrunk slope of one order's residuals against ``t`` and the lines it kept, or a zero slope*"""
    none = (0.0, np.zeros(len(t), dtype=bool))
    if len(t) < settings.slopeMinLines:
        return none
    fit = _robust_fit(np.column_stack([np.ones(len(t)), t]), residual, settings.clipSigma)
    if fit is None:
        return none
    coefficients, kept = fit
    nKept = int(kept.sum())
    span = np.ptp(t[kept]) if nKept else 0.0
    if nKept < settings.slopeMinLines or span < 2 * SLOPE_MIN_SPAN:
        return none
    return float(coefficients[1] * nKept / (nKept + SLOPE_SHRINK_LINES)), kept


def _fit_order_slopes(lines, t, globalResidual, kept, settings):
    """*each supported order's shrunk slope, and the lines kept once those orders are re-judged by their own fit*

    An order's own slope is fitted to the residuals of all its lines, including ones the global clip dropped, so an
    order with a shape error is not clipped out of the one fit that could correct it.
    """
    orders = lines["order"].to_numpy()
    slopes = {}
    newKept = kept.copy()
    for o in np.unique(orders):
        rows = np.flatnonzero(orders == o)
        slope, orderKept = _order_slope(t[rows], globalResidual[rows], settings)
        if slope != 0.0:
            slopes[o] = slope
            newKept[rows] = orderKept
    return slopes, newKept


def _order_offsets(lineOrders, residual):
    """*the part of each order's median offset from the global model that is statistically significant (px)*

    An order keeps ``z^2 / (z^2 + OFFSET_SHRINK_Z^2)`` of its median residual, with ``z`` the median in standard errors.
    The scatter is that of the lines about their own order's median, pooled over the orders, so a real offset between
    orders does not inflate it. An order with fewer than ``MIN_OFFSET_LINES`` lines has none.
    """
    orders = [o for o in np.unique(lineOrders) if (lineOrders == o).sum() >= MIN_OFFSET_LINES]
    if not orders:
        return {}
    medians = {o: float(np.median(residual[lineOrders == o])) for o in orders}
    within = np.concatenate([residual[lineOrders == o] - medians[o] for o in orders])
    sigma = max(MAD_TO_SIGMA * float(np.median(np.abs(within))), MIN_CLIP_SCALE_PX)
    offsets = {}
    for o in orders:
        standardError = sigma * MEDIAN_SE_FACTOR / math.sqrt((lineOrders == o).sum())
        z2 = (medians[o] / standardError) ** 2
        offsets[o] = medians[o] * z2 / (z2 + OFFSET_SHRINK_Z**2)
    return offsets


def _fit_with_offsets(design, y, lineOrders, clipSigma):
    """*fit the global model and the per-order offsets in turn; return the coefficients, kept lines and offsets*

    The first fit has no offsets and clips an order that sits far from the others. Each round then takes the offsets
    from the residuals of the last global fit and refits the global model to the lines with their offsets removed, so
    one order with a genuine offset neither biases the others nor is clipped away.

    **Return:**

    - ``(coefficients, kept, offsets)`` -- or *None* if the design matrix is rank-deficient
    """
    fit = _robust_fit(design, y, clipSigma)
    offsets = {}
    for _ in range(OFFSET_ITERATIONS):
        if fit is None:
            return None
        offsets = _order_offsets(lineOrders, y - design @ fit[0])
        offsetPx = np.array([offsets.get(o, 0.0) for o in lineOrders])
        fit = _robust_fit(design, y - offsetPx, clipSigma)
    if fit is None:
        return None
    return fit[0], fit[1], offsets


def _build_corrections(lines, coefficients, slopes, offsets, kept, geometry, normalise, settings):
    """*the ``order_sky_correction`` of every order in ``geometry``*"""
    nOrderTerms = settings.orderDeg + 1
    counts = lines.loc[kept, "order"].value_counts()
    corrections = {}
    for o, orderGeometry in geometry.items():
        m = float(normalise(o))
        offset = float(offsets.get(o, 0.0))
        shiftPx = float(sum(coefficients[i] * m**i for i in range(nOrderTerms))) + offset
        along = [float(c) for c in coefficients[nOrderTerms:]] or [0.0]
        slope = float(slopes.get(o, 0.0))
        along[0] += slope
        corrections[int(o)] = order_sky_correction(
            order=int(o),
            shiftPx=shiftPx,
            alongPx=tuple(along),
            perOrderSlopePx=slope,
            nLines=int(counts.get(o, 0)),
            geometry=orderGeometry,
            perOrderOffsetPx=offset,
            source="own-offset" if abs(offset) >= SOURCE_OFFSET_PX else "global",
        )
    return corrections


def _fallback(lines, reason, rmsPx=float("nan"), clipFraction=float("nan")):
    """*a result with no corrections, naming why the global fit cannot be used*"""
    return sky_refit_result({}, (), len(lines), rmsPx, clipFraction, lines.copy(), reason)


def _usable_lines(lines, geometry):
    """*the lines whose values are finite, with a positive pixel scale, and whose order has a usable geometry*"""
    isUsable = np.isfinite(lines["catalogueNm"].to_numpy(dtype=float)) & np.isfinite(
        lines["observedNm"].to_numpy(dtype=float)
    )
    if "pixelScaleNm" in lines.columns:
        scale = lines["pixelScaleNm"].to_numpy(dtype=float)
        isUsable &= np.isfinite(scale) & (scale > 0)
    isUsable &= lines["order"].isin(list(geometry)).to_numpy()
    return lines.loc[isUsable].reset_index(drop=True)


def fit_nir_sky_model(lines, geometry, settings) -> sky_refit_result:
    """*fit one detector-motion model to the matched sky lines of every order*

    The offsets are fitted in px, so orders with opposite dispersion signs share one model. Lines more than
    ``clipSigma`` robust sigmas from the model are clipped. An order with enough lines gets the significant part of its
    own offset from the model (``_order_offsets``), and then, if its lines span enough of it, a shrunk slope along the
    order. Every order with a usable geometry gets a correction, including an order with no lines, which takes the
    model's value.

    **Key Arguments:**

    - ``lines`` -- one row per matched line: ``order``, ``catalogueNm``, ``observedNm`` (the catalogue and the
      observed wavelength, nm) and optionally ``pixelScaleNm`` at the line. Rows with a non-finite value or a
      non-positive pixel scale are dropped, and so are the rows of an order whose geometry is not usable
    - ``geometry`` -- ``{order: order_geometry}`` of every order to correct
    - ``settings`` -- the ``sky_tuning_settings``

    **Return:**

    - ``result`` -- the ``sky_refit_result``. Its ``fallbackReason`` is set, and its corrections are empty, when there
      are too few lines, the model is rank-deficient or the lines scatter more than ``MAX_GLOBAL_RMS_PX``

    **Raises:**

    - ``ValueError`` -- if a line belongs to an order that has no geometry
    """
    lines = lines.reset_index(drop=True)
    unknown = sorted(set(lines["order"]) - set(geometry))
    if unknown:
        raise ValueError(f"order {unknown[0]} has sky lines but no geometry")
    geometry = {o: g for o, g in geometry.items() if _is_usable_geometry(g)}
    lines = _usable_lines(lines, geometry)
    if len(lines) < MIN_GLOBAL_LINES:
        return _fallback(lines, f"{len(lines)} sky lines matched (min {MIN_GLOBAL_LINES})")

    normalise = _normalised_orders(list(geometry))
    orderCentre = np.array([geometry[o].centreNm for o in lines["order"]])
    orderHalf = np.array([geometry[o].halfRangeNm for o in lines["order"]])
    t = (lines["observedNm"].to_numpy(dtype=float) - orderCentre) / orderHalf
    lineOrders = lines["order"].to_numpy()
    design = _global_design(normalise(lineOrders), t, settings)
    residual = _line_residuals_px(lines, geometry)

    fit = _fit_with_offsets(design, residual, lineOrders, settings.clipSigma)
    if fit is None:
        return _fallback(lines, "the global model is rank-deficient (too few distinct orders or positions)")
    coefficients, kept, offsets = fit
    offsetPx = np.array([offsets.get(o, 0.0) for o in lineOrders])
    slopes, kept = _fit_order_slopes(lines, t, residual - design @ coefficients - offsetPx, kept, settings)

    modelPx = design @ coefficients + offsetPx + np.array([slopes.get(o, 0.0) for o in lineOrders]) * t
    rmsPx = float(np.sqrt(np.mean((residual - modelPx)[kept] ** 2)))
    clipFraction = float(1.0 - kept.mean())
    if rmsPx > MAX_GLOBAL_RMS_PX:
        return _fallback(lines, f"fit rms {rmsPx:0.2f} px (max {MAX_GLOBAL_RMS_PX})", rmsPx, clipFraction)

    lines = lines.assign(residualPx=residual, modelPx=modelPx, clipped=~kept)
    corrections = _build_corrections(lines, coefficients, slopes, offsets, kept, geometry, normalise, settings)
    return sky_refit_result(
        corrections, tuple(float(c) for c in coefficients), int(kept.sum()), rmsPx, clipFraction, lines, None
    )


def apply_sky_corrections(extractedOrdersDF, corrections) -> pd.DataFrame:
    """*add each order's row-by-row correction to its ``wavelengthMean``, in a new data-frame*

    Orders without a correction are left as they are.

    **Key Arguments:**

    - ``extractedOrdersDF`` -- the extracted orders, with ``order``, ``wavelengthMean`` and optionally ``pixelScaleNm``
    - ``corrections`` -- ``{order: order_sky_correction}``

    **Return:**

    - ``shifted`` -- a copy of the data-frame with the corrected wavelengths
    """
    shifted = extractedOrdersDF.copy()
    for o, correction in corrections.items():
        rows = shifted["order"] == o
        wavelengths = pd.to_numeric(shifted.loc[rows, "wavelengthMean"], errors="coerce").to_numpy(dtype=float)
        scales = None
        if "pixelScaleNm" in shifted.columns:
            scales = pd.to_numeric(shifted.loc[rows, "pixelScaleNm"], errors="coerce").to_numpy(dtype=float)
        shifted.loc[rows, "wavelengthMean"] = wavelengths + correction.correction_nm(wavelengths, scales)
    return shifted


@dataclass(frozen=True, eq=False)
class nir_order_sky_lines:
    """*the detected sky-line centres and the catalogue lines of one order, found once for the iterated matching*

    The arrays are copied and made read-only.

    **Key Arguments:**

    - ``order`` -- the order
    - ``peakNm`` -- the detected line centres (nm) on the arc-solution wavelength scale
    - ``catalogueNm`` -- the calibration catalogue wavelengths (nm) in the order's range
    - ``pixelScaleNm`` -- ``pixelScaleNm`` at each detected centre
    - ``toleranceNm`` -- the match tolerance (nm) of the order
    """

    order: int
    peakNm: np.ndarray
    catalogueNm: np.ndarray
    pixelScaleNm: np.ndarray
    toleranceNm: float

    def __post_init__(self) -> None:
        for name in ("peakNm", "catalogueNm", "pixelScaleNm"):
            array = np.array(getattr(self, name), dtype=float)
            array.flags.writeable = False
            object.__setattr__(self, name, array)


def _empty_lines() -> pd.DataFrame:
    return pd.DataFrame({"order": [], "catalogueNm": [], "observedNm": []})


def _match_order_lines(orderLines, correction, toleranceNm=None) -> pd.DataFrame:
    """*the lines of one order matched once the peaks are moved by the correction, as rows for ``fit_nir_sky_model``*

    A peak with a non-finite centre or pixel scale takes no part. The tolerance is the order's own unless
    ``toleranceNm`` is given.
    """
    peaks, scales = orderLines.peakNm, orderLines.pixelScaleNm
    usable = np.flatnonzero(np.isfinite(peaks) & np.isfinite(scales) & (scales > 0))
    shiftNm = correction.correction_nm(peaks[usable], scales[usable])
    peakIndex, catalogueIndex = match_peaks_to_catalogue(
        peaks[usable] + shiftNm, orderLines.catalogueNm, orderLines.toleranceNm if toleranceNm is None else toleranceNm
    )
    peakIndex = usable[peakIndex]
    return pd.DataFrame(
        {
            "order": np.full(peakIndex.size, orderLines.order, dtype=int),
            "catalogueNm": orderLines.catalogueNm[catalogueIndex],
            "observedNm": peaks[peakIndex],
            "pixelScaleNm": scales[peakIndex],
        }
    )


def _rigid_correction(order, shiftPx, geometry) -> order_sky_correction:
    """*a correction that moves every row of the order by the same detector shift*"""
    return order_sky_correction(order, float(shiftPx), (0.0,), 0.0, 0, geometry)


class seed_match(NamedTuple):
    """The lines an order matched in the first round, how many each seed matched, and whether it keeps its own seed."""

    matched: pd.DataFrame
    globalLines: int
    ownLines: int | None
    keepsOwnSeed: bool


def _seed_scatter_nm(matched, correction) -> float:
    """*the rms (nm) of the matched lines about the seed's correction, infinite when nothing matched*"""
    if matched.empty:
        return float("inf")
    residualNm = (
        matched["catalogueNm"]
        - matched["observedNm"]
        - correction.correction_nm(matched["observedNm"].to_numpy(), matched["pixelScaleNm"].to_numpy())
    )
    return float(np.sqrt(np.mean(residualNm**2)))


def _seeded_match(orderLines, geometry, seedPx, ownSeedPx) -> seed_match:
    """*the order's lines matched about the seed, global or its own, that matches more of them within the cap*

    The two seeds are compared by the lines each matches within ``MAX_TIGHT_TOLERANCE_NM``, not within the order's tight
    tolerance. The broad lines of a wide slit have scattered centres, so the tight tolerance can match few or none
    about either seed, and one comparison there cannot tell an order with a genuine offset (its own seed right, the
    global one far off) from an order measured on aliased lines (the global seed right, its own seed off). The cap can.
    The tie goes to the lower scatter of the lines matched about the seed, then to the global seed. The chosen seed then
    matches the order's lines at its own tight tolerance.

    The order keeps its own seed (``keepsOwnSeed``) when that seed matches strictly more lines than the global seed
    within the cap, or when neither matches any: the lines of an order with a curved arc-solution error sit beyond the
    cap of every seed, and the global value would extrapolate several px. A tie at one line or more still goes to the
    global seed.

    **Return:**

    - the ``seed_match``, whose line counts are those within the cap
    """
    corrections = [
        _rigid_correction(orderLines.order, seed, geometry) for seed in (seedPx, ownSeedPx) if seed is not None
    ]
    atCap = [_match_order_lines(orderLines, c, MAX_TIGHT_TOLERANCE_NM) for c in corrections]
    scores = [(len(m), -_seed_scatter_nm(m, c)) for m, c in zip(atCap, corrections, strict=True)]
    chosen = 1 if len(scores) > 1 and scores[1] > scores[0] else 0
    counts = [len(m) for m in atCap]
    ownLines = counts[1] if len(counts) > 1 else None
    keepsOwn = ownLines is not None and (ownLines > counts[0] or ownLines == counts[0] == 0)
    return seed_match(_match_order_lines(orderLines, corrections[chosen]), counts[0], ownLines, keepsOwn)


def _kept_own_seed_corrections(orderLines, geometry, seedPx, ownSeedsPx):
    """*match the first round and set aside the orders that keep their own seed*

    An order that keeps its own seed (``_seeded_match``) but has fewer than ``MIN_OFFSET_LINES`` lines cannot carry an
    offset of its own and would be clipped as an outlier of the global model. It keeps its own seed as its shift,
    unshrunk, and its lines stay out of the global fit.

    **Return:**

    - ``frames`` -- the matched lines of the other orders
    - ``kept`` -- ``{order: order_sky_correction}`` of the orders that keep their own seed
    - ``seedLines`` -- ``{order: (lines matched about the global seed, about the own seed or None)}``
    """
    frames, kept, seedLines = [], {}, {}
    for o, lines in orderLines.items():
        match = _seeded_match(lines, geometry[o], seedPx, ownSeedsPx.get(o))
        seedLines[o] = (match.globalLines, match.ownLines)
        if match.keepsOwnSeed and len(match.matched) < MIN_OFFSET_LINES:
            rigid = _rigid_correction(o, ownSeedsPx[o], geometry[o])
            kept[o] = replace(rigid, nLines=len(match.matched), source="own-seed")
        else:
            frames.append(match.matched)
    return frames, kept, seedLines


def fit_with_iterated_matching(orderLines, geometry, seedPx, settings, iterations=MATCH_ITERATIONS, ownSeedsPx=None):
    """*match the sky lines of every order and fit the global model, matching again about each fit*

    The first match moves the peaks of each order by a rigid seed: ``seedPx``, or the order's own seed from
    ``ownSeedsPx`` when that matches more of its lines (``_seeded_match``). An order that wins with its own seed but has
    too few lines for an offset keeps that seed as its final shift (``_kept_own_seed_corrections``), and so does an
    order whose lines match about neither seed. Each later match of
    the other orders moves their peaks by the previous fit, which brings in the lines of an order whose shape error put
    them beyond the tolerance of the rigid seed. A fit that falls back at any iteration ends the loop with that
    fallback.

    **Key Arguments:**

    - ``orderLines`` -- ``{order: nir_order_sky_lines}`` of the orders that may contribute lines
    - ``geometry`` -- ``{order: order_geometry}`` of every order to correct
    - ``seedPx`` -- the rigid detector shift (px) of the first match, or *None* when no order could give one
    - ``settings`` -- the ``sky_tuning_settings``
    - ``iterations`` -- the number of match-and-fit rounds. Default ``MATCH_ITERATIONS``
    - ``ownSeedsPx`` -- ``{order: shift (px)}`` that each order measured for itself (``own_seeds_px``). Default *None*

    **Return:**

    - ``result`` -- the ``sky_refit_result`` of the last round. Its ``nLines`` and ``lines`` are those of the global
      fit, which leaves out the orders that kept their own seed

    **Raises:**

    - ``ValueError`` -- if an order with lines has no geometry
    """
    if seedPx is None:
        return _fallback(_empty_lines(), "no order gave a seed shift")
    missing = sorted(set(orderLines) - set(geometry))
    if missing:
        raise ValueError(f"order {missing[0]} has sky lines but no geometry")
    orderLines = {o: lines for o, lines in orderLines.items() if _is_usable_geometry(geometry[o])}
    ownSeedsPx = {} if ownSeedsPx is None else ownSeedsPx
    frames, kept, seedLines = _kept_own_seed_corrections(orderLines, geometry, seedPx, ownSeedsPx)
    result = None
    for _ in range(iterations):
        result = fit_nir_sky_model(
            pd.concat(frames, ignore_index=True) if frames else _empty_lines(), geometry, settings
        )
        if result.fallbackReason is not None:
            return result
        frames = [_match_order_lines(lines, result.corrections[o]) for o, lines in orderLines.items() if o not in kept]
    corrections = {**result.corrections, **kept}
    corrections = {
        o: replace(c, globalSeedLines=seedLines[o][0], ownSeedLines=seedLines[o][1]) if o in seedLines else c
        for o, c in corrections.items()
    }
    return replace(result, corrections=corrections)
