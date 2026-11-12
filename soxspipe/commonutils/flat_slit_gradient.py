"""The along-slit sky illumination taken from the master flat (DY-2147).

The QTH lamp lights the slit differently from the sky, and the science frame is divided by the flat. The sky in the
flat-fielded frame therefore carries a smooth gradient along the slit that a model of wavelength alone cannot follow.
For one order, this module fits a quadratic P(s) in slit position s to the flat's slit profile, inside the range
where both the sky and the flat are lit, and returns the ratio 1 / P(s) that the sky fit divides out. This module is
pure numpy and returns new objects without modifying its inputs.
"""

import numpy as np

from soxspipe.commonutils.slit_illumination_profiles import (
    PROFILE_BIN_ARCSEC,
    TOO_FEW_FLAT_BINS_REASON,
    chunk_normalised_flat_profile,
)

# THE LARGEST MEAN GRADIENT OF THE FITTED PROFILE ACROSS THE FIT RANGE, AS A FRACTION OF ITS MEAN LEVEL PER ARCSEC, THAT
# IS TRUSTED. THE MEASURED NIR GRADIENT IS ABOUT 0.0035 (DY-1990), SO A LARGER ONE MEANS A BAD FLAT
MAX_SLOPE_PER_ARCSEC = 0.015
# THE LARGEST LOCAL GRADIENT ANYWHERE IN THE FIT RANGE, IN THE SAME UNITS. THE MASTER FLAT'S REAL CURVATURE (THE
# ROLL-OFF TOWARD THE SLIT ENDS) MAKES THE GRADIENT AT THE ENDS 0.015 TO 0.024 IN SOXS NIR ORDERS 10 TO 22 OF THE
# 2026-01-27 AND 2026-01-12 FLATS WHILE THE MEAN GRADIENT IS ABOUT 0.004. THIS LIMIT IS SET ABOVE THAT, AND BELOW THE
# 0.029 TO 0.038 OF ORDERS 23 AND 24, WHICH HAVE THE WEAKEST FLAT. IT CATCHES A PROFILE THAT BENDS STEEPLY WITH NO NET
# SLOPE (DY-2147)
MAX_LOCAL_GRADIENT_PER_ARCSEC = 0.03
PROFILE_DEGREE = 2
NO_RANGE_REASON = "the order has no slit range to fit the flat over"
NO_SKY_REASON = "the order has no sky pixels"
NO_FLAT_REASON = "there is no master flat"
NOT_POSITIVE_REASON = "the fitted flat profile is not positive in the fit range"
UNSCALED = 1.0
# THE KEY OF THE PER-ORDER FACTOR K IN THE DETECTOR PARAMETERS OF AN ARM (DY-2150). AN ARM WITHOUT THE KEY KEEPS K = 1
SCALE_KEY = "flat-slit-gradient-scale"


def gradient_scale(detectorParams: dict, order: int) -> float:
    """*the factor k by which the master flat's gradient is scaled in one order*

    The factors are the ``flat-slit-gradient-scale`` table of the arm in ``soxs_detector_parameters.yaml``, which maps
    an order number to k. Every k in the table must be a number above 0 and at most 1, so the scaled profile stays
    positive (see ``flat_slit_ratio``). A table with any other value raises, whichever order is asked for.

    **Key Arguments:**

    - ``detectorParams`` -- the detector parameters of the arm, from ``detector_lookup``
    - ``order`` -- the order number

    **Return:**

    - ``k`` -- the table entry of the order, or 1 when the arm has no table or the table has no entry for the order

    **Usage:**

    ```python
    from soxspipe.commonutils.detector_lookup import detector_lookup
    from soxspipe.commonutils.flat_slit_gradient import gradient_scale
    detectorParams = detector_lookup(log=log, settings=settings).get("NIR")
    k = gradient_scale(detectorParams, 17)
    ```
    """
    table = detectorParams.get(SCALE_KEY, {})
    if not isinstance(table, dict):
        raise ValueError(f"{SCALE_KEY} must map an order to a factor k, not {type(table).__name__}")
    for tableOrder, k in table.items():
        if isinstance(k, bool) or not isinstance(k, (int, float)) or not 0 < k <= 1:
            raise ValueError(
                f"{SCALE_KEY} for order {tableOrder} is {k!r}, but k must be a number above 0 and at most 1"
            )
    return table.get(order, UNSCALED)


def flat_slit_ratio(log, order, flat, slit, wavelength, fitRange, isSky, gradientScale=UNSCALED):
    """*the ratio 1 / P(s) of one order, from its master flat, normalised to mean 1 over the sky pixels*

    P(s) is a quadratic fitted to the flat's slit profile (see ``chunk_normalised_flat_profile``) inside ``fitRange``
    only, so the roll-off at the slit ends does not reach the fit. A pixel outside ``fitRange`` takes the ratio of the
    nearest limit. The ratio is 1 for every pixel, with the reason logged as a warning, when there is no flat, when the
    order has no fit range or no sky pixels, when too few slit bins hold flat pixels, when the fitted profile is not
    positive anywhere in the fit range, or when its gradient is too steep (``MAX_SLOPE_PER_ARCSEC`` for the mean
    gradient across the fit range, ``MAX_LOCAL_GRADIENT_PER_ARCSEC`` anywhere in it). Those tests judge the flat as it
    is, before ``gradientScale``.

    ``gradientScale`` is the factor k that scales the fitted gradient (DY-2150). The profile used is
    ``P_k(s) = M + k (P(s) - M)``, where M is the mean of P over the fit range, so its gradient everywhere, and its
    deviation from the mean, are k times those of P. ``gradient_scale`` accepts only 0 < k <= 1, which keeps the profile
    positive.

    **Key Arguments:**

    - ``log`` -- logger
    - ``order`` -- the order number, for the log
    - ``flat`` -- the master flat value of each pixel of the order, or None when the frame has no master flat
    - ``slit`` -- the slit position of each pixel of the order in arcsec, the shape of ``flat``
    - ``wavelength`` -- the wavelength of each pixel of the order, the shape of ``flat``
    - ``fitRange`` -- the ``(lower, upper)`` slit range in arcsec to fit over, from ``flat_fit_ranges_from_result``, or
      None
    - ``isSky`` -- boolean array, the shape of ``slit``, True for the pixels the sky model is fitted to
    - ``gradientScale`` -- the factor k for this order, from ``gradient_scale``. 1 (the default) leaves the flat's
      gradient as fitted

    **Return:**

    - ``ratio`` -- the ratio of each pixel, a new array the shape of ``slit``

    **Usage:**

    ```python
    from soxspipe.commonutils.flat_slit_gradient import flat_slit_ratio
    ratio = flat_slit_ratio(log, 12, flat, slit, wavelength, (-4.9, 5.1), isSky)
    ```
    """
    slit = np.asarray(slit, dtype=float)
    isSky = np.asarray(isSky, dtype=bool) & np.isfinite(slit)
    ratio, slope, reason = _ratio_from_flat(flat, slit, wavelength, fitRange, isSky, gradientScale)
    if reason:
        log.warning(f"order {order}: the slit ratio from the master flat is 1, because {reason}")
        return np.ones(slit.shape)
    log.info(
        f"order {order}: the master flat gives a gradient of {slope * 100:+.2f} %/arcsec along the slit, "
        f"over {fitRange[0]:+.2f} to {fitRange[1]:+.2f} arcsec, scaled by k = {gradientScale:.2f}"
    )
    return ratio


def _ratio_from_flat(flat, slit, wavelength, fitRange, isSky, gradientScale):
    """*the normalised ratio, the mean gradient of the flat profile and the reason the ratio is unusable (or "")*"""
    if flat is None:
        return None, np.nan, NO_FLAT_REASON
    if fitRange is None:
        return None, np.nan, NO_RANGE_REASON
    if not isSky.any():
        return None, np.nan, NO_SKY_REASON
    flat = np.asarray(flat, dtype=float)
    wavelength = np.asarray(wavelength, dtype=float)
    fitLower, fitUpper = fitRange
    isFitted = np.isfinite(slit) & (slit >= fitLower) & (slit <= fitUpper)
    if not isFitted.any():
        return None, np.nan, TOO_FEW_FLAT_BINS_REASON
    centres, profile, reason = chunk_normalised_flat_profile(flat, slit, wavelength, isFitted, PROFILE_BIN_ARCSEC)
    if reason:
        return None, np.nan, reason
    isBin = np.isfinite(profile)
    coefficients = np.polyfit(centres[isBin], profile[isBin], PROFILE_DEGREE)
    minimum, maxGradient, meanLevel, meanGradient = _profile_extremes(coefficients, fitLower, fitUpper)
    if not (np.isfinite(minimum) and minimum > 0):
        return None, np.nan, NOT_POSITIVE_REASON
    meanSlope, steepestSlope = meanGradient / meanLevel, maxGradient / meanLevel
    if abs(meanSlope) > MAX_SLOPE_PER_ARCSEC:
        return None, meanSlope, _slope_reason("mean", abs(meanSlope), MAX_SLOPE_PER_ARCSEC)
    if steepestSlope > MAX_LOCAL_GRADIENT_PER_ARCSEC:
        return None, meanSlope, _slope_reason("local", steepestSlope, MAX_LOCAL_GRADIENT_PER_ARCSEC)
    # SCALING THE DEVIATION FROM THE MEAN LEVEL SCALES THE GRADIENT AT EVERY SLIT POSITION AND KEEPS THE MEAN LEVEL
    scaled = gradientScale * coefficients
    scaled[-1] += (1 - gradientScale) * meanLevel
    ratio = np.where(np.isfinite(slit), 1.0 / np.polyval(scaled, np.clip(slit, fitLower, fitUpper)), 1.0)
    return ratio / ratio[isSky].mean(), meanSlope, ""


def _slope_reason(kind, value, limit):
    """*the reason text for a flat profile whose mean or local slope is beyond its limit*"""
    return f"the {kind} flat slope of {value * 100:.2f} %/arcsec is beyond the limit of {limit * 100:.1f} %/arcsec"


def _profile_extremes(coefficients, lower, upper):
    """*the minimum, the largest gradient (size), the mean and the mean gradient of a quadratic over [lower, upper]*"""
    a, b, c = coefficients
    candidates = [lower, upper]
    # A QUADRATIC TERM TOO SMALL TO HAVE A VERTEX INSIDE THE RANGE (OR ONE THAT OVERFLOWS) LEAVES A LINE
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        vertex = -b / (2 * a)
    if np.isfinite(vertex) and lower < vertex < upper:
        candidates.append(vertex)
    minimum = np.polyval(coefficients, candidates).min()
    maxGradient = max(abs(2 * a * lower + b), abs(2 * a * upper + b))
    meanLevel = a * (lower**2 + lower * upper + upper**2) / 3 + b * (lower + upper) / 2 + c
    meanGradient = (np.polyval(coefficients, upper) - np.polyval(coefficients, lower)) / (upper - lower)
    return minimum, maxGradient, meanLevel, meanGradient
