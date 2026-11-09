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


def flat_slit_ratio(log, order, flat, slit, wavelength, fitRange, isSky):
    """*the ratio 1 / P(s) of one order, from its master flat, normalised to mean 1 over the sky pixels*

    P(s) is a quadratic fitted to the flat's slit profile (see ``chunk_normalised_flat_profile``) inside ``fitRange``
    only, so the roll-off at the slit ends does not reach the fit. A pixel outside ``fitRange`` takes the ratio of the
    nearest limit. The ratio is 1 for every pixel, with the reason logged as a warning, when there is no flat, when the
    order has no fit range or no sky pixels, when too few slit bins hold flat pixels, when the fitted profile is not
    positive anywhere in the fit range, or when its gradient is too steep (``MAX_SLOPE_PER_ARCSEC`` for the mean
    gradient across the fit range, ``MAX_LOCAL_GRADIENT_PER_ARCSEC`` anywhere in it).

    **Key Arguments:**

    - ``log`` -- logger
    - ``order`` -- the order number, for the log
    - ``flat`` -- the master flat value of each pixel of the order, or None when the frame has no master flat
    - ``slit`` -- the slit position of each pixel of the order in arcsec, the shape of ``flat``
    - ``wavelength`` -- the wavelength of each pixel of the order, the shape of ``flat``
    - ``fitRange`` -- the ``(lower, upper)`` slit range in arcsec to fit over, from ``flat_fit_ranges_from_result``, or
      None
    - ``isSky`` -- boolean array, the shape of ``slit``, True for the pixels the sky model is fitted to

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
    ratio, slope, reason = _ratio_from_flat(flat, slit, wavelength, fitRange, isSky)
    if reason:
        log.warning(f"order {order}: the slit ratio from the master flat is 1, because {reason}")
        return np.ones(slit.shape)
    log.info(
        f"order {order}: the master flat gives a gradient of {slope * 100:+.2f} %/arcsec along the slit, "
        f"over {fitRange[0]:+.2f} to {fitRange[1]:+.2f} arcsec"
    )
    return ratio


def _ratio_from_flat(flat, slit, wavelength, fitRange, isSky):
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
    ratio = np.where(np.isfinite(slit), 1.0 / np.polyval(coefficients, np.clip(slit, fitLower, fitUpper)), 1.0)
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
