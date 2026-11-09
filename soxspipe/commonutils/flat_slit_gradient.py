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
    normalised_flat_profile,
)

# THE LARGEST MEAN GRADIENT OF THE FLAT ACROSS THE KEPT RANGE THAT IS TRUSTED, AS A FRACTION PER ARCSEC. THE MEASURED
# NIR GRADIENT IS ABOUT 0.0035 (DY-1990), SO A LARGER ONE MEANS A BAD FLAT
MAX_SLOPE_PER_ARCSEC = 0.015
PROFILE_DEGREE = 2
# THE FITTED PROFILE IS CHECKED ON THIS MANY POINTS ACROSS THE KEPT RANGE
CHECK_POINTS = 101
NO_RANGE_REASON = "the order has no kept slit range"
NO_SKY_REASON = "the order has no sky pixels"
NOT_POSITIVE_REASON = "the fitted flat profile is not positive in the kept slit range"


def flat_slit_ratio(log, order, flat, slit, keepRange, isSky):
    """*the ratio 1 / P(s) of one order, from its master flat, normalised to mean 1 over the sky pixels*

    P(s) is a quadratic fitted to the flat's slit profile inside ``keepRange`` only, so the roll-off at the slit ends
    does not reach the fit. A pixel outside ``keepRange`` takes the ratio of the nearest limit. The ratio is 1 for
    every pixel, with the reason logged as a warning, when the order has no kept range or no sky pixels, when too few
    slit bins hold flat pixels, when the fitted profile is not positive, or when the mean gradient of the profile
    across the kept range is above ``MAX_SLOPE_PER_ARCSEC``.

    **Key Arguments:**

    - ``log`` -- logger
    - ``order`` -- the order number, for the log
    - ``flat`` -- the master flat value of each pixel of the order
    - ``slit`` -- the slit position of each pixel of the order in arcsec, the shape of ``flat``
    - ``keepRange`` -- the ``(lower, upper)`` both-lit slit range in arcsec from ``kept_slit_range``, or None
    - ``isSky`` -- boolean array, the shape of ``flat``, True for the pixels the sky model is fitted to

    **Return:**

    - ``ratio`` -- the ratio of each pixel, a new array the shape of ``flat``

    **Usage:**

    ```python
    from soxspipe.commonutils.flat_slit_gradient import flat_slit_ratio
    ratio = flat_slit_ratio(log, 12, flat, slit, (-4.9, 5.1), isSky)
    ```
    """
    flat = np.asarray(flat, dtype=float)
    slit = np.asarray(slit, dtype=float)
    isSky = np.asarray(isSky, dtype=bool) & np.isfinite(slit)
    unitRatio = np.ones(slit.shape)
    ratio, slope, reason = _ratio_from_flat(flat, slit, keepRange, isSky)
    if reason:
        log.warning(f"order {order}: the slit ratio from the master flat is 1, because {reason}")
        return unitRatio
    log.info(f"order {order}: the master flat gives a gradient of {slope * 100:+.2f} % per arcsec along the slit")
    return ratio


def _ratio_from_flat(flat, slit, keepRange, isSky):
    """*the normalised ratio, the mean gradient of the flat profile and the reason the ratio is unusable (or "")*"""
    if keepRange is None:
        return None, np.nan, NO_RANGE_REASON
    if not isSky.any():
        return None, np.nan, NO_SKY_REASON
    keepLower, keepUpper = keepRange
    isKept = np.isfinite(flat) & np.isfinite(slit) & (slit >= keepLower) & (slit <= keepUpper)
    if not isKept.any():
        return None, np.nan, TOO_FEW_FLAT_BINS_REASON
    centres, profile, reason = normalised_flat_profile(flat, slit, isKept, PROFILE_BIN_ARCSEC)
    if reason:
        return None, np.nan, reason
    isBin = np.isfinite(profile)
    coefficients = np.polyfit(centres[isBin], profile[isBin], PROFILE_DEGREE)
    check = np.polyval(coefficients, np.linspace(keepLower, keepUpper, CHECK_POINTS))
    if check.min() <= 0:
        return None, np.nan, NOT_POSITIVE_REASON
    slope = (check[-1] - check[0]) / (keepUpper - keepLower) / check.mean()
    if abs(slope) > MAX_SLOPE_PER_ARCSEC:
        return (
            None,
            slope,
            (
                f"the flat slope of {slope * 100:+.2f} % per arcsec is beyond the limit of "
                f"{MAX_SLOPE_PER_ARCSEC * 100:.1f} % per arcsec"
            ),
        )
    ratio = np.where(np.isfinite(slit), 1.0 / np.polyval(coefficients, np.clip(slit, keepLower, keepUpper)), 1.0)
    return ratio / ratio[isSky].mean(), slope, ""
