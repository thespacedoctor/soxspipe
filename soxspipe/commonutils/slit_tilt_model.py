"""Slit-tilt model: the part of the wavelength map's tilt along the slit that a frame's own skylines can measure.

The stare sky model S(λ) is a b-spline in wavelength only. When the 2D map puts a skyline at slightly different
wavelengths at different slit positions, the line lands on different map wavelengths along the slit, and the residual
of a pixel is close to S'(λ)·δλ, where δλ is the map error at that pixel (DY-2065). This module fits

    δλ(λ, s) = s · (c0 + c1 · u),    u = 2 (λ − λmin) / (λmax − λmin) − 1

to that residual by weighted least squares, so a tilt that changes sign along the order (a twist) is measured as
well as a uniform tilt. It is pure numpy on 1-D arrays, one value per pixel of one order, in any pixel order, and it
returns new objects without modifying its inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

# THE FEWEST USABLE PIXELS THAT MEASURE A TILT; AN ORDER WITH FEWER RAISES
MIN_FIT_PIXELS = 100
# THE NUMBER OF TILT COEFFICIENTS, C0 AND C1; A DESIGN OF LOWER RANK CANNOT MEASURE BOTH
N_COEFFICIENTS = 2


@dataclass(frozen=True)
class SlitTiltEstimate:
    """*the fitted slit tilt of one order*

    **Key Arguments:**

    - ``c0`` -- the mean tilt, in wavelength units per unit of slit position
    - ``c1`` -- the change of tilt from the order centre to its red end, in the same units
    - ``minWavelength`` -- the smallest wavelength of the order, which sets u = -1
    - ``maxWavelength`` -- the largest wavelength of the order, which sets u = +1
    """

    c0: float
    c1: float
    minWavelength: float
    maxWavelength: float


def estimate_slit_tilt(
    wavelength: np.ndarray,
    slitPosition: np.ndarray,
    flux: np.ndarray,
    skyModel: np.ndarray,
    skyDerivative: np.ndarray,
    error: np.ndarray,
    isFitted: np.ndarray,
) -> SlitTiltEstimate:
    """*fit the slit tilt of one order to the residual of its first sky fit*

    **Key Arguments:**

    - ``wavelength`` -- the map wavelength of each pixel
    - ``slitPosition`` -- the slit position of each pixel
    - ``flux`` -- the flux of each pixel
    - ``skyModel`` -- the first sky fit at each pixel
    - ``skyDerivative`` -- the wavelength derivative of the first sky fit at each pixel
    - ``error`` -- the flux error of each pixel
    - ``isFitted`` -- True for the pixels the first sky fit used

    **Return:**

    - ``estimate`` -- the fitted ``SlitTiltEstimate``

    **Raises:**

    - ``ValueError`` -- when the inputs differ in length, the order spans no wavelength range, fewer than
      ``MIN_FIT_PIXELS`` pixels are usable, or the pixels cannot measure both coefficients (for example a flat sky)

    **Usage:**

        estimate = estimate_slit_tilt(wavelength, slitPosition, flux, skyModel, skyDerivative, error, isFitted)
    """
    import numpy as np

    arrays = [np.asarray(a, dtype=float) for a in (wavelength, slitPosition, flux, skyModel, skyDerivative, error)]
    wavelength, s, flux, skyModel, d, error = arrays
    isFitted = np.asarray(isFitted, dtype=bool)
    if len({a.size for a in arrays} | {isFitted.size}) != 1:
        raise ValueError("the slit-tilt inputs must all have one value per pixel")
    isFiniteWavelength = np.isfinite(wavelength)
    if not isFiniteWavelength.any() or np.ptp(wavelength[isFiniteWavelength]) <= 0:
        raise ValueError("the order spans no wavelength range, so its slit tilt cannot be measured")
    minWavelength = float(wavelength[isFiniteWavelength].min())
    maxWavelength = float(wavelength[isFiniteWavelength].max())
    u = _order_coordinate(wavelength, minWavelength, maxWavelength)
    residual = flux - skyModel
    with np.errstate(divide="ignore"):
        weight = 1 / error
    isUsed = isFitted & np.isfinite(residual) & np.isfinite(d) & np.isfinite(weight) & np.isfinite(s) & np.isfinite(u)
    if isUsed.sum() < MIN_FIT_PIXELS:
        raise ValueError(f"only {isUsed.sum()} pixels are usable, fewer than the {MIN_FIT_PIXELS} a tilt needs")

    design = np.vstack([s * d, s * u * d]).T[isUsed] * weight[isUsed, None]
    coefficients, _, rank, _ = np.linalg.lstsq(design, residual[isUsed] * weight[isUsed], rcond=None)
    if rank < N_COEFFICIENTS:
        raise ValueError("the usable pixels cannot measure the slit tilt (no skyline flanks along the slit)")
    return SlitTiltEstimate(float(coefficients[0]), float(coefficients[1]), minWavelength, maxWavelength)


def corrected_wavelength(estimate: SlitTiltEstimate, wavelength: np.ndarray, slitPosition: np.ndarray) -> np.ndarray:
    """*the map wavelengths moved by the fitted slit tilt, so each skyline lands on one wavelength along the slit*

    **Key Arguments:**

    - ``estimate`` -- the ``SlitTiltEstimate`` of the order
    - ``wavelength`` -- the map wavelength of each pixel
    - ``slitPosition`` -- the slit position of each pixel

    **Return:**

    - ``corrected`` -- a new array of wavelengths, λ + s · (c0 + c1 · u)

    **Usage:**

        corrected = corrected_wavelength(estimate, wavelength, slitPosition)
    """
    import numpy as np

    wavelength = np.asarray(wavelength, dtype=float)
    u = _order_coordinate(wavelength, estimate.minWavelength, estimate.maxWavelength)
    return wavelength + np.asarray(slitPosition, dtype=float) * (estimate.c0 + estimate.c1 * u)


def _order_coordinate(wavelength: np.ndarray, minWavelength: float, maxWavelength: float) -> np.ndarray:
    """*the wavelength scaled to u = -1 at the blue end of the order and +1 at its red end*"""
    return 2 * (wavelength - minWavelength) / (maxWavelength - minWavelength) - 1
