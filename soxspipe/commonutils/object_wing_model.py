"""Object wing model: the object light that a wavelength-only sky fit absorbs.

The stare sky model S(λ) is a b-spline in wavelength only. In a bright standard, object wing light outside the object
mask enters the fit, so the model carries object light that is then subtracted from the object. This module estimates
that light from the frame's own spatial profile. It is pure numpy on 1-D arrays, one value per pixel of one order, in
any pixel order, and it returns new objects without modifying its inputs.

The estimate has five steps:

1. The residual is the flux minus the first sky fit. Its brightest slit bin gives the object centre.
2. The residual of the pixels near the centre gives the object spectrum G(λ), by a rolling median along wavelength.
3. The residual over G, binned along the slit, gives a non-parametric spatial profile P(s). Only the pixels the
   first sky fit used go into it, because a clipped fit did not see the brighter wing pixels.
4. The faintest 1 arcsec of the profile is taken to hold no object light. Its level P0 is the zero point.
5. If the zero point is negative and significant, the first sky fit holds object light of -P0·G(λ). The wing model is
   G(λ)·(P(s) − P0), which is to be subtracted from the sky pixels before the sky is refitted.

A negative zero point that is small against the sky or against its noise is not corrected.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

PROFILE_BIN_ARCSEC = 0.25
CORE_HALF_WIDTH_ARCSEC = 0.5
FAINT_RUN_ARCSEC = 1.0
CORE_MEDIAN_WINDOW = 201
MIN_CORE_PIXELS = 201
MIN_PIXELS_PER_PROFILE_BIN = 20
WING_BIAS_SIGNIFICANCE = 5


@dataclass(frozen=True, eq=False)
class object_wing_estimate:
    """*the object light found in a first sky fit, and the profile it was found from*

    **Key Arguments:**

    - ``isCorrected`` -- True if the sky fit holds object light that should be removed
    - ``reason`` -- why no correction is made. Empty if ``isCorrected``
    - ``centre`` -- the slit position of the object centre in arcsec
    - ``binCentres`` -- the centre of every profile bin in arcsec
    - ``profile`` -- the profile P in every bin, in units of the core spectrum, from the pixels the first sky fit
      used. NaN for empty bins
    - ``zeroPoint`` -- the level P0 of the faintest part of the profile
    - ``coreSpectrum`` -- the object spectrum G at every input pixel
    - ``biasLevel`` -- the object light in the sky fit at the median core spectrum, -P0 times that median
    - ``biasNoise`` -- the noise of ``biasLevel``
    - ``skyLevel`` -- the median sky fit over the usable pixels
    - ``slitLossFraction`` -- the estimated share of the slit's object light that the sky fit removes
    """

    isCorrected: bool
    reason: str
    centre: float
    binCentres: np.ndarray
    profile: np.ndarray
    zeroPoint: float
    coreSpectrum: np.ndarray
    biasLevel: float
    biasNoise: float
    skyLevel: float
    slitLossFraction: float


def _slit_bin_layout(slitPosition):
    """*bin the slit positions into ``PROFILE_BIN_ARCSEC`` bins that start at the lowest position*

    **Key Arguments:**

    - ``slitPosition`` -- the slit position of each pixel in arcsec. Must not be empty

    **Return:**

    - ``binIndex`` -- the bin number of each pixel
    - ``binCentres`` -- the centre of every bin from the lowest position to the highest, in arcsec
    """
    import numpy as np

    lowest = slitPosition.min()
    binIndex = np.floor((slitPosition - lowest) / PROFILE_BIN_ARCSEC + 0.5).astype(int)
    binCentres = lowest + PROFILE_BIN_ARCSEC * np.arange(binIndex.max() + 1)
    return binIndex, binCentres


def profile_peak(slitPosition, residual):
    """*find the slit position of the brightest bin in the median residual profile*

    Bins with fewer than ``MIN_PIXELS_PER_PROFILE_BIN`` pixels are ignored, as one noisy pixel at the slit end could
    otherwise win.

    **Key Arguments:**

    - ``slitPosition`` -- the slit position in arcsec of each usable pixel
    - ``residual`` -- the flux minus the first sky fit at each usable pixel

    **Return:**

    - ``centre`` -- the centre in arcsec of the bin with the highest median residual. NaN if no bin has enough pixels
    """
    import numpy as np
    import pandas as pd

    binIndex, binCentres = _slit_bin_layout(slitPosition)
    nBins = len(binCentres)
    counts = np.bincount(binIndex, minlength=nBins)
    medians = pd.Series(residual).groupby(binIndex).median().reindex(range(nBins)).to_numpy()
    medians[counts < MIN_PIXELS_PER_PROFILE_BIN] = np.nan
    if np.isnan(medians).all():
        return float("nan")
    return float(binCentres[np.nanargmax(medians)])


def core_spectrum(wavelength, slitPosition, residual, isUsable, centre):
    """*measure the object spectrum in the core of the slit profile and place it on every pixel*

    The core pixels are the usable pixels within ``CORE_HALF_WIDTH_ARCSEC`` of the object centre. They are sorted by
    wavelength, and a centred rolling median of their residual over ``CORE_MEDIAN_WINDOW`` pixels gives the spectrum.
    It is interpolated to the wavelength of every input pixel.

    **Key Arguments:**

    - ``wavelength`` -- the wavelength of each pixel, in any order
    - ``slitPosition`` -- the slit position of each pixel in arcsec
    - ``residual`` -- the flux minus the first sky fit at each pixel
    - ``isUsable`` -- True for the pixels that may be used
    - ``centre`` -- the slit position of the object centre in arcsec

    **Return:**

    - ``spectrum`` -- the core spectrum at the wavelength of every input pixel
    - ``nCorePixels`` -- the number of core pixels the spectrum is built from. The spectrum is all NaN if this is 0
    """
    import numpy as np
    import pandas as pd

    isCore = isUsable & (np.abs(slitPosition - centre) < CORE_HALF_WIDTH_ARCSEC)
    nCorePixels = int(np.count_nonzero(isCore))
    if nCorePixels == 0:
        return np.full(wavelength.shape, np.nan), nCorePixels
    order = np.argsort(wavelength[isCore], kind="stable")
    coreWavelength = wavelength[isCore][order]
    coreResidual = pd.Series(residual[isCore][order])
    coreMedian = coreResidual.rolling(CORE_MEDIAN_WINDOW, center=True, min_periods=1).median().to_numpy()
    return np.interp(wavelength, coreWavelength, coreMedian), nCorePixels


def slit_profile(slitPosition, residual, error, coreSpectrum):
    """*measure the spatial profile of the residual, in units of the core spectrum*

    In every bin the profile is the summed residual over the summed core spectrum, and its noise is the summed error
    in quadrature over the summed core spectrum. Pixels where the core spectrum is not positive are skipped. A bin
    with fewer than ``MIN_PIXELS_PER_PROFILE_BIN`` such pixels is NaN. The bins span all the pixels passed in.

    **Key Arguments:**

    - ``slitPosition`` -- the slit position in arcsec of each usable pixel
    - ``residual`` -- the flux minus the first sky fit at each usable pixel
    - ``error`` -- the flux error at each usable pixel
    - ``coreSpectrum`` -- the core spectrum at each usable pixel

    **Return:**

    - ``binCentres`` -- the centre of every slit bin in arcsec
    - ``profile`` -- the profile in every bin
    - ``profileNoise`` -- the noise of the profile in every bin
    """
    import numpy as np

    binIndex, binCentres = _slit_bin_layout(slitPosition)
    nBins = len(binCentres)
    isPositive = coreSpectrum > 0
    counts = np.bincount(binIndex[isPositive], minlength=nBins)
    coreSum = np.bincount(binIndex[isPositive], weights=coreSpectrum[isPositive], minlength=nBins)
    residualSum = np.bincount(binIndex[isPositive], weights=residual[isPositive], minlength=nBins)
    varianceSum = np.bincount(binIndex[isPositive], weights=error[isPositive] ** 2, minlength=nBins)
    isFilled = counts >= MIN_PIXELS_PER_PROFILE_BIN
    profile = np.full(nBins, np.nan)
    profileNoise = np.full(nBins, np.nan)
    profile[isFilled] = residualSum[isFilled] / coreSum[isFilled]
    profileNoise[isFilled] = np.sqrt(varianceSum[isFilled]) / coreSum[isFilled]
    return binCentres, profile, profileNoise


def faint_end_zero_point(binCentres, profile, profileNoise):
    """*find the faintest part of the slit profile and its level*

    Every run of ``FAINT_RUN_ARCSEC`` of consecutive bins that are all finite is a candidate. The zero point is the
    lowest mean of the profile over a run, and its noise is the quadrature sum of the bin noises over that run,
    divided by the number of bins.

    **Key Arguments:**

    - ``binCentres`` -- the centre of every slit bin in arcsec
    - ``profile`` -- the profile in every bin, with NaN for empty bins
    - ``profileNoise`` -- the noise of the profile in every bin

    **Return:**

    - ``zeroPoint`` -- the mean profile over the faintest run. NaN if no run is complete
    - ``zeroNoise`` -- the noise of ``zeroPoint``
    - ``runCentre`` -- the slit position in arcsec of the middle of the faintest run
    """
    import numpy as np
    from numpy.lib.stride_tricks import sliding_window_view

    nRun = round(FAINT_RUN_ARCSEC / PROFILE_BIN_ARCSEC)
    if len(profile) < nRun:
        return float("nan"), float("nan"), float("nan")
    windows = sliding_window_view(profile, nRun)
    isComplete = np.isfinite(windows).all(axis=1)
    if not isComplete.any():
        return float("nan"), float("nan"), float("nan")
    means = np.where(isComplete, windows.mean(axis=1), np.inf)
    start = int(np.argmin(means))
    noise = np.sqrt(np.sum(profileNoise[start : start + nRun] ** 2)) / nRun
    runCentre = binCentres[start : start + nRun].mean()
    return float(means[start]), float(noise), float(runCentre)


def _slit_loss_fraction(profile, zeroPoint):
    """*estimate the share of the slit's object light that the sky fit removes*

    **Key Arguments:**

    - ``profile`` -- the profile in every bin, with NaN for empty bins
    - ``zeroPoint`` -- the level of the faintest part of the profile

    **Return:**

    - ``slitLossFraction`` -- the removed light over the object light, summed over the finite bins. NaN if the profile
      holds no light above the zero point
    """
    import numpy as np

    finite = profile[np.isfinite(profile)]
    objectLight = float(np.sum(finite - zeroPoint))
    if not objectLight > 0:
        return float("nan")
    return float(-zeroPoint * len(finite) / objectLight)


def _no_correction_reason(zeroPoint, biasLevel, biasNoise, skyLevel, gateFraction):
    """*say why the zero point does not justify a correction, or return an empty string if it does*

    **Key Arguments:**

    - ``zeroPoint`` -- the level of the faintest part of the profile
    - ``biasLevel`` -- the object light in the sky fit
    - ``biasNoise`` -- the noise of ``biasLevel``
    - ``skyLevel`` -- the median sky fit
    - ``gateFraction`` -- the smallest bias, as a fraction of the sky level, that is worth correcting

    **Return:**

    - ``reason`` -- the first gate that fails, or an empty string
    """
    # THE COMPARISONS ARE WRITTEN SO THAT A NAN FAILS THEM
    if not zeroPoint < 0:
        return f"the zero point {zeroPoint:.3g} is not negative, so the sky fit holds no object light"
    if not biasLevel > gateFraction * max(skyLevel, 0.0):
        return f"the bias {biasLevel:.3g} is below the gate fraction {gateFraction:g} of the sky level {skyLevel:.3g}"
    if not biasLevel > WING_BIAS_SIGNIFICANCE * biasNoise:
        return (
            f"the bias {biasLevel:.3g} is below the significance limit of "
            f"{WING_BIAS_SIGNIFICANCE} times its noise {biasNoise:.3g}"
        )
    return ""


def _blank_estimate(reason, nPixels):
    """*build an estimate that makes no correction and holds no measurements*

    **Key Arguments:**

    - ``reason`` -- why no correction is made
    - ``nPixels`` -- the number of input pixels, for the length of the NaN core spectrum

    **Return:**

    - ``estimate`` -- an ``object_wing_estimate`` with ``isCorrected`` False and every measurement NaN or empty
    """
    import numpy as np

    nan = float("nan")
    return object_wing_estimate(
        False, reason, nan, np.array([]), np.array([]), nan, np.full(nPixels, nan), nan, nan, nan, nan
    )


def estimate_object_wings(wavelength, slitPosition, flux, skyModel, error, isUsable, gateFraction, isFitted=None):
    """*estimate the object light that the first sky fit absorbed, from the frame's own spatial profile*

    **Key Arguments:**

    - ``wavelength`` -- the wavelength of each pixel of one order, in any order
    - ``slitPosition`` -- the slit position of each pixel in arcsec
    - ``flux`` -- the flux of each pixel
    - ``skyModel`` -- the first sky fit evaluated at each pixel
    - ``error`` -- the flux error of each pixel
    - ``isUsable`` -- True for the pixels that may be used. Pixels with a value that is not finite are skipped too
    - ``gateFraction`` -- the smallest bias, as a fraction of the sky level, that is worth correcting
    - ``isFitted`` -- True for the pixels the first sky fit used. The profile and its zero point come from the usable
      pixels among these, so the wing model holds the light the fit saw. Default *None*, which is every usable pixel

    **Return:**

    - ``estimate`` -- the ``object_wing_estimate``. Its ``isCorrected`` is False, with a ``reason``, when no
      correction is justified

    **Usage:**

    ```python
    from soxspipe.commonutils.object_wing_model import estimate_object_wings, wing_model
    estimate = estimate_object_wings(wavelength, slitPosition, flux, skyModel, error, isUsable, 0.01)
    wings = wing_model(estimate, slitPosition)
    ```
    """
    import numpy as np

    arrays = [wavelength, slitPosition, flux, skyModel, error, isUsable]
    if isFitted is not None:
        arrays.append(isFitted)
    if any(array.shape != wavelength.shape for array in arrays):
        raise ValueError("the pixel arrays must all have the same shape")
    isUsable = isUsable & np.isfinite(wavelength) & np.isfinite(slitPosition) & np.isfinite(flux)
    isUsable = isUsable & np.isfinite(skyModel) & np.isfinite(error)
    estimate = _blank_estimate("there are no usable pixels", wavelength.size)
    if not isUsable.any():
        return estimate
    residual = flux - skyModel
    centre = profile_peak(slitPosition[isUsable], residual[isUsable])
    if np.isnan(centre):
        return replace(estimate, reason="no slit bin has enough usable pixels to find the object centre")
    spectrum, nCore = core_spectrum(wavelength, slitPosition, residual, isUsable, centre)
    estimate = replace(estimate, centre=centre, coreSpectrum=spectrum)
    if nCore < MIN_CORE_PIXELS:
        return replace(estimate, reason=f"only {nCore} core pixels, fewer than the {MIN_CORE_PIXELS} needed")
    # A CLIPPED FIT DID NOT SEE EVERY USABLE PIXEL, SO THE PROFILE COMES FROM THE PIXELS IT DID SEE
    isProfile = isUsable if isFitted is None else isUsable & isFitted
    binCentres, profile, profileNoise = slit_profile(
        slitPosition[isProfile], residual[isProfile], error[isProfile], spectrum[isProfile]
    )
    estimate = replace(estimate, binCentres=binCentres, profile=profile)
    zeroPoint, zeroNoise, _ = faint_end_zero_point(binCentres, profile, profileNoise)
    if np.isnan(zeroPoint):
        return replace(
            estimate, reason=f"no run of {FAINT_RUN_ARCSEC:g} arcsec of finite profile bins for the zero point"
        )
    medianSpectrum = float(np.median(spectrum[isUsable]))
    biasLevel = -zeroPoint * medianSpectrum
    biasNoise = zeroNoise * medianSpectrum
    skyLevel = float(np.median(skyModel[isUsable]))
    reason = _no_correction_reason(zeroPoint, biasLevel, biasNoise, skyLevel, gateFraction)
    # THE SLIT LOSS NEEDS THE WHOLE OBJECT, AND THE FIT DID NOT USE THE OBJECT CORE
    slitProfile = profile
    if isFitted is not None:
        _, slitProfile, _ = slit_profile(
            slitPosition[isUsable], residual[isUsable], error[isUsable], spectrum[isUsable]
        )
    return replace(
        estimate,
        isCorrected=reason == "",
        reason=reason,
        zeroPoint=zeroPoint,
        biasLevel=biasLevel,
        biasNoise=biasNoise,
        skyLevel=skyLevel,
        slitLossFraction=_slit_loss_fraction(slitProfile, zeroPoint),
    )


def wing_model(estimate, slitPosition):
    """*evaluate the object wing light that the sky fit holds, at every pixel*

    The wing model is G(λ)·Pc(s), where Pc is the profile minus its zero point, linearly interpolated over slit
    position from the finite bins. Empty bins are bridged, and Pc is constant beyond the first and last finite bin.
    The model is zero everywhere if the estimate is not corrected, and a value that is not finite becomes zero.

    **Key Arguments:**

    - ``estimate`` -- the ``object_wing_estimate`` from ``estimate_object_wings``
    - ``slitPosition`` -- the slit position in arcsec of the same pixels, in the same order, that the estimate used

    **Return:**

    - ``wings`` -- the wing light at each pixel, to be subtracted from the sky pixels before the sky is refitted
    """
    import numpy as np

    if slitPosition.shape != estimate.coreSpectrum.shape:
        raise ValueError("the slit positions must be those of the pixels the estimate was made from")
    if not estimate.isCorrected:
        return np.zeros(slitPosition.shape)
    isFinite = np.isfinite(estimate.profile)
    wingProfile = np.interp(
        slitPosition, estimate.binCentres[isFinite], estimate.profile[isFinite] - estimate.zeroPoint
    )
    return np.nan_to_num(estimate.coreSpectrum * wingProfile, nan=0.0, posinf=0.0, neginf=0.0)
