"""Object wing model: the object light that a wavelength-only sky fit absorbs.

The stare sky model S(λ) is a b-spline in wavelength only. In a bright standard, object wing light outside the object
mask enters the fit, so the model carries object light that is then subtracted from the object. This module estimates
that light from the frame's own spatial profile. It is pure numpy on 1-D arrays, one value per pixel of one order, in
any pixel order, and it returns new objects without modifying its inputs.

The estimate has five steps:

1. The residual is the flux minus the first sky fit. Its brightest slit bin gives the object centre.
2. The residual of the pixels near the centre gives the object spectrum G(λ), by a rolling median along wavelength.
3. The residual over G, binned along the slit, gives a non-parametric spatial profile P(s). Only the pixels the
   first sky fit used go into it, because a clipped fit did not see the brighter wing pixels. Of those, only the
   pixels with faint sky go into it, because a slit illumination error scales with the sky and would pass for
   object light.
4. The profile more than 3 arcsec from the object centre is taken to hold no object light. Its mean level P0, over
   both sides of the object, is the zero point (DY-1964). The faintest 1 arcsec, used before, picked a dim slit end
   or a low side of the slit and over-corrected.
5. If the zero point is negative and significant, the first sky fit holds object light of -P0·G(λ). The wing model is
   G(λ)·(P(s) − P0), which is to be subtracted from the sky pixels before the sky is refitted.

A negative zero point that is small against the noise of a faint-sky pixel, or against its own noise, is not
corrected. Nor is an order whose slit profile shows more than one significant source (DY-1933): on a crowded slit no
part of the slit is free of object light, so the far slit does not measure the sky.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

PROFILE_BIN_ARCSEC = 0.25
CORE_HALF_WIDTH_ARCSEC = 0.5
# THE PROFILE MORE THAN THIS FAR FROM THE OBJECT CENTRE IS TAKEN TO HOLD NO OBJECT LIGHT (DY-1964)
FAR_SLIT_ARCSEC = 3.0
# THE FAR BINS MUST COVER AT LEAST THIS MUCH SLIT FOR A ZERO POINT
FAR_SLIT_MIN_ARCSEC = 1.0
# THE PROFILE USES ONLY PIXELS WHOSE SKY MODEL IS AT OR BELOW THIS QUANTILE OF THE ORDER'S SKY, BECAUSE A SLIT
# ILLUMINATION ERROR SCALES WITH THE SKY AND THE OBJECT LIGHT DOES NOT (DY-1964)
FAINT_SKY_QUANTILE = 0.5
CORE_MEDIAN_WINDOW = 201
MIN_CORE_PIXELS = 201
MIN_PIXELS_PER_PROFILE_BIN = 20
WING_BIAS_SIGNIFICANCE = 5
# A SOURCE IS A PROFILE PEAK WHOSE PROMINENCE IS ABOVE THIS MANY TIMES ITS BIN NOISE (DY-1933)
SOURCE_SIGNIFICANCE = 5
# ... AND ABOVE THIS FRACTION OF THE PROFILE'S PEAK ABOVE ITS MEDIAN, SO A FAINT HALO BUMP NEXT TO A BRIGHT STAR IS NOT
# A SOURCE. THE MEDIAN, NOT THE MINIMUM, BECAUSE AN UNLIT SLIT END CAN SIT FAR BELOW THE REST OF THE PROFILE
MIN_SOURCE_FRACTION = 0.2
# ... AND AT LEAST THIS FAR FROM A BRIGHTER SOURCE, SO A NOISY CORE IS NOT SPLIT IN TWO
MIN_SOURCE_SEPARATION_ARCSEC = 1.5


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
    - ``zeroPoint`` -- the level P0 of the profile far from the object
    - ``coreSpectrum`` -- the object spectrum G at every input pixel
    - ``biasLevel`` -- the object light in the sky fit at the median core spectrum, -P0 times that median
    - ``biasNoise`` -- the noise of ``biasLevel``
    - ``skyLevel`` -- the median sky fit over the usable pixels
    - ``faintSkyNoise`` -- the median flux error of the usable pixels more than 3 arcsec from the object whose sky
      fit is at or below the median sky
    - ``slitLossFraction`` -- the estimated share of the slit's object light that the sky fit removes
    - ``sourceCentres`` -- the slit positions in arcsec of the significant sources in the profile of every usable
      pixel, most prominent first. Empty if the profile was not measured
    - ``isCrowdedSlit`` -- True if every other gate passed but the order is not corrected because the profile shows
      more than one source
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
    faintSkyNoise: float
    slitLossFraction: float
    sourceCentres: np.ndarray
    isCrowdedSlit: bool = False


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


def far_slit_zero_point(binCentres, profile, profileNoise, centre):
    """*find the level of the slit profile far from the object (DY-1964)*

    The zero point is the mean of the finite profile bins more than ``FAR_SLIT_ARCSEC`` from the object centre, on both
    sides of it. Its noise is the quadrature sum of the bin noises, divided by the number of bins. Fewer than
    ``FAR_SLIT_MIN_ARCSEC`` of such bins give no zero point.

    **Key Arguments:**

    - ``binCentres`` -- the centre of every slit bin in arcsec
    - ``profile`` -- the profile in every bin, with NaN for empty bins
    - ``profileNoise`` -- the noise of the profile in every bin
    - ``centre`` -- the slit position of the object centre in arcsec

    **Return:**

    - ``zeroPoint`` -- the mean profile of the far bins. NaN if there are too few
    - ``zeroNoise`` -- the noise of ``zeroPoint``
    """
    import numpy as np

    isFar = np.isfinite(profile) & (np.abs(binCentres - centre) > FAR_SLIT_ARCSEC)
    nFar = int(np.count_nonzero(isFar))
    if nFar < round(FAR_SLIT_MIN_ARCSEC / PROFILE_BIN_ARCSEC):
        return float("nan"), float("nan")
    noise = np.sqrt(np.sum(profileNoise[isFar] ** 2)) / nFar
    return float(np.mean(profile[isFar])), float(noise)


def profile_sources(binCentres, profile, profileNoise):
    """*find the significant sources in a slit profile (DY-1933)*

    A source is a local maximum of the profile whose prominence is above ``SOURCE_SIGNIFICANCE`` times the noise of
    its bin and above ``MIN_SOURCE_FRACTION`` of the profile's maximum minus its median. A source closer than
    ``MIN_SOURCE_SEPARATION_ARCSEC`` to a more prominent one is dropped. A bin that is not finite or has no positive
    noise is invalid. Invalid bins at the slit ends are dropped, and the profile and noise are linearly interpolated
    across invalid bins inside the slit, so an empty bin neither hides a peak nor makes one. A maximum in the first or
    last valid bin is not a peak.

    **Key Arguments:**

    - ``binCentres`` -- the centre of every slit bin in arcsec
    - ``profile`` -- the profile in every bin, with NaN for empty bins
    - ``profileNoise`` -- the noise of the profile in every bin

    **Return:**

    - ``sourceCentres`` -- the slit positions in arcsec of the sources, most prominent first
    """
    import numpy as np
    from scipy.signal import find_peaks

    isValid = np.isfinite(profile) & np.isfinite(profileNoise) & (profileNoise > 0)
    if not isValid.any():
        return np.array([])
    validIndex = np.flatnonzero(isValid)
    span = slice(validIndex[0], validIndex[-1] + 1)
    centres = binCentres[span]
    values = np.interp(centres, binCentres[isValid], profile[isValid])
    noise = np.interp(centres, binCentres[isValid], profileNoise[isValid])
    peaks, properties = find_peaks(values, prominence=0)
    prominences = properties["prominences"]
    isSource = (prominences > SOURCE_SIGNIFICANCE * noise[peaks]) & (
        prominences > MIN_SOURCE_FRACTION * (values.max() - np.median(values))
    )
    candidates = centres[peaks][isSource][np.argsort(-prominences[isSource], kind="stable")]
    sources = []
    for centre in candidates:
        if all(abs(centre - kept) >= MIN_SOURCE_SEPARATION_ARCSEC for kept in sources):
            sources.append(centre)
    return np.array(sources, dtype=float)


def _slit_loss_fraction(profile, zeroPoint):
    """*estimate the share of the slit's object light that the sky fit removes*

    **Key Arguments:**

    - ``profile`` -- the profile in every bin, with NaN for empty bins
    - ``zeroPoint`` -- the level of the profile far from the object

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


def _no_correction_reason(zeroPoint, biasLevel, biasNoise, faintSkyNoise, minLeakSigma):
    """*say why the zero point does not justify a correction, or return an empty string if it does*

    **Key Arguments:**

    - ``zeroPoint`` -- the level of the profile far from the object
    - ``biasLevel`` -- the object light in the sky fit
    - ``biasNoise`` -- the noise of ``biasLevel``
    - ``faintSkyNoise`` -- the typical flux error of a faint-sky pixel
    - ``minLeakSigma`` -- the smallest bias, in units of ``faintSkyNoise``, that is worth correcting

    **Return:**

    - ``reason`` -- the first gate that fails, or an empty string
    """
    # THE COMPARISONS ARE WRITTEN SO THAT A NAN FAILS THEM
    if not zeroPoint < 0:
        return f"the zero point {zeroPoint:.3g} is not negative, so the sky fit holds no object light"
    # A LEAK FAR BELOW THE NOISE OF ONE PIXEL CANNOT BE SEEN IN THE SKY-SUBTRACTED IMAGE, SO IT IS NOT WORTH A REFIT
    if not biasLevel > minLeakSigma * faintSkyNoise:
        return (
            f"the bias {biasLevel:.3g} is below {minLeakSigma:g} times the noise of a faint-sky pixel "
            f"{faintSkyNoise:.3g}"
        )
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
        isCorrected=False,
        reason=reason,
        centre=nan,
        binCentres=np.array([]),
        profile=np.array([]),
        zeroPoint=nan,
        coreSpectrum=np.full(nPixels, nan),
        biasLevel=nan,
        biasNoise=nan,
        skyLevel=nan,
        faintSkyNoise=nan,
        slitLossFraction=nan,
        sourceCentres=np.array([]),
    )


def estimate_object_wings(wavelength, slitPosition, flux, skyModel, error, isUsable, minLeakSigma, isFitted=None):
    """*estimate the object light that the first sky fit absorbed, from the frame's own spatial profile*

    **Key Arguments:**

    - ``wavelength`` -- the wavelength of each pixel of one order, in any order
    - ``slitPosition`` -- the slit position of each pixel in arcsec
    - ``flux`` -- the flux of each pixel
    - ``skyModel`` -- the first sky fit evaluated at each pixel
    - ``error`` -- the flux error of each pixel
    - ``isUsable`` -- True for the pixels that may be used. Pixels with a value that is not finite are skipped too
    - ``minLeakSigma`` -- the smallest bias, in units of the flux error of a faint-sky pixel, that is worth correcting
    - ``isFitted`` -- True for the pixels the first sky fit used. The profile and its zero point come from the usable
      pixels among these, so the wing model holds the light the fit saw. Default *None*, which is every usable pixel

    **Return:**

    - ``estimate`` -- the ``object_wing_estimate``. Its ``isCorrected`` is False, with a ``reason``, when no
      correction is justified

    **Usage:**

    ```python
    from soxspipe.commonutils.object_wing_model import estimate_object_wings, wing_model
    estimate = estimate_object_wings(wavelength, slitPosition, flux, skyModel, error, isUsable, 0.03)
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
    if not isProfile.any():
        return replace(estimate, reason="the first sky fit used none of the usable pixels")
    isFaintSky = isUsable & (skyModel <= np.quantile(skyModel[isUsable], FAINT_SKY_QUANTILE))
    isProfile = isProfile & isFaintSky
    if not isProfile.any():
        return replace(estimate, reason="the first sky fit used no usable pixel with faint sky")
    binCentres, profile, profileNoise = slit_profile(
        slitPosition[isProfile], residual[isProfile], error[isProfile], spectrum[isProfile]
    )
    estimate = replace(estimate, binCentres=binCentres, profile=profile)
    zeroPoint, zeroNoise = far_slit_zero_point(binCentres, profile, profileNoise, centre)
    if np.isnan(zeroPoint):
        return replace(
            estimate,
            reason=(
                f"fewer than {FAR_SLIT_MIN_ARCSEC:g} arcsec of finite profile bins more than "
                f"{FAR_SLIT_ARCSEC:g} arcsec from the object for the zero point"
            ),
        )
    medianSpectrum = float(np.median(spectrum[isUsable]))
    biasLevel = -zeroPoint * medianSpectrum
    biasNoise = zeroNoise * medianSpectrum
    skyLevel = float(np.median(skyModel[isUsable]))
    # THE NOISE OF THE SKY, NOT OF THE OBJECT: THE BRIGHT CORE WOULD MAKE THE LEAK GATE STRICTER
    isFarFaintSky = isFaintSky & (np.abs(slitPosition - centre) > FAR_SLIT_ARCSEC)
    faintSkyNoise = float(np.median(error[isFarFaintSky if isFarFaintSky.any() else isFaintSky]))
    reason = _no_correction_reason(zeroPoint, biasLevel, biasNoise, faintSkyNoise, minLeakSigma)
    # THE SLIT LOSS AND THE SOURCE COUNT NEED THE WHOLE OBJECT ON EVERY USABLE PIXEL, AND THE PROFILE ABOVE HAS
    # NEITHER THE OBJECT CORE (WHICH THE FIT DID NOT USE) NOR THE BRIGHT-SKY PIXELS
    slitBins, slitProfile, slitNoise = slit_profile(
        slitPosition[isUsable], residual[isUsable], error[isUsable], spectrum[isUsable]
    )
    sourceCentres = profile_sources(slitBins, slitProfile, slitNoise)
    isCrowdedSlit = reason == "" and len(sourceCentres) > 1
    if isCrowdedSlit:
        positions = ", ".join(f"{centre:+.2f}" for centre in sourceCentres)
        reason = (
            f"the slit profile shows {len(sourceCentres)} sources (at {positions} arcsec), so no part of the slit "
            "is free of object light"
        )
    return replace(
        estimate,
        isCorrected=reason == "",
        reason=reason,
        zeroPoint=zeroPoint,
        biasLevel=biasLevel,
        biasNoise=biasNoise,
        skyLevel=skyLevel,
        faintSkyNoise=faintSkyNoise,
        slitLossFraction=_slit_loss_fraction(slitProfile, zeroPoint),
        sourceCentres=sourceCentres,
        isCrowdedSlit=isCrowdedSlit,
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
