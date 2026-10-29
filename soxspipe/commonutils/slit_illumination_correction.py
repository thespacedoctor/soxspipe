"""Skyline-based multiplicative slit-illumination correction (DY-2054 prototype).

The stare sky model is a b-spline in wavelength only, so it assumes the sky is equally bright along the slit. A small
illumination gradient along the slit leaves a residual that scales with the sky, so it is largest at the bright
skylines. This module measures that gradient as the Kelson ``/fixillum`` step does (``longslit_reduce.pro``, lines 3709
to 3770), and returns the factor to multiply the sky model by:

1. The fractional residual of each pixel is (data minus object model minus sky model) over the sky model.
2. Only bright skyline pixels are kept: a sky model of at least ``MIN_SKY_SNR`` times the pixel error, a fractional
   residual of at most ``MAX_FRACTIONAL_RESIDUAL``, and a wavelength inside the central 90% of the order. The IDL
   uses a fixed 1000 counts. That is about 45 times the shot noise on X-Shooter NIR, so the same rule on SOXS, with
   its shorter exposures and its different gain, is a limit in units of the pixel error.
3. The mean fractional residual is taken in slit-position bins of ``BIN_WIDTH_ARCSEC``. A bin with
   ``MIN_PIXELS_PER_SLIT_BIN`` kept pixels or fewer is not used, and at least ``MIN_BINS`` bins must be used.
4. The bin means are interpolated onto every pixel, and held flat outside the outer bins used. The factor is one plus
   that curve, less its mean. The mean level of the sky belongs to the wavelength fit, so the factor only moves sky
   light along the slit.

The module is pure numpy on 1-D arrays, one value per pixel, in any pixel order. The pixels may be those of one order or
of every order of a detector. The caller supplies the object model by subtracting it from the data and by leaving the
object pixels out of ``isUsable``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

# THE SLIT-POSITION BIN WIDTH IN ARCSEC, FIXED AT THE SOXS NIR PIXEL SCALE OF ABOUT 0.25 ARCSEC (VIS IS ABOUT 0.21)
BIN_WIDTH_ARCSEC = 0.25
# A SKYLINE PIXEL HAS A SKY MODEL OF AT LEAST THIS MANY TIMES ITS ERROR, SO ITS FRACTIONAL NOISE IS AT MOST
# 1/MIN_SKY_SNR.
# ON THE SOXS NIR STANDARD (DY-2054) THE RMS GAIN IN ORDERS 11 TO 13 RISES WITH THIS GATE UP TO ABOUT 30, AND A LOWER
# GATE LETS THE NOISY FAINT-ORDER LINES SET A CURVE
MIN_SKY_SNR = 30.0
# PIXELS WITH A LARGER FRACTIONAL RESIDUAL ARE CRs, OBJECT LIGHT OR A MISFIT LINE, NOT ILLUMINATION
MAX_FRACTIONAL_RESIDUAL = 0.1
# THE WAVELENGTH RANGE LEFT OUT AT EACH END OF THE ORDER, WHERE THE SKY SPLINE IS NOT CONSTRAINED
WAVELENGTH_MARGIN_FRACTION = 0.05
# A SLIT BIN NEEDS MORE THAN THIS MANY KEPT PIXELS TO BE USED, AND THE CURVE IS INTERPOLATED ACROSS A BIN WITH
# FEWER. THE IDL LIMIT IS 20 PIXELS OF ABOUT 45 SIGMA, A BIN MEAN GOOD TO 0.3%. PIXELS OF MIN_SKY_SNR SIGMA NEED
# 20 * (45 / 30) ** 2 = 45 PIXELS FOR THE SAME NOISE, AND THE SPARSE BINS AT THE SLIT ENDS SCATTER BY SEVERAL PERCENT.
# 20 PER BIN LET A NOISY CURVE FROM THE 41 PIXELS PER BIN OF NIR ORDER 10 MAKE ITS RESIDUAL WORSE (DY-2054)
MIN_PIXELS_PER_SLIT_BIN = 50
# THE CURVE NEEDS AT LEAST THIS MANY USABLE BINS TO BE INTERPOLATED
MIN_BINS = 2


@dataclass(frozen=True, eq=False)
class illumination_correction:
    """*the multiplicative slit-illumination factor of a set of pixels*

    **Key Arguments:**

    - ``factor`` -- the factor to multiply the sky model of each input pixel by. All ones if ``isApplied`` is False
    - ``isApplied`` -- True if enough skyline pixels qualified to measure a curve
    - ``pixelCount`` -- the number of skyline pixels that qualified
    - ``binCentres`` -- the slit position in arcsec of the centre of each bin used
    - ``binValues`` -- the mean fractional residual of each bin used, less the mean over all bins
    - ``binCounts`` -- the number of pixels in each bin used
    """

    factor: np.ndarray
    isApplied: bool
    pixelCount: int
    binCentres: np.ndarray
    binValues: np.ndarray
    binCounts: np.ndarray


def _unapplied(pixelTotal, pixelCount=0):
    """*the result for pixels that are left as they are*"""
    import numpy as np

    empty = np.array([])
    return illumination_correction(np.ones(pixelTotal), False, pixelCount, empty, empty, empty)


def _skyline_pixels(wavelength, skyModel, dataMinusObject, error, isUsable, minSkySnr):
    """*the mask of the bright skyline pixels and the fractional residual of every pixel*"""
    import numpy as np

    with np.errstate(divide="ignore", invalid="ignore"):
        fractional = (dataMinusObject - skyModel) / skyModel
    usableWavelength = wavelength[isUsable & np.isfinite(wavelength)]
    if usableWavelength.size == 0:
        return np.zeros(wavelength.shape, dtype=bool), fractional
    lowest, highest = usableWavelength.min(), usableWavelength.max()
    margin = WAVELENGTH_MARGIN_FRACTION * (highest - lowest)
    isInsideRange = (wavelength >= lowest + margin) & (wavelength <= highest - margin)
    with np.errstate(invalid="ignore"):
        isBright = (skyModel > 0) & (error > 0) & (skyModel >= minSkySnr * error)
        isSmall = np.abs(fractional) <= MAX_FRACTIONAL_RESIDUAL
    isSkyline = isUsable & isInsideRange & isBright & isSmall & np.isfinite(fractional)
    return isSkyline, fractional


def measure_slit_illumination(
    slitPosition,
    wavelength,
    skyModel,
    dataMinusObject,
    error,
    isUsable,
    minSkySnr=MIN_SKY_SNR,
    binWidth=BIN_WIDTH_ARCSEC,
):
    """*the multiplicative factor along the slit that the bright skylines show the sky model is missing*

    **Key Arguments:**

    - ``slitPosition`` -- the slit position of each pixel in arcsec
    - ``wavelength`` -- the wavelength of each pixel
    - ``skyModel`` -- the fitted sky model at each pixel
    - ``dataMinusObject`` -- the flux of each pixel minus the object model, 0 where there is none
    - ``error`` -- the flux error of each pixel
    - ``isUsable`` -- False for the pixels that are not sky: object, bad pixels, slit ends
    - ``minSkySnr`` -- a skyline pixel has a sky model of at least this many times its error
    - ``binWidth`` -- the width of a slit bin in arcsec

    **Return:**

    - ``correction`` -- an ``illumination_correction``. Its factor is all ones when fewer than ``MIN_BINS`` bins have
      enough skyline pixels, or when no pixel is usable
        - the curve is held flat outside the outer bins used, so the slit ends take the value of the outermost bin

    **Raises:**

    - ``ValueError`` -- if the input arrays do not all have the same length

    **Usage:**

    ```python
    correction = measure_slit_illumination(
        slitPosition, wavelength, skyModel, flux - wings, error, isSky
    )
    skyModel = skyModel * correction.factor
    ```
    """
    import numpy as np

    slitPosition = np.asarray(slitPosition, dtype=float)
    wavelength, skyModel, dataMinusObject, error, isUsable = (
        np.asarray(values) for values in (wavelength, skyModel, dataMinusObject, error, isUsable)
    )
    if not (
        slitPosition.shape
        == wavelength.shape
        == skyModel.shape
        == dataMinusObject.shape
        == error.shape
        == isUsable.shape
    ):
        raise ValueError(
            "the slit position, wavelength, sky model, data, error and usable arrays must have the same length"
        )
    isUsable = isUsable.astype(bool) & np.isfinite(slitPosition)
    pixelTotal = slitPosition.size
    if not isUsable.any():
        return _unapplied(pixelTotal)

    isSkyline, fractional = _skyline_pixels(
        wavelength.astype(float),
        skyModel.astype(float),
        dataMinusObject.astype(float),
        error.astype(float),
        isUsable,
        minSkySnr,
    )
    pixelCount = int(isSkyline.sum())
    slitStart = slitPosition[isUsable].min()
    slitBins = int(np.ceil((slitPosition[isUsable].max() - slitStart) / binWidth)) or 1
    # A PIXEL AT THE FAR EDGE OF A SPAN THAT IS A WHOLE NUMBER OF BINS BELONGS TO THE LAST BIN
    binIndex = np.floor((slitPosition[isSkyline] - slitStart) / binWidth).astype(int)
    binIndex = np.clip(binIndex, 0, slitBins - 1)
    counts = np.bincount(binIndex, minlength=slitBins)
    sums = np.bincount(binIndex, weights=fractional[isSkyline], minlength=slitBins)
    isBinUsed = counts > MIN_PIXELS_PER_SLIT_BIN
    if isBinUsed.sum() < MIN_BINS:
        return _unapplied(pixelTotal, pixelCount)

    centres = slitStart + (np.arange(counts.size) + 0.5) * binWidth
    values = sums[isBinUsed] / counts[isBinUsed]
    # THE MEAN LEVEL IS SET BY THE WAVELENGTH FIT, SO THE CURVE ONLY REDISTRIBUTES THE SKY ALONG THE SLIT
    values = values - np.average(values, weights=counts[isBinUsed])
    curve = np.interp(slitPosition, centres[isBinUsed], values)
    factor = np.where(np.isfinite(slitPosition), 1.0 + curve, 1.0)
    return illumination_correction(factor, True, pixelCount, centres[isBinUsed], values, counts[isBinUsed])
