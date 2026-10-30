"""Skyline-based multiplicative slit-illumination correction (DY-2054 prototype).

The stare sky model is a b-spline in wavelength only, so it assumes the sky is equally bright along the slit. A small
illumination gradient along the slit leaves a residual that scales with the sky, so it is largest at the bright
skylines. This module measures that gradient as the Kelson ``/fixillum`` step does (``longslit_reduce.pro``, lines 3709
to 3770), and returns the factor to multiply the sky model by:

1. The fractional residual of each pixel is (data minus object model minus sky model) over the sky model.
2. Only bright skyline pixels are kept: a sky model of at least ``MIN_SKY_SNR`` times the pixel error, a fractional
   residual of at most ``MAX_FRACTIONAL_RESIDUAL``, and a wavelength inside the central 90% of the order. The IDL
   uses a fixed 1000 counts. That is about 45 times the shot noise on X-Shooter NIR, so the same rule on SOXS, with
   its shorter exposures and its different gain, is a limit in units of the pixel error. In an order where more than
   ``MAX_CONTINUUM_FRACTION`` of the usable pixels pass that floor, the bright pixels are the sky continuum, not
   skylines (the thermal NIR order 10), and the floor is ``CONTINUUM_SKY_SNR``.
3. The inverse-variance weighted mean fractional residual is taken in slit-position bins of ``BIN_WIDTH_ARCSEC``,
   with a weight of (sky model / error) squared per pixel. A bin is used when the error of its mean is at most
   ``MAX_BIN_ERROR``, and at least ``MIN_BINS`` bins must be used.
4. The bin means are interpolated onto every pixel, and held flat outside the outer bins used. The factor is one plus
   that curve, less its weighted mean. The mean level of the sky belongs to the wavelength fit, so the factor only
   moves sky light along the slit.
5. Inside an optional fill zone (the object and the margin around it, where no bin is measured) the curve is a
   quadratic in slit position, inverse-variance fitted to the bins within ``FILL_FIT_HALF_WIDTH_ARCSEC`` of the slit
   centre, instead of a straight line between the bins at the zone edges. It is used only between the outermost of
   those bins, and with fewer than ``MIN_FILL_BINS`` of them the curve stays the straight line.

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
# A SKYLINE PIXEL HAS A SKY MODEL OF AT LEAST THIS MANY TIMES ITS ERROR. THE BIN MEANS ARE INVERSE-VARIANCE WEIGHTED,
# SO FAINTER PIXELS ADD LITTLE NOISE. ON THE SOXS NIR STANDARDS AND 3FGL (DY-2054) FLOORS OF 5, 10 AND 20 GIVE THE SAME
# CURVES IN ORDERS 11 TO 13, AND 10 ALSO CORRECTS ORDERS 13, 15 AND 16, WHICH A FLOOR OF 30 LEAVES OUT
MIN_SKY_SNR = 10.0
# WHEN MORE THAN THIS FRACTION OF AN ORDER'S USABLE PIXELS PASS MIN_SKY_SNR, THE ORDER IS CONTINUUM-DOMINATED AND
# ITS FLOOR IS CONTINUUM_SKY_SNR. NIR ORDER 10 (1.86 TO 2.03 MICRON, THERMAL) HAS 0.67 TO 1.0 OF ITS PIXELS ABOVE
# 10 SIGMA ON THREE FRAMES AND EVERY OTHER ORDER AT MOST 0.39. AT 10 SIGMA ITS CURVE MADE THE RESIDUAL WORSE (DY-2054)
MAX_CONTINUUM_FRACTION = 0.5
CONTINUUM_SKY_SNR = 30.0
# PIXELS WITH A LARGER FRACTIONAL RESIDUAL ARE CRs, OBJECT LIGHT OR A MISFIT LINE, NOT ILLUMINATION
MAX_FRACTIONAL_RESIDUAL = 0.1
# THE WAVELENGTH RANGE LEFT OUT AT EACH END OF THE ORDER, WHERE THE SKY SPLINE IS NOT CONSTRAINED
WAVELENGTH_MARGIN_FRACTION = 0.05
# A SLIT BIN IS USED WHEN THE ERROR OF ITS WEIGHTED MEAN IS AT MOST THIS, AND THE CURVE IS INTERPOLATED ACROSS A BIN
# WITH A LARGER ERROR. THE IDL LIMIT OF 20 PIXELS OF ABOUT 45 SIGMA GIVES 1 / 45 / SQRT(20) = 0.5%
MAX_BIN_ERROR = 0.005
# THE CURVE NEEDS AT LEAST THIS MANY USABLE BINS TO BE INTERPOLATED
MIN_BINS = 2
# THE FILL-ZONE QUADRATIC IS FITTED TO THE BINS WITHIN THIS DISTANCE OF THE SLIT CENTRE, AND NEEDS THIS MANY OF THEM
FILL_FIT_HALF_WIDTH_ARCSEC = 5.0
MIN_FILL_BINS = 4


@dataclass(frozen=True, eq=False)
class illumination_correction:
    """*the multiplicative slit-illumination factor of a set of pixels*

    **Key Arguments:**

    - ``factor`` -- the factor to multiply the sky model of each input pixel by. All ones if ``isApplied`` is False
    - ``isApplied`` -- True if enough skyline pixels qualified to measure a curve
    - ``pixelCount`` -- the number of skyline pixels that qualified
    - ``binCentres`` -- the slit position in arcsec of the centre of each bin used
    - ``binValues`` -- the weighted mean fractional residual of each bin used, less the weighted mean over all bins
    - ``binCounts`` -- the number of pixels in each bin used
    - ``binErrors`` -- the error of the weighted mean of each bin used
    - ``skySnrFloor`` -- the signal-to-noise floor of a skyline pixel in this order: ``MIN_SKY_SNR``, or
      ``CONTINUUM_SKY_SNR`` for a continuum-dominated order
    """

    factor: np.ndarray
    isApplied: bool
    pixelCount: int
    binCentres: np.ndarray
    binValues: np.ndarray
    binCounts: np.ndarray
    binErrors: np.ndarray
    skySnrFloor: float


def _unapplied(pixelTotal, pixelCount=0, skySnrFloor=MIN_SKY_SNR):
    """*the result for pixels that are left as they are*"""
    import numpy as np

    empty = np.array([])
    return illumination_correction(np.ones(pixelTotal), False, pixelCount, empty, empty, empty, empty, skySnrFloor)


def _sky_snr_floor(skyModel, error, isUsable, minSkySnr):
    """*the skyline floor of an order: minSkySnr, or CONTINUUM_SKY_SNR when most usable pixels pass minSkySnr*"""
    import numpy as np

    with np.errstate(invalid="ignore"):
        isBright = (skyModel > 0) & (error > 0) & (skyModel >= minSkySnr * error)
    brightFraction = float(np.mean(isBright[isUsable]))
    return max(minSkySnr, CONTINUUM_SKY_SNR) if brightFraction > MAX_CONTINUUM_FRACTION else minSkySnr


def _fill_with_quadratic(curve, slitPosition, isFillZone, centres, values, errors):
    """*the curve with its fill-zone values replaced by a quadratic fitted to the bins near the slit centre*

    The quadratic is used only between the outermost bins it is fitted to, so a fill zone at a slit end keeps the
    interpolated curve, held flat, beyond them. It does not have to meet the interpolated curve at the zone edges, so
    the sky model can step there.
    """
    import numpy as np

    isNearCentre = np.abs(centres) < FILL_FIT_HALF_WIDTH_ARCSEC
    if isNearCentre.sum() < MIN_FILL_BINS:
        return curve
    fitCentres = centres[isNearCentre]
    coefficients = np.polyfit(fitCentres, values[isNearCentre], 2, w=1.0 / errors[isNearCentre])
    with np.errstate(invalid="ignore"):
        isFilled = isFillZone & (slitPosition >= fitCentres.min()) & (slitPosition <= fitCentres.max())
    return np.where(isFilled, np.polyval(coefficients, slitPosition), curve)


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
    isFillZone=None,
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
    - ``isFillZone`` -- True for the pixels whose factor comes from a quadratic fitted to the bins (the object and the
      margin around it). None, the default, interpolates the bins everywhere
    - ``minSkySnr`` -- a skyline pixel has a sky model of at least this many times its error
    - ``binWidth`` -- the width of a slit bin in arcsec

    **Return:**

    - ``correction`` -- an ``illumination_correction``. Its factor is all ones when fewer than ``MIN_BINS`` bins have
      a mean with an error of at most ``MAX_BIN_ERROR``, or when no pixel is usable
        - the curve is held flat outside the outer bins used, so the slit ends take the value of the outermost bin

    **Raises:**

    - ``ValueError`` -- if the input arrays do not all have the same length

    **Usage:**

    ```python
    correction = measure_slit_illumination(
        slitPosition, wavelength, skyModel, flux - wings, error, isSky, isFillZone=isNearObject
    )
    skyModel = skyModel * correction.factor
    ```
    """
    import numpy as np

    slitPosition = np.asarray(slitPosition, dtype=float)
    wavelength, skyModel, dataMinusObject, error, isUsable = (
        np.asarray(values) for values in (wavelength, skyModel, dataMinusObject, error, isUsable)
    )
    isFillZone = np.zeros(slitPosition.shape, dtype=bool) if isFillZone is None else np.asarray(isFillZone, dtype=bool)
    if not (
        slitPosition.shape
        == wavelength.shape
        == skyModel.shape
        == dataMinusObject.shape
        == error.shape
        == isUsable.shape
        == isFillZone.shape
    ):
        raise ValueError(
            "the slit position, wavelength, sky model, data, error, usable and fill-zone arrays must have the same "
            "length"
        )
    isUsable = isUsable.astype(bool) & np.isfinite(slitPosition)
    pixelTotal = slitPosition.size
    if not isUsable.any():
        return _unapplied(pixelTotal)

    skyModel, error = skyModel.astype(float), error.astype(float)
    floor = _sky_snr_floor(skyModel, error, isUsable, minSkySnr)
    isSkyline, fractional = _skyline_pixels(
        wavelength.astype(float), skyModel, dataMinusObject.astype(float), error, isUsable, floor
    )
    slitStart = slitPosition[isUsable].min()
    slitBins = int(np.ceil((slitPosition[isUsable].max() - slitStart) / binWidth)) or 1
    # A PIXEL AT THE FAR EDGE OF A SPAN THAT IS A WHOLE NUMBER OF BINS BELONGS TO THE LAST BIN
    binIndex = np.floor((slitPosition[isSkyline] - slitStart) / binWidth).astype(int)
    binIndex = np.clip(binIndex, 0, slitBins - 1)
    with np.errstate(over="ignore"):
        weights = (skyModel[isSkyline] / error[isSkyline]) ** 2
    # A WEIGHT THAT OVERFLOWS (A NEAR-ZERO ERROR) WOULD TURN ITS BIN MEAN, AND THE CURVE AROUND IT, INTO NAN
    isFiniteWeight = np.isfinite(weights)
    # DROP THOSE PIXELS FROM THE SKYLINE MASK ITSELF, SO fractional[isSkyline] STAYS ALIGNED WITH THE WEIGHTS
    isSkyline[isSkyline] = isFiniteWeight
    weights = weights[isFiniteWeight]
    pixelCount = int(isSkyline.sum())
    binIndex = binIndex[isFiniteWeight]
    counts = np.bincount(binIndex, minlength=slitBins)
    weightSums = np.bincount(binIndex, weights=weights, minlength=slitBins)
    weightedSums = np.bincount(binIndex, weights=weights * fractional[isSkyline], minlength=slitBins)
    with np.errstate(divide="ignore"):
        binErrors = 1.0 / np.sqrt(weightSums)
    isBinUsed = binErrors <= MAX_BIN_ERROR
    if isBinUsed.sum() < MIN_BINS:
        return _unapplied(pixelTotal, pixelCount, floor)

    centres = (slitStart + (np.arange(counts.size) + 0.5) * binWidth)[isBinUsed]
    values = weightedSums[isBinUsed] / weightSums[isBinUsed]
    # THE MEAN LEVEL IS SET BY THE WAVELENGTH FIT, SO THE CURVE ONLY REDISTRIBUTES THE SKY ALONG THE SLIT
    values = values - np.average(values, weights=weightSums[isBinUsed])
    curve = np.interp(slitPosition, centres, values)
    curve = _fill_with_quadratic(curve, slitPosition, isFillZone, centres, values, binErrors[isBinUsed])
    factor = np.where(np.isfinite(slitPosition), 1.0 + curve, 1.0)
    return illumination_correction(
        factor, True, pixelCount, centres, values, counts[isBinUsed], binErrors[isBinUsed], floor
    )
