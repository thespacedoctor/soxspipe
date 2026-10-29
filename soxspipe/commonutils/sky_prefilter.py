"""*the Kelson-style two-pass offset-bin pre-filter of the sky pixels (DY-2055)*

The filter is a port of the pre-fit outlier rejection of Kelson's ``longslit_reduce.pro``. The pixels of an order are
sorted by wavelength and cut into bins of equal wavelength width. Each bin flags the pixels that lie more than a
sigma limit from the bin median, in units of the robust sigma of the bin. The cut is made twice, the second time with
the bin edges shifted by half a bin, and only a pixel flagged in both passes is rejected. A narrow skyline that sits on
a bin edge is split between two bins in one pass, and so looks like an outlier there, but it sits in the middle of a
bin in the other pass.

The differences from the IDL are these:

- the IDL tests only the pixels more than 3 sigma from the object position, which soxspipe does not know at this
  stage, so the caller passes every pixel that is not already clipped
- the IDL drops the pixels below the shifted first edge from the second pass, and here they form a half-width first bin
- the IDL flags a pixel with ``>=`` and so flags every pixel of a bin whose sigma is zero, and here such a bin flags
  nothing
- the IDL tests a bin of any size, and here a bin with fewer than ``MIN_BIN_SAMPLES`` pixels is not tested
"""

# THE SIGMA LIMIT OF THE IDL FILTER
OFFSET_BIN_SIGMA = 3.0
# THE BIN WIDTH, IN DISPERSION PIXELS OF THE ORDER. THE IDL USES 5 TIMES THE MINIMUM BREAK-POINT SPACING, AND HERE THAT
# SPACING IS ONE DISPERSION PIXEL
OFFSET_BIN_WIDTH_PIXELS = 5.0
# A BIN WITH FEWER FINITE PIXELS THAN THIS HAS NO USABLE ROBUST SIGMA AND IS NOT TESTED
MIN_BIN_SAMPLES = 10

# THE ROBUST SIGMA OF GODDARD'S ROBUST_SIGMA.PRO, WHICH THE IDL CALLS. A MAD IS DIVIDED BY THIS TO GIVE A SIGMA
MAD_NORMALISATION = 0.6745
# A MEAN ABSOLUTE DEVIATION IS DIVIDED BY THIS TO GIVE A SIGMA WHEN THE MAD IS ZERO
MEAN_DEVIATION_NORMALISATION = 0.8
# THE BIWEIGHT CUTS OFF AT THIS MANY MAD-SIGMAS FROM THE MEDIAN
BIWEIGHT_CUTOFF = 6.0
# A SCALE BELOW THIS IS ZERO
SCALE_EPSILON = 1e-20
# THE BIWEIGHT NEEDS AT LEAST THIS MANY VALUES INSIDE ITS CUTOFF
BIWEIGHT_MIN_INLIERS = 3


def robust_sigma(values):
    """*the biweight robust sigma of a set of values, as in the IDL ``robust_sigma``*

    **Key Arguments:**

    - ``values`` -- array of values. Values that are not finite are ignored

    **Return:**

    - ``sigma`` -- the biweight midvariance scale; 0 when the values are all equal, or fewer than
      ``BIWEIGHT_MIN_INLIERS`` lie inside the biweight cutoff
    """
    import numpy as np

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0

    deviation = values - np.median(values)
    scale = np.median(np.abs(deviation)) / MAD_NORMALISATION
    if scale < SCALE_EPSILON:
        scale = np.mean(np.abs(deviation)) / MEAN_DEVIATION_NORMALISATION
        if scale < SCALE_EPSILON:
            return 0.0

    u2 = (deviation / (BIWEIGHT_CUTOFF * scale)) ** 2
    isInside = u2 <= 1.0
    if np.count_nonzero(isInside) < BIWEIGHT_MIN_INLIERS:
        return 0.0
    numerator = np.sum(deviation[isInside] ** 2 * (1.0 - u2[isInside]) ** 4)
    denominator = np.sum((1.0 - u2[isInside]) * (1.0 - 5.0 * u2[isInside]))
    variance = values.size * numerator / (denominator * (denominator - 1.0))
    return float(np.sqrt(variance)) if variance > 0 else 0.0


def order_dispersion(wavelength, x, y):
    """*the wavelength step of one dispersion pixel of an order*

    The step is the wavelength span of the order over the number of pixel intervals it covers along its longer
    detector axis. That axis is the dispersion axis, so the rule needs no knowledge of the detector orientation or of
    the binning: the pixel indices of a binned frame are binned pixels, and the step is per binned pixel.

    **Key Arguments:**

    - ``wavelength`` -- the wavelength of each pixel of the order
    - ``x``, ``y`` -- the detector column and row of each pixel

    **Return:**

    - ``step`` -- the wavelength step per dispersion pixel; 0 when the order spans fewer than two pixels or a zero
      wavelength range
    """
    import numpy as np

    wavelength = np.asarray(wavelength, dtype=float)
    isFinite = np.isfinite(wavelength)
    if not isFinite.any():
        return 0.0
    span = float(wavelength[isFinite].max() - wavelength[isFinite].min())
    pixelCount = max(np.unique(np.asarray(x)[isFinite]).size, np.unique(np.asarray(y)[isFinite]).size)
    if pixelCount < 2:
        return 0.0
    return span / (pixelCount - 1)


def _flag_bins(wavelength, flux, isUsable, binWidth, offset, sigma):
    """*flag the pixels of the usable set that lie more than ``sigma`` robust sigmas from their bin median*

    **Key Arguments:**

    - ``wavelength`` -- wavelength of each pixel
    - ``flux`` -- flux of each pixel
    - ``isUsable`` -- boolean array, True for a pixel with a finite wavelength and flux
    - ``binWidth`` -- the bin width, in wavelength units
    - ``offset`` -- the shift of the bin edges from the shortest usable wavelength
    - ``sigma`` -- the sigma limit

    **Return:**

    - ``flagged`` -- boolean array, True for a flagged pixel
    """
    import numpy as np

    flagged = np.zeros(wavelength.size, dtype=bool)
    usable = np.flatnonzero(isUsable)
    binIndex = np.floor((wavelength[usable] - (wavelength[usable].min() + offset)) / binWidth).astype(np.int64)
    byBin = np.argsort(binIndex, kind="stable")
    firstOfBin = np.flatnonzero(np.diff(binIndex[byBin], prepend=binIndex[byBin][0] - 1))
    for members in np.split(usable[byBin], firstOfBin[1:]):
        if members.size < MIN_BIN_SAMPLES:
            continue
        binFlux = flux[members]
        scatter = robust_sigma(binFlux)
        if scatter > 0:
            flagged[members] = np.abs(binFlux - np.median(binFlux)) > sigma * scatter
    return flagged


def offset_bin_pass_masks(wavelength, flux, binWidth, sigma=OFFSET_BIN_SIGMA):
    """*the pixels flagged by each of the two offset-bin passes*

    **Key Arguments:**

    - ``wavelength`` -- wavelength of each pixel, in any row order
    - ``flux`` -- flux of each pixel
    - ``binWidth`` -- the bin width, in wavelength units
    - ``sigma`` -- the sigma limit. Default *OFFSET_BIN_SIGMA*

    **Return:**

    - ``firstPass`` -- boolean array in the input row order, True for a pixel the first pass flags
    - ``secondPass`` -- the same for the pass whose bin edges are shifted by half a bin
        - a pixel with a non-finite wavelength or flux is never flagged, and does not set a bin statistic
        - nothing is flagged when ``binWidth`` is not a positive number
    """
    import numpy as np

    wavelength = np.asarray(wavelength, dtype=float)
    flux = np.asarray(flux, dtype=float)
    isUsable = np.isfinite(wavelength) & np.isfinite(flux)
    if not (np.isfinite(binWidth) and binWidth > 0) or not isUsable.any():
        return np.zeros(flux.size, dtype=bool), np.zeros(flux.size, dtype=bool)
    return (
        _flag_bins(wavelength, flux, isUsable, binWidth, 0.0, sigma),
        _flag_bins(wavelength, flux, isUsable, binWidth, 0.5 * binWidth, sigma),
    )


def two_pass_offset_bin_outliers(wavelength, flux, binWidth, sigma=OFFSET_BIN_SIGMA):
    """*the pixels flagged in both offset-bin passes*

    **Key Arguments:**

    - ``wavelength`` -- wavelength of each pixel, in any row order
    - ``flux`` -- flux of each pixel
    - ``binWidth`` -- the bin width, in wavelength units
    - ``sigma`` -- the sigma limit. Default *OFFSET_BIN_SIGMA*

    **Return:**

    - ``isRejected`` -- boolean array in the input row order, True for a pixel both passes flag
    """
    firstPass, secondPass = offset_bin_pass_masks(wavelength, flux, binWidth, sigma)
    return firstPass & secondPass
