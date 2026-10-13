"""Skyline slit profiles for the slit-illumination measurement.

This module holds the code that turns a flat-fielded frame and a skyline list into one raw-like sky profile along the
slit for each order, plus the slit-binning helpers it shares with the flat-edge fit in ``slit_illumination``. The
public function ``skyline_slit_profiles`` is also importable from ``soxspipe.commonutils.slit_illumination``.
"""

import warnings

import numpy as np
import pandas as pd

PROFILE_BIN_ARCSEC = 0.1
# THE INTERIOR USED TO NORMALISE A PROFILE IS THIS FRACTION OF THE ORDER'S SLIT RANGE, CENTRED ON THE RANGE
INTERIOR_HALF_FRACTION = 0.5
# A FLAT SLIT BIN NEEDS THIS MANY PIXELS TO GIVE A MEDIAN
MIN_FLAT_PIXELS_PER_BIN = 15
MIN_FLAT_BINS = 20
# SKYLINE WINDOWS, IN UNITS OF THE LOCAL WAVELENGTH STEP PER PIXEL. THE CORE HOLDS THE LINE, THE SIDE BANDS HOLD ONLY
# CONTINUUM. VIS LINES ARE ABOUT 3 PIXELS FWHM, SO THE SIDE BANDS START WELL OUTSIDE THE WINGS
CORE_HALF_WIDTH_PIXELS = 2.5
SIDE_BAND_INNER_PIXELS = 6.0
SIDE_BAND_OUTER_PIXELS = 12.0
MIN_CORE_PIXELS = 1
MIN_SIDE_PIXELS = 2
# A LINE IS USED ONLY IF ITS INTERIOR AMPLITUDE EXCEEDS THIS MANY TIMES ITS TYPICAL ERROR IN ONE INTERIOR SLIT BIN
MIN_LINE_SNR = 5.0
MIN_LINES_PER_BIN = 3
# THE LOCAL WAVELENGTH STEP IS THE MEDIAN OVER THIS MANY PIXELS NEAREST A LINE IN WAVELENGTH
LOCAL_STEP_PIXELS = 1000
# THE STANDARD ERROR OF A MEDIAN IS THIS MANY TIMES THE STANDARD ERROR OF THE MEAN
MEDIAN_ERROR_FACTOR = np.sqrt(np.pi / 2)
MAD_TO_SIGMA = 1.4826
PROFILE_COLUMNS = ["order", "slit", "profile", "profileErr", "nLines"]

NO_FLAT_LEVEL_REASON = "flat level is not positive"
TOO_FEW_FLAT_BINS_REASON = "too few populated flat bins"


def in_order_pixels(orderMap, slit, wavelength=None):
    """*the pixels that belong to a whole-numbered order and have a finite slit position (and wavelength)*

    Binned maps are block averages, so a pixel on an order boundary can carry a fractional order number. Those
    pixels belong to no order and are left out.
    """
    isInOrder = np.isfinite(orderMap) & np.isfinite(slit) & (orderMap > 0)
    isInOrder &= orderMap == np.round(np.where(np.isfinite(orderMap), orderMap, 0.0))
    if wavelength is not None:
        isInOrder &= np.isfinite(wavelength) & (wavelength > 0)
    return isInOrder


def _bin_geometry(orderSlit, binArcsec):
    """*the integer index of the first bin, the bin centres and the interior flags for an order's slit range*"""
    firstIndex = int(np.floor(np.min(orderSlit) / binArcsec))
    lastIndex = int(np.floor(np.max(orderSlit) / binArcsec))
    centres = (np.arange(firstIndex, lastIndex + 1) + 0.5) * binArcsec
    middle = 0.5 * (centres[0] + centres[-1])
    halfRange = 0.5 * (centres[-1] - centres[0])
    isInterior = np.abs(centres - middle) <= INTERIOR_HALF_FRACTION * halfRange
    return firstIndex, centres, isInterior


def _binned_medians(binIndex, values, nBins, minPixels):
    """*the median of ``values`` in each of ``nBins`` bins, NaN where a bin holds fewer than ``minPixels`` values*"""
    medians = np.full(nBins, np.nan)
    grouped = pd.Series(values).groupby(binIndex)
    counts = grouped.size()
    found = grouped.median()
    isEnough = counts.to_numpy() >= minPixels
    medians[found.index.to_numpy()[isEnough]] = found.to_numpy()[isEnough]
    return medians


def normalised_flat_profile(flatData, slit, isOrder, binArcsec):
    """*the order's flat in slit bins, normalised to its interior median*

    **Return:**

    - ``centres`` -- the slit bin centres in arcsec
    - ``profile`` -- the normalised flat profile, NaN in sparse bins and everywhere if the profile is unusable
    - ``reason`` -- why the profile is unusable, or an empty string
    """
    firstIndex, centres, isInterior = _bin_geometry(slit[isOrder], binArcsec)
    binIndex = np.floor(slit[isOrder] / binArcsec).astype(int) - firstIndex
    values = flatData[isOrder]
    isFinite = np.isfinite(values)
    profile = _binned_medians(binIndex[isFinite], values[isFinite], len(centres), MIN_FLAT_PIXELS_PER_BIN)
    if np.isfinite(profile).sum() < MIN_FLAT_BINS:
        return centres, np.full(len(centres), np.nan), TOO_FEW_FLAT_BINS_REASON
    level = np.nanmedian(profile[isInterior]) if np.isfinite(profile[isInterior]).any() else np.nan
    if not (np.isfinite(level) and level > 0):
        return centres, np.full(len(centres), np.nan), NO_FLAT_LEVEL_REASON
    return centres, profile / level, ""


def _nan_median(values, axis):
    """*the median along an axis ignoring NaN, NaN (without a warning) where everything is NaN*"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(values, axis=axis)


def _wavelength_step(wavelength):
    """*the wavelength change per pixel at every pixel, from the gradient of the wavelength image*"""
    gradientY, gradientX = np.gradient(wavelength)
    return np.hypot(gradientY, gradientX)


def _slice_bin_stats(binIndex, values, errs2, nBins):
    """*per slit bin: the median of ``values``, the mean of ``errs2`` and the pixel count*"""
    grouped = pd.DataFrame({"bin": binIndex, "value": values, "err2": errs2}).groupby("bin")
    stats = grouped.agg(median=("value", "median"), err2=("err2", "mean"), count=("value", "size"))
    median, err2, count = np.full(nBins, np.nan), np.full(nBins, np.nan), np.zeros(nBins, dtype=int)
    index = stats.index.to_numpy()
    median[index], err2[index], count[index] = stats["median"], stats["err2"], stats["count"]
    return median, err2, count


class _order_pixels:
    """The usable pixels of one order sorted by wavelength, with their slit bin and local wavelength step."""

    def __init__(self, wavelength, binIndex, values, errs2, step):
        sortIndex = np.argsort(wavelength, kind="stable")
        self.wavelength = wavelength[sortIndex]
        self.binIndex = binIndex[sortIndex]
        self.values = values[sortIndex]
        self.errs2 = errs2[sortIndex]
        self.step = step[sortIndex]

    def local_step(self, lineWavelength):
        """*the median wavelength step per pixel of the pixels nearest a line in wavelength*"""
        centre = np.searchsorted(self.wavelength, lineWavelength)
        first = max(centre - LOCAL_STEP_PIXELS // 2, 0)
        return _nan_median(self.step[first : first + LOCAL_STEP_PIXELS], axis=None)

    def window(self, lowWavelength, highWavelength):
        """*the slice of pixels with wavelength in ``[lowWavelength, highWavelength)``*"""
        first, last = np.searchsorted(self.wavelength, [lowWavelength, highWavelength])
        return slice(first, last)


def _line_difference(pixels, lineWavelength, step, nBins):
    """*per slit bin, the core median minus the side-band median of one line, with its error and NaN where sparse*"""
    core = pixels.window(lineWavelength - CORE_HALF_WIDTH_PIXELS * step, lineWavelength + CORE_HALF_WIDTH_PIXELS * step)
    inner, outer = SIDE_BAND_INNER_PIXELS * step, SIDE_BAND_OUTER_PIXELS * step
    below = pixels.window(lineWavelength - outer, lineWavelength - inner)
    above = pixels.window(lineWavelength + inner, lineWavelength + outer)
    side = np.concatenate([np.arange(below.start, below.stop), np.arange(above.start, above.stop)])
    coreStats = _slice_bin_stats(pixels.binIndex[core], pixels.values[core], pixels.errs2[core], nBins)
    sideStats = _slice_bin_stats(pixels.binIndex[side], pixels.values[side], pixels.errs2[side], nBins)
    coreMedian, coreErr2, coreCount = coreStats
    sideMedian, sideErr2, sideCount = sideStats
    isDense = (coreCount >= MIN_CORE_PIXELS) & (sideCount >= MIN_SIDE_PIXELS)
    with np.errstate(invalid="ignore", divide="ignore"):
        variance = MEDIAN_ERROR_FACTOR**2 * (coreErr2 / coreCount + sideErr2 / sideCount)
    difference = np.where(isDense, coreMedian - sideMedian, np.nan)
    return difference, np.where(isDense, np.sqrt(variance), np.nan)


def _clean_lines(pixels, lineWavelengths):
    """*the lines inside the order whose windows hold no neighbouring line, each with its local wavelength step*"""
    lines = np.sort(lineWavelengths[np.isfinite(lineWavelengths)])
    clean = []
    for lineWavelength in lines:
        step = pixels.local_step(lineWavelength)
        if not np.isfinite(step) or step <= 0:
            continue
        reach = SIDE_BAND_OUTER_PIXELS * step
        if lineWavelength - reach < pixels.wavelength[0] or lineWavelength + reach > pixels.wavelength[-1]:
            continue
        neighbours = np.abs(lines - lineWavelength) < reach + CORE_HALF_WIDTH_PIXELS * step
        if neighbours.sum() > 1:
            continue
        clean.append((lineWavelength, step))
    return clean


def _combine_lines(differences, errors, isInterior):
    """*normalise each line to its interior amplitude and combine the lines per slit bin*

    **Return:**

    - ``profile`` -- the median normalised line profile per bin
    - ``profileErr`` -- the larger of the line-to-line scatter and the propagated ERRS error of the median
    - ``nLines`` -- the number of lines in each bin

    or *None* if no line is bright enough.
    """
    amplitude = _nan_median(differences[:, isInterior], axis=1)
    amplitudeErr = _nan_median(errors[:, isInterior], axis=1)
    with np.errstate(invalid="ignore"):
        isGood = np.isfinite(amplitude) & (amplitude > 0) & (amplitude >= MIN_LINE_SNR * amplitudeErr)
    if not isGood.any():
        return None
    normalised = differences[isGood] / amplitude[isGood, None]
    normalisedErr = errors[isGood] / amplitude[isGood, None]
    nLines = np.isfinite(normalised).sum(axis=0)
    profile = _nan_median(normalised, axis=0)
    deviation = _nan_median(np.abs(normalised - profile), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        scatterErr = MEDIAN_ERROR_FACTOR * MAD_TO_SIGMA * deviation / np.sqrt(nLines)
        errsErr = MEDIAN_ERROR_FACTOR / np.sqrt(np.nansum(1.0 / normalisedErr**2, axis=0))
    profile = np.where(nLines >= MIN_LINES_PER_BIN, profile, np.nan)
    return profile, np.fmax(scatterErr, errsErr), nLines


def _order_skyline_profile(pixels, nBins, lineWavelengths, isInterior):
    """*the combined skyline slit profile of one order, or None with the reason there is none*

    **Return:**

    - ``combined`` -- the ``(profile, profileErr, nLines)`` of ``_combine_lines``, or *None*
    - ``note`` -- why there is no profile, or an empty string
    """
    firstWavelength, lastWavelength = pixels.wavelength[0], pixels.wavelength[-1]
    if not ((lineWavelengths >= firstWavelength) & (lineWavelengths <= lastWavelength)).any():
        return None, f"the skyline table has no lines between {firstWavelength:.0f} and {lastWavelength:.0f} nm"
    clean = _clean_lines(pixels, lineWavelengths)
    if not clean:
        return None, "every skyline in the order is blended with a neighbour or too close to the end of the order"
    differences = np.full((len(clean), nBins), np.nan)
    errors = np.full((len(clean), nBins), np.nan)
    for index, (lineWavelength, step) in enumerate(clean):
        differences[index], errors[index] = _line_difference(pixels, lineWavelength, step, nBins)
    combined = _combine_lines(differences, errors, isInterior)
    if combined is None:
        return None, f"every skyline in the order is too faint (under {MIN_LINE_SNR:.0f} sigma per slit bin)"
    return combined, ""


def skyline_slit_profiles(
    data,
    errs,
    isUsable,
    flatData,
    wavelength,
    slit,
    orderMap,
    lineWavelengths,
    binArcsec=PROFILE_BIN_ARCSEC,
):
    """*the raw-like sky slit profile of every order, from skyline cores minus their side bands*

    For each order, skyline and slit bin, the median of the side-band pixels is subtracted from the median of the
    core pixels. This removes the object and the sky continuum, so no object mask is needed. Each line is normalised
    to its interior amplitude, and the median over lines gives one profile per order. The frame is flat-fielded, so
    the profile is multiplied by the order's normalised flat profile to give the profile before flat-fielding. That
    flat profile uses the pixels with a positive finite wavelength, whereas ``flat_slit_edges`` does not require a
    wavelength, so the two pixel sets can differ at the order ends. The
    core and side-band windows are in units of the local wavelength step per pixel (``CORE_HALF_WIDTH_PIXELS``,
    ``SIDE_BAND_INNER_PIXELS`` and ``SIDE_BAND_OUTER_PIXELS``). A line is skipped when another line lies within
    its side bands, or when its side bands leave the order's wavelength range.

    **Key Arguments:**

    - ``data`` -- the flat-fielded, detrended frame
    - ``errs`` -- the frame's 1-sigma uncertainty
    - ``isUsable`` -- boolean image, True where ``data`` and ``errs`` can be used
    - ``flatData`` -- the master flat's flux image
    - ``wavelength`` -- the 2D map's WAVELENGTH image in nm
    - ``slit`` -- the 2D map's SLIT image in arcsec
    - ``orderMap`` -- the 2D map's ORDER image
    - ``lineWavelengths`` -- the skyline wavelengths in nm
    - ``binArcsec`` -- the slit bin width in arcsec. Default *PROFILE_BIN_ARCSEC*

    **Return:**

    - ``profiles`` -- dataframe with one row per order and populated slit bin and the columns ``order``, ``slit``,
      ``profile``, ``profileErr`` and ``nLines``. Empty if no order has a usable profile
    - ``notes`` -- ``{order: reason}`` for each order that has no profile rows

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination import skyline_slit_profiles
    profiles, notes = skyline_slit_profiles(
        frame.data, frame.uncertainty.array, isUsable, flat.data, wavelength, slit, orderMap, lineWavelengths
    )
    ```
    """
    isInOrder = in_order_pixels(orderMap, slit, wavelength)
    step = _wavelength_step(wavelength)
    lineWavelengths = np.asarray(lineWavelengths, dtype=float)
    lineWavelengths = lineWavelengths[np.isfinite(lineWavelengths)]
    pieces, notes = [], {}
    for order in np.unique(orderMap[isInOrder]):
        isOrder = isInOrder & (orderMap == order)
        centres, flatProfile, reason = normalised_flat_profile(flatData, slit, isOrder, binArcsec)
        isPixel = isOrder & isUsable & np.isfinite(data) & np.isfinite(errs)
        if reason or not isPixel.any():
            notes[int(order)] = (
                f"the flat profile is unusable: {reason}" if reason else "the order has no usable pixels"
            )
            continue
        firstIndex, _, isInterior = _bin_geometry(slit[isOrder], binArcsec)
        binIndex = np.floor(slit[isPixel] / binArcsec).astype(int) - firstIndex
        pixels = _order_pixels(wavelength[isPixel], binIndex, data[isPixel], errs[isPixel] ** 2, step[isPixel])
        combined, note = _order_skyline_profile(pixels, len(centres), lineWavelengths, isInterior)
        if combined is None:
            notes[int(order)] = note
            continue
        profile, profileErr, nLines = combined
        piece = pd.DataFrame(
            {
                "order": int(order),
                "slit": centres,
                "profile": profile * flatProfile,
                "profileErr": profileErr * flatProfile,
                "nLines": nLines,
            }
        )
        isValid = np.isfinite(piece["profile"]) & np.isfinite(piece["profileErr"]) & (piece["profileErr"] > 0)
        if not isValid.any():
            notes[int(order)] = f"no slit bin holds at least {MIN_LINES_PER_BIN} usable skylines"
            continue
        pieces.append(piece.loc[isValid])
    if not pieces:
        return pd.DataFrame(columns=PROFILE_COLUMNS), notes
    return pd.concat(pieces, ignore_index=True)[PROFILE_COLUMNS], notes
