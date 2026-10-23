"""Skyline and sky-continuum slit profiles for the slit-illumination measurement.

This module holds the code that turns a flat-fielded frame and a skyline list into one raw-like sky profile along the
slit for each order, plus the slit-binning helpers it shares with the flat-edge fit in ``slit_illumination``. For an
order the skyline list has no lines in, ``continuum_slit_profiles`` builds the profile from the sky continuum. The
public functions ``skyline_slit_profiles`` and ``continuum_slit_profiles`` are also importable from
``soxspipe.commonutils.slit_illumination``.
"""

import warnings
from typing import NamedTuple

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
# NLINES IS THE NUMBER OF SKYLINES BEHIND A BIN OF A SKYLINE PROFILE, AND THE NUMBER OF WAVELENGTH CHUNKS BEHIND A BIN
# OF A CONTINUUM PROFILE. THE OFFSET FIT DOES NOT USE IT
PROFILE_COLUMNS = ["order", "slit", "profile", "profileErr", "nLines"]

# THE CONTINUUM-EDGE PROFILE (ORDERS WITH NO SKYLINES) SPLITS AN ORDER INTO THIS MANY EQUAL-COUNT WAVELENGTH CHUNKS
CONTINUUM_CHUNKS = 20
# A SLIT BIN OF THE COMBINED CONTINUUM PROFILE NEEDS THIS MANY CHUNKS
CONTINUUM_MIN_CHUNKS_PER_BIN = 5
# BINS THIS FAR (BOTH SIDES) FROM THE OBJECT ARE LEFT OUT OF THE CONTINUUM PROFILE
CONTINUUM_OBJECT_MASK_ARCSEC = 2.5
# A SLIT BIN OF ONE CHUNK NEEDS THIS MANY PIXELS TO GIVE A MEDIAN
MIN_CHUNK_PIXELS_PER_BIN = 8
# THE OBJECT SEARCH KEEPS THIS FAR INSIDE THE FLAT EDGES, SO THAT A SLIT CLIFF IS NEVER TAKEN FOR THE OBJECT
OBJECT_SEARCH_MARGIN_ARCSEC = 1.0
MIN_OBJECT_SEARCH_BINS = 10
# THE OBJECT PEAK IS PICKED FROM THE PROFILE SMOOTHED OVER THIS MANY BINS
OBJECT_SMOOTHING_BINS = 5
OBJECT_SIGNIFICANCE = 5.0

NO_SKYLINES_NOTE_PREFIX = "the skyline table has no lines between"
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


def has_no_skylines_note(note):
    """*True if a ``skyline_slit_profiles`` note says the skyline table has no lines in the order*"""
    return isinstance(note, str) and note.startswith(NO_SKYLINES_NOTE_PREFIX)


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


class _order_setup(NamedTuple):
    """The per-order inputs that the skyline and the continuum profile builders share."""

    centres: np.ndarray
    flatProfile: np.ndarray
    isInterior: np.ndarray
    isPixel: np.ndarray
    binIndex: np.ndarray


def _prepare_order(data, errs, isUsable, flatData, slit, isOrder, binArcsec):
    """*the shared setup of one order, or None with the reason the order cannot be profiled*

    **Return:**

    - ``setup`` -- an ``_order_setup``, or *None*
    - ``note`` -- why there is no setup, or an empty string
    """
    centres, flatProfile, reason = normalised_flat_profile(flatData, slit, isOrder, binArcsec)
    isPixel = isOrder & isUsable & np.isfinite(data) & np.isfinite(errs)
    if reason or not isPixel.any():
        return None, f"the flat profile is unusable: {reason}" if reason else "the order has no usable pixels"
    firstIndex, _, isInterior = _bin_geometry(slit[isOrder], binArcsec)
    binIndex = np.floor(slit[isPixel] / binArcsec).astype(int) - firstIndex
    return _order_setup(centres, flatProfile, isInterior, isPixel, binIndex), ""


def _profile_piece(order, setup, profile, profileErr, counts, isKept=None):
    """*the profile rows (``PROFILE_COLUMNS``) of the finite bins of one combined profile, times the flat profile*

    ``isKept`` optionally limits the rows to some bins, as the object mask does.
    """
    piece = pd.DataFrame(
        {
            "order": int(order),
            "slit": setup.centres,
            "profile": profile * setup.flatProfile,
            "profileErr": profileErr * setup.flatProfile,
            "nLines": counts,
        }
    )
    isValid = np.isfinite(piece["profile"]) & np.isfinite(piece["profileErr"]) & (piece["profileErr"] > 0)
    if isKept is not None:
        isValid &= isKept
    return piece.loc[isValid][PROFILE_COLUMNS]


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
        return None, f"{NO_SKYLINES_NOTE_PREFIX} {firstWavelength:.0f} and {lastWavelength:.0f} nm"
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
        setup, note = _prepare_order(data, errs, isUsable, flatData, slit, isInOrder & (orderMap == order), binArcsec)
        if setup is None:
            notes[int(order)] = note
            continue
        isPixel = setup.isPixel
        pixels = _order_pixels(wavelength[isPixel], setup.binIndex, data[isPixel], errs[isPixel] ** 2, step[isPixel])
        combined, note = _order_skyline_profile(pixels, len(setup.centres), lineWavelengths, setup.isInterior)
        if combined is None:
            notes[int(order)] = note
            continue
        piece = _profile_piece(order, setup, *combined)
        if piece.empty:
            notes[int(order)] = f"no slit bin holds at least {MIN_LINES_PER_BIN} usable skylines"
            continue
        pieces.append(piece)
    if not pieces:
        return pd.DataFrame(columns=PROFILE_COLUMNS), notes
    return pd.concat(pieces, ignore_index=True)[PROFILE_COLUMNS], notes


def _chunk_bin_medians(binIndex, values, errs2, chunks, nBins):
    """*per chunk and slit bin: the median of ``values`` and the squared error of that median*

    A bin of a chunk that holds fewer than ``MIN_CHUNK_PIXELS_PER_BIN`` pixels is NaN. The error of the median is
    ``MEDIAN_ERROR_FACTOR`` times the standard error of the mean.
    """
    medians = np.full((len(chunks), nBins), np.nan)
    errors2 = np.full((len(chunks), nBins), np.nan)
    for index, members in enumerate(chunks):
        median, meanErr2, count = _slice_bin_stats(binIndex[members], values[members], errs2[members], nBins)
        isDense = count >= MIN_CHUNK_PIXELS_PER_BIN
        medians[index] = np.where(isDense, median, np.nan)
        errors2[index] = np.where(isDense, MEDIAN_ERROR_FACTOR**2 * meanErr2 / np.maximum(count, 1), np.nan)
    return medians, errors2


def _combine_chunks(medians, errors2, isNormBin, isSelected, minChunksPerBin):
    """*normalise each selected chunk to its median over ``isNormBin`` and combine the chunks per slit bin*

    **Return:**

    - ``profile`` -- the median normalised chunk profile, NaN in bins with fewer than ``minChunksPerBin`` chunks
    - ``profileErr`` -- the larger of the chunk-to-chunk scatter and the propagated ERRS error of the median
    - ``nChunks`` -- the number of chunks in each bin
    """
    selected, selectedErr2 = medians[isSelected], errors2[isSelected]
    level = _nan_median(selected[:, isNormBin], axis=1)
    with np.errstate(invalid="ignore"):
        isGood = np.isfinite(level) & (level > 0)
    normalised = selected[isGood] / level[isGood, None]
    normalisedErr2 = selectedErr2[isGood] / level[isGood, None] ** 2
    nChunks = np.isfinite(normalised).sum(axis=0)
    profile = _nan_median(normalised, axis=0)
    deviation = _nan_median(np.abs(normalised - profile), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        scatterErr = MEDIAN_ERROR_FACTOR * MAD_TO_SIGMA * deviation / np.sqrt(nChunks)
        errsErr = MEDIAN_ERROR_FACTOR / np.sqrt(np.nansum(1.0 / normalisedErr2, axis=0))
    profile = np.where(nChunks >= minChunksPerBin, profile, np.nan)
    return profile, np.fmax(scatterErr, errsErr), nChunks


def _find_object(centres, profile, lower, upper):
    """*the slit position of the object peak in a flat-fielded profile, or None if there is no significant peak*

    The search stays ``OBJECT_SEARCH_MARGIN_ARCSEC`` inside the flat's lit edges. The peak is picked from the profile
    smoothed over ``OBJECT_SMOOTHING_BINS`` bins, and it counts only above ``OBJECT_SIGNIFICANCE`` times the
    smoothed robust scatter of the profile.
    """
    isInside = (centres > lower + OBJECT_SEARCH_MARGIN_ARCSEC) & (centres < upper - OBJECT_SEARCH_MARGIN_ARCSEC)
    isInside &= np.isfinite(profile)
    if isInside.sum() < MIN_OBJECT_SEARCH_BINS:
        return None
    inside = profile[isInside]
    median = np.median(inside)
    scatter = MAD_TO_SIGMA * np.median(np.abs(inside - median))
    smooth = np.convolve(inside - median, np.ones(OBJECT_SMOOTHING_BINS) / OBJECT_SMOOTHING_BINS, mode="same")
    peak = int(np.argmax(smooth))
    significance = smooth[peak] / (scatter / np.sqrt(OBJECT_SMOOTHING_BINS)) if scatter > 0 else np.inf
    return float(centres[isInside][peak]) if significance >= OBJECT_SIGNIFICANCE else None


class continuum_chunking(NamedTuple):
    """How a continuum profile is built: the chunk count, the fewest chunks per bin and the object mask half-width."""

    nChunks: int
    minChunksPerBin: int
    objectMaskArcsec: float


def _order_continuum_profiles(pixels, setup, order, flatEdge, chunking):
    """*the full, blue-half and red-half continuum profiles of one order, or None with the reason there are none*

    ``pixels`` is the order's ``_order_pixels``, which is sorted by wavelength, so equal-count chunks of it are
    wavelength chunks.

    **Return:**

    - ``profiles`` -- ``(full, blue, red)`` dataframes of ``PROFILE_COLUMNS``, or *None*
    - ``note`` -- why there is no profile, or an empty string
    """
    nChunks, minChunksPerBin, objectMaskArcsec = chunking
    nPixels = len(pixels.values)
    if nPixels < nChunks * MIN_CHUNK_PIXELS_PER_BIN:
        return None, f"too few usable pixels ({nPixels}) for {nChunks} wavelength chunks"
    nBins = len(setup.centres)
    chunks = np.array_split(np.arange(nPixels), nChunks)
    medians, errors2 = _chunk_bin_medians(pixels.binIndex, pixels.values, pixels.errs2, chunks, nBins)
    allChunks = np.ones(nChunks, dtype=bool)
    # PASS 1 NORMALISES ON THE WHOLE INTERIOR TO FIND THE OBJECT. PASS 2 NORMALISES ON THE SKY INTERIOR ONLY
    firstPass = _combine_chunks(medians, errors2, setup.isInterior, allChunks, minChunksPerBin)
    objectCentre = _find_object(setup.centres, firstPass[0], flatEdge[0], flatEdge[1])
    isSky = np.ones(nBins, dtype=bool)
    if objectCentre is not None:
        isSky = np.abs(setup.centres - objectCentre) > objectMaskArcsec
    isSkyInterior = setup.isInterior & isSky
    firstHalf = np.arange(nChunks) < nChunks // 2
    profiles = tuple(
        _profile_piece(
            order,
            setup,
            *_combine_chunks(medians, errors2, isSkyInterior, isSelected, minChunksPerBin),
            isKept=isSky,
        )
        for isSelected in (allChunks, firstHalf, ~firstHalf)
    )
    return profiles, ""


def continuum_slit_profiles(
    data,
    errs,
    isUsable,
    flatData,
    wavelength,
    slit,
    orderMap,
    flatEdges,
    orders,
    binArcsec=PROFILE_BIN_ARCSEC,
    nChunks=CONTINUUM_CHUNKS,
    minChunksPerBin=CONTINUUM_MIN_CHUNKS_PER_BIN,
    objectMaskArcsec=CONTINUUM_OBJECT_MASK_ARCSEC,
):
    """*the raw-like sky-continuum slit profile of the given orders, for orders that have no skylines to use*

    Each order's usable pixels are split into ``nChunks`` equal-count wavelength chunks. Every chunk is median-binned
    along the slit and normalised to its sky interior, and the chunks are combined by median. The chunk-to-chunk
    scatter is the error, floored by the propagated ERRS error. The result is multiplied by the order's normalised
    flat profile, as the skyline profile is. The object is found as the strongest significant peak of the combined
    profile, and the bins within ``objectMaskArcsec`` of it are left out of the returned rows and out of the
    normalisation. The blue half and the red half of the chunks are combined the same way, so that the stability of
    the offset can be checked.

    **Key Arguments:**

    - ``data`` -- the flat-fielded, detrended frame
    - ``errs`` -- the frame's 1-sigma uncertainty
    - ``isUsable`` -- boolean image, True where ``data`` and ``errs`` can be used
    - ``flatData`` -- the master flat's flux image
    - ``wavelength`` -- the 2D map's WAVELENGTH image in nm
    - ``slit`` -- the 2D map's SLIT image in arcsec
    - ``orderMap`` -- the 2D map's ORDER image
    - ``flatEdges`` -- dataframe from ``flat_slit_edges``, whose ``lower`` and ``upper`` bound the object search
    - ``orders`` -- the order numbers to build profiles for
    - ``binArcsec`` -- the slit bin width in arcsec. Default *PROFILE_BIN_ARCSEC*
    - ``nChunks`` -- the number of wavelength chunks per order. Default *CONTINUUM_CHUNKS*
    - ``minChunksPerBin`` -- the fewest chunks a slit bin needs. Default *CONTINUUM_MIN_CHUNKS_PER_BIN*
    - ``objectMaskArcsec`` -- the half-width of the object mask in arcsec. Default *CONTINUUM_OBJECT_MASK_ARCSEC*

    **Return:**

    - ``profiles`` -- dataframe with the columns of ``PROFILE_COLUMNS``, where ``nLines`` holds the number of chunks in
      the bin (for a skyline profile it holds the number of skylines). Empty if no order has a profile
    - ``halfProfiles`` -- the ``(blue, red)`` dataframes of the same layout, built from the first and the second half
      of the chunks in wavelength
    - ``notes`` -- ``{order: reason}`` for each order that has no profile rows

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination_profiles import continuum_slit_profiles
    profiles, (blue, red), notes = continuum_slit_profiles(
        frame.data, frame.uncertainty.array, isUsable, flat.data, wavelength, slit, orderMap, flatEdges, [2, 3]
    )
    ```
    """
    edgeByOrder = {int(row["order"]): (row["lower"], row["upper"]) for _, row in flatEdges.iterrows()}
    isInOrder = in_order_pixels(orderMap, slit, wavelength)
    step = _wavelength_step(wavelength)
    chunking = continuum_chunking(nChunks, minChunksPerBin, objectMaskArcsec)
    pieces, notes = ([], [], []), {}
    for order in orders:
        setup, note = _prepare_order(data, errs, isUsable, flatData, slit, isInOrder & (orderMap == order), binArcsec)
        edge = edgeByOrder.get(int(order), (np.nan, np.nan))
        if setup is not None and not np.isfinite(edge).all():
            setup, note = None, "the order has no finite flat edges"
        built = None
        if setup is not None:
            isPixel = setup.isPixel
            pixels = _order_pixels(
                wavelength[isPixel], setup.binIndex, data[isPixel], errs[isPixel] ** 2, step[isPixel]
            )
            built, note = _order_continuum_profiles(pixels, setup, order, edge, chunking)
        if built is None or built[0].empty:
            notes[int(order)] = note or f"no slit bin holds at least {minChunksPerBin} usable wavelength chunks"
            continue
        for kind, piece in zip(pieces, built, strict=True):
            kind.append(piece)
    empty = pd.DataFrame(columns=PROFILE_COLUMNS)
    full, blue, red = (pd.concat(kind, ignore_index=True) if kind else empty for kind in pieces)
    return full, (blue, red), notes
