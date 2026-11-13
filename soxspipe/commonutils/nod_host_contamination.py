"""Host contamination of a nod pair: find extended emission that the A-B subtraction has removed from the object.

The A-B difference removes the sky, and with it any emission that is longer along the slit than the nod throw. A host
galaxy line that covers more of the slit than the throw sits at both nod positions, so part of it is subtracted from
the object, which can leave false dips below the continuum.

For a point source the negative trace of the A-B difference is the object at nod B, so it is the positive trace
scaled by a constant. A narrow feature in the sum of the negative trace and the scaled positive trace is signal that
the object does not explain. This module counts those features in each order:

1. The positive and the negative trace are the slit positions with the most extreme median flux.
2. The median of each detector line outside both trace apertures is subtracted. It is the sky residual left by a sky
   that changed between the nods, which is flat along the slit.
3. Each aperture gives one flux per detector line: the median of its pixels times their number, so an unmasked cosmic
   ray one or two pixels long hardly moves it.
4. A feature is a run of lines where both the residual and the negative trace itself rise above the limit, in units
   of the larger of the propagated and the local noise. Lines near a catalogue skyline and at the order ends are not
   searched.

It only flags the features; it changes no data. A host that is symmetric about the object leaves no residual at the
negative trace, so this check cannot see it.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from soxspipe.commonutils.nod_shift_and_add import slit_pixel_scale

# THE SCALE FROM A MEDIAN ABSOLUTE DEVIATION TO A STANDARD DEVIATION
MAD_TO_SIGMA = 1.4826
# THE VARIANCE OF THE MEDIAN OF GAUSSIAN SAMPLES IS ABOUT THIS FACTOR TIMES THE VARIANCE OF THEIR MEAN
MEDIAN_VARIANCE_FACTOR = np.pi / 2
# THE DEFAULT DETECTION LIMIT, IN UNITS OF THE LOCAL NOISE OF THE APERTURE RESIDUAL
DEFAULT_SIGMA_LIMIT = 5.0
# THE WINDOW (IN DETECTOR LINES) OF THE RUNNING MEDIAN THAT REMOVES THE SMOOTH PART OF THE RESIDUAL, AND OF THE
# RUNNING MAD THAT MEASURES ITS LOCAL NOISE. A NARROW LINE IS A FEW DETECTOR LINES WIDE
TREND_WINDOW_LINES = 51
# A FEATURE MUST SPAN AT LEAST THIS MANY DETECTOR LINES, SO A SINGLE HOT PIXEL OR COSMIC RAY DOES NOT COUNT
MIN_FEATURE_LINES = 2
# FLAGGED LINES THIS FAR APART OR CLOSER BELONG TO ONE FEATURE
MAX_FEATURE_GAP_LINES = 1
# DETECTOR LINES WITHIN THIS MANY LINES OF A CATALOGUE SKYLINE ARE NOT SEARCHED, AS SKY RESIDUALS ARE LARGE THERE
SKY_MASK_HALF_WIDTH_LINES = 8
# THE POSITIVE TRACE SCALED BY THIS RATIO MODELS THE NEGATIVE TRACE. A RATIO OUTSIDE THESE LIMITS MEANS THE TWO
# TRACES ARE NOT THE SAME OBJECT, SO THE ORDER IS SKIPPED
MIN_TRACE_RATIO = 0.25
MAX_TRACE_RATIO = 4.0
# FEATURES WHOSE WAVELENGTH RANGES ARE THIS CLOSE (NM) ARE ONE FEATURE WHEN FEATURES ARE MERGED ACROSS ORDERS
MERGE_GAP_NM = 0.2


@dataclass(frozen=True)
class host_feature:
    """A narrow feature at the negative trace that the object does not explain.

    ``order`` is the order it was found in (``None`` after a merge across orders). ``significance`` is the largest
    absolute residual in the feature, in units of the local noise.
    """

    order: int | None
    wavelengthMinNm: float
    wavelengthMaxNm: float
    significance: float


def find_host_contamination(
    difference, mapImages, aperturePixels, sigmaLimit=DEFAULT_SIGMA_LIMIT, skylineWavelengthsNm=()
):
    """*find the narrow features at the negative trace of an A-B difference that the object does not explain*

    **Key Arguments:**

    - ``difference`` -- the A-B difference, a CCDData frame with a standard-deviation uncertainty
    - ``mapImages`` -- the ``two_d_map_images`` of the frame's binning
    - ``aperturePixels`` -- the full aperture length in binned pixels along the slit
    - ``sigmaLimit`` -- a feature must exceed this many times the local noise. Default *5*
    - ``skylineWavelengthsNm`` -- wavelengths (nm) of bright skylines whose neighbourhood is not searched

    **Return:**

    - ``features`` -- a tuple of ``host_feature``, ordered by order and wavelength. Empty when none is found
    """
    data = np.asarray(difference.data, dtype=float)
    error = np.asarray(difference.uncertainty.array, dtype=float)
    isUsable = np.isfinite(data) & np.isfinite(error)
    if difference.mask is not None:
        isUsable &= ~np.asarray(difference.mask, dtype=bool)
    skylines = np.asarray(skylineWavelengthsNm, dtype=float)

    features = []
    for order in mapImages.orders:
        isOrder = (mapImages.order == order) & isUsable & np.isfinite(mapImages.wavelength)
        if not isOrder.any():
            continue
        rows, columns = np.nonzero(isOrder)
        stepRow = float(np.median(np.abs(mapImages.wavelengthRow[isOrder])))
        stepColumn = float(np.median(np.abs(mapImages.wavelengthCol[isOrder])))
        pixels = pd.DataFrame(
            {
                # A DETECTOR LINE ALONG THE SLIT: A ROW WHEN THE WAVELENGTH CHANGES FASTEST DOWN THE ROWS
                "line": rows if stepRow >= stepColumn else columns,
                "slit": mapImages.slit[isOrder],
                "wavelength": mapImages.wavelength[isOrder],
                "flux": data[isOrder],
                "variance": error[isOrder] ** 2,
            }
        )
        scales = (aperturePixels, slit_pixel_scale(mapImages, order), max(stepRow, stepColumn))
        features.extend(_order_features(order, pixels, scales, sigmaLimit, skylines))
    return tuple(features)


def _order_features(order, pixels, scales, sigmaLimit, skylines):
    """*the host features of one order, from its usable pixels and its (aperture pixels, arcsec/pixel, nm/line)*"""
    aperturePixels, slitScale, stepNm = scales
    apertureArcsec = aperturePixels * slitScale
    traces = _trace_positions(pixels, apertureArcsec, slitScale)
    if traces is None:
        return []
    pixels = _without_slit_background(pixels, traces, apertureArcsec)
    spectra = _paired_spectra(*(_aperture_spectrum(pixels, position, apertureArcsec) for position in traces))
    ratio = _trace_ratio(spectra) if len(spectra) else None
    if ratio is None:
        return []
    residual = spectra["flux_neg"] + ratio * spectra["flux_pos"]
    propagated = np.sqrt(spectra["variance_neg"] + ratio**2 * spectra["variance_pos"])
    significance = _significance(residual, propagated)
    # THE FEATURE MUST BE AT THE OPPOSITE NOD POSITION: A SPIKE IN THE POSITIVE TRACE ALONE ALSO LEAVES A RESIDUAL
    isAtNegativeTrace = _significance(spectra["flux_neg"], np.sqrt(spectra["variance_neg"])) > sigmaLimit
    significance[~isAtNegativeTrace] = 0.0
    wavelengths = spectra["wavelength"].to_numpy()
    significance[_near_skylines(wavelengths, skylines, stepNm)] = 0.0
    return [
        host_feature(
            order=int(order),
            wavelengthMinNm=float(wavelengths[first]),
            wavelengthMaxNm=float(wavelengths[last]),
            significance=float(np.max(significance[first : last + 1])),
        )
        for first, last in _flagged_runs(significance > sigmaLimit)
    ]


def _trace_positions(pixels, apertureArcsec, slitScale):
    """*the slit positions of the positive and the negative trace, or None when they are less than an aperture apart*

    The slit is cut into bins one pixel long, and each trace is the bin with the most extreme median flux.
    """
    slitBin = np.round(pixels["slit"] / slitScale)
    profile = pixels.groupby(slitBin)["flux"].median()
    if profile.empty or profile.max() <= 0 or profile.min() >= 0:
        return None
    positive = float(profile.idxmax() * slitScale)
    negative = float(profile.idxmin() * slitScale)
    if abs(positive - negative) < apertureArcsec:
        return None
    return positive, negative


def _without_slit_background(pixels, traces, apertureArcsec):
    """*the pixels less the median of each detector line away from both traces*

    The sky changes between the nod exposures, so the A-B difference keeps a sky residual that is flat along the slit,
    in the continuum and in the skylines. The median of the pixels outside both trace apertures measures
    it. A line with no such pixel keeps its flux.
    """
    isBackground = np.ones(len(pixels), dtype=bool)
    for position in traces:
        isBackground &= np.abs(pixels["slit"].to_numpy() - position) > apertureArcsec / 2
    background = pixels.loc[isBackground].groupby("line")["flux"].median()
    level = pixels["line"].map(background).fillna(0.0).to_numpy()
    return pixels.assign(flux=pixels["flux"].to_numpy() - level)


def _aperture_spectrum(pixels, position, apertureArcsec):
    """*the flux and its variance, and the mean wavelength, of each detector line of an aperture along the slit*

    The flux is the median of the aperture pixels times their number. Emission that fills the aperture moves the
    median as it moves the sum, but an unmasked cosmic ray one or two pixels long along the slit hardly moves it.
    A line with fewer aperture pixels than the median line, for example at the end of an order, is dropped.
    """
    inAperture = pixels.loc[np.abs(pixels["slit"] - position) <= apertureArcsec / 2]
    grouped = inAperture.groupby("line")
    count = grouped["flux"].count()
    spectrum = pd.DataFrame(
        {
            "wavelength": grouped["wavelength"].mean(),
            "flux": grouped["flux"].median() * count,
            "variance": MEDIAN_VARIANCE_FACTOR * grouped["variance"].sum(),
            "count": count,
        }
    )
    return spectrum.loc[spectrum["count"] >= spectrum["count"].median()].sort_values("wavelength")


def _paired_spectra(positive, negative):
    """*the positive trace spectrum with the negative trace spectrum interpolated onto its wavelengths*

    The two traces of a line sit at slightly different wavelengths when the lines of constant wavelength tilt across
    the slit, so the negative trace is resampled before the two are compared. Wavelengths outside the negative trace
    are dropped.
    """
    if len(positive) < 2 or len(negative) < 2:
        return pd.DataFrame(columns=["wavelength", "flux_pos", "variance_pos", "flux_neg", "variance_neg"])
    wavelengths = positive["wavelength"].to_numpy()
    isCovered = (wavelengths >= negative["wavelength"].iloc[0]) & (wavelengths <= negative["wavelength"].iloc[-1])
    wavelengths = wavelengths[isCovered]
    return pd.DataFrame(
        {
            "wavelength": wavelengths,
            "flux_pos": positive["flux"].to_numpy()[isCovered],
            "variance_pos": positive["variance"].to_numpy()[isCovered],
            "flux_neg": np.interp(wavelengths, negative["wavelength"], negative["flux"]),
            "variance_neg": np.interp(wavelengths, negative["wavelength"], negative["variance"]),
        }
    )


def _trace_ratio(spectra):
    """*the ratio of the negative trace to the positive trace, or None when it is not a usable constant*"""
    positiveLevel = float(np.median(spectra["flux_pos"]))
    negativeLevel = float(np.median(spectra["flux_neg"]))
    if positiveLevel <= 0 or negativeLevel >= 0:
        return None
    ratio = -negativeLevel / positiveLevel
    if not MIN_TRACE_RATIO <= ratio <= MAX_TRACE_RATIO:
        return None
    return ratio


def _significance(residual, propagated):
    """*the absolute residual above its running median, in units of the larger of the propagated and local noise*"""
    window = {"window": TREND_WINDOW_LINES, "center": True, "min_periods": 1}
    detrended = residual - residual.rolling(**window).median()
    localNoise = MAD_TO_SIGMA * detrended.abs().rolling(**window).median()
    noise = np.maximum(propagated.to_numpy(), localNoise.to_numpy())
    with np.errstate(divide="ignore", invalid="ignore"):
        significance = np.nan_to_num(np.abs(detrended.to_numpy()) / noise, nan=0.0, posinf=0.0)
    # THE RUNNING MEDIAN IS ONE-SIDED AT THE ENDS OF AN ORDER, WHERE THE FLUX ALSO FALLS OFF, SO THEY ARE NOT SEARCHED
    edge = TREND_WINDOW_LINES // 2
    significance[:edge] = 0.0
    significance[-edge:] = 0.0
    return significance


def _near_skylines(wavelengths, skylines, stepNm):
    """*a mask of the wavelengths that lie within the sky mask of a catalogue skyline*"""
    if not len(skylines):
        return np.zeros(len(wavelengths), dtype=bool)
    distance = np.min(np.abs(wavelengths[:, None] - skylines[None, :]), axis=1)
    return distance <= SKY_MASK_HALF_WIDTH_LINES * stepNm


def _flagged_runs(isFlagged):
    """*the (first, last) indices of each run of flagged detector lines long enough to be a feature*"""
    indices = np.flatnonzero(isFlagged)
    if not len(indices):
        return []
    breaks = np.flatnonzero(np.diff(indices) > MAX_FEATURE_GAP_LINES + 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [len(indices) - 1]))
    return [
        (int(indices[start]), int(indices[end]))
        for start, end in zip(starts, ends, strict=True)
        if indices[end] - indices[start] + 1 >= MIN_FEATURE_LINES
    ]


def merge_host_features(features):
    """*merge features whose wavelength ranges overlap, or nearly do, across orders and nodding cycles*

    **Key Arguments:**

    - ``features`` -- the ``host_feature`` objects to merge

    **Return:**

    - ``merged`` -- a tuple of ``host_feature`` with no order, sorted by wavelength. One per distinct feature
    """
    merged = []
    for feature in sorted(features, key=lambda item: item.wavelengthMinNm):
        if merged and feature.wavelengthMinNm <= merged[-1].wavelengthMaxNm + MERGE_GAP_NM:
            last = merged[-1]
            merged[-1] = host_feature(
                order=None,
                wavelengthMinNm=last.wavelengthMinNm,
                wavelengthMaxNm=max(last.wavelengthMaxNm, feature.wavelengthMaxNm),
                significance=max(last.significance, feature.significance),
            )
            continue
        merged.append(host_feature(None, feature.wavelengthMinNm, feature.wavelengthMaxNm, feature.significance))
    return tuple(merged)


def host_contamination_message(features):
    """*the warning that names the wavelength ranges of merged host features*

    **Key Arguments:**

    - ``features`` -- the merged ``host_feature`` objects, from ``merge_host_features``

    **Return:**

    - ``message`` -- the warning text, or *None* when there are no features
    """
    if not features:
        return None
    ranges = ", ".join(f"{feature.wavelengthMinNm:.2f}-{feature.wavelengthMaxNm:.2f} nm" for feature in features)
    return (
        f"Extended emission at the opposite nod position may have been subtracted from the object spectrum at "
        f"{ranges}. A host line longer along the slit than the nod throw can leave false dips there."
    )
