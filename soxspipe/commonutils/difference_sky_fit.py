"""Helpers of the sky fit on a nod A-B difference frame (DY-2057 prototype).

In an A-B difference most of the sky has cancelled. What is left is a small residual that is the same at every slit
position, plus a positive copy of the object at the A position and a negative copy at the B position. The helpers
here supply what the spline sky fit needs to work on such a frame: a two-sided object mask, the noise model, and the
rule that stops the fit from subtracting noise.
"""

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import median_filter

# THE SKY AND READ-NOISE VARIANCE OF AN A-B DIFFERENCE ARE EACH TWICE THAT OF ONE FRAME (KELSON, longslit_reduce.pro
# LINES 3224-3231). THE BIAS TERM IS DROPPED: BOTH FRAMES HAD THE SAME BIAS REMOVED
DIFFERENCE_VARIANCE_FACTOR = 2.0
# THE SHIFT-AND-ADD FRAME OF DY-1674 IS THE MEAN OF TWO INDEPENDENT A-B DIFFERENCES, (2 * VAR + 2 * VAR) / 4, SO ITS
# VARIANCE IS THAT OF ONE FRAME
COMBINED_VARIANCE_FACTOR = 1.0
# THE MODEL COUNTS AS SKY ONLY WHERE IT IS AT LEAST THIS MANY STANDARD ERRORS OF THE LOCAL WEIGHTED MEAN
MIN_MODEL_SIGNIFICANCE = 3.0
# THE NUMBER OF WAVELENGTH-SORTED SAMPLES OVER WHICH THE STANDARD ERROR OF THE LOCAL WEIGHTED MEAN IS TAKEN. IT IS THE
# SPACING OF THE STARTER KNOTS OF A NIR ORDER (300), MADE ODD SO THE WINDOW IS CENTRED, THE SMOOTHEST THE SPLINE CAN BE
GATE_WINDOW_SAMPLES = 301
# THE NUMBER OF WAVELENGTH-SORTED SAMPLES OF THE ROLLING MEDIAN THAT SMOOTHS THE SKY LEVEL BEHIND THE NOISE MODEL. THE
# MEDIAN IGNORES THE TRACES, WHICH FILL ONLY A FEW SLIT PIXELS OF EVERY COLUMN
NOISE_SKY_WINDOW = 501

# THE OBJECT MASK READS THE SIGNED SLIT PROFILE OF AN ORDER IN THIS MANY BINS
OBJECT_SLIT_BINS = 100
# THE BINS AT EACH END OF THE SLIT THAT ARE NOT EXAMINED, AS THE SLIT-EDGE CLIP HAS TAKEN THEM
OBJECT_EDGE_BINS = 5
# A BIN IS AN OBJECT SEED WHEN ITS MEAN SIGNAL-TO-NOISE DIFFERS FROM THE ORDER'S MEDIAN BY THIS MANY STANDARD ERRORS,
# AND THE OBJECT GROWS OVER THE NEIGHBOURING BINS THAT DIFFER BY MORE THAN OBJECT_EXTEND_SIGMA
OBJECT_SEED_SIGMA = 5.0
OBJECT_EXTEND_SIGMA = 2.0
# THE RANGE IS WIDENED BY THIS MANY BINS ON EACH SIDE
OBJECT_DILATE_BINS = 2
# A BIN NEEDS AT LEAST THIS MANY PIXELS TO BE JUDGED
OBJECT_MIN_BIN_SAMPLES = 20
# THE SCATTER OF THE SLIT PROFILE COMES FROM THE QUIETEST BINS: THE LOWEST DECILE OF THE ABSOLUTE DEVIATIONS, WHICH IS
# 0.1257 STANDARD DEVIATIONS OF A GAUSSIAN
OBJECT_QUIET_PERCENTILE = 10
OBJECT_QUIET_DECILE_TO_SIGMA = 0.1257
# THE RUNNING MEDIAN OF THE SLIT PROFILE THAT IS TAKEN AS ITS BASELINE IS THIS MANY BINS WIDE, ABOUT 60% OF THE SLIT
OBJECT_BASELINE_BINS = 63


@dataclass(frozen=True)
class difference_fit:
    """How the sky fit treats a difference frame.

    ``noiseSkyFrame`` is the frame whose flux is the sky level that sets the photon noise, as a CCDData object with the
    shape of the difference frame. ``subtract_sky`` refuses a ``difference_fit`` without one. A ``subtract_sky`` that
    was not built through ``__init__`` can instead carry its own ``noise_sky`` column on the pixels.
    ``varianceFactor`` multiplies the sky and read-noise variance of the noise model, 2 for one A-B difference.
    ``minSignificance`` is the significance below which the fitted model is dropped.
    """

    noiseSkyFrame: object = None
    varianceFactor: float = DIFFERENCE_VARIANCE_FACTOR
    minSignificance: float = MIN_MODEL_SIGNIFICANCE


def _binned_signal_to_noise(slitPosition, flux, error):
    """*the mean signal-to-noise of the usable pixels in each slit bin, and the bin edges and pixel counts*"""
    usable = np.isfinite(slitPosition) & np.isfinite(flux) & np.isfinite(error) & (error > 0)
    if not usable.any():
        return None
    edges = np.linspace(slitPosition[usable].min(), slitPosition[usable].max(), OBJECT_SLIT_BINS + 1)
    binIndex = np.clip(np.digitize(slitPosition[usable], edges) - 1, 0, OBJECT_SLIT_BINS - 1)
    counts = np.bincount(binIndex, minlength=OBJECT_SLIT_BINS)
    sums = np.bincount(binIndex, weights=(flux / error)[usable], minlength=OBJECT_SLIT_BINS)
    profile = np.full(OBJECT_SLIT_BINS, np.nan)
    hasSamples = counts >= OBJECT_MIN_BIN_SAMPLES
    profile[hasSamples] = sums[hasSamples] / counts[hasSamples]
    return profile, counts, edges


def _quiet_bin_level(profile, counts, isJudged):
    """*the profile minus its slow baseline, and the standard error of a bin, both read from the quietest bins*

    The two nod positions differ in flat-field and illumination, so the slit profile of a difference has a slope that
    is not an object. A bright object and its wings, in turn, can fill most of the slit. The baseline is therefore the
    running median of the bins over ``OBJECT_BASELINE_BINS``, wide enough to pass over the traces. What is left is
    the objects and the noise. The lowest decile of its absolute values sits in the bins with no object, and gives the
    scatter of a bin. The pure-noise standard error of a bin is ``1 / sqrt(count)``, which the measured scatter may
    only raise.
    """
    filled = np.where(isJudged, profile, np.median(profile[isJudged]))
    baseline = median_filter(filled, size=OBJECT_BASELINE_BINS, mode="nearest")
    deviation = np.where(isJudged, filled - baseline, 0.0)
    decile = np.percentile(np.abs(deviation[isJudged]), OBJECT_QUIET_PERCENTILE)
    standardError = np.fmax(1.0 / np.sqrt(np.fmax(counts, 1)), decile / OBJECT_QUIET_DECILE_TO_SIGMA)
    return deviation, standardError


def signed_object_slit_ranges(slitPosition, flux, error):
    """*the slit-position ranges of the positive and negative objects of a difference order*

    The signed slit profile is the mean signal-to-noise of the pixels in each slit bin. A uniform sky residual lifts
    every bin alike, so it is removed with the median of the bins. A bin whose profile departs from that median by
    ``OBJECT_SEED_SIGMA`` standard errors seeds an object of either sign, which grows over the neighbouring bins that
    depart by more than ``OBJECT_EXTEND_SIGMA`` and is widened by ``OBJECT_DILATE_BINS``.

    **Key Arguments:**

    - ``slitPosition`` -- slit position of each pixel (arcsec)
    - ``flux`` -- flux of each pixel
    - ``error`` -- flux error of each pixel

    **Return:**

    - ``objectRanges`` -- a ``[lower, upper]`` slit-position pair for each object, of either sign
    """
    binned = _binned_signal_to_noise(
        np.asarray(slitPosition, dtype=float), np.asarray(flux, dtype=float), np.asarray(error, dtype=float)
    )
    if binned is None:
        return []
    profile, counts, edges = binned
    isJudged = np.isfinite(profile)
    if isJudged.sum() <= 2 * OBJECT_EDGE_BINS:
        return []
    deviation, standardError = _quiet_bin_level(profile, counts, isJudged)
    significance = np.where(isJudged, np.abs(deviation) / standardError, 0.0)
    examined = np.zeros(OBJECT_SLIT_BINS, dtype=bool)
    examined[OBJECT_EDGE_BINS : OBJECT_SLIT_BINS - OBJECT_EDGE_BINS] = True
    isSeed = (significance > OBJECT_SEED_SIGMA) & examined
    if not isSeed.any():
        return []
    isStanding = significance > OBJECT_EXTEND_SIGMA
    isObject = np.zeros(OBJECT_SLIT_BINS, dtype=bool)
    for seed in np.flatnonzero(isSeed):
        for step in (-1, 1):
            position = seed
            while 0 <= position < OBJECT_SLIT_BINS and (position == seed or isStanding[position]):
                isObject[position] = True
                position += step
    isObject = np.convolve(isObject, np.ones(2 * OBJECT_DILATE_BINS + 1), mode="same") > 0
    return _runs_to_ranges(isObject, edges)


def _runs_to_ranges(isObject, edges):
    """*the ``[lower, upper]`` edges of every run of True bins*"""
    padded = np.concatenate(([False], isObject, [False]))
    changes = np.flatnonzero(padded[1:] != padded[:-1])
    return [[float(edges[start]), float(edges[stop])] for start, stop in zip(changes[::2], changes[1::2], strict=True)]


def gate_insignificant_model(model, weights, isUsed, minSignificance, window=GATE_WINDOW_SAMPLES):
    """*zero the sky model wherever it is within the noise of its own local mean*

    A spline fitted to a frame with almost no sky still follows the noise a little. At each sample the standard error
    of the local weighted mean is ``1 / sqrt(sum of the weights squared)`` over ``window`` neighbouring samples. The
    model is kept only where its size is at least ``minSignificance`` of those errors, so a skyline residual survives
    and a noise ripple is dropped. What is dropped is below the noise of one pixel, by construction.

    **Key Arguments:**

    - ``model`` -- the fitted sky model of each pixel, sorted by wavelength
    - ``weights`` -- the fit weights (1 / sigma) of each pixel
    - ``isUsed`` -- True for the pixels the fit used
    - ``minSignificance`` -- the model survives where it is at least this many standard errors from zero
    - ``window`` -- the number of neighbouring samples of the local mean, cut to the length of the order when the order
      is shorter. Default ``GATE_WINDOW_SAMPLES``

    **Return:**

    - ``gatedModel`` -- a new array. A model value that is NaN comes out as 0
    """
    model = np.asarray(model, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if model.size == 0:
        return model.copy()
    usableWeights = np.where(np.asarray(isUsed, dtype=bool) & np.isfinite(weights), weights, 0.0)
    # A WINDOW LONGER THAN THE ORDER WOULD MAKE THE "same" CONVOLUTION RETURN THE WINDOW'S LENGTH
    weightSum = np.convolve(usableWeights**2, np.ones(min(window, model.size)), mode="same")
    with np.errstate(divide="ignore"):
        standardError = 1.0 / np.sqrt(weightSum)
    return np.where(np.abs(model) >= minSignificance * standardError, model, 0.0)
