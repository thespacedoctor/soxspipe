#!/usr/bin/env python
"""
*sub-sample centres of sky emission lines, robust to the flat-topped, dipped profiles of wide slits*

Each detected peak is assigned to the line it sits on. Peaks whose line footprints overlap (for example the two
horns of a wide-slit line) are one line, and each line gets one baseline, one half-maximum level and one centre:
the midpoint of its outermost half-maximum crossings.

Author
: David Young

Date Created
: October 15, 2026
"""

import numpy as np

# HALF-WIDTH IN SAMPLES OF THE WINDOW AROUND A LINE THAT SUPPLIES ITS BASELINE
SKYLINE_BASELINE_HALF_WINDOW = 50
# LOW PERCENTILE OF THE WINDOW USED AS A FIRST BASELINE, BEFORE THE LINE ITSELF IS KNOWN
SKYLINE_BASELINE_PERCENTILE = 20
# A LINE'S FOOTPRINT IS THE RUN OF SAMPLES ABOVE THIS FRACTION OF ITS PEAK HEIGHT OVER THE FIRST BASELINE
SKYLINE_FOOTPRINT_FRACTION = 0.1
# FRACTION OF THE LINE HEIGHT ABOVE BASELINE AT WHICH THE CROSSINGS ARE MEASURED
SKYLINE_CENTRE_LEVEL = 0.5


def measure_line_centres(coordinates, skyVals, peaks):
    """*return the sorted half-maximum midpoints of the sky lines holding the detected peaks*

    **Key Arguments:**

    - ``coordinates`` -- ascending wavelength (or pixel) coordinate of each sky sample
    - ``skyVals`` -- the sky spectrum (any constant or slowly varying baseline)
    - ``peaks`` -- indices of the detected sky peaks

    **Return:**

    - ``centres`` -- one centre per line. A line whose half-maximum runs off the end of the order, or which does
      not rise above its baseline, is dropped.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    skyVals = np.asarray(skyVals, dtype=float)
    footprints = [_line_footprint(skyVals, int(peak)) for peak in peaks]
    lines = _merge_overlapping([footprint for footprint in footprints if footprint is not None])
    centres = [_half_maximum_centre(coordinates, skyVals, left, right) for left, right in lines]
    return np.sort([centre for centre in centres if np.isfinite(centre)])


def _window(length, index):
    """*return the slice bounds of the baseline window around ``index``*"""
    return max(0, index - SKYLINE_BASELINE_HALF_WINDOW), min(length, index + SKYLINE_BASELINE_HALF_WINDOW + 1)


def _line_footprint(skyVals, peak):
    """*return the first and last index of the line holding ``peak``, or None when the peak is not a line*"""
    if not np.isfinite(skyVals[peak]):
        return None
    lower, upper = _window(len(skyVals), peak)
    firstBaseline = np.nanpercentile(skyVals[lower:upper], SKYLINE_BASELINE_PERCENTILE)
    if not skyVals[peak] > firstBaseline:
        return None
    return _samples_above(skyVals, peak, firstBaseline + SKYLINE_FOOTPRINT_FRACTION * (skyVals[peak] - firstBaseline))


def _merge_overlapping(footprints):
    """*merge footprints that overlap or touch, so all the peaks of one line become one (left, right) pair*"""
    lines = []
    for left, right in sorted(footprints):
        if lines and left <= lines[-1][1] + 1:
            lines[-1] = (lines[-1][0], max(lines[-1][1], right))
        else:
            lines.append((left, right))
    return lines


def _half_maximum_centre(coordinates, skyVals, left, right):
    """*return the midpoint of the outermost half-maximum crossings of the line from ``left`` to ``right``, or NaN*"""
    top = left + int(np.nanargmax(skyVals[left : right + 1]))

    # THE BASELINE COMES FROM THE WINDOW WITH THE LINE ITSELF LEFT OUT, SO A WIDE LINE CANNOT RAISE IT
    lower, upper = _window(len(skyVals), top)
    outside = np.r_[skyVals[lower:left], skyVals[right + 1 : upper]]
    if not np.isfinite(outside).any():
        return np.nan
    baseline = np.nanmedian(outside)
    if not skyVals[top] > baseline:
        return np.nan
    level = baseline + SKYLINE_CENTRE_LEVEL * (skyVals[top] - baseline)

    # OUTERMOST SAMPLES ABOVE THE LEVEL, SO A DIP BELOW THE LEVEL INSIDE THE LINE DOES NOT SPLIT IT
    aboveLevel = np.flatnonzero(skyVals[left : right + 1] > level)
    first, _ = _samples_above(skyVals, left + aboveLevel[0], level)
    _, last = _samples_above(skyVals, left + aboveLevel[-1], level)
    if first == 0 or last == len(skyVals) - 1:
        return np.nan

    # LINEAR INTERPOLATION OF EACH CROSSING BETWEEN THE LAST SAMPLE ABOVE AND THE FIRST BELOW THE LEVEL
    if not np.isfinite(skyVals[[first - 1, last + 1]]).all():
        return np.nan
    leftCrossing = np.interp(level, skyVals[[first - 1, first]], coordinates[[first - 1, first]])
    rightCrossing = np.interp(level, skyVals[[last + 1, last]], coordinates[[last + 1, last]])
    return 0.5 * (leftCrossing + rightCrossing)


def _samples_above(skyVals, start, level):
    """*return the first and last index of the contiguous run of samples above ``level`` around ``start``*"""
    left = start
    while left > 0 and skyVals[left - 1] > level:
        left -= 1
    right = start
    while right < len(skyVals) - 1 and skyVals[right + 1] > level:
        right += 1
    return left, right
