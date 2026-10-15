#!/usr/bin/env python
"""
*sub-sample centres of sky emission lines, robust to the flat-topped, dipped profiles of wide slits*

Adjacent detected peaks are one line when the valley between them stays at or above ``SKYLINE_VALLEY_FRACTION`` of
the lower peak's height above the baseline (for example the two horns of a wide-slit line). Otherwise they are
separate lines: a valley that falls towards the baseline cannot be told apart from a resolved doublet. Each line
gets one baseline, one half-maximum level and one centre: the midpoint of its outermost half-maximum crossings.

``SKYLINE_VALLEY_FRACTION`` comes from real VIS sky (DY-1599). The dips between the horns of 5-arcsec lines in
orders r and i are 0.74-0.99 of the lower horn; the valleys between adjacent lines on a 1-arcsec frame are
-0.39-0.11. Peaks closer than ``find_peaks``'s ``distance`` never reach this module, so a blend closer than that is
one line whatever its valley.

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
# ADJACENT PEAKS ARE ONE LINE WHEN THE VALLEY BETWEEN THEM KEEPS THIS FRACTION OF THE LOWER PEAK'S HEIGHT
SKYLINE_VALLEY_FRACTION = 0.5
# THE SAMPLES ABOVE THIS FRACTION OF THE LINE HEIGHT ARE LEFT OUT OF THE LINE'S OWN BASELINE
SKYLINE_FOOTPRINT_FRACTION = 0.1
# FRACTION OF THE LINE HEIGHT ABOVE BASELINE AT WHICH THE CROSSINGS ARE MEASURED
SKYLINE_CENTRE_LEVEL = 0.5


def measure_line_centres(coordinates, skyVals, peaks):
    """*return the sorted half-maximum midpoints of the sky lines holding the detected peaks*

    **Key Arguments:**

    - ``coordinates`` -- ascending wavelength (or pixel) coordinate of each sky sample
    - ``skyVals`` -- the sky spectrum (any constant or slowly varying baseline), the same length as ``coordinates``
    - ``peaks`` -- indices of the detected sky peaks; indices outside the spectrum are ignored

    **Return:**

    - ``centres`` -- one centre per line. A line is dropped when its half maximum runs off the end of the order,
      when a sample next to one of its half-maximum crossings is NaN, or when it does not rise above its baseline.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    skyVals = np.asarray(skyVals, dtype=float)
    if len(coordinates) != len(skyVals):
        raise ValueError(f"coordinates ({len(coordinates)}) and skyVals ({len(skyVals)}) must be the same length")

    peaks = np.unique([int(peak) for peak in peaks if 0 <= peak < len(skyVals)])
    peaks = [peak for peak in peaks if np.isfinite(skyVals[peak])]
    centres = [_half_maximum_centre(coordinates, skyVals, line) for line in _group_peaks_into_lines(skyVals, peaks)]
    return np.sort([centre for centre in centres if np.isfinite(centre)])


def _first_baseline(skyVals, first, last):
    """*return a low percentile of the window around samples ``first`` to ``last``*"""
    lower = max(0, first - SKYLINE_BASELINE_HALF_WINDOW)
    upper = min(len(skyVals), last + SKYLINE_BASELINE_HALF_WINDOW + 1)
    return np.nanpercentile(skyVals[lower:upper], SKYLINE_BASELINE_PERCENTILE)


def _group_peaks_into_lines(skyVals, peaks):
    """*return the sorted peaks grouped into lines, joining neighbours whose valley stays high enough*"""
    lines = []
    for peak in peaks:
        if lines and _is_one_line(skyVals, lines[-1][-1], peak):
            lines[-1].append(peak)
        else:
            lines.append([peak])
    return lines


def _is_one_line(skyVals, left, right):
    """*return True when the valley between peaks ``left`` and ``right`` keeps the valley fraction of the lower peak*"""
    baseline = _first_baseline(skyVals, left, right)
    lowerHeight = min(skyVals[left], skyVals[right]) - baseline
    valley = np.min(skyVals[left : right + 1]) - baseline
    return bool(lowerHeight > 0 and valley >= SKYLINE_VALLEY_FRACTION * lowerHeight)


def _half_maximum_centre(coordinates, skyVals, line):
    """*return the midpoint of the outermost half-maximum crossings of the line holding peaks ``line``, or NaN*"""
    first, last = line[0], line[-1]
    top = first + int(np.nanargmax(skyVals[first : last + 1]))
    baseline = _line_baseline(skyVals, first, last, top)
    if not skyVals[top] > baseline:
        return np.nan
    level = baseline + SKYLINE_CENTRE_LEVEL * (skyVals[top] - baseline)

    # WALK OUT FROM THE OUTERMOST SAMPLES ABOVE THE LEVEL, SO A DIP INSIDE THE LINE DOES NOT SPLIT IT
    aboveLevel = np.flatnonzero(skyVals[first : last + 1] > level)
    left, _ = _samples_above(skyVals, first + aboveLevel[0], level)
    _, right = _samples_above(skyVals, first + aboveLevel[-1], level)
    if left == 0 or right == len(skyVals) - 1 or not np.isfinite(skyVals[[left - 1, right + 1]]).all():
        return np.nan

    # LINEAR INTERPOLATION OF EACH CROSSING BETWEEN THE LAST SAMPLE ABOVE AND THE FIRST BELOW THE LEVEL
    leftCrossing = np.interp(level, skyVals[[left - 1, left]], coordinates[[left - 1, left]])
    rightCrossing = np.interp(level, skyVals[[right + 1, right]], coordinates[[right + 1, right]])
    return 0.5 * (leftCrossing + rightCrossing)


def _line_baseline(skyVals, first, last, top):
    """*return the median of the window around ``top`` with the line's own footprint left out, or NaN*"""
    firstBaseline = _first_baseline(skyVals, first, last)
    footprintLevel = firstBaseline + SKYLINE_FOOTPRINT_FRACTION * (skyVals[top] - firstBaseline)
    left, _ = _samples_above(skyVals, first, footprintLevel)
    _, right = _samples_above(skyVals, last, footprintLevel)

    # THE LINE ITSELF IS LEFT OUT, SO A WIDE LINE CANNOT RAISE ITS OWN BASELINE
    lower = max(0, top - SKYLINE_BASELINE_HALF_WINDOW)
    upper = min(len(skyVals), top + SKYLINE_BASELINE_HALF_WINDOW + 1)
    outside = np.r_[skyVals[lower:left], skyVals[right + 1 : upper]]
    if not np.isfinite(outside).any():
        return np.nan
    return np.nanmedian(outside)


def _samples_above(skyVals, start, level):
    """*return the first and last index of the contiguous run of samples above ``level`` around ``start``*"""
    left = start
    while left > 0 and skyVals[left - 1] > level:
        left -= 1
    right = start
    while right < len(skyVals) - 1 and skyVals[right + 1] > level:
        right += 1
    return left, right
