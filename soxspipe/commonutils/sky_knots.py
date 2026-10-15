"""B-spline knot placement for the stare sky fit: the native wavelength step and uniform knots (DY-1287).

Kelson (2003) places the knots of a sky b-spline at a fixed spacing of about one native pixel in wavelength. The tilted
slit samples the sky more finely than that, so every knot interval holds many samples. This module is pure numpy and
returns new arrays without modifying its inputs.
"""

import numpy as np


def native_wavelength_step(wavelength, dispersionPixel, crossPixel):
    """*the median wavelength step between pixels that neighbour each other along the dispersion axis*

    **Key Arguments:**

    - ``wavelength`` -- array of the wavelength of each pixel, in any order
    - ``dispersionPixel`` -- array of the integer pixel position of each pixel along the dispersion axis
    - ``crossPixel`` -- array of the integer pixel position of each pixel along the cross-dispersion axis

    **Return:**

    - ``step`` -- the median absolute wavelength difference between two pixels on the same cross-dispersion
      line whose dispersion positions differ by exactly one

    **Raises:**

    - ``ValueError`` -- if no two pixels neighbour each other, or the median step is not finite and positive

    **Usage:**

    ```python
    from soxspipe.commonutils.sky_knots import native_wavelength_step
    step = native_wavelength_step(df["wavelength"].values, df["x"].values, df["y"].values)
    ```
    """
    wavelength = np.asarray(wavelength, dtype=float)
    dispersionPixel = np.asarray(dispersionPixel)
    crossPixel = np.asarray(crossPixel)
    ordering = np.lexsort((dispersionPixel, crossPixel))
    wavelength = wavelength[ordering]
    dispersionPixel = dispersionPixel[ordering]
    crossPixel = crossPixel[ordering]

    isNeighbour = (np.diff(crossPixel) == 0) & (np.diff(dispersionPixel) == 1)
    steps = np.abs(np.diff(wavelength))[isNeighbour]
    steps = steps[np.isfinite(steps)]
    if steps.size == 0:
        raise ValueError("no two pixels are neighbouring along the dispersion axis, so there is no wavelength step")
    step = float(np.median(steps))
    if not step > 0:
        raise ValueError(f"the native wavelength step must be positive, got {step}")
    return step


def uniform_interior_knots(wavelengthMin, wavelengthMax, knotSpacing):
    """*evenly spaced interior knots, starting one spacing above the bluest wavelength*

    **Key Arguments:**

    - ``wavelengthMin`` -- the bluest sample wavelength
    - ``wavelengthMax`` -- the reddest sample wavelength
    - ``knotSpacing`` -- the wavelength between two knots, finite and positive

    **Return:**

    - ``knots`` -- array of the knots strictly between ``wavelengthMin`` and ``wavelengthMax``

    **Raises:**

    - ``ValueError`` -- if ``knotSpacing`` is not finite and positive
    """
    if not (np.isfinite(knotSpacing) and knotSpacing > 0):
        raise ValueError(f"the knot spacing must be finite and positive, got {knotSpacing}")
    knotCount = int(np.ceil((wavelengthMax - wavelengthMin) / knotSpacing))
    knots = wavelengthMin + knotSpacing * np.arange(1, max(knotCount, 1))
    return knots[knots < wavelengthMax]


def knots_with_samples(knots, wavelength, minSamples=1):
    """*drop the knots that would leave a knot interval with fewer than ``minSamples`` samples*

    The knots are taken from blue to red. A knot is kept when at least ``minSamples`` samples lie between the
    last kept knot (or the blue end) and it. Then the reddest kept knots are dropped until the red end interval
    also holds ``minSamples`` samples. Dropping a knot merges its two intervals.

    **Key Arguments:**

    - ``knots`` -- the proposed interior knots, in any order
    - ``wavelength`` -- the sorted sample wavelengths
    - ``minSamples`` -- the fewest samples a knot interval may hold. Default *1*

    **Return:**

    - ``keptKnots`` -- a new sorted array of the kept knots
        - of two knots too close together, the bluer knot is kept
        - a sample exactly on a knot counts for neither neighbouring interval
    """
    wavelength = np.asarray(wavelength, dtype=float)
    proposedKnots = np.sort(np.asarray(knots, dtype=float))
    keptKnots = []
    if len(wavelength) == 0:
        return np.array(keptKnots, dtype=float)

    lastWavelength = wavelength[-1]
    firstSampleAbove = 0
    for knot in proposedKnots:
        if knot >= lastWavelength:
            break
        if np.searchsorted(wavelength, knot, side="left") - firstSampleAbove >= minSamples:
            keptKnots.append(knot)
            firstSampleAbove = np.searchsorted(wavelength, knot, side="right")

    while keptKnots and len(wavelength) - np.searchsorted(wavelength, keptKnots[-1], side="right") < minSamples:
        keptKnots.pop()
    return np.array(keptKnots, dtype=float)
