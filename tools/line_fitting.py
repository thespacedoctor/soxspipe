"""*find and fit isolated emission lines on a frame that has a wavelength, slit and order map*

This is the line measurement of the DY-2297 prototype `linewidth.py`, used by `measure_spec_res.py` on slit arcs and on
sky lines. It works on arrays only: the caller reads the frame and the `DISP_IMAGE` maps.

**Method**

1. Per order, build a median spectrum on half-pixel wavelength bins from the pixels inside the slit range.
2. Find the peaks of that spectrum. Keep the isolated ones, with no other peak within 8 pixels.
3. Fit a Gaussian to the pixels within 5 pixels of each kept peak. Leave out a line with a flagged or saturated pixel in
   the window, a failed fit, a width outside 0.8 to 8 pixels or a centre that moved more than 2 pixels.
4. The resolving power of a line is `centre / FWHM`, in nm.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

FWHM_PER_SIGMA = 2.3548
ARC_SLIT_RANGE_ARCSEC = (-3.0, 3.0)
MIN_ORDER_PIXELS = 500
PROFILE_BIN_FRACTION = 0.5
MIN_PEAK_DISTANCE_BINS = 6
ISOLATION_PIXELS = 8
FIT_WINDOW_PIXELS = 5
MIN_FIT_PIXELS = 20
MIN_FWHM_PIXELS = 0.8
MAX_FWHM_PIXELS = 8.0
MAX_CENTRE_DRIFT_PIXELS = 2.0
ARC_PROMINENCE = 30.0
MAD_TO_SIGMA = 1.4826
# THE KEYS ARE THE DETECTOR PARAMETER dispersion-axis OF THE ARM (VIS x, NIR y). THE VALUES ARE THE NUMPY AXIS OF THE
# DISP_IMAGE MAPS: THE WAVELENGTH RUNS ALONG AXIS 0 FOR VIS AND AXIS 1 FOR NIR (DY-2455). SO THE KEY "x" IS NUMPY AXIS 0
DISPERSION_AXIS_INDEX = {"x": 0, "y": 1}

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class LineFit:
    """*one fitted emission line: its order, centre, widths and resolving power*"""

    order: int
    centreNm: float
    fwhmNm: float
    fwhmPx: float
    resolvingPower: float


def gaussian(x: FloatArray, amplitude: float, centre: float, sigma: float, offset: float) -> FloatArray:
    """*a Gaussian with a constant offset*

    **Key Arguments:**

    - ``x`` -- the positions to evaluate at
    - ``amplitude`` -- the peak height above the offset
    - ``centre`` -- the position of the peak
    - ``sigma`` -- the standard deviation
    - ``offset`` -- the constant level

    **Return:**

    - ``values`` -- the Gaussian at ``x``

    **Usage:**

    ```python
    values = gaussian(wavelengths, 5000.0, 550.0, 0.04, 100.0)
    ```
    """
    return amplitude * np.exp(-0.5 * ((x - centre) / sigma) ** 2) + offset


def nm_per_pixel(wavelength: FloatArray, mask: npt.NDArray[np.bool_], axis: int) -> float:
    """*the median wavelength step between neighbouring pixels along the dispersion axis, inside a mask*

    **Key Arguments:**

    - ``wavelength`` -- the wavelength map
    - ``mask`` -- True for the pixels to use
    - ``axis`` -- the array axis the wavelength runs along

    **Return:**

    - ``nmStep`` -- the median absolute step between neighbours that are both in the mask, in nm

    **Usage:**

    ```python
    nmStep = nm_per_pixel(wavelength, mask, axis=1)
    ```
    """
    lower = tuple(slice(None, -1) if index == axis else slice(None) for index in range(2))
    upper = tuple(slice(1, None) if index == axis else slice(None) for index in range(2))
    both = mask[lower] & mask[upper]
    return float(np.median(np.abs(wavelength[upper] - wavelength[lower])[both]))


def median_profile(wavelengths: FloatArray, values: FloatArray, binWidth: float) -> tuple[FloatArray, FloatArray]:
    """*the median of the values in wavelength bins of ``binWidth``, and the bin edges*

    **Key Arguments:**

    - ``wavelengths`` -- the pixel wavelengths
    - ``values`` -- the pixel values
    - ``binWidth`` -- the bin width in nm

    **Return:**

    - ``profile`` -- the median per bin, with empty bins filled with the overall median
    - ``edges`` -- the left edge of each bin

    **Usage:**

    ```python
    profile, edges = median_profile(wavelengths, values, 0.5 * nmStep)
    ```
    """
    edges = np.arange(wavelengths.min(), wavelengths.max(), binWidth)
    binIndex = np.digitize(wavelengths, edges)
    medians = pd.Series(values).groupby(binIndex).median()
    profile = medians.reindex(range(1, len(edges) + 1)).to_numpy(dtype=float)
    return np.nan_to_num(profile, nan=np.nanmedian(profile)), edges


def isolated_centres(
    profile: FloatArray, edges: FloatArray, binWidth: float, nmStep: float, prominence: float
) -> list[float]:
    """*the wavelengths of the profile's peaks that have no other peak within ``ISOLATION_PIXELS`` pixels*

    **Key Arguments:**

    - ``profile`` -- the median spectrum
    - ``edges`` -- the left edge of each profile bin
    - ``binWidth`` -- the bin width in nm
    - ``nmStep`` -- the nm per pixel
    - ``prominence`` -- the peak prominence, in units of the profile's noise

    **Return:**

    - ``centres`` -- the wavelengths of the isolated peaks

    **Usage:**

    ```python
    centres = isolated_centres(profile, edges, binWidth, nmStep, 30.0)
    ```
    """
    base = np.nanmedian(profile)
    noise = MAD_TO_SIGMA * np.nanmedian(np.abs(profile - base))
    peaks, _ = find_peaks(profile, prominence=prominence * noise, distance=MIN_PEAK_DISTANCE_BINS)
    centres = edges[peaks] + binWidth / 2
    return [
        float(centre)
        for index, centre in enumerate(centres)
        if not np.any(np.abs(np.delete(centres, index) - centre) < ISOLATION_PIXELS * nmStep)
    ]


def fit_line(
    wavelengths: FloatArray,
    values: FloatArray,
    isGood: npt.NDArray[np.bool_],
    centre: float,
    nmStep: float,
    order: int,
) -> LineFit | None:
    """*fit a Gaussian to the pixels around ``centre``; None when the line cannot be used*

    **Key Arguments:**

    - ``wavelengths`` -- the sorted pixel wavelengths of the order
    - ``values`` -- the pixel values
    - ``isGood`` -- True where the pixel is usable
    - ``centre`` -- the wavelength of the peak
    - ``nmStep`` -- the nm per pixel
    - ``order`` -- the order number, recorded on the result

    **Return:**

    - ``line`` -- the fitted line, or *None* when the window holds a flagged pixel, too few pixels, a failed fit, a
      negative amplitude, a width outside 0.8 to 8 pixels or a centre that moved more than 2 pixels

    **Usage:**

    ```python
    line = fit_line(wavelengths, values, isGood, 550.0, 0.03, 1)
    ```
    """
    window = np.abs(wavelengths - centre) < FIT_WINDOW_PIXELS * nmStep
    if window.sum() < MIN_FIT_PIXELS or not isGood[window].all():
        return None
    windowWavelengths, windowValues = wavelengths[window], values[window]
    start = [windowValues.max() - np.median(windowValues), centre, nmStep, np.median(windowValues)]
    try:
        params, _ = curve_fit(gaussian, windowWavelengths, windowValues, p0=start)
    except (RuntimeError, ValueError):
        return None
    fwhm = FWHM_PER_SIGMA * abs(params[2])
    isWidthOk = MIN_FWHM_PIXELS * nmStep < fwhm < MAX_FWHM_PIXELS * nmStep
    if not isWidthOk or params[0] <= 0 or abs(params[1] - centre) > MAX_CENTRE_DRIFT_PIXELS * nmStep:
        return None
    return LineFit(order, float(params[1]), float(fwhm), float(fwhm / nmStep), float(params[1] / fwhm))


def measure_order(
    data: FloatArray,
    isGood: npt.NDArray[np.bool_],
    wavelength: FloatArray,
    mask: npt.NDArray[np.bool_],
    axis: int,
    order: int,
    prominence: float,
) -> list[LineFit]:
    """*the fitted lines of one order*

    **Key Arguments:**

    - ``data`` -- the frame
    - ``isGood`` -- True where the pixel is usable
    - ``wavelength`` -- the wavelength map
    - ``mask`` -- True for the pixels of this order inside the slit range
    - ``axis`` -- the array axis the wavelength runs along
    - ``order`` -- the order number
    - ``prominence`` -- the peak prominence, in units of the profile's noise

    **Return:**

    - ``lines`` -- the lines kept

    **Usage:**

    ```python
    lines = measure_order(data, isGood, wavelength, mask, 1, 3, 30.0)
    ```
    """
    nmStep = nm_per_pixel(wavelength, mask, axis)
    sortIndex = np.argsort(wavelength[mask])
    sortedWavelengths = wavelength[mask][sortIndex]
    sortedValues = data[mask][sortIndex]
    sortedGood = isGood[mask][sortIndex]
    binWidth = nmStep * PROFILE_BIN_FRACTION
    profile, edges = median_profile(sortedWavelengths, sortedValues, binWidth)
    centres = isolated_centres(profile, edges, binWidth, nmStep, prominence)
    fitted = (fit_line(sortedWavelengths, sortedValues, sortedGood, centre, nmStep, order) for centre in centres)
    return [line for line in fitted if line is not None]


def measure_lines(
    data: FloatArray,
    good: npt.NDArray[np.bool_],
    wavelength: FloatArray,
    slit: FloatArray,
    order: npt.NDArray[np.integer],
    dispersionAxis: str,
    slitRange: tuple[float, float] = ARC_SLIT_RANGE_ARCSEC,
    prominence: float = ARC_PROMINENCE,
) -> list[LineFit]:
    """*the isolated, unsaturated emission lines of a frame, with their resolving power*

    **Key Arguments:**

    - ``data`` -- the frame, trimmed to the shape of the maps
    - ``good`` -- True where the pixel is usable (not flagged, not saturated)
    - ``wavelength`` -- the ``WAVELENGTH`` map of the ``DISP_IMAGE``, in nm
    - ``slit`` -- the ``SLIT`` map, in arcsec
    - ``order`` -- the ``ORDER`` map; 0 is no order
    - ``dispersionAxis`` -- ``x`` or ``y``, the ``dispersion-axis`` of the arm in the detector parameters. ``x`` is
      numpy axis 0 of the maps (VIS) and ``y`` is numpy axis 1 (NIR)
    - ``slitRange`` -- the open interval of slit positions to use, in arcsec. Default ``ARC_SLIT_RANGE_ARCSEC``, within
      3 arcsec of the centre
    - ``prominence`` -- the peak prominence, in units of the profile's noise. Default ``ARC_PROMINENCE``, 30

    **Return:**

    - ``lines`` -- one ``LineFit`` per kept line

    **Usage:**

    ```python
    lines = measure_lines(data, good, wavelength, slit, order, "x")
    ```
    """
    axis = DISPERSION_AXIS_INDEX[dispersionAxis]
    inSlit = (slit > slitRange[0]) & (slit < slitRange[1]) & np.isfinite(data) & (wavelength > 0)
    lines: list[LineFit] = []
    for orderNumber in np.unique(order[order > 0]):
        mask = (order == orderNumber) & inSlit
        if mask.sum() >= MIN_ORDER_PIXELS:
            lines.extend(measure_order(data, good, wavelength, mask, axis, int(orderNumber), prominence))
    return lines


def median_resolving_power(lines: list[LineFit]) -> float:
    """*the median resolving power over fitted lines*

    **Key Arguments:**

    - ``lines`` -- the fitted lines, at least one

    **Return:**

    - ``resolvingPower`` -- the median of ``centre / FWHM``

    **Usage:**

    ```python
    resolvingPower = median_resolving_power(lines)
    ```
    """
    return float(np.median([line.resolvingPower for line in lines]))
