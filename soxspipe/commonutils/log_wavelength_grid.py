#!/usr/bin/env python
"""
*The fixed log-wavelength grid of a SOXS arm and the flux-conserving resampler that moves a spectrum onto it*

Every point of the grid is :math:`\\lambda_k = \\lambda_0 \\exp(k \\Delta v / c)` with :math:`\\lambda_0 = 1000` nm,
so every bin has the same width in km/s and the grid is the same for every frame. A grid point is the geometric centre
of its bin: bin ``k`` spans the edges ``k - 1/2`` and ``k + 1/2`` (edge ``m`` sits at
:math:`\\lambda_0 \\exp((m - 1/2) \\Delta v / c)`).

Every column the merge resamples is a density (flat-fielded electrons per pixel, the flat has already removed the pixel
scale), so a bin takes the weighted mean of the pixels that overlap it. The weights are the overlaps, normalized to 1
over the bin. The variance of a bin is the sum of the squared weights times the variances of the pixels (Sum w^2
sigma^2), which ignores the covariance between neighbouring bins. A bin that no pixel overlaps has no data: the
resampler returns 0 for it, reports it in ``OverlapWeights.hasData``, and never interpolates across it.

The module holds no recipe state: every function takes and returns numbers and arrays.

:Author:
    David Young

:Date Created:
    November 27, 2026
"""

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

# THE WAVELENGTH OF GRID POINT 0 (nm), THE SAME FOR EVERY ARM
LOG_GRID_ORIGIN_NM = 1000.0
SPEED_OF_LIGHT_KMS = 299792.458

# THE YAML KEY THAT HOLDS THE STEP OF AN ARM PER UNBINNED PIXEL (km/s)
LOG_STEP_KEY = "log-wavelength-step-kms"

# ROUNDING NOISE TOLERATED WHEN A WAVELENGTH THAT SITS ON AN EDGE IS CONVERTED TO AN EDGE INDEX (IN BINS)
EDGE_INDEX_TOLERANCE = 1.0e-9


def grid_step_kms(detectorParams: Mapping[str, Any], binning: int) -> float:
    """*the log-wavelength step of the grid of an arm (km/s), for the spectral binning of a frame*

    **Key Arguments:**

    - ``detectorParams`` -- the detector parameters of the arm, holding ``log-wavelength-step-kms`` per unbinned pixel
    - ``binning`` -- the detector binning along the dispersion axis, 1 or more

    **Return:**

    - ``stepKms`` -- the step per bin in km/s

    **Raises:**

    - ``ValueError`` -- if the detector parameters hold no step, or a step that is not a positive finite number
    """
    step = detectorParams.get(LOG_STEP_KEY)
    isNumber = isinstance(step, (int, float)) and not isinstance(step, bool)
    if not isNumber or not math.isfinite(step) or step <= 0:
        raise ValueError(
            f"The detector parameters need a positive number for '{LOG_STEP_KEY}' (km/s per unbinned pixel), "
            f"found {step!r}."
        )
    return float(step) * binning


def _velocity_fraction(stepKms: float) -> float:
    """*the step as a fraction of the speed of light*

    **Key Arguments:**

    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``fraction`` -- the step divided by the speed of light
    """
    return stepKms / SPEED_OF_LIGHT_KMS


def grid_wavelength(index: Any, stepKms: float) -> np.ndarray:
    """*the wavelength (nm) of grid point ``index``*

    **Key Arguments:**

    - ``index`` -- the integer grid index, or an array of them
    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``wavelength`` -- ``1000 * exp(index * step / c)`` in nm
    """
    return LOG_GRID_ORIGIN_NM * np.exp(np.asarray(index, dtype=float) * _velocity_fraction(stepKms))


def bin_edges(firstIndex: int, lastIndex: int, stepKms: float) -> np.ndarray:
    """*the edges (nm) of the bins ``firstIndex`` to ``lastIndex``, both included*

    **Key Arguments:**

    - ``firstIndex`` -- the grid index of the first bin
    - ``lastIndex`` -- the grid index of the last bin
    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``edges`` -- the ``lastIndex - firstIndex + 2`` edges in nm, rising
    """
    return grid_wavelength(np.arange(firstIndex, lastIndex + 2) - 0.5, stepKms)


def _edge_position(wavelength: Any, stepKms: float) -> np.ndarray:
    """*the wavelength as a fractional edge index: an integer where the wavelength is on a bin edge*

    **Key Arguments:**

    - ``wavelength`` -- the wavelength in nm, or an array of them
    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``position`` -- the fractional edge index
    """
    return np.log(np.asarray(wavelength, dtype=float) / LOG_GRID_ORIGIN_NM) / _velocity_fraction(stepKms) + 0.5


def snap_to_edge_index(wavelength: float, stepKms: float) -> int:
    """*the index of the bin edge nearest to a wavelength, measured in log wavelength*

    Edge ``m`` is the lower edge of bin ``m`` and the upper edge of bin ``m - 1``.

    **Key Arguments:**

    - ``wavelength`` -- the wavelength in nm
    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``edgeIndex`` -- the index of the nearest edge
    """
    return int(np.rint(_edge_position(wavelength, stepKms)))


def covered_bin_range(lowerWavelength: float, upperWavelength: float, stepKms: float) -> tuple[int, int]:
    """*the first and last bin that lie wholly between two wavelengths*

    A bin that the wavelength range covers only in part is left out.

    **Key Arguments:**

    - ``lowerWavelength`` -- the shortest wavelength covered (nm)
    - ``upperWavelength`` -- the longest wavelength covered (nm)
    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``firstIndex``, ``lastIndex`` -- the grid indices of the first and last covered bin. ``lastIndex`` is below
      ``firstIndex`` when no bin is covered
    """
    firstIndex = math.ceil(_edge_position(lowerWavelength, stepKms) - EDGE_INDEX_TOLERANCE)
    lastIndex = math.floor(_edge_position(upperWavelength, stepKms) + EDGE_INDEX_TOLERANCE) - 1
    return int(firstIndex), int(lastIndex)


@dataclass(frozen=True)
class OverlapWeights:
    """*the normalized overlap of each pixel with each bin it touches*

    Each entry ``i`` says that pixel ``pixelIndex[i]`` overlaps bin ``binIndex[i]`` with the normalized weight
    ``weight[i]``. The weights of one bin add up to 1.

    **Key Arguments:**

    - ``binCount`` -- the number of bins
    - ``binIndex`` -- the bin of each entry
    - ``pixelIndex`` -- the pixel of each entry
    - ``weight`` -- the normalized overlap of each entry
    - ``hasData`` -- *True* for each bin that at least one kept pixel overlaps
    """

    binCount: int
    binIndex: np.ndarray
    pixelIndex: np.ndarray
    weight: np.ndarray
    hasData: np.ndarray

    def mean(self, values: np.ndarray) -> np.ndarray:
        """*the weighted mean of a per-pixel density in each bin, 0 where the bin has no data*

        **Key Arguments:**

        - ``values`` -- the value of each pixel

        **Return:**

        - ``resampled`` -- the value of each bin
        """
        values = np.asarray(values, dtype=float)
        return np.bincount(self.binIndex, weights=self.weight * values[self.pixelIndex], minlength=self.binCount)

    def propagate_variance(self, variance: np.ndarray) -> np.ndarray:
        """*the variance of each bin, the sum of the squared weights times the pixel variances*

        **Key Arguments:**

        - ``variance`` -- the variance of each pixel

        **Return:**

        - ``resampled`` -- the variance of each bin, 0 where the bin has no data
        """
        variance = np.asarray(variance, dtype=float)
        return np.bincount(self.binIndex, weights=self.weight**2 * variance[self.pixelIndex], minlength=self.binCount)


def usable_pixel_mask(pixelLower: np.ndarray, pixelUpper: np.ndarray, keep: np.ndarray | None = None) -> np.ndarray:
    """*flag the pixels that can contribute: finite edges, a positive width, and kept*

    **Key Arguments:**

    - ``pixelLower`` -- the lower edge of each pixel (nm)
    - ``pixelUpper`` -- the upper edge of each pixel (nm)
    - ``keep`` -- *True* for each pixel to use (default: all of them)

    **Return:**

    - ``isUsable`` -- boolean array, *True* for each usable pixel
    """
    isUsable = np.isfinite(pixelLower) & np.isfinite(pixelUpper) & (pixelUpper > pixelLower)
    if keep is not None:
        isUsable &= np.asarray(keep, dtype=bool)
    return isUsable


def overlap_weights(
    pixelLower: np.ndarray, pixelUpper: np.ndarray, binEdges: np.ndarray, keep: np.ndarray | None = None
) -> OverlapWeights:
    """*the normalized overlap weights that resample pixels onto bins*

    A pixel with a non-finite edge, or with no width, is dropped. A dropped pixel leaves a hole: the bins that only it
    overlapped have no data, and the other pixels do not fill it.

    **Key Arguments:**

    - ``pixelLower`` -- the lower edge of each pixel (nm)
    - ``pixelUpper`` -- the upper edge of each pixel (nm)
    - ``binEdges`` -- the rising edges of the bins (nm), one more than there are bins
    - ``keep`` -- *True* for each pixel to use (default: all of them)

    **Return:**

    - ``weights`` -- the ``OverlapWeights`` of the pixels on the bins
    """
    pixelLower = np.asarray(pixelLower, dtype=float)
    pixelUpper = np.asarray(pixelUpper, dtype=float)
    binEdges = np.asarray(binEdges, dtype=float)
    binCount = len(binEdges) - 1

    isUsable = usable_pixel_mask(pixelLower, pixelUpper, keep)
    lower = np.where(isUsable, pixelLower, 0.0)
    upper = np.where(isUsable, pixelUpper, 0.0)

    # THE BINS EACH PIXEL TOUCHES: FROM THE ONE HOLDING ITS LOWER EDGE TO THE ONE HOLDING ITS UPPER EDGE
    firstBin = np.clip(np.searchsorted(binEdges, lower, side="right") - 1, 0, binCount - 1)
    lastBin = np.clip(np.searchsorted(binEdges, upper, side="left") - 1, 0, binCount - 1)
    touched = np.where(isUsable, np.maximum(lastBin - firstBin + 1, 0), 0)
    pixelIndex = np.repeat(np.arange(len(lower)), touched)
    entryStart = np.cumsum(touched) - touched
    binIndex = np.repeat(firstBin, touched) + np.arange(touched.sum()) - np.repeat(entryStart, touched)

    overlap = np.minimum(upper[pixelIndex], binEdges[binIndex + 1]) - np.maximum(lower[pixelIndex], binEdges[binIndex])
    hasOverlap = overlap > 0
    pixelIndex, binIndex, overlap = pixelIndex[hasOverlap], binIndex[hasOverlap], overlap[hasOverlap]

    binTotal = np.bincount(binIndex, weights=overlap, minlength=binCount)
    return OverlapWeights(
        binCount=binCount,
        binIndex=binIndex,
        pixelIndex=pixelIndex,
        weight=overlap / binTotal[binIndex],
        hasData=binTotal > 0,
    )
