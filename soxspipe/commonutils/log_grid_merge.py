#!/usr/bin/env python
"""
*Merge the extracted echelle orders of a SOXS arm onto its fixed log-wavelength grid*

Each order is cut at its join with the next order, and the join is snapped to the nearest bin edge, so no bin takes
flux from two orders and there is no overlap-weighted sum between orders. Within its cut, an order gives every bin
that its pixels cover wholly. A bin that the order covers only in part at an end is dropped. A bin that no pixel
overlaps, such as a bin in the gap between two orders, has no data: it holds 0 and is flagged, and nothing is
interpolated across it. See ``log_wavelength_grid`` for the grid and the resampling weights.

The module holds no recipe state: every function takes and returns numbers and arrays.

:Author:
    David Young

:Date Created:
    November 27, 2026
"""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from soxspipe.commonutils.log_wavelength_grid import (
    bin_edges,
    covered_bin_range,
    grid_wavelength,
    overlap_weights,
    snap_to_edge_index,
)


@dataclass(frozen=True)
class OrderPixels:
    """*the extracted pixels of one order*

    **Key Arguments:**

    - ``wavelength`` -- the wavelength at the centre of each pixel (nm)
    - ``pixelScale`` -- the width of each pixel in wavelength (nm)
    - ``keep`` -- *True* for each pixel to use. A pixel left out leaves a hole, which other pixels do not fill
    - ``values`` -- the density columns to resample, by name, one value per pixel
    - ``variance`` -- the variance of the flux of each pixel
    """

    wavelength: np.ndarray
    pixelScale: np.ndarray
    keep: np.ndarray
    values: Mapping[str, np.ndarray]
    variance: np.ndarray


@dataclass(frozen=True)
class MergedOnGrid:
    """*the merged spectrum on the log-wavelength grid*

    **Key Arguments:**

    - ``firstIndex`` -- the grid index of the first bin
    - ``wavelength`` -- the grid point of each bin (nm)
    - ``hasData`` -- *True* for each bin that at least one pixel overlaps
    - ``values`` -- the resampled density columns by name, 0 in a bin with no data
    - ``variance`` -- the resampled variance, 0 in a bin with no data
    - ``joinWavelengths`` -- the wavelength (nm) of each join after it was snapped to a bin edge, from red to blue
    """

    firstIndex: int
    wavelength: np.ndarray
    hasData: np.ndarray
    values: Mapping[str, np.ndarray]
    variance: np.ndarray
    joinWavelengths: tuple[float, ...]


def _order_bin_range(
    order: OrderPixels, lowerCutEdge: int | None, upperCutEdge: int | None, stepKms: float
) -> tuple[int, int]:
    """*the first and last bin an order gives: the bins it covers wholly and that lie between its two cuts*

    **Key Arguments:**

    - ``order`` -- the pixels of the order
    - ``lowerCutEdge`` -- the edge index of the cut on the blue side, or *None* for the bluest order
    - ``upperCutEdge`` -- the edge index of the cut on the red side, or *None* for the reddest order
    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``firstIndex``, ``lastIndex`` -- the grid indices of the first and last bin. ``lastIndex`` is below
      ``firstIndex`` when the order gives no bin
    """
    halfWidth = np.asarray(order.pixelScale, dtype=float) / 2.0
    lower = np.asarray(order.wavelength, dtype=float) - halfWidth
    upper = np.asarray(order.wavelength, dtype=float) + halfWidth
    isUsable = np.isfinite(lower) & np.isfinite(upper) & (upper > lower)
    if not isUsable.any():
        return 0, -1
    firstIndex, lastIndex = covered_bin_range(lower[isUsable].min(), upper[isUsable].max(), stepKms)
    # EDGE m IS THE LOWER EDGE OF BIN m, SO A BIN LIES BELOW THE CUT AT EDGE m WHEN ITS INDEX IS BELOW m
    if lowerCutEdge is not None:
        firstIndex = max(firstIndex, lowerCutEdge)
    if upperCutEdge is not None:
        lastIndex = min(lastIndex, upperCutEdge - 1)
    return firstIndex, lastIndex


def merge_orders_on_log_grid(
    orders: Sequence[OrderPixels], joinWavelengths: Sequence[float], stepKms: float
) -> MergedOnGrid:
    """*resample the orders onto the log-wavelength grid, each order cut at its joins*

    **Key Arguments:**

    - ``orders`` -- the orders from the reddest to the bluest
    - ``joinWavelengths`` -- the wavelength of the join of each pair of neighbouring orders (nm), red pair first
    - ``stepKms`` -- the step of the grid in km/s

    **Return:**

    - ``merged`` -- the ``MergedOnGrid`` spectrum. Its grid is unbroken from the first to the last bin that holds data

    **Raises:**

    - ``ValueError`` -- if there is not one join for each pair of neighbouring orders, or if no bin of the grid lies
      wholly inside the pixels of an order
    """
    if len(joinWavelengths) != max(len(orders) - 1, 0):
        raise ValueError(
            f"Expected one join for each pair of neighbouring orders ({max(len(orders) - 1, 0)}), "
            f"got {len(joinWavelengths)}."
        )
    joinEdges = [snap_to_edge_index(wavelength, stepKms) for wavelength in joinWavelengths]

    segments = []
    for position, order in enumerate(orders):
        upperCut = joinEdges[position - 1] if position > 0 else None
        lowerCut = joinEdges[position] if position < len(joinEdges) else None
        firstIndex, lastIndex = _order_bin_range(order, lowerCut, upperCut, stepKms)
        if lastIndex >= firstIndex:
            segments.append((order, firstIndex, lastIndex))
    if not segments:
        raise ValueError("No bin of the log-wavelength grid lies inside the extracted wavelength range.")

    firstIndex = min(segment[1] for segment in segments)
    lastIndex = max(segment[2] for segment in segments)
    binCount = lastIndex - firstIndex + 1
    hasData = np.zeros(binCount, dtype=bool)
    variance = np.zeros(binCount)
    values = {name: np.zeros(binCount) for name in segments[0][0].values}

    for order, segmentFirst, segmentLast in segments:
        halfWidth = np.asarray(order.pixelScale, dtype=float) / 2.0
        wavelength = np.asarray(order.wavelength, dtype=float)
        weights = overlap_weights(
            wavelength - halfWidth, wavelength + halfWidth, bin_edges(segmentFirst, segmentLast, stepKms), order.keep
        )
        window = slice(segmentFirst - firstIndex, segmentLast - firstIndex + 1)
        fill = weights.hasData & ~hasData[window]
        for name, column in order.values.items():
            values[name][window][fill] = weights.mean(column)[fill]
        variance[window][fill] = weights.propagate_variance(order.variance)[fill]
        hasData[window] |= fill

    return MergedOnGrid(
        firstIndex=firstIndex,
        wavelength=grid_wavelength(firstIndex + np.arange(binCount), stepKms),
        hasData=hasData,
        values=values,
        variance=variance,
        joinWavelengths=tuple(float(grid_wavelength(edge - 0.5, stepKms)) for edge in joinEdges),
    )
