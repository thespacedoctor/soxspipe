"""Slit coverage: the pixels of the slit that not every frame of a nod position covers (DY-2290).

The frames of one nod position are aligned to each other along the slit (see ``slit_alignment``). The shift of a frame
can leave slit pixels that the frame does not reach, so a stack of the frames holds fewer frames there. This module
finds those pixels, so the sky fit and the extraction can leave them out.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from soxspipe.commonutils.slit_alignment import (
    _frame_names,
    _is_source_on_order,
    _order_shifts,
    slit_alignment_result,
    two_d_map_images,
)
from soxspipe.commonutils.slit_alignment_sources import _order_sources


def alignment_coverage_mask(frameShifts: Sequence[Mapping[int, float]], mapImages: two_d_map_images) -> np.ndarray:
    """*the pixels of the slit that not every frame of one nod position covers*

    The frames of one nod position are aligned to each other along the slit, so a frame shifted by ``s`` arcsec takes
    its value in order ``o`` from where the slit coordinate is ``s`` arcsec higher at the same wavelength (see
    ``shift_frame_along_slit``). A pixel is uncovered when that source lies off the detector or outside the order. A
    pixel is covered by every frame when it is covered at the lowest and at the highest shift, so those two set the
    uncovered ends of an order. This is the slit range that the stack does not reach in every frame, not a flag for a
    bad pixel. Pixels outside every order are never uncovered, and a shift longer than the slit leaves the whole order
    uncovered.

    **Key Arguments:**

    - ``frameShifts`` -- for each frame of the position, its ``{order: shift in arcsec}`` from the alignment to the
      other frames. The reference frame has no shift. An empty sequence is one frame that needs no alignment
    - ``mapImages`` -- the ``two_d_map_images`` of the frame's binning

    **Return:**

    - ``uncovered`` -- boolean image, True for a pixel that a frame does not cover

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_coverage import alignment_coverage_mask, orders_without_coverage

    uncovered = alignment_coverage_mask(frame_order_shifts(result, names), result.mapImages)
    emptyOrders = orders_without_coverage(uncovered, result.mapImages)
    ```
    """
    frames = list(frameShifts) or [{}]
    lowest, highest = {}, {}
    for order in mapImages.orders:
        totals = [shifts.get(order, 0.0) for shifts in frames]
        lowest[order], highest[order] = min(totals), max(totals)
    return _uncovered_by_shifts(lowest, mapImages) | _uncovered_by_shifts(highest, mapImages)


def orders_without_coverage(uncovered: np.ndarray, mapImages: two_d_map_images) -> list[int]:
    """*the orders whose every pixel is uncovered, so they have no slit range that every frame covers*

    **Key Arguments:**

    - ``uncovered`` -- the boolean image of ``alignment_coverage_mask``
    - ``mapImages`` -- the ``two_d_map_images`` of the frame's binning

    **Return:**

    - ``orders`` -- the order numbers that have no covered pixel

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_coverage import orders_without_coverage

    emptyOrders = orders_without_coverage(uncovered, mapImages)
    ```
    """
    return [order for order in mapImages.orders if uncovered[mapImages.order == order].all()]


def frame_order_shifts(
    result: slit_alignment_result, frameNames: Sequence[str] | None = None
) -> list[dict[int, float]]:
    """*the alignment shift of each order of each frame of an aligned group*

    **Key Arguments:**

    - ``result`` -- the ``slit_alignment_result`` of the group
    - ``frameNames`` -- the names the group was aligned with. Default *None*, for the names it was given by default

    **Return:**

    - ``shifts`` -- for each frame in the order of ``result.frames``, a ``{order: shift in arcsec}`` dictionary

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_coverage import frame_order_shifts

    shiftsOfEachFrame = frame_order_shifts(results[0], namesA)
    ```
    """
    orders = result.mapImages.orders if result.mapImages is not None else ()
    names = _frame_names(result.frames, frameNames)
    return [_order_shifts(result.shifts, name, orders) for name in names]


def _uncovered_by_shifts(orderShifts: Mapping[int, float], mapImages: two_d_map_images) -> np.ndarray:
    """*the pixels with no source on the order for one slit shift per order*"""
    uncovered = np.zeros(mapImages.order.shape, dtype=bool)
    for order, shift in orderShifts.items():
        if shift == 0:
            continue
        targets, sources = _order_sources(mapImages, order, shift)
        uncovered[targets] = ~_is_source_on_order(mapImages, order, sources)
    return uncovered
