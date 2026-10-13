"""The both-lit overlap mask: the slit-end pixels where the science sky and the master flat are not both lit.

Spatial flexure shifts the science slit image along the slit by an offset δ per order (``slit_illumination``). The
master flat keeps its own lit edges, so at one slit end the science is lit where the flat falls away, and at the other
the flat is lit where the science is dark. Neither end can be fixed by a coordinate shift, so it is excluded from the
sky model instead. For each order the kept range is [max(l, l+δ) + m, min(u, u+δ) − m], where l and u are the flat's
lit edges and m is a margin. This module is pure numpy and returns new objects without modifying its inputs.
"""

import math
from numbers import Real

import numpy as np

from soxspipe.commonutils.slit_illumination import STATUS_BORROWED, STATUS_FITTED
from soxspipe.commonutils.slit_illumination_profiles import in_order_pixels

__all__ = ["overlap_exclude_mask"]


def _usable_offsets(offsets):
    """*``{order: offset}`` for every order whose offset was fitted or borrowed and is finite*"""
    isUsable = offsets["status"].isin([STATUS_FITTED, STATUS_BORROWED]) & np.isfinite(
        offsets["offset"].to_numpy(dtype=float)
    )
    usable = offsets.loc[isUsable]
    return dict(zip(usable["order"].astype(int), usable["offset"].astype(float), strict=True))


def _usable_flat_edges(flatEdges):
    """*``{order: (lower, upper)}`` for every order whose flat edges are finite*"""
    lower = flatEdges["lower"].to_numpy(dtype=float)
    upper = flatEdges["upper"].to_numpy(dtype=float)
    isUsable = np.isfinite(lower) & np.isfinite(upper)
    orders = flatEdges["order"].to_numpy()[isUsable]
    return {int(order): (lo, up) for order, lo, up in zip(orders, lower[isUsable], upper[isUsable], strict=True)}


def overlap_exclude_mask(log, slit, orderMap, flatEdges, offsets, marginArcsec):
    """*the pixels of each order that lie outside the range where both the science sky and the master flat are lit*

    An order is masked only if it has a usable offset (fitted or borrowed) and finite flat edges. Any other order is
    left unmasked, as is an order whose kept range is empty, which is logged as a warning. Only pixels that belong to
    a whole-numbered order and have a finite slit position can be excluded, so inter-order pixels never are. A pixel
    exactly on a limit of the kept range is kept.

    **Key Arguments:**

    - ``log`` -- logger
    - ``slit`` -- the 2D map's SLIT image in arcsec, at the frame's binning
    - ``orderMap`` -- the 2D map's ORDER image, the same shape as ``slit``
    - ``flatEdges`` -- the per-order flat edges from ``slit_illumination.flat_slit_edges``
    - ``offsets`` -- the per-order offsets from ``slit_illumination.fit_order_offsets``
    - ``marginArcsec`` -- how far inside the both-lit range the mask edges sit, in arcsec; finite and not negative

    **Return:**

    - ``excludeMask`` -- boolean image, the shape of ``slit``, True for a pixel to exclude from the sky model

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination_mask import overlap_exclude_mask
    excludeMask = overlap_exclude_mask(
        log, slit, orderMap, result.flatEdges, result.offsets, result.overlapMarginArcsec
    )
    ```
    """
    log.debug("starting the ``overlap_exclude_mask`` function")
    slit = np.asarray(slit, dtype=float)
    orderMap = np.asarray(orderMap, dtype=float)
    if slit.shape != orderMap.shape:
        raise ValueError(f"image shapes differ: SLIT {slit.shape}, ORDER {orderMap.shape}")
    isNumber = isinstance(marginArcsec, Real) and not isinstance(marginArcsec, bool)
    if not (isNumber and math.isfinite(marginArcsec) and marginArcsec >= 0):
        raise ValueError(f"the overlap margin must be a finite number of arcsec, not negative: {marginArcsec!r}")

    excludeMask = np.zeros(slit.shape, dtype=bool)
    isInOrder = in_order_pixels(orderMap, slit)
    flatByOrder = _usable_flat_edges(flatEdges)
    for order, offset in _usable_offsets(offsets).items():
        if order not in flatByOrder:
            continue
        lower, upper = flatByOrder[order]
        keepLower = max(lower, lower + offset) + marginArcsec
        keepUpper = min(upper, upper + offset) - marginArcsec
        if not keepLower < keepUpper:
            log.warning(
                f"order {order} is not masked: its both-lit slit range "
                f"[{keepLower:.2f}, {keepUpper:.2f}] arcsec is empty"
            )
            continue
        isOrder = isInOrder & (orderMap == order)
        # THE NAN SLIT OF A NON-ORDER PIXEL MAKES BOTH COMPARISONS FALSE, AND isOrder EXCLUDES IT ANYWAY
        with np.errstate(invalid="ignore"):
            excludeMask |= isOrder & ((slit < keepLower) | (slit > keepUpper))

    log.debug(f"the overlap mask excludes {np.count_nonzero(excludeMask)} pixels")
    log.debug("completed the ``overlap_exclude_mask`` function")
    return excludeMask
