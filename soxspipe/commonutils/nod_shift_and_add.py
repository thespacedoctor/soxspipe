"""Shift-and-add of a nod pair: the B-A difference moves onto the A-B difference along the slit, and the two average.

The B-A difference is exactly the negative of the A-B difference, so its positive trace is the object at nod B, and
its negative trace is the object at nod A. Resampling it by the B minus A slit position puts its positive trace on the
positive trace of A-B. The mean of the two frames then holds the object once, with the noise of two exposures
averaged, so one extraction replaces the two the pair needed before.
"""

import copy
from dataclasses import dataclass

import numpy as np

from soxspipe.commonutils.slit_alignment import WHOLE_FRAME

NO_PLAN_MESSAGE = (
    "The slit position of the object could not be measured at both nod positions, so the A-B and B-A differences "
    "are extracted separately instead of being combined."
)


@dataclass(frozen=True)
class shift_and_add_plan:
    """The slit shift of each order from the B-A difference to the A-B difference, and the map it applies through."""

    orderShifts: dict
    mapImages: object


def combine_nod_differences(abDifference, baDifference, orderShifts, mapImages):
    """*average the A-B difference with the B-A difference shifted along the slit onto the A-B positive trace*

    The variance of the mean of two independent frames is ``(var_AB + var_shifted) / 4``. A pixel is flagged when it
    is flagged in either frame. The header of the A-B difference is kept.

    **Key Arguments:**

    - ``abDifference`` -- the A-B difference, with the A header and a standard-deviation uncertainty
    - ``baDifference`` -- the B-A difference, with a standard-deviation uncertainty
    - ``orderShifts`` -- the slit shift in arcsec of each order, from the B-A frame to the A-B frame: the slit
      position of the object at nod B minus its position at nod A
    - ``mapImages`` -- the ``two_d_map_images`` of the frames' binning

    **Return:**

    - ``combined`` -- the combined frame, a new float32 CCDData frame. The inputs are not changed
    """
    from astropy.nddata import CCDData, StdDevUncertainty

    from soxspipe.commonutils.slit_alignment import shift_frame_along_slit

    if abDifference.uncertainty is None or baDifference.uncertainty is None:
        raise ValueError("both differences need an uncertainty to be combined")
    shifted = shift_frame_along_slit(baDifference, orderShifts, mapImages)

    data = (np.asarray(abDifference.data, dtype=np.float32) + shifted.data) / 2
    variance = (np.asarray(abDifference.uncertainty.array, dtype=np.float32) ** 2 + shifted.uncertainty.array**2) / 4
    mask = shifted.mask.copy()
    if abDifference.mask is not None:
        mask |= np.asarray(abDifference.mask) != 0
    return CCDData(
        data,
        unit=abDifference.unit,
        meta=copy.deepcopy(abDifference.meta),
        mask=mask,
        uncertainty=StdDevUncertainty(np.sqrt(variance), unit=abDifference.uncertainty.unit),
    )


def reference_slit_positions(shifts):
    """*the slit position of the object in the reference frame of a group, for each order*

    **Key Arguments:**

    - ``shifts`` -- the shift table of a ``slit_alignment_result``

    **Return:**

    - ``positions`` -- a dictionary from order (``WHOLE_FRAME`` for a position that belongs to the whole frame) to the
      position in arcsec. An order with no reliable trace is left out
    """
    reference = shifts.loc[shifts["is_reference"] & shifts["position_arcsec"].notna()]
    return {
        int(order): float(position)
        for order, position in zip(reference["order"], reference["position_arcsec"], strict=True)
    }


def build_shift_and_add_plan(log, resultA, resultB):
    """*the slit shifts that put the B-A difference onto the A-B difference, from the measured alignment positions*

    The stack of each group is aligned to the group's reference frame, so the object sits at that frame's position
    in the A stack and in the B stack. The shift of an order is the B position minus the A position. An order that
    lacks a position in either group takes the median of the measured shifts, and the log names it.

    **Key Arguments:**

    - ``log`` -- logger
    - ``resultA`` -- the ``slit_alignment_result`` of the A frames, which were traced on A-B
    - ``resultB`` -- the ``slit_alignment_result`` of the B frames, which were traced on B-A

    **Return:**

    - ``plan`` -- a ``shift_and_add_plan``, or *None* (after a warning) when no order has a position in both groups
    """
    mapImages = resultA.mapImages if resultA.mapImages is not None else resultB.mapImages
    if mapImages is None:
        log.warning(NO_PLAN_MESSAGE)
        return None
    positionsA = reference_slit_positions(resultA.shifts)
    positionsB = reference_slit_positions(resultB.shifts)
    orderShifts = {}
    for order in mapImages.orders:
        positionA = positionsA.get(order, positionsA.get(WHOLE_FRAME, np.nan))
        positionB = positionsB.get(order, positionsB.get(WHOLE_FRAME, np.nan))
        orderShifts[order] = float(positionB - positionA)
    measured = [shift for shift in orderShifts.values() if np.isfinite(shift)]
    if not measured:
        log.warning(NO_PLAN_MESSAGE)
        return None
    medianShift = float(np.median(measured))
    unmeasured = [order for order, shift in orderShifts.items() if not np.isfinite(shift)]
    if unmeasured:
        log.warning(
            f"No object position at both nod positions in orders {', '.join(str(order) for order in unmeasured)}, "
            f"so their shift between the nod positions is the median of the other orders ({medianShift:+.3f} arcsec)."
        )
    return shift_and_add_plan(
        {order: shift if np.isfinite(shift) else medianShift for order, shift in orderShifts.items()}, mapImages
    )
