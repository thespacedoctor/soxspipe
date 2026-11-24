"""Robust set-level offset test for multi-pinhole arc-line sets.

In a multi-pinhole frame each arc line is measured once per pinhole, and those measurements form a set. A misidentified
or blended line sits off the model by the same amount in every pinhole, with little scatter inside the set, so it
survives clips that judge each pinhole on its own. This module finds such sets by comparing each set's median signed
residual with the robust spread of set offsets in the same order.

The spread has a floor (``minSigma``, in pixels). Without it, the decision about one set depends on how tightly the
other correct sets in its order scatter: when they improve, the MAD shrinks and sets that were kept before are dropped
(DY-2170). With the floor, a set in an order that scatters less than the floor is judged against a fixed cut of
``sigmaLimit * minSigma`` pixels.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # TYPING-ONLY IMPORTS. NUMPY IS IMPORTED INSIDE THE FUNCTION THAT USES IT, AS IN create_dispersion_map
    import numpy as np
    import numpy.typing as npt
    import pandas as pd

# THE FITTED_LINES COLUMN THAT MARKS THE LINES OF A SET DROPPED BY THIS TEST
SET_OFFSET_COLUMN = "set_offset_clipped"
# DEFAULT CUT, IN MAD-SIGMA OF THE SET OFFSETS IN AN ORDER (DY-1526 RULING)
OFFSET_SET_CLIP_SIGMA = 5.0
# DEFAULT FLOOR ON THE MAD-SIGMA OF THE SET OFFSETS, IN PIXELS. 0 KEEPS THE MAD-ONLY RULE OF DY-1526; THE SETTINGS FILES
# SET THE NIR VALUE (DY-2170)
OFFSET_SET_MIN_SIGMA = 0.0
# FEWER SETS THAN THIS GIVE NO MEANINGFUL MAD, SO THE ORDER IS SKIPPED
MIN_SETS_PER_ORDER = 5
# A LINE MEASURED IN ONE PINHOLE ONLY IS NOT A SET
MIN_PINHOLES_PER_SET = 2
MAD_TO_SIGMA = 1.4826
SET_KEYS = ["wavelength", "order"]
RESIDUAL_COLUMNS = ["residuals_x", "residuals_y"]


def flag_offset_pinhole_sets(
    lines: pd.DataFrame,
    sigmaLimit: float = OFFSET_SET_CLIP_SIGMA,
    minSets: int = MIN_SETS_PER_ORDER,
    minSigma: float = OFFSET_SET_MIN_SIGMA,
) -> npt.NDArray[np.bool_]:
    """*flag the multi-pinhole sets whose signed offset from the model is an outlier within their order*

    For each set (one arc line, all its pinholes) the median signed residual is taken separately in x and in y. Within
    each order, a set is flagged when its median offset lies more than ``sigmaLimit`` MAD-sigma from the median of
    that order's set offsets, in either axis. The MAD-sigma is never taken below ``minSigma`` pixels, so in an order
    whose sets scatter less than that the cut is ``sigmaLimit * minSigma`` pixels, whatever the other sets do. Orders
    with fewer than ``minSets`` sets, or with no spread and no floor, are skipped.
    Lines measured in only one pinhole are never flagged, so single-pinhole fits are unchanged.

    **Key Arguments:**

    - ``lines`` -- dataframe with ``wavelength``, ``order``, ``residuals_x`` and ``residuals_y`` columns, one row per
      measured pinhole
    - ``sigmaLimit`` -- the cut, in MAD-sigma of the set offsets in the order. Default *OFFSET_SET_CLIP_SIGMA*
    - ``minSets`` -- the smallest number of sets an order needs to be tested. Default *MIN_SETS_PER_ORDER*
    - ``minSigma`` -- the floor on the MAD-sigma, in pixels. Default *OFFSET_SET_MIN_SIGMA* (no floor)

    **Return:**

    - ``flagged`` -- boolean array, parallel to the rows of ``lines``, *True* for every pinhole of a flagged set

    **Usage:**

    ```python
    from soxspipe.commonutils.pinhole_set_offsets import flag_offset_pinhole_sets
    flagged = flag_offset_pinhole_sets(orderPixelTable, sigmaLimit=5, minSigma=0.08)
    ```
    """
    import numpy as np

    if not (np.isfinite(minSigma) and minSigma >= 0):
        raise ValueError(f"minSigma must be a finite number, zero or positive, not {minSigma}")

    if lines.empty:
        return np.zeros(0, dtype=bool)

    grouped = lines.groupby(SET_KEYS, sort=True)
    setOffsets = grouped[RESIDUAL_COLUMNS].median()
    setOffsets = setOffsets.loc[grouped.size() >= MIN_PINHOLES_PER_SET]

    flaggedSets = []
    for _, orderSets in setOffsets.groupby(level="order", sort=True):
        if len(orderSets) < minSets:
            continue
        isOutlier = np.zeros(len(orderSets), dtype=bool)
        for column in RESIDUAL_COLUMNS:
            offsets = orderSets[column].to_numpy(dtype=float)
            deviations = np.abs(offsets - np.median(offsets))
            sigma = max(MAD_TO_SIGMA * np.median(deviations), minSigma)
            # NO SPREAD AND NO FLOOR GIVES NO SCALE TO JUDGE AN OUTLIER AGAINST. A NaN OFFSET MAKES THE MEDIAN NaN, AND
            # max() THEN KEEPS THE NaN, SO THE AXIS IS SKIPPED AS BEFORE THE FLOOR
            if not sigma > 0:
                continue
            isOutlier |= deviations > sigmaLimit * sigma
        flaggedSets.extend(orderSets.index[isOutlier])

    return lines.set_index(SET_KEYS).index.isin(flaggedSets)
