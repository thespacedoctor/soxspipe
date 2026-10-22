"""Robust set-level offset test for multi-pinhole arc-line sets.

In a multi-pinhole frame each arc line is measured once per pinhole, and those measurements form a set. A misidentified
or blended line sits off the model by the same amount in every pinhole, with little scatter inside the set, so it
survives clips that judge each pinhole on its own. This module finds such sets by comparing each set's median signed
residual with the robust spread of set offsets in the same order.
"""

import numpy as np

# DEFAULT CUT, IN MAD-SIGMA OF THE SET OFFSETS IN AN ORDER (DY-1526 RULING)
OFFSET_SET_CLIP_SIGMA = 5.0
# FEWER SETS THAN THIS GIVE NO MEANINGFUL MAD, SO THE ORDER IS SKIPPED
MIN_SETS_PER_ORDER = 5
# A LINE MEASURED IN ONE PINHOLE ONLY IS NOT A SET
MIN_PINHOLES_PER_SET = 2
MAD_TO_SIGMA = 1.4826
SET_KEYS = ["wavelength", "order"]
RESIDUAL_COLUMNS = ["residuals_x", "residuals_y"]


def flag_offset_pinhole_sets(lines, sigmaLimit=OFFSET_SET_CLIP_SIGMA, minSets=MIN_SETS_PER_ORDER):
    """*flag the multi-pinhole sets whose signed offset from the model is an outlier within their order*

    For each set (one arc line, all its pinholes) the median signed residual is taken separately in x and in y. Within
    each order, a set is flagged when its median offset lies more than ``sigmaLimit`` MAD-sigma from the median of
    that order's set offsets, in either axis. Orders with fewer than ``minSets`` sets, or with no spread, are skipped.
    Lines measured in only one pinhole are never flagged, so single-pinhole fits are unchanged.

    **Key Arguments:**

    - ``lines`` -- dataframe with ``wavelength``, ``order``, ``residuals_x`` and ``residuals_y`` columns, one row per
      measured pinhole
    - ``sigmaLimit`` -- the cut, in MAD-sigma of the set offsets in the order. Default *OFFSET_SET_CLIP_SIGMA*
    - ``minSets`` -- the smallest number of sets an order needs to be tested. Default *MIN_SETS_PER_ORDER*

    **Return:**

    - ``flagged`` -- boolean array, parallel to the rows of ``lines``, *True* for every pinhole of a flagged set

    **Usage:**

    ```python
    from soxspipe.commonutils.pinhole_set_offsets import flag_offset_pinhole_sets
    flagged = flag_offset_pinhole_sets(orderPixelTable, sigmaLimit=5)
    ```
    """
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
            sigma = MAD_TO_SIGMA * np.median(deviations)
            # NO SPREAD (OR NaN) GIVES NO SCALE TO JUDGE AN OUTLIER AGAINST
            if not sigma > 0:
                continue
            isOutlier |= deviations > sigmaLimit * sigma
        flaggedSets.extend(orderSets.index[isOutlier])

    return lines.set_index(SET_KEYS).index.isin(flaggedSets)
