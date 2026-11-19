"""The rule for the master-flat pixels that cannot flat-field or unflatten a frame (DY-2279).

A master-flat pixel is unusable when its flux is zero, negative or not finite. Every step that divides by the
master flat, or multiplies a frame back by it, uses this rule, so the steps cannot disagree about which pixels to
drop. Dividing by such a pixel gives ``inf``, a huge value or a sign flip, and a ``RuntimeWarning`` as well.
"""

from __future__ import annotations


def unusable_flat_pixel_mask(flatData):
    """*return the boolean image of the master-flat pixels that cannot be used*

    **Key Arguments:**

    - ``flatData`` -- the master-flat flux array

    **Return:**

    - ``unusable`` -- boolean array, True where the flux is zero, negative, NaN or infinite

    **Usage:**

    ```python
    unusable = unusable_flat_pixel_mask(master_flat.data)
    ```
    """
    import numpy as np

    return ~(np.isfinite(flatData) & (flatData > 0))


def neutralise_unusable_flat(flat, unusable):
    """*return a copy of the master flat that is safe to divide by or multiply with*

    The unusable pixels hold 1 in the data and 0 in the uncertainty, so arithmetic with the copy raises no
    ``RuntimeWarning``. The results at those pixels mean nothing: pass the result to ``flag_unusable_pixels``.

    **Key Arguments:**

    - ``flat`` -- the master flat (CCDData). It is not changed
    - ``unusable`` -- boolean image of the unusable pixels, from ``unusable_flat_pixel_mask``

    **Return:**

    - ``neutralFlat`` -- a copy of ``flat`` with the unusable pixels neutralised

    **Usage:**

    ```python
    neutralFlat = neutralise_unusable_flat(master_flat, unusable)
    ```
    """
    neutralFlat = flat.copy()
    neutralFlat.data[unusable] = 1.0
    if neutralFlat.uncertainty is not None:
        neutralFlat.uncertainty.array[unusable] = 0.0
    return neutralFlat


def flag_unusable_pixels(frame, unusable):
    """*return a copy of a frame with the unusable pixels set to NaN and masked*

    The data and the uncertainty become NaN, and the mask is set, where ``unusable`` is True. Integer data or an
    integer uncertainty array is promoted to float first, because an integer array cannot hold NaN. ``frame`` is
    not changed.

    **Key Arguments:**

    - ``frame`` -- the frame (CCDData) to flag
    - ``unusable`` -- boolean image of the unusable pixels, from ``unusable_flat_pixel_mask``

    **Return:**

    - ``flaggedFrame`` -- a copy of ``frame``, flagged

    **Usage:**

    ```python
    flaggedFrame = flag_unusable_pixels(frame, unusable)
    ```
    """
    import numpy as np

    flaggedFrame = frame.copy()
    flaggedFrame.data = flaggedFrame.data.astype(np.result_type(flaggedFrame.data.dtype, np.float32), copy=False)
    flaggedFrame.data[unusable] = np.nan
    if flaggedFrame.uncertainty is not None:
        uncertaintyArray = flaggedFrame.uncertainty.array
        uncertaintyArray = uncertaintyArray.astype(np.result_type(uncertaintyArray.dtype, np.float32), copy=False)
        uncertaintyArray[unusable] = np.nan
        flaggedFrame.uncertainty.array = uncertaintyArray
    flaggedFrame.mask = unusable.copy() if flaggedFrame.mask is None else (flaggedFrame.mask | unusable)
    return flaggedFrame
