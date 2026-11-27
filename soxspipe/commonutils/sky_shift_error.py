"""*the pure parts of the uncertainty of an order's sky shift (DY-2359)*

The ``SKY SHIFT O<n> ERR`` QC is the uncertainty (nm) of the wavelength shift that the skyline tuning applied to order
n. Its formula depends on where the shift came from (DY-2307):

- a shift measured from the sky lines of the order has the standard error of the median of the kept line offsets,
  ``scatter * sqrt(pi / 2) / sqrt(n)``;
- a shift predicted by a model has the model's rms in px, times the dispersion of the order. The model's rms for the
  VIS r order is in ``r-from-i-transfer``, and for the u and g orders in ``ug-wavelength-transfer``;
- an order that was not shifted has the fixed per-arm uncertainty of an untuned wavelength solution, in the
  ``sky-shift-uncertainty-nm`` block of the detector parameters.

The pieces that need the classification of a shift (``measured``, ``model`` or ``none``) live in the extraction, which
holds it. This module imports nothing from a recipe.

Author
: David Young

Date Created
: November 27, 2026
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

# THE STANDARD ERROR OF THE MEDIAN OF A NORMAL SAMPLE IS THIS FACTOR TIMES THE STANDARD ERROR OF THE MEAN (1.2533)
MEDIAN_ERROR_FACTOR = math.sqrt(math.pi / 2)

UNCERTAINTY_KEY = "sky-shift-uncertainty-nm"
# THE TWO WAYS A SOLUTION CAN BE LEFT UNTUNED: TUNING SWITCHED OFF OR AN ORDER THAT FAILED (untuned), AND TUNING
# SKIPPED BECAUSE AN ARC IS ATTACHED TO THE SCIENCE FRAME (attached-arc, DY-1665)
UNTUNED_KEY = "untuned"
ATTACHED_ARC_KEY = "attached-arc"


def median_standard_error_nm(scatterNm: float, nLines: int) -> float:
    """*the standard error of the median of the kept sky-line offsets of one order*

    **Key Arguments:**

    - ``scatterNm`` -- the standard deviation (nm) of the kept sky-line offsets
    - ``nLines`` -- the number of kept sky lines, at least 1

    **Return:**

    - ``errorNm`` -- ``scatterNm * sqrt(pi / 2) / sqrt(nLines)`` (nm)

    **Usage:**

    ```python
    from soxspipe.commonutils.sky_shift_error import median_standard_error_nm
    errorNm = median_standard_error_nm(stats["scatterNm"], stats["nLines"])
    ```
    """
    return MEDIAN_ERROR_FACTOR * scatterNm / math.sqrt(nLines)


def uncertainty_of_untuned_order_nm(detectorParams: Mapping[str, Any], kind: str) -> float | None:
    """*the fixed uncertainty (nm) of an order that the skyline tuning did not shift*

    **Key Arguments:**

    - ``detectorParams`` -- the arm's detector parameters
    - ``kind`` -- ``UNTUNED_KEY`` or ``ATTACHED_ARC_KEY``, the key of ``sky-shift-uncertainty-nm``

    **Return:**

    - ``errorNm`` -- the per-arm value (nm), or *None* if the arm has no value under ``kind``. An arm leaves a key out
      when its value was not measured

    **Raises:**

    - ``ValueError`` -- if the arm has a value that is not a finite positive number

    **Usage:**

    ```python
    from soxspipe.commonutils.sky_shift_error import UNTUNED_KEY, uncertainty_of_untuned_order_nm
    errorNm = uncertainty_of_untuned_order_nm(self.detectorParams, UNTUNED_KEY)
    ```
    """
    block = detectorParams.get(UNCERTAINTY_KEY) or {}
    if kind not in block:
        return None
    try:
        errorNm = float(block[kind])
    except (TypeError, ValueError) as e:
        raise ValueError(f"the {UNCERTAINTY_KEY} detector parameters need a numeric `{kind}` value ({e})") from e
    if not math.isfinite(errorNm) or errorNm <= 0:
        raise ValueError(f"the {UNCERTAINTY_KEY} detector parameters need a finite positive `{kind}` value")
    return errorNm
