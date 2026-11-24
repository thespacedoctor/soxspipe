"""Read the number of integrations (NDIT) a raw frame averages (DY-2317).

A raw NIR file is the mean of its ``ESO DET NDIT`` integrations, so the variance of its pixels is the variance of one
integration divided by NDIT. VIS and UVB frames carry no NDIT card and are one integration.
"""

import math
from numbers import Real

import numpy as np

# THE LARGEST NDIT A CARD MAY HOLD. THE ARCHIVE MAXIMUM IS 600 (AN ABORTED TEST DARK), SO 10000 IS FAR ABOVE ANY REAL
# VALUE AND FAR BELOW THE POINT WHERE int(value) STOPS FITTING A NUMPY INTEGER
MAX_NDIT = 10000


def _is_whole_number(value):
    """*True for a finite integer or integer-valued float, and False for anything else, a boolean included*"""
    if isinstance(value, bool) or not isinstance(value, Real):
        return False
    return math.isfinite(value) and value == int(value)


def read_ndit(header, kw, log):
    """*read the number of integrations a raw frame averages*

    **Key Arguments:**

    - ``header`` -- the frame header, or any mapping keyed by the FITS card name
    - ``kw`` -- the keyword lookup (``keyword_lookup(...).get``)
    - ``log`` -- logger

    **Return:**

    - ``ndit`` -- the NDIT as an integer from 1 to ``MAX_NDIT`` (10000, about 16 times the largest archive value of 600,
      and far below the point where a huge card value would no longer fit a NumPy integer). A missing card gives 1
      with no message; a card that is not a whole number in that range gives 1 and a warning

    **Usage:**

    ```python
    from soxspipe.commonutils.ndit import read_ndit
    ndit = read_ndit(frame.header, kw, log)
    ```
    """
    card = kw("DET_NDIT")
    if card not in header:
        return 1

    value = header[card]
    if not _is_whole_number(value) or not 1 <= value <= MAX_NDIT:
        log.warning(
            f"The {card} card holds {value!r}, which is not a whole number from 1 to {MAX_NDIT}; NDIT is taken as 1."
        )
        return 1

    return int(value)


def scale_by_sqrt_ndit(values: np.ndarray, ndit: int) -> np.ndarray:
    """*scale a noise weight or a signal-to-noise by the square root of NDIT*

    A frame is the mean of NDIT integrations, so its variance is that of one integration divided by NDIT. Every
    1/sigma or signal-to-noise is then sqrt(NDIT) times the NDIT 1 value (DY-2323).

    **Key Arguments:**

    - ``values`` -- numpy array of weights or signal-to-noise values
    - ``ndit`` -- the number of integrations the frame averages (``read_ndit``)

    **Return:**

    - ``scaled`` -- ``values`` itself for NDIT 1, so the NDIT 1 result is bit-identical to the unscaled one, else
      ``values`` times sqrt(NDIT) as a new array
    """
    return values if ndit == 1 else values * np.sqrt(ndit)
