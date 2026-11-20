"""Read the number of integrations (NDIT) a raw frame averages (DY-2317).

A raw NIR file is the mean of its ``ESO DET NDIT`` integrations, so the variance of its pixels is the variance of one
integration divided by NDIT. VIS and UVB frames carry no NDIT card and are one integration.
"""

import math
from numbers import Real


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

    - ``ndit`` -- the NDIT as an integer of at least 1. A missing card gives 1 with no message; a card that is not a
      whole number of at least 1 gives 1 and a warning

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
    if not _is_whole_number(value) or value < 1:
        log.warning(f"The {card} card holds {value!r}, which is not a whole number of at least 1; NDIT is taken as 1.")
        return 1

    return int(value)
