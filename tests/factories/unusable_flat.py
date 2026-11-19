"""Factory for a master flat that holds every kind of unusable pixel (DY-2279)."""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

USABLE_FLAT_VALUE = 2.0
# COLUMN 0 IS USABLE; COLUMNS 1 TO 5 ARE A ZERO, A NEGATIVE, A NAN, AN INF AND A -INF
UNUSABLE_FLAT_COLUMNS = (1, 2, 3, 4, 5)
FLAT_SHAPE = (2, 6)


def unusable_flat_pixel_mask() -> np.ndarray:
    """Return the boolean image of the pixels of ``flat_with_unusable_pixels`` that are unusable."""
    unusable = np.zeros(FLAT_SHAPE, dtype=bool)
    unusable[:, list(UNUSABLE_FLAT_COLUMNS)] = True
    return unusable


def flat_with_unusable_pixels() -> CCDData:
    """Return a 2 x 6 master flat whose unusable pixels lie both on an order (row 0) and off it (row 1).

    Row 1 carries the flat's own mask, as the inter-order pixels of a real master flat do.
    """
    row = [USABLE_FLAT_VALUE, 0.0, -0.5, np.nan, np.inf, -np.inf]
    return CCDData(
        np.array([row, row], dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full(FLAT_SHAPE, 0.1, dtype=np.float32)),
        mask=np.array([[False] * 6, [True] * 6]),
    )
