"""Constants shared by the multi-pinhole set-offset tests (DY-1526)."""

from __future__ import annotations

import numpy as np

PINHOLES = 9
# ELEVEN SET OFFSETS SPREAD EVENLY OVER +-0.3 PX: A MAD-SIGMA OF ABOUT 0.27 PX
NOISE_OFFSETS = np.linspace(-0.3, 0.3, 11)
