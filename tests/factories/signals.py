"""Factories for deterministic synthetic signals."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def synthetic_signal(
    *,
    size: int = 64,
    seed: int = 7,
    noise: float = 0.01,
) -> NDArray[np.float64]:
    """Return a seeded Gaussian signal on a linear continuum."""
    pixels = np.arange(size, dtype=float)
    centre = (size - 1) / 2
    profile = np.exp(-0.5 * ((pixels - centre) / 3.0) ** 2)
    rng = np.random.default_rng(seed)
    signal = 1.0 + (0.001 * pixels) + profile + rng.normal(0.0, noise, size)
    signal.setflags(write=False)
    return signal
