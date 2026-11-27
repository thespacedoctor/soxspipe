"""Factory for wavelengths that sit on the log-wavelength grid a stack reads (DY-2349)."""

from __future__ import annotations

from soxspipe.commonutils.log_wavelength_grid import grid_index, grid_wavelength

# THE STEP OF THE UNBINNED VIS GRID (KM/S), WHICH THE FACTORY HEADERS AND DETECTOR PARAMETERS OF THE STACK TESTS GIVE
VIS_GRID_STEP_KMS = 17.83


def grid_point_nm(wavelength: float, stepKms: float = VIS_GRID_STEP_KMS) -> float:
    """Return the wavelength (nm) of the grid point of the bin that holds ``wavelength``."""
    return float(grid_wavelength(grid_index(wavelength, stepKms), stepKms))
