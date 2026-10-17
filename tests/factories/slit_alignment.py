"""Factory for synthetic nodding frames whose objects sit at known slit positions, for slit-alignment tests."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from .frames import instrument_header
from .slit_illumination import write_map_images

ROWS = 600
BAND_COLUMNS = 70
ORDER_PITCH = 80
# THE SLIT COORDINATE CHANGES BY THIS MANY ARCSEC PER COLUMN, AND ITS ISO-LINES TILT BY THIS MANY COLUMNS PER ROW
PIXEL_ARCSEC = 0.2
SLIT_TILT = 0.0137
# THE WAVELENGTH CHANGES BY THIS MANY NM PER ROW, AND ITS ISO-LINES TILT BY THIS MANY ROWS PER COLUMN
NM_PER_ROW = 0.02
WAVELENGTH_TILT = 0.05
FIRST_WAVELENGTH_NM = 500.0
ORDER_WAVELENGTH_STEP_NM = 100.0
OBJECT_SIGMA_ARCSEC = 0.3
OBJECT_PEAK = 400.0
SKY_LEVEL = 50.0
SKY_LINE_PEAK = 300.0
SKY_LINE_PERIOD_ROWS = 40
SKY_LINE_SIGMA_ROWS = 1.5
READ_NOISE = 3.0
# THE ESO KEYWORD THE NODDING SPLIT READS, POSITIVE FOR AN A FRAME
CUMULATIVE_OFFSET = "ESO SEQ CUMOFF Y"
# THE HALF-LENGTH OF THE LIT SLIT IN THE SYNTHETIC MAP, IN ARCSEC
FOOTPRINT_HALF = 5.4

VIS_ORDERS = (1, 2, 3, 4)
NIR_ORDERS = tuple(range(11, 17))


@dataclass(frozen=True)
class alignment_map:
    """A synthetic 2D dispersion map: orders laid side by side as tilted column bands."""

    path: Path
    wavelength: np.ndarray
    slit: np.ndarray
    orderMap: np.ndarray
    orders: tuple[int, ...]
    isTransposed: bool = False

    @property
    def shape(self) -> tuple[int, int]:
        """Return the image shape."""
        return self.wavelength.shape


def _geometry(orders: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the wavelength, slit and order images of orders laid side by side as tilted column bands.

    The slit coordinate is ``PIXEL_ARCSEC * (column - centre - SLIT_TILT * (row - middle))`` and the wavelength is
    ``first + NM_PER_ROW * (row - middle + WAVELENGTH_TILT * (column - centre))``, so neither set of iso-lines is
    aligned with the pixel grid. Pixels outside the lit slit of an order are NaN.
    """
    shape = (ROWS, ORDER_PITCH * len(orders))
    wavelength, slit, orderMap = (np.full(shape, np.nan) for _ in range(3))
    rows = np.arange(ROWS)[:, None] - 0.5 * (ROWS - 1)
    columns = np.arange(BAND_COLUMNS)[None, :] - 0.5 * (BAND_COLUMNS - 1)
    for index, order in enumerate(orders):
        local = PIXEL_ARCSEC * (columns - SLIT_TILT * rows)
        band = FIRST_WAVELENGTH_NM + ORDER_WAVELENGTH_STEP_NM * index + NM_PER_ROW * (rows + WAVELENGTH_TILT * columns)
        inFootprint = np.abs(local) <= FOOTPRINT_HALF
        block = (slice(None), slice(index * ORDER_PITCH, index * ORDER_PITCH + BAND_COLUMNS))
        wavelength[block] = np.where(inFootprint, band, np.nan)
        slit[block] = np.where(inFootprint, local, np.nan)
        orderMap[block] = np.where(inFootprint, float(order), np.nan)
    return wavelength, slit, orderMap


def synthetic_alignment_map(
    workspacePath: Path,
    *,
    orders: tuple[int, ...] = VIS_ORDERS,
    isTransposed: bool = False,
) -> alignment_map:
    """Write a tilted 2D map for ``orders`` and return it with its images.

    The dispersion runs along the rows, as on the SOXS VIS arm. ``isTransposed`` swaps the axes, as on the NIR arm.
    """
    images = _geometry(orders)
    if isTransposed:
        images = tuple(image.T.copy() for image in images)
    path = write_map_images(workspacePath / "alignment" / "twoDMap.fits", *images)
    return alignment_map(path, *images, tuple(orders), isTransposed)


def zero_slit_samples(scene: alignment_map, *, step: int = 20) -> pd.DataFrame:
    """Return what ``detect_continuum.create_pixel_arrays`` returns: the pixel of slit coordinate 0 along each order."""
    slit, orderMap, wavelength = (
        (image.T if scene.isTransposed else image) for image in (scene.slit, scene.orderMap, scene.wavelength)
    )
    records = []
    for order in scene.orders:
        for row in range(30, ROWS - 30, step):
            columns = np.flatnonzero(orderMap[row] == order)
            column = float(np.interp(0.0, slit[row][columns], columns))
            records.append(
                {
                    "order": float(order),
                    "wavelength": float(wavelength[row, int(round(column))]),
                    "slit_position": 0.0,
                    "fit_x": row if scene.isTransposed else column,
                    "fit_y": column if scene.isTransposed else row,
                }
            )
    return pd.DataFrame(records)


def _sky(wavelength: np.ndarray) -> np.ndarray:
    """Return the sky counts, a function of wavelength alone: a continuum and evenly spaced Gaussian lines."""
    phase = (np.nan_to_num(wavelength) / NM_PER_ROW) % SKY_LINE_PERIOD_ROWS
    distance = np.minimum(phase, SKY_LINE_PERIOD_ROWS - phase)
    return SKY_LEVEL + SKY_LINE_PEAK * np.exp(-0.5 * (distance / SKY_LINE_SIGMA_ROWS) ** 2)


def synthetic_nod_frame(
    scene: alignment_map,
    *,
    objectSlit: float | Mapping[int, float],
    mjd: float,
    offset: float = 3.0,
    arm: str = "VIS",
    objectPeak: float = OBJECT_PEAK,
    hasSky: bool = True,
    badPixels: tuple[tuple[int, int], ...] = (),
    seed: int | None = None,
    headerOverrides: Mapping[str, object] | None = None,
) -> CCDData:
    """Return a prepared frame: sky plus a Gaussian object trace at ``objectSlit`` arcsec in each order.

    ``objectSlit`` is one slit position for every order, or a mapping from order to slit position (an order left
    out of the mapping carries no object). The sky depends on wavelength alone, so it cancels in a difference
    of two frames, and ``hasSky=False`` leaves it out. ``badPixels`` are ``(row, column)`` pixels flagged in the
    mask. ``seed`` adds white noise.
    """
    positions = (
        {order: float(objectSlit) for order in scene.orders}
        if not isinstance(objectSlit, Mapping)
        else {order: float(slit) for order, slit in objectSlit.items()}
    )
    data = _sky(scene.wavelength) if hasSky else np.zeros(scene.shape)
    for order, slitPosition in positions.items():
        isOrder = scene.orderMap == order
        profile = np.exp(-0.5 * ((np.where(isOrder, scene.slit, np.nan) - slitPosition) / OBJECT_SIGMA_ARCSEC) ** 2)
        data = data + objectPeak * np.nan_to_num(profile)
    data = np.where(np.isfinite(scene.orderMap), data, 0.0)
    sigma = np.sqrt(READ_NOISE**2 + np.clip(data, 0.0, None))
    if seed is not None:
        data = data + np.random.default_rng(seed).normal(0.0, 1.0, data.shape) * sigma
    mask = np.zeros(data.shape, dtype=bool)
    for row, column in badPixels:
        mask[row, column] = True
    header = instrument_header(
        arm=arm,
        overrides={"EXPTIME": 60.0, "MJD-OBS": mjd, CUMULATIVE_OFFSET: offset, **(headerOverrides or {})},
    )
    return CCDData(
        data,
        unit=u.electron,
        meta=header,
        mask=mask,
        uncertainty=StdDevUncertainty(sigma, unit=u.electron),
    )


__all__ = [
    "NIR_ORDERS",
    "PIXEL_ARCSEC",
    "VIS_ORDERS",
    "alignment_map",
    "synthetic_alignment_map",
    "synthetic_nod_frame",
]
