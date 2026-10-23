"""Factory for synthetic echelle frames with skylines, a flat and a 2D map, for slit-illumination tests."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from scipy.special import ndtr

SCENE_ROWS = 600
SCENE_COLUMNS = 70
SCENE_ORDER_PITCH = 80
# THE SLIT CENTRE DRIFTS ACROSS THE COLUMNS BY THIS MANY PIXELS PER ROW, AS A TILTED SLIT DOES
SCENE_TILT_PIXELS_PER_ROW = 0.0137
SCENE_SLIT_PIXEL_ARCSEC = 0.2
SCENE_NM_PER_PIXEL = 0.02
SCENE_LINE_SIGMA_PIXELS = 1.3
SCENE_LINE_SPACING_NM = 0.8
SCENE_LINES_PER_ORDER = 14
SCENE_ORDER_WAVELENGTH_STEP_NM = 100.0
SCENE_FIRST_WAVELENGTH_NM = 500.0
SCENE_READ_NOISE = 3.3


@dataclass(frozen=True)
class illumination_scene:
    """A synthetic detrended frame, master flat and 2D map with known slit offsets."""

    frame: CCDData
    masterFlat: CCDData
    twoDMapPath: Path
    skylinesDF: pd.DataFrame
    wavelength: np.ndarray
    slit: np.ndarray
    orderMap: np.ndarray
    orders: tuple[int, ...]


def top_hat(slit: np.ndarray, amplitude: float, lower: float, upper: float, width: float) -> np.ndarray:
    """Return an erf top-hat: ``amplitude`` between ``lower`` and ``upper`` with Gaussian-blurred edges."""
    return amplitude * (ndtr((slit - lower) / width) - ndtr((slit - upper) / width))


def _order_geometry(orders: tuple[int, ...], footprintHalf: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the wavelength, slit and order images of orders laid side by side as column bands."""
    shape = (SCENE_ROWS, SCENE_ORDER_PITCH * len(orders))
    wavelength = np.full(shape, np.nan)
    slit = np.full(shape, np.nan)
    orderMap = np.full(shape, np.nan)
    rows = np.arange(SCENE_ROWS)
    drift = SCENE_TILT_PIXELS_PER_ROW * (rows - 0.5 * (SCENE_ROWS - 1))
    for index, order in enumerate(orders):
        columns = np.arange(SCENE_COLUMNS)
        local = (columns[None, :] - 0.5 * (SCENE_COLUMNS - 1) + drift[:, None]) * SCENE_SLIT_PIXEL_ARCSEC
        start = index * SCENE_ORDER_PITCH
        inFootprint = np.abs(local) <= footprintHalf
        firstWavelength = SCENE_FIRST_WAVELENGTH_NM + SCENE_ORDER_WAVELENGTH_STEP_NM * index
        band = np.broadcast_to((firstWavelength + rows * SCENE_NM_PER_PIXEL)[:, None], local.shape)
        block = (slice(None), slice(start, start + SCENE_COLUMNS))
        wavelength[block] = np.where(inFootprint, band, np.nan)
        slit[block] = np.where(inFootprint, local, np.nan)
        orderMap[block] = np.where(inFootprint, float(order), np.nan)
    return wavelength, slit, orderMap


def _skylines(orders: tuple[int, ...]) -> pd.DataFrame:
    """Return an evenly spaced skyline table, with a different brightness for every line."""
    wavelengths, fluxes = [], []
    for index in range(len(orders)):
        firstWavelength = SCENE_FIRST_WAVELENGTH_NM + SCENE_ORDER_WAVELENGTH_STEP_NM * index
        for lineIndex in range(SCENE_LINES_PER_ORDER):
            wavelengths.append(firstWavelength + SCENE_LINE_SPACING_NM * (lineIndex + 0.5))
            fluxes.append(1.0 + 0.5 * np.sin(lineIndex + index))
    return pd.DataFrame({"WAVELENGTH": wavelengths, "FLUX": fluxes})


def _sky_spectrum(wavelength: np.ndarray, skylinesDF: pd.DataFrame, lineAmplitude: float, continuum: float):
    """Return the sky counts per pixel: Gaussian lines on a flat continuum."""
    sigmaNm = SCENE_LINE_SIGMA_PIXELS * SCENE_NM_PER_PIXEL
    sky = np.full(wavelength.shape, continuum)
    for lineWavelength, flux in zip(skylinesDF["WAVELENGTH"], skylinesDF["FLUX"], strict=True):
        sky += lineAmplitude * flux * np.exp(-0.5 * ((wavelength - lineWavelength) / sigmaNm) ** 2)
    return sky


def write_map_images(path: Path, wavelength: np.ndarray, slit: np.ndarray, orderMap: np.ndarray) -> Path:
    """Write the three map images in the extension layout of a 2D dispersion map."""
    hdus = fits.HDUList(
        [
            fits.PrimaryHDU(wavelength.astype("float32"), header=fits.Header({"EXTNAME": "WAVELENGTH"})),
            fits.ImageHDU(slit.astype("float32"), name="SLIT"),
            fits.ImageHDU(orderMap.astype("float32"), name="ORDER"),
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    hdus.writeto(path)
    return path


def synthetic_illumination_scene(
    workspacePath: Path,
    *,
    orders: tuple[int, ...] = (11, 12),
    offsets: Mapping[int, float] | None = None,
    flatEdges: tuple[float, float, float] = (-5.0, 5.0, 0.15),
    scienceWidth: float | None = None,
    footprintHalf: float = 5.4,
    objectSlit: float | None = None,
    objectCounts: float = 0.0,
    continuum: float = 20.0,
    lineAmplitude: float = 600.0,
    seed: int = 3,
    hasNoise: bool = True,
    slitName: str | None = "SLIT1.0",
) -> illumination_scene:
    """Return a detrended frame whose lit slit is shifted by ``offsets[order]`` arcsec from the master flat.

    The raw counts are (sky lines + continuum + optional object) x the shifted top-hat x a blaze. The master flat
    is the unshifted top-hat x the same blaze. The frame is raw / flat, with StdDev uncertainty, and is masked
    where the flat is below half its slit-profile maximum, as the master flat's own quality mask does. The frame
    header names the slit (``slitName``, or no slit keyword if *None*) as the VIS slit keyword does.
    """
    shifts = dict(offsets or {})
    lower, upper, flatWidth = flatEdges
    sciWidth = flatWidth if scienceWidth is None else scienceWidth
    wavelength, slit, orderMap = _order_geometry(orders, footprintHalf)
    skylinesDF = _skylines(orders)
    sky = _sky_spectrum(wavelength, skylinesDF, lineAmplitude, continuum)
    if objectSlit is not None:
        sky = sky + objectCounts * np.exp(-0.5 * ((slit - objectSlit) / 0.5) ** 2)
    blaze = 1.0 + 0.2 * np.sin(wavelength / 3.0)
    illumination = np.zeros(slit.shape)
    for order in orders:
        delta = float(shifts.get(order, 0.0))
        isOrder = orderMap == order
        illumination[isOrder] = top_hat(slit[isOrder], 1.0, lower + delta, upper + delta, sciWidth)
    flatProfile = top_hat(np.nan_to_num(slit), 1.0, lower, upper, flatWidth)
    flat = np.where(np.isfinite(slit), flatProfile * blaze, np.nan)
    raw = np.nan_to_num(sky * illumination * blaze)
    sigma = np.sqrt(SCENE_READ_NOISE**2 + np.clip(raw, 0, None))
    if hasNoise:
        raw = raw + np.random.default_rng(seed).normal(0.0, 1.0, raw.shape) * sigma
    isMasked = ~np.isfinite(slit) | (flatProfile < 0.5)
    safeFlat = np.where(isMasked, 1.0, flat)
    data = np.where(isMasked, 0.0, raw / safeFlat)
    errs = np.where(isMasked, 1.0, sigma / safeFlat)
    frame = CCDData(
        data,
        unit=u.electron,
        mask=isMasked,
        uncertainty=StdDevUncertainty(errs, unit=u.electron),
        meta={} if slitName is None else {"ESO INS VISE NAME": slitName},
    )
    masterFlat = CCDData(np.nan_to_num(flat), unit=u.dimensionless_unscaled, mask=~np.isfinite(slit))
    mapPath = write_map_images(workspacePath / "scene" / "twoDMap.fits", wavelength, slit, orderMap)
    return illumination_scene(frame, masterFlat, mapPath, skylinesDF, wavelength, slit, orderMap, tuple(orders))
