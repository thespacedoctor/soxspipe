"""*synthetic arc and sky frames, with the maps of a DISP_IMAGE, for the line-measurement tests of tools/*

The wavelength map is linear along the dispersion and tilted by `ROW_TILT` pixels per row, so the pixels of one
wavelength bin come from many columns, as they do on a real frame. The lines are Gaussian in wavelength, so they follow
that map and their FWHM in pixels is known.

The arrays have the layout of the real `DISP_IMAGE` maps (DY-2455): the wavelength runs along numpy axis 0 for VIS and
along numpy axis 1 for NIR. On the 2026-01-12 maps the median step between neighbouring pixels is 0.0335 nm along axis 0
and 0.0051 nm along axis 1 for VIS, and 0.0088 nm along axis 0 and 0.0498 nm along axis 1 for NIR.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import numpy.typing as npt
from astropy.io import fits

ROWS = 60
COLUMNS = 2000
NM_PER_PIXEL = 0.03
START_NM = 500.0
BACKGROUND = 100.0
AMPLITUDE = 5000.0
FWHM_PER_SIGMA = 2.3548
LINE_COLUMNS = list(range(150, 1900, 150))
ROW_TILT = 0.37
DISPERSION_NUMPY_AXIS = {"VIS": 0, "NIR": 1}


class Scene:
    """*a synthetic frame with the maps that a DISP_IMAGE carries, in the real layout of the arm*

    **Key Arguments:**

    - ``arm`` -- ``VIS`` (wavelength along numpy axis 0) or ``NIR`` (along axis 1). Default ``VIS``
    - ``fwhmPx`` -- the FWHM of every line, in pixels. Default 3
    - ``lineColumns`` -- the columns of the line centres at the middle row. Default every 150 columns
    - ``minSlit`` -- the lines are present only where ``|slit|`` is at least this, in arcsec. Default 0
    - ``maxSlit`` -- the lines are present only where ``|slit|`` is at most this, in arcsec. Default 100
    - ``side`` -- 1 to light the positive side of the slit only, -1 for the negative side only, 0 for both. Default 0
    - ``seed`` -- the seed of the noise. Default 0
    """

    def __init__(
        self,
        arm: str = "VIS",
        fwhmPx: float = 3.0,
        lineColumns: list[int] | None = None,
        minSlit: float = 0.0,
        maxSlit: float = 100.0,
        side: int = 0,
        seed: int = 0,
    ) -> None:
        lineColumns = LINE_COLUMNS if lineColumns is None else lineColumns
        columns = np.arange(COLUMNS)
        rows = np.arange(ROWS)
        self.wavelength = START_NM + NM_PER_PIXEL * (columns[None, :] + ROW_TILT * rows[:, None])
        slitAxis = np.linspace(-6.0, 6.0, ROWS)
        self.slit = np.tile(slitAxis[:, None], (1, COLUMNS))
        self.order = np.ones((ROWS, COLUMNS), dtype=int)
        self.fwhmPx = fwhmPx
        self.centres = [START_NM + NM_PER_PIXEL * (column + ROW_TILT * ROWS / 2) for column in lineColumns]
        rng = np.random.default_rng(seed)
        data = BACKGROUND + rng.normal(0.0, 1.0, (ROWS, COLUMNS))
        sigmaNm = fwhmPx * NM_PER_PIXEL / FWHM_PER_SIGMA
        isLit = (np.abs(self.slit) >= minSlit) & (np.abs(self.slit) <= maxSlit) & (self.slit * side >= 0)
        for centre in self.centres:
            data += AMPLITUDE * np.exp(-0.5 * ((self.wavelength - centre) / sigmaNm) ** 2) * isLit
        self.data = data
        self.good = np.ones((ROWS, COLUMNS), dtype=bool)
        self.arm = arm
        self.isTransposed = DISPERSION_NUMPY_AXIS[arm] == 0
        if self.isTransposed:
            for name in ("wavelength", "slit", "order", "data", "good"):
                setattr(self, name, getattr(self, name).T.copy())

    def slit_row(self, slitArcsec: float) -> int:
        """*the slit-position index (the ``row`` of ``pixel_at``) nearest to a slit position in arcsec*"""
        slitAxis = self.slit[0, :] if self.isTransposed else self.slit[:, 0]
        return int(np.argmin(np.abs(slitAxis - slitArcsec)))

    def pixel_at(self, centreNm: float, row: int) -> tuple[int, int]:
        """*the array index of the pixel, along a slit position ``row``, whose wavelength is closest to a wavelength*"""
        column = int(round((centreNm - START_NM) / NM_PER_PIXEL - ROW_TILT * row))
        return (column, row) if self.isTransposed else (row, column)

    def expected_resolving_power(self, centreNm: float) -> float:
        """*the resolving power of a line at a wavelength: centre over FWHM in nm*"""
        return centreNm / (self.fwhmPx * NM_PER_PIXEL)


def write_raw_arc(path: Path, scene: Scene, header: fits.Header) -> Path:
    """*write a scene as a raw frame with the given header*"""
    fits.PrimaryHDU(scene.data.astype(np.float32), header=header).writeto(path)
    return path


def write_disp_image(path: Path, scene: Scene, mjd: float = 60000.1) -> Path:
    """*write a scene's maps as a DISP_IMAGE of the scene's arm: ``WAVELENGTH`` primary, then ``SLIT`` and ``ORDER``*"""
    header = fits.Header()
    header["ESO SEQ ARM"] = scene.arm
    header["MJD-OBS"] = mjd
    primary = fits.PrimaryHDU(scene.wavelength.astype(np.float32), header=header)
    primary.header["EXTNAME"] = "WAVELENGTH"
    slitHdu = fits.ImageHDU(scene.slit.astype(np.float32), name="SLIT")
    orderHdu = fits.ImageHDU(scene.order.astype(np.int32), name="ORDER")
    fits.HDUList([primary, slitHdu, orderHdu]).writeto(path)
    return path


def write_sky_pair(path: Path, scene: Scene, header: fits.Header, badPixel: tuple[int, int] | None = None) -> Path:
    """*write a scene as a ``<stem>_SKYSUB.fits`` and ``<stem>_SKYMODEL.fits`` pair and return the SKYSUB path*

    The sky model holds every line of the scene and the sky-subtracted frame is flat, so the lines are found only
    when the model is added back. ``badPixel`` is flagged in the ``QUAL`` extension.
    """
    if not path.name.endswith("_SKYSUB.fits"):
        raise ValueError("the path must end in _SKYSUB.fits")
    model: npt.NDArray[np.float64] = scene.data - BACKGROUND
    quality = np.zeros(scene.data.shape, dtype=np.int16)
    if badPixel is not None:
        quality[badPixel] = 1
    primary = fits.PrimaryHDU(np.full(scene.data.shape, BACKGROUND, dtype=np.float32), header=header)
    fits.HDUList([primary, fits.ImageHDU(quality, name="QUAL")]).writeto(path)
    modelPath = path.with_name(path.name.replace("_SKYSUB.fits", "_SKYMODEL.fits"))
    fits.PrimaryHDU(model.astype(np.float32)).writeto(modelPath)
    return path


def raw_header(
    arm: str = "VIS",
    slit: str = "SLIT1.0",
    mjd: float = 60000.4,
    tech: str = "ECHELLE,SLIT",
    dprType: str = "LAMP,WAVE",
    binning: tuple[int, int] = (1, 1),
) -> fits.Header:
    """*the header of a raw frame with the cards the SPEC_RES script reads*"""
    header = fits.Header()
    header["ESO SEQ ARM"] = arm
    header["ESO INS VISE NAME" if arm == "VIS" else "ESO INS NISE NAME"] = slit
    header["ESO DPR TECH"] = tech
    header["ESO DPR TYPE"] = dprType
    header["ESO DET BINX"], header["ESO DET BINY"] = binning
    header["MJD-OBS"] = mjd
    return header


def sky_header(slit: str = "SLIT1.0", mjd: float = 60000.45, arm: str = "VIS") -> fits.Header:
    """*the header of a reduced stare frame: the arm, slit and MJD of the science frame*"""
    return raw_header(arm=arm, slit=slit, mjd=mjd, tech="ECHELLE,SLIT,STARE", dprType="OBJECT")
