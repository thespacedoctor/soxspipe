"""Shared synthetic SOXS NIR scene for the Horne-extraction global sky-tuning tests (DY-1968).

Eight 20 nm orders (11 to 18) at 0.05 nm per sample. Every order has the same signed dispersion, -0.05 nm/px, and holds
ten catalogue lines. The sky lines sit where a smooth detector motion puts them: ``TRUE_A0 + TRUE_A1 * m`` px, with
``m`` the order normalised to [-1, 1].
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits

import soxspipe
from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.nir_sky_refit import read_sky_tuning_settings
from soxspipe.commonutils.ug_wavelength_transfer import detector_px_to_nm_shift
from tests.factories import dispersion_map_fits
from tests.unit._horne_vis_scene import KEYWORDS, night_header

ORDERS = list(range(11, 19))
PIXEL_NM = 0.05
DISPERSION = -0.05
ORDER_WIDTH_NM = 20.0
LINES_PER_ORDER = 10
LINE_HEIGHT = 200.0
LINE_SIGMA_NM = 0.08
# THE DETECTOR MOTION (px) OF THE SKY LINES: TRUE_A0 + TRUE_A1 * m
TRUE_A0 = -3.0
TRUE_A1 = -0.4
ORDER_MID = 14.5
# THE TWO CATALOGUE COMPONENTS OF AN UNRESOLVED OH LAMBDA-DOUBLET ARE THIS FAR APART (nm)
BLEND_GAP_NM = 0.03
ORDER_HALF = 3.5


def order_start_nm(order: int) -> float:
    return 1000.0 + 30.0 * (order - 11)


def true_shift_px(order: int, offsets: dict | None = None) -> float:
    """*the detector motion (px) of an order's sky lines, plus the genuine offset of its arc solution if it has one*"""
    return TRUE_A0 + TRUE_A1 * (order - ORDER_MID) / ORDER_HALF + (offsets or {}).get(order, 0.0)


def true_shift_nm(order: int, offsets: dict | None = None) -> float:
    return detector_px_to_nm_shift(true_shift_px(order, offsets), DISPERSION)


def catalogue_nm(order: int) -> np.ndarray:
    return order_start_nm(order) + 1.5 + 1.8 * np.arange(LINES_PER_ORDER)


def dispersion_map(tmp_path: Path, arcRmsPx: float | None) -> str:
    coefficients = pd.DataFrame(
        {
            "axis": ["x", "y"],
            "order_deg": [0, 0],
            "wavelength_deg": [1, 1],
            "slit_deg": [0, 0],
            "c000": [100.0, 30000.0],
            "c010": [0.0, DISPERSION**-1],
        }
    )
    path = str(dispersion_map_fits(tmp_path / "dispersion_map.fits", coefficients=coefficients))
    if arcRmsPx is not None:
        fits.setval(path, "ESO QC Y RES MEDIAN", value=arcRmsPx)
    return path


def sky_flux(wavelength: np.ndarray, centres: np.ndarray) -> np.ndarray:
    lines = sum(np.exp(-0.5 * ((wavelength - c) / LINE_SIGMA_NM) ** 2) for c in centres)
    return 1.0 + LINE_HEIGHT * lines


def sky_line_nm(order: int, blendedFirstLine: bool = False) -> np.ndarray:
    """*the wavelengths of the sky lines the order holds; with ``blendedFirstLine`` the first is an unresolved blend*"""
    lines = catalogue_nm(order).copy()
    if blendedFirstLine:
        lines[0] += BLEND_GAP_NM / 2
    return lines


def build_orders(
    linelessOrders=(), linesOnly=None, offsets=None, lineCounts=None, blendedFirstLine=False
) -> pd.DataFrame:
    """*the extracted NIR orders; ``linelessOrders`` have a flat sky, ``linesOnly`` keeps lines in only those orders,
    ``offsets`` is ``{order: genuine offset of the arc solution (px)}``,
    ``lineCounts`` is ``{order: lines in its sky}``*"""
    frames = []
    for order in ORDERS:
        wavelength = order_start_nm(order) + PIXEL_NM * np.arange(int(ORDER_WIDTH_NM / PIXEL_NM))
        hasLines = order not in linelessOrders and (linesOnly is None or order in linesOnly)
        observed = sky_line_nm(order, blendedFirstLine)[: (lineCounts or {}).get(order, LINES_PER_ORDER)]
        observed = observed - true_shift_nm(order, offsets)
        skyFlux = sky_flux(wavelength, observed) if hasLines else np.ones_like(wavelength)
        frames.append(
            pd.DataFrame({"order": order, "wavelengthMean": wavelength, "pixelScaleNm": PIXEL_NM, "skyFlux": skyFlux})
        )
    return pd.concat(frames, ignore_index=True)


def catalogue_with_decoys(decoyOrders, offsetNm: float = 0.47):
    """*a ``_get_local_skylines_for_order`` double; ``decoyOrders`` list a decoy ``offsetNm`` past each line*"""

    def lines(wmin, wmax, order, column):
        real = catalogue_nm(order)
        real = np.concatenate([real, real + offsetNm]) if order in decoyOrders else real
        inRange = real[(real >= wmin) & (real <= wmax)]
        return inRange, inRange

    return lines


def nir_parameters() -> dict:
    path = Path(soxspipe.__file__).parent / "resources" / "soxs_detector_parameters.yaml"
    return yaml.safe_load(path.read_text())["NIR"]


def build_extractor(
    log: object,
    tmp_path: Path,
    sky_tuning: dict | None = None,
    inst: str = "SOXS",
    arm: str = "NIR",
    arcRmsPx: float | None = 0.2,
    decoyOrders=(),
) -> horne_extraction:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.debug = False
    # THE DISPERSION DIRECTION OF THIS SYNTHETIC MAP, WHICH IS WHAT axisB NAMES
    extractor.axisB = "y"
    extractor.inst = inst
    extractor.arm = arm
    extractor.dispersionMap = dispersion_map(tmp_path, arcRmsPx)
    extractor.detectorParams = nir_parameters()
    extractor.skySubtractedFrame = SimpleNamespace(header=fits.Header(night_header()))
    extractor.scienceHeaders = [extractor.skySubtractedFrame.header]
    extractor.kw = KEYWORDS.__getitem__
    extractor.qc = pd.DataFrame()
    extractor.recipeName = "soxs-stare"
    extractor.dateObs = "2026-07-05T05:26:29"
    extractor.nirSkyTuning = read_sky_tuning_settings(sky_tuning)
    extractor._get_local_skylines_for_order = catalogue_with_decoys(decoyOrders)
    return extractor


def applied_shift(before: pd.DataFrame, after: pd.DataFrame, order: int) -> float:
    mask = before["order"] == order
    return float((after.loc[mask, "wavelengthMean"] - before.loc[mask, "wavelengthMean"]).mean())


def catalogue_table(blendedFirstLine: bool = False) -> pd.DataFrame:
    """*the isolated-line catalogue table of the scene*

    It lists every line of every order, and with ``blendedFirstLine`` both components of the first line's doublet.
    """
    tables = []
    for order in ORDERS:
        lines = list(catalogue_nm(order))
        if blendedFirstLine:
            lines.append(lines[0] + BLEND_GAP_NM)
        tables.append(pd.DataFrame({"WAVELENGTH": lines, "FLUX": 300.0, "ISOLATED": True}))
    return pd.concat(tables, ignore_index=True)


def build_catalogue_extractor(
    log: object, tmp_path: Path, sky_tuning: dict | None = None, blendedFirstLine: bool = False
) -> horne_extraction:
    """*an extractor that runs the real ``_get_local_skylines_for_order`` against the scene's catalogue table*"""
    extractor = build_extractor(log, tmp_path, sky_tuning=sky_tuning)
    del extractor._get_local_skylines_for_order
    extractor.skylinesDF = catalogue_table(blendedFirstLine)
    return extractor
