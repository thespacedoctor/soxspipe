"""Shared synthetic VIS scene for the Horne-extraction sky-tuning tests (DY-1608, DY-1645)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import yaml
from astropy.io import fits

import soxspipe
from soxspipe.commonutils.horne_extraction import horne_extraction
from tests.factories import dispersion_map_fits

PIXEL_NM = 0.05
# SYNTHETIC DISPERSION MAP: y = 30000 + (30 o^2 - 150 o + 140) * wavelength, SO dy/dwavelength IS +20 px/nm IN
# ORDERS 1 AND 4 (r, i: +0.05 nm/px) AND -40 px/nm IN ORDERS 2 AND 3 (u, g: -0.025 nm/px)
RI_DISPERSION = 0.05
UG_DISPERSION = -0.025
# SKY LINES APPEAR THIS FAR RED OF THEIR CATALOGUE WAVELENGTHS: +0.6 px IN r AND +0.8 px IN i
R_LINE_OFFSET_NM = 0.03
I_LINE_OFFSET_NM = 0.04
SCIENCE_ABSROT = -19.76611
KEYWORDS = {"ABSROT": "ESO ADA ABSROT END", "PAR_ANG_START": "ESO TEL PARANG START"}
CATALOGUE = 501.0 + 1.5 * np.arange(12)


def dispersion_map(tmp_path: Path) -> str:
    coefficients = pd.DataFrame(
        {
            "axis": ["x", "y"],
            "order_deg": [2, 2],
            "wavelength_deg": [1, 1],
            "slit_deg": [0, 0],
            "c000": [100.0, 30000.0],
            "c010": [0.0, 140.0],
            "c100": [0.0, 0.0],
            "c110": [0.0, -150.0],
            "c200": [0.0, 0.0],
            "c210": [0.0, 30.0],
        }
    )
    return str(dispersion_map_fits(tmp_path / "dispersion_map.fits", coefficients=coefficients))


def sky_from_lines(wavelength: np.ndarray, centres: np.ndarray) -> np.ndarray:
    lines = sum(np.exp(-0.5 * ((wavelength - c) / 0.08) ** 2) for c in centres)
    return 1.0 + 20.0 * lines


def sky_spectrum(wavelength: np.ndarray, offsetNm: float) -> np.ndarray:
    return sky_from_lines(wavelength, CATALOGUE + offsetNm)


def build_orders() -> pd.DataFrame:
    wavelength = np.arange(500.0, 520.0, PIXEL_NM)
    frames = []
    for order, offset in [(1, R_LINE_OFFSET_NM), (2, None), (3, None), (4, I_LINE_OFFSET_NM)]:
        skyFlux = np.ones_like(wavelength) if offset is None else sky_spectrum(wavelength, offset)
        frames.append(
            pd.DataFrame({"order": order, "wavelengthMean": wavelength, "pixelScaleNm": PIXEL_NM, "skyFlux": skyFlux})
        )
    return pd.concat(frames, ignore_index=True)


def vis_parameters() -> dict:
    path = Path(soxspipe.__file__).parent / "resources" / "soxs_detector_parameters.yaml"
    return yaml.safe_load(path.read_text())["VIS"]


def build_extractor(
    log: object,
    tmp_path: Path,
    scienceHeader: dict,
    detectorParams: dict | None = None,
    inst: str = "SOXS",
) -> horne_extraction:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.debug = False
    extractor.axisB = "y"
    extractor.inst = inst
    extractor.dispersionMap = dispersion_map(tmp_path)
    extractor.detectorParams = vis_parameters() if detectorParams is None else detectorParams
    extractor.skySubtractedFrame = SimpleNamespace(header=fits.Header(scienceHeader))
    extractor.kw = KEYWORDS.__getitem__
    extractor.qc = pd.DataFrame()
    extractor.recipeName = "soxs-stare"
    extractor.dateObs = "2026-07-05T05:26:29"
    extractor._get_local_skylines_for_order = lambda *args: (CATALOGUE, CATALOGUE)
    return extractor


def night_header() -> dict:
    return {KEYWORDS["ABSROT"]: SCIENCE_ABSROT, KEYWORDS["PAR_ANG_START"]: -85.867}


def qc_by_name(extractor: horne_extraction, name: str) -> pd.Series:
    rows = extractor.qc.loc[extractor.qc["qc_name"] == name]
    assert len(rows) == 1, name
    return rows.iloc[0]


def applied_shift(before: pd.DataFrame, after: pd.DataFrame, order: int) -> float:
    mask = before["order"] == order
    return float((after.loc[mask, "wavelengthMean"] - before.loc[mask, "wavelengthMean"]).mean())
