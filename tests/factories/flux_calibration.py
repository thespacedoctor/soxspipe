"""A synthetic flux calibration run, for tests of the calibration and of what feeds it."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.flux_calibration import flux_calibration
from soxspipe.commonutils.response_basis import RESPONSE_BASIS, RESPONSE_BASIS_KEYWORD

WAVES = [500.0, 600.0]
RESPONSE = [7.0, 8.0]
EXPOSURE_TIME = 10.0
KEYWORDS = {
    "DPR_CATG": "ESO DPR CATG",
    "DPR_TYPE": "ESO DPR TYPE",
    "DET_READ_SPEED": "ESO DET READ CURID",
    "CONAD": "ESO DET OUT1 CONAD",
    "GAIN": "ESO DET OUT1 GAIN",
    "RON": "ESO DET OUT1 RON",
    "SEQ_ARM": "ESO SEQ ARM",
    "PRO_TYPE": "ESO PRO TYPE",
    "PRO_CATG": "ESO PRO CATG",
    "PRO_AIRM_MEAN": "ESO PRO AIRM MEAN",
    "WIN_BINY": "ESO DET BINY",
}


def write_response(tmp_path: Path, wave: list[float], response: list[float]) -> Path:
    """Write a RESP product holding a RESPONSE curve with the current basis card."""
    curve = fits.table_to_hdu(Table({"WAVE": wave, "RESPONSE": response}))
    curve.name = "RESPONSE"
    primary = fits.PrimaryHDU(header=fits.Header({RESPONSE_BASIS_KEYWORD: RESPONSE_BASIS}))
    path = tmp_path / "SYNTHETIC_RESP.fits"
    fits.HDUList([primary, curve]).writeto(path)
    return path


def calibrate_spectrum(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    spectrum: object,
    *,
    arm: str = "NIR",
    binning: int = 1,
    responseWave: list[float] | None = None,
    extinction: np.ndarray | None = None,
) -> Table:
    """Run ``calibrate`` on ``spectrum`` and return the table it hands to the writer."""
    import soxspipe.commonutils.toolkit as toolkit

    # A RESPONSE COVERING OTHER WAVELENGTHS IS FLAT, SO ONLY ITS RANGE MATTERS
    responsePath = (
        write_response(tmp_path, WAVES, RESPONSE)
        if responseWave is None
        else write_response(tmp_path, responseWave, [7.0, 7.0])
    )
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = spectrum
    calibration.exptime = EXPOSURE_TIME
    calibration.airmass = 1.4
    calibration.extinctionPath = "extinction.fits"
    calibration.arm = arm
    calibration.debug = False
    calibration.header = fits.Header({"DATE-OBS": "2024-01-02T03:04:05", "ESO DET BINY": binning})
    calibration.recipeName = "soxs-stare"
    calibration.sofName = "synthetic"
    calibration.productDir = str(tmp_path)
    calibration.products = pd.DataFrame()
    calibration.kw = KEYWORDS.__getitem__
    captured: list[object] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.append(kwargs["tables"][0]),
    )
    if extinction is not None:
        monkeypatch.setattr(toolkit, "extinction_correction_factor", lambda wavelengths, path, airmass: extinction)
    calibration.calibrate()
    return captured[0]
