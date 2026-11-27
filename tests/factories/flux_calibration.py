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


def bare_flux_calibration(**attributes: object) -> flux_calibration:
    """Return a ``flux_calibration`` built without its constructor, with ``attributes`` set on it.

    The constructor needs a workspace, so tests set the attributes ``calibrate`` reads. The defaults are the ones
    most tests share: no settings, no QC table, an empty products table and debug off.
    """
    calibration = object.__new__(flux_calibration)
    defaults = {"settings": {}, "debug": False, "qc": None, "products": pd.DataFrame()}
    for name, value in {**defaults, **attributes}.items():
        setattr(calibration, name, value)
    return calibration


def calibrate_and_capture_write(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    spectrum: object,
    *,
    arm: str = "NIR",
    binning: int = 1,
    responseWave: list[float] | None = None,
    extinction: np.ndarray | None = None,
    header: fits.Header | None = None,
    settings: dict[str, object] | None = None,
) -> dict[str, object]:
    """Run ``calibrate`` on ``spectrum`` and return the keyword arguments it hands to the writer."""
    import soxspipe.commonutils.toolkit as toolkit

    # A RESPONSE COVERING OTHER WAVELENGTHS IS FLAT, SO ONLY ITS RANGE MATTERS
    responsePath = (
        write_response(tmp_path, WAVES, RESPONSE)
        if responseWave is None
        else write_response(tmp_path, responseWave, [7.0, 7.0])
    )
    calibration = bare_flux_calibration(
        log=log,
        settings=settings if settings is not None else {},
        responseFunction=str(responsePath),
        extractedSpectrum=spectrum,
        exptime=EXPOSURE_TIME,
        airmass=1.4,
        extinctionPath="extinction.fits",
        arm=arm,
        header=header or fits.Header({"DATE-OBS": "2024-01-02T03:04:05", "ESO DET BINY": binning}),
        recipeName="soxs-stare",
        sofName="synthetic",
        productDir=str(tmp_path),
        kw=KEYWORDS.__getitem__,
    )
    captured: list[dict[str, object]] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk", lambda **kwargs: captured.append(kwargs)
    )
    if extinction is not None:
        monkeypatch.setattr(toolkit, "extinction_correction_factor", lambda wavelengths, path, airmass: extinction)
    calibration.calibrate()
    return captured[0]


def calibrate_spectrum(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    spectrum: object,
    **options: object,
) -> Table:
    """Run ``calibrate`` on ``spectrum`` and return the table it hands to the writer."""
    return calibrate_and_capture_write(tmp_path, log, monkeypatch, spectrum, **options)["tables"][0]
