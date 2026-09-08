"""Factories for deterministic FITS headers and CCDData frames."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from astropy import units as u
from astropy.io.fits import Header
from astropy.nddata import CCDData, StdDevUncertainty

_KEYWORDS = {
    "soxs": {
        "DATE_OBS": "DATE-OBS",
        "DET_READ_SPEED": "ESO DET READ CURID",
        "DPR_CATG": "ESO DPR CATG",
        "DPR_TECH": "ESO DPR TECH",
        "DPR_TYPE": "ESO DPR TYPE",
        "EXPTIME": "EXPTIME",
        "INSTRUME": "INSTRUME",
        "OBJECT": "OBJECT",
        "OBS_ID": "ESO OBS ID",
        "SEQ_ARM": "ESO SEQ ARM",
    },
    "xsh": {
        "DATE_OBS": "DATE-OBS",
        "DET_READ_SPEED": "ESO DET READ SPEED",
        "DPR_CATG": "ESO DPR CATG",
        "DPR_TECH": "ESO DPR TECH",
        "DPR_TYPE": "ESO DPR TYPE",
        "EXPTIME": "EXPTIME",
        "INSTRUME": "INSTRUME",
        "OBJECT": "OBJECT",
        "OBS_ID": "ESO OBS ID",
        "SEQ_ARM": "ESO SEQ ARM",
    },
}


def instrument_header(
    *,
    instrument: str = "soxs",
    arm: str = "VIS",
    overrides: Mapping[str, object] | None = None,
) -> Header:
    """Return a fresh header using the instrument's resolved FITS keywords."""
    normalizedInstrument = instrument.lower()
    keywordMap = _KEYWORDS[normalizedInstrument]
    semanticValues: dict[str, object] = {
        "DATE_OBS": "2024-01-02T03:04:05.678",
        "DET_READ_SPEED": 2 if normalizedInstrument == "soxs" else "400k",
        "DPR_CATG": "CALIB",
        "DPR_TECH": "IMAGE",
        "DPR_TYPE": "BIAS",
        "EXPTIME": 0.0,
        "INSTRUME": normalizedInstrument.upper(),
        "OBJECT": "SYNTHETIC TARGET",
        "OBS_ID": 42,
        "SEQ_ARM": arm.upper(),
    }
    semanticValues = {**semanticValues, **(overrides or {})}
    return Header(
        {keywordMap.get(key, key): value for key, value in semanticValues.items()}
    )


def synthetic_ccd(
    *,
    shape: tuple[int, int] = (32, 32),
    seed: int = 7,
    instrument: str = "soxs",
    arm: str = "VIS",
    prepared: bool = False,
    headerOverrides: Mapping[str, object] | None = None,
) -> CCDData:
    """Return a fresh deterministic frame with optional mask and uncertainty."""
    rng = np.random.default_rng(seed)
    data = rng.normal(loc=100.0, scale=2.0, size=shape)
    header = instrument_header(
        instrument=instrument,
        arm=arm,
        overrides=headerOverrides,
    )
    if not prepared:
        return CCDData(data, unit=u.electron, meta=header)

    mask = np.zeros(shape, dtype=bool)
    mask[0, 0] = True
    uncertainty = StdDevUncertainty(np.full(shape, 1.5), unit=u.electron)
    return CCDData(
        data,
        unit=u.electron,
        meta=header,
        mask=mask,
        uncertainty=uncertainty,
    )
