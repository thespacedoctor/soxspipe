"""Factory for a Phase 3 `SCIENCE.SPECTRUM` product written through the real Phase 3 table writer."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Column, Table

from soxspipe.commonutils.phase3 import PER_ARM_MODE, write_fits_table_to_disk
from soxspipe.commonutils.phase3_check import FAILURE, check_phase3

from .settings import pipeline_settings

PRODUCT_NAME = "SOXS_20260102T030405_STARE_ph3.fits"
SPECTRUM_LENGTH = 6
WAVE_START_NM = 600.0
WAVE_STEP_NM = 0.5
FLUX_UNIT = "erg / (s cm2 Angstrom)"
AIR_WAVELENGTH_UCD = "em.wl;obs.atmos"
SPEC_UTYPE = "spec:Data.SpectralAxis.Value"
FLUX_UTYPE = "spec:Data.FluxAxis.Value"
ERR_UTYPE = "spec:Data.FluxAxis.Accuracy.StatError"
ESO_FLUX_UTYPE = "eso:Data.FluxAxis.Value"
ESO_ERR_UTYPE = "eso:Data.FluxAxis.Accuracy.StatError"
MAIN_FLUX_UCD = "phot.flux.density;em.wl;meta.main"
MAIN_ERR_UCD = "stat.error;phot.flux.density;em.ql;meta.main"
OTHER_FLUX_UCD = "phot.flux.density;em.wl"
OTHER_ERR_UCD = "stat.error;phot.flux.density;em.ql"

# KEYWORDS OF THE PRIMARY HEADER, WITH A VALUE FOR EVERY ONE THE PHASE 3 STANDARD LISTS FOR A SPECTRUM
PRIMARY_CARDS: dict[str, object] = {
    "PRODCATG": "SCIENCE.SPECTRUM",
    "TELESCOP": "ESO-NTT",
    "INSTRUME": "SOXS",
    "OBJECT": "SYNTHETIC TARGET",
    "RA": 150.1,
    "DEC": -27.5,
    "RADESYS": "FK5",
    "EQUINOX": 2000.0,
    "EXPTIME": 600.0,
    "TEXPTIME": 600.0,
    "MJD-OBS": 61000.1,
    "MJD-END": 61000.108,
    "PROG_ID": "115.ABCD.001",
    "OBID1": 42,
    "NCOMBINE": 2,
    "OBSTECH": "ECHELLE,SLIT,STARE",
    "FLUXCAL": "ABSOLUTE",
    "PROV1": "SOXS.2026-01-02T03:04:05.678",
    "PROV2": "SOXS.2026-01-02T03:14:05.678",
    "SPECSYS": "BARYCENT",
    "CONTNORM": False,
    "TOT_FLUX": False,
    "FLUXERR": -2.0,
    "WAVELMIN": WAVE_START_NM,
    "WAVELMAX": WAVE_START_NM + WAVE_STEP_NM * (SPECTRUM_LENGTH - 1),
    "SPEC_BIN": WAVE_STEP_NM,
    "SNR": 12.5,
    "SPEC_RES": 5500.0,
    "SPEC_ERR": 0.002,
    "SPEC_SYE": 0.003,
    "RA_ERR": 0.0001,
    "DEC_ERR": 0.0001,
    "LAMRMS": 0.004,
    "LAMNLIN": 120,
}

# KEYWORDS OF THE TABLE HEADER THAT THE WRITER DOES NOT SET
EXTENSION_CARDS: dict[str, object] = {
    "NELEM": SPECTRUM_LENGTH,
    "VOCLASS": "SPECTRUM v2.0",
    "VOPUB": "ESO/SAF",
    "TITLE": "SYNTHETIC TARGET spectrum",
    "OBJECT": "SYNTHETIC TARGET",
    "RA": 150.1,
    "DEC": -27.5,
    "APERTURE": 0.0003,
    "TELAPSE": 691.2,
    "TMID": 61000.104,
    "SPEC_VAL": 600.0 + WAVE_STEP_NM * (SPECTRUM_LENGTH - 1) / 2,
    "SPEC_BW": WAVE_STEP_NM * (SPECTRUM_LENGTH - 1),
}

# THE COLUMNS IN ORDER: NAME, UNIT, TUCD, TUTYP
COLUMN_LAYOUT = (
    ("WAVE", "nm", AIR_WAVELENGTH_UCD, SPEC_UTYPE),
    ("FLUX", FLUX_UNIT, MAIN_FLUX_UCD, FLUX_UTYPE),
    ("ERR", FLUX_UNIT, MAIN_ERR_UCD, ERR_UTYPE),
    ("QUAL", "", "meta.code.qual", "eso:Data.FluxAxis.Quality"),
    ("FLUX_COUNTS", "adu", OTHER_FLUX_UCD, ESO_FLUX_UTYPE),
    ("ERR_COUNTS", "adu", OTHER_ERR_UCD, ESO_ERR_UTYPE),
)


def _default_arrays() -> dict[str, np.ndarray]:
    wave = WAVE_START_NM + WAVE_STEP_NM * np.arange(SPECTRUM_LENGTH)
    return {
        "WAVE": wave,
        "FLUX": np.full(SPECTRUM_LENGTH, 1.0e-16),
        "ERR": np.full(SPECTRUM_LENGTH, 1.0e-17),
        "QUAL": np.zeros(SPECTRUM_LENGTH, dtype=np.int32),
        "FLUX_COUNTS": np.full(SPECTRUM_LENGTH, 100.0),
        "ERR_COUNTS": np.full(SPECTRUM_LENGTH, 10.0),
    }


def _apply_overrides(cards: Mapping[str, object], overrides: Mapping[str, object] | None) -> dict[str, object]:
    """Return the cards with each override applied, an override of ``None`` removing the card."""
    merged = {**cards, **(overrides or {})}
    return {keyword: value for keyword, value in merged.items() if value is not None}


def _spectrum_table(
    columns: Mapping[str, np.ndarray] | None,
    dropColumns: tuple[str, ...],
    columnCards: Mapping[str, tuple[str, str, str]] | None,
) -> Table:
    """Build the one-row table, with the unit, UCD and UType of each column written as its header cards."""
    arrays = {**_default_arrays(), **(columns or {})}
    table = Table()
    tableCards: dict[str, object] = {}
    kept = [layout for layout in COLUMN_LAYOUT if layout[0] not in dropColumns]
    kept += [(name, "", "", "") for name in arrays if name not in {layout[0] for layout in COLUMN_LAYOUT}]
    for index, (name, unit, ucd, utype) in enumerate(kept, start=1):
        unit, ucd, utype = (columnCards or {}).get(name, (unit, ucd, utype))
        table[name] = Column(arrays[name] if arrays[name].ndim == 2 else [arrays[name]], unit=unit)
        tableCards[f"TCOMM{index}"] = f"{name} column"
        tableCards[f"TUCD{index}"] = ucd
        tableCards[f"TUTYP{index}"] = utype
        if name == "WAVE":
            tableCards[f"TDMIN{index}"] = float(np.min(arrays[name]))
            tableCards[f"TDMAX{index}"] = float(np.max(arrays[name]))
    table.meta.update(tableCards)
    return table


def phase3_spectrum_file(
    destination: Path,
    log: object,
    *,
    primary: Mapping[str, object] | None = None,
    extension: Mapping[str, object] | None = None,
    columns: Mapping[str, np.ndarray] | None = None,
    dropColumns: tuple[str, ...] = (),
    columnCards: Mapping[str, tuple[str, str, str]] | None = None,
    mode: str = PER_ARM_MODE,
) -> Path:
    """Write a good one-row `SCIENCE.SPECTRUM` product with the real Phase 3 table writer.

    ``primary`` and ``extension`` change header cards, and an override of ``None`` removes the card.
    ``columns`` replaces or adds column arrays, ``dropColumns`` leaves columns out, and ``columnCards`` replaces the
    ``(unit, TUCD, TUTYP)`` of a column.
    """
    header = fits.Header(_apply_overrides(PRIMARY_CARDS, primary))
    table = _spectrum_table(columns, dropColumns, columnCards)
    tableCards = _apply_overrides({**EXTENSION_CARDS, **table.meta}, extension)
    table.meta.clear()
    table.meta.update(tableCards)
    write_fits_table_to_disk(
        log=log,
        settings=pipeline_settings(destination.parent),
        header=header,
        tables=[table],
        filePath=destination,
        phase3=mode,
        extname="SPECTRUM",
    )
    return destination


def finding_messages(findings: list, severity: str) -> list[str]:
    """Return the messages of the findings that have the given severity."""
    return [finding.message for finding in findings if finding.severity == severity]


def failure_messages(tmp_path: Path, log: object, **factoryArguments: object) -> list[str]:
    """Write a product with ``phase3_spectrum_file(**factoryArguments)`` and return the failure messages."""
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, **factoryArguments)
    return finding_messages(check_phase3(product), FAILURE)


def rewrite_product(path: Path, change: Callable[[fits.HDUList], fits.HDUList | None]) -> Path:
    """Read a written product, let ``change`` break it, and write it back with fresh checksums."""
    with fits.open(path) as hdus:
        broken = change(hdus) or hdus
        broken.writeto(path, checksum=True, overwrite=True)
    return path


def add_cards(
    path: Path, primary: Mapping[str, object] | None = None, extension: Mapping[str, object] | None = None
) -> Path:
    """Add cards to the primary and extension headers of a written product."""

    def add(hdus: fits.HDUList) -> None:
        hdus[0].header.update(primary or {})
        hdus[1].header.update(extension or {})

    return rewrite_product(path, add)
