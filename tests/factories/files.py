"""Factories for temporary FITS and SOF files."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from .dataframes import dispersion_table
from .frames import synthetic_ccd


def dispersion_map_fits(
    destination: Path,
    *,
    coefficients: pd.DataFrame | None = None,
) -> Path:
    """Write a deterministic dispersion-map coefficient table."""
    table = dispersion_table() if coefficients is None else coefficients.copy()
    Table.from_pandas(table).write(destination, format="fits")
    return destination


def order_table_fits(
    destination: Path,
    *,
    polynomials: pd.DataFrame,
    metadata: pd.DataFrame,
) -> Path:
    """Write polynomial and metadata extensions in the order-table layout."""
    hdus = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(Table.from_pandas(polynomials.copy())),
            fits.BinTableHDU(Table.from_pandas(metadata.copy())),
        ]
    )
    hdus.writeto(destination)
    return destination


def raw_fits(
    destination: Path,
    *,
    shape: tuple[int, int] = (32, 32),
    seed: int = 7,
    instrument: str = "soxs",
) -> Path:
    """Write a deterministic primary-HDU raw frame and return its path."""
    frame = synthetic_ccd(shape=shape, seed=seed, instrument=instrument)
    fits.PrimaryHDU(data=np.asarray(frame.data), header=frame.header).writeto(
        destination
    )
    return destination


def prepared_fits(
    destination: Path,
    *,
    shape: tuple[int, int] = (32, 32),
    seed: int = 7,
    instrument: str = "soxs",
) -> Path:
    """Write a prepared FLUX/ERRS/QUAL FITS layout and return its path."""
    frame = synthetic_ccd(
        shape=shape,
        seed=seed,
        instrument=instrument,
        prepared=True,
    )
    hdus = frame.to_hdu(
        hdu_mask="QUAL",
        hdu_uncertainty="ERRS",
        hdu_flags=None,
    )
    hdus[0].name = "FLUX"
    hdus.writeto(destination)
    return destination


def sof_file(destination: Path, members: Iterable[tuple[Path, str]]) -> Path:
    """Write a deterministic set-of-files inventory and return its path."""
    rows = [f"{path} {category}" for path, category in members]
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return destination
