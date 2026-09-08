"""Factories for temporary FITS and SOF files."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
from astropy.io import fits

from .frames import synthetic_ccd


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
