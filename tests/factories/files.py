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
    fits.PrimaryHDU(data=np.asarray(frame.data), header=frame.header).writeto(destination)
    return destination


def harvestable_raw_fits(destination: Path, *, includeDprType: bool = True) -> Path:
    """Write a raw bias frame carrying every header card the data organiser indexes."""
    raw_fits(destination)
    with fits.open(destination, mode="update") as hdus:
        header = hdus[0].header
        header["MJD-OBS"] = 60311.75
        header["ESO DET3 EXPO TIME"] = 8.0
        header["ESO DET BINX"] = 1
        header["ESO DET BINY"] = 2
        header["ESO TPL ID"] = "SOXS_cal_bias"
        header["ESO INS ACFW ID"] = "g"
        header["ESO INS VISE NAME"] = "SLIT_1.0"
        header["ESO DET3 CAM NAME"] = "VIS"
        header["ESO ADA ABSROT END"] = 12.5
        if not includeDprType:
            del header["ESO DPR TYPE"]
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


LZW_MAGIC = b"\x1f\x9d"
LZW_BLOCK_MODE = 0x80
LZW_MAX_BITS = 16
LZW_FIRST_FREE_CODE = 257
LZW_CODES_PER_GROUP = 8


def lzw_compress(data: bytes) -> bytes:
    """Compress bytes into the Unix ``compress`` (``.Z``) LZW format.

    The output matches ``compress -b16`` for inputs that never trigger its
    dictionary reset: block mode is declared but no CLEAR code is ever sent,
    which every ``.Z`` reader accepts. It lets tests build ``.fits.Z`` fixtures
    without the ``compress`` command, which CI runners do not have.
    """
    header = LZW_MAGIC + bytes([LZW_BLOCK_MODE | LZW_MAX_BITS])
    if not data:
        return header
    return header + _pack_lzw_codes(_lzw_codes(data))


def _lzw_codes(data: bytes) -> list[tuple[int, int]]:
    """Return the LZW ``(code, bitWidth)`` sequence for ``data``."""
    maxMaxCode = 1 << LZW_MAX_BITS
    nBits = 9
    maxCode = (1 << nBits) - 1
    freeCode = LZW_FIRST_FREE_CODE
    table: dict[int, int] = {}
    codes: list[tuple[int, int]] = []

    def emit(code: int) -> None:
        nonlocal nBits, maxCode
        codes.append((code, nBits))
        # THE WIDTH GROWS AFTER THE CODE THAT FILLS THE CURRENT RANGE, AS IN compress(1)
        if freeCode > maxCode:
            nBits += 1
            maxCode = maxMaxCode if nBits == LZW_MAX_BITS else (1 << nBits) - 1

    prefix = data[0]
    for byte in data[1:]:
        key = (prefix << 8) | byte
        if key in table:
            prefix = table[key]
            continue
        emit(prefix)
        prefix = byte
        if freeCode < maxMaxCode:
            table[key] = freeCode
            freeCode += 1
    emit(prefix)
    return codes


def _pack_lzw_codes(codes: list[tuple[int, int]]) -> bytes:
    """Pack LZW codes least-significant bit first, in compress(1) groups.

    Codes are written in groups of eight. When the code width grows, the
    unfinished group is padded to its full size, ``bitWidth`` bytes, because
    readers skip to the next group boundary at each width change.
    """
    packed = bytearray()
    group = 0
    groupBits = 0
    groupCount = 0
    groupWidth = codes[0][1]
    for code, width in codes:
        if width != groupWidth:
            if groupCount:
                packed += group.to_bytes(groupWidth, "little")
            group, groupBits, groupCount, groupWidth = 0, 0, 0, width
        group |= code << groupBits
        groupBits += width
        groupCount += 1
        if groupCount == LZW_CODES_PER_GROUP:
            packed += group.to_bytes(width, "little")
            group, groupBits, groupCount = 0, 0, 0
    if groupCount:
        packed += group.to_bytes((groupBits + 7) // 8, "little")
    return bytes(packed)


def lzw_compressed_fits(source: Path, destination: Path) -> Path:
    """Write an LZW-compressed (``.fits.Z``) copy of ``source`` and return its path."""
    destination.write_bytes(lzw_compress(source.read_bytes()))
    return destination


def sof_file(destination: Path, members: Iterable[tuple[Path, str]]) -> Path:
    """Write a deterministic set-of-files inventory and return its path."""
    rows = [f"{path} {category}" for path, category in members]
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return destination
