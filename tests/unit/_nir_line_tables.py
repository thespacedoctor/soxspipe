"""Helpers shared by the packaged SOXS NIR arc line table tests (DY-1970)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from astropy.table import Table

import soxspipe

CALIBRATIONS_DIR = Path(soxspipe.__file__).parent / "resources" / "static_calibrations" / "soxs"
MPH_TABLE = "SOXS_MPH_TAB_NIR_ArHgNeXe.fits"
SPH_TABLE = "ArHgNeXe_clean_within_2.0pixel.fits"

PINHOLES_PER_LINE = 9
WAVELENGTH_ATOL_NM = 5e-4
# DY-1970: ORDER 10 IS THE ONLY ORDER WHOSE ROWS CHANGED
EDITED_ORDER = 10
NEI_NM = 1957.7118
HGI_NM = 1970.018
XEI_NM = 1878.8146
# ROW COUNTS OF THE ORDER 10 BLOCK BEFORE THE EDIT (DEVELOP b9ab9446)
SPH_O10_ROWS_BEFORE = 184
MPH_O10_ROWS_BEFORE = 72
INTEGER_COLUMNS = ("order", "slit_index")
COLUMNS = ["ion", "wavelength", "order", "slit_index", "slit_position", "detector_x", "detector_y"]


def read_table(name: str) -> Table:
    return Table.read(CALIBRATIONS_DIR / name)


def ions(table: Table) -> np.ndarray:
    return np.array([i.decode().strip() if isinstance(i, bytes) else str(i).strip() for i in table["ion"]])


def set_rows(table: Table, ion: str, wavelength: float) -> Table:
    """Return the order 10 rows of one line, in file order."""
    isLine = (
        (np.asarray(table["order"]) == EDITED_ORDER)
        & (ions(table) == ion)
        & np.isclose(np.asarray(table["wavelength"], dtype=float), wavelength, rtol=0, atol=WAVELENGTH_ATOL_NM)
    )
    return table[isLine]


def order_checksums(table: Table) -> dict[int, str]:
    """Return one SHA-256 per order, over every column of its rows in file order."""
    checksums = {}
    for order in sorted({int(o) for o in table["order"]}):
        rows = table[np.asarray(table["order"]) == order]
        digest = hashlib.sha256()
        digest.update("|".join(ions(rows)).encode())
        for column in COLUMNS[1:]:
            dtype = ">i8" if column in INTEGER_COLUMNS else ">f8"
            digest.update(np.asarray(rows[column]).astype(dtype).tobytes())
        checksums[order] = digest.hexdigest()
    return checksums
