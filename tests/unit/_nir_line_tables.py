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
SLIT_INDEX_AT_CENTRE = 4
WAVELENGTH_ATOL_NM = 5e-4
POSITION_ATOL_PX = 1e-6
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
    """*read one packaged SOXS NIR line table*

    **Key Arguments:**

    - ``name`` -- file name of the table in the packaged SOXS static calibrations directory

    **Return:**

    - ``table`` -- the table, with every row in file order
    """
    return Table.read(CALIBRATIONS_DIR / name)


def ions(table: Table) -> np.ndarray:
    """*return the ion labels of a table as stripped strings*

    **Key Arguments:**

    - ``table`` -- a line table, or a selection of its rows

    **Return:**

    - ``labels`` -- one stripped ion label per row, in row order
    """
    return np.array([i.decode().strip() if isinstance(i, bytes) else str(i).strip() for i in table["ion"]])


def set_mask(table: Table, order: int, ion: str, wavelength: float, atol: float = WAVELENGTH_ATOL_NM) -> np.ndarray:
    """*flag the rows of one line in one order*

    **Key Arguments:**

    - ``table`` -- a line table
    - ``order`` -- the echelle order of the line
    - ``ion`` -- the ion label of the line
    - ``wavelength`` -- the wavelength of the line in nm
    - ``atol`` -- the largest distance in nm between a row wavelength and ``wavelength``

    **Return:**

    - ``mask`` -- boolean array, true for each row of the line
    """
    return (
        (np.asarray(table["order"]) == order)
        & (ions(table) == ion)
        & np.isclose(np.asarray(table["wavelength"], dtype=float), wavelength, rtol=0, atol=atol)
    )


def set_rows(table: Table, order: int, ion: str, wavelength: float, atol: float = WAVELENGTH_ATOL_NM) -> Table:
    """*select the rows of one line in one order*

    **Key Arguments:**

    - ``table`` -- a line table
    - ``order`` -- the echelle order of the line
    - ``ion`` -- the ion label of the line
    - ``wavelength`` -- the wavelength of the line in nm
    - ``atol`` -- the largest distance in nm between a row wavelength and ``wavelength``

    **Return:**

    - ``rows`` -- the rows of the line, in file order
    """
    return table[set_mask(table, order, ion, wavelength, atol)]


def positions_by_slit(rows: Table) -> dict[int, tuple[float, float]]:
    """*map the slit index of each row to its detector position*

    **Key Arguments:**

    - ``rows`` -- the rows of one line

    **Return:**

    - ``positions`` -- ``{slit_index: (detector_x, detector_y)}``
    """
    columns = zip(rows["slit_index"], rows["detector_x"], rows["detector_y"], strict=True)
    return {int(s): (float(x), float(y)) for s, x, y in columns}


def rows_checksum(rows: Table) -> str:
    """*hash every column of some rows, in row order*

    **Key Arguments:**

    - ``rows`` -- the rows to hash

    **Return:**

    - ``checksum`` -- the SHA-256 hex digest of the ion labels and of every numeric column
    """
    digest = hashlib.sha256()
    digest.update("|".join(ions(rows)).encode())
    for column in COLUMNS[1:]:
        dtype = ">i8" if column in INTEGER_COLUMNS else ">f8"
        digest.update(np.asarray(rows[column]).astype(dtype).tobytes())
    return digest.hexdigest()


def order_checksums(table: Table) -> dict[int, str]:
    """*hash the rows of every order of a table*

    **Key Arguments:**

    - ``table`` -- a line table

    **Return:**

    - ``checksums`` -- ``{order: checksum}``, one ``rows_checksum`` per order, over its rows in file order
    """
    orders = sorted({int(o) for o in table["order"]})
    return {order: rows_checksum(table[np.asarray(table["order"]) == order]) for order in orders}


def untouched_rows_checksum(table: Table, order: int, touchedSets: list[tuple[str, float]]) -> str:
    """*hash the rows of one order that the edit did not touch*

    **Key Arguments:**

    - ``table`` -- a line table
    - ``order`` -- the echelle order that was edited
    - ``touchedSets`` -- ``(ion, wavelength)`` of each line of that order that the edit added or moved

    **Return:**

    - ``checksum`` -- the ``rows_checksum`` of the rows of ``order`` outside ``touchedSets``, in file order
    """
    keep = np.asarray(table["order"]) == order
    for ion, wavelength in touchedSets:
        keep &= ~set_mask(table, order, ion, wavelength)
    return rows_checksum(table[keep])
