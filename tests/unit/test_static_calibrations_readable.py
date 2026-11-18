"""Every packaged static calibration FITS file can be opened and read in full."""

from __future__ import annotations

from pathlib import Path

import pytest
from astropy.io import fits

from soxspipe.commonutils.getpackagepath import getpackagepath

pytestmark = pytest.mark.unit

STATIC_CALIBRATIONS_DIR = Path(getpackagepath()) / "resources" / "static_calibrations"
STATIC_FITS_FILES = sorted(STATIC_CALIBRATIONS_DIR.rglob("*.fits"))


def _relative_name(path: Path) -> str:
    return path.relative_to(STATIC_CALIBRATIONS_DIR).as_posix()


def test_static_calibrations_are_packaged() -> None:
    """The directory scan finds files, so the readability test below cannot pass on an empty list."""
    assert STATIC_FITS_FILES


@pytest.mark.parametrize("fitsPath", STATIC_FITS_FILES, ids=_relative_name)
def test_static_calibration_reads_every_hdu(fitsPath: Path) -> None:
    """Reading the data of each HDU raises nothing, and a failure names the file."""
    try:
        with fits.open(fitsPath, memmap=False) as hduList:
            for hdu in hduList:
                _ = hdu.data
    except Exception as err:
        pytest.fail(f"{_relative_name(fitsPath)} cannot be read: {type(err).__name__}: {err}")
