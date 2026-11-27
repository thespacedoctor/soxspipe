"""Factory for a `DISP_TAB` file that holds chosen `ESO QC` cards."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from astropy.io import fits


def disp_tab_with_qcs(directory: Path, **qcCards: Any) -> str:
    """Write ``DISP_TAB_VIS.fits`` with the given ``ESO QC`` cards in its primary header, and return its path.

    A card is given by its QC name, for example ``**{"WAVE RMS": 0.0036}``. No card is written for a name not given, so
    a map from before a QC existed is a call without it.
    """
    header = fits.Header({f"HIERARCH ESO QC {name}": value for name, value in qcCards.items()})
    path = directory / "DISP_TAB_VIS.fits"
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU.from_columns([])]).writeto(path, overwrite=True)
    return str(path)
