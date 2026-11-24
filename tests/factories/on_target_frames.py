"""Factory for the raw-frame headers and frame lists that the on-target frame accounting reads (DY-2356)."""

from __future__ import annotations

from typing import Any

from astropy.io import fits
from astropy.table import MaskedColumn, Table

TPL_START = "2026-01-02T03:04:05"
PROG_ID = "115.ABCD.001"
OBS_ID = 42
VIS_EXPTIME_SECONDS = 1800.0
NIR_DIT_SECONDS = 900.0
NIR_NDIT = 2


def raw_science_header(
    arm: str = "VIS",
    *,
    mjd: float = 61000.1,
    arcfile: str | None = "SOXS.2026-01-02T03:04:05.678.fits",
    exptime: float | None = None,
    dit: float | None = NIR_DIT_SECONDS,
    ndit: int | None = NIR_NDIT,
    obsId: int | None = OBS_ID,
    progId: str | None = PROG_ID,
    tplStart: str | None = TPL_START,
    **extra: Any,
) -> fits.Header:
    """Return the header of one raw science frame, with the cards the on-target frame list reads.

    A NIR frame (``arm="NIR"``) has ``EXPTIME`` equal to one DIT, as on the real instrument, and the DIT and NDIT cards
    the time on target comes from. A VIS frame has the whole exposure in ``EXPTIME``. A card given as ``None`` is left
    out. ``extra`` adds or overrides any other card.
    """
    isNir = arm.upper() == "NIR"
    if exptime is None:
        exptime = dit if isNir else VIS_EXPTIME_SECONDS
    cards: dict[str, Any] = {
        "ARCFILE": arcfile,
        "MJD-OBS": mjd,
        "EXPTIME": exptime,
        "HIERARCH ESO SEQ ARM": arm,
        "HIERARCH ESO OBS ID": obsId,
        "HIERARCH ESO OBS PROG ID": progId,
        "HIERARCH ESO TPL START": tplStart,
    }
    if isNir:
        cards["HIERARCH ESO DET SEQ1 DIT"] = dit
        cards["HIERARCH ESO DET NDIT"] = ndit
    cards.update(extra)
    return fits.Header({keyword: value for keyword, value in cards.items() if value is not None})


def raw_frame_summary(names: list[str], arm: str = "VIS", imageType: str = "OBJECT") -> Table:
    """Return the raw-frame summary of a recipe whose raw files are named ``names`` (``_pre`` is added to each).

    The names differ from the ``ARCFILE`` cards of the frames, as they do for a renamed or simulated file. Every row
    is a raw row: it has no ``ESO PRO TYPE``, which reaches pandas as NaN.
    """
    proType = MaskedColumn([""] * len(names), name="ESO PRO TYPE", mask=[True] * len(names))
    return Table(
        {
            "filename": [name.replace(".fits", "_pre.fits") for name in names],
            "file": [f"/raw/{name}" for name in names],
            "TYPE": [imageType] * len(names),
            "ARM": [arm] * len(names),
            "ESO PRO TYPE": proType,
            "ESO PRO CATG": MaskedColumn([""] * len(names), name="ESO PRO CATG", mask=[True] * len(names)),
        },
        masked=True,
    )
