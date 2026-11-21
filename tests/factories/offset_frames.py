"""Factories for the headers of one offset template run, with the real ESO card names."""

from __future__ import annotations

import math
from collections.abc import Sequence

FIXED_SKY_TECHNIQUE = "ECHELLE,SLIT,OFFSET"
GENERIC_TECHNIQUE = "ECHELLE,SLIT,GENOFFSET"

# THE POINTING THE ACQUISITION LEFT THE TELESCOPE AT, IN DEGREES
BASE_RA_DEG = 244.0
BASE_DEC_DEG = -30.0
ARCSEC_PER_DEGREE = 3600.0

# ONE TEMPLATE RUN LASTS MINUTES; FRAMES OF A RUN ARE SPACED THIS FAR APART IN MJD
FRAME_SPACING_DAYS = 0.004
FIRST_MJD = 60000.0


def pointing_cards(eastArcsec: float, northArcsec: float) -> dict[str, float]:
    """Return the ``RA`` and ``DEC`` pointing cards of a telescope placed (east, north) arcsec from the base."""
    cosDec = math.cos(math.radians(BASE_DEC_DEG))
    return {
        "RA": BASE_RA_DEG + eastArcsec / cosDec / ARCSEC_PER_DEGREE,
        "DEC": BASE_DEC_DEG + northArcsec / ARCSEC_PER_DEGREE,
    }


def fixed_sky_run(
    steps: Sequence[tuple[float, float]],
    *,
    start: str = "2025-05-14T04:00:00",
    jitter: float = 1.0,
    acquisition: tuple[float, float] = (-9.6, -3.9),
    firstExpno: int = 1,
    firstMjd: float = FIRST_MJD,
    dprType: str = "OBJECT",
    technique: str = FIXED_SKY_TECHNIQUE,
) -> list[dict[str, object]]:
    """Return the headers of one FixedSkyOffset run, one per ``SEQ FIXOFF`` step.

    The telescope starts at ``acquisition`` (arcsec east and north of the base). Each frame carries the step the
    telescope made to reach it, as the template writes it, and a pointing that is the sum of the steps so far.
    """
    headers: list[dict[str, object]] = []
    east, north = acquisition
    for index, (stepRa, stepDec) in enumerate(steps):
        east, north = east + stepRa, north + stepDec
        headers.append(
            {
                **pointing_cards(east, north),
                "ESO TPL START": start,
                "ESO TPL EXPNO": firstExpno + index,
                "ESO SEQ JITTER WIDTH": jitter,
                "HIERARCH ESO SEQ FIXOFF RA": stepRa,
                "HIERARCH ESO SEQ FIXOFF DEC": stepDec,
                "ESO DPR TYPE": dprType,
                "ESO DPR TECH": technique,
                "MJD-OBS": firstMjd + index * FRAME_SPACING_DAYS,
            }
        )
    return headers


def generic_offset_run(
    pointings: Sequence[tuple[float, float]],
    dprTypes: Sequence[str],
    *,
    start: str = "2025-12-01T02:43:49",
    firstExpno: int = 1,
    firstMjd: float = FIRST_MJD,
) -> list[dict[str, object]]:
    """Return the headers of one GenericOffset run: a pointing (arcsec east, north of the base) and ``DPR TYPE`` each.

    These frames carry no ``SEQ FIXOFF`` or ``SEQ JITTER WIDTH`` card.
    """
    return [
        {
            **pointing_cards(east, north),
            "ESO TPL START": start,
            "ESO TPL EXPNO": firstExpno + index,
            "ESO DPR TYPE": dprType,
            "ESO DPR TECH": GENERIC_TECHNIQUE,
            "MJD-OBS": firstMjd + index * FRAME_SPACING_DAYS,
        }
        for index, ((east, north), dprType) in enumerate(zip(pointings, dprTypes, strict=True))
    ]
