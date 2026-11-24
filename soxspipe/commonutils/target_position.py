"""*the target position of a merged spectrum, decoded from the telescope target cards (DY-2357)*

``ESO TEL TARG ALPHA`` holds the right ascension as a packed ``HHMMSS.TTT`` number and ``ESO TEL TARG DELTA`` the
declination as a packed ``+-DDMMSS.TTT`` number. The ``RA`` and ``DEC`` cards of a raw frame are where the telescope
pointed, which is not the target, so the product never uses them. The target position goes in as it is: J2000, with no
proper motion.

The decoders split the digits with integer arithmetic, so the seconds pick up no float error.

This module imports nothing from a recipe.

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

HOURS_TO_DEGREES = 15
MILLI = 1000
SECONDS_PER_MINUTE = 60
MINUTES_PER_HOUR = 60
HOURS_PER_DAY = 24
MAX_DECLINATION = 90
MILLISECONDS_PER_HOUR = MINUTES_PER_HOUR * SECONDS_PER_MINUTE * MILLI

RA_CARD = "RA"
DEC_CARD = "DEC"
RADESYS_CARD = "RADESYS"
EQUINOX_CARD = "EQUINOX"
TARGET_FRAME = "FK5"
TARGET_EQUINOX = 2000.0
POSITION_COMMENT = "[deg] Spectroscopic target position (J2000)"


def _packed_number(value: Any) -> float | None:
    """*a header value as a finite float, or None for a bool, a non-number or a value that is not finite*"""
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _split_packed(number: float) -> tuple[int, int, int]:
    """*split the size of a packed number into whole units, whole minutes, and milliseconds of seconds*

    ``3019.8`` gives ``(0, 30, 19800)``. The size is rounded to the nearest thousandth first.
    """
    thousandths = round(abs(number) * MILLI)
    wholeDigits, millis = divmod(thousandths, MILLI)
    unitsAndMinutes, seconds = divmod(wholeDigits, 100)
    units, minutes = divmod(unitsAndMinutes, 100)
    return units, minutes, seconds * MILLI + millis


def _is_valid_sexagesimal(minutes: int, milliseconds: int) -> bool:
    """*return True when the minutes and the seconds (in thousandths) of a packed value are each under 60*"""
    return minutes < MINUTES_PER_HOUR and milliseconds < SECONDS_PER_MINUTE * MILLI


def decode_packed_ra(value: Any) -> float | None:
    """*decode a packed* ``HHMMSS.TTT`` *right ascension to degrees*

    **Key Arguments:**

    - ``value`` -- the card value, a number or a numeric string

    **Return:**

    - ``degrees`` -- the right ascension in degrees, or None when the value is not a finite number, is negative, or has
      hours of 24 or more, minutes of 60 or more or seconds of 60 or more

    **Usage:**

    ```python
    from soxspipe.commonutils.target_position import decode_packed_ra
    ra = decode_packed_ra(header["ESO TEL TARG ALPHA"])
    ```
    """
    number = _packed_number(value)
    if number is None or number < 0:
        return None
    hours, minutes, milliseconds = _split_packed(number)
    if hours >= HOURS_PER_DAY or not _is_valid_sexagesimal(minutes, milliseconds):
        return None
    totalMilliseconds = (hours * MINUTES_PER_HOUR + minutes) * SECONDS_PER_MINUTE * MILLI + milliseconds
    return totalMilliseconds * HOURS_TO_DEGREES / MILLISECONDS_PER_HOUR


def decode_packed_dec(value: Any) -> float | None:
    """*decode a packed* ``+-DDMMSS.TTT`` *declination to degrees*

    The sign belongs to the whole value, so ``-000530.0`` is -0.091667 deg.

    **Key Arguments:**

    - ``value`` -- the card value, a number or a numeric string

    **Return:**

    - ``degrees`` -- the declination in degrees, or None when the value is not a finite number, or has degrees over 90,
      minutes of 60 or more or seconds of 60 or more

    **Usage:**

    ```python
    from soxspipe.commonutils.target_position import decode_packed_dec
    dec = decode_packed_dec(header["ESO TEL TARG DELTA"])
    ```
    """
    number = _packed_number(value)
    if number is None:
        return None
    degrees, minutes, milliseconds = _split_packed(number)
    if not _is_valid_sexagesimal(minutes, milliseconds):
        return None
    totalMilliseconds = (degrees * MINUTES_PER_HOUR + minutes) * SECONDS_PER_MINUTE * MILLI + milliseconds
    if totalMilliseconds > MAX_DECLINATION * MILLISECONDS_PER_HOUR:
        return None
    size = totalMilliseconds / MILLISECONDS_PER_HOUR
    return -size if number < 0 and size else size


def target_position_cards(entries: Sequence[Any], log: Any) -> list[tuple[str, Any, str]]:
    """*the target position cards of an on-target list: the position of the earliest frame that has one*

    **Key Arguments:**

    - ``entries`` -- the on-target list, earliest frame first. Each entry has ``targetRa`` and ``targetDec`` in degrees,
      or None
    - ``log`` -- logger

    **Return:**

    - ``cards`` -- ``(keyword, value, comment)`` tuples for ``RA``, ``DEC``, ``RADESYS`` and ``EQUINOX``. Empty, with
      one warning, when no frame has both coordinates

    **Usage:**

    ```python
    for keyword, value, comment in target_position_cards(entries, self.log):
        set_header_card(header, keyword, value, comment)
    ```
    """
    for entry in entries:
        if entry.targetRa is not None and entry.targetDec is not None:
            return [
                (RA_CARD, entry.targetRa, POSITION_COMMENT),
                (DEC_CARD, entry.targetDec, POSITION_COMMENT),
                (RADESYS_CARD, TARGET_FRAME, "Reference frame of RA and DEC"),
                (EQUINOX_CARD, TARGET_EQUINOX, "[yr] Equinox of RA and DEC"),
            ]
    log.warning(
        "No on-target frame has both ESO TEL TARG ALPHA and ESO TEL TARG DELTA, so the spectrum has no RA and DEC "
        "and the arm fails."
    )
    return []
