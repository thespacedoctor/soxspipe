"""*the on-target frames of a recipe and the frame-accounting cards they give a merged spectrum (DY-2356)*

An on-target frame is a frame whose light is in the product: a stare frame, a nod A or B frame, or an offset ON frame.
OFF frames, rejected frames and nod cycles or positions that were skipped are not on target. Each recipe builds the list
from the raw headers when it reads the frames, because the recipe edits those headers later. The list gives the
``EXPTIME``, ``TEXPTIME``, ``NCOMBINE``, ``MJD-OBS``, ``MJD-END``, ``OBIDi``, ``PROG_ID`` and ``PROVi`` cards of the
merged spectrum. The cards go on a copy of the product header. They never go on a frame header: flux calibration and
the response divide the counts by the exposure time of one frame, and that header is where they read it.

This module imports nothing from a recipe.

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

import copy
import math
import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

from soxspipe.commonutils.fits_headers import set_header_card
from soxspipe.commonutils.ndit import read_ndit

SECONDS_PER_DAY = 86400.0
FITS_SUFFIX = ".fits"
NIR_ARM = "NIR"
MULTI_PROGRAMME = "MULTI"

# THE CARDS A PRODUCT HEADER GETS. THESE ARE PHASE 3 KEYWORDS, SO THEY DO NOT CHANGE WITH THE INSTRUMENT
EXPTIME_CARD = "EXPTIME"
TEXPTIME_CARD = "TEXPTIME"
NCOMBINE_CARD = "NCOMBINE"
MJD_OBS_CARD = "MJD-OBS"
MJD_END_CARD = "MJD-END"
PROG_ID_CARD = "PROG_ID"

# THE CARDS AN EARLIER STAMP OR A RAW FRAME MAY HAVE LEFT IN A HEADER, WHICH A NEW STAMP REPLACES
STALE_CARD_PATTERN = re.compile(r"PROV\d+|OBID\d+|PROGID\d+")
STALE_CARDS = (PROG_ID_CARD, NCOMBINE_CARD, TEXPTIME_CARD, MJD_END_CARD)

KeywordLookup = Callable[..., str]


@dataclass(frozen=True)
class on_target_frame:
    """One on-target frame, as its raw header gave it.

    **Key Arguments:**

    - ``sourceName`` -- the name of the file the recipe read the frame from, which identifies the frame in the recipe
    - ``name`` -- the name the frame goes by in the product: its ``ARCFILE``, else its ``ORIGFILE``, else ``sourceName``
    - ``arcfile`` -- the ``ARCFILE`` card, or *None* when the frame has none
    - ``mjdObs`` -- the start of the exposure (MJD), or *None* when the frame has no readable ``MJD-OBS``
    - ``exposureSeconds`` -- the time on target in seconds, or *None* when the frame has no readable time
    - ``obsId`` -- the observation block ID, or *None*
    - ``progId`` -- the programme ID, or *None*
    - ``tplStart`` -- the start time of the template the frame belongs to, or *None*
    """

    sourceName: str
    name: str
    arcfile: str | None
    mjdObs: float | None
    exposureSeconds: float | None
    obsId: int | None
    progId: str | None
    tplStart: str | None

    @property
    def provenance(self) -> str:
        """*the name of the frame in a* ``PROVi`` *card: its name without the* ``.fits`` *suffix*"""
        return self.name.removesuffix(FITS_SUFFIX)

    @property
    def mjdEnd(self) -> float | None:  # noqa: N802
        """*the end of the exposure (MJD), or None when the start or the time on target is unknown*"""
        if self.mjdObs is None or self.exposureSeconds is None:
            return None
        return self.mjdObs + self.exposureSeconds / SECONDS_PER_DAY


def _text(header: Any, keyword: str) -> str | None:
    """*the stripped text of a header card, or None when the card is missing or empty*"""
    text = str(header.get(keyword, "")).strip()
    return text or None


def _number(header: Any, keyword: str) -> float | None:
    """*a header card as a finite float, or None when it is missing or not a finite number*"""
    try:
        value = float(header[keyword])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _whole_number(header: Any, keyword: str) -> int | None:
    """*a header card as an integer, or None when it is missing or not a whole number*"""
    value = _number(header, keyword)
    return int(value) if value is not None and value == int(value) else None


def _time_on_target(header: Any, kw: KeywordLookup, arm: str, sourceName: str, log: Any) -> float | None:
    """*the time on target of one frame in seconds: DIT times NDIT on the NIR arm, else* ``EXPTIME``

    A frame with no readable positive time gives None and a warning that names the card.
    """
    if str(arm).upper() == NIR_ARM:
        card = kw("DET_DIT")
        dit = _number(header, card)
        seconds = None if dit is None else dit * read_ndit(header, kw, log)
    else:
        card = kw("EXPTIME")
        seconds = _number(header, card)
    if seconds is None or seconds <= 0:
        log.warning(f"{sourceName} has no readable {card} card, so its time on target is unknown.")
        return None
    return seconds


def read_on_target_frame(header: Any, kw: KeywordLookup, sourceName: str, arm: str, log: Any) -> on_target_frame:
    """*read the on-target entry of one raw frame from its header*

    A missing card gives None for that field and a warning. A frame with no ``ARCFILE`` takes its ``ORIGFILE``, else
    ``sourceName``, as its name, with a warning.

    **Key Arguments:**

    - ``header`` -- the header of the raw frame, before the recipe edits it
    - ``kw`` -- the keyword lookup (``keyword_lookup(...).get``)
    - ``sourceName`` -- the name of the file the recipe read the frame from
    - ``arm`` -- the arm of the frame. ``NIR`` takes its time from DIT and NDIT, every other arm from ``EXPTIME``
    - ``log`` -- logger

    **Return:**

    - ``entry`` -- the ``on_target_frame``

    **Usage:**

    ```python
    from soxspipe.commonutils.on_target_frames import read_on_target_frame
    entry = read_on_target_frame(frame.header, self.kw, "SOXS.2026-01-02T03:04:05.678.fits", self.arm, self.log)
    ```
    """
    arcfile = _text(header, kw("ARCFILE"))
    name = arcfile or _text(header, kw("ORIGFILE")) or sourceName
    if arcfile is None:
        log.warning(f"{sourceName} has no {kw('ARCFILE')} card, so its provenance card names it {name}.")

    obsId = _whole_number(header, kw("OBS_ID"))
    progId = _text(header, kw("OBS_PROG_ID"))
    tplStart = _text(header, kw("TPL_START"))
    missing = [
        kw(alias)
        for alias, value in (("OBS_ID", obsId), ("OBS_PROG_ID", progId), ("TPL_START", tplStart))
        if value is None
    ]
    if missing:
        log.warning(f"{sourceName} has no readable {', '.join(missing)} card, so the product leaves it out.")

    return on_target_frame(
        sourceName=sourceName,
        name=name,
        arcfile=arcfile,
        mjdObs=_number(header, kw("MJDOBS")),
        exposureSeconds=_time_on_target(header, kw, arm, sourceName, log),
        obsId=obsId,
        progId=progId,
        tplStart=tplStart,
    )


def build_on_target_frames(
    headers: Sequence[Any], sourceNames: Sequence[str], kw: KeywordLookup, arm: str, log: Any
) -> tuple[on_target_frame, ...]:
    """*build the on-target list of the frames a recipe has read, sorted by observation time*

    **Key Arguments:**

    - ``headers`` -- the raw headers of the frames, before the recipe edits them
    - ``sourceNames`` -- the file name of each frame, in the order of ``headers``
    - ``kw`` -- the keyword lookup
    - ``arm`` -- the arm of the frames
    - ``log`` -- logger

    **Return:**

    - ``entries`` -- one ``on_target_frame`` per frame, earliest first. A frame with no ``MJD-OBS`` goes last

    **Raises:**

    - ``ValueError`` -- when ``headers`` and ``sourceNames`` differ in length, or two frames have one file name

    **Usage:**

    ```python
    entries = build_on_target_frames([f.header for f in frames], names, self.kw, self.arm, self.log)
    ```
    """
    if len(headers) != len(sourceNames):
        raise ValueError(
            f"The on-target frame list needs one file name for each header, but it was given {len(headers)} headers "
            f"and {len(sourceNames)} file names."
        )
    duplicates = sorted({name for name in sourceNames if list(sourceNames).count(name) > 1})
    if duplicates:
        raise ValueError(f"Two frames have the same file name, so they cannot be told apart: {', '.join(duplicates)}")
    entries = [
        read_on_target_frame(header, kw, sourceName, arm, log)
        for header, sourceName in zip(headers, sourceNames, strict=True)
    ]
    return tuple(sorted(entries, key=lambda entry: (entry.mjdObs is None, entry.mjdObs or 0.0)))


def frames_named(entries: Sequence[on_target_frame], sourceNames: Sequence[str]) -> tuple[on_target_frame, ...]:
    """*select the entries of the frames read from the given files*

    **Key Arguments:**

    - ``entries`` -- the on-target list
    - ``sourceNames`` -- the file names of the frames to select

    **Return:**

    - ``selected`` -- the entries of those frames, in the order of ``entries``

    **Raises:**

    - ``ValueError`` -- when a file name has no entry

    **Usage:**

    ```python
    positionA = frames_named(self.onTargetFrames, allFrameANames)
    ```
    """
    wanted = set(sourceNames)
    unknown = sorted(wanted - {entry.sourceName for entry in entries})
    if unknown:
        raise ValueError(f"The on-target frame list has no entry for: {', '.join(unknown)}")
    return tuple(entry for entry in entries if entry.sourceName in wanted)


def without_frames(
    entries: Sequence[on_target_frame], sourceNames: Sequence[str], reason: str, log: Any
) -> tuple[on_target_frame, ...]:
    """*leave frames out of the on-target list, with one warning that names them and says why*

    **Key Arguments:**

    - ``entries`` -- the on-target list
    - ``sourceNames`` -- the file names of the frames to leave out
    - ``reason`` -- why they are not on target, as the end of a sentence
    - ``log`` -- logger

    **Return:**

    - ``kept`` -- the entries that remain, as a new tuple. ``entries`` itself when no frame is left out

    **Raises:**

    - ``ValueError`` -- when a file name has no entry

    **Usage:**

    ```python
    self.onTargetFrames = without_frames(self.onTargetFrames, namesB, "no object trace was found in them", self.log)
    ```
    """
    dropped = frames_named(entries, sourceNames)
    if not dropped:
        return tuple(entries)
    log.warning(
        f"{len(dropped)} frame{'s' if len(dropped) > 1 else ''} left out of the on-target list because "
        f"{reason}: {', '.join(entry.name for entry in dropped)}."
    )
    return tuple(entry for entry in entries if entry not in dropped)


def _distinct(values: Sequence[Any]) -> list[Any]:
    """*the values that are not None, once each, in order of first appearance*"""
    return list(dict.fromkeys(value for value in values if value is not None))


def _time_cards(entries: Sequence[on_target_frame]) -> list[tuple[str, Any, str]]:
    """*the exposure and the start and end time cards of a list, or none when a time is unknown*"""
    cards = []
    started = [entry for entry in entries if entry.mjdObs is not None]
    if started:
        cards.append((MJD_OBS_CARD, min(entry.mjdObs for entry in started), "[d] Start of the first frame"))
    times = [entry.exposureSeconds for entry in entries]
    if None in times:
        return cards
    total = float(sum(times))
    cards[:0] = [
        (EXPTIME_CARD, total, "[s] Total integration time on target"),
        (TEXPTIME_CARD, total, "[s] Total integration time on target"),
    ]
    ends = [entry.mjdEnd for entry in started]
    if ends:
        cards.append((MJD_END_CARD, max(ends), "[d] End of the last frame to end"))
    return cards


def frame_accounting_cards(entries: Sequence[on_target_frame]) -> list[tuple[str, Any, str]]:
    """*the frame-accounting cards of an on-target list*

    The exposure cards (``EXPTIME``, ``TEXPTIME`` and ``MJD-END``) are left out when any frame has no time on target,
    because the sum would be too small. ``PROV<i>`` and ``OBID<i>`` have no zero padding, as ESO writes them.

    **Key Arguments:**

    - ``entries`` -- the on-target list

    **Return:**

    - ``cards`` -- ``(keyword, value, comment)`` tuples. Empty for an empty list

    **Usage:**

    ```python
    for keyword, value, comment in frame_accounting_cards(entries):
        set_header_card(header, keyword, value, comment)
    ```
    """
    if not entries:
        return []
    cards = _time_cards(entries)
    cards.append((NCOMBINE_CARD, len(entries), "Number of frames combined"))
    cards.extend(
        (f"OBID{index}", obsId, "Observation block ID")
        for index, obsId in enumerate(_distinct([entry.obsId for entry in entries]), start=1)
    )
    programmes = _distinct([entry.progId for entry in entries])
    if len(programmes) == 1:
        cards.append((PROG_ID_CARD, programmes[0], "ESO programme ID"))
    elif programmes:
        cards.append((PROG_ID_CARD, MULTI_PROGRAMME, "Several ESO programmes, see PROGIDi"))
        cards.extend((f"PROGID{index}", progId, "ESO programme ID") for index, progId in enumerate(programmes, start=1))
    cards.extend(
        (f"PROV{index}", entry.provenance, "Originating raw frame") for index, entry in enumerate(entries, start=1)
    )
    return cards


def with_frame_accounting(header: Any, entries: Sequence[on_target_frame] | None) -> Any:
    """*give a header the frame-accounting cards of an on-target list, on a copy*

    The cards an earlier stamp left (``PROV<i>``, ``OBID<i>``, ``PROGID<i>``, ``PROG_ID``, ``NCOMBINE``, ``TEXPTIME``
    and ``MJD-END``) go first. The raw ``EXPTIME`` goes too when the list has no total, so a single frame's time cannot
    pass for one.

    **Key Arguments:**

    - ``header`` -- the header to stamp. It is not changed
    - ``entries`` -- the on-target list

    **Return:**

    - ``stamped`` -- a stamped copy of ``header``, or ``header`` itself when ``entries`` is None or empty

    **Usage:**

    ```python
    header = with_frame_accounting(self.masterHeaderFrame.header, self.onTargetFrames)
    ```
    """
    if not entries:
        return header
    stamped = copy.deepcopy(header)
    for keyword in [key for key in stamped if STALE_CARD_PATTERN.fullmatch(key) or key in STALE_CARDS]:
        stamped.pop(keyword, None)
    cards = frame_accounting_cards(entries)
    if EXPTIME_CARD not in {keyword for keyword, _, _ in cards}:
        stamped.pop(EXPTIME_CARD, None)
    for keyword, value, comment in cards:
        set_header_card(stamped, keyword, value, comment)
    return stamped
