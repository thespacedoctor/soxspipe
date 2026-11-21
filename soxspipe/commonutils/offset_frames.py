"""*the ON and OFF frames of an offset SOF: the classification per template run and the balance check (DY-2354)*

An offset SOF holds the exposures of one or more template runs (one ``TPL START`` each). Each run is classified on its
own, by where the telescope pointed relative to the run's first exposure. The rule depends on the template:

- ``ECHELLE,SLIT,OFFSET`` (FixedSkyOffset and SpecphotStdOffset): a frame is ON when it is within the jitter width plus
  1 arcsec of the run's first exposure. The running sum of the ``SEQ FIXOFF`` steps must give the same positions.
- ``ECHELLE,SLIT,GENOFFSET`` (GenericOffset): ``DPR TYPE`` says SKY (OFF) or OBJECT (ON), and an OBJECT frame must be
  within 1 arcsec of the run's first OBJECT frame.

A frame that cannot be classified fails the SOF. Frames are never dropped silently.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import Any

ARCSEC_PER_DEGREE = 3600.0

FIXED_SKY_TECHNIQUE = "ECHELLE,SLIT,OFFSET"
GENERIC_TECHNIQUE = "ECHELLE,SLIT,GENOFFSET"

# THE DPR TYPE OF A SKY FRAME OF EITHER TEMPLATE. A DPR TYPE WITH THIS WORD IS OFF
SKY_TYPE = "SKY"

# THE DPR TYPE WORDS OF A FRAME THAT IS ON THE TARGET
ON_TYPE_WORDS = {"OBJECT", "STD"}

# A FIXEDSKYOFFSET FRAME IS ON THE OBJECT WHEN IT IS WITHIN THE JITTER WIDTH OF THE RUN'S FIRST EXPOSURE, PLUS THIS
# MARGIN
ON_POSITION_MARGIN_ARCSEC = 1.0

# THE POINTING OF A FRAME AND THE SUM OF THE SEQ FIXOFF STEPS UP TO IT MUST AGREE TO THIS DISTANCE
FIXOFF_CROSS_CHECK_ARCSEC = 1.0

# A GENERICOFFSET OBJECT FRAME MUST BE THIS CLOSE TO THE RUN'S FIRST OBJECT FRAME. THE TEMPLATE HAS NO JITTER PARAMETER
GENERIC_ON_TOLERANCE_ARCSEC = 1.0


class OffsetClassificationError(ValueError):
    """*the frames of an offset SOF cannot be classified as ON and OFF frames*"""


class frame_label(Enum):
    """The label of a classified frame. A frame that could not be classified has no label (None)."""

    ON = "ON"
    OFF = "OFF"


@dataclass(frozen=True)
class _offset_exposure:
    """One frame of a template run, with the header cards that the classification reads."""

    frame: Any
    filename: str
    start: str
    expno: int
    ra: float
    dec: float
    technique: str
    dprType: str | None = None
    jitter: float | None = None
    fixoffRa: float | None = None
    fixoffDec: float | None = None


def _read_card(
    header: Any, keyword: str, filename: str, problems: list[str], cast: Callable[[Any], Any] = float
) -> Any:
    """*read a header card, recording a problem when it is missing or unreadable*

    **Key Arguments:**

    - ``header`` -- the frame's header
    - ``keyword`` -- the header keyword to read
    - ``filename`` -- the frame's file name, which names the frame in a problem
    - ``problems`` -- the list that a problem message is added to
    - ``cast`` -- the type to convert the value to. Default *float*

    **Return:**

    - ``value`` -- the converted value, or None when the card is missing or unreadable
    """
    try:
        return cast(header[keyword])
    except KeyError:
        problems.append(f"{filename} has no {keyword} card.")
    except (TypeError, ValueError):
        problems.append(f"{filename} has an unreadable {keyword} card: {header[keyword]!r}.")
    return None


def _read_exposure(frame: Any, filename: str, kw: Callable[[str], str], problems: list[str]) -> _offset_exposure | None:
    """*read the cards that the classification of one frame needs*

    The cards depend on the template. FixedSkyOffset reads the jitter width and the ``SEQ FIXOFF`` steps, and
    GenericOffset reads ``DPR TYPE``.

    **Key Arguments:**

    - ``frame`` -- the science frame
    - ``filename`` -- the frame's file name
    - ``kw`` -- the keyword lookup
    - ``problems`` -- the list that a problem message is added to

    **Return:**

    - ``exposure`` -- the frame and its cards, or None when a card is missing or the template is not an offset one.
      Each reason is added to ``problems``
    """
    header = frame.header
    problemCount = len(problems)
    start = _read_card(header, kw("TPL_START"), filename, problems, str)
    expno = _read_card(header, kw("TPL_EXPNO"), filename, problems, int)
    ra = _read_card(header, kw("RA"), filename, problems)
    dec = _read_card(header, kw("DEC"), filename, problems)
    technique = _read_card(header, kw("DPR_TECH"), filename, problems, str)

    dprType = jitter = fixoffRa = fixoffDec = None
    if technique == FIXED_SKY_TECHNIQUE:
        jitter = _read_card(header, kw("NOD_JITTER_BOX"), filename, problems)
        fixoffRa = _read_card(header, kw("OFFSET_RA"), filename, problems)
        fixoffDec = _read_card(header, kw("OFFSET_DEC"), filename, problems)
    elif technique == GENERIC_TECHNIQUE:
        dprType = _read_card(header, kw("DPR_TYPE"), filename, problems, str)
    elif technique is not None:
        problems.append(
            f"{filename} has DPR TECH {technique}, which is neither {FIXED_SKY_TECHNIQUE} nor {GENERIC_TECHNIQUE}."
        )

    if len(problems) > problemCount:
        return None
    return _offset_exposure(
        frame=frame,
        filename=filename,
        start=start,
        expno=expno,
        ra=ra,
        dec=dec,
        technique=technique,
        dprType=dprType,
        jitter=jitter,
        fixoffRa=fixoffRa,
        fixoffDec=fixoffDec,
    )


def _sky_offset_arcsec(exposure: _offset_exposure, reference: _offset_exposure) -> tuple[float, float]:
    """*the (east, north) offset in arcsec of a frame's pointing from the reference frame's pointing*

    **Key Arguments:**

    - ``exposure`` -- the frame
    - ``reference`` -- the frame the offset is measured from

    **Return:**

    - ``east``, ``north`` -- the offset in arcsec, east on the sky (RA times cos(DEC)) and north
    """
    north = (exposure.dec - reference.dec) * ARCSEC_PER_DEGREE
    deltaRa = (exposure.ra - reference.ra + 180.0) % 360.0 - 180.0
    east = deltaRa * math.cos(math.radians(reference.dec)) * ARCSEC_PER_DEGREE
    return east, north


def _describe(exposure: _offset_exposure, start: str) -> str:
    """*name a frame, with its place in its template run, for an error message*

    **Key Arguments:**

    - ``exposure`` -- the frame
    - ``start`` -- the ``TPL START`` of the frame's run

    **Return:**

    - ``description`` -- the frame's file name, its exposure number and its run
    """
    return f"{exposure.filename} (exposure {exposure.expno} of the template run that started at {start})"


def _classify_fixed_sky(start: str, run: list[_offset_exposure], problems: list[str]) -> list[frame_label | None]:
    """*the ON frames of a FixedSkyOffset run: those near the run's first exposure*

    The run is cross-checked against the running sum of its ``SEQ FIXOFF`` steps, which is the position of each frame
    if the telescope did what the template commanded. A frame whose two positions differ by more than 1 arcsec is a
    problem.

    **Key Arguments:**

    - ``start`` -- the ``TPL START`` of the run
    - ``run`` -- the frames of the run, in exposure order
    - ``problems`` -- the list that a problem message is added to

    **Return:**

    - ``labels`` -- one ``frame_label`` (ON or OFF) per frame of ``run``
    """
    reference = run[0]
    labels = []
    stepEast, stepNorth = 0.0, 0.0
    referenceEast, referenceNorth = reference.fixoffRa, reference.fixoffDec
    for exposure in run:
        stepEast, stepNorth = stepEast + exposure.fixoffRa, stepNorth + exposure.fixoffDec
        east, north = _sky_offset_arcsec(exposure, reference)
        commandedEast, commandedNorth = stepEast - referenceEast, stepNorth - referenceNorth
        if math.hypot(east - commandedEast, north - commandedNorth) > FIXOFF_CROSS_CHECK_ARCSEC:
            problems.append(
                f"{_describe(exposure, start)} is at ({east:+.1f}, {north:+.1f}) arcsec from the first exposure by "
                f"its pointing, but at ({commandedEast:+.1f}, {commandedNorth:+.1f}) arcsec by the sum of the "
                "SEQ FIXOFF steps. The two must agree to within 1 arcsec."
            )
        isOn = math.hypot(east, north) <= exposure.jitter + ON_POSITION_MARGIN_ARCSEC
        labels.append(frame_label.ON if isOn else frame_label.OFF)
    return labels


def _classify_generic(start: str, run: list[_offset_exposure], problems: list[str]) -> list[frame_label | None]:
    """*label the frames of a GenericOffset run from ``DPR TYPE``, and check the ON frames stay put*

    A frame whose ``DPR TYPE`` has the word SKY is OFF, wherever the telescope pointed. An OBJECT frame is ON when it
    is within 1 arcsec of the run's first OBJECT frame. A frame of any other type, or an OBJECT frame farther away, is a
    problem and gets no label.

    **Key Arguments:**

    - ``start`` -- the ``TPL START`` of the run
    - ``run`` -- the frames of the run, in exposure order
    - ``problems`` -- the list that a problem message is added to

    **Return:**

    - ``labels`` -- one ``frame_label`` (ON or OFF) per frame of ``run``, or None for a frame that could not be
      classified
    """
    labels = []
    reference = None
    for exposure in run:
        words = set(exposure.dprType.upper().split(","))
        if SKY_TYPE in words:
            labels.append(frame_label.OFF)
        elif not ON_TYPE_WORDS & words:
            problems.append(
                f"{_describe(exposure, start)} has DPR TYPE {exposure.dprType}, which is not OBJECT or SKY."
            )
            labels.append(None)
        else:
            reference = reference or exposure
            east, north = _sky_offset_arcsec(exposure, reference)
            distance = math.hypot(east, north)
            if distance > GENERIC_ON_TOLERANCE_ARCSEC:
                problems.append(
                    f"{_describe(exposure, start)} is an OBJECT frame {distance:.1f} arcsec from the first OBJECT "
                    f"frame of its run, {reference.filename}. An OBJECT frame must be within "
                    f"{GENERIC_ON_TOLERANCE_ARCSEC:g} arcsec of it."
                )
            labels.append(frame_label.ON if distance <= GENERIC_ON_TOLERANCE_ARCSEC else None)
    return labels


def _check_balance(
    start: str,
    onFrames: list[_offset_exposure],
    offFrames: list[_offset_exposure],
    hasProblems: bool,
    problems: list[str],
) -> None:
    """*record a problem when a run has no ON frame, no OFF frame or unequal numbers of ON and OFF frames*

    A run that already has a problem can be unequal only because of it, so the unequal counts are then left unsaid. A
    missing ON or OFF frame is always said.

    **Key Arguments:**

    - ``start`` -- the ``TPL START`` of the run
    - ``onFrames`` -- the run's ON frames
    - ``offFrames`` -- the run's OFF frames
    - ``hasProblems`` -- whether the run already has a classification problem
    - ``problems`` -- the list that a problem message is added to
    """
    if not onFrames:
        problems.append(
            f"The template run that started at {start} has no ON frame and {len(offFrames)} OFF frames. "
            "An offset reduction needs frames on the target."
        )
    elif not offFrames:
        problems.append(
            f"The template run that started at {start} has {len(onFrames)} ON frames and no OFF frame. "
            "An offset reduction needs sky frames. Use soxs-stare for data taken without a sky position."
        )
    elif len(onFrames) != len(offFrames) and not hasProblems:
        problems.append(
            f"Found {len(onFrames)} ON frames and {len(offFrames)} OFF frames in the template run that started at "
            f"{start}. The number of ON and OFF frames must be the same for offset reductions."
        )


def _check_one_template(start: str, run: list[_offset_exposure], problems: list[str]) -> bool:
    """*record a problem when the frames of one run were taken with more than one offset template*

    **Key Arguments:**

    - ``start`` -- the ``TPL START`` of the run
    - ``run`` -- the frames of the run
    - ``problems`` -- the list that a problem message is added to

    **Return:**

    - ``isOneTemplate`` -- whether the run holds frames of one template only
    """
    techniques = sorted({exposure.technique for exposure in run})
    if len(techniques) > 1:
        parts = [
            f"{technique} ({', '.join(e.filename for e in run if e.technique == technique)})"
            for technique in techniques
        ]
        problems.append(
            f"The template run that started at {start} mixes offset templates: {' and '.join(parts)}. "
            "Each run must hold frames of one template."
        )
    return len(techniques) == 1


def _split_run(
    start: str, run: list[_offset_exposure], problems: list[str]
) -> tuple[list[_offset_exposure], list[_offset_exposure]]:
    """*split one template run into its ON and OFF frames, and record the run's classification problems*

    **Key Arguments:**

    - ``start`` -- the ``TPL START`` of the run
    - ``run`` -- the frames of the run, in exposure order
    - ``problems`` -- the list that a problem message is added to

    **Return:**

    - ``onFrames``, ``offFrames`` -- the run's ``_offset_exposure`` objects, in exposure order
    """
    if not _check_one_template(start, run, problems):
        return [], []
    problemCount = len(problems)
    classify = _classify_fixed_sky if run[0].technique == FIXED_SKY_TECHNIQUE else _classify_generic
    labels = classify(start, run, problems)
    onFrames = [exposure for exposure, label in zip(run, labels, strict=True) if label is frame_label.ON]
    offFrames = [exposure for exposure, label in zip(run, labels, strict=True) if label is frame_label.OFF]
    _check_balance(start, onFrames, offFrames, len(problems) > problemCount, problems)
    return onFrames, offFrames


def classify_offset_frames(
    log: Any, allObjectFrames: Sequence[Any], allFilenames: Sequence[str], kw: Callable[[str], str]
) -> tuple[list[Any], list[Any], list[str], list[str]]:
    """*split the science frames of an offset SOF into ON and OFF frames, one template run at a time*

    **Key Arguments:**

    - ``log`` -- logger
    - ``allObjectFrames`` -- the science frames
    - ``allFilenames`` -- the science frames' file names, in the same order
    - ``kw`` -- the keyword lookup

    **Return:**

    - ``allFrameON``, ``allFrameOFF``, ``allFrameONNames``, ``allFrameOFFNames`` -- the ON and OFF frames and their
      file names. Runs follow each other in the order they started and, within a run, frames follow their exposure
      numbers, so the nth ON frame and the nth OFF frame are always from the same run

    **Raises:**

    - ``OffsetClassificationError`` -- when a frame cannot be classified. The message names every such frame

    **Usage**

    ```python
    from soxspipe.commonutils.offset_frames import classify_offset_frames
    allFrameON, allFrameOFF, allFrameONNames, allFrameOFFNames = classify_offset_frames(
        log=log,
        allObjectFrames=frames,
        allFilenames=filenames,
        kw=kw,
    )
    ```
    """
    problems = []
    exposures = [
        _read_exposure(frame, filename, kw, problems)
        for frame, filename in zip(allObjectFrames, allFilenames, strict=True)
    ]
    runs = {}
    for exposure in filter(None, exposures):
        runs.setdefault(exposure.start, []).append(exposure)

    allOn, allOff = [], []
    for start in sorted(runs):
        run = sorted(runs[start], key=lambda exposure: exposure.expno)
        onFrames, offFrames = _split_run(start, run, problems)
        allOn.extend(onFrames)
        allOff.extend(offFrames)

    if problems:
        error = " ".join(problems)
        log.error(error)
        raise OffsetClassificationError(error)

    plural = "S" if len(runs) > 1 else ""
    log.print(f"# PROCESSING {len(allOn)} ON-OFF OFFSET CYCLES IN {len(runs)} TEMPLATE RUN{plural}")

    return (
        [exposure.frame for exposure in allOn],
        [exposure.frame for exposure in allOff],
        [exposure.filename for exposure in allOn],
        [exposure.filename for exposure in allOff],
    )
