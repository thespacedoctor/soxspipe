"""*the airmass and exposure time of a sequence of frames, for the flux calibration of the product they make (DY-1695)*

A stacked product collects light over every frame that contributes to it, so its extinction correction uses the mean
airmass of those frames, not the airmass at the end of the first one. Each frame contributes the midpoint of its
``ESO TEL AIRM START`` and ``ESO TEL AIRM END``, and frames are weighted by how long each exposure lasted.

The same frames must also share one exposure time, because flux calibration and the response curve divide the counts
of the product by it (DY-1714).
"""

import math

# AN AIRMASS IS AT LEAST 1 (THE ZENITH); LOWER VALUES ARE SENTINELS SUCH AS THE DATA ORGANISER'S -99.99
MINIMUM_AIRMASS = 1.0

# TWO EXPOSURE TIMES THIS CLOSE (IN SECONDS) ARE THE SAME EXPOSURE; THEY DIFFER BY FLOAT NOISE
EXPOSURE_TIME_TOLERANCE_SECONDS = 1e-3


def _finite_float(header, keyword):
    """*the value of a header card as a finite float, or None when it is missing or not a finite number*"""
    try:
        value = float(header[keyword])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _airmass_value(header, keyword):
    """*one airmass card as a float, or None when it is missing, not a finite number, or below 1*"""
    value = _finite_float(header, keyword)
    return value if value is not None and value >= MINIMUM_AIRMASS else None


def _frame_airmass(header, kw, log, frameNumber):
    """*the midpoint of one frame's start and end airmass, or the one of the two it has, or None*"""
    startKeyword, endKeyword = kw("AIRMASS_START"), kw("AIRMASS_END")
    start = _airmass_value(header, startKeyword)
    end = _airmass_value(header, endKeyword)
    if start is not None and end is not None:
        return (start + end) / 2
    if start is None and end is None:
        log.warning(f"Frame {frameNumber} has no readable {startKeyword} or {endKeyword}, so its airmass is not used")
        return None
    missing, present, value = (startKeyword, endKeyword, end) if start is None else (endKeyword, startKeyword, start)
    log.warning(f"Frame {frameNumber} has no readable {missing}, so its {present} ({value}) is used as its airmass")
    return value


def _exposure_duration(header, kw):
    """*a frame's exposure duration in seconds, or None if it cannot be read*

    ``ESO DET SEQ1 EXPTIME`` comes first: on SOXS NIR ``EXPTIME`` is one DIT, and the sequence of NDIT reads lasts
    much longer (DY-1666). VIS frames lack it, and their ``EXPTIME`` is the whole exposure. The X-Shooter keyword
    file has no alias for it.
    """
    try:
        duration = _finite_float(header, kw("SEQ_EXPTIME"))
    except LookupError:
        # THE X-SHOOTER KEYWORD FILE HAS NO SEQ_EXPTIME ALIAS, SO EXPTIME IS THE DURATION
        duration = None
    if duration is None:
        duration = _finite_float(header, kw("EXPTIME"))
    if duration is None or duration <= 0:
        return None
    return duration


def mean_airmass(headers, kw, log):
    """*the exposure-time-weighted mean airmass of the frames that contribute to a product*

    A frame without a readable start airmass uses its end airmass, and the reverse, with a warning. A frame with
    neither is left out. When any frame used lacks a readable positive exposure duration, every frame gets equal
    weight.

    **Key Arguments:**

    - ``headers`` -- the FITS headers of the raw frames that contribute to the product
    - ``kw`` -- the keyword lookup, mapping an alias such as ``AIRMASS_END`` to the header keyword
    - ``log`` -- logger

    **Return:**

    - ``airmass`` -- the mean airmass, or *None* when no frame has a readable airmass

    **Usage:**

    ```python
    from soxspipe.commonutils.sequence_airmass import mean_airmass
    airmass = mean_airmass([frame.header for frame in allObjectFrames], self.kw, self.log)
    ```
    """
    framesUsed = []
    for frameNumber, header in enumerate(headers, start=1):
        airmass = _frame_airmass(header, kw, log, frameNumber)
        if airmass is not None:
            framesUsed.append((airmass, _exposure_duration(header, kw)))
    if not framesUsed:
        log.warning("No contributing frame has a readable airmass, so the product has no airmass")
        return None

    weights = [duration for _, duration in framesUsed]
    if any(weight is None for weight in weights):
        log.warning(
            "A contributing frame has no readable exposure time, so every frame gets equal weight in the airmass"
        )
        weights = [1.0] * len(framesUsed)
    return sum(airmass * weight for (airmass, _), weight in zip(framesUsed, weights, strict=True)) / sum(weights)


def _exposure_time(header, keyword):
    """*one frame's exposure time as a float, or None when it is missing, not a finite number, or not positive*"""
    value = _finite_float(header, keyword)
    return value if value is not None and value > 0 else None


def _frame_label(header, frameNumber):
    """*the raw file name of a frame, or its position in the sequence when its header has none*"""
    for keyword in ("ARCFILE", "ORIGFILE"):
        name = str(header.get(keyword, "")).strip()
        if name:
            return name
    return f"frame {frameNumber}"


def require_equal_exposure_times(headers, kw):
    """*raise an error unless every frame that contributes to a product has the same exposure time*

    Stacked frames are mean-combined and nod and offset extractions median-combined, so the counts of the product
    are per single exposure. Flux calibration and the response curve divide them by one frame's ``EXPTIME``, which
    is right only when every exposure lasted that long. ``EXPTIME`` is the divisor, so it is the value compared:
    on SOXS NIR it is one DIT.

    **Key Arguments:**

    - ``headers`` -- the FITS headers of the frames that contribute to the product
    - ``kw`` -- the keyword lookup, mapping an alias such as ``EXPTIME`` to the header keyword

    **Return:**

    - None

    **Usage:**

    ```python
    from soxspipe.commonutils.sequence_airmass import require_equal_exposure_times
    require_equal_exposure_times([frame.header for frame in allObjectFrames], self.kw)
    ```
    """
    exptimes = [_exposure_time(header, kw("EXPTIME")) for header in headers]
    if not exptimes:
        return
    labels = [_frame_label(header, frameNumber) for frameNumber, header in enumerate(headers, start=1)]
    unreadable = [label for label, exptime in zip(labels, exptimes, strict=True) if exptime is None]
    if unreadable:
        raise ValueError(
            f"Flux calibration needs the exposure time of every frame, but there is no readable {kw('EXPTIME')} "
            f"in: {', '.join(unreadable)}"
        )
    # COMPARE THE WHOLE SPREAD, SO FRAMES EACH CLOSE TO THE FIRST CANNOT DIFFER FROM EACH OTHER BY MORE
    if max(exptimes) - min(exptimes) <= EXPOSURE_TIME_TOLERANCE_SECONDS:
        return
    frames = ", ".join(f"{label} ({exptime:g} s)" for label, exptime in zip(labels, exptimes, strict=True))
    raise ValueError(
        "Flux calibration divides the counts by one exposure time, so every frame of the sequence must have the "
        f"same exposure time, but they differ: {frames}"
    )
