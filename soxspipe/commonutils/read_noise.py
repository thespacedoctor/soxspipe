"""Resolve the read noise a UVB/VIS reduction uses from the master bias its SOF names (DY-1458).

The master bias records ``ESO QC RAW RON``, the single-frame read noise measured from a pair of raw frames, and
``ESO PRO DATANCOM``, the number of frames it combined. The measured value replaces the nominal header value. A master
bias that lacks either is an old product, and the resolver falls back to the nominal value in the raw-frame header,
then to the detector default.
"""

import math
from dataclasses import dataclass

# THE HEADER CARD OF THE READ NOISE MEASURED FROM A RAW-FRAME PAIR
MEASURED_RON_CARD = "ESO QC RAW RON"
# THE CARD OF THE FIRST RAW FRAME A PRODUCT WAS BUILT FROM; THE OTHERS COUNT UP FROM 1
RAW_FRAME_CARD_TEMPLATE = "ESO PRO REC1 RAW{index} NAME"
# A MEASURED VALUE THIS FAR (AS A FRACTION) FROM THE NOMINAL ONE IS WORTH A WARNING
NOMINAL_DISAGREEMENT_FRACTION = 0.15

SOURCE_MASTER_BIAS = "master bias"
SOURCE_RAW_HEADER = "raw-frame header"
SOURCE_DETECTOR_DEFAULT = "detector default"


@dataclass(frozen=True)
class ReadNoise:
    """*the read noise of a reduction and where it came from*

    **Key Arguments:**

    - ``singleFrame`` -- the read noise (e-) of one frame
    - ``nCombined`` -- the number of frames the master bias combined, or *None* when there is no master bias header
    - ``skyFit`` -- the read noise (e-) of a frame with the master bias subtracted: ``singleFrame * sqrt(1 + 1/N)``,
      or ``singleFrame`` when there is no master bias header
    - ``source`` -- where ``singleFrame`` came from: the master bias, the raw-frame header or the detector default
    """

    singleFrame: float
    nCombined: int | None
    skyFit: float
    source: str


def resolve_read_noise(masterBiasHeader, nominalRon, defaultRon, kw, log):
    """*resolve the read noise (e-) of a UVB/VIS reduction from its master bias header*

    **Key Arguments:**

    - ``masterBiasHeader`` -- the primary header of the master bias, or *None* when it cannot be read
    - ``nominalRon`` -- the read noise (e-) in the raw-frame headers, or *None* when they carry none
    - ``defaultRon`` -- the detector-parameter read noise (e-)
    - ``kw`` -- the keyword lookup (``keyword_lookup(...).get``)
    - ``log`` -- logger

    **Return:**

    - ``readNoise`` -- a ``ReadNoise``. The measured value is used when the master bias holds a finite, positive
      ``ESO QC RAW RON`` and an ``ESO PRO DATANCOM``; otherwise the nominal value, then the default, with a warning.
      A measured value more than 15% from the nominal one also logs a warning but is still used

    **Usage:**

    ```python
    from soxspipe.commonutils.read_noise import resolve_read_noise
    readNoise = resolve_read_noise(masterBiasHeader, 3.3, 3.8, kw, log)
    ```
    """
    nominal = _positive_number(nominalRon)
    referenceRon = nominal if nominal is not None else float(defaultRon)

    if masterBiasHeader is None:
        singleFrame, source = _fallback(nominal, defaultRon)
        log.warning(
            f"The master bias header could not be read; the read noise is the {source} value of {singleFrame} e-."
        )
        return ReadNoise(singleFrame, None, singleFrame, source)

    measured = _positive_number(masterBiasHeader.get(MEASURED_RON_CARD))
    datancom = _positive_integer(masterBiasHeader.get(kw("DATANCOM")))
    nCombined = datancom if datancom is not None else _count_raw_frames(masterBiasHeader)

    if measured is not None and datancom is not None:
        if abs(measured - referenceRon) / referenceRon > NOMINAL_DISAGREEMENT_FRACTION:
            log.warning(
                f"The master bias read noise of {measured:0.2f} e- is more than "
                f"{NOMINAL_DISAGREEMENT_FRACTION:.0%} from the nominal {referenceRon:0.2f} e-; "
                "using the measured value."
            )
        return ReadNoise(measured, datancom, _sky_fit_ron(measured, datancom), SOURCE_MASTER_BIAS)

    if nCombined is None:
        nCombined = 1
        log.warning("The master bias does not say how many frames it combined; counting it as 1.")

    singleFrame, source = _fallback(nominal, defaultRon)
    reason = _unusable_reason(measured, datancom, kw)
    log.warning(
        f"The master bias has no usable measured read noise ({reason}); using the {source} value of {singleFrame} e-."
    )
    return ReadNoise(singleFrame, nCombined, _sky_fit_ron(singleFrame, nCombined), source)


def _fallback(nominal, defaultRon):
    """*the nominal read noise and its source, else the detector default and its source*"""
    if nominal is not None:
        return nominal, SOURCE_RAW_HEADER
    return float(defaultRon), SOURCE_DETECTOR_DEFAULT


def _unusable_reason(measured, datancom, kw):
    """*say why a master bias read noise could not be used*"""
    if measured is None:
        return f"{MEASURED_RON_CARD} is missing or not a positive number"
    return f"it lacks {kw('DATANCOM')}, so it is an old product"


def _sky_fit_ron(singleFrame, nCombined):
    """*the read noise (e-) of a frame with a master bias of ``nCombined`` frames subtracted*"""
    return singleFrame * math.sqrt(1.0 + 1.0 / nCombined)


def _positive_number(value):
    """*the value as a finite, positive float, else None; a boolean is not a number here*"""
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0 else None


def _positive_integer(value):
    """*the value as a positive integer, else None*"""
    number = _positive_number(value)
    if number is None or number != int(number):
        return None
    return int(number)


def _count_raw_frames(header):
    """*the number of consecutive ``ESO PRO REC1 RAW<i> NAME`` cards from 1, or None when there are none*"""
    count = 0
    while RAW_FRAME_CARD_TEMPLATE.format(index=count + 1) in header:
        count += 1
    return count or None
