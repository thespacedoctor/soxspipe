"""Resolve the read noise a reduction uses (DY-1458, DY-1581, DY-2341).

The master bias (UVB/VIS) or master dark (X-Shooter NIR) records ``ESO QC RAW RON``, the single-frame read noise
measured from a pair of raw frames, and ``ESO PRO DATANCOM``, the number of frames it combined. The measured value
replaces the nominal header value. A master frame that lacks either is an old product, and the resolver falls back to
the nominal value in the raw-frame header, then to the detector default.

A SOXS NIR frame takes its read noise from its own up-the-ramp sample count instead (``ESO DET NDSAMPLES``), looked up
in a table of the detector parameters. The master dark's ``ESO QC RAW RON`` is a QC value there and sets nothing.
"""

import math
from dataclasses import dataclass

# THE HEADER CARD OF THE READ NOISE MEASURED FROM A RAW-FRAME PAIR
MEASURED_RON_CARD = "ESO QC RAW RON"
# THE CARD OF THE FIRST RAW FRAME A PRODUCT WAS BUILT FROM; THE OTHERS COUNT UP FROM 1
RAW_FRAME_CARD_TEMPLATE = "ESO PRO REC1 RAW{index} NAME"
# THE WORKSPACE DIRECTORY HOLDING THE RAW FRAMES, IN ONE DATED FOLDER EACH
RAW_DIRECTORY_NAME = "raw"
# THE MOST RAW FRAMES OF ONE NAME THAT A MASTER DARK'S SAMPLE COUNT IS LOOKED FOR IN
MAX_RAW_CANDIDATES = 10
# WHAT READING A CORRUPT FITS FILE CAN RAISE
# IMPORTERROR: ASTROPY NEEDS AN OPTIONAL PACKAGE TO DECOMPRESS A .Z FILE
FITS_READ_ERRORS = (OSError, ValueError, EOFError, ImportError)
# A MEASURED VALUE THIS FAR (AS A FRACTION) FROM THE NOMINAL ONE IS WORTH A WARNING
NOMINAL_DISAGREEMENT_FRACTION = 0.15
# A NIR MASTER DARK QC RAW RON THIS FAR (AS A FRACTION) FROM THE TABLE VALUE AT ITS SAMPLE COUNT IS WORTH A WARNING
# (DY-2341)
DARK_RON_TABLE_TOLERANCE = 0.20

SOURCE_MASTER_BIAS = "master bias"
SOURCE_MASTER_DARK = "master dark"
SOURCE_RAW_HEADER = "raw-frame header"
SOURCE_DETECTOR_DEFAULT = "detector default"
SOURCE_NDSAMPLES_TABLE = "NDSAMPLES table"


@dataclass(frozen=True)
class ReadNoise:
    """*the read noise of a reduction and where it came from*

    **Key Arguments:**

    - ``singleFrame`` -- the read noise (e-) of one frame
    - ``nCombined`` -- the number of frames the master frame combined, or *None* when there is no master header
    - ``skyFit`` -- the read noise (e-) of a frame with the master frame subtracted: ``singleFrame * sqrt(1 + 1/N)``,
      or ``singleFrame`` when there is no master header. For a SOXS NIR frame it is
      ``sqrt(singleFrame**2 + masterDarkRon**2 / N)``, the master dark's noise looked up at its own sample count
    - ``source`` -- where ``singleFrame`` came from: the master bias, the master dark, the raw-frame header, the NIR
      sample-count table or the detector default
    """

    singleFrame: float
    nCombined: int | None
    skyFit: float
    source: str


def resolve_read_noise(masterHeader, nominalRon, defaultRon, kw, log, masterLabel=SOURCE_MASTER_BIAS):
    """*resolve the read noise (e-) of a reduction from the header of its master bias (UVB/VIS) or master dark (NIR)*

    **Key Arguments:**

    - ``masterHeader`` -- the primary header of the master frame, or *None* when it cannot be read
    - ``nominalRon`` -- the read noise (e-) in the raw-frame headers, or *None* when they carry none. When *None*, a
      measured value is compared with ``defaultRon``
    - ``defaultRon`` -- the detector-parameter read noise (e-)
    - ``kw`` -- the keyword lookup (``keyword_lookup(...).get``)
    - ``log`` -- logger
    - ``masterLabel`` -- what the master frame is called in warnings and in ``ReadNoise.source``: ``"master bias"``
      (default) or ``"master dark"``

    **Return:**

    - ``readNoise`` -- a ``ReadNoise``. The measured value is used when the master frame holds a finite, positive
      ``ESO QC RAW RON`` and an ``ESO PRO DATANCOM``; otherwise the nominal value, then the default, with a warning.
      A measured value more than 15% from the nominal one also logs a warning but is still used

    **Usage:**

    ```python
    from soxspipe.commonutils.read_noise import resolve_read_noise
    readNoise = resolve_read_noise(masterBiasHeader, 3.3, 3.8, kw, log)
    nirReadNoise = resolve_read_noise(masterDarkHeader, None, 7.3, kw, log, masterLabel="master dark")
    ```
    """
    nominal = _positive_number(nominalRon)
    referenceRon = nominal if nominal is not None else float(defaultRon)

    if masterHeader is None:
        singleFrame, source = _fallback(nominal, defaultRon)
        log.warning(
            f"The {masterLabel} header could not be read; the read noise is the {source} value of {singleFrame} e-."
        )
        return ReadNoise(singleFrame, None, singleFrame, source)

    measured = _positive_number(masterHeader.get(MEASURED_RON_CARD))
    datancom = _positive_integer(masterHeader.get(kw("DATANCOM")))
    nCombined = datancom if datancom is not None else _count_raw_frames(masterHeader)

    if measured is not None and datancom is not None:
        if abs(measured - referenceRon) / referenceRon > NOMINAL_DISAGREEMENT_FRACTION:
            log.info(
                f"The {masterLabel} read noise of {measured:0.2f} e- is more than "
                f"{NOMINAL_DISAGREEMENT_FRACTION:.0%} from the nominal {referenceRon:0.2f} e-; "
                "using the measured value."
            )
        return ReadNoise(measured, datancom, _sky_fit_ron(measured, datancom), masterLabel)

    if nCombined is None:
        nCombined = 1
        log.warning(f"The {masterLabel} does not say how many frames it combined; counting it as 1.")

    singleFrame, source = _fallback(nominal, defaultRon)
    reason = _unusable_reason(measured, datancom, kw)
    log.info(
        f"The {masterLabel} has no usable measured read noise ({reason}); using the {source} value of {singleFrame} e-."
    )
    return ReadNoise(singleFrame, nCombined, _sky_fit_ron(singleFrame, nCombined), source)


def _fallback(nominal, defaultRon):
    """*the nominal read noise and its source, else the detector default and its source*"""
    if nominal is not None:
        return nominal, SOURCE_RAW_HEADER
    return float(defaultRon), SOURCE_DETECTOR_DEFAULT


def _unusable_reason(measured, datancom, kw):
    """*say why a master frame read noise could not be used*"""
    if measured is None:
        return f"{MEASURED_RON_CARD} is missing or not a positive number"
    return f"it lacks {kw('DATANCOM')}, so it is an old product"


def _sky_fit_ron(singleFrame, nCombined):
    """*the read noise (e-) of a frame with a master frame of ``nCombined`` frames subtracted*"""
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


def table_read_noise(ndsamples, table):
    """*the NIR read noise (e-) at a sample count, interpolated linearly in log(NDSAMPLES) between the table points*

    **Key Arguments:**

    - ``ndsamples`` -- the ``ESO DET NDSAMPLES`` value of a frame
    - ``table`` -- the measured read noise (e-) keyed by sample count, from the detector parameters

    **Return:**

    - ``ron`` -- the read noise (e-), or *None* when ``ndsamples`` is not a whole number inside the range of the table

    **Usage:**

    ```python
    from soxspipe.commonutils.read_noise import table_read_noise
    ron = table_read_noise(97, detectorParams["read-noise-by-ndsamples"])
    ```
    """
    import numpy as np

    if _ndsamples_problem(ndsamples, table) is not None:
        return None
    counts = sorted(table)
    return float(np.interp(math.log(ndsamples), np.log(counts), [table[count] for count in counts]))


def read_noise_from_ndsamples(ndsamples, table, fallbackRon, frameLabel, log):
    """*the NIR read noise (e-) of a frame from its sample count, else the fallback with a warning*

    **Key Arguments:**

    - ``ndsamples`` -- the ``ESO DET NDSAMPLES`` value of the frame, or *None* when it has no such card
    - ``table`` -- the measured read noise (e-) keyed by sample count
    - ``fallbackRon`` -- the read noise (e-) for a frame whose count the table cannot serve
    - ``frameLabel`` -- the frame name, for the warning
    - ``log`` -- logger

    **Return:**

    - ``ron`` -- the table value from ``table_read_noise``, or ``fallbackRon`` with a warning that names the frame and
      the reason

    **Usage:**

    ```python
    from soxspipe.commonutils.read_noise import read_noise_from_ndsamples
    ron = read_noise_from_ndsamples(10, table, 7.55, "SOXS.2026-01-31T19:16:38.123.fits", log)
    ```
    """
    ron = table_read_noise(ndsamples, table)
    if ron is not None:
        return ron

    log.warning(
        f"The frame {frameLabel} has {_ndsamples_problem(ndsamples, table)}; "
        f"the read noise is the detector value of {float(fallbackRon)} e-."
    )
    return float(fallbackRon)


def resolve_nir_read_noise(rawFrames, masterDarkPath, table, fallbackRon, workspaceRoot, kw, log):
    """*resolve the read noise of a SOXS NIR reduction from the sample counts of its frames*

    The single-frame read noise is the table value at the sample count of the raw frames. The sky-fit noise adds
    the master dark's own noise when a master dark is subtracted: ``sqrt(ron_frame**2 + ron_dark**2 / N)``, where
    ``ron_dark`` is the table value at the dark's sample count and ``N`` the number of frames the dark combined.

    **Key Arguments:**

    - ``rawFrames`` -- a list of ``(name, ndsamples)`` pairs, one for each raw frame; ``ndsamples`` is *None* for a
      frame with no ``ESO DET NDSAMPLES`` card
    - ``masterDarkPath`` -- the path to the master dark that is subtracted, or *None* when none is
    - ``table`` -- the measured read noise (e-) keyed by sample count, from the NIR detector parameters
    - ``fallbackRon`` -- the detector-parameter read noise (e-), for a count the table cannot serve
    - ``workspaceRoot`` -- the workspace root, where the dated folders of raw frames are looked for when a master dark
      lacks the sample-count card
    - ``kw`` -- the keyword lookup (``keyword_lookup(...).get``)
    - ``log`` -- logger

    **Return:**

    - ``readNoise`` -- a ``ReadNoise``, or *None* when there are no raw frames. Raw frames with different sample
      counts give the largest read noise and a warning. A frame or dark whose count the table cannot serve uses
      ``fallbackRon`` and a warning that names it

    **Usage:**

    ```python
    from soxspipe.commonutils.read_noise import resolve_nir_read_noise
    readNoise = resolve_nir_read_noise([("raw_1.fits", 128)], "master_dark.fits", table, 7.55, workspaceRoot, kw, log)
    ```
    """
    import os

    from astropy.io import fits

    if not rawFrames:
        return None

    frameRon = _raw_frames_read_noise(rawFrames, table, fallbackRon, log)
    isFromTable = all(table_read_noise(ndsamples, table) is not None for _, ndsamples in rawFrames)
    source = SOURCE_NDSAMPLES_TABLE if isFromTable else SOURCE_DETECTOR_DEFAULT

    if masterDarkPath is None:
        return ReadNoise(frameRon, None, frameRon, source)

    try:
        masterHeader = fits.getheader(masterDarkPath)
    except FITS_READ_ERRORS as e:
        log.warning(f"The master dark header could not be read ({e}); the sky fit adds no master dark noise.")
        return ReadNoise(frameRon, None, frameRon, source)

    darkRon = read_noise_from_ndsamples(
        _master_ndsamples(masterHeader, workspaceRoot, kw, log),
        table,
        fallbackRon,
        os.path.basename(masterDarkPath),
        log,
    )
    nCombined = _combined_frame_count(masterHeader, kw, log)
    skyFit = math.sqrt(frameRon**2 + darkRon**2 / nCombined)
    return ReadNoise(frameRon, nCombined, skyFit, source)


def _raw_frames_read_noise(rawFrames, table, fallbackRon, log):
    """*the largest read noise (e-) of the raw frames, each looked up from its own sample count*"""
    frameRons = [read_noise_from_ndsamples(count, table, fallbackRon, name, log) for name, count in rawFrames]
    usableCounts = sorted({int(count) for _, count in rawFrames if table_read_noise(count, table) is not None})
    if len(usableCounts) > 1:
        log.warning(
            f"The raw frames have different sample counts ({usableCounts}); "
            f"the read noise is the largest, {max(frameRons):0.2f} e-."
        )
    return max(frameRons)


def _master_ndsamples(masterHeader, workspaceRoot, kw, log):
    """*the sample count of a master dark: its own card, else that of a raw frame it lists, else None*

    The raw frame is the first one the master dark lists. Its name comes from file content, so it is used only when
    it is a plain file name, and only files that resolve inside ``<workspaceRoot>/raw`` are read.
    """
    card = kw("DET_NDSAMPLES")
    ndsamples = masterHeader.get(card)
    if ndsamples is not None:
        return ndsamples

    rawName = masterHeader.get(RAW_FRAME_CARD_TEMPLATE.format(index=1))
    if not rawName:
        return None
    if not _is_plain_file_name(rawName):
        log.warning(
            f"The master dark lists the raw frame {rawName!r}, which is not a plain file name; "
            "its sample count is not looked up there."
        )
        return None

    for rawPath in _raw_frame_candidates(rawName, workspaceRoot):
        ndsamples = _raw_frame_ndsamples(rawPath, rawName, card, log)
        if ndsamples is not None:
            return ndsamples
    return None


def _is_plain_file_name(name):
    """*True for a string that is one file name: no directory part, backslash or NUL, and not ``.`` or ``..``*"""
    import os

    return (
        isinstance(name, str)
        and name not in ("", ".", "..")
        and not any(character in name for character in ("/", "\\", "\0"))
        and name == os.path.basename(name)
        and not os.path.isabs(name)
    )


def _raw_frame_candidates(rawName, workspaceRoot):
    """*the files named ``rawName`` in a dated folder of ``<workspaceRoot>/raw``, at most ``MAX_RAW_CANDIDATES``*

    An uncompressed file comes before a ``.Z`` one, then by path. A file that resolves outside the raw directory (for
    example through a symbolic link) is left out.
    """
    import glob
    import os

    rawRoot = os.path.realpath(os.path.join(str(workspaceRoot), RAW_DIRECTORY_NAME))
    candidates = []
    for name in (rawName, rawName + ".Z"):
        for path in sorted(glob.glob(os.path.join(glob.escape(rawRoot), "*", glob.escape(name)))):
            resolved = os.path.realpath(path)
            if _is_inside(resolved, rawRoot) and os.path.isfile(resolved):
                # THE CHECKED PATH IS THE ONE READ, SO A LINK SWAPPED AFTER THE CHECK IS NOT FOLLOWED
                candidates.append(resolved)
    return candidates[:MAX_RAW_CANDIDATES]


def _is_inside(path, root):
    """*True when ``path`` lies inside ``root``; paths on different drives are not inside*"""
    import os

    try:
        return os.path.commonpath([root, path]) == root
    except ValueError:
        return False


def _raw_frame_ndsamples(rawPath, rawName, card, log):
    """*the ``card`` value in the header of a raw frame, or None when it has none or cannot be read*"""
    from astropy.io import fits

    try:
        return fits.getheader(rawPath).get(card)
    except FITS_READ_ERRORS as e:
        log.warning(f"The raw frame {rawName} listed by the master dark could not be read ({e}).")
        return None


def _combined_frame_count(masterHeader, kw, log):
    """*the number of frames a master frame combined, counting a master that does not say as 1 with a warning*"""
    nCombined = _positive_integer(masterHeader.get(kw("DATANCOM"))) or _count_raw_frames(masterHeader)
    if nCombined is None:
        log.warning("The master dark does not say how many frames it combined; counting it as 1.")
        return 1
    return nCombined


def _ndsamples_problem(ndsamples, table):
    """*say why a sample count cannot be looked up in the table, or None when it can*"""
    import numpy as np

    if ndsamples is None or ndsamples is np.ma.masked:
        return "no ESO DET NDSAMPLES card"
    isNumber = isinstance(ndsamples, (int, float, np.integer, np.floating)) and not isinstance(ndsamples, bool)
    if not isNumber or not math.isfinite(ndsamples) or ndsamples != int(ndsamples):
        return f"an ESO DET NDSAMPLES value ({ndsamples!r}) that is not a whole number"
    lowest, highest = min(table), max(table)
    if not lowest <= ndsamples <= highest:
        return f"an ESO DET NDSAMPLES of {int(ndsamples)}, outside the measured range of {lowest} to {highest}"
    return None
