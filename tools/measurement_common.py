"""*shared helpers of the developer measurement scripts*

The scripts `measure_log_wavelength_step.py`, `measure_spec_res.py` and `measure_sky_shift_uncertainty.py` measure
instrument values that ship in `soxspipe/resources/soxs_detector_parameters.yaml`. They are developer tools, not part
of a reduction and not in the package CLI. This module holds what the three share: path expansion, the arm, the
spectral binning and the other required cards of a product, the yaml block layout with its provenance comment, and the
exit-status convention of the other tools (0 done, 2 the input could not be used).
"""

from __future__ import annotations

import functools
import glob
import math
import re
import sys
from collections.abc import Callable, Iterable, Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from types import MappingProxyType
from typing import Any

import yaml
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[1]
DETECTOR_PARAMETERS_PATH = REPO_ROOT / "soxspipe" / "resources" / "soxs_detector_parameters.yaml"
KEYWORDS_PATH = REPO_ROOT / "soxspipe" / "resources" / "soxs_keywords.yaml"
VERSION_PATH = REPO_ROOT / "soxspipe" / "__version__.py"
ARM_CARD = "ESO SEQ ARM"
MJD_CARD = "MJD-OBS"
# THE VALUE OF dispersion-axis IN soxs_detector_parameters.yaml ("x" OR "y") -> THE NUMPY AXIS THE WAVELENGTH RUNS ALONG
# IN THE ARRAYS OF A DISP_IMAGE. HERE "x" MEANS NUMPY AXIS 0, WHICH IS FITS NAXIS2, THE OPPOSITE OF THE FITS x AXIS
# (NAXIS1, NUMPY AXIS 1). VERIFIED ON THE 2026-01-12 MAPS: THE WAVELENGTH STEP IS 0.0335 nm ALONG AXIS 0 FOR VIS AND
# 0.0498 nm ALONG AXIS 1 FOR NIR (DY-2455)
DISP_IMAGE_AXIS_BY_DISPERSION_AXIS_VALUE = {"x": 0, "y": 1}
# THE ESO DET BIN CARD THAT BINS EACH NUMPY AXIS OF A RAW FRAME: BINX BINS NAXIS1 (AXIS 1), BINY BINS NAXIS2 (AXIS 0).
# VIS 1X1 RAW FRAMES ARE 858 BY 4096 (NAXIS1 BY NAXIS2) AND 2X1 FRAMES ARE 429 BY 4096 (BINX 2, BINY 1)
ROW_AXIS = 0
COLUMN_AXIS = 1
BINNING_CARD_BY_NUMPY_AXIS = {ROW_AXIS: "ESO DET BINY", COLUMN_AXIS: "ESO DET BINX"}
# THE NIR DETECTOR DOES NOT BIN, SO ITS RAW FRAMES AND PRODUCTS CARRY NO BINNING CARD (DY-2456)
UNBINNED_ARMS = frozenset({"NIR"})
EXIT_INPUT_ERROR = 2
UNKNOWN_VERSION = "unknown"
YAML_INDENT = "  "


class ToolError(Exception):
    """*an input problem the measurement scripts report and stop on*"""


def load_yaml(path: Path) -> dict[str, Any]:
    """*read a yaml file*

    **Key Arguments:**

    - ``path`` -- the yaml file to read

    **Return:**

    - ``content`` -- the parsed mapping

    **Usage:**

    ```python
    parameters = load_yaml(DETECTOR_PARAMETERS_PATH)
    ```
    """
    with open(path, encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def freeze(value: Any) -> Any:
    """*a read-only copy of a parsed yaml value: mappings become ``MappingProxyType``, lists become tuples*

    **Key Arguments:**

    - ``value`` -- a mapping, list or scalar

    **Return:**

    - ``frozen`` -- the read-only copy, frozen all the way down

    **Usage:**

    ```python
    parameters = freeze({"VIS": {"binning": [1, 1]}})
    ```
    """
    if isinstance(value, Mapping):
        return MappingProxyType({key: freeze(item) for key, item in value.items()})
    if isinstance(value, list | tuple):
        return tuple(freeze(item) for item in value)
    return value


def frozen_mapping[Key, Value](mapping: Mapping[Key, Value]) -> Mapping[Key, Value]:
    """*a read-only copy of a mapping, for the fields of a frozen dataclass*

    **Key Arguments:**

    - ``mapping`` -- the mapping to copy

    **Return:**

    - ``copy`` -- a ``MappingProxyType`` over a copy of the mapping

    **Usage:**

    ```python
    perOrder = frozen_mapping({1: 0.3})
    ```
    """
    return MappingProxyType(dict(mapping))


@functools.cache
def detector_parameters() -> Mapping[str, Any]:
    """*the packaged detector parameters, read once and read-only*

    **Return:**

    - ``parameters`` -- the parsed ``soxs_detector_parameters.yaml``

    **Usage:**

    ```python
    scienceRows = detector_parameters()["VIS"]["science-pixels"]["rows"]
    ```
    """
    return freeze(load_yaml(DETECTOR_PARAMETERS_PATH))


@functools.cache
def keyword_aliases() -> Mapping[str, str]:
    """*the packaged SOXS keyword aliases, read once and read-only*

    **Return:**

    - ``aliases`` -- the parsed ``soxs_keywords.yaml``, alias to header keyword

    **Usage:**

    ```python
    keyword = keyword_aliases()["MJDOBS"]
    ```
    """
    return freeze(load_yaml(KEYWORDS_PATH))


def expand_paths(patterns: Iterable[str]) -> list[Path]:
    """*expand glob patterns, with a leading ``~``, into one sorted list of unique paths*

    **Key Arguments:**

    - ``patterns`` -- file paths or glob patterns

    **Return:**

    - ``paths`` -- the matching paths, sorted and without duplicates

    **Raises:**

    - ``ToolError`` -- if a pattern matches no file

    **Usage:**

    ```python
    paths = expand_paths(["~/data/*_EXTRACTED_ORDERS.fits"])
    ```
    """
    paths: set[Path] = set()
    for pattern in patterns:
        matches = glob.glob(str(Path(pattern).expanduser()))
        if not matches:
            raise ToolError(f"no file matches {pattern!r}")
        paths.update(Path(match) for match in matches)
    return sorted(paths)


def read_header(path: Path) -> fits.Header:
    """*the primary header of a FITS file*

    **Key Arguments:**

    - ``path`` -- the FITS file

    **Return:**

    - ``header`` -- the primary header

    **Raises:**

    - ``ToolError`` -- if the file cannot be read as FITS

    **Usage:**

    ```python
    header = read_header(Path("frame.fits"))
    ```
    """
    try:
        return fits.getheader(path)
    except OSError as error:
        raise ToolError(f"cannot read {path}: {error}") from error


def _where(path: Path | None) -> str:
    """*how an error names the source of a header: the file, or ``the header`` when there is no file*"""
    return f"{path}" if path is not None else "the header"


def required_card(header: fits.Header, card: str, path: Path | None = None) -> Any:
    """*the value of a header card that the measurement needs*

    **Key Arguments:**

    - ``header`` -- the FITS header
    - ``card`` -- the card name, for example ``ESO DET BINX``
    - ``path`` -- the file the header came from, named in the error. Default *None*

    **Return:**

    - ``value`` -- the card value

    **Raises:**

    - ``ToolError`` -- if the card is missing, naming the file and the card

    **Usage:**

    ```python
    slit = required_card(header, "ESO INS VISE NAME", path)
    ```
    """
    if card not in header:
        raise ToolError(f"{_where(path)} has no {card} card")
    return header[card]


def read_mjd_obs(header: fits.Header, path: Path | None = None) -> float:
    """*the observation MJD of a frame, from its ``MJD-OBS`` card*

    **Key Arguments:**

    - ``header`` -- the FITS header
    - ``path`` -- the file the header came from, named in the error. Default *None*

    **Return:**

    - ``mjd`` -- the observation MJD

    **Raises:**

    - ``ToolError`` -- if the card is missing or not a number

    **Usage:**

    ```python
    mjd = read_mjd_obs(header, path)
    ```
    """
    value = required_card(header, MJD_CARD, path)
    try:
        return float(value)
    except (TypeError, ValueError) as error:
        raise ToolError(f"{_where(path)} has a {MJD_CARD} card that is not a number: {value!r}") from error


def arm_of(header: fits.Header, path: Path | None = None) -> str:
    """*the arm of a product, ``VIS`` or ``NIR``, from the sequencer arm card*

    **Key Arguments:**

    - ``header`` -- the FITS header of the product
    - ``path`` -- the file the header came from, named in the error. Default *None*

    **Return:**

    - ``arm`` -- the arm, as it is keyed in ``soxs_detector_parameters.yaml``

    **Raises:**

    - ``ToolError`` -- if the card is missing or the arm has no detector parameters

    **Usage:**

    ```python
    arm = arm_of(header, path)
    ```
    """
    arm = str(required_card(header, ARM_CARD, path)).strip().upper()
    if arm not in detector_parameters():
        raise ToolError(f"{ARM_CARD} is {arm!r}, which has no entry in soxs_detector_parameters.yaml")
    return arm


def dispersion_axis(arm: str) -> str:
    """*the ``dispersion-axis`` value of the arm in the detector parameters, ``x`` or ``y``*

    This is not a numpy axis, and for the arrays of a ``DISP_IMAGE`` it is the opposite of the FITS axis name. Use
    ``disp_image_axis`` for the numpy axis.

    **Key Arguments:**

    - ``arm`` -- ``VIS`` or ``NIR``

    **Return:**

    - ``axis`` -- ``x`` or ``y``

    **Usage:**

    ```python
    axis = dispersion_axis("NIR")
    ```
    """
    return detector_parameters()[arm]["dispersion-axis"]


def disp_image_axis(arm: str) -> int:
    """*the numpy axis that the wavelength runs along in the ``DISP_IMAGE`` arrays and the raw frame of an arm*

    The detector parameter ``dispersion-axis`` of VIS is ``x``, and that is numpy axis 0 (FITS ``NAXIS2``). The one of
    NIR is ``y``, and that is numpy axis 1 (FITS ``NAXIS1``).

    **Key Arguments:**

    - ``arm`` -- ``VIS`` or ``NIR``

    **Return:**

    - ``axis`` -- 0 or 1

    **Usage:**

    ```python
    axis = disp_image_axis("VIS")
    ```
    """
    return DISP_IMAGE_AXIS_BY_DISPERSION_AXIS_VALUE[dispersion_axis(arm)]


def binning_along(header: fits.Header, numpyAxis: int, arm: str, path: Path | None = None) -> int:
    """*the binning factor along a numpy axis of a raw frame, from its ``ESO DET BINX`` or ``BINY`` card*

    ``BINX`` bins numpy axis 1 (``NAXIS1``) and ``BINY`` bins numpy axis 0 (``NAXIS2``).

    **Key Arguments:**

    - ``header`` -- the FITS header of the frame
    - ``numpyAxis`` -- 0 or 1
    - ``arm`` -- the arm of the frame. A card that is missing on an arm in ``UNBINNED_ARMS`` (NIR) counts as 1
    - ``path`` -- the file the header came from, named in the error. Default *None*

    **Return:**

    - ``binning`` -- the number of detector pixels in one binned pixel along the axis

    **Raises:**

    - ``ToolError`` -- if the card is missing (except on an unbinned arm), not a whole number, or below 1

    **Usage:**

    ```python
    binColumns = binning_along(header, COLUMN_AXIS, "VIS", path)
    ```
    """
    card = BINNING_CARD_BY_NUMPY_AXIS[numpyAxis]
    if card not in header and arm in UNBINNED_ARMS:
        return 1
    value = required_card(header, card, path)
    notWhole = ToolError(f"{_where(path)} has a {card} card that is not a whole number: {value!r}")
    if isinstance(value, bool):
        raise notWhole
    try:
        binning = int(value)
    except (TypeError, ValueError) as error:
        raise notWhole from error
    if binning != value and not isinstance(value, str):
        raise notWhole
    if binning < 1:
        raise ToolError(f"{_where(path)} has a {card} card below 1: {value!r}")
    return binning


def spectral_binning(header: fits.Header, arm: str, path: Path | None = None) -> int:
    """*the binning factor along the dispersion of a product: ``BINY`` for VIS, ``BINX`` for NIR*

    The card follows the numpy axis the wavelength runs along (``disp_image_axis``).

    **Key Arguments:**

    - ``header`` -- the FITS header of the product
    - ``arm`` -- the arm of the product, which sets the dispersion axis
    - ``path`` -- the file the header came from, named in the error. Default *None*

    **Return:**

    - ``binning`` -- the number of detector pixels in one binned pixel along the dispersion

    **Raises:**

    - ``ToolError`` -- if the binning card is missing (a NIR frame without one counts as unbinned), not a whole number,
      or below 1

    **Usage:**

    ```python
    binning = spectral_binning(header, "VIS", path)
    ```
    """
    return binning_along(header, disp_image_axis(arm), arm, path)


def soxspipe_version(path: Path = VERSION_PATH) -> str:
    """*the soxspipe version of this checkout, read without importing the package*

    **Key Arguments:**

    - ``path`` -- the generated ``__version__.py``. Default the checkout's own

    **Return:**

    - ``version`` -- the version string, or ``unknown`` when the file is missing or unreadable

    **Usage:**

    ```python
    version = soxspipe_version()
    ```
    """
    try:
        match = re.search(r'__version__\s*=\s*"([^"]+)"', path.read_text(encoding="utf-8"))
    except OSError:
        return UNKNOWN_VERSION
    return match.group(1) if match else UNKNOWN_VERSION


def night_count(mjds: Iterable[float]) -> int:
    """*the number of distinct UTC dates among observation MJDs*

    **Key Arguments:**

    - ``mjds`` -- the observation MJDs

    **Return:**

    - ``count`` -- the number of distinct integer MJDs

    **Usage:**

    ```python
    nights = night_count([60000.1, 60001.4])
    ```
    """
    return len({math.floor(mjd) for mjd in mjds})


def provenance_comments(subject: str, mjds: Sequence[float], version: str, today: str) -> list[str]:
    """*the comment line that records how a printed yaml value was measured*

    **Key Arguments:**

    - ``subject`` -- what was counted, for example ``FRAMES``
    - ``mjds`` -- the observation MJD of each input
    - ``version`` -- the soxspipe version of the code that measured
    - ``today`` -- the measurement date, ISO format

    **Return:**

    - ``lines`` -- the uppercase yaml comment lines, house style

    **Usage:**

    ```python
    lines = provenance_comments("FRAMES", mjds, "0.19.1", "2026-11-22")
    ```
    """
    first, last = math.floor(min(mjds)), math.floor(max(mjds))
    return [
        f"# MEASURED {today} WITH soxspipe {version} FROM {len(mjds)} {subject}"
        f" ON {night_count(mjds)} NIGHTS (MJD {first} TO {last})"
    ]


def arm_block(arm: str, comments: Iterable[str], body: Iterable[str]) -> list[str]:
    """*the lines of one arm's yaml block: the arm key, its comments, then its body, indented*

    **Key Arguments:**

    - ``arm`` -- the arm key, ``VIS`` or ``NIR``
    - ``comments`` -- the comment lines, each starting with ``#``
    - ``body`` -- the yaml lines of the value, already indented relative to the arm key

    **Return:**

    - ``lines`` -- the block, with the comments and body indented under the arm

    **Usage:**

    ```python
    lines = arm_block("VIS", ["# NOTE"], ["key: 1"])
    ```
    """
    return [f"{arm}:", *(f"{YAML_INDENT}{line}" for line in [*comments, *body])]


def join_blocks(lines: Iterable[str]) -> str:
    """*the yaml text of several blocks: the lines joined, with a final newline*

    **Key Arguments:**

    - ``lines`` -- the lines of every block

    **Return:**

    - ``text`` -- the text to print

    **Usage:**

    ```python
    print(join_blocks(arm_block("VIS", [], ["key: 1"])), end="")
    ```
    """
    return "\n".join(lines) + "\n"


def ensure_repo_on_path() -> None:
    """*put this checkout first on ``sys.path``, so that ``import soxspipe`` reads this code, not another install*

    **Usage:**

    ```python
    ensure_repo_on_path()
    from soxspipe.commonutils.ug_wavelength_transfer import attached_arc_offsets
    ```
    """
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))


def today_iso() -> str:
    """*today's UTC date, ISO format*

    **Return:**

    - ``today`` -- the date as ``YYYY-MM-DD``

    **Usage:**

    ```python
    today = today_iso()
    ```
    """
    return datetime.now(UTC).date().isoformat()


def run_main(function: Callable[[Sequence[str]], int], argv: Sequence[str]) -> int:
    """*run a script's main function and turn a ``ToolError`` into a message and exit status 2*

    **Key Arguments:**

    - ``function`` -- the script's ``main``
    - ``argv`` -- the command-line arguments

    **Return:**

    - ``status`` -- the function's status, or 2 after a ``ToolError``

    **Usage:**

    ```python
    sys.exit(run_main(main, sys.argv[1:]))
    ```
    """
    try:
        return function(argv)
    except ToolError as error:
        print(f"error: {error}", file=sys.stderr)
        return EXIT_INPUT_ERROR
