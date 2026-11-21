"""*shared helpers of the developer measurement scripts*

The scripts `measure_log_wavelength_step.py`, `measure_spec_res.py` and `measure_sky_shift_uncertainty.py` measure
instrument values that ship in `soxspipe/resources/soxs_detector_parameters.yaml`. They are developer tools, not part
of a reduction and not in the package CLI. This module holds what the three share: path expansion, the arm and the
spectral binning of a product, the provenance comment that goes above each printed yaml block, and the exit-status
convention of the other tools (0 done, 2 the input could not be used).
"""

from __future__ import annotations

import glob
import math
import re
import sys
from collections.abc import Callable, Iterable, Sequence
from datetime import UTC, datetime
from pathlib import Path

import yaml
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[1]
DETECTOR_PARAMETERS_PATH = REPO_ROOT / "soxspipe" / "resources" / "soxs_detector_parameters.yaml"
KEYWORDS_PATH = REPO_ROOT / "soxspipe" / "resources" / "soxs_keywords.yaml"
ARM_CARD = "ESO SEQ ARM"
SPECTRAL_BINNING_CARDS = {"x": "ESO DET BINX", "y": "ESO DET BINY"}
VERSION_PATH = REPO_ROOT / "soxspipe" / "__version__.py"
EXIT_INPUT_ERROR = 2
UNKNOWN_VERSION = "unknown"


class ToolError(Exception):
    """*an input problem the measurement scripts report and stop on*"""


def load_yaml(path: Path) -> dict:
    """*read a packaged yaml resource*

    **Key Arguments:**

    - ``path`` -- the yaml file to read

    **Return:**

    - ``content`` -- the parsed mapping
    """
    with open(path, encoding="utf-8") as stream:
        return yaml.safe_load(stream)


def expand_paths(patterns: Iterable[str]) -> list[Path]:
    """*expand glob patterns, with a leading ``~``, into one sorted list of unique paths*

    **Key Arguments:**

    - ``patterns`` -- file paths or glob patterns

    **Return:**

    - ``paths`` -- the matching paths, sorted and without duplicates

    **Raises:**

    - ``ToolError`` -- if a pattern matches no file
    """
    paths: set[Path] = set()
    for pattern in patterns:
        matches = glob.glob(str(Path(pattern).expanduser()))
        if not matches:
            raise ToolError(f"no file matches {pattern!r}")
        paths.update(Path(match) for match in matches)
    return sorted(paths)


def read_header(path: Path):
    """*the primary header of a FITS file*

    **Raises:**

    - ``ToolError`` -- if the file cannot be read as FITS
    """
    try:
        return fits.getheader(path)
    except OSError as error:
        raise ToolError(f"cannot read {path}: {error}") from error


def arm_of(header) -> str:
    """*the arm of a product, ``VIS`` or ``NIR``, from the sequencer arm card*

    **Key Arguments:**

    - ``header`` -- the FITS header of the product

    **Return:**

    - ``arm`` -- the arm, as it is keyed in ``soxs_detector_parameters.yaml``

    **Raises:**

    - ``ToolError`` -- if the card is missing or the arm has no detector parameters
    """
    arm = header.get(ARM_CARD)
    if arm is None:
        raise ToolError(f"the header has no {ARM_CARD} card")
    arm = str(arm).strip().upper()
    if arm not in load_yaml(DETECTOR_PARAMETERS_PATH):
        raise ToolError(f"{ARM_CARD} is {arm!r}, which has no entry in soxs_detector_parameters.yaml")
    return arm


def dispersion_axis(arm: str) -> str:
    """*the detector axis the arm disperses along, ``x`` or ``y``, from the detector parameters*"""
    return load_yaml(DETECTOR_PARAMETERS_PATH)[arm]["dispersion-axis"]


def spectral_binning(header, arm: str) -> int:
    """*the binning factor along the dispersion axis of a product*

    **Key Arguments:**

    - ``header`` -- the FITS header of the product
    - ``arm`` -- the arm of the product, which sets the dispersion axis

    **Return:**

    - ``binning`` -- the number of detector pixels in one binned pixel along the dispersion

    **Raises:**

    - ``ToolError`` -- if the binning card is missing or below 1
    """
    card = SPECTRAL_BINNING_CARDS[dispersion_axis(arm)]
    binning = header.get(card)
    if binning is None or binning < 1:
        raise ToolError(f"the header has no usable {card} card (value {binning!r})")
    return int(binning)


def soxspipe_version(path: Path = VERSION_PATH) -> str:
    """*the soxspipe version of this checkout, read without importing the package*

    **Key Arguments:**

    - ``path`` -- the generated ``__version__.py``. Default the checkout's own

    **Return:**

    - ``version`` -- the version string, or ``unknown`` when the file is missing or unreadable
    """
    try:
        match = re.search(r'__version__\s*=\s*"([^"]+)"', path.read_text(encoding="utf-8"))
    except OSError:
        return UNKNOWN_VERSION
    return match.group(1) if match else UNKNOWN_VERSION


def night_count(mjds: Iterable[float]) -> int:
    """*the number of distinct UTC dates among observation MJDs*"""
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
    """
    first, last = math.floor(min(mjds)), math.floor(max(mjds))
    return [
        f"# MEASURED {today} WITH soxspipe {version} FROM {len(mjds)} {subject}"
        f" ON {night_count(mjds)} NIGHTS (MJD {first} TO {last})"
    ]


def ensure_repo_on_path() -> None:
    """*put this checkout first on ``sys.path``, so that ``import soxspipe`` reads this code, not another install*"""
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))


def today_iso() -> str:
    """*today's UTC date, ISO format*"""
    return datetime.now(UTC).date().isoformat()


def run_main(function: Callable[[Sequence[str]], int], argv: Sequence[str]) -> int:
    """*run a script's main function and turn a ``ToolError`` into a message and exit status 2*"""
    try:
        return function(argv)
    except ToolError as error:
        print(f"error: {error}", file=sys.stderr)
        return EXIT_INPUT_ERROR
