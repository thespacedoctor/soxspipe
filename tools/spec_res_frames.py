"""*read slit arcs, reduced stare frames and ``DISP_IMAGE`` files, and measure their lines, for `measure_spec_res.py`*

A frame is described by a small record (`ArcRef`, `SkyRef`, `MapRef`) read from its header. `pair_frames_with_maps`
pairs frames with the `DISP_IMAGE` of the same arm that is nearest in time, and `measure_arc` and `measure_sky` run the
line measurement of `line_fitting.py` through that image's wavelength, slit and order maps.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import measurement_common as common
import numpy as np
import numpy.typing as npt
from astropy.io import fits
from line_fitting import ARC_PROMINENCE, ARC_SLIT_RANGE_ARCSEC, LineFit, measure_lines
from measurement_common import ToolError

SKY_SLIT_RANGE_ARCSEC = (2.0, 5.0)
SKY_PROMINENCE = 20.0
SLIT_CARDS = {"VIS": "ESO INS VISE NAME", "NIR": "ESO INS NISE NAME"}
DPR_TECH_CARD = "ESO DPR TECH"
DPR_TYPE_CARD = "ESO DPR TYPE"
SLIT_ARC_TECH = "ECHELLE,SLIT"
SLIT_ARC_TYPES = ("LAMP,WAVE", "WAVE,LAMP")
SKYSUB_SUFFIX = "_SKYSUB.fits"
SKYMODEL_SUFFIX = "_SKYMODEL.fits"


class _TimedFrame(Protocol):
    """*what pairing needs of a frame: where it is, its arm and its observation time*"""

    path: Path
    arm: str
    mjd: float


@dataclass(frozen=True)
class ArcRef:
    """*a raw slit arc: where it is, its arm, slit, observation MJD and binning*"""

    path: Path
    arm: str
    slit: str
    mjd: float
    binX: int
    binY: int


@dataclass(frozen=True)
class SkyRef:
    """*a reduced stare frame: its sky-subtracted frame, its sky model, arm, slit and observation MJD*"""

    path: Path
    modelPath: Path
    arm: str
    slit: str
    mjd: float


@dataclass(frozen=True)
class MapRef:
    """*a ``DISP_IMAGE``: where it is, its arm and the MJD of the arc it was solved from*"""

    path: Path
    arm: str
    mjd: float


def is_slit_arc(header: fits.Header) -> bool:
    """*whether a raw frame header is a slit arc lamp frame (``ECHELLE,SLIT`` with ``LAMP,WAVE`` or ``WAVE,LAMP``)*

    **Key Arguments:**

    - ``header`` -- the raw frame header

    **Return:**

    - ``isSlitArc`` -- *True* for a slit arc

    **Usage:**

    ```python
    if is_slit_arc(header):
        ...
    ```
    """
    return header.get(DPR_TECH_CARD) == SLIT_ARC_TECH and header.get(DPR_TYPE_CARD) in SLIT_ARC_TYPES


def science_slices(arm: str, binX: int, binY: int) -> tuple[slice, slice]:
    """*the row and column slices of the science pixels of a raw frame, for its binning*

    **Key Arguments:**

    - ``arm`` -- ``VIS`` or ``NIR``
    - ``binX`` -- the binning along x
    - ``binY`` -- the binning along y

    **Return:**

    - ``rowSlice``, ``columnSlice`` -- the science-pixel limits of the detector parameters, divided by the binning

    **Usage:**

    ```python
    rows, columns = science_slices("VIS", 1, 1)
    ```
    """
    pixels = common.detector_parameters()[arm]["science-pixels"]
    rows, columns = pixels["rows"], pixels["columns"]
    return (
        slice(rows["start"] // binY, rows["end"] // binY),
        slice(columns["start"] // binX, columns["end"] // binX),
    )


def _identity(header: fits.Header, path: Path) -> tuple[str, str, float]:
    """*the arm, slit name and observation MJD of a frame, raising a ``ToolError`` that names a missing card*"""
    arm = common.arm_of(header, path)
    slit = str(common.required_card(header, SLIT_CARDS[arm], path)).strip()
    return arm, slit, common.read_mjd_obs(header, path)


def read_arc_ref(path: Path) -> ArcRef | None:
    """*read the header of a raw frame; None when it is not a slit arc*

    **Key Arguments:**

    - ``path`` -- the raw frame

    **Return:**

    - ``arc`` -- the slit arc's identity, or *None* for any other frame

    **Raises:**

    - ``ToolError`` -- if the file cannot be read or a card is missing

    **Usage:**

    ```python
    arc = read_arc_ref(Path("raw.fits"))
    ```
    """
    header = common.read_header(path)
    if not is_slit_arc(header):
        return None
    arm, slit, mjd = _identity(header, path)
    binX, binY = common.binning_along(header, 1, arm, path), common.binning_along(header, 0, arm, path)
    return ArcRef(Path(path), arm, slit, mjd, binX, binY)


def read_sky_ref(path: Path) -> SkyRef | None:
    """*read the header of a ``_SKYSUB.fits`` product and find its sky model; None when the model is missing*

    **Key Arguments:**

    - ``path`` -- the sky-subtracted frame, named ``<stem>_SKYSUB.fits``

    **Return:**

    - ``sky`` -- the frame's identity with the path of its ``_SKYMODEL.fits``, or *None* when that file does not exist

    **Raises:**

    - ``ToolError`` -- if the name does not end in ``_SKYSUB.fits``, or the file cannot be read or lacks a card

    **Usage:**

    ```python
    sky = read_sky_ref(Path("stare_SKYSUB.fits"))
    ```
    """
    if not path.name.endswith(SKYSUB_SUFFIX):
        raise ToolError(f"{path} is not a {SKYSUB_SUFFIX} product")
    arm, slit, mjd = _identity(common.read_header(path), path)
    modelPath = path.with_name(path.name[: -len(SKYSUB_SUFFIX)] + SKYMODEL_SUFFIX)
    return SkyRef(path, modelPath, arm, slit, mjd) if modelPath.exists() else None


def read_map_ref(path: Path) -> MapRef:
    """*read the arm and MJD of a ``DISP_IMAGE`` from its primary header*

    **Key Arguments:**

    - ``path`` -- the ``DISP_IMAGE`` file

    **Return:**

    - ``mapRef`` -- the arm and MJD

    **Raises:**

    - ``ToolError`` -- if the file cannot be read or a card is missing

    **Usage:**

    ```python
    mapRef = read_map_ref(Path("night_DISP_IMAGE.fits"))
    ```
    """
    header = common.read_header(path)
    return MapRef(Path(path), common.arm_of(header, path), common.read_mjd_obs(header, path))


def pair_frames_with_maps[Frame: _TimedFrame](
    frames: Sequence[Frame], maps: Sequence[MapRef], maxGapDays: float
) -> tuple[list[tuple[Frame, MapRef]], list[Frame]]:
    """*pair each frame with the ``DISP_IMAGE`` of its arm that is nearest in time*

    **Key Arguments:**

    - ``frames`` -- the arcs or sky frames, each with an ``arm`` and an ``mjd``
    - ``maps`` -- the ``DISP_IMAGE`` products
    - ``maxGapDays`` -- the largest time gap that still counts as the same night's solution

    **Return:**

    - ``pairs`` -- ``(frame, map)`` for each frame with a map in range
    - ``unpaired`` -- the frames with no map of their arm within ``maxGapDays``

    **Usage:**

    ```python
    pairs, unpaired = pair_frames_with_maps(arcs, maps, maxGapDays=1.0)
    ```
    """
    pairs, unpaired = [], []
    for frame in frames:
        candidates = [m for m in maps if m.arm == frame.arm and abs(m.mjd - frame.mjd) <= maxGapDays]
        if candidates:
            pairs.append((frame, min(candidates, key=lambda m: abs(m.mjd - frame.mjd))))
        else:
            unpaired.append(frame)
    return pairs, unpaired


def _measure_through_maps(
    framePath: Path,
    data: npt.NDArray[np.float64],
    isGood: npt.NDArray[np.bool_],
    arm: str,
    mapRef: MapRef,
    slitRange: tuple[float, float],
    prominence: float,
) -> list[LineFit]:
    """*measure a trimmed frame through the ``WAVELENGTH``, ``SLIT`` and ``ORDER`` maps of a ``DISP_IMAGE``*"""
    with fits.open(mapRef.path) as maps:
        wavelength = np.array(maps["WAVELENGTH"].data, dtype=float)
        slit = np.array(maps["SLIT"].data, dtype=float)
        order = np.array(maps["ORDER"].data)
    if data.shape != wavelength.shape:
        raise ToolError(
            f"the trimmed frame {framePath} has shape {data.shape}, "
            f"but the maps of {mapRef.path} have shape {wavelength.shape}"
        )
    return measure_lines(
        data, isGood, wavelength, slit, order, common.disp_image_axis(arm), slitRange=slitRange, prominence=prominence
    )


def measure_arc(arc: ArcRef, mapRef: MapRef, saturation: float) -> list[LineFit]:
    """*measure the lines of one raw slit arc through the maps of its ``DISP_IMAGE``*

    **Key Arguments:**

    - ``arc`` -- the raw slit arc
    - ``mapRef`` -- the ``DISP_IMAGE`` of the same arm
    - ``saturation`` -- the level in ADU at and above which a pixel is not used

    **Return:**

    - ``lines`` -- the fitted lines

    **Raises:**

    - ``ToolError`` -- if the trimmed frame and the maps have different shapes

    **Usage:**

    ```python
    lines = measure_arc(arc, mapRef, saturation=40000.0)
    ```
    """
    rows, columns = science_slices(arc.arm, arc.binX, arc.binY)
    data = fits.getdata(arc.path, 0).astype(float)[rows, columns]
    return _measure_through_maps(
        arc.path, data, data < saturation, arc.arm, mapRef, ARC_SLIT_RANGE_ARCSEC, ARC_PROMINENCE
    )


def measure_sky(sky: SkyRef, mapRef: MapRef) -> list[LineFit]:
    """*measure the night-sky lines of a reduced stare frame through the maps of its ``DISP_IMAGE``*

    The data is the sky-subtracted frame plus the sky model, so the sky lines are back in. Pixels that ``QUAL`` flags
    are left out. Only the pixels with ``2 < slit < 5`` arcsec are used, on the positive side of the slit only, as the
    DY-2297 prototype did, so the numbers stay comparable with it. The negative side is not used.

    **Key Arguments:**

    - ``sky`` -- the stare frame
    - ``mapRef`` -- the ``DISP_IMAGE`` of the same arm

    **Return:**

    - ``lines`` -- the fitted lines

    **Raises:**

    - ``ToolError`` -- if the frame and the maps have different shapes

    **Usage:**

    ```python
    lines = measure_sky(sky, mapRef)
    ```
    """
    with fits.open(sky.path) as subtracted:
        data = subtracted[0].data.astype(float) + fits.getdata(sky.modelPath, 0).astype(float)
        isGood = subtracted["QUAL"].data == 0
    return _measure_through_maps(sky.path, data, isGood, sky.arm, mapRef, SKY_SLIT_RANGE_ARCSEC, SKY_PROMINENCE)
