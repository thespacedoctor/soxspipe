#!/usr/bin/env python
"""*measure the nominal SPEC_RES of each SOXS arm and slit from slit arcs*

This is a developer tool for DY-2377. It is not part of a reduction and not in the soxspipe command line.

The Phase 3 keyword `SPEC_RES` is the resolving power through the science slit. It comes from a nominal table per arm
and slit, `spec-res` in `soxs_detector_parameters.yaml`, keyed `SLIT0.5` to `SLIT5.0`. The value is the median arc-line
resolving power over the arm, measured from slit arcs over several nights (DY-2306). `R PIN MEDIAN` is never used,
because the pinhole is a 0.5 arcsec aperture and overstates R at 1.0 arcsec by about 1.8 times. The method is that of
the DY-2297 prototype `linewidth.py`.

**Method**

1. Pair each raw slit arc (`ESO DPR TECH = ECHELLE,SLIT`, `LAMP,WAVE`) with the `DISP_IMAGE` of the same arm that is
   nearest in time, within `--max-gap-days`. The `DISP_IMAGE` holds the `WAVELENGTH`, `SLIT` and `ORDER` maps of the
   night's multi-pinhole solution.
2. Trim the raw arc to the science pixels and keep the pixels within 3 arcsec of the slit centre.
3. Per order, build a median spectrum on half-pixel wavelength bins and find the emission lines. Keep isolated lines
   (no neighbour within 8 pixels) with no flagged or saturated pixel in the fit window.
4. Fit a Gaussian to the pixels of each line. The resolving power of a line is `centre / FWHM`, in nm.
5. Per arm and slit, take the median R over all kept lines of all frames. Slit `SLIT5.0` is not measured. Its value is
   the `SLIT1.0` value times 1.0 / 5.0, and the printed comment says so.

**Usage**

```bash
python tools/measure_spec_res.py --arcs "path/to/raw/*SLIT*.fits" --disp-images "path/to/*_DISP_IMAGE.fits"
python tools/measure_spec_res.py --arcs a.fits b.fits --disp-images night1_IMAGE.fits --saturation 40000
```

It prints the yaml block to paste into `soxspipe/resources/soxs_detector_parameters.yaml`:

```yaml
VIS:
  # MEASURED 2026-11-22 WITH soxspipe 0.19.1 FROM 30 SLIT ARCS ON 5 NIGHTS (MJD 60000 TO 60040)
  spec-res:
    SLIT0.5: 8965
    SLIT1.0: 4906
    SLIT1.5: 3271
    SLIT5.0: 981
```

Raw arcs must be readable by astropy, so uncompressed or gzip. Exit status is 0 when the block was printed and 2 when
the input could not be used.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple

import measurement_common as common
import numpy as np
import pandas as pd
from astropy.io import fits
from measurement_common import ToolError, today_iso
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

FWHM_PER_SIGMA = 2.3548
SLIT_HALF_RANGE_ARCSEC = 3.0
MIN_ORDER_PIXELS = 500
PROFILE_BIN_FRACTION = 0.5
MIN_PEAK_DISTANCE_BINS = 6
ISOLATION_PIXELS = 8
FIT_WINDOW_PIXELS = 5
MIN_FIT_PIXELS = 20
MIN_FWHM_PIXELS = 0.8
MAX_FWHM_PIXELS = 8.0
MAX_CENTRE_DRIFT_PIXELS = 2.0
DEFAULT_PROMINENCE = 30.0
DEFAULT_SATURATION = 40000.0
DEFAULT_MAX_GAP_DAYS = 1.0
MAD_TO_SIGMA = 1.4826
DERIVED_SLIT = "SLIT5.0"
REFERENCE_SLIT = "SLIT1.0"
SLIT_CARDS = {"VIS": "ESO INS VISE NAME", "NIR": "ESO INS NISE NAME"}
DPR_TECH_CARD = "ESO DPR TECH"
DPR_TYPE_CARD = "ESO DPR TYPE"
SLIT_ARC_TECH = "ECHELLE,SLIT"
SLIT_ARC_TYPES = ("LAMP,WAVE", "WAVE,LAMP")
SLIT_WIDTH_PATTERN = re.compile(r"^SLIT(\d+(?:\.\d+)?)$")
DISPERSION_AXIS_INDEX = {"x": 1, "y": 0}


@dataclass(frozen=True)
class LineFit:
    """*one fitted emission line: its order, centre, widths and resolving power*"""

    order: int
    centreNm: float
    fwhmNm: float
    fwhmPx: float
    resolvingPower: float


def gaussian(x, amplitude, centre, sigma, offset):
    """*a Gaussian with a constant offset*"""
    return amplitude * np.exp(-0.5 * ((x - centre) / sigma) ** 2) + offset


def nm_per_pixel(wavelength, mask, axis: int) -> float:
    """*the median wavelength step between neighbouring pixels along the dispersion axis, inside a mask*"""
    lower = tuple(slice(None, -1) if index == axis else slice(None) for index in range(2))
    upper = tuple(slice(1, None) if index == axis else slice(None) for index in range(2))
    both = mask[lower] & mask[upper]
    return float(np.median(np.abs(wavelength[upper] - wavelength[lower])[both]))


def median_profile(wavelengths, values, binWidth: float):
    """*the median of the values in wavelength bins of ``binWidth``, and the bin edges*

    **Return:**

    - ``profile`` -- the median per bin, with empty bins filled with the overall median
    - ``edges`` -- the left edge of each bin
    """
    edges = np.arange(wavelengths.min(), wavelengths.max(), binWidth)
    binIndex = np.digitize(wavelengths, edges)
    medians = pd.Series(values).groupby(binIndex).median()
    profile = medians.reindex(range(1, len(edges) + 1)).to_numpy(dtype=float)
    return np.nan_to_num(profile, nan=np.nanmedian(profile)), edges


def isolated_centres(profile, edges, binWidth: float, nmStep: float, prominence: float):
    """*the wavelengths of the profile's peaks that have no other peak within ``ISOLATION_PIXELS`` pixels*"""
    base = np.nanmedian(profile)
    noise = MAD_TO_SIGMA * np.nanmedian(np.abs(profile - base))
    peaks, _ = find_peaks(profile, prominence=prominence * noise, distance=MIN_PEAK_DISTANCE_BINS)
    centres = edges[peaks] + binWidth / 2
    return [
        centre
        for index, centre in enumerate(centres)
        if not np.any(np.abs(np.delete(centres, index) - centre) < ISOLATION_PIXELS * nmStep)
    ]


def fit_line(wavelengths, values, good, centre: float, nmStep: float, order: int) -> LineFit | None:
    """*fit a Gaussian to the pixels around ``centre``; None when the line cannot be used*

    A line is not used when the window holds a flagged or saturated pixel, too few pixels, a failed fit, a negative
    amplitude, a width outside 0.8 to 8 pixels or a centre that moved more than 2 pixels.
    """
    window = np.abs(wavelengths - centre) < FIT_WINDOW_PIXELS * nmStep
    if window.sum() < MIN_FIT_PIXELS or not good[window].all():
        return None
    w, f = wavelengths[window], values[window]
    start = [f.max() - np.median(f), centre, nmStep, np.median(f)]
    try:
        params, _ = curve_fit(gaussian, w, f, p0=start)
    except (RuntimeError, ValueError):
        return None
    fwhm = FWHM_PER_SIGMA * abs(params[2])
    isWidthOk = MIN_FWHM_PIXELS * nmStep < fwhm < MAX_FWHM_PIXELS * nmStep
    if not isWidthOk or params[0] <= 0 or abs(params[1] - centre) > MAX_CENTRE_DRIFT_PIXELS * nmStep:
        return None
    return LineFit(order, float(params[1]), float(fwhm), float(fwhm / nmStep), float(params[1] / fwhm))


def measure_order(data, good, wavelength, mask, axis: int, order: int, prominence: float) -> list[LineFit]:
    """*the fitted lines of one order*"""
    nmStep = nm_per_pixel(wavelength, mask, axis)
    sortIndex = np.argsort(wavelength[mask])
    w = wavelength[mask][sortIndex]
    f = data[mask][sortIndex]
    ok = good[mask][sortIndex]
    binWidth = nmStep * PROFILE_BIN_FRACTION
    profile, edges = median_profile(w, f, binWidth)
    centres = isolated_centres(profile, edges, binWidth, nmStep, prominence)
    fitted = (fit_line(w, f, ok, centre, nmStep, order) for centre in centres)
    return [line for line in fitted if line is not None]


def measure_lines(
    data,
    good,
    wavelength,
    slit,
    order,
    dispersionAxis: str,
    slitRange: tuple[float, float] = (-SLIT_HALF_RANGE_ARCSEC, SLIT_HALF_RANGE_ARCSEC),
    prominence: float = DEFAULT_PROMINENCE,
) -> list[LineFit]:
    """*the isolated, unsaturated emission lines of an arc frame, with their resolving power*

    **Key Arguments:**

    - ``data`` -- the arc frame, trimmed to the shape of the maps
    - ``good`` -- True where the pixel is usable (not flagged, not saturated)
    - ``wavelength`` -- the ``WAVELENGTH`` map of the ``DISP_IMAGE``, in nm
    - ``slit`` -- the ``SLIT`` map, in arcsec
    - ``order`` -- the ``ORDER`` map; 0 is no order
    - ``dispersionAxis`` -- ``x`` or ``y``, the detector axis the wavelength runs along
    - ``slitRange`` -- the slit positions to use, in arcsec. Default within 3 arcsec of the centre
    - ``prominence`` -- the peak prominence, in units of the profile's noise. Default 30

    **Return:**

    - ``lines`` -- one ``LineFit`` per kept line
    """
    axis = DISPERSION_AXIS_INDEX[dispersionAxis]
    inSlit = (slit > slitRange[0]) & (slit < slitRange[1]) & np.isfinite(data) & (wavelength > 0)
    lines: list[LineFit] = []
    for orderNumber in np.unique(order[order > 0]):
        mask = (order == orderNumber) & inSlit
        if mask.sum() >= MIN_ORDER_PIXELS:
            lines.extend(measure_order(data, good, wavelength, mask, axis, int(orderNumber), prominence))
    return lines


def median_resolving_power(lines: Sequence[LineFit]) -> float:
    """*the median resolving power over fitted lines*"""
    return float(np.median([line.resolvingPower for line in lines]))


def slit_width(slitName: str) -> float:
    """*the width in arcsec of a slit name such as ``SLIT1.5``*

    **Raises:**

    - ``ToolError`` -- if the name is not ``SLIT`` followed by a width
    """
    match = SLIT_WIDTH_PATTERN.match(slitName)
    if match is None:
        raise ToolError(f"the slit name {slitName!r} is not of the form SLIT<width>")
    return float(match.group(1))


def derive_wide_slit(valuesBySlit: Mapping[str, int]) -> dict[str, int]:
    """*add the derived 5.0 arcsec value: the 1.0 arcsec value scaled by 1 / width*

    The measurement at 5.0 arcsec is unreliable (few usable lines), so the table derives it (DY-2306).

    **Key Arguments:**

    - ``valuesBySlit`` -- ``{slit name: integer R}`` of the measured slits

    **Return:**

    - ``table`` -- a copy with ``SLIT5.0`` added, or unchanged when ``SLIT1.0`` was not measured
    """
    table = dict(valuesBySlit)
    if REFERENCE_SLIT in table:
        table[DERIVED_SLIT] = round(table[REFERENCE_SLIT] * slit_width(REFERENCE_SLIT) / slit_width(DERIVED_SLIT))
    return table


class ArcRef(NamedTuple):
    """*a raw slit arc: where it is, its arm, slit, observation MJD and binning*"""

    path: Path
    arm: str
    slit: str
    mjd: float
    binX: int
    binY: int


class MapRef(NamedTuple):
    """*a ``DISP_IMAGE``: where it is, its arm and the MJD of the arc it was solved from*"""

    path: Path
    arm: str
    mjd: float


def is_slit_arc(header) -> bool:
    """*whether a raw frame header is a slit arc lamp frame (``ECHELLE,SLIT`` with ``LAMP,WAVE`` or ``WAVE,LAMP``)*"""
    return header.get(DPR_TECH_CARD) == SLIT_ARC_TECH and header.get(DPR_TYPE_CARD) in SLIT_ARC_TYPES


def science_slices(arm: str, binX: int, binY: int) -> tuple[slice, slice]:
    """*the row and column slices of the science pixels of a raw frame, for its binning*

    **Key Arguments:**

    - ``arm`` -- ``VIS`` or ``NIR``
    - ``binX`` -- the binning along x
    - ``binY`` -- the binning along y

    **Return:**

    - ``rowSlice``, ``columnSlice`` -- the science-pixel limits of the detector parameters, divided by the binning
    """
    pixels = common.load_yaml(common.DETECTOR_PARAMETERS_PATH)[arm]["science-pixels"]
    rows, columns = pixels["rows"], pixels["columns"]
    return (
        slice(rows["start"] // binY, rows["end"] // binY),
        slice(columns["start"] // binX, columns["end"] // binX),
    )


def read_arc_ref(path: Path) -> ArcRef | None:
    """*read the header of a raw frame; None when it is not a slit arc*

    **Raises:**

    - ``ToolError`` -- if the file cannot be read or a card is missing
    """
    header = common.read_header(path)
    if not is_slit_arc(header):
        return None
    arm = common.arm_of(header)
    slit = header.get(SLIT_CARDS[arm])
    if slit is None or "MJD-OBS" not in header:
        raise ToolError(f"{path} has no {SLIT_CARDS[arm]} or MJD-OBS card")
    return ArcRef(
        Path(path),
        arm,
        str(slit).strip(),
        float(header["MJD-OBS"]),
        int(header["ESO DET BINX"]),
        int(header["ESO DET BINY"]),
    )


def read_map_ref(path: Path) -> MapRef:
    """*read the arm and MJD of a ``DISP_IMAGE`` from its primary header*

    **Raises:**

    - ``ToolError`` -- if the file cannot be read or a card is missing
    """
    header = common.read_header(path)
    if "MJD-OBS" not in header:
        raise ToolError(f"{path} has no MJD-OBS card")
    return MapRef(Path(path), common.arm_of(header), float(header["MJD-OBS"]))


def pair_arcs_with_maps(
    arcs: Sequence[ArcRef], maps: Sequence[MapRef], maxGapDays: float
) -> tuple[list[tuple[ArcRef, MapRef]], list[ArcRef]]:
    """*pair each slit arc with the ``DISP_IMAGE`` of its arm that is nearest in time*

    **Key Arguments:**

    - ``arcs`` -- the raw slit arcs
    - ``maps`` -- the ``DISP_IMAGE`` products
    - ``maxGapDays`` -- the largest time gap that still counts as the same night's solution

    **Return:**

    - ``pairs`` -- ``(arc, map)`` for each arc with a map in range
    - ``unpaired`` -- the arcs with no map of their arm within ``maxGapDays``
    """
    pairs, unpaired = [], []
    for arc in arcs:
        candidates = [mapRef for mapRef in maps if mapRef.arm == arc.arm and abs(mapRef.mjd - arc.mjd) <= maxGapDays]
        if candidates:
            pairs.append((arc, min(candidates, key=lambda mapRef: abs(mapRef.mjd - arc.mjd))))
        else:
            unpaired.append(arc)
    return pairs, unpaired


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
    """
    rows, columns = science_slices(arc.arm, arc.binX, arc.binY)
    data = fits.getdata(arc.path, 0).astype(float)[rows, columns]
    with fits.open(mapRef.path) as maps:
        wavelength = np.array(maps["WAVELENGTH"].data, dtype=float)
        slit = np.array(maps["SLIT"].data, dtype=float)
        order = np.array(maps["ORDER"].data)
    if data.shape != wavelength.shape:
        raise ToolError(
            f"the trimmed arc {arc.path} has shape {data.shape}, "
            f"but the maps of {mapRef.path} have shape {wavelength.shape}"
        )
    return measure_lines(data, data < saturation, wavelength, slit, order, common.dispersion_axis(arc.arm))


class ArcMeasurement(NamedTuple):
    """*the lines fitted on one slit arc*"""

    arm: str
    slit: str
    mjd: float
    lines: list[LineFit]


@dataclass(frozen=True)
class SlitSummary:
    """*the nominal resolving power of one arm and slit, with what it was measured from*"""

    resolvingPower: int
    lineCount: int
    frameCount: int
    mjds: list[float]
    nightMedians: dict[int, int]


def summarise(measurements: Sequence[ArcMeasurement]) -> dict[str, dict[str, SlitSummary]]:
    """*the median resolving power per arm and slit over every line of every frame*

    **Key Arguments:**

    - ``measurements`` -- the lines of each slit arc. A frame with no line adds nothing.

    **Return:**

    - ``summary`` -- ``{arm: {slit: SlitSummary}}``. ``nightMedians`` is the median R of each UTC date, for the comment.

    **Raises:**

    - ``ToolError`` -- if no frame found a line
    """
    groups: dict[tuple[str, str], list[ArcMeasurement]] = defaultdict(list)
    for measurement in measurements:
        if measurement.lines:
            groups[(measurement.arm, measurement.slit)].append(measurement)
    if not groups:
        raise ToolError("no arc line was measured on any frame")
    summary: dict[str, dict[str, SlitSummary]] = defaultdict(dict)
    for (arm, slit), frames in groups.items():
        allLines = [line for frame in frames for line in frame.lines]
        nights: dict[int, list[LineFit]] = defaultdict(list)
        for frame in frames:
            nights[math.floor(frame.mjd)].extend(frame.lines)
        summary[arm][slit] = SlitSummary(
            resolvingPower=round(median_resolving_power(allLines)),
            lineCount=len(allLines),
            frameCount=len(frames),
            mjds=[frame.mjd for frame in frames],
            nightMedians={night: round(median_resolving_power(lines)) for night, lines in sorted(nights.items())},
        )
    return dict(summary)


def format_yaml(summary: Mapping[str, Mapping[str, SlitSummary]], version: str, today: str) -> str:
    """*the yaml block to paste into ``soxs_detector_parameters.yaml``*

    **Key Arguments:**

    - ``summary`` -- the measured arms and slits
    - ``version`` -- the soxspipe version, for the provenance comment
    - ``today`` -- the measurement date, ISO format

    **Return:**

    - ``text`` -- one block per arm: provenance comments, then ``spec-res`` keyed by slit, ``SLIT5.0`` derived
    """
    lines = []
    for arm, slits in summary.items():
        ordered = sorted(slits, key=slit_width)
        mjds = [mjd for slit in ordered for mjd in slits[slit].mjds]
        lines.append(f"{arm}:")
        lines.extend(f"  {line}" for line in common.provenance_comments("SLIT ARCS", mjds, version, today))
        for slit in ordered:
            lines.append(f"  # {slit}: {_describe(slits[slit])}")
        table = derive_wide_slit({slit: slits[slit].resolvingPower for slit in ordered})
        if DERIVED_SLIT in table:
            lines.append(f"  # {DERIVED_SLIT} IS DERIVED, NOT MEASURED: THE {REFERENCE_SLIT} VALUE TIMES 1.0 / 5.0")
        else:
            lines.append(f"  # {DERIVED_SLIT} NOT DERIVED: NO {REFERENCE_SLIT} MEASUREMENT")
        lines.append("  spec-res:")
        lines.extend(f"    {slit}: {table[slit]}" for slit in sorted(table, key=slit_width))
    return "\n".join(lines) + "\n"


def _describe(slitSummary: SlitSummary) -> str:
    """*the comment text of one slit: line and frame counts and the median R of each night*"""
    nights = ", ".join(f"{power} (MJD {night})" for night, power in slitSummary.nightMedians.items())
    return f"{slitSummary.lineCount} LINES FROM {slitSummary.frameCount} FRAMES; MEDIAN R PER NIGHT {nights}"


def collect_inputs(arcPatterns: Sequence[str], mapPatterns: Sequence[str]) -> tuple[list[ArcRef], list[MapRef]]:
    """*read the slit arcs and ``DISP_IMAGE`` files named by the patterns*

    A frame that is not a slit arc, and a ``SLIT5.0`` arc (its value is derived), are named on stderr and left out.

    **Raises:**

    - ``ToolError`` -- if a pattern matches nothing or a file cannot be read
    """
    arcs = []
    for path in common.expand_paths(arcPatterns):
        arc = read_arc_ref(path)
        if arc is None:
            print(f"skipped {path}: not a slit arc", file=sys.stderr)
        elif arc.slit == DERIVED_SLIT:
            print(f"skipped {path}: {DERIVED_SLIT} is derived from {REFERENCE_SLIT}, not measured", file=sys.stderr)
        else:
            arcs.append(arc)
    maps = [read_map_ref(path) for path in common.expand_paths(mapPatterns)]
    return arcs, maps


def build_parser() -> argparse.ArgumentParser:
    """*the command-line parser of the script*"""
    parser = argparse.ArgumentParser(
        description="Measure the nominal SPEC_RES of each SOXS arm and slit from slit arcs."
    )
    parser.add_argument("--arcs", nargs="+", required=True, help="raw slit arc FITS files or glob patterns")
    parser.add_argument("--disp-images", nargs="+", required=True, help="DISP_IMAGE files or glob patterns")
    parser.add_argument("--max-gap-days", type=float, default=DEFAULT_MAX_GAP_DAYS, help="largest arc to map gap")
    parser.add_argument(
        "--saturation", type=float, default=DEFAULT_SATURATION, help="pixels at or above this ADU are bad"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """*measure the lines of every paired slit arc and print the yaml block*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments. Default the process arguments

    **Return:**

    - ``status`` -- 0 when the block was printed

    **Raises:**

    - ``ToolError`` -- if no arc could be paired with a ``DISP_IMAGE`` or no line was measured
    """
    arguments = build_parser().parse_args(argv)
    arcs, maps = collect_inputs(arguments.arcs, arguments.disp_images)
    pairs, unpaired = pair_arcs_with_maps(arcs, maps, arguments.max_gap_days)
    for arc in unpaired:
        print(
            f"skipped {arc.path}: no DISP_IMAGE of the {arc.arm} arm within {arguments.max_gap_days} days",
            file=sys.stderr,
        )
    if not pairs:
        raise ToolError("no slit arc was paired with a DISP_IMAGE")
    measurements = [
        ArcMeasurement(arc.arm, arc.slit, arc.mjd, measure_arc(arc, mapRef, arguments.saturation))
        for arc, mapRef in pairs
    ]
    print(format_yaml(summarise(measurements), common.soxspipe_version(), today_iso()), end="")
    return 0


if __name__ == "__main__":
    sys.exit(common.run_main(main, sys.argv[1:]))
