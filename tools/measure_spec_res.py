#!/usr/bin/env python
"""*measure the nominal SPEC_RES of each SOXS arm and slit from slit arcs, with a sky-line cross-check*

This is a developer tool for DY-2377. It is not part of a reduction and not in the soxspipe command line.

The Phase 3 keyword `SPEC_RES` is the resolving power through the science slit. It comes from a nominal table per arm
and slit, `spec-res` in `soxs_detector_parameters.yaml`, keyed `SLIT0.5` to `SLIT5.0`. The value is the median arc-line
resolving power over the arm, measured from slit arcs over several nights (DY-2306). `R PIN MEDIAN` is never used,
because the pinhole is a 0.5 arcsec aperture and overstates R at 1.0 arcsec by about 1.8 times. The method is that of
the DY-2297 prototype `linewidth.py`, in `line_fitting.py`.

**Method**

1. Pair each raw slit arc (`ESO DPR TECH = ECHELLE,SLIT`, `LAMP,WAVE`) with the `DISP_IMAGE` of the same arm that is
   nearest in time, within `--max-gap-days`. The `DISP_IMAGE` holds the `WAVELENGTH`, `SLIT` and `ORDER` maps of the
   night's multi-pinhole solution.
2. Trim the raw arc to the science pixels. Find the isolated emission lines within 3 arcsec of the slit centre and fit a
   Gaussian to each. The resolving power of a line is `centre / FWHM`, in nm.
3. Per arm and slit, take the median R over all kept lines of all frames. Slit `SLIT5.0` is not measured. Its value is
   the `SLIT1.0` value times 1.0 / 5.0, and the printed comment says so.

**Sky-line cross-check (optional)**

With `--sky-products`, the script also measures night-sky lines in reduced stare frames. The data is the `_SKYSUB.fits`
frame plus its `_SKYMODEL.fits` frame, with the pixels that `QUAL` flags left out. Only the pixels with `2 < slit < 5`
arcsec are used, on the positive side of the slit only, as in the DY-2297 prototype, so the numbers stay comparable
with it. The `DISP_IMAGE` pairing is the same. The median sky-line R per arm and slit is printed as a yaml comment
beside the table, with its ratio to the arc value. It never changes the table.

**Usage**

```bash
python tools/measure_spec_res.py --arcs "path/to/raw/*SLIT*.fits" --disp-images "path/to/*_DISP_IMAGE.fits"
python tools/measure_spec_res.py --arcs a.fits b.fits --disp-images night1_IMAGE.fits --saturation 40000
python tools/measure_spec_res.py --arcs "raw/*.fits" --disp-images "maps/*.fits" --sky-products "stare/*_SKYSUB.fits"
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

import measurement_common as common
from line_fitting import LineFit, median_resolving_power
from measurement_common import ToolError, today_iso
from spec_res_frames import (
    SKYMODEL_SUFFIX,
    ArcRef,
    MapRef,
    SkyRef,
    measure_arc,
    measure_sky,
    pair_frames_with_maps,
    read_arc_ref,
    read_map_ref,
    read_sky_ref,
)

DEFAULT_SATURATION = 40000.0
DEFAULT_MAX_GAP_DAYS = 1.0
DERIVED_SLIT = "SLIT5.0"
REFERENCE_SLIT = "SLIT1.0"
SLIT_WIDTH_PATTERN = re.compile(r"^SLIT(\d+(?:\.\d+)?)$")


@dataclass(frozen=True)
class FrameMeasurement:
    """*the lines fitted on one frame*"""

    arm: str
    slit: str
    mjd: float
    lines: tuple[LineFit, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "lines", tuple(self.lines))


@dataclass(frozen=True)
class SlitSummary:
    """*the resolving power of one arm and slit, with what it was measured from*"""

    resolvingPower: int
    lineCount: int
    frameCount: int
    mjds: tuple[float, ...]
    nightMedians: Mapping[int, int]

    def __post_init__(self) -> None:
        object.__setattr__(self, "mjds", tuple(self.mjds))
        object.__setattr__(self, "nightMedians", common.frozen_mapping(self.nightMedians))


Summary = dict[str, dict[str, SlitSummary]]


def slit_width(slitName: str) -> float:
    """*the width in arcsec of a slit name such as ``SLIT1.5``*

    **Key Arguments:**

    - ``slitName`` -- the slit name, as the raw header holds it

    **Return:**

    - ``width`` -- the width in arcsec

    **Raises:**

    - ``ToolError`` -- if the name is not ``SLIT`` followed by a width

    **Usage:**

    ```python
    width = slit_width("SLIT1.5")
    ```
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

    **Usage:**

    ```python
    table = derive_wide_slit({"SLIT1.0": 4906})
    ```
    """
    table = dict(valuesBySlit)
    if REFERENCE_SLIT in table:
        table[DERIVED_SLIT] = round(table[REFERENCE_SLIT] * slit_width(REFERENCE_SLIT) / slit_width(DERIVED_SLIT))
    return table


def summarise(measurements: Sequence[FrameMeasurement]) -> Summary:
    """*the median resolving power per arm and slit over every line of every frame*

    **Key Arguments:**

    - ``measurements`` -- the lines of each frame. A frame with no line adds nothing.

    **Return:**

    - ``summary`` -- ``{arm: {slit: SlitSummary}}``, empty when no frame found a line. ``nightMedians`` is the median R
      of each UTC date, for the comment.

    **Usage:**

    ```python
    summary = summarise(measurements)
    ```
    """
    groups: dict[tuple[str, str], list[FrameMeasurement]] = defaultdict(list)
    for measurement in measurements:
        if measurement.lines:
            groups[(measurement.arm, measurement.slit)].append(measurement)
    summary: Summary = defaultdict(dict)
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


def _describe(slitSummary: SlitSummary) -> str:
    """*the comment text of one slit: line and frame counts and the median R of each night*"""
    nights = ", ".join(f"{power} (MJD {night})" for night, power in slitSummary.nightMedians.items())
    return f"{slitSummary.lineCount} LINES FROM {slitSummary.frameCount} FRAMES; MEDIAN R PER NIGHT {nights}"


def _sky_comments(arcSlits: Mapping[str, SlitSummary], skySlits: Mapping[str, SlitSummary]) -> list[str]:
    """*the comment lines of the sky-line cross-check of one arm, with the ratio to the arc value where there is one*"""
    lines = []
    for slit in sorted(skySlits, key=slit_width):
        sky = skySlits[slit]
        ratio = f", SKY/ARC {sky.resolvingPower / arcSlits[slit].resolvingPower:.2f}" if slit in arcSlits else ""
        lines.append(
            f"# SKY-LINE CROSS-CHECK, NOT IN THE TABLE: {slit} R {sky.resolvingPower} "
            f"({sky.lineCount} LINES FROM {sky.frameCount} FRAMES{ratio})"
        )
    return lines


def format_yaml(summary: Summary, version: str, today: str, skySummary: Summary | None = None) -> str:
    """*the yaml block to paste into ``soxs_detector_parameters.yaml``*

    **Key Arguments:**

    - ``summary`` -- the arc-measured arms and slits
    - ``version`` -- the soxspipe version, for the provenance comment
    - ``today`` -- the measurement date, ISO format
    - ``skySummary`` -- the sky-line measurement, printed as comments only. Default *None*

    **Return:**

    - ``text`` -- one block per arm: provenance comments, then ``spec-res`` keyed by slit, ``SLIT5.0`` derived

    **Usage:**

    ```python
    text = format_yaml(summary, "0.19.1", "2026-11-22")
    ```
    """
    lines: list[str] = []
    for arm, slits in summary.items():
        ordered = sorted(slits, key=slit_width)
        mjds = [mjd for slit in ordered for mjd in slits[slit].mjds]
        comments = common.provenance_comments("SLIT ARCS", mjds, version, today)
        comments += [f"# {slit}: {_describe(slits[slit])}" for slit in ordered]
        table = derive_wide_slit({slit: slits[slit].resolvingPower for slit in ordered})
        if DERIVED_SLIT in table:
            comments.append(
                f"# {DERIVED_SLIT} IS DERIVED, NOT MEASURED: THE {REFERENCE_SLIT} VALUE TIMES "
                f"{slit_width(REFERENCE_SLIT)} / {slit_width(DERIVED_SLIT)}"
            )
        else:
            comments.append(f"# {DERIVED_SLIT} NOT DERIVED: NO {REFERENCE_SLIT} MEASUREMENT")
        comments += _sky_comments(slits, (skySummary or {}).get(arm, {}))
        body = ["spec-res:", *(f"{common.YAML_INDENT}{slit}: {table[slit]}" for slit in sorted(table, key=slit_width))]
        lines += common.arm_block(arm, comments, body)
    return common.join_blocks(lines)


def _skip_derived_slit(path: Path, slit: str) -> bool:
    """*whether a frame is of the derived slit, which is named on stderr and not measured*"""
    if slit != DERIVED_SLIT:
        return False
    print(f"skipped {path}: {DERIVED_SLIT} is derived from {REFERENCE_SLIT}, not measured", file=sys.stderr)
    return True


def collect_arcs(patterns: Sequence[str]) -> list[ArcRef]:
    """*read the slit arcs named by the patterns, naming each frame it leaves out on stderr*

    A frame that is not a slit arc and a ``SLIT5.0`` arc are left out.

    **Key Arguments:**

    - ``patterns`` -- raw frame paths or glob patterns

    **Return:**

    - ``arcs`` -- the arcs to measure

    **Raises:**

    - ``ToolError`` -- if a pattern matches nothing or a file cannot be read

    **Usage:**

    ```python
    arcs = collect_arcs(["raw/*.fits"])
    ```
    """
    arcs = []
    for path in common.expand_paths(patterns):
        arc = read_arc_ref(path)
        if arc is None:
            print(f"skipped {path}: not a slit arc", file=sys.stderr)
        elif not _skip_derived_slit(path, arc.slit):
            arcs.append(arc)
    return arcs


def collect_skies(patterns: Sequence[str]) -> list[SkyRef]:
    """*read the stare frames named by the patterns, naming each frame it leaves out on stderr*

    A frame with no ``_SKYMODEL.fits`` beside it and a ``SLIT5.0`` frame are left out.

    **Key Arguments:**

    - ``patterns`` -- ``_SKYSUB.fits`` paths or glob patterns

    **Return:**

    - ``skies`` -- the frames to measure

    **Raises:**

    - ``ToolError`` -- if a pattern matches nothing or a file cannot be read

    **Usage:**

    ```python
    skies = collect_skies(["stare/*_SKYSUB.fits"])
    ```
    """
    skies = []
    for path in common.expand_paths(patterns):
        sky = read_sky_ref(path)
        if sky is None:
            print(f"skipped {path}: its {SKYMODEL_SUFFIX} file is missing", file=sys.stderr)
        elif not _skip_derived_slit(path, sky.slit):
            skies.append(sky)
    return skies


def _pair_and_report[Frame: ArcRef | SkyRef](
    frames: Sequence[Frame], maps: Sequence[MapRef], maxGapDays: float
) -> list[tuple[Frame, MapRef]]:
    """*pair frames with maps and name each frame left without one on stderr*"""
    pairs, unpaired = pair_frames_with_maps(frames, maps, maxGapDays)
    for frame in unpaired:
        print(f"skipped {frame.path}: no DISP_IMAGE of the {frame.arm} arm within {maxGapDays} days", file=sys.stderr)
    return pairs


def _sky_cross_check(arguments: argparse.Namespace, maps: Sequence[MapRef], arcArms: set[str]) -> Summary | None:
    """*the sky-line summary, or None when no sky input was given or the cross-check failed*

    The cross-check never stops the arc table. A ``ToolError`` on the sky path is named on stderr, and the table prints
    without the cross-check. A sky arm with no arc arm has no table to sit beside, so it is named on stderr too.
    """
    if not arguments.sky_products:
        return None
    try:
        skyPairs = _pair_and_report(collect_skies(arguments.sky_products), maps, arguments.max_gap_days)
        skySummary = summarise(
            [FrameMeasurement(sky.arm, sky.slit, sky.mjd, measure_sky(sky, mapRef)) for sky, mapRef in skyPairs]
        )
    except ToolError as error:
        print(f"skipped the sky-line cross-check: {error}", file=sys.stderr)
        return None
    for arm in sorted(set(skySummary) - arcArms):
        print(f"sky-line cross-check for the {arm} arm not recorded: no arc measurement for that arm", file=sys.stderr)
    return skySummary


def build_parser() -> argparse.ArgumentParser:
    """*the command-line parser of the script*

    **Return:**

    - ``parser`` -- the argument parser

    **Usage:**

    ```python
    arguments = build_parser().parse_args(argv)
    ```
    """
    parser = argparse.ArgumentParser(description="Measure the nominal SPEC_RES of each SOXS arm and slit.")
    parser.add_argument("--arcs", nargs="+", required=True, help="raw slit arc FITS files or glob patterns")
    parser.add_argument("--disp-images", nargs="+", required=True, help="DISP_IMAGE files or glob patterns")
    parser.add_argument(
        "--sky-products", nargs="+", help="_SKYSUB.fits files or patterns, for the sky-line cross-check"
    )
    parser.add_argument("--max-gap-days", type=float, default=DEFAULT_MAX_GAP_DAYS, help="largest frame to map gap")
    parser.add_argument("--saturation", type=float, default=DEFAULT_SATURATION, help="pixels at or above this are bad")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """*measure the lines of every paired slit arc and print the yaml block*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments. Default the process arguments

    **Return:**

    - ``status`` -- 0 when the block was printed

    **Raises:**

    - ``ToolError`` -- if no arc could be paired with a ``DISP_IMAGE`` or no line was measured

    **Usage:**

    ```python
    status = main(["--arcs", "raw/*.fits", "--disp-images", "maps/*.fits"])
    ```
    """
    arguments = build_parser().parse_args(argv)
    maps = [read_map_ref(path) for path in common.expand_paths(arguments.disp_images)]
    arcPairs = _pair_and_report(collect_arcs(arguments.arcs), maps, arguments.max_gap_days)
    if not arcPairs:
        raise ToolError("no slit arc was paired with a DISP_IMAGE")
    summary = summarise(
        [
            FrameMeasurement(arc.arm, arc.slit, arc.mjd, measure_arc(arc, mapRef, arguments.saturation))
            for arc, mapRef in arcPairs
        ]
    )
    if not summary:
        raise ToolError("no arc line was measured on any frame")
    skySummary = _sky_cross_check(arguments, maps, set(summary))
    print(format_yaml(summary, common.soxspipe_version(), today_iso(), skySummary), end="")
    return 0


if __name__ == "__main__":
    sys.exit(common.run_main(main, sys.argv[1:]))
