#!/usr/bin/env python
"""*measure the sky-shift uncertainty of orders that are not tuned to the sky lines, per SOXS arm*

This is a developer tool for DY-2377. It is not part of a reduction and not in the soxspipe command line.

The Phase 3 keyword `SPEC_SYE` is the largest, over orders, of the uncertainty of the sky shift applied to an order
(DY-2307). An order that is not tuned has no measured uncertainty, so each arm carries two fixed values in
`soxs_detector_parameters.yaml`:

- `untuned`: skyline tuning is switched off, or an order failed. The arc solution is used as it is, and the arc may be
  hours from the science frame.
- `attached-arc`: tuning is skipped because an arc is attached to the science frame (7 minutes and less than 5 degrees
  of rotator angle, DY-1665).

**What is measured**

A sky shift is the wavelength correction that tuning applied to an order, so for a frame that was tuned it is also the
error the order would have had without tuning. The script reads the `ESO QC SKY SHIFT O<n>` cards (nm) of reduced
merged spectra (`_EXTRACTED_MERGED` today, `_SPECTRUM` after DY-2346) and takes the root mean square about zero, not
the scatter about the mean, because an untuned order keeps the mean shift as an error. It sorts each tuned frame by its
arc: the frame's `DISP_TAB` is the one named by its `ESO PRO REC1 CAL<i> NAME` card, and the arc is attached when
`attached_arc_offsets` of `ug_wavelength_transfer` says so. The shifts of frames with an unattached arc give
`untuned`, and the shifts of frames with an attached arc give `attached-arc`.

The printed value is the largest over orders of the per-order RMS, because `SPEC_SYE` is a maximum over orders. The RMS
pooled over all orders and frames is printed as a yaml comment, with the per-order values.

The shifts include the noise of the sky-line measurement, so both values are upper limits. A frame that skipped tuning
(`SKY TUNING = skipped: attached arc`) holds no shift. The `attached-arc` value needs frames that were tuned although
their arc is attached, which means products from a build before DY-1665 (soxspipe 9d0d450c) or a run with the skip
removed. The script names each product it leaves out.

**Usage**

```bash
# merged spectra: _EXTRACTED_MERGED today, _SPECTRUM after DY-2346
python tools/measure_sky_shift_uncertainty.py --products "path/to/merged_spectra/*.fits" --arcs "path/to/*_SSOL.fits"
python tools/measure_sky_shift_uncertainty.py --products a.fits b.fits --arcs "arcs/*.fits"
```

`--arcs` are the multi-pinhole `DISP_TAB` products. A product name in the header can be cut short, so a product matches
the `DISP_TAB` whose file name starts with it. It prints the yaml block to paste into
`soxspipe/resources/soxs_detector_parameters.yaml`:

```yaml
VIS:
  # MEASURED 2026-11-22 WITH soxspipe 0.19.1 FROM 40 FRAMES ON 5 NIGHTS (MJD 60000 TO 60040)
  sky-shift-uncertainty-nm:
    untuned: 0.1234
    attached-arc: 0.0123
```

Exit status is 0 when the block was printed and 2 when the input could not be used.
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
from astropy.io import fits
from measurement_common import ToolError, today_iso

UNTUNED = "untuned"
ATTACHED_ARC = "attached-arc"
KINDS = (UNTUNED, ATTACHED_ARC)
YAML_KEY = "sky-shift-uncertainty-nm"
VALUE_DECIMALS = 4
SKY_SHIFT_PATTERN = re.compile(r"^ESO QC SKY SHIFT O(\d+)$")
TUNING_CARD = "ESO QC SKY TUNING"
CALIBRATION_CATG_PATTERN = re.compile(r"^ESO PRO REC1 CAL(\d+) CATG$")
DISP_TAB_CATG_PREFIX = "DISP_TAB"
VALUE_COMMENT = "# VALUE: LARGEST PER-ORDER RMS ABOUT ZERO OF THE SKY SHIFTS (SPEC_SYE IS A MAXIMUM OVER ORDERS)"


@dataclass(frozen=True)
class ShiftSample:
    """*the sky shifts (nm) of one tuned frame, with its arm, MJD and arc kind*"""

    arm: str
    kind: str
    mjd: float
    shifts: Mapping[int, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "shifts", common.frozen_mapping(self.shifts))


@dataclass(frozen=True)
class Uncertainty:
    """*the uncertainty of one arm and kind: the worst-order value, the pooled RMS and the shifts they came from*"""

    valueNm: float
    pooledNm: float
    shiftCount: int
    frameCount: int
    perOrderNm: Mapping[int, float]
    mjds: tuple[float, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "perOrderNm", common.frozen_mapping(self.perOrderNm))
        object.__setattr__(self, "mjds", tuple(self.mjds))


def read_sky_shifts(header: fits.Header) -> dict[int, float]:
    """*the ``SKY SHIFT O<n>`` QCs of a product, in nm, by order; the ``PX`` and ``ERR`` cards are not read*

    **Key Arguments:**

    - ``header`` -- the product's primary header

    **Return:**

    - ``shifts`` -- ``{order: shift in nm}``, finite values only

    **Usage:**

    ```python
    shifts = read_sky_shifts(header)
    ```
    """
    shifts = {}
    for key in header:
        match = SKY_SHIFT_PATTERN.match(key)
        if match is not None and isinstance(header[key], int | float) and math.isfinite(header[key]):
            shifts[int(match.group(1))] = float(header[key])
    return shifts


def rms_about_zero(values: Sequence[float]) -> float:
    """*the root mean square of the values, about zero*

    **Key Arguments:**

    - ``values`` -- the shifts, at least one

    **Return:**

    - ``rms`` -- ``sqrt(mean(value^2))``

    **Raises:**

    - ``ToolError`` -- if there is no shift

    **Usage:**

    ```python
    rms = rms_about_zero([0.1, -0.2])
    ```
    """
    if not values:
        raise ToolError("there is no shift to take the RMS of")
    return math.sqrt(sum(value * value for value in values) / len(values))


def find_arc(header: fits.Header, candidates: Sequence[Path]) -> Path | None:
    """*the ``DISP_TAB`` file that a product's header names, matched by the start of the file name*

    The header card is cut to fit the 80-column card, so the name can be a prefix of the file name.

    **Key Arguments:**

    - ``header`` -- the product header, with ``ESO PRO REC1 CAL<i> CATG`` and ``NAME`` cards
    - ``candidates`` -- the ``DISP_TAB`` files to choose from

    **Return:**

    - ``path`` -- the matching file, or *None* when the header names no ``DISP_TAB`` or no candidate matches

    **Raises:**

    - ``ToolError`` -- if the name matches more than one candidate

    **Usage:**

    ```python
    arcPath = find_arc(header, arcPaths)
    ```
    """
    for key in header:
        match = CALIBRATION_CATG_PATTERN.match(key)
        if match is None or not str(header[key]).startswith(DISP_TAB_CATG_PREFIX):
            continue
        name = str(header.get(f"ESO PRO REC1 CAL{match.group(1)} NAME", ""))
        matches = [path for path in candidates if name and path.name.startswith(name)]
        if len(matches) > 1:
            raise ToolError(f"the arc name {name!r} matches more than one file: {[str(path) for path in matches]}")
        return matches[0] if matches else None
    return None


def classify_frame(scienceHeader: fits.Header, arcHeader: fits.Header) -> str:
    """*whether the arc of a frame is attached, by the rule that skips tuning in the extraction*

    **Key Arguments:**

    - ``scienceHeader`` -- the header of the reduced product, standing for the science frame
    - ``arcHeader`` -- the header of the arc's ``DISP_TAB``

    **Return:**

    - ``kind`` -- ``attached-arc`` when the arc is within the attached-arc limits, otherwise ``untuned``

    **Usage:**

    ```python
    kind = classify_frame(productHeader, arcHeader)
    ```
    """
    common.ensure_repo_on_path()
    from soxspipe.commonutils.ug_wavelength_transfer import attached_arc_offsets, within_attached_arc_limits

    offsets = attached_arc_offsets([scienceHeader], arcHeader, common.keyword_aliases().__getitem__)
    return ATTACHED_ARC if within_attached_arc_limits(offsets) else UNTUNED


def _kind_uncertainty(frames: Sequence[ShiftSample]) -> Uncertainty:
    """*the worst-order and pooled uncertainty of the frames of one arm and kind*"""
    byOrder: dict[int, list[float]] = defaultdict(list)
    for frame in frames:
        for order, shift in frame.shifts.items():
            byOrder[order].append(shift)
    perOrder = {order: rms_about_zero(values) for order, values in sorted(byOrder.items())}
    everyShift = [shift for values in byOrder.values() for shift in values]
    return Uncertainty(
        valueNm=max(perOrder.values()),
        pooledNm=rms_about_zero(everyShift),
        shiftCount=len(everyShift),
        frameCount=len(frames),
        perOrderNm=perOrder,
        mjds=[frame.mjd for frame in frames],
    )


def uncertainty_by_arm(samples: Sequence[ShiftSample]) -> dict[str, dict[str, Uncertainty]]:
    """*the uncertainty of each arm and kind from the sky shifts of the frames*

    **Key Arguments:**

    - ``samples`` -- the tuned frames

    **Return:**

    - ``result`` -- ``{arm: {kind: Uncertainty}}``. The value is the largest per-order RMS about zero, and the pooled
      RMS over all orders and frames is kept beside it. A kind with no frame is absent.

    **Usage:**

    ```python
    result = uncertainty_by_arm(samples)
    ```
    """
    groups: dict[tuple[str, str], list[ShiftSample]] = defaultdict(list)
    for sample in samples:
        groups[(sample.arm, sample.kind)].append(sample)
    result: dict[str, dict[str, Uncertainty]] = defaultdict(dict)
    for (arm, kind), frames in groups.items():
        result[arm][kind] = _kind_uncertainty(frames)
    return dict(result)


def _describe(uncertainty: Uncertainty | None) -> str:
    """*the comment text of one kind: the counts, the pooled RMS and the RMS of each order*"""
    if uncertainty is None:
        return "NO FRAME; NOT MEASURED"
    perOrder = ", ".join(f"O{order} {rms:.{VALUE_DECIMALS}f}" for order, rms in uncertainty.perOrderNm.items())
    return (
        f"{uncertainty.shiftCount} SHIFTS FROM {uncertainty.frameCount} FRAMES; "
        f"POOLED RMS {uncertainty.pooledNm:.{VALUE_DECIMALS}f}; PER ORDER {perOrder} nm"
    )


def format_yaml(result: Mapping[str, Mapping[str, Uncertainty]], version: str, today: str) -> str:
    """*the yaml block to paste into ``soxs_detector_parameters.yaml``*

    **Key Arguments:**

    - ``result`` -- the uncertainties by arm and kind
    - ``version`` -- the soxspipe version, for the provenance comment
    - ``today`` -- the measurement date, ISO format

    **Return:**

    - ``text`` -- one block per arm. A kind with no frame gets a comment and no key.

    **Usage:**

    ```python
    text = format_yaml(result, "0.19.1", "2026-11-22")
    ```
    """
    lines: list[str] = []
    for arm, kinds in result.items():
        mjds = [mjd for uncertainty in kinds.values() for mjd in uncertainty.mjds]
        comments = [*common.provenance_comments("FRAMES", mjds, version, today), VALUE_COMMENT]
        comments += [f"# {kind.upper()}: {_describe(kinds.get(kind))}" for kind in KINDS]
        body = [f"{YAML_KEY}:"]
        body += [
            f"{common.YAML_INDENT}{kind}: {round(kinds[kind].valueNm, VALUE_DECIMALS)}"
            for kind in KINDS
            if kind in kinds
        ]
        lines += common.arm_block(arm, comments, body)
    return common.join_blocks(lines)


def _read_or_skip(path: Path) -> fits.Header | None:
    """*the header of a file, or None after naming the unreadable file on stderr*"""
    try:
        return common.read_header(path)
    except ToolError as error:
        print(f"skipped {path}: {error}", file=sys.stderr)
        return None


def build_samples(productPaths: Sequence[Path], arcPaths: Sequence[Path]) -> list[ShiftSample]:
    """*read the sky shifts of each product and sort the frame by its arc*

    A product that cannot be read, a product with no shift and a product whose named arc is not among ``arcPaths`` are
    named on stderr and left out.

    **Key Arguments:**

    - ``productPaths`` -- the reduced spectra
    - ``arcPaths`` -- the ``DISP_TAB`` files

    **Return:**

    - ``samples`` -- one per tuned frame

    **Raises:**

    - ``ToolError`` -- if no product gives a sky shift, or a card a measured product needs is missing

    **Usage:**

    ```python
    samples = build_samples(productPaths, arcPaths)
    ```
    """
    samples = []
    arcHeaders: dict[Path, fits.Header | None] = {}
    for path in productPaths:
        header = _read_or_skip(path)
        if header is None:
            continue
        shifts = read_sky_shifts(header)
        if not shifts:
            print(f"skipped {path}: no sky shifts (SKY TUNING = {header.get(TUNING_CARD)!r})", file=sys.stderr)
            continue
        arcPath = find_arc(header, arcPaths)
        if arcPath is None:
            print(f"skipped {path}: the DISP_TAB it names is not among the --arcs", file=sys.stderr)
            continue
        if arcPath not in arcHeaders:
            arcHeaders[arcPath] = _read_or_skip(arcPath)
        if arcHeaders[arcPath] is None:
            continue
        kind = classify_frame(header, arcHeaders[arcPath])
        samples.append(ShiftSample(common.arm_of(header, path), kind, common.read_mjd_obs(header, path), shifts))
    if not samples:
        raise ToolError("no sky shift was found in any product")
    return samples


def build_parser() -> argparse.ArgumentParser:
    """*the command-line parser of the script*

    **Return:**

    - ``parser`` -- the argument parser

    **Usage:**

    ```python
    arguments = build_parser().parse_args(argv)
    ```
    """
    parser = argparse.ArgumentParser(
        description="Measure the per-arm sky-shift uncertainty of untuned and attached-arc orders."
    )
    parser.add_argument("--products", nargs="+", required=True, help="merged spectra with SKY SHIFT QCs, or patterns")
    parser.add_argument("--arcs", nargs="+", required=True, help="multi-pinhole DISP_TAB files, or patterns")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """*measure the sky-shift uncertainty of each arm and print the yaml block*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments. Default the process arguments

    **Return:**

    - ``status`` -- 0 when the block was printed

    **Usage:**

    ```python
    status = main(["--products", "merged/*.fits", "--arcs", "arcs/*.fits"])
    ```
    """
    arguments = build_parser().parse_args(argv)
    samples = build_samples(common.expand_paths(arguments.products), common.expand_paths(arguments.arcs))
    print(format_yaml(uncertainty_by_arm(samples), common.soxspipe_version(), today_iso()), end="")
    return 0


if __name__ == "__main__":
    sys.exit(common.run_main(main, sys.argv[1:]))
