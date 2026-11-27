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
merged spectra (`_SPECTRUM.fits`) and takes the root mean square about zero, not
the scatter about the mean, because an untuned order keeps the mean shift as an error. It sorts each tuned frame by its
arc: the frame's `DISP_TAB` is the one named by its `ESO PRO REC1 CAL<i> NAME` card, and the arc is attached when
`attached_arc_offsets` of `ug_wavelength_transfer` says so. The shifts of frames with an unattached arc give
`untuned`, and the shifts of frames with an attached arc give `attached-arc`.

The printed value is the largest over orders of the per-order RMS, because `SPEC_SYE` is a maximum over orders. The RMS
pooled over all orders and frames is printed as a yaml comment, with the per-order values.

An order that was neither measured nor modelled is left out of the RMS, because the pipeline records it as a shift of
exactly 0.000 nm. The script names each order it leaves out on stderr. Only an order with a recorded shift of exactly 0
can be left out, and a non-zero shift always counts.

A product from a build with DY-2461 has an `ESO QC SKY SHIFT O<n> SRC` card for every order that has a shift card. It is
`measured`, `model` or `none`, and `none` marks an order whose 0.000 nm is a placeholder. The script reads that card
for each order that has one, so a measured zero counts. Only an order with no `SRC` card, which means every order of a
product from an older build, goes through the inference below. These cards mark such an order (a measured shift of
exactly zero has none of them, and counts):

- VIS r (order 1): `ESO QC R SHIFT SOURCE = none`. `measured`, `combined` and `predicted` all count.
- VIS u and g (orders 2 and 3): `ESO QC UG SHIFT SOURCE = none`. `model` counts.
- VIS i (order 4) has no source card. The pipeline applies no i shift from fewer than `MIN_SHIFT_LINES` sky lines or a
  sky-line SNR under `MIN_SKY_LINE_SNR`. The card holds the SNR rounded to 0.1, so a card at the limit may come from a
  refused SNR just under it. The script therefore leaves i out when `ESO QC SKY LINE SNR O4` is missing or at most that
  limit, or when the comment of `UG SHIFT SOURCE` says that order 4 kept fewer lines than the limit. A product with no
  `UG SHIFT SOURCE` card (a build before DY-1608) has no such evidence, so its i shift counts.
- Any arm, NIR orders included: a shift of exactly 0 with no `ESO QC SKY SHIFT O<n> PX` card. The pipeline writes the
  PX card whenever it knows the dispersion of the order, and the NIR global fit writes a bare 0 for an order that has
  no correction.

**Known gaps** of the inference for products with no `SRC` card. The script does not leave these zeros out, because no
card marks them:

- The NIR per-order fallback (`SKY TUNING = skylines-per-order`) writes 0.0 with a PX card for an order that matched no
  sky line. The script cannot tell it from a measured zero.
- VIS i is missed when r fails its guards first and i kept fewer than 2 lines with an SNR above the limit. The comment
  of `UG SHIFT SOURCE` then names order 1, not order 4.
- A VIS arm with no r coefficients writes no `R SHIFT SOURCE` card, so an unmeasured r there is not marked either.
- The opposite case: a measured i shift below 0.0005 nm, with a true SNR from 5.40 to 5.449, reads as 0.000 with an
  SNR card of 5.4, so the script leaves it out. Leaving out a near-zero shift only raises the RMS.

The shifts include the noise of the sky-line measurement, so both values are upper limits. A frame that skipped tuning
(`SKY TUNING = skipped: attached arc`) holds no shift. The `attached-arc` value needs frames that were tuned although
their arc is attached, which means products from a build before DY-1665 (soxspipe 9d0d450c) or a run with the skip
removed. The script names each product it leaves out.

**Usage**

```bash
# merged spectra: <SOF stem>_SPECTRUM.fits
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
SKY_SHIFT_CARD_FORMAT = "ESO QC SKY SHIFT O{order}"
SKY_SHIFT_PATTERN = re.compile(r"^" + SKY_SHIFT_CARD_FORMAT.format(order=r"(\d+)") + r"$")
TUNING_CARD = "ESO QC SKY TUNING"
CALIBRATION_CATG_PATTERN = re.compile(r"^ESO PRO REC1 CAL(\d+) CATG$")
DISP_TAB_CATG_PREFIX = "DISP_TAB"
VIS_ARM = "VIS"
R_SOURCE_CARD = "ESO QC R SHIFT SOURCE"
UG_SOURCE_CARD = "ESO QC UG SHIFT SOURCE"
SOURCE_NONE = "none"
SHIFT_SOURCE_PATTERN = re.compile(r"^" + SKY_SHIFT_CARD_FORMAT.format(order=r"(\d+)") + r" SRC$")
# THE REASON THE PIPELINE GIVES IN THE COMMENT OF THE UG SHIFT SOURCE CARD, E.G. "order 4: 0 sky lines (min 10)"
# (``_order_guard_reason`` IN ug_wavelength_transfer.py)
SKY_LINES_REASON_PATTERN = re.compile(r"^order (\d+): (\d+) sky lines")
VALUE_COMMENT = "# VALUE: LARGEST PER-ORDER RMS ABOUT ZERO OF THE SKY SHIFTS (SPEC_SYE IS A MAXIMUM OVER ORDERS)"


def shift_card(order: int) -> str:
    """*the keyword of the ``SKY SHIFT`` QC (nm) of an order*"""
    return SKY_SHIFT_CARD_FORMAT.format(order=order)


def shift_px_card(order: int) -> str:
    """*the keyword of the ``SKY SHIFT ... PX`` QC (px) of an order*"""
    return f"{shift_card(order)} PX"


def shift_source_card(order: int) -> str:
    """*the keyword of the ``SKY SHIFT ... SRC`` QC of an order: ``measured``, ``model`` or ``none`` (DY-2461)*"""
    return f"{shift_card(order)} SRC"


def line_snr_card(order: int) -> str:
    """*the keyword of the ``SKY LINE SNR`` QC of an order*"""
    return f"ESO QC SKY LINE SNR O{order}"


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


def read_shift_sources(header: fits.Header) -> dict[int, tuple[str, str]]:
    """*the ``SKY SHIFT O<n> SRC`` QCs of a product, by order, as ``(source, reason)``*

    **Key Arguments:**

    - ``header`` -- the product's primary header

    **Return:**

    - ``sources`` -- ``{order: (source, reason)}``; empty for a product from a build before DY-2461

    **Usage:**

    ```python
    sources = read_shift_sources(header)
    ```
    """
    sources = {}
    for key in header:
        match = SHIFT_SOURCE_PATTERN.match(key)
        if match is not None and isinstance(header[key], str):
            sources[int(match.group(1))] = (header[key], header.comments[key])
    return sources


def _unmeasured_vis_orders(header: fits.Header) -> dict[int, str]:
    """*the VIS orders that the source and sky-line cards of a product mark as not measured and not modelled*

    The result may name an order whose shift is not 0. ``unmeasured_orders`` drops those.
    """
    common.ensure_repo_on_path()
    from soxspipe.commonutils.ug_wavelength_transfer import (
        MIN_SHIFT_LINES,
        MIN_SKY_LINE_SNR,
        PREDICTED_ORDERS,
        VIS_I_ORDER,
        VIS_R_ORDER,
    )

    leftOut: dict[int, str] = {}
    if header.get(R_SOURCE_CARD) == SOURCE_NONE:
        leftOut[VIS_R_ORDER] = f"{R_SOURCE_CARD} is {SOURCE_NONE}"
    ugSource = header.get(UG_SOURCE_CARD)
    if ugSource == SOURCE_NONE:
        leftOut.update(dict.fromkeys(PREDICTED_ORDERS, f"{UG_SOURCE_CARD} is {SOURCE_NONE}"))
    if ugSource is None:
        return leftOut
    snr = header.get(line_snr_card(VIS_I_ORDER))
    if not isinstance(snr, int | float) or not math.isfinite(snr):
        leftOut[VIS_I_ORDER] = "no finite sky-line SNR card, so no sky shift was measured"
    elif snr <= MIN_SKY_LINE_SNR:
        # THE CARD IS ROUNDED TO 0.1, SO A CARD AT THE LIMIT CAN HOLD AN SNR THE PIPELINE REFUSED
        leftOut[VIS_I_ORDER] = f"sky-line SNR {snr:.1f} is not above {MIN_SKY_LINE_SNR}"
    else:
        match = SKY_LINES_REASON_PATTERN.match(header.comments[UG_SOURCE_CARD])
        if match and int(match.group(1)) == VIS_I_ORDER and int(match.group(2)) < MIN_SHIFT_LINES:
            leftOut[VIS_I_ORDER] = f"{match.group(2)} sky lines kept (min {MIN_SHIFT_LINES})"
    return leftOut


def unmeasured_orders(header: fits.Header, arm: str, shifts: Mapping[int, float]) -> dict[int, str]:
    """*the orders of a product whose shift is a recorded zero, not a measured or modelled shift*

    A measured or modelled shift of exactly 0.000 nm has none of the marks, so it is not named. An order with a
    ``SKY SHIFT O<n> SRC`` card is decided by that card alone. See the module docstring for the cards of the older
    inference.

    **Key Arguments:**

    - ``header`` -- the product's primary header
    - ``arm`` -- the arm of the product, ``VIS`` or ``NIR`` (``common.arm_of``)
    - ``shifts`` -- ``{order: shift in nm}`` of the product (``read_sky_shifts``)

    **Return:**

    - ``leftOut`` -- ``{order: reason}`` for each order of ``shifts`` that is not a measurement

    **Usage:**

    ```python
    leftOut = unmeasured_orders(header, "VIS", read_sky_shifts(header))
    ```
    """
    sources = read_shift_sources(header)
    leftOut = {
        order: f"{shift_source_card(order)} is {source} ({reason})"
        for order, (source, reason) in sources.items()
        if source == SOURCE_NONE
    }
    unsourced = {order: shift for order, shift in shifts.items() if order not in sources}
    inferred = _unmeasured_vis_orders(header) if arm == VIS_ARM else {}
    leftOut.update({order: reason for order, reason in inferred.items() if order not in sources})
    for order in unsourced:
        if shift_px_card(order) not in header:
            leftOut.setdefault(order, f"a shift of 0 and no {shift_px_card(order)} card")
    # ONLY A RECORDED ZERO CAN BE A NON-MEASUREMENT: A NON-ZERO SHIFT IS ALWAYS A MEASUREMENT
    return {order: reason for order, reason in sorted(leftOut.items()) if shifts.get(order) == 0}


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


def _measured_shifts(path: Path, header: fits.Header, arm: str, shifts: Mapping[int, float]) -> dict[int, float]:
    """*the shifts of a product without its unmeasured orders, each of which is named on stderr*"""
    common.ensure_repo_on_path()
    from soxspipe.commonutils.toolkit import order_label

    leftOut = unmeasured_orders(header, arm, shifts)
    for order, reason in leftOut.items():
        print(f"left out order {order_label(order, arm)} of {path}: {reason}", file=sys.stderr)
    return {order: shift for order, shift in shifts.items() if order not in leftOut}


def build_samples(productPaths: Sequence[Path], arcPaths: Sequence[Path]) -> list[ShiftSample]:
    """*read the sky shifts of each product and sort the frame by its arc*

    A product that cannot be read, a product with no shift and a product whose named arc is not among ``arcPaths`` are
    named on stderr and left out. So is each order that is not a measurement (``unmeasured_orders``), and a product
    with no measured order.

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
        arm = common.arm_of(header, path)
        shifts = _measured_shifts(path, header, arm, shifts)
        if not shifts:
            print(f"skipped {path}: no sky shift was measured", file=sys.stderr)
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
        samples.append(ShiftSample(arm, kind, common.read_mjd_obs(header, path), shifts))
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
