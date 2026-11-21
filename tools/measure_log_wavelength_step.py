#!/usr/bin/env python
"""*measure the log-wavelength step of each SOXS arm from extracted-order tables*

This is a developer tool for DY-2377. It is not part of a reduction and not in the soxspipe command line.

Each arm's merged spectrum sits on one fixed grid that is uniform in log wavelength, `lambda_k = lambda_0 exp(k dv / c)`
(ADR 0003). The step `dv` is a property of the instrument and lives in `soxs_detector_parameters.yaml`. It is the
smallest per-order median native step, multiplied by the spectral binning factor when a frame is reduced. This script
measures that value so it can replace the provisional 18 km/s (VIS) and 11 km/s (NIR).

**Method**

1. Read each `<stem>_EXTRACTED_ORDERS*.fits` table. Its `order` and `wavelengthMean` columns hold one row per native
   bin of an order, and `wavelengthMean` is the bin centre.
2. Per order, take the median of `c ln(lambda_{i+1} / lambda_i)` over adjacent bins, in km/s. Rows of every input file
   are pooled per order, so each order's median is over all frames.
3. Per arm, take the smallest per-order median and divide it by the spectral binning factor of the frames (the
   `ESO DET BINX` card for VIS, `ESO DET BINY` for NIR). The printed value is per unbinned pixel, the way the yaml holds
   it. Use unbinned (1x1) frames where you can, because every input of an arm must share one binning.

**Usage**

```bash
python tools/measure_log_wavelength_step.py "path/to/*_EXTRACTED_ORDERS.fits"
python tools/measure_log_wavelength_step.py a.fits b.fits
```

It prints the yaml block to paste into `soxspipe/resources/soxs_detector_parameters.yaml`:

```yaml
VIS:
  # MEASURED 2026-11-22 WITH soxspipe 0.19.1 FROM 20 FRAMES ON 5 NIGHTS (MJD 60000 TO 60040)
  log-wavelength-step-kms: 17.93
```

Exit status is 0 when the block was printed and 2 when the input could not be used.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import measurement_common as common
import numpy as np
from astropy.io import fits
from measurement_common import ToolError, today_iso

LIGHT_SPEED_KMS = 299792.458
MIN_ORDER_ROWS = 3
STEP_DECIMALS = 2
YAML_KEY = "log-wavelength-step-kms"
ORDER_COLUMN = "order"
WAVELENGTH_COLUMN = "wavelengthMean"
EXTRACTED_ORDERS_HDU = 1


def native_steps_kms(wavelengths) -> np.ndarray:
    """*the velocity step between each pair of adjacent wavelength bins, in km/s*

    **Key Arguments:**

    - ``wavelengths`` -- the bin wavelengths of one order of one frame, in any order, with repeats and non-finite
      values allowed

    **Return:**

    - ``steps`` -- ``c ln(lambda_{i+1} / lambda_i)`` over the sorted distinct bins

    **Raises:**

    - ``ToolError`` -- if fewer than two distinct positive wavelengths remain
    """
    values = np.asarray(wavelengths, dtype=float)
    distinct = np.unique(values[np.isfinite(values) & (values > 0)])
    if distinct.size < 2:
        raise ToolError("an order needs at least two distinct wavelengths to have a step")
    return LIGHT_SPEED_KMS * np.log(distinct[1:] / distinct[:-1])


def native_step_kms(wavelengths) -> float:
    """*the median velocity step between adjacent wavelength bins of one order, in km/s*"""
    return float(np.median(native_steps_kms(wavelengths)))


def per_order_step_samples(orders, wavelengths) -> dict[int, np.ndarray]:
    """*the bin-to-bin steps of each order of one frame, in km/s*

    **Key Arguments:**

    - ``orders`` -- the order number of each row
    - ``wavelengths`` -- the bin wavelength of each row

    **Return:**

    - ``samples`` -- ``{order: steps}``. An order with fewer than ``MIN_ORDER_ROWS`` usable rows is left out.
    """
    orderArray = np.asarray(orders)
    wavelengthArray = np.asarray(wavelengths, dtype=float)
    samples = {}
    for order in np.unique(orderArray):
        rows = wavelengthArray[orderArray == order]
        if np.count_nonzero(np.isfinite(rows) & (rows > 0)) < MIN_ORDER_ROWS:
            continue
        samples[int(order)] = native_steps_kms(rows)
    return samples


def per_order_steps_kms(orders, wavelengths) -> dict[int, float]:
    """*the median native step of each order, in km/s*"""
    return {order: float(np.median(steps)) for order, steps in per_order_step_samples(orders, wavelengths).items()}


def unbinned_step_kms(perOrderSteps: Mapping[int, float], binning: int) -> float:
    """*the log-wavelength step per unbinned pixel: the smallest order step over the binning factor*

    The grid step of a reduction is this value times its binning factor, so the yaml holds the unbinned one.

    **Key Arguments:**

    - ``perOrderSteps`` -- ``{order: median native step in km/s}`` measured on frames with the given binning
    - ``binning`` -- the spectral binning factor of those frames

    **Return:**

    - ``stepKms`` -- the smallest per-order step divided by ``binning``

    **Raises:**

    - ``ToolError`` -- if no order has a step
    """
    if not perOrderSteps:
        raise ToolError("no order has a measured step")
    return min(perOrderSteps.values()) / binning


@dataclass(frozen=True)
class Product:
    """*the columns and header values of one extracted-orders table that the measurement reads*"""

    path: Path
    arm: str
    binning: int
    mjd: float
    orders: np.ndarray
    wavelengths: np.ndarray


@dataclass(frozen=True)
class ArmResult:
    """*the measured step of one arm, with the per-order values and inputs it came from*"""

    stepKms: float
    perOrderKms: dict[int, float]
    binning: int
    mjds: list[float]


def read_product(path: Path) -> Product:
    """*read the arm, binning, MJD and the order and wavelength columns of an extracted-orders table*

    **Key Arguments:**

    - ``path`` -- a ``_EXTRACTED_ORDERS`` FITS file: the cards in the primary header, the rows in extension 1

    **Return:**

    - ``product`` -- the values read

    **Raises:**

    - ``ToolError`` -- if the file cannot be read, a card is missing or a column is missing
    """
    try:
        with fits.open(path) as hdus:
            header = hdus[0].header
            table = hdus[EXTRACTED_ORDERS_HDU].data
            names = table.columns.names
            if ORDER_COLUMN not in names or WAVELENGTH_COLUMN not in names:
                raise ToolError(f"{path} has no {ORDER_COLUMN} and {WAVELENGTH_COLUMN} columns")
            arm = common.arm_of(header)
            return Product(
                path=Path(path),
                arm=arm,
                binning=common.spectral_binning(header, arm),
                mjd=float(header["MJD-OBS"]),
                orders=np.array(table[ORDER_COLUMN]),
                wavelengths=np.array(table[WAVELENGTH_COLUMN], dtype=float),
            )
    except (OSError, IndexError, KeyError) as error:
        raise ToolError(f"cannot read an extracted-orders table from {path}: {error}") from error


def measure_arms(products: Sequence[Product]) -> dict[str, ArmResult]:
    """*the log-wavelength step of each arm, from the frames of that arm*

    The steps of one order are taken frame by frame and then pooled, because the native bins of two frames differ.

    **Key Arguments:**

    - ``products`` -- the extracted-orders tables, of any arms

    **Return:**

    - ``results`` -- ``{arm: ArmResult}``, with steps per unbinned pixel

    **Raises:**

    - ``ToolError`` -- if the frames of one arm have different binning, or no order of an arm has a step
    """
    results = {}
    for arm in sorted({product.arm for product in products}):
        armProducts = [product for product in products if product.arm == arm]
        binnings = {product.binning for product in armProducts}
        if len(binnings) > 1:
            raise ToolError(f"the {arm} frames have different spectral binning {sorted(binnings)}; pass one binning")
        binning = binnings.pop()
        pooled: dict[int, list[np.ndarray]] = {}
        for product in armProducts:
            for order, steps in per_order_step_samples(product.orders, product.wavelengths).items():
                pooled.setdefault(order, []).append(steps)
        perOrder = {order: float(np.median(np.concatenate(steps))) / binning for order, steps in pooled.items()}
        results[arm] = ArmResult(
            stepKms=unbinned_step_kms(perOrder, 1),
            perOrderKms=perOrder,
            binning=binning,
            mjds=[product.mjd for product in armProducts],
        )
    return results


def format_yaml(results: Mapping[str, ArmResult], version: str, today: str) -> str:
    """*the yaml block to paste into ``soxs_detector_parameters.yaml``*

    **Key Arguments:**

    - ``results`` -- the measured arms
    - ``version`` -- the soxspipe version, for the provenance comment
    - ``today`` -- the measurement date, ISO format

    **Return:**

    - ``text`` -- one block per arm, with provenance comments and the ``log-wavelength-step-kms`` key
    """
    lines = []
    for arm, result in results.items():
        perOrder = ", ".join(f"ORDER {order} {step:.2f}" for order, step in sorted(result.perOrderKms.items()))
        lines.append(f"{arm}:")
        lines.extend(f"  {line}" for line in common.provenance_comments("FRAMES", result.mjds, version, today))
        lines.append(f"  # SMALLEST PER-ORDER MEDIAN STEP PER UNBINNED PIXEL; {perOrder} km/s")
        lines.append(f"  {YAML_KEY}: {round(result.stepKms, STEP_DECIMALS)}")
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    """*the command-line parser of the script*"""
    parser = argparse.ArgumentParser(
        description="Measure the log-wavelength step of each SOXS arm from extracted-order tables.",
    )
    parser.add_argument("tables", nargs="+", help="_EXTRACTED_ORDERS.fits files or glob patterns")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """*measure the step of each arm and print the yaml block*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments. Default the process arguments

    **Return:**

    - ``status`` -- 0 when the block was printed
    """
    arguments = build_parser().parse_args(argv)
    products = [read_product(path) for path in common.expand_paths(arguments.tables)]
    results = measure_arms(products)
    print(format_yaml(results, common.soxspipe_version(), today_iso()), end="")
    return 0


if __name__ == "__main__":
    sys.exit(common.run_main(main, sys.argv[1:]))
