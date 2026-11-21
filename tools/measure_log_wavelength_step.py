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
2. Per order and frame, take the steps `c ln(lambda_{i+1} / lambda_i)` between adjacent bins, in km/s. Pool the steps of
   each order over all frames and take their median. The bins of two frames differ, so wavelengths are never pooled.
3. Per arm, take the smallest per-order median and divide it by the spectral binning factor of the frames (the
   `ESO DET BINX` card for VIS, `ESO DET BINY` for NIR). Use unbinned (1x1) frames where you can, because every input
   of an arm must share one binning.

The printed value is per unbinned pixel; the reduction multiplies it by BINX (VIS) or BINY (NIR).

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
import numpy.typing as npt
from astropy.io import fits
from measurement_common import ToolError, today_iso

LIGHT_SPEED_KMS = 299792.458
MIN_ORDER_ROWS = 3
STEP_DECIMALS = 2
YAML_KEY = "log-wavelength-step-kms"
BINNING_COMMENT = "# PER UNBINNED PIXEL; THE REDUCTION MULTIPLIES IT BY BINX (VIS) OR BINY (NIR)"
ORDER_COLUMN = "order"
WAVELENGTH_COLUMN = "wavelengthMean"
EXTRACTED_ORDERS_HDU = 1

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class Product:
    """*the columns and header values of one extracted-orders table that the measurement reads*"""

    path: Path
    arm: str
    binning: int
    mjd: float
    orders: npt.NDArray[np.integer]
    wavelengths: FloatArray


@dataclass(frozen=True)
class ArmResult:
    """*the measured step of one arm, with the per-order values and inputs it came from*"""

    stepKms: float
    perOrderKms: Mapping[int, float]
    binning: int
    mjds: tuple[float, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "perOrderKms", common.frozen_mapping(self.perOrderKms))
        object.__setattr__(self, "mjds", tuple(self.mjds))


def native_steps_kms(wavelengths: npt.ArrayLike) -> FloatArray:
    """*the velocity step between each pair of adjacent wavelength bins, in km/s*

    **Key Arguments:**

    - ``wavelengths`` -- the bin wavelengths of one order of one frame, in any order, with repeats and non-finite
      values allowed

    **Return:**

    - ``steps`` -- ``c ln(lambda_{i+1} / lambda_i)`` over the sorted distinct bins

    **Raises:**

    - ``ToolError`` -- if fewer than two distinct positive wavelengths remain

    **Usage:**

    ```python
    steps = native_steps_kms(orderWavelengths)
    ```
    """
    values = np.asarray(wavelengths, dtype=float)
    distinct = np.unique(values[np.isfinite(values) & (values > 0)])
    if distinct.size < 2:
        raise ToolError("an order needs at least two distinct wavelengths to have a step")
    return LIGHT_SPEED_KMS * np.log(distinct[1:] / distinct[:-1])


def per_order_step_samples(orders: npt.ArrayLike, wavelengths: npt.ArrayLike) -> dict[int, FloatArray]:
    """*the bin-to-bin steps of each order of one frame, in km/s*

    **Key Arguments:**

    - ``orders`` -- the order number of each row
    - ``wavelengths`` -- the bin wavelength of each row

    **Return:**

    - ``samples`` -- ``{order: steps}``. An order with fewer than ``MIN_ORDER_ROWS`` usable rows is left out.

    **Usage:**

    ```python
    samples = per_order_step_samples(product.orders, product.wavelengths)
    ```
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


def read_product(path: Path) -> Product:
    """*read the arm, binning, MJD and the order and wavelength columns of an extracted-orders table*

    **Key Arguments:**

    - ``path`` -- a ``_EXTRACTED_ORDERS`` FITS file: the cards in the primary header, the rows in extension 1

    **Return:**

    - ``product`` -- the values read

    **Raises:**

    - ``ToolError`` -- if the file cannot be read, a card is missing or a column is missing

    **Usage:**

    ```python
    product = read_product(Path("stare_EXTRACTED_ORDERS.fits"))
    ```
    """
    try:
        with fits.open(path) as hdus:
            header = hdus[0].header
            table = hdus[EXTRACTED_ORDERS_HDU].data
            names = table.columns.names
            if ORDER_COLUMN not in names or WAVELENGTH_COLUMN not in names:
                raise ToolError(f"{path} has no {ORDER_COLUMN} and {WAVELENGTH_COLUMN} columns")
            arm = common.arm_of(header, path)
            return Product(
                path=Path(path),
                arm=arm,
                binning=common.spectral_binning(header, arm, path),
                mjd=common.read_mjd_obs(header, path),
                orders=np.array(table[ORDER_COLUMN]),
                wavelengths=np.array(table[WAVELENGTH_COLUMN], dtype=float),
            )
    except (OSError, IndexError) as error:
        raise ToolError(f"cannot read an extracted-orders table from {path}: {error}") from error


def _arm_binning(armProducts: Sequence[Product], arm: str) -> int:
    """*the one spectral binning of an arm's frames, or a ``ToolError`` when they differ*"""
    binnings = {product.binning for product in armProducts}
    if len(binnings) > 1:
        raise ToolError(f"the {arm} frames have different spectral binning {sorted(binnings)}; pass one binning")
    return binnings.pop()


def measure_arms(products: Sequence[Product]) -> dict[str, ArmResult]:
    """*the log-wavelength step of each arm, from the frames of that arm*

    The steps of one order are taken frame by frame and then pooled, because the native bins of two frames differ.

    **Key Arguments:**

    - ``products`` -- the extracted-orders tables, of any arms

    **Return:**

    - ``results`` -- ``{arm: ArmResult}``, with steps per unbinned pixel

    **Raises:**

    - ``ToolError`` -- if the frames of one arm have different binning, or no order of an arm has a step

    **Usage:**

    ```python
    results = measure_arms([read_product(path) for path in paths])
    ```
    """
    results = {}
    for arm in sorted({product.arm for product in products}):
        armProducts = [product for product in products if product.arm == arm]
        binning = _arm_binning(armProducts, arm)
        pooled: dict[int, list[FloatArray]] = {}
        for product in armProducts:
            for order, steps in per_order_step_samples(product.orders, product.wavelengths).items():
                pooled.setdefault(order, []).append(steps)
        if not pooled:
            raise ToolError(f"no order of the {arm} arm has a measured step")
        perOrder = {order: float(np.median(np.concatenate(steps))) / binning for order, steps in pooled.items()}
        results[arm] = ArmResult(min(perOrder.values()), perOrder, binning, [product.mjd for product in armProducts])
    return results


def format_yaml(results: Mapping[str, ArmResult], version: str, today: str) -> str:
    """*the yaml block to paste into ``soxs_detector_parameters.yaml``*

    **Key Arguments:**

    - ``results`` -- the measured arms
    - ``version`` -- the soxspipe version, for the provenance comment
    - ``today`` -- the measurement date, ISO format

    **Return:**

    - ``text`` -- one block per arm, with provenance comments and the ``log-wavelength-step-kms`` key

    **Usage:**

    ```python
    text = format_yaml(results, "0.19.1", "2026-11-22")
    ```
    """
    lines: list[str] = []
    for arm, result in results.items():
        perOrder = ", ".join(f"ORDER {order} {step:.2f}" for order, step in sorted(result.perOrderKms.items()))
        comments = common.provenance_comments("FRAMES", result.mjds, version, today)
        comments += [BINNING_COMMENT, f"# SMALLEST PER-ORDER MEDIAN STEP; {perOrder} km/s"]
        lines += common.arm_block(arm, comments, [f"{YAML_KEY}: {round(result.stepKms, STEP_DECIMALS)}"])
    return common.join_blocks(lines)


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

    **Usage:**

    ```python
    status = main(["path/to/*_EXTRACTED_ORDERS.fits"])
    ```
    """
    arguments = build_parser().parse_args(argv)
    products = [read_product(path) for path in common.expand_paths(arguments.tables)]
    print(format_yaml(measure_arms(products), common.soxspipe_version(), today_iso()), end="")
    return 0


if __name__ == "__main__":
    sys.exit(common.run_main(main, sys.argv[1:]))
