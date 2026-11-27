#!/usr/bin/env python
"""
*check an ESO Phase 3 product against the Science Data Products standard (version 8)*

ESO has no offline validator, so the ``phase3_check`` modules are the one place that defines what a compliant
`SCIENCE.SPECTRUM` product is. The tests, the real-data gate and the ``soxspipe check-phase3 <file>...`` command all
call ``check_phase3``, which opens one file, runs every registered check and returns a list of ``Finding`` records.
A ``FAILURE`` is a rule that ESO rejects the file for. A ``WARNING`` is a recommended or optional item that is
missing. ``soxspipe check-phase3`` exits non-zero when any file has a failure.

**Add a check:**

Write a function that takes the ``Phase3Product`` and returns a list of findings, and register it with
``register_check``. Use ``fail`` and ``warn`` to make the findings, and return an empty list when the product is fine.
Put it in a module named ``phase3_check_<topic>.py`` next to this one, and add that module to ``CHECK_MODULES``
below. A module imports ``register_check``, ``fail`` and ``warn`` from ``phase3_check_framework``, not from this
module, which imports the check modules. A check of the table data needs ``needsTable=True``, and a check of its column
cards only needs ``needsTableHeader=True``: the runner then skips it when the file has no such table, and the
structure check reports that fault once. A check that raises is reported as a failure and does not stop the others.

```python
from soxspipe.commonutils.phase3_check_framework import Finding, Phase3Product, fail, register_check, warn


@register_check
def check_exptime_is_positive(product: Phase3Product) -> list[Finding]:
    exptime = product.primary.get("EXPTIME")
    if exptime is not None and exptime <= 0:
        return [fail(f"EXPTIME is {exptime}, it must be greater than 0")]
    return []


@register_check(needsTable=True)
def check_flux_is_not_all_zero(product: Phase3Product) -> list[Finding]:
    if not product.table["FLUX"][0].any():
        return [warn("every FLUX value is 0")]
    return []
```

A keyword that the standard lists as mandatory, recommended or optional goes into ``KEYWORD_RULES`` in
``phase3_check_headers`` instead. The presence check then reports it with the right severity, and the type check reads
the kind of value.

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

from dataclasses import replace
from importlib import import_module
from pathlib import Path

from soxspipe.commonutils import phase3_check_framework
from soxspipe.commonutils.phase3_check_framework import (
    FAILURE,
    WARNING,
    Finding,
    fail,
    is_phase3_name,
    open_product,
)

# THE MODULES THAT HOLD THE CHECKS, IN THE ORDER THEIR CHECKS RUN. A MODULE REGISTERS ITS CHECKS WHEN IT IS IMPORTED
CHECK_MODULES = (
    "soxspipe.commonutils.phase3_check_headers",
    "soxspipe.commonutils.phase3_check_provenance",
    "soxspipe.commonutils.phase3_check_target_flux",
    "soxspipe.commonutils.phase3_check_wavelength",
    "soxspipe.commonutils.phase3_check_file",
    "soxspipe.commonutils.phase3_check_table",
)

# SAID ONCE FOR A FILE THAT IS NOT NAMED AS A PHASE 3 PRODUCT, SO ITS OUTPUT DOES NOT READ AS A BROKEN FILE
NOT_PHASE3_NOTE = "note: not a `_ph3` product; length rules skipped"

for moduleName in CHECK_MODULES:
    import_module(moduleName)


def check_phase3(filePath: str | Path) -> list[Finding]:
    """*run every registered check on one product*

    **Key Arguments:**

    - ``filePath`` -- the path to the product

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, failures and warnings, empty for a compliant product. A file that
      cannot be opened gives one failure

    **Usage:**

    ```python
    from soxspipe.commonutils.phase3_check import FAILURE, check_phase3
    failures = [f for f in check_phase3("SOXS_20260102T030405_STARE_ph3.fits") if f.severity == FAILURE]
    ```
    """
    try:
        product = open_product(filePath)
    except (OSError, ValueError) as error:
        return [replace(fail(f"the file cannot be opened as FITS: {error}"), check="open_file")]

    findings = []
    try:
        for check in phase3_check_framework.CHECKS.values():
            if (check.needsTable and product.table is None) or (check.needsTableHeader and not product.isBinTable):
                continue
            try:
                found = check.function(product)
            except Exception as error:
                found = [fail(f"the check crashed: {type(error).__name__}: {error}")]
            findings.extend(replace(finding, check=check.name) for finding in found)
    finally:
        product.hdus.close()
    return findings


def _count(number: int, noun: str) -> str:
    return f"{number} {noun}" if number == 1 else f"{number} {noun}s"


def check_phase3_command(filePaths: list[str]) -> int:
    """*check each file, print the result of each and return the exit status for the command line*

    **Key Arguments:**

    - ``filePaths`` -- the paths of the products to check

    **Return:**

    - ``status`` -- 1 when any file has a failure, otherwise 0. Warnings do not change it

    **Usage:**

    ```python
    from soxspipe.commonutils.phase3_check import check_phase3_command
    status = check_phase3_command(["SOXS_20260102T030405_STARE_ph3.fits"])
    ```
    """
    status = 0
    for filePath in filePaths:
        findings = check_phase3(filePath)
        failures = [finding for finding in findings if finding.severity == FAILURE]
        warnings = [finding for finding in findings if finding.severity == WARNING]
        verdict = "FAIL" if failures else "PASS"
        print(f"{filePath}: {verdict} ({_count(len(failures), 'failure')}, {_count(len(warnings), 'warning')})")
        if not is_phase3_name(Path(filePath).name):
            print(f"  {NOT_PHASE3_NOTE}")
        for finding in failures + warnings:
            print(f"  {finding.severity} [{finding.check}] {finding.message}")
        if failures:
            status = 1
    return status
