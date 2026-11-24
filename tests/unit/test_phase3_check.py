"""Unit tests for the check registry and runner of `soxspipe.commonutils.phase3_check` (DY-2345)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils import phase3_check, phase3_check_framework
from soxspipe.commonutils.phase3_check import FAILURE, WARNING, check_phase3
from tests.factories.phase3_products import (
    PRODUCT_NAME,
    SPECTRUM_LENGTH,
    phase3_spectrum_file,
)

pytestmark = pytest.mark.unit


def test_a_file_that_is_not_fits_gives_one_failure(tmp_path: Path) -> None:
    # ARRANGE
    notFits = tmp_path / PRODUCT_NAME
    notFits.write_text("this is not a FITS file")

    # ACT
    findings = check_phase3(notFits)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]
    assert "cannot be opened" in findings[0].message


@pytest.mark.parametrize("name", ["", "."])
def test_a_path_that_is_not_a_file_gives_one_failure(
    tmp_path: Path, name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    monkeypatch.chdir(tmp_path)

    # ACT
    findings = check_phase3(name)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]


def test_a_truncated_fits_file_gives_one_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    product.write_bytes(product.read_bytes()[:100])

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]


def test_a_missing_file_gives_one_failure(tmp_path: Path) -> None:
    # ARRANGE
    missing = tmp_path / PRODUCT_NAME

    # ACT
    findings = check_phase3(missing)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]


@pytest.fixture
def isolated_registry(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Give a test its own copy of the check registry, so a check it registers is gone afterwards."""
    registry = dict(phase3_check_framework.CHECKS)
    monkeypatch.setattr(phase3_check_framework, "CHECKS", registry)
    return registry


def test_a_registered_check_runs_and_its_finding_carries_the_check_name(
    tmp_path: Path, log: object, isolated_registry: dict
) -> None:
    # ARRANGE
    @phase3_check.register_check
    def check_object_is_not_a_test_target(product: phase3_check.Phase3Product) -> list:
        return [phase3_check.fail(f"{product.primary['OBJECT']} is a test target")]

    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert findings == [
        phase3_check.Finding(FAILURE, "SYNTHETIC TARGET is a test target", "check_object_is_not_a_test_target")
    ]


def test_a_registered_check_may_return_a_warning(tmp_path: Path, log: object, isolated_registry: dict) -> None:
    # ARRANGE
    @phase3_check.register_check
    def check_for_a_doubt(product: phase3_check.Phase3Product) -> list:
        return [phase3_check.warn("doubtful")]

    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert [(finding.severity, finding.message) for finding in findings] == [(WARNING, "doubtful")]


def test_registering_two_checks_with_one_name_is_an_error(isolated_registry: dict) -> None:
    # ARRANGE
    def check_something(product: phase3_check.Phase3Product) -> list:
        return []

    phase3_check.register_check(check_something)

    # ACT, ASSERT
    with pytest.raises(ValueError, match="check_something"):
        phase3_check.register_check(check_something)


def test_a_check_that_raises_is_a_failure_and_the_other_checks_still_run(
    tmp_path: Path, log: object, isolated_registry: dict
) -> None:
    # ARRANGE
    @phase3_check.register_check
    def check_that_crashes(product: phase3_check.Phase3Product) -> list:
        raise KeyError("EXPTIME")

    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"PRODCATG": None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    crash = [finding for finding in findings if finding.check == "check_that_crashes"]
    assert len(crash) == 1
    assert crash[0].severity == FAILURE
    assert "KeyError" in crash[0].message
    assert any("PRODCATG" in finding.message for finding in findings)


def test_a_table_check_is_skipped_when_the_file_has_no_table(tmp_path: Path, isolated_registry: dict) -> None:
    # ARRANGE
    calls: list[str] = []

    @phase3_check.register_check(needsTable=True)
    def check_that_reads_the_table(product: phase3_check.Phase3Product) -> list:
        calls.append("ran")
        return []

    emptyFile = tmp_path / PRODUCT_NAME
    fits.HDUList([fits.PrimaryHDU()]).writeto(emptyFile)

    # ACT
    findings = check_phase3(emptyFile)

    # ASSERT
    assert calls == []
    assert any("exactly one extension" in finding.message for finding in findings)


@pytest.mark.parametrize("extension", [None, "image"])
def test_no_check_crashes_on_a_file_without_a_binary_table(tmp_path: Path, extension: str | None) -> None:
    # ARRANGE
    product = tmp_path / PRODUCT_NAME
    hdus = [fits.PrimaryHDU()] + ([fits.ImageHDU(np.zeros(3), name="SPECTRUM")] if extension else [])
    fits.HDUList(hdus).writeto(product)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert not [finding for finding in findings if "crashed" in finding.message]
    assert not [finding for finding in findings if finding.check in ("check_column_names", "check_column_descriptors")]


def test_the_module_docstring_examples_register_and_run(tmp_path: Path, log: object, isolated_registry: dict) -> None:
    # ARRANGE
    namespace: dict = {}
    exampleStart = phase3_check.__doc__.index("```python") + len("```python")
    exec(phase3_check.__doc__[exampleStart : phase3_check.__doc__.index("```", exampleStart)], namespace)  # noqa: S102
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, columns={"FLUX": np.zeros(SPECTRUM_LENGTH)})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert [(finding.check, finding.severity) for finding in findings] == [("check_flux_is_not_all_zero", WARNING)]
