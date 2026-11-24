"""Unit tests for the `soxspipe check-phase3` command and its usage text (DY-2345)."""

from __future__ import annotations

import os
import subprocess
import sys
from importlib import import_module
from pathlib import Path

import pytest
from docopt import docopt as parse_docopt

from soxspipe.commonutils import phase3_check
from tests.factories.phase3_products import (
    PRODUCT_NAME,
    phase3_spectrum_file,
)

pytestmark = pytest.mark.unit


def test_the_command_prints_each_file_and_returns_0_when_there_is_no_failure(
    tmp_path: Path, log: object, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"LAMRMS": None})

    # ACT
    status = phase3_check.check_phase3_command([str(product)])

    # ASSERT
    output = capsys.readouterr().out
    assert status == 0
    assert f"{product}: PASS (0 failures, 1 warning)" in output
    assert "warning [check_mandatory_and_recommended_keywords] optional keyword LAMRMS is missing" in output


def test_the_command_returns_1_when_any_file_has_a_failure(
    tmp_path: Path, log: object, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    good = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    bad = phase3_spectrum_file(tmp_path / "SOXS_20260102T030406_STARE_ph3.fits", log, primary={"PRODCATG": None})

    # ACT
    status = phase3_check.check_phase3_command([str(good), str(bad)])

    # ASSERT
    output = capsys.readouterr().out
    assert status == 1
    assert f"{good}: PASS (0 failures, 0 warnings)" in output
    assert f"{bad}: FAIL (1 failure, 0 warnings)" in output
    assert "failure [check_mandatory_and_recommended_keywords] mandatory keyword PRODCATG is missing" in output


def test_the_command_returns_1_for_a_file_that_does_not_exist(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE, ACT
    status = phase3_check.check_phase3_command([str(tmp_path / "missing.fits")])

    # ASSERT
    assert status == 1
    assert "FAIL" in capsys.readouterr().out


def _run_soxspipe(monkeypatch: pytest.MonkeyPatch, *arguments: str) -> int:
    clUtils = import_module("soxspipe.cl_utils")
    monkeypatch.setattr(sys, "argv", ["soxspipe", *arguments])
    with pytest.raises(SystemExit) as exit_info:
        clUtils.main()
    return exit_info.value.code


def test_soxspipe_check_phase3_exits_0_for_a_good_product_and_writes_nothing(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    products = tmp_path / "products"
    products.mkdir()
    product = phase3_spectrum_file(products / PRODUCT_NAME, log)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    monkeypatch.chdir(workspace)

    # ACT
    code = _run_soxspipe(monkeypatch, "check-phase3", str(product))

    # ASSERT
    assert code == 0
    assert "PASS" in capsys.readouterr().out
    assert list(workspace.iterdir()) == []


def test_soxspipe_check_phase3_exits_1_when_any_of_several_files_fails(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    good = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    bad = phase3_spectrum_file(tmp_path / "SOXS_20260102T030406_STARE_ph3.fits", log, primary={"PRODCATG": None})
    monkeypatch.chdir(tmp_path)

    # ACT
    code = _run_soxspipe(monkeypatch, "check-phase3", str(good), str(bad))

    # ASSERT
    assert code == 1
    output = capsys.readouterr().out
    assert "PASS" in output
    assert "FAIL" in output


def test_the_usage_text_that_fundamentals_reads_has_only_names_it_can_unpack() -> None:
    # ARRANGE
    clUtils = import_module("soxspipe.cl_utils")

    # ACT
    arguments = parse_docopt(clUtils.TOOLS_DOC, argv=["prep"])

    # ASSERT
    unpackedNames = [
        argument.replace("-", "") + "Flag" if argument[0] == "-" else argument.replace("<", "").replace(">", "")
        for argument in arguments
    ]
    assert all(name.isidentifier() for name in unpackedNames)
    assert "check-phase3" not in arguments


def test_the_full_usage_text_still_has_the_check_phase3_command() -> None:
    # ARRANGE
    clUtils = import_module("soxspipe.cl_utils")

    # ACT
    arguments = parse_docopt(clUtils.__doc__, argv=["check-phase3", "a.fits", "b.fits"])

    # ASSERT
    assert arguments["check-phase3"] is True
    assert arguments["<file>"] == ["a.fits", "b.fits"]


def test_the_real_command_line_exits_1_for_a_bad_product_and_0_for_a_good_one(tmp_path: Path, log: object) -> None:
    # ARRANGE
    good = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    bad = phase3_spectrum_file(tmp_path / "SOXS_20260102T030406_STARE_ph3.fits", log, primary={"PRODCATG": None})
    driver = "import sys; sys.argv = ['soxspipe', *sys.argv[1:]]; from soxspipe.cl_utils import main; main()"
    packageRoot = str(Path(phase3_check.__file__).resolve().parents[2])

    def run(*files: Path) -> subprocess.CompletedProcess:
        return subprocess.run(  # noqa: S603
            [sys.executable, "-W", "ignore", "-c", driver, "check-phase3", *map(str, files)],
            cwd=tmp_path,
            env={**os.environ, "PYTHONPATH": packageRoot},
            capture_output=True,
            text=True,
            check=False,
        )

    # ACT
    goodRun, badRun = run(good), run(good, bad)

    # ASSERT
    assert (goodRun.returncode, badRun.returncode) == (0, 1)
    assert "PASS" in goodRun.stdout
    assert "mandatory keyword PRODCATG is missing" in badRun.stdout
    assert goodRun.stderr == ""
