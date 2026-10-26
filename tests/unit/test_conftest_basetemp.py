"""*tests for the per-run default basetemp that keeps parallel runs apart*"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from tests.conftest import default_basetemp, remove_stale_basetemps

pytestmark = pytest.mark.unit

ROOT = Path("/work/checkout-one")


@pytest.fixture
def system_temp(tmp_path, monkeypatch) -> Path:
    """*point the system temporary directory at an empty directory*

    **Return:**

    - ``systemTemp`` -- the empty directory
    """
    systemTemp = tmp_path / "systmp"
    systemTemp.mkdir()
    monkeypatch.setattr(tempfile, "gettempdir", lambda: str(systemTemp))

    return systemTemp


DEAD_PID = 424242
LIVE_PID = 424243


@pytest.fixture(autouse=True)
def fake_process_table(monkeypatch):
    """*make `os.kill` report DEAD_PID as gone and every other process as alive*"""

    def kill(processId, signal):
        if processId == DEAD_PID:
            raise ProcessLookupError

    monkeypatch.setattr(os, "kill", kill)


def test_two_checkouts_get_different_basetemps():
    first = default_basetemp(Path("/work/checkout-one"), 100)
    second = default_basetemp(Path("/work/checkout-two"), 100)

    assert first != second


def test_two_runs_of_one_checkout_get_different_basetemps():
    assert default_basetemp(ROOT, 100) != default_basetemp(ROOT, 101)


def test_one_run_of_one_checkout_always_gets_the_same_basetemp():
    assert default_basetemp(ROOT, 100) == default_basetemp(ROOT, 100)


def test_the_basetemp_defaults_to_the_current_process():
    assert default_basetemp(ROOT).name.endswith(f"-{os.getpid()}")


def test_the_basetemp_sits_in_the_system_temp_directory_not_in_the_checkout():
    basetemp = default_basetemp(ROOT, 100)

    assert basetemp.parent == Path(tempfile.gettempdir())
    assert not basetemp.is_relative_to(ROOT)
    assert basetemp.name.startswith("soxspipe-pytest-")


def test_a_basetemp_of_a_dead_run_in_this_checkout_is_removed(system_temp):
    stale = default_basetemp(ROOT, DEAD_PID)
    (stale / "test_x0").mkdir(parents=True)

    remove_stale_basetemps(ROOT)

    assert not stale.exists()


def test_a_basetemp_of_a_live_run_in_this_checkout_is_kept(system_temp):
    live = default_basetemp(ROOT, LIVE_PID)
    live.mkdir()

    remove_stale_basetemps(ROOT)

    assert live.exists()


def test_the_basetemp_of_this_very_run_is_kept(system_temp):
    own = default_basetemp(ROOT, os.getpid())
    own.mkdir()

    remove_stale_basetemps(ROOT)

    assert own.exists()


def test_a_folder_whose_name_does_not_parse_is_kept(system_temp):
    unparsed = [
        default_basetemp(ROOT, 1).with_name(default_basetemp(ROOT, 1).name.rsplit("-", 1)[0] + "-notapid"),
        default_basetemp(ROOT, 1).with_name(default_basetemp(ROOT, 1).name.rsplit("-", 1)[0] + "-12x"),
        default_basetemp(ROOT, 1).with_name(default_basetemp(ROOT, 1).name.rsplit("-", 1)[0]),
    ]
    for folder in unparsed:
        folder.mkdir()

    remove_stale_basetemps(ROOT)

    assert all(folder.exists() for folder in unparsed)


def test_a_dead_runs_folder_of_another_checkout_is_kept(system_temp):
    other = default_basetemp(Path("/work/checkout-two"), DEAD_PID)
    other.mkdir()

    remove_stale_basetemps(ROOT)

    assert other.exists()


def test_a_file_with_a_stale_looking_name_is_kept(system_temp):
    lookalike = default_basetemp(ROOT, DEAD_PID)
    lookalike.write_text("not a directory", encoding="utf-8")

    remove_stale_basetemps(ROOT)

    assert lookalike.exists()


def test_a_process_owned_by_someone_else_counts_as_alive(system_temp, monkeypatch):
    foreign = default_basetemp(ROOT, 424242)
    foreign.mkdir()

    def denied(processId, signal):
        raise PermissionError

    monkeypatch.setattr(os, "kill", denied)

    remove_stale_basetemps(ROOT)

    assert foreign.exists()


@pytest.mark.parametrize("suffix", ["9" * 30, str(2**63)])
def test_a_pid_beyond_the_platform_range_keeps_the_folder_instead_of_aborting(system_temp, monkeypatch, suffix):
    huge = default_basetemp(ROOT, 1).with_name(default_basetemp(ROOT, 1).name.rsplit("-", 1)[0] + "-" + suffix)
    huge.mkdir()

    def overflow(processId, signal):
        raise OverflowError("Python int too large to convert to C int")

    monkeypatch.setattr(os, "kill", overflow)

    remove_stale_basetemps(ROOT)

    assert huge.exists()


def test_a_value_error_from_the_process_check_keeps_the_folder(system_temp, monkeypatch):
    folder = default_basetemp(ROOT, DEAD_PID)
    folder.mkdir()

    def invalid(processId, signal):
        raise ValueError("bad pid")

    monkeypatch.setattr(os, "kill", invalid)

    remove_stale_basetemps(ROOT)

    assert folder.exists()
