"""Unit contracts for joining a path under a workspace root (DY-1702)."""

from __future__ import annotations

from pathlib import Path

import pytest

from soxspipe.commonutils.workspace_paths import path_under_root

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("rootDir", [".", "./", ""])
def test_path_is_unchanged_when_the_root_is_the_current_directory(rootDir: str) -> None:
    assert path_under_root("./raw/2024-01-01/a.fits", rootDir) == "./raw/2024-01-01/a.fits"


def test_relative_path_is_joined_under_the_root_and_tidied() -> None:
    assert path_under_root("./sessions/s1/reduced", "/data/workspace") == "/data/workspace/sessions/s1/reduced"


def test_absolute_path_is_unchanged_whatever_the_root() -> None:
    assert path_under_root("/elsewhere/a.fits", "/data/workspace") == "/elsewhere/a.fits"


def test_leading_tilde_of_the_root_is_expanded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))

    assert path_under_root("raw/a.fits", "~/workspace") == str(tmp_path / "workspace" / "raw" / "a.fits")
