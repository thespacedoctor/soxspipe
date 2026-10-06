"""Unit contracts for workspace paths used by the data organiser."""

import pytest

from soxspipe.commonutils.data_organiser import data_organiser

pytestmark = pytest.mark.unit


def test_tilde_root_dir_is_expanded(log, tmp_path, monkeypatch) -> None:
    """Expand the documented ``~/...`` workspace form before deriving paths."""
    homePath = tmp_path / "custom-home"
    monkeypatch.setenv("HOME", str(homePath))
    (homePath / "dy598-workspace").mkdir(parents=True)

    organiser = data_organiser(
        log=log,
        rootDir="~/dy598-workspace",
        dbConnect=False,
    )

    expectedRoot = homePath / "dy598-workspace"
    assert organiser.rootDir == str(expectedRoot)
    assert organiser.rawDir == str(expectedRoot / "raw")
    assert organiser.miscDir == str(expectedRoot / "misc")
    assert organiser.sessionsDir == str(expectedRoot / "sessions")
    assert all(
        "~" not in path
        for path in (
            organiser.rootDir,
            organiser.rawDir,
            organiser.miscDir,
            organiser.sessionsDir,
        )
    )


def test_absolute_root_dir_is_unchanged(log, tmp_path) -> None:
    """Preserve an absolute workspace path and its derived directories."""
    rootPath = tmp_path / "absolute-workspace"
    rootPath.mkdir()

    organiser = data_organiser(
        log=log,
        rootDir=str(rootPath),
        dbConnect=False,
    )

    assert organiser.rootDir == str(rootPath)
    assert organiser.rawDir == str(rootPath / "raw")
    assert organiser.miscDir == str(rootPath / "misc")
    assert organiser.sessionsDir == str(rootPath / "sessions")
