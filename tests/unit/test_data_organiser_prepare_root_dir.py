"""Unit contracts for `data_organiser.prepare` called through the Python API from outside the workspace (DY-1702)."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import pytest

from soxspipe.commonutils.data_organiser import data_organiser
from tests.factories import harvestable_raw_fits, workspace_organiser

pytestmark = pytest.mark.unit


def _prepare(organiser: data_organiser) -> None:
    try:
        organiser.prepare(report=False)
    finally:
        organiser.close()


def test_prepare_files_frames_from_the_workspace_when_the_cwd_is_empty(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    emptyCwd = tmp_path / "empty"
    emptyCwd.mkdir()
    monkeypatch.chdir(emptyCwd)

    # ACT
    _prepare(organiser)

    # ASSERT
    assert list((root / "raw").rglob("frame1.fits"))
    assert list(emptyCwd.iterdir()) == []


def test_prepare_deletes_a_duplicate_from_the_workspace_and_leaves_the_cwd_copy(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: PREPARE ONCE SO THE FRAME IS KNOWN, THEN DROP A DUPLICATE IN THE WORKSPACE ROOT
    firstOrganiser = workspace_organiser(tmp_path, log=log)
    root = Path(firstOrganiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    _prepare(firstOrganiser)
    filed = next((root / "raw").rglob("frame1.fits"))
    duplicate = root / "frame1.fits"
    shutil.copy(filed, duplicate)
    bystander = tmp_path / "bystander"
    bystander.mkdir()
    shutil.copy(filed, bystander / "frame1.fits")
    monkeypatch.chdir(bystander)

    # ACT
    _prepare(data_organiser(log=log, rootDir=str(root)))

    # ASSERT
    assert not duplicate.exists()
    assert (bystander / "frame1.fits").exists()
    assert filed.exists()
