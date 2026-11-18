"""Unit contracts for `data_organiser.prepare` called through the Python API from outside the workspace (DY-1702)."""

from __future__ import annotations

import shutil
import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd
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


def _link_frame_into_workspace(root: Path, target: Path) -> Path:
    link = root / "linked.fits"
    link.symlink_to(target)
    return link


def test_prepare_never_moves_a_linked_frame_that_lies_outside_the_workspace(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    outside = tmp_path / "outside"
    outside.mkdir()
    victim = harvestable_raw_fits(outside / "victim.fits")
    _link_frame_into_workspace(root, victim)
    foreignCwd = tmp_path / "foreign"
    foreignCwd.mkdir()
    monkeypatch.chdir(foreignCwd)

    # ACT
    _prepare(organiser)

    # ASSERT
    assert victim.is_file()
    assert not list((root / "raw").rglob("victim.fits"))


def test_sync_never_deletes_the_target_of_a_linked_duplicate_outside_the_workspace(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: AN INDEX FROM AN EARLIER RELEASE HOLDS A ROW WHOSE NAME IS THE `..` PATH OF A LINK TARGET
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    outside = tmp_path / "outside"
    outside.mkdir()
    victim = harvestable_raw_fits(outside / "victim.fits")
    _link_frame_into_workspace(root, victim)
    organiser.conn = sqlite3.connect(":memory:")
    linkedName = "./../outside/victim.fits"
    row = {"file": linkedName, "eso dpr tech": "ECHELLE,SLIT,STARE", "mjd-date": "2024-01-02"}
    pd.DataFrame([{**row, "filepath": f"./raw/2024-01-02/{linkedName}"}]).to_sql(
        "raw_frames", organiser.conn, index=False
    )
    scans = iter([(pd.DataFrame([{**row, "filepath": None}]), [linkedName], 1), (None, None, 0)])
    monkeypatch.setattr(organiser, "_create_directory_table", lambda **_: next(scans))
    monkeypatch.setattr(organiser, "_populate_raw_frames_extra_columns", lambda table: table.copy())

    # ACT
    organiser._sync_raw_frames(skipSqlSync=True)

    # ASSERT
    assert victim.is_file()


def test_prepare_files_the_good_frame_and_moves_an_empty_fits_file_to_misc(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    (root / "empty.fits").write_bytes(b"")
    foreignCwd = tmp_path / "foreign"
    foreignCwd.mkdir()
    monkeypatch.chdir(foreignCwd)

    # ACT
    _prepare(organiser)

    # ASSERT
    assert list((root / "raw").rglob("frame1.fits"))
    assert (root / "misc" / "empty.fits").exists()
