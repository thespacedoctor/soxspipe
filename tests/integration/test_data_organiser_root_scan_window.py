"""Integration contracts for the raw-frame sync when the workspace root holds more frames than one scan window.

The sync reads the workspace root in windows of `ROOT_SCAN_WINDOW` entries. A frame it cannot file stays in the root,
so a later window must move on to the entries not yet examined, and the sync must end (DY-2270). A `raw/<date>/`
folder is created only for a frame that is about to be moved into it (DY-2278).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

from soxspipe.commonutils.data_organiser import data_organiser
from tests.factories import harvestable_raw_fits, workspace_layout, workspace_organiser

pytestmark = pytest.mark.integration

NIGHT = "2024-01-02"
WINDOW = 2
# A HANGING SYNC SCANS THE SAME WINDOW FOR EVER; FAR MORE PASSES THAN ANY TEST HERE NEEDS FAILS IT INSTEAD
PASS_CAP = 20
DATA_ORGANISER_MODULE = sys.modules[data_organiser.__module__]


class _SyncDidNotEndError(AssertionError):
    pass


def _cap_scan_passes(organiser: data_organiser, monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """Fail the sync after `PASS_CAP` scans, and return the list that records each scan's remaining-file count."""
    scan = organiser._create_directory_table
    remainingCounts: list[int] = []

    def capped_scan(**kwargs: Any) -> tuple[Any, Any, Any]:
        if len(remainingCounts) >= PASS_CAP:
            raise _SyncDidNotEndError(f"the sync scanned the workspace root {PASS_CAP} times without ending")
        rawFrames, fitsPaths, remainingFiles = scan(**kwargs)
        remainingCounts.append(remainingFiles)
        return rawFrames, fitsPaths, remainingFiles

    monkeypatch.setattr(organiser, "_create_directory_table", capped_scan)
    return remainingCounts


def _windowed_organiser(tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch) -> tuple[data_organiser, Path]:
    monkeypatch.setattr(DATA_ORGANISER_MODULE, "ROOT_SCAN_WINDOW", WINDOW)
    organiser = workspace_organiser(tmp_path, log=log)
    return organiser, Path(organiser.rootDir)


def _add_held_frame(root: Path, name: str) -> None:
    """Write a root frame whose `raw/<date>/` destination is already taken, so the sync holds it (DY-2262)."""
    harvestable_raw_fits(root / name)
    taken = root / "raw" / NIGHT / name
    taken.parent.mkdir(parents=True, exist_ok=True)
    taken.write_bytes(b"already here")


def _add_stuck_frames(root: Path, names: list[str]) -> None:
    """Write root frames the sync cannot file: alternately one that lacks DPR TYPE and one that is held."""
    for index, name in enumerate(names):
        if index % 2 == 0:
            harvestable_raw_fits(root / name, includeDprType=False)
        else:
            _add_held_frame(root, name)


def _sync(organiser: data_organiser) -> list[str]:
    """Run the raw-frame sync and return the `filepath` of every `raw_frames` row it leaves."""
    organiser.conn, _ = organiser._get_or_create_db_connection()
    organiser._select_instrument()
    try:
        organiser._sync_raw_frames(skipSqlSync=True)
        return [row[0] for row in organiser.conn.execute("SELECT filepath FROM raw_frames").fetchall()]
    finally:
        organiser.close()


def _root_fits_names(root: Path) -> list[str]:
    return sorted(entry.name for entry in root.iterdir() if entry.name.endswith(".fits"))


def _left_in_root_warnings(log: Any) -> list[str]:
    return [
        message for level, message in log.messages if level == "warning" and "left in the workspace root" in message
    ]


def test_sync_ends_when_more_root_frames_than_one_window_are_all_skipped(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: A MISSING-DPR FRAME, A DY-2262 HELD FRAME AND ANOTHER MISSING-DPR FRAME, WITH A WINDOW OF 2
    organiser, root = _windowed_organiser(tmp_path, log, monkeypatch)
    stuck = ["a_stuck.fits", "b_stuck.fits", "c_stuck.fits"]
    _add_stuck_frames(root, stuck)
    remainingCounts = _cap_scan_passes(organiser, monkeypatch)

    # ACT
    _sync(organiser)

    # ASSERT
    assert remainingCounts[0] == 1
    assert _root_fits_names(root) == stuck
    assert (root / "raw" / NIGHT / "b_stuck.fits").read_bytes() == b"already here"


def test_sync_files_the_fileable_frames_that_sort_before_and_after_the_stuck_frames(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: WINDOW 2, THREE STUCK FRAMES BETWEEN TWO FILEABLE FRAMES IN SORTED NAME ORDER
    organiser, root = _windowed_organiser(tmp_path, log, monkeypatch)
    harvestable_raw_fits(root / "a_fileable.fits")
    stuck = ["b_stuck.fits", "c_stuck.fits", "d_stuck.fits"]
    _add_stuck_frames(root, stuck)
    harvestable_raw_fits(root / "e_fileable.fits")
    _cap_scan_passes(organiser, monkeypatch)

    # ACT
    _sync(organiser)

    # ASSERT
    assert _root_fits_names(root) == stuck
    assert (root / "raw" / NIGHT / "a_fileable.fits").is_file()
    assert (root / "raw" / NIGHT / "e_fileable.fits").is_file()


def test_sync_warns_once_with_the_number_of_frames_left_in_the_root(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: FIVE STUCK FRAMES TAKE THREE WINDOWS OF 2
    organiser, root = _windowed_organiser(tmp_path, log, monkeypatch)
    _add_stuck_frames(root, [f"{letter}_stuck.fits" for letter in "abcde"])
    _cap_scan_passes(organiser, monkeypatch)

    # ACT
    _sync(organiser)

    # ASSERT
    [warning] = _left_in_root_warnings(log)
    assert "5 FITS frames" in warning


def test_sync_of_one_window_holding_a_held_frame_warns_only_about_the_hold(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: NO MORE FRAMES THAN ONE WINDOW, SO THE SYNC BEHAVES AS IT ALWAYS HAS
    organiser, root = _windowed_organiser(tmp_path, log, monkeypatch)
    harvestable_raw_fits(root / "a_fileable.fits")
    _add_held_frame(root, "b_held.fits")
    remainingCounts = _cap_scan_passes(organiser, monkeypatch)

    # ACT
    _sync(organiser)

    # ASSERT
    assert remainingCounts == [0]
    assert _left_in_root_warnings(log) == []
    holdWarnings = [message for level, message in log.messages if level == "warning" and "already exists" in message]
    assert len(holdWarnings) == 1
    assert _root_fits_names(root) == ["b_held.fits"]


def test_sync_creates_no_dated_folder_for_a_dated_archive_frame_it_does_not_move(log: Any, tmp_path: Path) -> None:
    # ARRANGE: THE RAW DIRECTORY IS AN EXTERNAL ARCHIVE, AND THE ROOT LINK POINTS AT A FRAME IN ITS DATED FOLDER
    layout = workspace_layout(tmp_path)
    archive = (tmp_path / "archive" / "raw").resolve()
    (archive / NIGHT).mkdir(parents=True)
    frame = harvestable_raw_fits(archive / NIGHT / "x.fits")
    (layout.root / "x.fits").symlink_to(frame)
    organiser = data_organiser(log=log, rootDir=str(layout.root), dbConnect=False)
    organiser.rawDir = str(archive)
    workspaceRawBefore = sorted(path.name for path in (layout.root / "raw").iterdir())

    # ACT
    rowPaths = _sync(organiser)

    # ASSERT
    assert rowPaths == []
    assert sorted(path.name for path in (layout.root / "raw").iterdir()) == workspaceRawBefore
    assert not (layout.root / "raw" / NIGHT).exists()
    assert sorted(path.name for path in archive.iterdir()) == [NIGHT]
    assert frame.is_file()


def test_sync_creates_the_dated_folder_of_a_root_frame_it_moves(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: NO `raw/<date>/` FOLDER EXISTS BEFORE THE SYNC
    organiser, root = _windowed_organiser(tmp_path, log, monkeypatch)
    harvestable_raw_fits(root / "a_fileable.fits")
    assert not (root / "raw" / NIGHT).exists()

    # ACT
    rowPaths = _sync(organiser)

    # ASSERT
    assert (root / "raw" / NIGHT / "a_fileable.fits").is_file()
    assert _root_fits_names(root) == []
    assert rowPaths == [f"./raw/{NIGHT}/a_fileable.fits"]
