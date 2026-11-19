"""Unit contracts for filing a root frame without replacing a file that already sits at its destination (DY-2262)."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.commonutils.data_organiser import data_organiser
from tests.factories import harvestable_raw_fits, workspace_organiser

pytestmark = pytest.mark.unit

NIGHT = "2024-01-02"


def _frames_table(*names: str) -> pd.DataFrame:
    return pd.DataFrame({"file": list(names), "filepath": [f"./raw/{NIGHT}/{name}" for name in names]})


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_a_frame_whose_destination_is_taken_is_dropped_with_a_warning_naming_both_paths(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    harvestable_raw_fits(root / "taken.fits")
    harvestable_raw_fits(root / "free.fits")
    archived = root / "raw" / NIGHT / "taken.fits"
    archived.parent.mkdir(parents=True)
    archived.write_bytes(b"archive")

    # ACT
    kept = organiser._without_overwriting_frames(_frames_table("taken.fits", "free.fits"))

    # ASSERT
    assert kept["file"].tolist() == ["free.fits"]
    assert archived.read_bytes() == b"archive"
    assert (root / "taken.fits").is_file()
    [warning] = _warnings(log)
    assert str(root / "taken.fits") in warning
    assert str(archived.resolve()) in warning


def test_a_frame_whose_destination_is_a_dangling_link_is_dropped(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    dangling = root / "raw" / NIGHT / "frame1.fits"
    dangling.parent.mkdir(parents=True)
    dangling.symlink_to(tmp_path / "missing.fits")

    # ACT
    kept = organiser._without_overwriting_frames(_frames_table("frame1.fits"))

    # ASSERT
    assert kept.empty
    assert not (tmp_path / "missing.fits").exists()


def test_a_link_to_a_frame_already_in_its_dated_folder_is_kept(log: Any, tmp_path: Path) -> None:
    # ARRANGE: THE ARCHIVE LINK THAT `prep` MAKES IN THE ROOT TO INDEX AN ARCHIVE FRAME IS ALREADY IN PLACE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    (root / "raw" / NIGHT).mkdir(parents=True)
    archived = harvestable_raw_fits(root / "raw" / NIGHT / "archived.fits")
    (root / "archived.fits").symlink_to(archived)

    # ACT
    kept = organiser._without_overwriting_frames(_frames_table("archived.fits"))

    # ASSERT
    assert kept["file"].tolist() == ["archived.fits"]
    assert _warnings(log) == []


def test_a_known_frame_is_removed_by_the_duplicate_check_before_the_move(log: Any, tmp_path: Path) -> None:
    # ARRANGE: PREPARE ONCE SO THE FRAME IS INDEXED, THEN DROP A COPY OF IT IN THE WORKSPACE ROOT
    firstOrganiser = workspace_organiser(tmp_path, log=log)
    root = Path(firstOrganiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    firstOrganiser.prepare(report=False)
    firstOrganiser.close()
    filed = next((root / "raw").rglob("frame1.fits"))
    filedBytes = filed.read_bytes()
    shutil.copy(filed, root / "frame1.fits")

    # ACT
    secondOrganiser = data_organiser(log=log, rootDir=str(root))
    secondOrganiser.prepare(report=False)
    secondOrganiser.close()

    # ASSERT
    assert filed.read_bytes() == filedBytes
    assert not (root / "frame1.fits").exists()
    assert _warnings(log) == []


def _indexed_files(organiser: data_organiser) -> list[str]:
    return [row[0] for row in organiser.conn.execute("SELECT file FROM raw_frames")]


def test_prepare_leaves_both_frames_in_place_and_does_not_index_the_skipped_frame(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    rootFrame = harvestable_raw_fits(root / "frame1.fits")
    rootBytes = rootFrame.read_bytes()
    archived = root / "raw" / NIGHT / "frame1.fits"
    archived.parent.mkdir(parents=True)
    archived.write_bytes(b"archive")
    (root / "misc" / "frame1.fits").write_bytes(b"misc")

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    indexed = _indexed_files(organiser)
    organiser.close()
    assert archived.read_bytes() == b"archive"
    assert rootFrame.read_bytes() == rootBytes
    assert (root / "misc" / "frame1.fits").read_bytes() == b"misc"
    assert indexed == []


def test_every_prepare_warns_again_while_the_skipped_frame_stays_in_the_root(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    firstOrganiser = workspace_organiser(tmp_path, log=log)
    root = Path(firstOrganiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    archived = root / "raw" / NIGHT / "frame1.fits"
    archived.parent.mkdir(parents=True)
    harvestable_raw_fits(archived)
    archivedBytes = archived.read_bytes()
    firstOrganiser.prepare(report=False)
    firstOrganiser.close()

    # ACT
    secondOrganiser = data_organiser(log=log, rootDir=str(root))
    secondOrganiser.prepare(report=False)
    secondOrganiser.close()

    # ASSERT
    assert len(_warnings(log)) == 2
    assert (root / "frame1.fits").is_file()
    assert archived.read_bytes() == archivedBytes


def test_a_destination_that_appears_after_the_filter_is_not_replaced_by_the_move(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE FILTER SEES A FREE DESTINATION, AS IT WOULD IF THE FILE ARRIVED BETWEEN THE FILTER AND THE MOVE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    rootFrame = harvestable_raw_fits(root / "frame1.fits")
    rootBytes = rootFrame.read_bytes()
    archived = root / "raw" / NIGHT / "frame1.fits"
    archived.parent.mkdir(parents=True)
    archived.write_bytes(b"archive")
    monkeypatch.setattr(organiser, "_without_overwriting_frames", lambda frames: frames)

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    indexed = _indexed_files(organiser)
    organiser.close()
    assert archived.read_bytes() == b"archive"
    assert rootFrame.read_bytes() == rootBytes
    assert indexed == []
    [warning] = _warnings(log)
    assert str(root / "frame1.fits") in warning
    assert str(archived.resolve()) in warning
