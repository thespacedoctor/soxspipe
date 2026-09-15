"""Unit tests for the verified real-data ZIP snapshot contract."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
import zipfile
from pathlib import Path

import pytest

from tests.real_data.snapshot import _manifest_hash, prepare_snapshot

pytestmark = pytest.mark.unit


@pytest.fixture
def cache_root(tmp_path: Path) -> Path:
    """Return a snapshot cache root that cannot collide with the autouse XDG cache directory."""
    return tmp_path / "snapshot-cache"


def _cached_snapshot_path(snapshot_cache_root: Path, manifest_path: Path) -> Path:
    # PREPARE_SNAPSHOT NAMES THE CACHED SNAPSHOT AFTER THE MANIFEST DIGEST
    return snapshot_cache_root / _manifest_hash(manifest_path)


def _write_zip(path: Path, members: dict[str, bytes], symlinks: tuple[str, ...] = ()) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        for name, content in members.items():
            info = zipfile.ZipInfo(name)
            if name in symlinks:
                info.external_attr = (stat.S_IFLNK | 0o777) << 16
            archive.writestr(info, content)


def _write_manifest(path: Path, archive_url: str, members: dict[str, bytes]) -> None:
    files = [
        {"path": name, "sha256": hashlib.sha256(content).hexdigest()}
        for name, content in sorted(members.items())
    ]
    path.write_text(json.dumps({"dataset_version": "test-v1", "archive_url": archive_url, "files": files}), encoding="utf-8")


def test_prepare_snapshot_verifies_archive_and_copies_disposable_workspace(cache_root: Path, tmp_path: Path) -> None:
    members = {"raw/frame.fits": b"fits", "sofs/reduce.sof": b"raw/frame.fits RAW"}
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, members)
    _write_manifest(manifest_path, "https://example.test/dataset.zip", members)

    workspace_path = prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")

    assert workspace_path == tmp_path / "workspace"
    assert (workspace_path / "raw/frame.fits").read_bytes() == b"fits"
    assert not cache_root.is_symlink()


def test_prepare_snapshot_cache_root_holds_only_the_cached_snapshot(cache_root: Path, tmp_path: Path) -> None:
    members = {"frame.fits": b"fits"}
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, members)
    _write_manifest(manifest_path, "https://example.test/dataset.zip", members)

    prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")

    # THE AUTOUSE ISOLATED_RUNTIME FIXTURE OWNS XDG_CACHE_HOME, SO THE SNAPSHOT CACHE MUST LIVE ELSEWHERE
    assert cache_root != Path(os.environ["XDG_CACHE_HOME"])
    assert list(cache_root.iterdir()) == [_cached_snapshot_path(cache_root, manifest_path)]


def test_prepare_snapshot_accepts_only_the_root_directory_metadata_entry(cache_root: Path, tmp_path: Path) -> None:
    members = {"raw/frame.fits": b"fits"}
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("/", b"")
        archive.writestr("raw/frame.fits", b"fits")
    _write_manifest(manifest_path, "https://example.test/dataset.zip", members)

    workspace_path = prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")

    assert workspace_path.joinpath("raw/frame.fits").read_bytes() == b"fits"


@pytest.mark.parametrize("root_content", [b"unexpected", b""])
def test_prepare_snapshot_rejects_invalid_root_directory_metadata(cache_root: Path, tmp_path: Path, root_content: bytes) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    root_info = zipfile.ZipInfo("/")
    if not root_content:
        root_info.external_attr = (stat.S_IFLNK | 0o777) << 16
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(root_info, root_content)
        archive.writestr("frame.fits", b"fits")
    _write_manifest(manifest_path, "https://example.test/dataset.zip", {"frame.fits": b"fits"})

    with pytest.raises(ValueError, match="root directory metadata"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


def test_prepare_snapshot_rejects_duplicate_root_directory_metadata(cache_root: Path, tmp_path: Path) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    with pytest.warns(UserWarning, match="Duplicate name"):
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr("/", b"")
            archive.writestr("/", b"")
            archive.writestr("frame.fits", b"fits")
    _write_manifest(manifest_path, "https://example.test/dataset.zip", {"frame.fits": b"fits"})

    with pytest.raises(ValueError, match="duplicate root"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


def test_prepare_snapshot_rejects_non_root_absolute_directory_member(cache_root: Path, tmp_path: Path) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("/escape/", b"")
        archive.writestr("frame.fits", b"fits")
    _write_manifest(manifest_path, "https://example.test/dataset.zip", {"frame.fits": b"fits"})

    with pytest.raises(ValueError, match="unsafe"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


@pytest.mark.parametrize("member_name", ["/escape", "../escape", "raw\\..\\escape", "C:/escape"])
def test_prepare_snapshot_rejects_unsafe_archive_members_before_extraction(cache_root: Path, tmp_path: Path, member_name: str) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, {member_name: b"unsafe"})
    _write_manifest(manifest_path, "https://example.test/dataset.zip", {member_name: b"unsafe"})

    with pytest.raises(ValueError, match="unsafe"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")

    assert not (tmp_path / "workspace").exists()


def test_prepare_snapshot_rejects_symlink_and_extra_members(cache_root: Path, tmp_path: Path) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, {"frame.fits": b"fits", "link": b"frame.fits"}, ("link",))
    _write_manifest(manifest_path, "https://example.test/dataset.zip", {"frame.fits": b"fits"})

    with pytest.raises(ValueError, match="symlink"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


def test_prepare_snapshot_rejects_missing_extra_and_duplicate_archive_members(cache_root: Path, tmp_path: Path) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_manifest(manifest_path, "https://example.test/dataset.zip", {"frame.fits": b"fits"})

    _write_zip(archive_path, {"frame.fits": b"fits", "extra.fits": b"extra"})
    with pytest.raises(ValueError, match="checksum mismatch"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")

    _write_zip(archive_path, {"different.fits": b"fits"})
    with pytest.raises(ValueError, match="checksum mismatch"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")

    with pytest.warns(UserWarning, match="Duplicate name"):
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr("frame.fits", b"fits")
            archive.writestr("frame.fits", b"fits")
    with pytest.raises(ValueError, match="duplicate"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


def test_prepare_snapshot_rejects_digest_mismatch_and_tampered_cache(cache_root: Path, tmp_path: Path) -> None:
    members = {"frame.fits": b"fits"}
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, members)
    _write_manifest(manifest_path, "https://example.test/dataset.zip", members)
    workspace_path = prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")
    cache_path = _cached_snapshot_path(cache_root, manifest_path)
    cached_file = cache_path / "frame.fits"
    cached_file.chmod(0o644)
    cached_file.write_bytes(b"changed")
    shutil.rmtree(workspace_path)

    with pytest.raises(ValueError, match="inventory or checksum"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


def test_prepare_snapshot_rejects_a_directory_symlink_in_the_cached_snapshot(cache_root: Path, tmp_path: Path) -> None:
    members = {"frame.fits": b"fits"}
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, members)
    _write_manifest(manifest_path, "https://example.test/dataset.zip", members)
    workspace_path = prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")
    cache_path = _cached_snapshot_path(cache_root, manifest_path)
    cache_path.chmod(0o755)
    (cache_path / "linked-directory").symlink_to(tmp_path, target_is_directory=True)
    shutil.rmtree(workspace_path)

    with pytest.raises(ValueError, match="symlink"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


def test_prepare_snapshot_rejects_manifest_with_duplicate_paths(cache_root: Path, tmp_path: Path) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, {"frame.fits": b"fits"})
    manifest_path.write_text(
        json.dumps({"dataset_version": "test-v1", "archive_url": "https://example.test/dataset.zip", "files": [{"path": "frame.fits", "sha256": "0" * 64}, {"path": "frame.fits", "sha256": "0" * 64}]}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")


def test_prepare_snapshot_rejects_non_https_and_unsorted_manifest_paths(cache_root: Path, tmp_path: Path) -> None:
    archive_path = tmp_path / "dataset.zip"
    manifest_path = tmp_path / "manifest.json"
    _write_zip(archive_path, {"a.fits": b"a", "b.fits": b"b"})
    manifest_path.write_text(
        json.dumps({"dataset_version": "test-v1", "archive_url": "http://example.test/dataset.zip", "files": []}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="HTTPS"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")

    manifest_path.write_text(
        json.dumps({"dataset_version": "test-v1", "archive_url": "https://example.test/dataset.zip", "files": [{"path": "b.fits", "sha256": hashlib.sha256(b"b").hexdigest()}, {"path": "a.fits", "sha256": hashlib.sha256(b"a").hexdigest()}]}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="sorted"):
        prepare_snapshot(archive_path, manifest_path, cache_root, tmp_path / "workspace")
