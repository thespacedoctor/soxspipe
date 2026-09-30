"""Verified, disposable real-data snapshots for opt-in integration tests."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import stat
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath


@dataclass(frozen=True)
class Dataset_manifest:
    """The immutable inventory accepted from a repository-owned manifest."""

    dataset_version: str
    archive_url: str
    files: dict[str, str]


def _safe_path(path: str) -> PurePosixPath:
    candidate = PurePosixPath(path)
    parts = path.rstrip("/").split("/")
    if (
        not path
        or "\\" in path
        or "\x00" in path
        or candidate.is_absolute()
        or any(part in {"", ".", ".."} for part in parts)
        or (candidate.parts and candidate.parts[0].endswith(":"))
    ):
        raise ValueError(f"unsafe archive path: {path!r}")
    return candidate


def _load_manifest(manifest_path: Path) -> Dataset_manifest:
    try:
        contents = json.loads(manifest_path.read_text(encoding="utf-8"))
        dataset_version = contents["dataset_version"]
        archive_url = contents["archive_url"]
        entries = contents["files"]
    except (OSError, KeyError, TypeError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid real-data manifest: {manifest_path}") from error
    if not isinstance(dataset_version, str) or not dataset_version:
        raise ValueError("invalid real-data manifest dataset_version")
    if not isinstance(archive_url, str) or not archive_url.startswith("https://"):
        raise ValueError("real-data manifest archive_url must use HTTPS")
    if not isinstance(entries, list) or not entries:
        raise ValueError("real-data manifest files must be a non-empty list")
    files: dict[str, str] = {}
    for entry in entries:
        if not isinstance(entry, dict) or set(entry) != {"path", "sha256"}:
            raise ValueError("invalid real-data manifest file entry")
        path = entry["path"]
        digest = entry["sha256"]
        if not isinstance(path, str) or not isinstance(digest, str):
            raise ValueError("invalid real-data manifest file entry")
        normalized_path = _safe_path(path).as_posix()
        if normalized_path in files:
            raise ValueError(f"duplicate manifest path: {normalized_path}")
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            raise ValueError(f"invalid SHA-256 for {normalized_path}")
        files[normalized_path] = digest
    if list(files) != sorted(files):
        raise ValueError("real-data manifest files must be sorted")
    return Dataset_manifest(dataset_version, archive_url, files)


def _manifest_hash(manifest_path: Path) -> str:
    return hashlib.sha256(manifest_path.read_bytes()).hexdigest()


def _is_regular_member(member: zipfile.ZipInfo) -> bool:
    mode = member.external_attr >> 16
    return stat.S_IFMT(mode) in {0, stat.S_IFREG}


def _is_root_directory_metadata(member: zipfile.ZipInfo) -> bool:
    if member.filename != "/":
        return False
    mode = member.external_attr >> 16
    if (
        member.file_size != 0
        or not member.is_dir()
        or stat.S_IFMT(mode) == stat.S_IFLNK
    ):
        raise ValueError("invalid root directory metadata entry")
    return True


def _sha256_stream(stream: object) -> str:
    digest = hashlib.sha256()
    while chunk := stream.read(1024 * 1024):
        digest.update(chunk)
    return digest.hexdigest()


def _validate_archive(archive_path: Path, manifest: Dataset_manifest) -> None:
    try:
        with zipfile.ZipFile(archive_path) as archive:
            member_names: set[str] = set()
            file_names: set[str] = set()
            root_directory_seen = False
            for member in archive.infolist():
                if _is_root_directory_metadata(member):
                    if root_directory_seen:
                        raise ValueError("duplicate root directory metadata entry")
                    root_directory_seen = True
                    continue
                normalized_path = _safe_path(member.filename).as_posix()
                if normalized_path in member_names:
                    raise ValueError(f"duplicate archive member: {normalized_path}")
                member_names.add(normalized_path)
                mode = member.external_attr >> 16
                if stat.S_IFMT(mode) == stat.S_IFLNK:
                    raise ValueError(f"archive symlink is not allowed: {normalized_path}")
                if member.is_dir():
                    continue
                if not _is_regular_member(member):
                    raise ValueError(f"archive member is not a regular file: {normalized_path}")
                file_names.add(normalized_path)
                with archive.open(member) as source:
                    digest = _sha256_stream(source)
                if manifest.files.get(normalized_path) != digest:
                    raise ValueError(f"archive checksum mismatch: {normalized_path}")
    except zipfile.BadZipFile as error:
        raise ValueError(f"invalid ZIP archive: {archive_path}") from error
    if file_names != set(manifest.files):
        raise ValueError("archive inventory has missing or extra files")


def _validate_snapshot(snapshot_path: Path, manifest: Dataset_manifest) -> None:
    if snapshot_path.is_symlink() or not snapshot_path.is_dir():
        raise ValueError(f"invalid snapshot directory: {snapshot_path}")
    files: dict[str, str] = {}
    for path in snapshot_path.rglob("*"):
        relative_path = path.relative_to(snapshot_path).as_posix()
        _safe_path(relative_path)
        if path.is_symlink():
            raise ValueError(f"snapshot symlink is not allowed: {relative_path}")
        if path.is_dir():
            continue
        if not path.is_file():
            raise ValueError(f"snapshot member is not a regular file: {relative_path}")
        with path.open("rb") as source:
            files[relative_path] = _sha256_stream(source)
    if files != manifest.files:
        raise ValueError("snapshot inventory or checksum mismatch")


def _make_read_only(path: Path) -> None:
    for member in sorted(path.rglob("*"), reverse=True):
        if not member.is_symlink():
            member.chmod(0o555 if member.is_dir() else 0o444)
    path.chmod(0o555)


def _extract_archive(archive_path: Path, manifest: Dataset_manifest, destination: Path) -> None:
    resolved_destination = destination.resolve()
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            if _is_root_directory_metadata(member):
                continue
            if member.is_dir():
                continue
            relative_path = _safe_path(member.filename)
            target = destination.joinpath(*relative_path.parts)
            if not target.resolve().is_relative_to(resolved_destination):
                raise ValueError(f"unsafe archive path: {member.filename!r}")
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, target.open("xb") as output:
                shutil.copyfileobj(source, output)
    _validate_snapshot(destination, manifest)


def prepare_snapshot(archive_path: Path, manifest_path: Path, cache_root: Path, workspace_path: Path) -> Path:
    """Verify an archive, cache its immutable contents, and create a verified workspace copy."""
    manifest = _load_manifest(manifest_path)
    _validate_archive(archive_path, manifest)
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_path = cache_root / _manifest_hash(manifest_path)
    if cache_path.exists():
        _validate_snapshot(cache_path, manifest)
    else:
        staging_path = Path(tempfile.mkdtemp(prefix="snapshot-", dir=cache_root))
        try:
            _extract_archive(archive_path, manifest, staging_path)
            _make_read_only(staging_path)
            os.replace(staging_path, cache_path)
        except Exception:
            shutil.rmtree(staging_path, ignore_errors=True)
            raise
    if workspace_path.exists():
        _validate_snapshot(workspace_path, manifest)
        raise ValueError(f"workspace already exists: {workspace_path}")
    shutil.copytree(cache_path, workspace_path)
    workspace_path.chmod(0o755)
    for member in workspace_path.rglob("*"):
        if not member.is_symlink():
            member.chmod(0o755 if member.is_dir() else 0o644)
    _validate_snapshot(workspace_path, manifest)
    return workspace_path


def main() -> None:
    """Create a verified, disposable real-data workspace from an archive."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--cache-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--workspace", required=True, type=Path)
    arguments = parser.parse_args()
    print(
        prepare_snapshot(
            arguments.archive,
            arguments.manifest,
            arguments.cache_root,
            arguments.workspace,
        )
    )


if __name__ == "__main__":
    main()
