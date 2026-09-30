"""Security contracts for paths owned by reduction workspaces."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

from soxspipe.commonutils import data_organiser
from soxspipe.commonutils.data_organiser import (
    _UnsafePathError,
    _validate_owned_path,
)
from tests.factories import workspace_organiser

pytestmark = pytest.mark.integration


@pytest.mark.parametrize(
    "sessionId", ["../escape", "nested/session", "x" * 17, "my-session"]
)
def test_session_create_rejects_invalid_identifier_before_writing(
    tmp_path: Path,
    log: object,
    sessionId: str,
) -> None:
    """Reject identifiers outside the documented session-name grammar."""
    organiser = workspace_organiser(tmp_path, log=log)

    with pytest.raises(ValueError, match="Session ID"):
        organiser.session_create(sessionId=sessionId)

    assert list(Path(organiser.sessionsDir).iterdir()) == []


@pytest.mark.parametrize("sessionId", ["../escape", "my-session"])
def test_constructor_rejects_traversal_in_active_session_file(
    tmp_path: Path,
    log: object,
    sessionId: str,
) -> None:
    """Do not construct a session path from an untrusted active-session file."""
    workspacePath = tmp_path / "workspace"
    sessionsPath = workspacePath / "sessions"
    sessionsPath.mkdir(parents=True)
    (sessionsPath / ".sessionid").write_text(sessionId, encoding="utf-8")

    with pytest.raises(ValueError, match="Session ID"):
        data_organiser(log=log, rootDir=str(workspacePath), dbConnect=False)


@pytest.mark.parametrize("isBroken", [False, True])
def test_constructor_rejects_active_session_file_symlink_outside_workspace(
    tmp_path: Path,
    log: object,
    isBroken: bool,
) -> None:
    """Do not read or later overwrite an externally targeted control file."""
    workspacePath = tmp_path / "workspace"
    sessionsPath = workspacePath / "sessions"
    sessionsPath.mkdir(parents=True)
    outsidePath = tmp_path / "outside-session-id"
    if not isBroken:
        outsidePath.write_text("base", encoding="utf-8")
    (sessionsPath / ".sessionid").symlink_to(outsidePath)

    with pytest.raises(ValueError, match="session ID path"):
        data_organiser(log=log, rootDir=str(workspacePath), dbConnect=False)

    if not isBroken:
        assert outsidePath.read_text(encoding="utf-8") == "base"


def test_constructor_rejects_sessions_directory_symlink_outside_workspace(
    tmp_path: Path,
    log: object,
) -> None:
    """Anchor the complete session tree to the configured workspace root."""
    workspacePath = tmp_path / "workspace"
    workspacePath.mkdir()
    outsidePath = tmp_path / "outside-sessions"
    outsidePath.mkdir()
    (workspacePath / "sessions").symlink_to(outsidePath, target_is_directory=True)

    with pytest.raises(ValueError, match="sessions directory"):
        data_organiser(log=log, rootDir=str(workspacePath), dbConnect=False)


@pytest.mark.parametrize(
    ("workspaceName", "label", "targetIsDirectory"),
    [
        ("raw", "raw directory", True),
        ("misc", "misc directory", True),
        ("soxspipe.db", "database path", False),
    ],
)
def test_constructor_rejects_workspace_asset_symlink_outside_workspace(
    tmp_path: Path,
    log: object,
    workspaceName: str,
    label: str,
    targetIsDirectory: bool,
) -> None:
    """Reject symlinked workspace paths before a later operation can use them."""
    workspacePath = tmp_path / "workspace"
    workspacePath.mkdir()
    outsidePath = tmp_path / f"outside-{workspaceName}"
    if targetIsDirectory:
        outsidePath.mkdir()
    else:
        sqlite3.connect(outsidePath).close()
    (workspacePath / workspaceName).symlink_to(
        outsidePath, target_is_directory=targetIsDirectory
    )

    with pytest.raises(_UnsafePathError, match=label):
        data_organiser(log=log, rootDir=str(workspacePath), dbConnect=False)

    assert outsidePath.exists()


def test_session_switch_rejects_directory_symlink_outside_workspace(
    tmp_path: Path,
    log: object,
) -> None:
    """Do not activate a session whose directory resolves outside sessions."""
    organiser = workspace_organiser(tmp_path, log=log)
    sessionsPath = Path(organiser.sessionsDir)
    outsidePath = tmp_path / "outside-session"
    outsidePath.mkdir()
    (sessionsPath / "external").symlink_to(outsidePath, target_is_directory=True)
    (sessionsPath / ".sessionid").write_text("base", encoding="utf-8")

    with pytest.raises(ValueError, match="session path"):
        organiser.session_switch("external")

    assert (sessionsPath / ".sessionid").read_text(encoding="utf-8") == "base"


@pytest.mark.parametrize("sofName", ["../escape.sof", "/tmp/escape.sof"])
def test_sof_writer_rejects_non_leaf_output_names(
    tmp_path: Path,
    log: object,
    sofName: str,
) -> None:
    """Keep generated SOFs directly inside the active session SOF directory."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.sessionPath = str(Path(organiser.sessionsDir) / "base")
    Path(organiser.sessionPath).mkdir()
    organiser.conn = sqlite3.connect(":memory:")
    pd.DataFrame(
        [{"sof": sofName, "filepath": "./raw/a.fits", "tag": "BIAS_VIS", "complete": 1}]
    ).to_sql("sof_map_base", organiser.conn, index=False)

    with pytest.raises(_UnsafePathError, match="SOF filename"):
        organiser._write_sof_files()


def test_owned_path_validation_normalizes_symlink_loop(
    tmp_path: Path,
) -> None:
    """Report a symlink-resolution loop as an unsafe workspace path."""
    ownerPath = tmp_path / "workspace"
    ownerPath.mkdir()
    firstLink = ownerPath / "first"
    secondLink = ownerPath / "second"
    firstLink.symlink_to(secondLink)
    secondLink.symlink_to(firstLink)

    with pytest.raises(_UnsafePathError, match="Unsafe test path"):
        _validate_owned_path(firstLink, ownerPath, "test path")


def test_owned_path_validation_normalizes_owner_symlink_loop(
    tmp_path: Path,
) -> None:
    """Report a symlink-resolution loop in the owning workspace path."""
    firstLink = tmp_path / "owner-one"
    secondLink = tmp_path / "owner-two"
    firstLink.symlink_to(secondLink)
    secondLink.symlink_to(firstLink)

    with pytest.raises(_UnsafePathError, match="Unsafe test owner"):
        _validate_owned_path(tmp_path / "candidate", firstLink, "test owner")


def test_sof_writer_rejects_destination_symlink_outside_session(
    tmp_path: Path,
    log: object,
) -> None:
    """Reject an existing SOF symlink that redirects outside the session."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.sessionPath = str(Path(organiser.sessionsDir) / "base")
    sofDirectory = Path(organiser.sessionPath) / "sof"
    sofDirectory.mkdir(parents=True)
    outsidePath = tmp_path / "outside.sof"
    outsidePath.write_text("outside", encoding="utf-8")
    (sofDirectory / "linked.sof").symlink_to(outsidePath)
    organiser.conn = sqlite3.connect(":memory:")
    pd.DataFrame(
        [
            {
                "sof": "linked.sof",
                "filepath": "./raw/a.fits",
                "tag": "BIAS_VIS",
                "complete": 1,
            }
        ]
    ).to_sql("sof_map_base", organiser.conn, index=False)

    with pytest.raises(ValueError, match="SOF path"):
        organiser._write_sof_files()

    assert outsidePath.read_text(encoding="utf-8") == "outside"


def test_sof_writer_rejects_directory_symlink_outside_session(
    tmp_path: Path,
    log: object,
) -> None:
    """Do not generate SOFs through a directory symlink outside the session."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.sessionPath = str(Path(organiser.sessionsDir) / "base")
    Path(organiser.sessionPath).mkdir()
    outsidePath = tmp_path / "outside-sofs"
    outsidePath.mkdir()
    (Path(organiser.sessionPath) / "sof").symlink_to(
        outsidePath, target_is_directory=True
    )
    organiser.conn = sqlite3.connect(":memory:")
    pd.DataFrame(
        [
            {
                "sof": "science.sof",
                "filepath": "./raw/a.fits",
                "tag": "OBJECT_VIS",
                "complete": 1,
            }
        ]
    ).to_sql("sof_map_base", organiser.conn, index=False)

    with pytest.raises(ValueError, match="SOF directory"):
        organiser._write_sof_files()

    assert list(outsidePath.iterdir()) == []


def test_database_connection_rejects_symlink_outside_workspace(
    tmp_path: Path,
    log: object,
) -> None:
    """Do not open a workspace database symlinked to an external file."""
    organiser = workspace_organiser(tmp_path, log=log)
    outsidePath = tmp_path / "outside.db"
    sqlite3.connect(outsidePath).close()
    Path(organiser.rootDbPath).symlink_to(outsidePath)

    with pytest.raises(ValueError, match="database path"):
        organiser._get_or_create_db_connection()


def test_reducer_selection_rejects_database_symlink_outside_workspace(
    tmp_path: Path,
    log: object,
) -> None:
    """Do not follow a reducer database symlink beyond its owning workspace."""
    from soxspipe.commonutils.reducer import reducer

    workspacePath = tmp_path / "workspace"
    workspacePath.mkdir()
    outsidePath = tmp_path / "outside.db"
    sqlite3.connect(outsidePath).close()
    databasePath = workspacePath / "soxspipe.db"
    databasePath.symlink_to(outsidePath)
    collection = reducer.__new__(reducer)
    collection.log = log
    collection.workspaceDirectory = str(workspacePath)
    collection.sessionDB = str(databasePath)
    collection.pathToSettings = False

    with pytest.raises(ValueError, match="database path"):
        collection.select_sof_files_to_process(reductionTarget="all")
