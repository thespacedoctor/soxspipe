"""Unit contracts for reading and classifying preserved session statuses (DY-218)."""

from __future__ import annotations

import hashlib
import sqlite3
import struct
import tempfile
from pathlib import Path

import pytest

from soxspipe.commonutils.session_status_restore import (
    classify_session_sofs,
    format_session_summary,
    open_backup_read_only,
    read_session_snapshots,
    session_status_snapshot,
)
from soxspipe.commonutils.sql_identifiers import UnsafeSqlIdentifierError

pytestmark = pytest.mark.unit

BACKUP_PATH = "/ws/backups/soxspipe_refresh_20260930T101010000000Z.db"


def test_summary_names_every_count_and_the_backup_when_some_statuses_did_not_survive() -> None:
    # ACT
    summary = format_session_summary("science", 212, 2, 3, BACKUP_PATH)

    # ASSERT
    assert summary == (
        "session 'science': 212 statuses restored, 2 changed (frames differ, requeued), "
        "3 dropped (product no longer exists); "
        f"the previous statuses are kept in {BACKUP_PATH}"
    )


def test_summary_omits_the_backup_path_when_every_status_was_restored() -> None:
    # ACT
    summary = format_session_summary("science", 212, 0, 0, BACKUP_PATH)

    # ASSERT
    assert summary == (
        "session 'science': 212 statuses restored, 0 changed (frames differ, requeued), "
        "0 dropped (product no longer exists)"
    )
    assert BACKUP_PATH not in summary


def _snapshot(statuses, frameSets):
    return session_status_snapshot(sessionId="science", statuses=statuses, frameSets=frameSets)


def test_classify_restores_a_sof_whose_frames_are_unchanged() -> None:
    # ARRANGE
    snapshot = _snapshot({"a.sof": "fail"}, {"a.sof": frozenset({"f1", "f2"})})

    # ACT
    result = classify_session_sofs(snapshot, {"a.sof"}, {"a.sof"}, {"a.sof": frozenset({"f2", "f1"})})

    # ASSERT
    assert result.toRestore == {"a.sof": "fail"}
    assert result.pending == frozenset()
    assert result.changed == frozenset()
    assert result.dropped == frozenset()


def test_classify_leaves_a_sof_pending_when_its_frames_differ() -> None:
    # ARRANGE
    snapshot = _snapshot({"a.sof": "pass"}, {"a.sof": frozenset({"f1", "cal_old"})})

    # ACT
    result = classify_session_sofs(snapshot, {"a.sof"}, {"a.sof"}, {"a.sof": frozenset({"f1", "cal_new"})})

    # ASSERT
    assert result.toRestore == {}
    assert result.pending == frozenset({"a.sof"})


def test_classify_drops_a_sof_the_rebuilt_database_no_longer_has() -> None:
    # ARRANGE
    snapshot = _snapshot({"gone.sof": "pass"}, {"gone.sof": frozenset({"f1"})})

    # ACT
    result = classify_session_sofs(snapshot, {"gone.sof"}, set(), {})

    # ASSERT
    assert result.dropped == frozenset({"gone.sof"})
    assert result.toRestore == {}
    assert result.pending == frozenset()


def test_classify_counts_every_sof_changed_when_the_backup_has_no_frame_sets() -> None:
    # ARRANGE
    snapshot = _snapshot({"a.sof": "pass"}, None)

    # ACT
    result = classify_session_sofs(snapshot, {"a.sof"}, {"a.sof"}, {"a.sof": frozenset({"f1"})})

    # ASSERT
    assert result.changed == frozenset({"a.sof"})
    assert result.toRestore == {}
    assert result.pending == frozenset()


def test_classify_restores_a_sof_with_no_frames_before_or_after() -> None:
    # ARRANGE
    snapshot = _snapshot({"a.sof": "pass"}, {})

    # ACT
    result = classify_session_sofs(snapshot, {"a.sof"}, {"a.sof"}, {})

    # ASSERT
    assert result.toRestore == {"a.sof": "pass"}


def test_classify_considers_only_the_candidate_sofs() -> None:
    # ARRANGE
    snapshot = _snapshot(
        {"a.sof": "pass", "b.sof": "fail"},
        {"a.sof": frozenset({"f1"}), "b.sof": frozenset({"f2"})},
    )
    rebuiltFrameSets = {"a.sof": frozenset({"f1"}), "b.sof": frozenset({"f2"})}

    # ACT
    result = classify_session_sofs(snapshot, {"b.sof"}, {"a.sof", "b.sof"}, rebuiltFrameSets)

    # ASSERT
    assert result.toRestore == {"b.sof": "fail"}


def _write_corrupt_but_openable_database(path):
    """Write a SQLite file that opens but fails `PRAGMA quick_check` (a bad freelist count in the header)."""
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA page_size=4096")
    connection.execute("CREATE TABLE product_frames (file TEXT, sof TEXT)")
    connection.executemany(
        "INSERT INTO product_frames VALUES (?, 'a.sof')", [(f"f{i}" + "x" * 50,) for i in range(500)]
    )
    connection.commit()
    connection.execute("VACUUM")
    connection.close()
    # CORRUPT THE FREELIST PAGE COUNT IN THE FILE HEADER (BYTES 36-39)
    with open(path, "r+b") as databaseFile:
        databaseFile.seek(36)
        databaseFile.write(struct.pack(">I", 999))


def _write_backup(path, statuses, frameSets):
    """Write a small preserved database.

    ``statuses`` maps a session ID to a list of ``(sof, status)`` rows and ``frameSets`` maps a session ID to a
    list of ``(sof, file)`` rows. A session missing from ``frameSets`` gets no ``sof_map_<id>`` table.
    """
    connection = sqlite3.connect(path)
    connection.execute("CREATE TABLE product_frames (file TEXT, sof TEXT)")
    for sessionId, rows in statuses.items():
        connection.execute(f"ALTER TABLE product_frames ADD status_{sessionId} TEXT")
        for index, (sof, status) in enumerate(rows):
            connection.execute(
                f"INSERT INTO product_frames (file, sof, status_{sessionId}) VALUES (?, ?, ?)",  # noqa: S608
                (f"product-{sessionId}-{index}.fits", sof, status),
            )
    for sessionId, rows in frameSets.items():
        connection.execute(f"CREATE TABLE sof_map_{sessionId} (sof TEXT, file TEXT)")
        connection.executemany(f"INSERT INTO sof_map_{sessionId} VALUES (?, ?)", rows)  # noqa: S608
    connection.commit()
    connection.close()


def _debug_messages(log):
    return [message for level, message in log.messages if level == "debug"]


def test_read_snapshots_returns_pass_and_fail_statuses_and_frame_sets_per_session(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(
        backup,
        {
            "science": [("a.sof", "pass"), ("b.sof", "fail"), ("c.sof", None)],
            "archive": [("a.sof", "fail")],
        },
        {
            "science": [("a.sof", "f1"), ("a.sof", "f2"), ("b.sof", "f3")],
            "archive": [("a.sof", "f1")],
        },
    )

    # ACT
    snapshots = read_session_snapshots(backup, ["science", "archive"], log)

    # ASSERT
    assert set(snapshots) == {"science", "archive"}
    assert snapshots["science"].sessionId == "science"
    assert snapshots["science"].statuses == {"a.sof": "pass", "b.sof": "fail"}
    assert snapshots["science"].frameSets == {"a.sof": frozenset({"f1", "f2"}), "b.sof": frozenset({"f3"})}
    assert snapshots["archive"].statuses == {"a.sof": "fail"}


def test_read_snapshots_skips_a_session_the_backup_has_no_status_column_for(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass")]}, {"science": [("a.sof", "f1")]})

    # ACT
    snapshots = read_session_snapshots(backup, ["science", "late"], log)

    # ASSERT
    assert set(snapshots) == {"science"}
    assert any("late" in message for message in _debug_messages(log))


def test_read_snapshots_marks_frame_sets_unknown_when_the_backup_has_no_sof_map_table(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass")]}, {})

    # ACT
    snapshots = read_session_snapshots(backup, ["science"], log)

    # ASSERT
    assert snapshots["science"].statuses == {"a.sof": "pass"}
    assert snapshots["science"].frameSets is None
    assert any("sof_map_science" in message for message in _debug_messages(log))


def test_read_snapshots_lets_fail_win_when_a_sof_has_both_statuses(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass"), ("a.sof", "fail")]}, {"science": []})

    # ACT
    snapshots = read_session_snapshots(backup, ["science"], log)

    # ASSERT
    assert snapshots["science"].statuses == {"a.sof": "fail"}
    assert any("a.sof" in message for message in _debug_messages(log))


def test_read_snapshots_skips_status_values_that_are_neither_pass_nor_fail(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass"), ("b.sof", "junk")]}, {"science": []})

    # ACT
    snapshots = read_session_snapshots(backup, ["science"], log)

    # ASSERT
    assert snapshots["science"].statuses == {"a.sof": "pass"}
    assert any("junk" in message for message in _debug_messages(log))


def test_read_snapshots_raises_when_the_backup_fails_its_quick_check(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_corrupt_but_openable_database(backup)

    # ACT / ASSERT
    with pytest.raises(sqlite3.DatabaseError, match="quick check"):
        read_session_snapshots(backup, ["science"], log)


def test_read_snapshots_raises_when_the_backup_does_not_exist(tmp_path, log) -> None:
    # ACT / ASSERT
    with pytest.raises(sqlite3.OperationalError):
        read_session_snapshots(tmp_path / "missing.db", ["science"], log)
    assert not (tmp_path / "missing.db").exists()


def test_read_snapshots_leaves_the_backup_and_its_directory_unchanged(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass")]}, {"science": [("a.sof", "f1")]})
    digestBefore = hashlib.sha256(backup.read_bytes()).hexdigest()
    listingBefore = sorted(path.name for path in tmp_path.iterdir())

    # ACT
    read_session_snapshots(backup, ["science"], log)

    # ASSERT
    assert hashlib.sha256(backup.read_bytes()).hexdigest() == digestBefore
    assert sorted(path.name for path in tmp_path.iterdir()) == listingBefore


@pytest.fixture
def scratch_root(tmp_path, monkeypatch):
    """Send every temporary directory to a known folder, so a leaked scratch directory is visible."""
    root = tmp_path / "scratch-root"
    root.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(root))
    return root


def test_open_backup_read_only_reads_from_a_copy_outside_the_backups_directory(tmp_path, scratch_root) -> None:
    # ARRANGE
    backup = tmp_path / "backups" / "backup.db"
    backup.parent.mkdir()
    _write_backup(backup, {"science": [("a.sof", "pass")]}, {})

    # ACT
    with open_backup_read_only(backup) as connection:
        databasePath = connection.execute("PRAGMA database_list;").fetchone()[2]
        rows = connection.execute("SELECT sof FROM product_frames;").fetchall()
        scratchWhileOpen = list(scratch_root.iterdir())

    # ASSERT
    assert rows == [("a.sof",)]
    assert len(scratchWhileOpen) == 1
    assert scratchWhileOpen[0] in Path(databasePath).parents


def test_open_backup_read_only_removes_the_scratch_directory_after_a_successful_read(tmp_path, scratch_root) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass")]}, {})

    # ACT
    with open_backup_read_only(backup) as connection:
        connection.execute("SELECT 1;")

    # ASSERT
    assert list(scratch_root.iterdir()) == []


def test_open_backup_read_only_removes_the_scratch_directory_when_the_quick_check_fails(tmp_path, scratch_root) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_corrupt_but_openable_database(backup)

    # ACT / ASSERT
    with pytest.raises(sqlite3.DatabaseError, match="quick check"), open_backup_read_only(backup):
        pytest.fail("the body must not run for a database that fails its quick check")
    assert list(scratch_root.iterdir()) == []


def test_open_backup_read_only_removes_the_scratch_directory_when_the_caller_raises(tmp_path, scratch_root) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass")]}, {})

    # ACT / ASSERT
    with pytest.raises(RuntimeError, match="caller failed"), open_backup_read_only(backup):
        raise RuntimeError("caller failed")
    assert list(scratch_root.iterdir()) == []


def test_open_backup_read_only_raises_an_operational_error_for_a_missing_file(tmp_path, scratch_root) -> None:
    # ACT / ASSERT
    with pytest.raises(sqlite3.OperationalError), open_backup_read_only(tmp_path / "missing.db"):
        pytest.fail("the body must not run for a missing database")
    assert list(scratch_root.iterdir()) == []
    assert not (tmp_path / "missing.db").exists()


# THE BACKUP IS UNTRUSTED INPUT, SO ITS VALUES MUST NOT REACH THE LOG OR AN EXCEPTION AT FULL LENGTH

HOSTILE_LENGTH = 10_000


def test_read_snapshots_logs_only_a_bounded_slice_of_a_hostile_status_value(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    hostileStatus = "S" * HOSTILE_LENGTH
    _write_backup(backup, {"science": [("a.sof", hostileStatus)]}, {"science": []})

    # ACT
    read_session_snapshots(backup, ["science"], log)

    # ASSERT
    messages = _debug_messages(log)
    assert any("a.sof" in message for message in messages)
    assert all(len(message) < 500 for message in messages)
    assert all("S" * 200 not in message for message in messages)


def test_read_snapshots_logs_only_a_bounded_slice_of_a_hostile_sof_name(tmp_path, log) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    hostileSof = "H" * HOSTILE_LENGTH
    _write_backup(
        backup, {"science": [(hostileSof, "pass"), (hostileSof, "fail"), (hostileSof, "junk")]}, {"science": []}
    )

    # ACT
    snapshots = read_session_snapshots(backup, ["science"], log)

    # ASSERT
    assert snapshots["science"].statuses == {hostileSof: "fail"}
    messages = _debug_messages(log)
    assert messages
    assert all(len(message) < 500 for message in messages)
    assert all("H" * 200 not in message for message in messages)


def test_quick_check_failure_message_is_bounded_however_long_the_sqlite_output(tmp_path, log, monkeypatch) -> None:
    # ARRANGE
    hostileOutput = "*** in database main ***\n" + "P" * HOSTILE_LENGTH
    backup = tmp_path / "backup.db"
    backup.write_bytes(b"")

    class _FakeConnection:
        def execute(self, sqlQuery):
            return self

        def fetchall(self):
            return [(hostileOutput,)]

        def close(self):
            pass

    monkeypatch.setattr(sqlite3, "connect", lambda *args, **kwargs: _FakeConnection())

    # ACT
    with pytest.raises(sqlite3.DatabaseError, match="quick check") as caught:
        read_session_snapshots(backup, ["science"], log)

    # ASSERT
    assert len(str(caught.value)) < 500
    assert "P" * 200 not in str(caught.value)


# THE SNAPSHOT AND RESULT ARE FROZEN DATACLASSES, SO THE MAPPINGS THEY HOLD MUST BE READ-ONLY TOO


def test_snapshot_statuses_cannot_be_mutated() -> None:
    # ARRANGE
    source = {"a.sof": "pass"}
    snapshot = _snapshot(source, {"a.sof": frozenset({"f1"})})

    # ACT / ASSERT
    with pytest.raises(TypeError):
        snapshot.statuses["a.sof"] = "fail"
    source["b.sof"] = "fail"
    assert dict(snapshot.statuses) == {"a.sof": "pass"}


def test_snapshot_frame_sets_cannot_be_mutated() -> None:
    # ARRANGE
    snapshot = _snapshot({"a.sof": "pass"}, {"a.sof": frozenset({"f1"})})

    # ACT / ASSERT
    with pytest.raises(TypeError):
        snapshot.frameSets["b.sof"] = frozenset()
    with pytest.raises(AttributeError):
        snapshot.frameSets["a.sof"].add("f2")


def test_snapshot_keeps_unknown_frame_sets_as_none() -> None:
    # ACT
    snapshot = _snapshot({"a.sof": "pass"}, None)

    # ASSERT
    assert snapshot.frameSets is None


def test_snapshot_freezes_inner_frame_sets_given_as_plain_sets() -> None:
    # ARRANGE
    inner = {"f1"}

    # ACT
    snapshot = _snapshot({"a.sof": "pass"}, {"a.sof": inner})
    inner.add("f2")

    # ASSERT
    assert snapshot.frameSets["a.sof"] == frozenset({"f1"})


def test_classify_result_to_restore_cannot_be_mutated() -> None:
    # ARRANGE
    snapshot = _snapshot({"a.sof": "pass"}, {"a.sof": frozenset({"f1"})})
    result = classify_session_sofs(snapshot, {"a.sof"}, {"a.sof"}, {"a.sof": frozenset({"f1"})})

    # ACT / ASSERT
    with pytest.raises(TypeError):
        result.toRestore["b.sof"] = "fail"


@pytest.mark.parametrize("sessionId", ["x; DROP TABLE product_frames", 'a"b', "a-b", "a" * 60])
def test_read_snapshots_rejects_a_session_id_that_is_not_a_safe_identifier(tmp_path, log, sessionId) -> None:
    # ARRANGE
    backup = tmp_path / "backup.db"
    _write_backup(backup, {"science": [("a.sof", "pass")]}, {"science": []})

    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        read_session_snapshots(backup, [sessionId], log)
