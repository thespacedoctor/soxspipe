"""Unit tests for safe and deterministic FITS decompression commands."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from soxspipe.commonutils import uncompress

pytestmark = pytest.mark.unit


class _CompletedProcess:
    def __init__(self, stderr: bytes = b"", returncode: int = 0) -> None:
        self.stderr = stderr
        self.returncode = returncode

    def communicate(self) -> tuple[bytes, bytes]:
        return b"", self.stderr


def test_uncompress_passes_sorted_paths_without_a_shell(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    secondPath = tmp_path / "b frame.fits.Z"
    firstPath = tmp_path / "a;frame.fits.Z"
    secondPath.touch()
    firstPath.touch()
    calls: list[tuple[object, dict[str, object]]] = []

    def record_popen(command: object, **kwargs: object) -> _CompletedProcess:
        calls.append((command, kwargs))
        return _CompletedProcess()

    monkeypatch.setattr(subprocess, "Popen", record_popen)

    result = uncompress(log=log, directory=tmp_path)

    assert result is None
    assert calls == [
        (
            ["uncompress", "-f", str(firstPath), str(secondPath)],
            {"stdout": subprocess.PIPE, "stderr": subprocess.PIPE},
        )
    ]


def test_uncompress_ignores_non_archives(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "raw.fits").touch()
    (tmp_path / "notes.Z").touch()

    def fail_popen(*args: object, **kwargs: object) -> None:
        raise AssertionError("No decompression process should be started")

    monkeypatch.setattr(subprocess, "Popen", fail_popen)

    assert uncompress(log=log, directory=tmp_path) is None


def test_uncompress_reports_a_missing_executable(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (tmp_path / "raw.fits.Z").touch()

    def missing_popen(*args: object, **kwargs: object) -> None:
        raise FileNotFoundError("uncompress")

    monkeypatch.setattr(subprocess, "Popen", missing_popen)

    with pytest.raises(SystemExit) as error:
        uncompress(log=log, directory=tmp_path)

    assert error.value.code == 0
    assert "The uncompress command was not found" in capsys.readouterr().out


def test_uncompress_reports_command_not_found_stderr(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (tmp_path / "raw.fits.Z").touch()
    monkeypatch.setattr(
        subprocess,
        "Popen",
        lambda *args, **kwargs: _CompletedProcess(
            b"uncompress: command not found",
            returncode=127,
        ),
    )

    with pytest.raises(SystemExit) as error:
        uncompress(log=log, directory=tmp_path)

    assert error.value.code == 0
    assert "uncompress: command not found" in capsys.readouterr().out


def test_uncompress_logs_a_nonzero_process_result_without_misclassifying_it(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    archivePath = tmp_path / "uncompress-data.fits.Z"
    archivePath.touch()
    monkeypatch.setattr(
        subprocess,
        "Popen",
        lambda *args, **kwargs: _CompletedProcess(
            f"{archivePath}: corrupt input".encode("ascii"),
            returncode=1,
        ),
    )

    assert uncompress(log=log, directory=tmp_path) is None
    assert "The uncompress command was not found" not in capsys.readouterr().out
    assert (
        "error",
        f"Could not uncompress .Z files (exit code 1): {archivePath}: corrupt input",
    ) in log.messages


def test_uncompress_does_not_report_success_for_a_nonzero_empty_result(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    (tmp_path / "raw.fits.Z").touch()
    monkeypatch.setattr(
        subprocess,
        "Popen",
        lambda *args, **kwargs: _CompletedProcess(returncode=1),
    )

    assert uncompress(log=log, directory=tmp_path) is None
    assert "Decompressed" not in capsys.readouterr().out


def test_uncompress_logs_other_process_errors(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "raw.fits.Z").touch()

    def denied_popen(*args: object, **kwargs: object) -> None:
        raise PermissionError("denied")

    monkeypatch.setattr(subprocess, "Popen", denied_popen)

    assert uncompress(log=log, directory=tmp_path) is None
    assert (
        "error",
        "Could not uncompress .Z files: denied",
    ) in log.messages


CURSOR_UP_AND_CLEAR = "\x1b[1A\x1b[2K"


def _make_archives(directory: Path, count: int) -> list[str]:
    paths = []
    for index in range(count):
        path = directory / f"frame_{index:03d}.fits.Z"
        path.touch()
        paths.append(str(path))
    return paths


def test_uncompress_splits_archives_into_batches_of_25_in_sorted_order(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    expectedPaths = _make_archives(tmp_path, 60)
    commands: list[list[str]] = []

    def record_popen(command: list[str], **kwargs: object) -> _CompletedProcess:
        commands.append(command)
        return _CompletedProcess()

    monkeypatch.setattr(subprocess, "Popen", record_popen)

    uncompress(log=log, directory=tmp_path)

    assert commands == [
        ["uncompress", "-f", *expectedPaths[0:25]],
        ["uncompress", "-f", *expectedPaths[25:50]],
        ["uncompress", "-f", *expectedPaths[50:60]],
    ]


def test_uncompress_does_not_start_an_empty_batch_for_an_exact_multiple_of_25(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _make_archives(tmp_path, 50)
    batchSizes: list[int] = []

    def record_popen(command: list[str], **kwargs: object) -> _CompletedProcess:
        batchSizes.append(len(command) - 2)
        return _CompletedProcess()

    monkeypatch.setattr(subprocess, "Popen", record_popen)

    uncompress(log=log, directory=tmp_path)

    assert batchSizes == [25, 25]


def test_uncompress_overwrites_the_progress_line_after_the_first_batch(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _make_archives(tmp_path, 30)
    monkeypatch.setattr(
        subprocess, "Popen", lambda *args, **kwargs: _CompletedProcess()
    )

    uncompress(log=log, directory=tmp_path)

    assert capsys.readouterr().out == (
        "Decompressed 25/30 fits.Z files (83.3%)\n"
        + CURSOR_UP_AND_CLEAR
        + "Decompressed 30/30 fits.Z files (100.0%)\n"
    )


def test_uncompress_carries_on_with_later_batches_after_an_os_error(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _make_archives(tmp_path, 30)
    callCount = 0

    def flaky_popen(command: list[str], **kwargs: object) -> _CompletedProcess:
        nonlocal callCount
        callCount += 1
        if callCount == 1:
            raise PermissionError("denied")
        return _CompletedProcess()

    monkeypatch.setattr(subprocess, "Popen", flaky_popen)

    uncompress(log=log, directory=tmp_path)

    assert callCount == 2
    assert ("error", "Could not uncompress .Z files: denied") in log.messages
    assert capsys.readouterr().out == (
        CURSOR_UP_AND_CLEAR + "Decompressed 30/30 fits.Z files (100.0%)\n"
    )


def test_uncompress_carries_on_with_later_batches_after_a_nonzero_exit(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    _make_archives(tmp_path, 30)
    results = iter([_CompletedProcess(b"bad", returncode=2), _CompletedProcess()])
    monkeypatch.setattr(subprocess, "Popen", lambda *args, **kwargs: next(results))

    uncompress(log=log, directory=tmp_path)

    assert ("error", "Could not uncompress .Z files (exit code 2): bad") in log.messages
    assert capsys.readouterr().out == (
        CURSOR_UP_AND_CLEAR + "Decompressed 30/30 fits.Z files (100.0%)\n"
    )


def test_uncompress_skips_directories_named_like_archives(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    (tmp_path / "folder.fits.Z").mkdir()
    monkeypatch.setattr(
        subprocess,
        "Popen",
        lambda *args, **kwargs: pytest.fail("no process expected"),
    )

    assert uncompress(log=log, directory=tmp_path) is None
