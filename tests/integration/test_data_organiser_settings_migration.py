"""`soxspipe prep` brings the active session's `soxspipe.yaml` up to date before it loads the settings (DY-2475)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.data_organiser import data_organiser
from tests.factories import harvestable_raw_fits, workspace_organiser

pytestmark = pytest.mark.integration

SPECTRAL_DEG = ("soxs-spat-solution", "nir", "order-deg")


def _session_yaml(root: Path) -> Path:
    return root / "sessions" / "base" / "soxspipe.yaml"


def _prepared_workspace(tmp_path: Path, log: Any) -> Path:
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    organiser.prepare(report=False)
    organiser.close()
    return root


def _rewrite_yaml(path: Path, edit: Any) -> None:
    settings = yaml.safe_load(path.read_text())
    edit(settings)
    path.write_text(yaml.safe_dump(settings, sort_keys=False, default_flow_style=None))


def _migration_lines(log: Any) -> list[str]:
    return [
        message
        for level, message in log.messages
        if level != "debug" and ("settings migration" in message.lower() or "superseded and" in message.lower())
    ]


def test_prepare_rewrites_a_superseded_default_and_the_loaded_settings_are_the_new_value(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE
    root = _prepared_workspace(tmp_path, log)
    _rewrite_yaml(_session_yaml(root), lambda s: s[SPECTRAL_DEG[0]][SPECTRAL_DEG[1]].update({"order-deg": [4, 4]}))
    log.messages.clear()
    organiser = data_organiser(log=log, rootDir=str(root))

    # ACT
    organiser.prepare(report=False)
    organiser.close()

    # ASSERT
    assert yaml.safe_load(_session_yaml(root).read_text())["soxs-spat-solution"]["nir"]["order-deg"] == [5, 4]
    assert organiser.settings["soxs-spat-solution"]["nir"]["order-deg"] == [5, 4]
    assert any("soxs-spat-solution/nir/order-deg" in line and "DY-1989" in line for line in _migration_lines(log))


def test_prepare_adds_a_missing_key_before_the_settings_are_loaded(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    root = _prepared_workspace(tmp_path, log)
    _rewrite_yaml(_session_yaml(root), lambda s: s.pop("nir-max-detector-temp-k"))
    organiser = data_organiser(log=log, rootDir=str(root))

    # ACT
    organiser.prepare(report=False)
    organiser.close()

    # ASSERT
    assert organiser.settings["nir-max-detector-temp-k"] == 80


def test_prepare_with_refresh_migrates_the_session_yaml_too(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    root = _prepared_workspace(tmp_path, log)
    _rewrite_yaml(_session_yaml(root), lambda s: s["soxs-spat-solution"]["nir"].update({"order-deg": [4, 4]}))
    organiser = data_organiser(log=log, rootDir=str(root))

    # ACT
    organiser.prepare(refresh=True, report=False)
    organiser.close()

    # ASSERT
    assert yaml.safe_load(_session_yaml(root).read_text())["soxs-spat-solution"]["nir"]["order-deg"] == [5, 4]


def test_prepare_leaves_an_up_to_date_session_yaml_alone_and_logs_no_migration(log: Any, tmp_path: Path) -> None:
    # ARRANGE: THE FIRST PREPARE CREATES THE SESSION, WHICH IS ALREADY CURRENT
    root = _prepared_workspace(tmp_path, log)
    assert _migration_lines(log) == []
    before = _session_yaml(root).read_bytes()
    modifiedBefore = _session_yaml(root).stat().st_mtime_ns
    organiser = data_organiser(log=log, rootDir=str(root))

    # ACT
    organiser.prepare(report=False)
    organiser.close()

    # ASSERT
    assert _session_yaml(root).read_bytes() == before
    assert _session_yaml(root).stat().st_mtime_ns == modifiedBefore
    assert _migration_lines(log) == []


def _organiser_with_session_yaml(tmp_path: Path, log: Any, settings: dict[str, Any]) -> tuple[data_organiser, Path]:
    organiser = workspace_organiser(tmp_path, log=log)
    sessionPath = Path(organiser.rootDir) / "sessions" / "base"
    sessionPath.mkdir(parents=True)
    organiser.sessionPath = str(sessionPath)
    yamlPath = sessionPath / "soxspipe.yaml"
    yamlPath.write_text(yaml.safe_dump(settings, sort_keys=False))
    return organiser, yamlPath


def test_a_session_without_a_known_instrument_is_not_migrated_and_a_warning_says_so(log: Any, tmp_path: Path) -> None:
    # ARRANGE: NO RAW FRAME HAS SET `instrument`, SO THE DEFAULTS TO MIGRATE TO ARE UNKNOWN
    organiser, yamlPath = _organiser_with_session_yaml(tmp_path, log, {"version": 1})
    before = yamlPath.read_bytes()

    # ACT
    organiser._migrate_session_settings()

    # ASSERT
    assert yamlPath.read_bytes() == before
    assert any(level == "warning" and "instrument" in message for level, message in log.messages)


def test_an_xshooter_session_is_migrated_with_the_xshooter_defaults(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    organiser, yamlPath = _organiser_with_session_yaml(tmp_path, log, {"version": 1})
    organiser.instrument = "XSH"
    xshDefaults = yaml.safe_load((Path(soxspipe.__file__).parent / "xsh_default_settings.yaml").read_text())

    # ACT
    organiser._migrate_session_settings()

    # ASSERT
    del xshDefaults["logging settings"]["handlers"]["file"]["filename"]
    assert yaml.safe_load(yamlPath.read_text()) == xshDefaults


@pytest.mark.parametrize("content", ["", "- a\n- list\n", "key: [unclosed\n"])
def test_prepare_carries_on_with_a_warning_when_the_session_yaml_cannot_be_migrated(
    log: Any, tmp_path: Path, content: str
) -> None:
    # ARRANGE
    organiser, yamlPath = _organiser_with_session_yaml(tmp_path, log, {"version": 1})
    organiser.instrument = "SOXS"
    yamlPath.write_text(content)

    # ACT
    organiser._migrate_session_settings()

    # ASSERT
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert str(yamlPath) in warnings[0]
    assert yamlPath.read_text() == content
