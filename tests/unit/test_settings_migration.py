"""Contracts for migrating the `soxspipe.yaml` of an existing session to the current defaults (DY-2475)."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.settings_migration import (
    apply_migration,
    migrate_session_settings,
    superseded_defaults,
)
from soxspipe.commonutils.settings_migration import superseded_default as entry

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
SPECTRAL_DEG_PATH = "soxs-spat-solution/nir/order-deg"


def _shipped_defaults(instrument: str) -> dict[str, Any]:
    return yaml.safe_load((PACKAGE_DIR / f"{instrument}_default_settings.yaml").read_text())


def _value_at(settings: dict[str, Any], path: str) -> Any:
    value: Any = settings
    for part in path.split("/"):
        value = value[part]
    return value


def _set_at(settings: dict[str, Any], path: str, value: Any) -> None:
    *parents, leaf = path.split("/")
    node = settings
    for part in parents:
        node = node[part]
    node[leaf] = value


def _delete_at(settings: dict[str, Any], path: str) -> None:
    *parents, leaf = path.split("/")
    node = settings
    for part in parents:
        node = node[part]
    del node[leaf]


def _write_yaml(path: Path, settings: dict[str, Any]) -> Path:
    path.write_text(yaml.safe_dump(settings, sort_keys=False, default_flow_style=None))
    return path


def _info_lines(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "info"]


def _print_lines(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "print"]


# THE PURE MIGRATION


def test_value_equal_to_an_old_default_becomes_the_new_default() -> None:
    # ARRANGE
    defaults = {"recipe": {"nir": {"order-deg": [5, 4]}}}
    settings = {"recipe": {"nir": {"order-deg": [4, 4]}}}
    registry = [entry("recipe/nir/order-deg", [4, 4], [5, 4], "DY-1989")]

    # ACT
    result = apply_migration(settings, defaults, registry)

    # ASSERT
    assert result.settings == {"recipe": {"nir": {"order-deg": [5, 4]}}}
    assert [(c.kind, c.path, c.old, c.new, c.ticket) for c in result.changes] == [
        ("superseded", "recipe/nir/order-deg", [4, 4], [5, 4], "DY-1989")
    ]


def test_value_that_differs_from_the_old_and_new_defaults_is_not_changed() -> None:
    # ARRANGE
    defaults = {"recipe": {"order-deg": [5, 4]}}
    settings = {"recipe": {"order-deg": [3, 3]}}
    registry = [entry("recipe/order-deg", [4, 4], [5, 4], "DY-1989")]

    # ACT
    result = apply_migration(settings, defaults, registry)

    # ASSERT
    assert result.settings == settings
    assert result.changes == ()


def test_up_to_date_settings_give_no_changes() -> None:
    # ARRANGE
    defaults = {"recipe": {"order-deg": [5, 4]}, "flag": True}
    registry = [entry("recipe/order-deg", [4, 4], [5, 4], "DY-1989")]

    # ACT
    result = apply_migration(copy.deepcopy(defaults), defaults, registry)

    # ASSERT
    assert result.settings == defaults
    assert result.changes == ()


def test_registry_path_missing_from_the_settings_is_not_added_as_a_superseded_value() -> None:
    # ARRANGE: THE MISSING-KEY STEP ADDS IT AT THE CURRENT DEFAULT, NOT THE SUPERSEDED STEP
    defaults = {"recipe": {"order-deg": [5, 4]}}
    registry = [entry("recipe/order-deg", [4, 4], [5, 4], "DY-1989")]

    # ACT
    result = apply_migration({"recipe": {}}, defaults, registry)

    # ASSERT
    assert [c.kind for c in result.changes] == ["missing"]


def test_a_boolean_is_not_taken_for_the_integer_old_default() -> None:
    # ARRANGE
    defaults = {"limit": 5}
    registry = [entry("limit", 1, 5, "DY-1")]

    # ACT
    result = apply_migration({"limit": True}, defaults, registry)

    # ASSERT
    assert result.settings == {"limit": True}
    assert result.changes == ()


def test_a_chain_of_superseded_defaults_ends_at_the_latest_default() -> None:
    # ARRANGE
    defaults = {"limit": 9}
    registry = [entry("limit", 1, 5, "DY-1"), entry("limit", 5, 9, "DY-2")]

    # ACT
    result = apply_migration({"limit": 1}, defaults, registry)

    # ASSERT
    assert result.settings == {"limit": 9}
    assert [(c.old, c.new, c.ticket) for c in result.changes] == [(1, 5, "DY-1"), (5, 9, "DY-2")]


def test_missing_top_level_key_is_added_at_its_default() -> None:
    # ARRANGE
    defaults = {"keep": 1, "nir-max-detector-temp-k": 80}

    # ACT
    result = apply_migration({"keep": 1}, defaults, [])

    # ASSERT
    assert result.settings == defaults
    assert [(c.kind, c.path, c.old, c.new, c.ticket) for c in result.changes] == [
        ("missing", "nir-max-detector-temp-k", "missing", 80, None)
    ]


def test_missing_nested_key_is_added_beside_the_keys_that_are_present() -> None:
    # ARRANGE
    defaults = {"stare": {"nir": {"sky-tuning": {"clip-sigma": 3.0, "global-fit": True}}}}
    settings = {"stare": {"nir": {"sky-tuning": {"clip-sigma": 2.0}}}}

    # ACT
    result = apply_migration(settings, defaults, [])

    # ASSERT
    assert result.settings == {"stare": {"nir": {"sky-tuning": {"clip-sigma": 2.0, "global-fit": True}}}}
    assert [c.path for c in result.changes] == ["stare/nir/sky-tuning/global-fit"]


def test_missing_subtree_is_added_as_one_change() -> None:
    # ARRANGE
    defaults = {"stare": {"nir": {"sky-tuning": {"clip-sigma": 3.0, "global-fit": True}}}}

    # ACT
    result = apply_migration({"stare": {}}, defaults, [])

    # ASSERT
    assert result.settings == defaults
    assert [c.path for c in result.changes] == ["stare/nir"]


def test_key_the_defaults_no_longer_have_is_left_in_place() -> None:
    # ARRANGE
    defaults = {"stare": {"clip-sigma": 3.0}}
    settings = {"stare": {"clip-sigma": 3.0, "object_wing_gate_fraction": 0.1}}

    # ACT
    result = apply_migration(settings, defaults, [])

    # ASSERT
    assert result.settings == settings
    assert result.changes == ()


def test_migration_does_not_change_its_inputs_or_share_default_values() -> None:
    # ARRANGE
    defaults = {"stare": {"nir": {"orders": [1, 2]}}}
    settings: dict[str, Any] = {"stare": {}}
    registry = [entry("x", 1, 2, "DY-1")]

    # ACT
    result = apply_migration(settings, defaults, registry)
    result.settings["stare"]["nir"]["orders"].append(3)

    # ASSERT
    assert settings == {"stare": {}}
    assert defaults == {"stare": {"nir": {"orders": [1, 2]}}}


def test_a_value_that_is_not_a_mapping_where_the_defaults_have_one_is_left_alone() -> None:
    # ARRANGE
    defaults = {"stare": {"clip-sigma": 3.0}}

    # ACT
    result = apply_migration({"stare": None}, defaults, [])

    # ASSERT
    assert result.settings == {"stare": None}
    assert result.changes == ()


# THE REGISTRY


@pytest.mark.parametrize("instrument", ["soxs", "xsh"])
def test_every_registry_path_exists_in_the_current_defaults(instrument: str) -> None:
    # ARRANGE
    defaults = _shipped_defaults(instrument)

    # ACT
    registry = superseded_defaults(instrument)

    # ASSERT
    for item in registry:
        assert _value_at(defaults, item.path) is not None, item.path


@pytest.mark.parametrize("instrument", ["soxs", "xsh"])
def test_the_last_registry_value_of_each_path_is_the_current_default(instrument: str) -> None:
    # ARRANGE
    defaults = _shipped_defaults(instrument)
    registry = superseded_defaults(instrument)

    # ACT
    latest = {item.path: item.new for item in registry}

    # ASSERT
    for path, new in latest.items():
        assert _value_at(defaults, path) == new, path


def test_every_registry_entry_names_a_ticket_and_changes_the_value() -> None:
    # ACT
    registry = superseded_defaults("soxs")

    # ASSERT
    assert registry
    for item in registry:
        assert item.ticket.startswith("DY-"), item.path
        assert item.old != item.new, item.path


def test_the_soxs_registry_holds_the_p3215_values() -> None:
    # ACT
    registry = {item.path: (item.old, item.new) for item in superseded_defaults("soxs")}

    # ASSERT
    assert registry[SPECTRAL_DEG_PATH] == ([4, 4], [5, 4])
    assert registry["soxs-spat-solution/nir/wavelength-deg"] == ([5, 4], [6, 4])
    assert registry["soxs-disp-solution/nir/qc-acceptable-ranges/detlines-frac"] == ([0.9, 1.01], [0.85, 1.01])
    assert registry["soxs-mflat/nir/qc-acceptable-ranges/ordexp90"] == ([30000, 60000], [30000, 70000])
    assert registry["soxs-stare-std/vis/sky-subtraction/min_points_per_knot"] == (155, 77)


def test_the_xsh_registry_may_be_empty() -> None:
    # ACT
    registry = superseded_defaults("xsh")

    # ASSERT
    assert isinstance(registry, tuple)


# THE FILE MIGRATION


def test_every_soxs_registry_entry_is_rewritten_in_a_yaml_file(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    registry = superseded_defaults("soxs")
    for item in registry:
        _set_at(settings, item.path, copy.deepcopy(item.old))
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    migrated = yaml.safe_load(yamlPath.read_text())
    assert migrated == _shipped_defaults("soxs")


def test_a_rewritten_value_is_logged_with_its_path_old_value_new_value_and_ticket(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    _set_at(settings, SPECTRAL_DEG_PATH, [4, 4])
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    logged = [line for line in _info_lines(log) if SPECTRAL_DEG_PATH in line]
    assert len(logged) == 1
    assert "[4, 4]" in logged[0]
    assert "[5, 4]" in logged[0]
    assert "DY-1989" in logged[0]


def test_a_missing_key_is_logged_as_missing_with_its_default(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    _delete_at(settings, "nir-max-detector-temp-k")
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    assert yaml.safe_load(yamlPath.read_text())["nir-max-detector-temp-k"] == 80
    logged = [line for line in _info_lines(log) if "nir-max-detector-temp-k" in line]
    assert len(logged) == 1
    assert "missing" in logged[0]
    assert "80" in logged[0]


def test_one_summary_line_gives_the_counts(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    _set_at(settings, SPECTRAL_DEG_PATH, [4, 4])
    _delete_at(settings, "nir-max-detector-temp-k")
    _delete_at(settings, "soxs-stare/measure-telluric-shift")
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    summary = _print_lines(log)
    assert len(summary) == 1
    assert "1 superseded" in summary[0]
    assert "2 missing" in summary[0]


def test_a_user_value_and_a_per_workspace_value_are_not_changed(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    _set_at(settings, SPECTRAL_DEG_PATH, [7, 7])
    _set_at(settings, "logging settings/handlers/file/filename", "/work/space/soxspipe.log")
    _delete_at(settings, "nir-max-detector-temp-k")
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    migrated = yaml.safe_load(yamlPath.read_text())
    assert _value_at(migrated, SPECTRAL_DEG_PATH) == [7, 7]
    assert _value_at(migrated, "logging settings/handlers/file/filename") == "/work/space/soxspipe.log"


def test_a_key_the_defaults_no_longer_have_stays_in_the_file(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    settings["soxs-nod"]["horne-extraction-profile-global-clipping-sigma"] = 25
    _delete_at(settings, "nir-max-detector-temp-k")
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    migrated = yaml.safe_load(yamlPath.read_text())
    assert migrated["soxs-nod"]["horne-extraction-profile-global-clipping-sigma"] == 25


def test_an_up_to_date_yaml_is_not_rewritten_and_logs_no_changes(log: Any, tmp_path: Path) -> None:
    # ARRANGE: A HAND-WRITTEN FILE WITH A COMMENT, SO A REWRITE WOULD SHOW AS A CHANGED FILE
    settings = _shipped_defaults("soxs")
    yamlPath = tmp_path / "soxspipe.yaml"
    yamlPath.write_text("# KEPT WHEN NOTHING CHANGES\n" + yaml.safe_dump(settings, sort_keys=False))
    before = yamlPath.read_bytes()
    modifiedBefore = yamlPath.stat().st_mtime_ns

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    assert yamlPath.read_bytes() == before
    assert yamlPath.stat().st_mtime_ns == modifiedBefore
    assert _info_lines(log) == []
    assert _print_lines(log) == []


def test_migrating_twice_changes_nothing_the_second_time(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    _set_at(settings, SPECTRAL_DEG_PATH, [4, 4])
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)
    migrate_session_settings(yamlPath, "soxs", log)
    afterFirst = yamlPath.read_bytes()
    log.messages.clear()

    # ACT
    migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    assert yamlPath.read_bytes() == afterFirst
    assert _info_lines(log) == []


def test_xshooter_missing_keys_are_added_from_the_xsh_defaults(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    defaults = _shipped_defaults("xsh")
    settings = copy.deepcopy(defaults)
    topKey = next(key for key, value in defaults.items() if not isinstance(value, dict))
    del settings[topKey]
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", settings)

    # ACT
    migrate_session_settings(yamlPath, "xsh", log)

    # ASSERT
    assert yaml.safe_load(yamlPath.read_text()) == defaults
    assert any(topKey in line for line in _info_lines(log))


def test_the_instrument_name_of_the_data_organiser_selects_the_defaults(log: Any, tmp_path: Path) -> None:
    # ARRANGE: THE DATA ORGANISER NAMES THE INSTRUMENTS `SOXS` AND `XSH`
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", {"version": 1})

    # ACT
    migrate_session_settings(yamlPath, "XSH", log)

    # ASSERT
    assert yaml.safe_load(yamlPath.read_text()) == _shipped_defaults("xsh")


def test_a_yaml_that_is_not_a_mapping_is_refused_and_left_alone(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    yamlPath = tmp_path / "soxspipe.yaml"
    yamlPath.write_text("- not\n- a mapping\n")

    # ACT
    with pytest.raises(ValueError, match="not a mapping"):
        migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    assert yamlPath.read_text() == "- not\n- a mapping\n"


def test_an_unknown_instrument_is_refused(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    yamlPath = _write_yaml(tmp_path / "soxspipe.yaml", {"version": 1})

    # ACT
    with pytest.raises(ValueError, match="instrument"):
        migrate_session_settings(yamlPath, "hubble", log)


def test_a_registry_path_whose_parent_is_missing_is_skipped() -> None:
    # ARRANGE
    registry = [entry("recipe/nir/order-deg", [4, 4], [5, 4], "DY-1989")]

    # ACT
    result = apply_migration({"recipe": {}}, {"recipe": {}}, registry)

    # ASSERT
    assert result.settings == {"recipe": {}}
    assert result.changes == ()


def test_a_failed_write_leaves_the_original_file_and_no_temporary_file(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    settings = _shipped_defaults("soxs")
    _delete_at(settings, "nir-max-detector-temp-k")
    sessionDir = tmp_path / "session"
    sessionDir.mkdir()
    yamlPath = _write_yaml(sessionDir / "soxspipe.yaml", settings)
    before = yamlPath.read_bytes()

    def failing_dump(*args: Any, **kwargs: Any) -> None:
        raise OSError("disk full")

    monkeypatch.setattr(yaml, "safe_dump", failing_dump)

    # ACT
    with pytest.raises(OSError, match="disk full"):
        migrate_session_settings(yamlPath, "soxs", log)

    # ASSERT
    assert yamlPath.read_bytes() == before
    assert [path.name for path in sessionDir.iterdir()] == ["soxspipe.yaml"]
