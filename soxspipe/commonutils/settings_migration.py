#!/usr/bin/env python
"""
*bring the `soxspipe.yaml` of an existing session up to date with the current default settings*

Author
: David Young

Date Created
: November 27, 2026
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, NamedTuple

PACKAGE_DIR = Path(__file__).resolve().parent.parent
REGISTRY_PATH = PACKAGE_DIR / "resources" / "settings_migrations.yaml"
INSTRUMENTS = ("soxs", "xsh")
PATH_SEPARATOR = "/"
SUPERSEDED_KIND = "superseded"
# THE KIND OF A KEY ADDED AT ITS DEFAULT, AND THE `old` VALUE LOGGED FOR IT
MISSING = "missing"
# VALUES THAT BELONG TO ONE WORKSPACE, SO THE DEFAULTS MUST NEVER FILL THEM IN
PER_WORKSPACE_PATHS = frozenset({"logging settings/handlers/file/filename"})


class superseded_default(NamedTuple):
    """One default value that a later release replaced."""

    path: str
    old: Any
    new: Any
    ticket: str


class settings_change(NamedTuple):
    """One edit made to a settings mapping: a superseded value rewritten or a missing key added."""

    kind: str
    path: str
    old: Any
    new: Any
    ticket: str | None


class migration_result(NamedTuple):
    """A migrated settings mapping, the edits that produced it, and the paths left alone for not being a mapping."""

    settings: dict[str, Any]
    changes: tuple[settings_change, ...]
    not_mappings: tuple[str, ...] = ()


def _instrument_key(instrument: str) -> str:
    """Return `soxs` or `xsh` for an instrument name as the data organiser spells it."""
    key = str(instrument).lower()
    if key not in INSTRUMENTS:
        raise ValueError(f"unknown instrument `{instrument}`: expected one of {INSTRUMENTS}")
    return key


def superseded_defaults(instrument: str) -> tuple[superseded_default, ...]:
    """*return the superseded default values of an instrument, oldest first*

    **Key Arguments:**

    - ``instrument`` -- `soxs` or `xsh` (any letter case)

    **Return:**

    - the `superseded_default` entries of the packaged `resources/settings_migrations.yaml`

    **Usage:**

    ```python
    from soxspipe.commonutils.settings_migration import superseded_defaults

    entries = superseded_defaults("soxs")
    ```
    """
    import yaml

    key = _instrument_key(instrument)
    with open(REGISTRY_PATH, encoding="utf-8") as stream:
        registry = yaml.safe_load(stream)
    return tuple(superseded_default(item["path"], item["old"], item["new"], item["ticket"]) for item in registry[key])


def _same_value(left: Any, right: Any) -> bool:
    """Return True when two settings values are equal and of the same YAML types (`True` is not `1`)."""
    import json

    return json.dumps(left, sort_keys=True) == json.dumps(right, sort_keys=True)


def _parent_and_key(settings: Mapping[str, Any], path: str) -> tuple[Any, str] | None:
    """Return the mapping that holds the key at `path` and the key, or None when the path is not in `settings`."""
    *parents, leaf = path.split(PATH_SEPARATOR)
    node: Any = settings
    for part in parents:
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    if not isinstance(node, dict) or leaf not in node:
        return None
    return node, leaf


def _rewrite_superseded(settings: dict[str, Any], superseded: Sequence[superseded_default]) -> list[settings_change]:
    """Replace each value that equals a superseded default with the new default, in registry order."""
    import copy

    changes = []
    for item in superseded:
        found = _parent_and_key(settings, item.path)
        if found is None:
            continue
        node, key = found
        if _same_value(node[key], item.old):
            node[key] = copy.deepcopy(item.new)
            changes.append(settings_change(SUPERSEDED_KIND, item.path, item.old, item.new, item.ticket))
    return changes


def _without_per_workspace(default: Any, path: str) -> Any:
    """Return a copy of a default value without the per-workspace paths beneath `path`."""
    import copy

    if not isinstance(default, dict):
        return copy.deepcopy(default)
    return {
        key: _without_per_workspace(value, f"{path}{PATH_SEPARATOR}{key}")
        for key, value in default.items()
        if f"{path}{PATH_SEPARATOR}{key}" not in PER_WORKSPACE_PATHS
    }


def _add_missing(
    settings: dict[str, Any], defaults: Mapping[str, Any], prefix: str = ""
) -> tuple[list[settings_change], list[str]]:
    """Add every key of `defaults` that `settings` lacks, at any depth, with a copy of its default value.

    Per-workspace paths are never added. Returns the changes and the paths where the settings hold a value that
    is not a mapping while the defaults hold one, which are left as they are.
    """
    changes: list[settings_change] = []
    notMappings: list[str] = []
    for key, default in defaults.items():
        path = f"{prefix}{key}"
        if path in PER_WORKSPACE_PATHS:
            continue
        if key not in settings:
            value = _without_per_workspace(default, path)
            settings[key] = value
            changes.append(settings_change(MISSING, path, MISSING, value, None))
        elif isinstance(default, dict):
            if isinstance(settings[key], dict):
                nested, nestedNotMappings = _add_missing(settings[key], default, f"{path}{PATH_SEPARATOR}")
                changes.extend(nested)
                notMappings.extend(nestedNotMappings)
            else:
                notMappings.append(path)
    return changes, notMappings


def apply_migration(
    settings: Mapping[str, Any], defaults: Mapping[str, Any], superseded: Sequence[superseded_default]
) -> migration_result:
    """*migrate a settings mapping to the current defaults*

    A value equal to a superseded default becomes the new default. A key that the defaults have and the settings
    lack is added at its default. Every other key and value stays as it is, including keys the defaults no longer
    have and values the user chose. Neither input is changed.

    **Key Arguments:**

    - ``settings`` -- the settings of the session
    - ``defaults`` -- the current default settings of the instrument
    - ``superseded`` -- the superseded defaults of the instrument (see `superseded_defaults`)

    **Return:**

    - a `migration_result` with the migrated settings, the `settings_change` of every edit, and the paths whose
      value is not a mapping where the defaults have one

    **Usage:**

    ```python
    from soxspipe.commonutils.settings_migration import apply_migration, superseded_defaults

    result = apply_migration(settings, defaults, superseded_defaults("soxs"))
    ```
    """
    import copy

    migrated = copy.deepcopy(dict(settings))
    changes = _rewrite_superseded(migrated, superseded)
    added, notMappings = _add_missing(migrated, defaults)
    changes.extend(added)
    return migration_result(migrated, tuple(changes), tuple(notMappings))


def _write_yaml_atomically(yamlPath: Path, settings: Mapping[str, Any]) -> None:
    """Write `settings` to `yamlPath` through a temporary file, so a failed write leaves the old file."""
    import shutil
    import tempfile

    import yaml

    handle, tempName = tempfile.mkstemp(dir=yamlPath.parent, prefix=f".{yamlPath.name}.", suffix=".tmp")
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            yaml.safe_dump(settings, stream, sort_keys=False, default_flow_style=None, width=1000, allow_unicode=True)
            stream.flush()
            os.fsync(stream.fileno())
        shutil.copymode(yamlPath, tempName)
        os.replace(tempName, yamlPath)
    except BaseException:
        Path(tempName).unlink(missing_ok=True)
        raise


def _log_result(log: Any, result: migration_result, yamlPath: Path) -> None:
    """Log each change at info level, a warning for each value left alone, and one summary line with the counts."""
    for change in result.changes:
        ticket = f" ({change.ticket})" if change.ticket else ""
        log.info(f"settings migration: `{change.path}` was {change.old}, now {change.new}{ticket}")
    for path in result.not_mappings:
        log.warning(
            f"settings migration: `{path}` in `{yamlPath}` is not a mapping, so its default keys were not added"
        )
    superseded = sum(change.kind == SUPERSEDED_KIND for change in result.changes)
    missing = len(result.changes) - superseded
    summary = f"`{yamlPath}`: {superseded} superseded and {missing} missing settings changed"
    if result.changes:
        log.print(f"Brought {summary}")
    else:
        log.debug(f"Settings migration found nothing to change in {summary}")


def migrate_session_settings(yamlPath: str | os.PathLike[str], instrument: str, log: Any) -> migration_result:
    """*bring a session's `soxspipe.yaml` up to date with the current defaults*

    Rewrites the values that equal a superseded default (see `superseded_defaults`) and adds the keys that the
    instrument's default settings file has and the yaml lacks. Other keys and values are not changed. Comments in the
    yaml do not survive a rewrite, and no backup is written. The file is not touched when nothing changes.

    **Key Arguments:**

    - ``yamlPath`` -- path to the session `soxspipe.yaml`
    - ``instrument`` -- `soxs` or `xsh` (any letter case), selecting the default settings file and the registry
    - ``log`` -- logger with `debug`, `info`, `warning` and `print` methods

    **Return:**

    - a `migration_result` with the migrated settings and the changes made

    **Raises:**

    - `ValueError` when the instrument is unknown or the yaml is not a mapping. The file is left as it is.

    **Usage:**

    ```python
    from soxspipe.commonutils.settings_migration import migrate_session_settings

    migrate_session_settings("sessions/base/soxspipe.yaml", "soxs", log)
    ```
    """
    import yaml

    key = _instrument_key(instrument)
    # RESOLVE THE PATH, SO THE REWRITE REPLACES THE TARGET OF A SYMLINK AND NOT THE LINK
    yamlPath = Path(yamlPath).resolve()
    with open(yamlPath, encoding="utf-8") as stream:
        settings = yaml.safe_load(stream)
    if not isinstance(settings, dict):
        raise ValueError(f"`{yamlPath}` is not a mapping of settings")
    with open(PACKAGE_DIR / f"{key}_default_settings.yaml", encoding="utf-8") as stream:
        defaults = yaml.safe_load(stream)

    result = apply_migration(settings, defaults, superseded_defaults(key))
    if result.changes:
        _write_yaml_atomically(yamlPath, result.settings)
    _log_result(log, result, yamlPath)
    return result
