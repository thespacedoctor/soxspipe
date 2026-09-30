"""Compatibility contracts for the retained legacy test helper."""

from __future__ import annotations

from pathlib import Path

import pytest

from soxspipe.utKit import utKit

pytestmark = pytest.mark.unit


def test_utkit_builds_legacy_paths_and_optional_database_configuration(
    tmp_path: Path,
) -> None:
    withoutDatabase = utKit(str(tmp_path))
    withDatabase = utKit(str(tmp_path), dbConn=True)

    assert withoutDatabase.pathToInputDir == f"{tmp_path}/input/"
    assert withoutDatabase.pathToOutputDir == f"{tmp_path}/output/"
    assert withoutDatabase.dbConfig is False
    assert "dryx_unit_testing" in withDatabase.dbConfig
    assert Path(withDatabase.get_project_root()).name == "soxspipe"


def test_utkit_no_longer_defines_dead_refresh_database_method() -> None:
    # DY-32: REFRESH_DATABASE() POINTED AT A SETTINGS FILE THAT NEVER SHIPPED
    # (TEST_SETTINGS.YAML, SINGULAR) AND RAN A MYSQL SCRIPT RUNNER AGAINST A
    # PIPELINE WHOSE PERSISTED STATE IS SQLITE. IT HAD NO REACHABLE CALLER.
    #
    # CHECK THE SOXSPIPE OVERRIDE'S OWN __DICT__, NOT HASATTR(): THE PARENT
    # FUNDAMENTALS.UTKIT CLASS STILL DEFINES REFRESH_DATABASE, SO HASATTR()
    # WOULD REPORT TRUE VIA INHERITANCE REGARDLESS OF THIS DELETION.
    assert "refresh_database" not in utKit.__dict__
