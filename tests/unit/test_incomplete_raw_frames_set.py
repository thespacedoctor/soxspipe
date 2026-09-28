"""Unit tests for the `missing calibrations` column of `get_incomplete_raw_frames_set`.

DY-265: the table of science sets that cannot be reduced named the sets but
not the calibrations they lack. These tests run the real query against a real
in-memory SQLite database. Plain tables stand in for the `cal_<type>` views,
because the completeness test only reads their `sof` and `upstream_status`
columns.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pytest
import yaml

import soxspipe
from soxspipe.commonutils import data_organiser as data_organiser_module
from soxspipe.commonutils.missing_calibrations import CALIBRATION_DESCRIPTIONS

pytestmark = pytest.mark.unit

SOF_MAP_PATH = Path(soxspipe.__file__).parent / "resources" / "soxs_sof_map.yaml"
CALIBRATION_TYPES = ("mbias", "mdark", "disp_solution", "order_centres", "mflat", "spat_solution", "std_flux")
DISPLAY_COLUMNS = [
    "eso seq arm",
    "mjd-obs",
    "eso dpr tech",
    "eso dpr type",
    "slit",
    "eso obs name",
    "eso obs id",
]


def _connection(calibrationTypes: tuple[str, ...] = CALIBRATION_TYPES) -> sqlite3.Connection:
    connection = sqlite3.connect(":memory:")
    connection.execute(
        'CREATE TABLE raw_frame_sets ("eso seq arm" TEXT, "mjd-obs" REAL, "eso dpr tech" TEXT, '
        '"eso dpr type" TEXT, "slit" TEXT, "eso obs name" TEXT, "eso obs id" INTEGER, '
        "sof TEXT, recipe TEXT, complete INTEGER)"
    )
    for calibrationType in calibrationTypes:
        connection.execute(f"CREATE TABLE cal_{calibrationType} (sof TEXT, upstream_status TEXT)")  # noqa: S608
    return connection


def _add_set(
    connection: sqlite3.Connection,
    sof: str,
    arm: str = "VIS",
    recipe: str = "nod_obj",
    complete: int = 0,
    obsId: int = 1,
    mjd: float = 60000.04,
    slit: str | None = "SLIT1.0",
) -> None:
    connection.execute(
        "INSERT INTO raw_frame_sets VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (arm, mjd, "ECHELLE,SLIT,NODDING", "OBJECT", slit, "target", obsId, sof, recipe, complete),
    )


def _add_calibration(connection: sqlite3.Connection, calibrationType: str, sof: str, status: str | None) -> None:
    connection.execute(f"INSERT INTO cal_{calibrationType} VALUES (?, ?)", (sof, status))  # noqa: S608


def _organiser(connection: sqlite3.Connection, sofMapLookup: dict[str, Any] | None) -> Any:
    organiser = data_organiser_module.__new__(data_organiser_module)
    organiser.conn = connection
    if sofMapLookup is not None:
        organiser.sofMapLookup = sofMapLookup
    return organiser


def _real_sof_map() -> dict[str, Any]:
    with open(SOF_MAP_PATH) as stream:
        return yaml.safe_load(stream)


def test_names_only_the_calibrations_absent_from_the_cal_tables() -> None:
    """A VIS nod set with a master bias but no flats is told what it lacks, in yaml order."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof")
    _add_calibration(connection, "mbias", "NOD_A.sof", "pass")
    _add_calibration(connection, "std_flux", "NOD_A.sof", "pass")

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert result["missing calibrations"].tolist() == ["master flat, spatial solution"]
    assert list(result.columns) == [*DISPLAY_COLUMNS, "missing calibrations"]


def test_uses_the_arm_specific_calibration_list() -> None:
    """A NIR nod set never reports a master bias, because only VIS requires one."""
    connection = _connection()
    _add_set(connection, "NOD_NIR.sof", arm="NIR")

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert result["missing calibrations"].tolist() == ["master flat, spatial solution, flux standard"]


def test_treats_a_failed_upstream_calibration_as_missing_and_a_null_status_as_present() -> None:
    """A `fail` row does not satisfy the completeness test, and a NULL status does."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof")
    _add_calibration(connection, "mbias", "NOD_A.sof", "pass")
    _add_calibration(connection, "mflat", "NOD_A.sof", "fail")
    _add_calibration(connection, "spat_solution", "NOD_A.sof", None)
    _add_calibration(connection, "std_flux", "NOD_A.sof", "pass")

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert result["missing calibrations"].tolist() == ["master flat"]


def test_ignores_calibration_rows_that_belong_to_another_sof() -> None:
    """A calibration present for one set does not hide its absence from another."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof", obsId=1)
    _add_set(connection, "NOD_B.sof", obsId=2)
    for calibrationType in ("mbias", "mflat", "spat_solution", "std_flux"):
        _add_calibration(connection, calibrationType, "NOD_A.sof", "pass")

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    byObsId = dict(zip(result["eso obs id"], result["missing calibrations"], strict=True))
    assert byObsId == {
        1: "unknown",
        2: "master bias, master flat, spatial solution, flux standard",
    }


def test_merges_sofs_sharing_a_display_row_into_one_ordered_deduplicated_description() -> None:
    """Two sofs that print as one row give one union, in yaml order, with no repeats."""
    connection = _connection()
    _add_set(connection, "NOD_B.sof")
    _add_set(connection, "NOD_A.sof")
    # B LACKS MFLAT AND STD_FLUX; A LACKS SPAT_SOLUTION AND STD_FLUX
    for calibrationType in ("mbias", "spat_solution"):
        _add_calibration(connection, calibrationType, "NOD_B.sof", "pass")
    for calibrationType in ("mbias", "mflat"):
        _add_calibration(connection, calibrationType, "NOD_A.sof", "pass")

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert len(result) == 1
    assert result["missing calibrations"].tolist() == ["master flat, spatial solution, flux standard"]


def test_keeps_sets_with_different_display_values_on_separate_rows() -> None:
    """Grouping only merges sofs whose seven display columns are identical."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof", obsId=1)
    _add_set(connection, "NOD_B.sof", obsId=2)

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert sorted(result["eso obs id"]) == [1, 2]


def test_returns_an_empty_frame_with_the_new_column_when_nothing_is_incomplete() -> None:
    """A workspace with only complete sets yields no rows and no column error."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof", complete=1)
    _add_set(connection, "MBIAS.sof", recipe="mbias")

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert result.empty
    assert list(result.columns) == [*DISPLAY_COLUMNS, "missing calibrations"]


def test_falls_back_to_the_raw_name_for_an_unknown_calibration_type() -> None:
    """A calibration type with no description is reported by its raw name."""
    connection = _connection(("mbias", "novel_cal"))
    _add_set(connection, "NOD_A.sof")
    sofMapLookup = {"science": {"recipe": "nod_obj", "calibrations": {"VIS": ["mbias", "novel_cal"]}}}
    _add_calibration(connection, "mbias", "NOD_A.sof", "pass")

    result = _organiser(connection, sofMapLookup).get_incomplete_raw_frames_set()
    connection.close()

    assert "novel_cal" not in CALIBRATION_DESCRIPTIONS
    assert result["missing calibrations"].tolist() == ["novel_cal"]


@pytest.mark.parametrize(
    "sofMapLookup",
    [
        None,
        {},
        {"science": {"recipe": "stare_obj", "calibrations": {"VIS": ["mbias"]}}},
        {"science": {"recipe": "nod_obj", "calibrations": {"NIR": ["mdark"]}}},
        {"science": {"recipe": "nod_obj"}},
    ],
    ids=["no-attribute", "empty-map", "other-recipe", "other-arm", "no-calibrations"],
)
def test_reports_unknown_when_no_requirement_can_be_found(sofMapLookup: dict[str, Any] | None) -> None:
    """A missing map, recipe, arm or calibrations entry degrades to `unknown` rather than raising."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof")

    result = _organiser(connection, sofMapLookup).get_incomplete_raw_frames_set()
    connection.close()

    assert result["missing calibrations"].tolist() == ["unknown"]


def test_keeps_a_set_whose_display_column_is_null() -> None:
    """A NULL `slit` must not drop the set from the table, or the user is never told it is stuck."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof", slit=None)

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert len(result) == 1
    assert result["slit"].isna().tolist() == [True]


def test_groups_sofs_by_the_mjd_rounded_to_one_decimal() -> None:
    """Sofs whose `mjd-obs` round to the same tenth share a row, and those that do not stay apart."""
    connection = _connection()
    _add_set(connection, "NOD_A.sof", mjd=60000.04)
    _add_set(connection, "NOD_B.sof", mjd=60000.02)
    _add_set(connection, "NOD_C.sof", mjd=60000.06)

    result = _organiser(connection, _real_sof_map()).get_incomplete_raw_frames_set()
    connection.close()

    assert sorted(result["mjd-obs"]) == [60000.0, 60000.1]


def test_describes_every_calibration_type_the_sof_maps_can_require() -> None:
    """Every calibration type in the shipped sof maps has a label, so none falls back to a raw name."""
    requiredTypes = set()
    for mapName in ("soxs_sof_map.yaml", "xsh_sof_map.yaml"):
        with open(SOF_MAP_PATH.with_name(mapName)) as stream:
            sofMap = yaml.safe_load(stream) or {}
        for filters in sofMap.values():
            for calibrationTypes in (filters.get("calibrations") or {}).values():
                requiredTypes.update(calibrationTypes)

    assert requiredTypes
    assert requiredTypes <= set(CALIBRATION_DESCRIPTIONS)
