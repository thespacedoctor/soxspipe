"""Public finalization contracts for completed recipe runs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import base_recipe
from tests.factories import product_table, qc_table

pytestmark = pytest.mark.unit


@dataclass
class RecordingCursor:
    """Database cursor boundary that records recipe status statements."""

    statements: list[str]
    params: list[tuple[Any, ...] | None]
    wasClosed: bool = False

    def execute(self, statement: str, params: tuple[Any, ...] | None = None) -> None:
        """Record an executed database statement and its bound parameters."""
        self.statements.append(statement)
        self.params.append(params)

    def close(self) -> None:
        """Mark the cursor closed."""
        self.wasClosed = True


@dataclass
class RecordingConnection:
    """Database connection boundary that creates recording cursors."""

    statements: list[str] = field(default_factory=list)
    params: list[tuple[Any, ...] | None] = field(default_factory=list)
    cursors: list[RecordingCursor] = field(default_factory=list)
    wasClosed: bool = False

    def cursor(self) -> RecordingCursor:
        """Return a cursor that records status updates."""
        cursor = RecordingCursor(self.statements, self.params)
        self.cursors.append(cursor)
        return cursor

    def close(self) -> None:
        """Mark the connection closed."""
        self.wasClosed = True


def _recipe(*, log: Any, tmp_path: Path, qc: pd.DataFrame, status: str) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.qc = qc
    recipe.status = status
    recipe.currentSession = "20240102"
    recipe.sofName = "synthetic"
    recipe.workspaceRootPath = str(tmp_path)
    recipe.outDir = str(tmp_path / "intermediate")
    Path(recipe.outDir).mkdir()
    return recipe


def test_clean_up_marks_successful_recipe_complete_and_removes_intermediates(
    log: Any,
    tmp_path: Path,
) -> None:
    """A passing recipe commits its session status and removes intermediates."""
    recipe = _recipe(
        log=log,
        tmp_path=tmp_path,
        qc=pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]}),
        status="fail",
    )
    connection = RecordingConnection()
    recipe.conn = connection

    result = recipe.clean_up()

    assert result is None
    assert connection.statements == [
        "update product_frames set status_20240102 = 'pass' where sof = ?"
    ]
    assert connection.params == [("synthetic.sof",)]
    assert connection.cursors[0].wasClosed
    assert connection.wasClosed
    assert not Path(recipe.outDir).exists()
    assert not hasattr(recipe, "conn")


def test_clean_up_refreshes_session_when_failed_qc_reverses_passing_status(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A failed QC refreshes an earlier passing session as failed."""
    recipe = _recipe(
        log=log,
        tmp_path=tmp_path,
        qc=pd.DataFrame({"qc_flag": ["fail"], "qc_name": ["TRACE RMS"]}),
        status="pass",
    )
    connection = RecordingConnection()
    recipe.conn = connection
    organiserCalls: list[tuple[str, object]] = []

    class RecordingOrganiser:
        def __init__(self, **kwargs: object) -> None:
            organiserCalls.append(("init", kwargs))

        def session_refresh(self, *, failure: bool) -> None:
            organiserCalls.append(("session_refresh", failure))

        def close(self) -> None:
            organiserCalls.append(("close", None))

    monkeypatch.setattr("soxspipe.commonutils.data_organiser", RecordingOrganiser)

    recipe.clean_up()

    assert connection.statements == []
    assert connection.wasClosed
    assert organiserCalls == [
        ("init", {"log": log, "rootDir": str(tmp_path)}),
        ("session_refresh", True),
        ("close", None),
    ]
    assert ("error", "\nRecipe marked as failed in the database as the following QC values are outside of the acceptable limits: TRACE RMS.") in log.messages
    assert not Path(recipe.outDir).exists()


def test_clean_up_records_an_explicit_failure_reason(log: Any, tmp_path: Path) -> None:
    """An explicit recipe failure reason is stored before intermediates are removed."""
    recipe = _recipe(
        log=log,
        tmp_path=tmp_path,
        qc=pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]}),
        status="pass",
    )
    connection = RecordingConnection()
    recipe.conn = connection

    recipe.clean_up(forceFail="synthetic calibration mismatch")

    assert connection.statements == [
        "update product_frames set error_message = ? where sof = ?"
    ]
    assert connection.params == [("synthetic calibration mismatch", "synthetic.sof")]
    assert connection.wasClosed
    assert not Path(recipe.outDir).exists()
    assert (
        "error",
        "\nRecipe marked as failed in the database. synthetic calibration mismatch",
    ) in log.messages


def test_constructor_initializes_a_sof_recipe_in_an_isolated_workspace(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Initialize the public base-recipe state without touching a user session database."""
    import soxspipe.commonutils as commonutils
    import soxspipe.commonutils.toolkit as toolkit

    productPath = tmp_path / "reduced" / "synthetic.fits"
    organiserCalls: list[dict[str, object]] = []

    class IsolatedOrganiser:
        """Session lookup boundary with no external workspace state."""

        def __init__(self, **kwargs: object) -> None:
            organiserCalls.append(kwargs)

        def session_list(self, *, silent: bool) -> tuple[bool, list[str]]:
            assert silent is True
            return False, []

        def close(self) -> None:
            return None

    monkeypatch.setattr(commonutils, "data_organiser", IsolatedOrganiser)
    monkeypatch.setattr(
        toolkit,
        "predict_product_path",
        lambda *_: (str(productPath), "2024-01-02"),
    )
    monkeypatch.setattr(toolkit, "add_recipe_logger", lambda receivedLog, _: receivedLog)
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: "calibrations")
    monkeypatch.setattr(
        toolkit,
        "utility_setup",
        lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")),
    )

    recipe = base_recipe(
        log=log,
        settings={"workspace-root-dir": str(tmp_path), "instrument": "soxs"},
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
        command="soxspipe mbias synthetic.sof",
        turnOffMP=True,
    )

    assert recipe.sofName == "synthetic"
    assert recipe.productPath == str(productPath)
    assert recipe.startNightDate == "2024-01-02"
    assert recipe.calibrationRootPath == "calibrations"
    assert recipe.conn is None
    assert recipe.qc.empty
    assert recipe.products.empty
    assert recipe.kw("INSTRUME") == "INSTRUME"
    assert organiserCalls == [
        {"log": log, "rootDir": str(tmp_path), "dbConnect": False}
    ]


def test_input_frame_validation_records_consistent_vis_detector_metadata(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Accept a uniform VIS inventory and derive its detector measurement units."""
    import importlib

    from astropy.table import Table

    recipeModule = importlib.import_module("soxspipe.recipes.base_recipe")
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.kw = lambda name: name
    recipe.settings = {"instrument": "soxs"}
    recipe.recipeName = "soxs-stare"
    recipe.get_recipe_settings = lambda: {"frame-clipping-sigma": 3.0}
    summary = Table(
        {
            "SEQ_ARM": ["VIS"],
            "INSTRUME": ["SOXS"],
            "PRO_CATG": ["SCIENCE_VIS"],
            "CDELT1": [1],
            "CDELT2": [1],
            "DET_READ_SPEED": ["fast"],
            "CONAD": [1.1],
            "GAIN": [1.0],
            "DPR_TYPE": ["OBJECT"],
            "SLIT_VIS": ["SLIT1.0"],
            "RON": [3.5],
            "PRO_TYPE": [None],
            "DPR_TECH": ["ECHELLE,SLIT"],
            "PRO_TECH": [None],
            "DPR_CATG": ["SCIENCE"],
        }
    )

    class UniformFrames:
        """One-frame inventory exposing the ImageFileCollection query seam."""

        def __init__(self) -> None:
            self.summary = summary

        def files_filtered(self, *, include_path: bool) -> list[str]:
            assert include_path is True
            return ["science.fits"]

        def values(self, keyword: str, *, unique: bool) -> list[object]:
            assert unique is True
            return list(summary[keyword])

    class DetectorLookup:
        """Deterministic detector settings source for validation behavior."""

        def __init__(self, **_: object) -> None:
            return None

        def get(self, arm: str) -> dict[str, object]:
            assert arm == "VIS"
            return {"dispersion-axis": "x", "gain": 2.0, "ron": 4.0}

    recipe.inputFrames = UniformFrames()
    monkeypatch.setattr(recipeModule, "detector_lookup", DetectorLookup)

    imageTypes, imageTechnology, imageCategories = recipe._verify_input_frames_basics()

    assert recipe.arm == "VIS"
    assert recipe.inst == "SOXS"
    assert recipe.axisA == "x"
    assert recipe.axisB == "y"
    assert recipe.detectorParams["binning"] == [1, 1]
    assert recipe.detectorParams["gain"].value == pytest.approx(1.1)
    assert recipe.detectorParams["ron"].value == pytest.approx(3.5)
    assert set(imageTypes) == {"OBJECT"}
    assert set(imageTechnology) == {"ECHELLE,SLIT"}
    assert set(imageCategories) == {"SCIENCE", "SCIENCE_VIS"}


def test_report_output_formats_whole_number_qc_values_to_four_decimals(
    log: Any,
) -> None:
    """A qc_value column of whole numbers is numpy int64, and still prints formatted."""
    # ARRANGE
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.verbose = False
    recipe.conn = None
    recipe.sofName = "synthetic"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.recipeName = "soxs-mbias"
    recipe.inst = "XSHOOTER"
    recipe.recipeSettings = {}
    recipe.qc = qc_table().assign(qc_name="N ORDERS", qc_value=5)
    recipe.products = product_table()

    # ACT
    recipe.report_output()

    # ASSERT
    assert recipe.qc["qc_value"].dtype == "int64"
    assert any("5.0000" in message for level, message in log.messages)


def test_report_output_deduplicates_qc_and_hides_paths_in_console_mode(
    log: Any,
) -> None:
    """Console reporting retains the last QC value and does not expose file paths."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.verbose = False
    recipe.conn = None
    recipe.sofName = "synthetic"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.recipeName = "soxs-mbias"
    recipe.inst = "XSHOOTER"
    recipe.recipeSettings = {}
    recipe.qc = pd.concat(
        [
            qc_table(),
            qc_table().assign(qc_value=4.1, reduction_date_utc="2024-01-02T05:06:07"),
        ],
        ignore_index=True,
    )
    recipe.products = product_table()

    reported = recipe.report_output()

    assert reported["qc_value"].tolist() == [4.1]
    assert reported["sof_name"].tolist() == ["synthetic.sof"]
    assert "file_path" not in recipe.products.columns
    assert any("SOXS-MBIAS QC METRICS" in message for level, message in log.messages)
