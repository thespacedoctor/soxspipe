"""Characterization tests for recipe dispatch and bulk reduction bookkeeping."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytestmark = pytest.mark.unit


class RecordingCursor:
    """Record the SQL statements issued by bulk reduction bookkeeping."""

    def __init__(self, changedRows: int = 0) -> None:
        self.changedRows = changedRows
        self.executed: list[str] = []
        self.executedMany: list[tuple[str, list[tuple[str, str]]]] = []
        self.isClosed = False

    def execute(self, query: str) -> None:
        self.executed.append(query)

    def executemany(self, query: str, rows: list[tuple[str, str]]) -> None:
        self.executedMany.append((query, rows))

    def fetchone(self) -> tuple[int]:
        return (self.changedRows,)

    def close(self) -> None:
        self.isClosed = True


class RecordingConnection:
    """Provide the subset of the SQLite connection contract used by the reducer."""

    def __init__(self, changedRows: int = 0) -> None:
        self.recordingCursor = RecordingCursor(changedRows=changedRows)
        self.isClosed = False

    def cursor(self) -> RecordingCursor:
        return self.recordingCursor

    def close(self) -> None:
        self.isClosed = True


@pytest.mark.parametrize(
    ("recipeToken", "recipeAttribute"),
    [
        ("mbias", "soxs_mbias"),
        ("mdark", "soxs_mdark"),
        ("disp_solution", "soxs_disp_solution"),
        ("order_centres", "soxs_order_centres"),
        ("mflat", "soxs_mflat"),
        ("spat_solution", "soxs_spatial_solution"),
        ("stare_obj", "soxs_stare"),
        ("nod_std", "soxs_nod"),
        ("offset_obj", "soxs_offset"),
    ],
)
def test_run_recipe_dispatches_supported_tokens_with_unchanged_arguments(
    monkeypatch: pytest.MonkeyPatch,
    log: Any,
    recipeToken: str,
    recipeAttribute: str,
) -> None:
    """Map each reducer token to its public recipe constructor."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    recipesModule = import_module("soxspipe.recipes")
    constructorCalls: list[dict[str, object]] = []
    expectedQc = pd.DataFrame([{"qc_name": "signal", "qc_flag": "pass"}])

    class RecordingRecipe:
        def __init__(self, **kwargs: object) -> None:
            constructorCalls.append(kwargs)

        def produce_product(self) -> tuple[str, pd.DataFrame]:
            return "/workspace/product.fits", expectedQc

    monkeypatch.setattr(recipesModule, recipeAttribute, RecordingRecipe)

    result = reducerModule.run_recipe(
        log=log,
        recipe=recipeToken,
        sof="/workspace/input.sof",
        settings={"workspace-root-dir": "/workspace"},
        overwrite=True,
        command="soxspipe recipe input.sof",
        verbose=True,
        turnOffMP=True,
    )

    assert result == ("/workspace/product.fits", expectedQc)
    assert constructorCalls == [
        {
            "log": log,
            "settings": {"workspace-root-dir": "/workspace"},
            "inputFrames": "/workspace/input.sof",
            "overwrite": True,
            "command": "soxspipe recipe input.sof",
            "verbose": True,
            "turnOffMP": True,
        }
    ]


def test_run_recipe_preserves_unknown_token_failure(log: Any) -> None:
    """Keep the current failure type for unsupported recipe tokens."""
    reducerModule = import_module("soxspipe.commonutils.reducer")

    with pytest.raises(UnboundLocalError):
        reducerModule.run_recipe(
            log=log,
            recipe="unsupported",
            sof="input.sof",
            settings={},
            overwrite=False,
        )


def test_run_recipe_bulk_records_status_qc_errors_and_refreshes_changed_sofs(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Translate worker outcomes into session status, error, QC, and refresh writes."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    commonutilsModule = import_module("soxspipe.commonutils")
    fundamentalsModule = import_module("fundamentals")
    workerCalls: list[dict[str, object]] = []
    organiserCalls: list[tuple[str, object]] = []
    workspacePath = tmp_path / "workspace"
    temporaryPath = workspacePath / "tmp"
    temporaryPath.mkdir(parents=True)

    def fake_run_recipe(**kwargs: object) -> tuple[str, pd.DataFrame]:
        workerCalls.append(kwargs)
        sofName = Path(str(kwargs["sof"])).name
        if sofName == "good.sof":
            return "good.fits", pd.DataFrame([{"qc_name": "signal", "qc_flag": "pass"}])
        if sofName == "bad-qc.sof":
            return "bad.fits", pd.DataFrame(
                [
                    {"qc_name": "signal", "qc_flag": "fail"},
                    {"qc_name": "noise", "qc_flag": "fail"},
                ]
            )
        if sofName == "previous-fail.sof":
            raise FileExistsError("this recipe previously failed")
        if sofName == "previous-pass.sof":
            raise FileExistsError("the product of this recipe exists")
        if sofName == "other-exists.sof":
            raise FileExistsError("unclassified existing output")
        raise RuntimeError("worker exploded")

    def fake_multiprocess(**kwargs: object) -> list[dict[str, object]]:
        assert kwargs["poolSize"] == 3
        assert kwargs["timeout"] == 36000
        assert kwargs["wrapperTurnOffMP"] is True
        assert kwargs["turnOffMP"] is False
        wrapper = kwargs["function"]
        return [
            wrapper(
                inputDict=inputDict,
                log=kwargs["log"],
                recipe=kwargs["recipe"],
                settings=kwargs["settings"],
                overwrite=kwargs["overwrite"],
                workspaceDirectory=kwargs["workspaceDirectory"],
                wrapperTurnOffMP=kwargs["wrapperTurnOffMP"],
            )
            for inputDict in kwargs["inputArray"]
        ]

    class RecordingOrganiser:
        def __init__(self, **kwargs: object) -> None:
            organiserCalls.append(("init", kwargs))

        def session_refresh(self, **kwargs: object) -> bool:
            organiserCalls.append(("refresh", kwargs))
            return False

        def _dataframe_to_sqlite(
            self, dataframe: pd.DataFrame, tableName: str, replace: bool
        ) -> None:
            organiserCalls.append(("sqlite", (dataframe.copy(), tableName, replace)))

        def close(self) -> None:
            organiserCalls.append(("close", None))

    monkeypatch.setattr(reducerModule, "run_recipe", fake_run_recipe)
    monkeypatch.setattr(fundamentalsModule, "fmultiprocess", fake_multiprocess)
    monkeypatch.setattr(commonutilsModule, "data_organiser", RecordingOrganiser)
    connection = RecordingConnection(changedRows=1)
    sofNames = [
        "good.sof",
        "bad-qc.sof",
        "previous-fail.sof",
        "previous-pass.sof",
        "other-exists.sof",
        "exception.sof",
    ]
    sofPaths = [str(workspacePath / "sof" / sofName) for sofName in sofNames]

    result = reducerModule.run_recipe_bulk(
        log=log,
        recipe="mflat",
        sofList=sofPaths,
        commandList=[f"soxspipe mflat {sofName}" for sofName in sofNames],
        settings={"workspace-root-dir": str(workspacePath)},
        overwrite=True,
        workspaceDirectory=str(workspacePath),
        conn=connection,
        sessionId="night1",
    )

    assert result is None
    assert all(call["turnOffMP"] is True for call in workerCalls)
    assert connection.isClosed is True
    assert connection.recordingCursor.isClosed is True
    assert any(
        "status_night1 = 'pass'" in query
        for query in connection.recordingCursor.executed
    )
    assert any(
        "status_night1 = 'fail'" in query
        for query in connection.recordingCursor.executed
    )
    assert connection.recordingCursor.executedMany == [
        (
            "update product_frames set error_message = ? where sof = ?",
            [
                (
                    "The following QC values are outside of acceptable limits: signal, noise.",
                    "bad-qc.sof",
                ),
                ("unclassified existing output", "other-exists.sof"),
                ("worker exploded", "exception.sof"),
            ],
        )
    ]
    assert ("refresh", {"failure": True}) in organiserCalls
    sqliteWrites = [call for call in organiserCalls if call[0] == "sqlite"]
    assert len(sqliteWrites) == 1
    writtenQc, tableName, replace = sqliteWrites[0][1]
    assert tableName == "quality_control"
    assert replace is False
    assert writtenQc["qc_name"].tolist() == [
        "signal",
        "signal",
        "signal",
        "noise",
    ]
    assert not temporaryPath.exists()
    output = capsys.readouterr().out
    assert "pool size of 3" in output
    assert "Number of successful mflat reductions: 1" in output
    assert "Number of failed mflat reductions: 3" in output
    assert "Number of pre-existing mflat reductions: 2" in output


def test_run_recipe_bulk_skips_refresh_and_qc_when_workers_return_nothing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """Avoid organizer writes when a bulk dispatch has no worker results."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    commonutilsModule = import_module("soxspipe.commonutils")
    fundamentalsModule = import_module("fundamentals")
    organiserWasCreated = False

    def fail_if_created(**kwargs: object) -> None:
        nonlocal organiserWasCreated
        organiserWasCreated = True

    monkeypatch.setattr(fundamentalsModule, "fmultiprocess", lambda **kwargs: [])
    monkeypatch.setattr(commonutilsModule, "data_organiser", fail_if_created)
    connection = RecordingConnection(changedRows=0)

    reducerModule.run_recipe_bulk(
        log=log,
        recipe="mbias",
        sofList=[],
        commandList=[],
        settings={},
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        conn=connection,
        sessionId="night1",
    )

    assert organiserWasCreated is False
    assert connection.recordingCursor.executedMany == []
