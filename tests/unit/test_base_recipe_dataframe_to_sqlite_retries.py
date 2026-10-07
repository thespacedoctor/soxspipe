"""Retry contract for `base_recipe._dataframe_to_sqlite` (DY-603).

The recipe-side insert retries `DataFrame.to_sql` up to seven times. If every
attempt fails it must re-raise the original error, so a lost `quality_control`
write fails the recipe instead of vanishing silently.
"""

from __future__ import annotations

import sqlite3
import time
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes.base_recipe import base_recipe

pytestmark = pytest.mark.unit

MAX_ATTEMPTS = 7


def _recipe(log: Any, connection: sqlite3.Connection) -> base_recipe:
    """Build a `base_recipe` carrying only the attributes the method reads."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.conn = connection
    return recipe


@pytest.fixture
def sleeps(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    """Record every `time.sleep` call instead of waiting."""
    calls: list[float] = []
    monkeypatch.setattr(time, "sleep", calls.append)
    return calls


@pytest.fixture
def to_sql_calls(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    """Count `DataFrame.to_sql` calls, failing the first `failures` of them with a locked-database error."""
    state: dict[str, Any] = {"attempts": 0, "failures": 0}
    originalToSql = pd.DataFrame.to_sql

    def counting_to_sql(self: pd.DataFrame, *args: Any, **kwargs: Any) -> Any:
        state["attempts"] += 1
        if state["attempts"] <= state["failures"]:
            raise sqlite3.OperationalError("database is locked")
        return originalToSql(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_sql", counting_to_sql)
    return state


def test_insert_that_fails_every_attempt_reraises_the_original_error(log, sleeps, to_sql_calls) -> None:
    # ARRANGE
    connection = sqlite3.connect(":memory:")
    recipe = _recipe(log, connection)
    to_sql_calls["failures"] = MAX_ATTEMPTS

    # ACT
    with pytest.raises(sqlite3.OperationalError, match="database is locked"):
        recipe._dataframe_to_sqlite(pd.DataFrame({"qc_name": ["RON"]}), "quality_control")

    # ASSERT
    assert to_sql_calls["attempts"] == MAX_ATTEMPTS
    # SLEEPS FALL BETWEEN ATTEMPTS ONLY, SO SEVEN ATTEMPTS MEAN SIX SLEEPS
    assert len(sleeps) == MAX_ATTEMPTS - 1
    connection.close()


def test_insert_on_a_closed_connection_raises_instead_of_returning(log, sleeps) -> None:
    # ARRANGE
    connection = sqlite3.connect(":memory:")
    connection.close()
    recipe = _recipe(log, connection)

    # ACT / ASSERT
    with pytest.raises(sqlite3.ProgrammingError):
        recipe._dataframe_to_sqlite(pd.DataFrame({"qc_name": ["RON"]}), "quality_control")
    assert len(sleeps) == MAX_ATTEMPTS - 1


def test_insert_that_succeeds_after_failures_writes_rows_and_stops_retrying(log, sleeps, to_sql_calls) -> None:
    # ARRANGE
    connection = sqlite3.connect(":memory:")
    connection.execute("create table quality_control (qc_name text, qc_value text)")
    recipe = _recipe(log, connection)
    to_sql_calls["failures"] = 3

    # ACT
    recipe._dataframe_to_sqlite(pd.DataFrame({"qc_name": ["RON"], "qc_value": ["5.0"]}), "quality_control")

    # ASSERT
    assert connection.execute("select qc_name, qc_value from quality_control").fetchall() == [("RON", "5.0")]
    assert to_sql_calls["attempts"] == 4
    assert len(sleeps) == 3
    connection.close()
