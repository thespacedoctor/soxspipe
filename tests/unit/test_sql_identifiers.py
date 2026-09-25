"""Unit tests for `soxspipe.commonutils.sql_identifiers`.

`validate_sql_identifier` is the boundary check for values that end up
interpolated into SQL text as a table or column name, because SQLite cannot
bind an identifier as a query parameter. These tests characterize the
grammar it enforces, and downstream tests in this file characterize the
callers that bind everything else -- the sof name and the failure message --
as parameters instead of interpolating them.
"""

from __future__ import annotations

import sqlite3
from typing import Any

import pytest

from soxspipe.commonutils import data_organiser as data_organiser_module
from soxspipe.commonutils.data_organiser import _validate_session_id
from soxspipe.commonutils.sql_identifiers import (
    UnsafeSqlIdentifierError,
    validate_sql_identifier,
)
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# `validate_sql_identifier` -- THE GRAMMAR ITSELF.
# ---------------------------------------------------------------------------


def test_validate_sql_identifier_returns_a_normal_name_unchanged() -> None:
    """A plain, short, alphanumeric-plus-underscore name is returned as-is."""
    # ARRANGE / ACT
    result = validate_sql_identifier("status_base", "column name")

    # ASSERT
    assert result == "status_base"


def test_validate_sql_identifier_rejects_an_empty_string() -> None:
    """An empty string never matches the grammar's `{0,63}` tail on a mandatory head character."""
    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        validate_sql_identifier("", "column name")


def test_validate_sql_identifier_rejects_a_non_string() -> None:
    """A non-`str` identifier is rejected outright, regardless of its value."""
    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        validate_sql_identifier(123, "table name")


@pytest.mark.parametrize(
    "hostileIdentifier",
    [
        "table'; drop table product_frames; --",
        "table;drop",
        "table name",
    ],
)
def test_validate_sql_identifier_rejects_names_containing_sql_metacharacters(
    hostileIdentifier: str,
) -> None:
    """A quote, semicolon or space in the candidate name is rejected."""
    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        validate_sql_identifier(hostileIdentifier, "table name")


def test_validate_sql_identifier_rejects_a_name_longer_than_64_characters() -> None:
    """The grammar caps identifiers at 64 characters: one head plus 63 tail characters."""
    # ARRANGE
    tooLong = "a" * 65

    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        validate_sql_identifier(tooLong, "table name")


def test_validate_sql_identifier_accepts_a_name_of_exactly_64_characters() -> None:
    """64 characters is the last accepted length, one below the rejection boundary above."""
    # ARRANGE
    atTheLimit = "a" * 64

    # ACT
    result = validate_sql_identifier(atTheLimit, "table name")

    # ASSERT
    assert result == atTheLimit


def test_validate_sql_identifier_rejects_a_digit_leading_name() -> None:
    """The head character must be a letter or an underscore.

    This is why the callers validate the composed `status_<session>` column name
    rather than the session name alone, which is digit-leading by default.
    """
    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        validate_sql_identifier("1table", "table name")


def test_validate_sql_identifier_rejects_a_non_ascii_name() -> None:
    """The grammar is ASCII-only: an accented letter is not a permitted identifier character."""
    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        validate_sql_identifier("café", "table name")


def test_validate_sql_identifier_error_names_the_label_and_the_rejected_value() -> None:
    """The raised error names both the caller-supplied label and the rejected value."""
    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError, match="table name") as excInfo:
        validate_sql_identifier("bad name", "table name")

    assert "bad name" in str(excInfo.value)


# ---------------------------------------------------------------------------
# `_validate_session_id` -- THE SESSION-ID GRAMMAR, NARROWER THAN THE GENERAL
# SQL-IDENTIFIER GRAMMAR ABOVE, AND A CHARACTERIZATION OF HOW IT COMPOSES WITH
# `validate_sql_identifier` ONCE PREFIXED WITH ``status_``.
# ---------------------------------------------------------------------------


def test_validate_session_id_rejects_a_hyphen() -> None:
    """A hyphenated session id is outside the documented grammar."""
    # ACT / ASSERT
    with pytest.raises(ValueError, match="Session ID") as excinfo:
        _validate_session_id("my-session")

    assert "_-" not in str(excinfo.value)


@pytest.mark.parametrize("sessionId", ["x" * 17, ""])
def test_validate_session_id_rejects_length_boundary_violations(
    sessionId: str,
) -> None:
    """A session id past the 16-character limit, or empty, is rejected directly by the grammar."""
    # ACT / ASSERT
    with pytest.raises(ValueError, match="Session ID"):
        _validate_session_id(sessionId)


@pytest.mark.parametrize(
    "sessionId",
    [
        "base",
        "my_supernova",
        "20260920t143000",
        "_",
        "9",
        "A" * 16,
        "z9_" * 5 + "x",
    ],
)
def test_every_accepted_session_id_composes_a_valid_status_column(
    sessionId: str,
) -> None:
    """Every session id the grammar accepts also composes a valid `status_<id>` column name."""
    # ACT / ASSERT
    assert _validate_session_id(sessionId) == sessionId
    assert (
        validate_sql_identifier(f"status_{sessionId}", "status column")
        == f"status_{sessionId}"
    )


# ---------------------------------------------------------------------------
# `base_recipe.clean_up` -- HOSTILE VALUES REACH THE DATABASE AS BOUND
# PARAMETERS, NOT AS INTERPOLATED SQL TEXT.
# ---------------------------------------------------------------------------


class _RecordingCursor:
    """A stub cursor that records every `execute()` call verbatim."""

    def __init__(self, executedQueries: list[tuple[str, tuple[Any, ...] | None]]) -> None:
        self.executedQueries = executedQueries

    def execute(self, sqlQuery: str, params: tuple[Any, ...] | None = None) -> None:
        self.executedQueries.append((sqlQuery, params))

    def close(self) -> None:
        return None


class _RecordingConnection:
    """A stub DB connection that hands out recording cursors."""

    def __init__(self) -> None:
        self.executedQueries: list[tuple[str, tuple[Any, ...] | None]] = []

    def cursor(self) -> _RecordingCursor:
        return _RecordingCursor(self.executedQueries)

    def close(self) -> None:
        return None


def _clean_up_recipe(
    log: Any,
    tmp_path: Any,
    *,
    sofName: str,
    status: str = "fail",
    currentSession: str = "base",
) -> tuple[base_recipe, _RecordingConnection]:
    """Build a recipe carrying only the attributes `clean_up` reads, plus its stub connection."""
    import pandas as pd

    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]})
    recipe.status = status
    recipe.currentSession = currentSession
    recipe.sofName = sofName
    conn = _RecordingConnection()
    recipe.conn = conn
    recipe.outDir = str(tmp_path / "already-removed")
    return recipe, conn


def test_clean_up_binds_an_injection_shaped_sof_name_as_a_parameter(
    log: Any,
    tmp_path: Any,
) -> None:
    """A sof name shaped like a SQL injection payload never reaches the SQL text itself."""
    # ARRANGE
    hostileSofName = "x' OR 1=1 --"
    recipe, conn = _clean_up_recipe(log, tmp_path, sofName=hostileSofName, status="fail")

    # ACT
    recipe.clean_up(forceFail=False)

    # ASSERT
    assert len(conn.executedQueries) == 1
    sqlQuery, params = conn.executedQueries[0]
    assert hostileSofName not in sqlQuery
    assert params == (f"{hostileSofName}.sof",)


def test_clean_up_binds_an_injection_shaped_force_fail_message_as_a_parameter(
    log: Any,
    tmp_path: Any,
) -> None:
    """A `forceFail` message containing a single quote never reaches the SQL text itself."""
    # ARRANGE
    hostileMessage = "boom: it's broken"
    recipe, conn = _clean_up_recipe(log, tmp_path, sofName="synthetic", status="fail")

    # ACT
    recipe.clean_up(forceFail=hostileMessage)

    # ASSERT
    assert len(conn.executedQueries) == 1
    sqlQuery, params = conn.executedQueries[0]
    assert hostileMessage not in sqlQuery
    assert params == (hostileMessage, "synthetic.sof")


def test_clean_up_accepts_the_default_digit_leading_session_id(
    log: Any,
    tmp_path: Any,
) -> None:
    """A digit-leading session id, the shape `data_organiser` generates by default, is accepted.

    `validate_sql_identifier`'s grammar requires a leading letter or
    underscore, so `clean_up` validates the composed `status_<session>`
    column name -- guaranteed letter-led by its literal `status_` prefix --
    rather than the bare session fragment, which is not.
    """
    # ARRANGE
    defaultShapedSessionId = "20260920t143000"
    recipe, conn = _clean_up_recipe(
        log,
        tmp_path,
        sofName="synthetic",
        status="pass",
        currentSession=defaultShapedSessionId,
    )

    # ACT
    recipe.clean_up(forceFail=False)

    # ASSERT
    assert len(conn.executedQueries) == 1
    sqlQuery, params = conn.executedQueries[0]
    assert sqlQuery == f"update product_frames set status_{defaultShapedSessionId} = 'pass' where sof = ?"  # noqa: S608
    assert params == ("synthetic.sof",)


# ---------------------------------------------------------------------------
# `_dataframe_to_sqlite` -- A HOSTILE TABLE NAME IS REJECTED BEFORE IT REACHES
# SQL TEXT, IN BOTH `base_recipe` AND `data_organiser`.
# ---------------------------------------------------------------------------


def test_base_recipe_dataframe_to_sqlite_rejects_a_hostile_table_name(log: Any) -> None:
    """A table name shaped like a SQL injection payload is rejected before any query runs."""
    # ARRANGE
    import pandas as pd

    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.conn = sqlite3.connect(":memory:")
    dataframe = pd.DataFrame({"qc_name": ["RON"]})

    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        recipe._dataframe_to_sqlite(dataframe, "quality_control; drop table product_frames; --")

    recipe.conn.close()


def test_data_organiser_dataframe_to_sqlite_rejects_a_hostile_table_name(log: Any) -> None:
    """A table name shaped like a SQL injection payload is rejected before any query runs."""
    # ARRANGE
    import pandas as pd

    organiser = data_organiser_module.__new__(data_organiser_module)
    organiser.log = log
    organiser.conn = sqlite3.connect(":memory:")
    dataframe = pd.DataFrame({"file": ["synthetic.fits"]})

    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        organiser._dataframe_to_sqlite(dataframe, "safe_rows; drop table raw_frames; --")

    organiser.conn.close()


# ---------------------------------------------------------------------------
# END-TO-END: A HOSTILE SOF NAME UPDATES ZERO ROWS AGAINST A REAL SQLITE
# `product_frames` TABLE, NOT EVERY ROW.
# ---------------------------------------------------------------------------


def test_clean_up_with_hostile_sof_name_updates_no_rows_in_a_real_database(
    log: Any,
    tmp_path: Any,
) -> None:
    """An injection-shaped sof name matches no real row, leaving every existing row untouched."""
    # ARRANGE
    import pandas as pd

    # A FILE-BACKED DATABASE, NOT `:memory:`, SO THE ROWS SURVIVE THE CONNECTION
    # THAT `clean_up` CLOSES AND CAN BE RE-OPENED TO VERIFY THEM.
    # `autocommit=True` MATCHES THE CONNECTION `base_recipe.__init__` OPENS IN
    # PRODUCTION -- WITHOUT IT, SQLITE'S IMPLICIT TRANSACTION IS NEVER
    # COMMITTED AND THE INJECTED UPDATE WOULD ROLL BACK ON CLOSE REGARDLESS OF
    # WHETHER THE VALUE WAS BOUND OR INTERPOLATED, MASKING THE BEHAVIOUR UNDER TEST.
    databasePath = tmp_path / "verify.db"
    connection = sqlite3.connect(str(databasePath), autocommit=True)
    connection.execute("create table product_frames (sof text, status_base text)")
    connection.executemany(
        "insert into product_frames values (?, ?)",
        [("one.sof", "fail"), ("two.sof", "fail")],
    )
    connection.commit()

    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]})
    recipe.status = "fail"
    recipe.currentSession = "base"
    # AN INJECTION PAYLOAD THAT, IF STILL INTERPOLATED, CLOSES THE STRING
    # LITERAL EARLY, FORCES AN ALWAYS-TRUE CONDITION, AND COMMENTS OUT THE
    # TRAILING `.sof'` SO EVERY ROW WOULD BE UPDATED TO 'pass'.
    recipe.sofName = "x' OR '1'='1' --"
    recipe.conn = connection
    recipe.outDir = str(tmp_path / "already-removed")

    # ACT
    recipe.clean_up(forceFail=False)

    # ASSERT
    verifyConnection = sqlite3.connect(str(databasePath))
    statuses = verifyConnection.execute(
        "select sof, status_base from product_frames order by sof"
    ).fetchall()
    verifyConnection.close()
    # NEITHER ROW MATCHES THE BOUND PAYLOAD, SO BOTH STAY 'fail'.
    assert statuses == [("one.sof", "fail"), ("two.sof", "fail")]


def test_clean_up_with_hostile_force_fail_message_writes_no_error_in_a_real_database(
    log: Any,
    tmp_path: Any,
) -> None:
    """An injection-shaped failure message matches no real row, so no row gains an error message."""
    # ARRANGE
    import pandas as pd

    # FILE-BACKED AND `autocommit=True`, FOR THE SAME REASONS AS THE SOF-NAME
    # TEST ABOVE.
    databasePath = tmp_path / "verify-force-fail.db"
    connection = sqlite3.connect(str(databasePath), autocommit=True)
    connection.execute("create table product_frames (sof text, status_base text, error_message text)")
    connection.executemany(
        "insert into product_frames values (?, ?, ?)",
        [("one.sof", "pass", None), ("two.sof", "pass", None)],
    )
    connection.commit()

    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]})
    # 'fail' RATHER THAN 'pass', SO `clean_up` TAKES THE `error_message` BRANCH
    # WITHOUT ALSO TAKING THE `data_organiser` SESSION-REFRESH BRANCH.
    recipe.status = "fail"
    recipe.currentSession = "base"
    recipe.sofName = "one"
    recipe.conn = connection
    recipe.outDir = str(tmp_path / "already-removed")
    # A MESSAGE THAT, IF STILL INTERPOLATED, WOULD CLOSE THE STRING LITERAL,
    # SUBSTITUTE AN ALWAYS-TRUE `where` CLAUSE, AND COMMENT OUT THE REAL ONE,
    # STAMPING EVERY ROW. THE TRAILING `--` MATTERS: WITHOUT IT THE INJECTED
    # STATEMENT IS A SYNTAX ERROR, AND THIS TEST WOULD FAIL ON THE OLD CODE FOR
    # THE WRONG REASON.
    hostileMessage = "boom' where sof like '%' --"

    # ACT
    recipe.clean_up(forceFail=hostileMessage)

    # ASSERT
    verifyConnection = sqlite3.connect(str(databasePath))
    rows = verifyConnection.execute(
        "select sof, error_message from product_frames order by sof"
    ).fetchall()
    verifyConnection.close()
    # ONLY THE MATCHING ROW IS STAMPED, AND THE PAYLOAD IS STORED VERBATIM AS A
    # VALUE RATHER THAN EXECUTED AS SQL.
    assert rows == [("one.sof", hostileMessage), ("two.sof", None)]


# ---------------------------------------------------------------------------
# `data_organiser._qc_acceptable_range_queries` -- WORKSPACE-SETTINGS RECIPE
# AND QC-RANGE KEYS REACH THE DATABASE AS BOUND PARAMETERS, NOT INTERPOLATED
# SQL TEXT.
# ---------------------------------------------------------------------------


def _qc_range_organiser(*, sessionId: str = "base", settings: dict[str, Any]) -> Any:
    """Build a `data_organiser` carrying only the attributes the method reads."""
    organiser = data_organiser_module.__new__(data_organiser_module)
    organiser.sessionId = sessionId
    organiser.settings = settings
    return organiser


def test_qc_acceptable_range_queries_binds_a_hostile_recipe_key_as_a_parameter() -> None:
    """A workspace-settings recipe key shaped like a SQL injection payload is bound, not interpolated."""
    # ARRANGE
    hostileRecipe = 'soxs-x" OR 1=1 --'
    organiser = _qc_range_organiser(
        settings={hostileRecipe: {"qc-acceptable-ranges": {"ron": [1.0, 2.0]}}}
    )

    # ACT
    queries = organiser._qc_acceptable_range_queries()

    # ASSERT
    matching = [q for q in queries if q[1] and q[1][2] == hostileRecipe]
    assert len(matching) == 1
    sqlQuery, params = matching[0]
    assert hostileRecipe not in sqlQuery
    assert params == (1.0, 2.0, hostileRecipe, "RON")


def test_qc_acceptable_range_queries_updates_no_row_for_a_hostile_recipe_key_in_a_real_database() -> None:
    """An injection-shaped recipe key matches no real row, leaving every existing row untouched."""
    # ARRANGE
    hostileRecipe = 'soxs-x" OR 1=1 --'
    organiser = _qc_range_organiser(
        settings={hostileRecipe: {"qc-acceptable-ranges": {"ron": [1.0, 2.0]}}}
    )
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "create table quality_control (soxspipe_recipe text, qc_name text, "
        "qc_order text, qc_value_min real, qc_value_max real, qc_flag text, "
        "qc_value text, sof_name text)"
    )
    connection.executemany(
        "insert into quality_control values (?, ?, ?, ?, ?, ?, ?, ?)",
        [
            ("soxs-mbias", "RON", "-1", None, None, "pass", "1.5", "one.sof"),
            ("soxs-mflat", "RON", "-1", None, None, "pass", "1.5", "two.sof"),
        ],
    )
    connection.execute("create table product_frames (sof text, status_base text)")
    connection.commit()

    # ACT
    for sqlQuery, sqlParams in organiser._qc_acceptable_range_queries():
        connection.execute(sqlQuery, sqlParams)
    connection.commit()

    # ASSERT
    rows = connection.execute(
        "select soxspipe_recipe, qc_value_min, qc_value_max from quality_control order by soxspipe_recipe"
    ).fetchall()
    connection.close()
    # NEITHER REAL ROW MATCHES THE HOSTILE RECIPE, SO NEITHER GAINS THE
    # HOSTILE RANGE. IF THE PAYLOAD WERE STILL INTERPOLATED, THE `OR 1=1`
    # WOULD HAVE SET THE RANGE ON BOTH.
    assert rows == [
        ("soxs-mbias", None, None),
        ("soxs-mflat", None, None),
    ]


def test_qc_acceptable_range_queries_updates_no_row_for_a_hostile_qc_key_in_a_real_database() -> None:
    """An injection-shaped QC-range key matches no real row, leaving every existing row untouched.

    The recipe name is the other value read from the same workspace settings
    and is covered by the test above; this covers the QC-range key itself
    (the `kk` in `v["qc-acceptable-ranges"].items()`), which becomes `qc_name`.
    """
    # ARRANGE
    # `/*` RATHER THAN `--`: `qc_name` IS BUILT VIA `kk.upper().replace("-", " ")`,
    # WHICH TURNS A TRAILING `--` COMMENT MARKER INTO TWO SPACES AND ONLY BREAKS
    # THE STILL-VULNERABLE SQL WITH A SYNTAX ERROR RATHER THAN SILENTLY
    # CORRUPTING BOTH ROWS. `/*` SURVIVES THE REPLACEMENT AND REPRODUCES THE
    # REAL DOUBLE-ROW CORRUPTION AGAINST THE OLD, UNPARAMETERIZED QUERY.
    hostileQcKey = 'ron" OR 1=1 /*'
    organiser = _qc_range_organiser(
        settings={"soxs-mbias": {"qc-acceptable-ranges": {hostileQcKey: [1.0, 2.0]}}}
    )
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "create table quality_control (soxspipe_recipe text, qc_name text, "
        "qc_order text, qc_value_min real, qc_value_max real, qc_flag text, "
        "qc_value text, sof_name text)"
    )
    connection.executemany(
        "insert into quality_control values (?, ?, ?, ?, ?, ?, ?, ?)",
        [
            ("soxs-mbias", "RON", "-1", None, None, "pass", "1.5", "one.sof"),
            ("soxs-mbias", "FLUX", "-1", None, None, "pass", "1.5", "two.sof"),
        ],
    )
    connection.execute("create table product_frames (sof text, status_base text)")
    connection.commit()

    # ACT
    for sqlQuery, sqlParams in organiser._qc_acceptable_range_queries():
        connection.execute(sqlQuery, sqlParams)
    connection.commit()

    # ASSERT
    rows = connection.execute(
        "select soxspipe_recipe, qc_name, qc_value_min, qc_value_max from quality_control order by qc_name"
    ).fetchall()
    connection.close()
    # NEITHER REAL ROW MATCHES THE HOSTILE QC NAME, SO NEITHER GAINS THE
    # HOSTILE RANGE. IF THE PAYLOAD WERE STILL INTERPOLATED, THE `OR 1=1 /*`
    # WOULD HAVE SET THE RANGE ON BOTH -- CONFIRMED BY RECONSTRUCTING THE OLD
    # UNPARAMETERIZED QUERY WITH THIS EXACT PAYLOAD AGAINST THIS SAME FIXTURE.
    assert rows == [
        ("soxs-mbias", "FLUX", None, None),
        ("soxs-mbias", "RON", None, None),
    ]
