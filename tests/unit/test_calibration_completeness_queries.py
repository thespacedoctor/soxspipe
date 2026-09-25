"""Unit tests for `data_organiser`'s calibration-completeness query builders.

DY-254: `build_sof_files`'s calibration-completeness cluster interpolated
`recipe`, `arm` and `ttype` as SQL text values, and a calibration type (`ct`)
as a `cal_<type>` table-name identifier, several times over. It also joined
a list of product SOFs unbound into an `IN (...)` clause. The cluster was
extracted into three small query-builder methods -- mirroring
`_qc_acceptable_range_queries` in `tests/unit/test_sql_identifiers.py` -- so
each concern is independently testable here.
"""

from __future__ import annotations

import sqlite3
from typing import Any

import pytest

from soxspipe.commonutils import data_organiser as data_organiser_module
from soxspipe.commonutils.sql_identifiers import UnsafeSqlIdentifierError

pytestmark = pytest.mark.unit


def _organiser() -> Any:
    return data_organiser_module.__new__(data_organiser_module)


# ---------------------------------------------------------------------------
# `_calibration_completeness_select_query` -- `recipe`, `arm` AND `ttype`
# REACH SQLITE AS BOUND PARAMETERS; `calType` ENTRIES ARE VALIDATED
# IDENTIFIERS.
# ---------------------------------------------------------------------------


def test_select_query_rejects_a_hostile_calibration_type() -> None:
    """A calibration type shaped like a SQL injection payload is rejected before any query runs."""
    organiser = _organiser()

    with pytest.raises(UnsafeSqlIdentifierError):
        organiser._calibration_completeness_select_query(
            "mbias", "UVB", "BIAS", ["bias; drop table product_frames; --"]
        )


def test_select_query_binds_recipe_and_arm_as_parameters() -> None:
    """`recipe` and `arm` are bound placeholders, never interpolated into the SQL text."""
    organiser = _organiser()
    hostileRecipe = "mflat' OR '1'='1' -- "
    hostileArm = "UVB' OR '1'='1' -- "

    sqlQuery, sqlParams = organiser._calibration_completeness_select_query(
        hostileRecipe, hostileArm, "FLAT", ["bias"]
    )

    assert hostileRecipe not in sqlQuery
    assert hostileArm not in sqlQuery
    assert sqlParams == (hostileRecipe, hostileArm)


def test_select_query_adds_a_bound_dpr_type_filter_only_for_std_types() -> None:
    """The extra `eso dpr type` filter, and its parameter, only appear for a `STD` type."""
    organiser = _organiser()

    sqlQuery, sqlParams = organiser._calibration_completeness_select_query(
        "mbias", "UVB", "STD,FLUX", ["bias"]
    )

    assert '"eso dpr type" = ?' in sqlQuery
    assert sqlParams == ("mbias", "UVB", "STD,FLUX")

    sqlQueryNonStd, sqlParamsNonStd = organiser._calibration_completeness_select_query(
        "mbias", "UVB", "BIAS", ["bias"]
    )

    assert '"eso dpr type"' not in sqlQueryNonStd
    assert sqlParamsNonStd == ("mbias", "UVB")


def test_select_query_matches_only_the_real_row_against_a_real_database() -> None:
    """A hostile recipe/arm pair matches no real row, rather than widening the query.

    IF STILL INTERPOLATED, THE PAYLOADS BELOW CLOSE THEIR STRING LITERALS,
    INJECT AN ALWAYS-TRUE CONDITION, AND COMMENT OUT THE REST OF THE QUERY
    (INCLUDING THE `EXISTS (...)` CALIBRATION CHECK), MATCHING BOTH ROWS
    INSTEAD OF NEITHER.
    """
    connection = sqlite3.connect(":memory:")
    connection.execute(
        'CREATE TABLE product_frames_plus (sof TEXT, recipe TEXT, complete INTEGER, '
        '"eso seq arm" TEXT, "eso dpr type" TEXT)'
    )
    connection.execute("CREATE TABLE cal_bias (sof TEXT, upstream_status TEXT)")
    connection.executemany(
        "INSERT INTO product_frames_plus VALUES (?, ?, ?, ?, ?)",
        [
            ("target.sof", "mflat", 0, "UVB", "FLAT"),
            ("other.sof", "mflat", 0, "VIS", "FLAT"),
        ],
    )
    connection.executemany(
        "INSERT INTO cal_bias VALUES (?, ?)",
        [("target.sof", "pass"), ("other.sof", "pass")],
    )
    connection.commit()
    organiser = _organiser()

    hostileRecipe = "mflat' OR '1'='1' -- "
    sqlQuery, sqlParams = organiser._calibration_completeness_select_query(
        hostileRecipe, "UVB", "FLAT", ["bias"]
    )
    matches = connection.execute(sqlQuery, sqlParams).fetchall()
    connection.close()

    assert matches == []


# ---------------------------------------------------------------------------
# `_calibration_completeness_update_query` -- `containerSofs` REACH SQLITE AS
# BOUND PARAMETERS, NEVER JOINED INTO THE SQL TEXT.
# ---------------------------------------------------------------------------


def test_update_query_matches_no_row_for_an_empty_container_list() -> None:
    """An empty `containerSofs` list produces a clause that matches no row."""
    organiser = _organiser()

    sqlQuery, sqlParams = organiser._calibration_completeness_update_query([])

    assert sqlParams == ()
    assert "sof in (NULL)" in sqlQuery


def test_update_query_binds_one_placeholder_per_sof_in_call_order() -> None:
    """Each `containerSofs` entry gets its own `?` placeholder, bound in the same order it was given."""
    organiser = _organiser()

    sqlQuery, sqlParams = organiser._calibration_completeness_update_query(
        ["first.sof", "second.sof", "third.sof"]
    )

    assert sqlQuery.count("?") == 3
    assert "sof in (?, ?, ?)" in sqlQuery
    assert sqlParams == ("first.sof", "second.sof", "third.sof")


def test_update_query_marks_only_the_hostile_sof_itself_complete_in_a_real_database() -> None:
    """A hostile SOF name in `containerSofs` updates only the row that literally matches it.

    IF STILL JOINED AS RAW TEXT, THIS PAYLOAD CLOSES THE STRING LITERAL,
    INJECTS AN ALWAYS-TRUE CONDITION, AND MARKS EVERY ROW COMPLETE -- NOT
    JUST THE ONE (NON-EXISTENT) ROW ACTUALLY REQUESTED.
    """
    connection = sqlite3.connect(":memory:")
    connection.execute("CREATE TABLE product_frames (sof TEXT, complete INTEGER)")
    connection.executemany(
        "INSERT INTO product_frames VALUES (?, ?)",
        [("one.sof", 0), ("two.sof", 0)],
    )
    connection.commit()
    organiser = _organiser()
    # THE OLD CODE WRAPPED EACH ENTRY IN DOUBLE QUOTES (`sof in ("...","...")`).
    # THIS PAYLOAD CLOSES THAT DOUBLE-QUOTED STRING, INJECTS AN ALWAYS-TRUE
    # CONDITION, AND REOPENS A DOUBLE-QUOTED STRING TO BALANCE THE TEMPLATE'S
    # OWN TRAILING `")`, MARKING BOTH REAL ROWS COMPLETE INSTEAD OF NEITHER.
    hostileSof = 'x") OR complete < 1 OR ("x'

    sqlQuery, sqlParams = organiser._calibration_completeness_update_query([hostileSof])
    connection.execute(sqlQuery, sqlParams)
    connection.commit()

    rows = connection.execute(
        "select sof, complete from product_frames order by sof"
    ).fetchall()
    connection.close()

    assert rows == [("one.sof", 0), ("two.sof", 0)]


# ---------------------------------------------------------------------------
# `_calibration_raw_frames_query` -- A CALIBRATION TYPE IS A VALIDATED
# IDENTIFIER, NEVER AN UNCHECKED TABLE-NAME INTERPOLATION.
# ---------------------------------------------------------------------------


def test_raw_frames_query_rejects_a_hostile_calibration_type() -> None:
    """A calibration type shaped like a SQL injection payload is rejected before any query runs."""
    organiser = _organiser()

    with pytest.raises(UnsafeSqlIdentifierError):
        organiser._calibration_raw_frames_query("bias; drop table product_frames; --")


def test_raw_frames_query_reads_from_the_composed_calibration_table_in_a_real_database() -> None:
    """A valid calibration type composes the expected `cal_<type>` table reference."""
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "CREATE TABLE product_frames (sof TEXT, complete INTEGER)"
    )
    connection.execute(
        "CREATE TABLE cal_bias (sof TEXT, file TEXT, upstream_tag TEXT, filepath TEXT)"
    )
    connection.execute("INSERT INTO product_frames VALUES ('one.sof', -1)")
    connection.execute(
        "INSERT INTO cal_bias VALUES ('one.sof', 'one.fits', 'MASTER_BIAS', './cal/one.fits')"
    )
    connection.commit()
    organiser = _organiser()

    sqlQuery = organiser._calibration_raw_frames_query("bias")
    rows = connection.execute(sqlQuery).fetchall()
    connection.close()

    assert sqlQuery.count("cal_bias") == 5
    assert "cal_bad" not in sqlQuery
    assert rows == [("one.fits", "MASTER_BIAS", "one.sof", "./cal/one.fits", -1)]
