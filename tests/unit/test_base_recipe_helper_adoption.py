"""Characterization of the base recipe QC rows the shared helpers replace.

These tests pin what `base_recipe` does today at the five QC-table row
builders and the timestamp sites feeding them, so the helper adoption in
`toolkit.append_qc` can be proved to change nothing. They assert three things
the helper API could silently break: the column order of the resulting table,
the shared reduction timestamp across rows built in one call, and the rendered
timestamp format.

The timestamp assertions are deliberately structural rather than frozen. The
base revision mints its timestamp from a function-local `datetime.utcnow()`
and the converted revision from `toolkit.utcnow_string()`, so a test that
patches one module's name passes against only one of the two. Matching the
format and requiring one shared value pins the contract under both.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.recipes import base_recipe
from tests.factories import qc_table, synthetic_ccd

pytestmark = pytest.mark.unit

# THE TIMESTAMP FORMAT EVERY QC ROW IN THIS MODULE RENDERS: SECONDS
# RESOLUTION, NO FRACTIONAL PART AND NO UTC OFFSET MARKER.
QC_TIMESTAMP_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$")

# THE COLUMN CONTRACT `base_recipe.__init__` CREATES `self.qc` WITH. EVERY
# INLINE ROW BUILDER APPENDS TO A TABLE THAT ALREADY CARRIES THESE COLUMNS, SO
# THE ROW DICTIONARY'S OWN KEY ORDER NEVER REACHES THE TABLE.
QC_COLUMNS = [
    "soxspipe_recipe",
    "qc_name",
    "qc_value",
    "qc_unit",
    "qc_order",
    "qc_comment",
    "obs_date_utc",
    "reduction_date_utc",
    "to_header",
]


def _recipe(log: Any) -> base_recipe:
    """Build a base recipe with only the attributes the QC methods read."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.recipeName = "soxs-mbias"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.qc = pd.DataFrame({column: [] for column in QC_COLUMNS})
    return recipe


def test_init_builds_the_qc_table_in_the_expected_column_order() -> None:
    """`__init__` fixes the column order every QC row in this module appends into.

    Read from the source rather than from a constructed recipe: `__init__`
    needs a workspace and a settings file, and the only thing under test here
    is the order of the keys in its `self.qc` literal. That order is why the
    inline rows can adopt `append_qc` even though the helper builds its own
    keys in a different order -- `pd.concat` aligns to the columns the table
    already has.
    """
    # ARRANGE
    import ast
    import inspect

    source = inspect.getsource(base_recipe.__init__)
    tree = ast.parse(inspect.cleandoc(source))

    # ACT
    qcAssignments = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and isinstance(node.targets[0], ast.Attribute)
        and node.targets[0].attr == "qc"
    ]

    # ASSERT
    # ASSERT THE COUNT BEFORE INDEXING, SO A LATER `__init__` THAT BUILDS ITS
    # QC TABLE SOME OTHER WAY FAILS HERE WITH A READABLE MESSAGE INSTEAD OF AN
    # `IndexError` FROM THE LINE BELOW.
    assert len(qcAssignments) == 1
    assert [key.value for key in qcAssignments[0].value.args[0].keys] == QC_COLUMNS


def test_flag_poor_data_shares_one_timestamp_across_both_temperature_rows(
    log: Any,
) -> None:
    """Both SOXS temperature rows record the same reduction timestamp."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.inst = "SOXS"
    recipe.detectorTemp = 82.5
    recipe.cptemp = 11.2
    recipe.recipeSettings = {}

    # ACT
    recipe.flag_poor_data()

    # ASSERT
    rows = recipe.qc.set_index("qc_name")
    detector = rows.loc["DETECTOR TEMP"]
    commonPath = rows.loc["CPATH TEMP"]
    assert detector["reduction_date_utc"] == commonPath["reduction_date_utc"]
    assert QC_TIMESTAMP_PATTERN.match(detector["reduction_date_utc"])


def _counting_clock(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Replace the shared clock with one that returns a new value every call.

    Comparing two rendered second-resolution timestamps cannot tell "one
    timestamp shared by both rows" from "two timestamps minted a few
    microseconds apart", because both render the same string. This stub makes
    the two cases distinguishable: every call returns a different value, so
    rows built from one call still match and rows built from two no longer do.

    Returns the list of values handed out, so a test can assert how many times
    the clock was read.
    """
    from soxspipe.commonutils import toolkit

    issued: list[str] = []

    def _clock(microseconds: bool = False) -> str:
        issued.append(f"2024-01-02T03:04:{len(issued):02d}")
        return issued[-1]

    monkeypatch.setattr(toolkit, "utcnow_string", _clock)
    return issued


def test_flag_poor_data_reads_the_clock_once_for_both_temperature_rows(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both temperature rows come from a single clock read, not one each."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.inst = "SOXS"
    recipe.detectorTemp = 82.5
    recipe.cptemp = 11.2
    recipe.recipeSettings = {}
    issued = _counting_clock(monkeypatch)

    # ACT
    recipe.flag_poor_data()

    # ASSERT
    assert len(issued) == 1
    assert list(recipe.qc["reduction_date_utc"]) == [issued[0], issued[0]]


def test_qc_ron_reads_the_clock_once_for_both_noise_rows(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both RON rows come from a single clock read, not one each."""
    # ARRANGE
    recipe = _recipe(log)
    issued = _counting_clock(monkeypatch)

    # ACT
    recipe.qc_ron(frameType="MBIAS", rawRon=2.5, masterRon=0.8)

    # ASSERT
    assert len(issued) == 1
    assert list(recipe.qc["reduction_date_utc"]) == [issued[0], issued[0]]


def test_flag_poor_data_temperature_rows_keep_their_values_and_columns(
    log: Any,
) -> None:
    """The two temperature rows are appended verbatim, in the table's column order."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.inst = "SOXS"
    recipe.detectorTemp = 82.5
    recipe.cptemp = 11.2
    recipe.recipeSettings = {}

    # ACT
    recipe.flag_poor_data()

    # ASSERT
    # `flag_poor_data` ADDS `qc_flag` OF ITS OWN, AFTER THE ROWS ARE BUILT.
    assert list(recipe.qc.columns) == [*QC_COLUMNS, "qc_flag"]
    rows = recipe.qc.set_index("qc_name")
    detector = rows.loc["DETECTOR TEMP"]
    assert detector["soxspipe_recipe"] == "soxs-mbias"
    assert detector["qc_value"] == pytest.approx(82.5)
    assert detector["qc_comment"] == "[K] temp of detector"
    assert detector["qc_unit"] == "kelvin"
    assert detector["obs_date_utc"] == "2024-01-02T03:04:05"
    assert bool(detector["to_header"]) is False
    commonPath = rows.loc["CPATH TEMP"]
    assert commonPath["qc_value"] == pytest.approx(11.2)
    assert commonPath["qc_comment"] == "[C] temp of common path"
    assert commonPath["qc_unit"] == "celsius"
    assert bool(commonPath["to_header"]) is False
    # NEITHER ROW SETS `qc_order`, SO `flag_poor_data` BACKFILLS IT TO "-1".
    assert set(recipe.qc["qc_order"]) == {"-1"}


def test_flag_poor_data_adds_no_rows_for_a_non_soxs_instrument(log: Any) -> None:
    """The temperature rows are SOXS-only; X-shooter gets neither."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.inst = "XSHOOTER"
    recipe.recipeSettings = {}
    recipe.qc = qc_table()

    # ACT
    recipe.flag_poor_data()

    # ASSERT
    assert list(recipe.qc["qc_name"]) == ["RON"]


def test_qc_ron_shares_one_timestamp_across_both_noise_rows(log: Any) -> None:
    """The raw and master RON rows record the same reduction timestamp."""
    # ARRANGE
    recipe = _recipe(log)

    # ACT
    recipe.qc_ron(
        frameType="MBIAS",
        frameName="master bias",
        rawRon=2.5,
        masterRon=0.8,
    )

    # ASSERT
    rows = recipe.qc.set_index("qc_name")
    rawRow = rows.loc["RAW RON"]
    masterRow = rows.loc["MASTER RON"]
    assert rawRow["reduction_date_utc"] == masterRow["reduction_date_utc"]
    assert QC_TIMESTAMP_PATTERN.match(rawRow["reduction_date_utc"])


def test_qc_ron_rows_keep_their_values_and_columns(log: Any) -> None:
    """Both RON rows are appended verbatim, in the table's column order."""
    # ARRANGE
    recipe = _recipe(log)

    # ACT
    rawRon, masterRon = recipe.qc_ron(
        frameType="MBIAS",
        frameName="master bias",
        rawRon=2.5,
        masterRon=0.8,
    )

    # ASSERT
    assert rawRon == pytest.approx(2.5, rel=1e-12)
    assert masterRon == pytest.approx(0.8, rel=1e-12)
    assert list(recipe.qc.columns) == QC_COLUMNS
    rows = recipe.qc.set_index("qc_name")
    rawRow = rows.loc["RAW RON"]
    assert rawRow["soxspipe_recipe"] == "soxs-mbias"
    assert rawRow["qc_value"] == pytest.approx(2.5, rel=1e-12)
    # THE LEADING "M" IS STRIPPED FROM THE FRAME TYPE FOR THE SINGLE-FRAME ROW.
    assert rawRow["qc_comment"] == "[e-] RON in single BIAS"
    assert rawRow["qc_unit"] == "electrons"
    assert rawRow["obs_date_utc"] == "2024-01-02T03:04:05"
    assert bool(rawRow["to_header"]) is True
    masterRow = rows.loc["MASTER RON"]
    assert masterRow["qc_value"] == pytest.approx(0.8, rel=1e-12)
    assert masterRow["qc_comment"] == "[e-] Combined RON in MBIAS"
    assert masterRow["qc_unit"] == "electrons"
    assert bool(masterRow["to_header"]) is True


def test_qc_ron_appends_only_the_raw_row_when_no_master_noise_is_supplied(
    log: Any,
) -> None:
    """Without a master frame or master RON, only the raw row is appended."""
    # ARRANGE
    recipe = _recipe(log)

    # ACT
    rawRon, masterRon = recipe.qc_ron(frameType="MBIAS", rawRon=2.5)

    # ASSERT
    assert rawRon == pytest.approx(2.5, rel=1e-12)
    assert masterRon is None
    assert list(recipe.qc["qc_name"]) == ["RAW RON"]


def test_qc_median_flux_level_row_keeps_its_values_and_columns(log: Any) -> None:
    """The median-flux row is appended verbatim, in the table's column order."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.kw = lambda keyword: keyword
    frame = synthetic_ccd(
        shape=(2, 2),
        prepared=True,
        headerOverrides={"EXPTIME": 2.0},
    )
    frame.data = np.array([[1.0, 3.0], [5.0, 1000.0]])
    frame.mask = np.array([[False, False], [False, True]])

    # ACT
    medianFlux = recipe.qc_median_flux_level(
        frame,
        frameType="MDARK",
        frameName="synthetic dark",
    )

    # ASSERT
    assert medianFlux == pytest.approx(3.0, rel=1e-12)
    assert list(recipe.qc.columns) == QC_COLUMNS
    row = recipe.qc.iloc[-1]
    assert row["soxspipe_recipe"] == "soxs-mbias"
    assert row["qc_name"] == "MDARK MEDIAN"
    assert row["qc_value"] == pytest.approx(3.0, rel=1e-12)
    assert row["qc_comment"] == "[e-] Median flux level of synthetic dark"
    assert row["qc_unit"] == "electrons"
    assert row["obs_date_utc"] == "2024-01-02T03:04:05"
    assert bool(row["to_header"]) is True
    assert QC_TIMESTAMP_PATTERN.match(row["reduction_date_utc"])


def test_qc_median_flux_level_uppercases_the_frame_type_in_the_metric_name(
    log: Any,
) -> None:
    """A lower-case frame type still produces an upper-case QC name."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.kw = lambda keyword: keyword
    frame = synthetic_ccd(
        shape=(2, 2),
        prepared=True,
        headerOverrides={"EXPTIME": 2.0},
    )
    frame.data = np.array([[1.0, 3.0], [5.0, 7.0]])
    frame.mask = np.zeros((2, 2), dtype=bool)

    # ACT
    recipe.qc_median_flux_level(frame, frameType="mflat", frameName="synthetic flat")

    # ASSERT
    assert recipe.qc.iloc[-1]["qc_name"] == "MFLAT MEDIAN"


def test_appending_every_qc_row_leaves_the_column_order_untouched(log: Any) -> None:
    """Three separate QC methods append to one table without reordering it."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.inst = "SOXS"
    recipe.detectorTemp = 82.5
    recipe.cptemp = 11.2
    recipe.recipeSettings = {}
    recipe.kw = lambda keyword: keyword
    frame = synthetic_ccd(
        shape=(2, 2),
        prepared=True,
        headerOverrides={"EXPTIME": 2.0},
    )
    frame.data = np.array([[1.0, 3.0], [5.0, 7.0]])
    frame.mask = np.zeros((2, 2), dtype=bool)

    # ACT
    recipe.qc_ron(frameType="MBIAS", rawRon=2.5, masterRon=0.8)
    recipe.qc_median_flux_level(frame, frameType="MBIAS", frameName="master bias")
    recipe.flag_poor_data()

    # ASSERT
    assert list(recipe.qc.columns) == [*QC_COLUMNS, "qc_flag"]
    assert list(recipe.qc["qc_name"]) == [
        "RAW RON",
        "MASTER RON",
        "MBIAS MEDIAN",
        "DETECTOR TEMP",
        "CPATH TEMP",
    ]
    assert pd.api.types.is_numeric_dtype(recipe.qc["qc_value"])
