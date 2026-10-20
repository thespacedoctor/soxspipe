"""Unit contracts for the late-frame rejoin and session regroup helpers (DY-263, DY-274)."""

from __future__ import annotations

import sqlite3

import pandas as pd
import pytest

from soxspipe.commonutils.late_frame_rejoin import (
    collect_sof_name_pins,
    delete_stale_files,
    find_rejoined_sofs,
    format_rejoin_summary,
    partition_labels,
    pin_rejoined_sof_names,
    read_set_members,
    session_map_tables,
)

pytestmark = pytest.mark.unit

OLD_SOF = "20240102T060411_MFLAT.sof"
NEW_SOF = "20240102T060416_MFLAT.sof"
MEMBERS = frozenset({"a.fits", "b.fits"})


def _groups(recipe="mflat", filepaths=("a.fits", "b.fits", "c.fits")):
    return pd.DataFrame([{"sof": NEW_SOF, "recipe": recipe, "filepaths": list(filepaths)}])


def test_pin_gives_the_old_name_to_a_group_that_holds_every_old_member() -> None:
    # ARRANGE
    groups = _groups()

    # ACT
    pinned = pin_rejoined_sof_names(groups, {OLD_SOF: ("mflat", MEMBERS)})

    # ASSERT
    assert pinned["sof"].tolist() == [OLD_SOF]


def test_pin_does_not_change_the_frame_it_is_given() -> None:
    # ARRANGE
    groups = _groups()

    # ACT
    pin_rejoined_sof_names(groups, {OLD_SOF: ("mflat", MEMBERS)})

    # ASSERT
    assert groups["sof"].tolist() == [NEW_SOF]


def test_pin_keeps_the_new_name_when_the_group_lacks_an_old_member() -> None:
    # ARRANGE
    groups = _groups(filepaths=("a.fits", "c.fits"))

    # ACT
    pinned = pin_rejoined_sof_names(groups, {OLD_SOF: ("mflat", MEMBERS)})

    # ASSERT
    assert pinned["sof"].tolist() == [NEW_SOF]


def test_pin_keeps_the_new_name_when_the_recipe_differs() -> None:
    # ARRANGE
    groups = _groups(recipe="order_centres")

    # ACT
    pinned = pin_rejoined_sof_names(groups, {OLD_SOF: ("mflat", MEMBERS)})

    # ASSERT
    assert pinned["sof"].tolist() == [NEW_SOF]


def test_pin_passes_an_empty_frame_through() -> None:
    # ARRANGE
    groups = pd.DataFrame()

    # ACT
    pinned = pin_rejoined_sof_names(groups, {OLD_SOF: ("mflat", MEMBERS)})

    # ASSERT
    assert pinned.empty


def test_summary_names_every_regrouped_sof_in_order() -> None:
    # ARRANGE
    sofs = [NEW_SOF, OLD_SOF]

    # ACT
    summary = format_rejoin_summary(sofs)

    # ASSERT
    assert summary.startswith("2 SOF(S)")
    assert summary.index(OLD_SOF) < summary.index(NEW_SOF)


def _database_with_a_late_frame(arm, tech, dprType="LAMP,FLAT"):
    """A database holding one processed member and one late frame of the same set, both with the given arm and tech."""
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "CREATE TABLE raw_frames_valid (filepath TEXT, processed INTEGER, `eso tpl start` TEXT, "
        "`eso seq arm` TEXT, `eso dpr tech` TEXT, `eso dpr type` TEXT)"
    )
    connection.execute("CREATE TABLE sof_map_test (file TEXT, tag TEXT, sof TEXT, filepath TEXT, complete INT)")
    connection.execute("CREATE TABLE product_frames (sof TEXT, recipe TEXT)")
    connection.executemany(
        "INSERT INTO raw_frames_valid VALUES (?, ?, '2024-01-02T06:04:00', ?, ?, ?)",
        [("a.fits", 1, arm, tech, dprType), ("b.fits", 0, arm, tech, dprType)],
    )
    connection.execute("INSERT INTO sof_map_test VALUES ('a.fits', 'T', ?, 'a.fits', 1)", (OLD_SOF,))
    connection.execute("INSERT INTO product_frames VALUES (?, 'mflat')", (OLD_SOF,))
    return connection


@pytest.mark.parametrize("arm, tech", [(None, "IMAGE"), ("VIS", None), ("VIS", "IMAGE"), ("NIR", None)])
def test_a_null_arm_or_tech_does_not_hide_a_late_frame_from_its_set(arm, tech) -> None:
    # ARRANGE
    connection = _database_with_a_late_frame(arm, tech)

    # ACT
    rejoined = find_rejoined_sofs(connection, "sof_map_test", ["eso seq arm", "eso dpr tech", "eso dpr type"])

    # ASSERT
    connection.close()
    assert rejoined == {OLD_SOF: ("mflat", frozenset({"a.fits"}))}


def test_delete_stale_files_removes_files_and_warns_about_one_it_cannot_remove(tmp_path, log) -> None:
    # ARRANGE
    removable = tmp_path / "product.fits"
    removable.write_text("stale")
    unremovable = tmp_path / "a_directory"
    unremovable.mkdir()
    missing = tmp_path / "gone.fits"

    # ACT
    delete_stale_files([removable, unremovable, missing], log)

    # ASSERT
    warnings = [message for level, message in log.messages if level == "warning"]
    assert not removable.exists()
    assert len(warnings) == 1 and "a_directory" in warnings[0]


def test_a_nir_image_frame_with_a_null_type_is_hidden_like_any_other_off_frame() -> None:
    # ARRANGE
    connection = _database_with_a_late_frame("NIR", "IMAGE", dprType=None)

    # ACT
    rejoined = find_rejoined_sofs(connection, "sof_map_test", ["eso seq arm", "eso dpr tech", "eso dpr type"])

    # ASSERT
    connection.close()
    assert rejoined == {}


# SESSION REGROUP HELPERS (DY-274)


def _workspace_with_two_sessions():
    """A database where `alpha` holds a set of two raw frames and a calibration row, and `beta` holds one raw frame."""
    connection = sqlite3.connect(":memory:")
    connection.execute("CREATE TABLE raw_frames (filepath TEXT)")
    connection.executemany("INSERT INTO raw_frames VALUES (?)", [("a.fits",), ("b.fits",), ("c.fits",)])
    connection.execute("CREATE TABLE product_frames (sof TEXT, recipe TEXT)")
    connection.executemany("INSERT INTO product_frames VALUES (?, ?)", [(OLD_SOF, "mflat"), ("NO_RAW.sof", "mflat")])
    for table in ("sof_map_alpha", "sof_map_beta", "z_sof_map", "sof_map_not_a_session_name_that_is_far_too_long"):
        connection.execute(f"CREATE TABLE {table} (file TEXT, tag TEXT, sof TEXT, filepath TEXT, complete INT)")  # noqa: S608
    connection.executemany(
        "INSERT INTO sof_map_alpha (sof, filepath) VALUES (?, ?)",
        [(OLD_SOF, "a.fits"), (OLD_SOF, "b.fits"), (OLD_SOF, "./reduced/master_bias.fits"), ("NO_RAW.sof", "x.fits")],
    )
    connection.execute("INSERT INTO sof_map_beta (sof, filepath) VALUES (?, ?)", (OLD_SOF, "c.fits"))
    return connection


def test_session_map_tables_lists_only_the_tables_that_sessions_own() -> None:
    # ARRANGE
    connection = _workspace_with_two_sessions()

    # ACT
    tables = session_map_tables(connection)

    # ASSERT
    connection.close()
    assert tables == ["sof_map_alpha", "sof_map_beta"]


def test_read_set_members_leaves_out_calibration_rows_and_sets_without_a_recipe() -> None:
    # ARRANGE
    connection = _workspace_with_two_sessions()

    # ACT
    members = read_set_members(connection, "sof_map_alpha")

    # ASSERT
    connection.close()
    assert members == {OLD_SOF: ("mflat", frozenset({"a.fits", "b.fits"}))}


def test_collect_name_pins_joins_the_frames_a_set_has_in_each_session() -> None:
    # ARRANGE
    connection = _workspace_with_two_sessions()

    # ACT
    pins = collect_sof_name_pins(connection, "sof_map_beta", {})

    # ASSERT
    connection.close()
    assert pins == {OLD_SOF: ("mflat", frozenset({"a.fits", "b.fits", "c.fits"}))}


def test_collect_name_pins_keeps_the_pinned_names_first() -> None:
    # ARRANGE
    connection = _workspace_with_two_sessions()
    pinned = {NEW_SOF: ("mflat", frozenset({"z.fits"}))}

    # ACT
    pins = collect_sof_name_pins(connection, "sof_map_alpha", pinned)

    # ASSERT
    connection.close()
    assert list(pins) == [NEW_SOF, OLD_SOF]


def test_partition_labels_name_the_known_set_of_each_frame_of_one_recipe() -> None:
    # ARRANGE
    pins = {
        OLD_SOF: ("mflat", frozenset({"a.fits", "b.fits"})),
        NEW_SOF: ("order_centres", frozenset({"c.fits"})),
    }

    # ACT
    labels = partition_labels(pins, "mflat")

    # ASSERT
    assert labels == {"a.fits": OLD_SOF, "b.fits": OLD_SOF}


def test_partition_labels_leave_out_the_excluded_sets() -> None:
    # ARRANGE
    pins = {OLD_SOF: ("mflat", MEMBERS), NEW_SOF: ("mflat", frozenset({"c.fits"}))}

    # ACT
    labels = partition_labels(pins, "mflat", excludedSofs=[OLD_SOF])

    # ASSERT
    assert labels == {"c.fits": NEW_SOF}


def test_find_rejoined_sofs_does_not_read_the_frames_when_the_sof_map_is_empty() -> None:
    # ARRANGE: THERE IS NO `raw_frames_valid` TABLE, SO READING THE FRAMES WOULD RAISE
    connection = sqlite3.connect(":memory:")
    connection.execute("CREATE TABLE sof_map_test (file TEXT, tag TEXT, sof TEXT, filepath TEXT, complete INT)")

    # ACT
    rejoined = find_rejoined_sofs(connection, "sof_map_test", ["eso seq arm"])

    # ASSERT
    connection.close()
    assert rejoined == {}
