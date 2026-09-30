"""Unit contracts for the late-frame rejoin helpers (DY-263)."""

from __future__ import annotations

import sqlite3

import pandas as pd
import pytest

from soxspipe.commonutils.late_frame_rejoin import (
    delete_stale_files,
    find_rejoined_sofs,
    format_rejoin_summary,
    pin_rejoined_sof_names,
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
