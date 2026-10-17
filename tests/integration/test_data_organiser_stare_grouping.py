"""Integration contracts for grouping stare science frames into one SOF per OB and arm (DY-1611)."""

from __future__ import annotations

import sqlite3

import pandas as pd
import pytest

from soxspipe.commonutils.late_frame_rejoin import find_rejoined_sofs
from tests.factories import raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration

STARE_TECH = "ECHELLE,SLIT,STARE"


def _stare_frames(arm, obsId, count, firstSecond=0, exptime=300.0, tplStart="2024-01-02T03:00:00"):
    """Return ``count`` stare object exposures of one OB on one arm."""
    template = raw_frame_table().iloc[0].to_dict()
    rows = []
    for index in range(count):
        second = firstSecond + index
        name = f"stare-{arm}-{obsId}-{second}.fits"
        rows.append(
            {
                **template,
                "eso seq arm": arm,
                "eso dpr catg": "SCIENCE",
                "eso dpr tech": STARE_TECH,
                "eso dpr type": "OBJECT",
                "eso obs id": obsId,
                "eso obs name": f"OB {obsId}",
                "exptime": exptime,
                "slit": "1.0",
                "object": "SN2024abc",
                "template": "SOXS_slt_obs_Stare",
                "eso tpl name": "SOXS_slt_obs_Stare",
                "eso tpl nexp": count,
                "eso tpl expno": index + 1,
                "eso tpl start": tplStart,
                "file": name,
                "filepath": f"./raw/2024-01-01/{name}",
                "date-obs": f"2024-01-02T03:05:{second:02d}.000",
                "mjd-obs": 60311.2118 + second * 0.0001,
                "set_first_file": name,
            }
        )
    return pd.DataFrame(rows)


def _stare_groups(tmp_path, log, frames):
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.instrument = "SOXS"
    organiser.conn = sqlite3.connect(":memory:")
    frames.to_sql("raw_frames_valid", organiser.conn, index=False)
    _, rawGroups = organiser.get_raw_frames_and_groups(
        ttype="OBJECT",
        tech=[STARE_TECH],
        recipe="stare_obj",
        recipeOrder=0,
        filterName="science_stare",
        unprocessedOnly=True,
    )
    organiser.conn.close()
    return rawGroups


def _members_by_arm(rawGroups):
    return {
        row["eso seq arm"]: sorted(path.split("/")[-1] for path in row["filepaths"]) for _, row in rawGroups.iterrows()
    }


def test_one_stare_ob_of_three_frames_per_arm_gives_one_sof_per_arm(tmp_path, log) -> None:
    # ARRANGE
    frames = pd.concat([_stare_frames("VIS", "101", 3), _stare_frames("NIR", "101", 3)], ignore_index=True)

    # ACT
    rawGroups = _stare_groups(tmp_path, log, frames)

    # ASSERT
    assert len(rawGroups.index) == 2
    assert rawGroups["counts"].tolist() == [3, 3]
    assert _members_by_arm(rawGroups) == {
        "NIR": ["stare-NIR-101-0.fits", "stare-NIR-101-1.fits", "stare-NIR-101-2.fits"],
        "VIS": ["stare-VIS-101-0.fits", "stare-VIS-101-1.fits", "stare-VIS-101-2.fits"],
    }
    # THE SOF IS NAMED FROM THE EARLIEST FRAME OF THE OB
    assert set(rawGroups["sof"]) == {
        "20240102T030500_VIS_1X1_FAST_STARE_OBJ_1_0_300_0S_SOXS_SN2024abc.sof",
        "20240102T030500_NIR_1X1_FAST_STARE_OBJ_1_0_300_0S_SOXS_SN2024abc.sof",
    }


def test_two_stare_obs_give_two_sofs(tmp_path, log) -> None:
    # ARRANGE
    frames = pd.concat(
        [_stare_frames("VIS", "101", 2), _stare_frames("VIS", "102", 2, firstSecond=30)], ignore_index=True
    )

    # ACT
    rawGroups = _stare_groups(tmp_path, log, frames)

    # ASSERT
    assert sorted(rawGroups["eso obs id"].tolist()) == ["101", "102"]
    assert rawGroups["counts"].tolist() == [2, 2]


def test_stare_frames_with_a_different_exposure_time_in_one_ob_give_separate_sofs(tmp_path, log) -> None:
    # ARRANGE
    frames = pd.concat(
        [_stare_frames("VIS", "101", 2), _stare_frames("VIS", "101", 1, firstSecond=30, exptime=600.0)],
        ignore_index=True,
    )

    # ACT
    rawGroups = _stare_groups(tmp_path, log, frames)

    # ASSERT
    assert sorted(zip(rawGroups["exptime"], rawGroups["counts"], strict=True)) == [(300.0, 2), (600.0, 1)]


def test_a_late_stare_frame_rejoins_the_sof_of_its_ob(tmp_path, log) -> None:
    # ARRANGE: TWO FRAMES OF ONE OB ARE ALREADY IN A SOF; THE THIRD ARRIVES LATE
    organiser = workspace_organiser(tmp_path, log=log)
    frames = _stare_frames("VIS", "101", 3)
    frames["eso tpl nexp"] = 3
    frames["processed"] = [1, 1, 0]
    sof = "20240102T030500_VIS_1X1_FAST_STARE_OBJ_1_0_300_0S_SOXS_SN2024abc.sof"
    connection = sqlite3.connect(":memory:")
    frames.to_sql("raw_frames_valid", connection, index=False)
    pd.DataFrame({"sof": [sof, sof], "filepath": frames["filepath"].tolist()[:2]}).to_sql(
        "sof_map_science", connection, index=False
    )
    pd.DataFrame({"sof": [sof], "recipe": ["soxs-stare"]}).to_sql("product_frames", connection, index=False)
    groupingKeys = [key for key in organiser.filterKeywords if key not in organiser.proKeywords]

    # ACT
    rejoined = find_rejoined_sofs(connection, "sof_map_science", groupingKeys)

    # ASSERT
    assert rejoined == {sof: ("soxs-stare", frozenset(frames["filepath"].tolist()[:2]))}


def test_a_late_stare_frame_from_a_later_template_of_the_ob_rejoins_its_sof(tmp_path, log) -> None:
    """Stare frames group by OB, not by template run, so the rejoin must not demand the same template start."""
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    frames = pd.concat(
        [
            _stare_frames("VIS", "101", 2, tplStart="2024-01-02T03:00:00"),
            _stare_frames("VIS", "101", 1, firstSecond=30, tplStart="2024-01-02T03:30:00"),
        ],
        ignore_index=True,
    )
    frames["eso tpl nexp"] = 3
    frames["processed"] = [1, 1, 0]
    sof = "20240102T030500_VIS_1X1_FAST_STARE_OBJ_1_0_300_0S_SOXS_SN2024abc.sof"
    connection = sqlite3.connect(":memory:")
    frames.to_sql("raw_frames_valid", connection, index=False)
    pd.DataFrame({"sof": [sof, sof], "filepath": frames["filepath"].tolist()[:2]}).to_sql(
        "sof_map_science", connection, index=False
    )
    pd.DataFrame({"sof": [sof], "recipe": ["soxs-stare"]}).to_sql("product_frames", connection, index=False)
    groupingKeys = [key for key in organiser.filterKeywords if key not in organiser.proKeywords]

    # ACT
    rejoined = find_rejoined_sofs(connection, "sof_map_science", groupingKeys)

    # ASSERT
    assert rejoined == {sof: ("soxs-stare", frozenset(frames["filepath"].tolist()[:2]))}
