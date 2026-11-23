"""The data organiser routes GenericOffset frames and offset SKY frames to `soxs-offset` (DY-2355)."""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.commonutils.data_organiser import data_organiser
from tests.factories import insert_frames, raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration

FIXED_SKY_TECH = "ECHELLE,SLIT,OFFSET"
GENERIC_TECH = "ECHELLE,SLIT,GENOFFSET"
STARE_TECH = "ECHELLE,SLIT,STARE"
NOD_TECH = "ECHELLE,SLIT,NODDING"
INSTRUMENT = "SOXS"
TPL_START = "2025-12-01T02:43:49"


def _insert_frame(
    connection: sqlite3.Connection, name: str, *, arm: str = "VIS", dprType: str = "SKY", tech: str = FIXED_SKY_TECH
) -> None:
    """Seed one raw-frame row with the columns the ignore rules read."""
    connection.execute(
        """INSERT INTO raw_frames (
            instrume, file, "mjd-obs", "mjd-date", "date-obs", "night start mjd",
            "eso dpr catg", "eso dpr type", "eso dpr tech", exptime, object, filepath,
            "eso seq arm", slit
        ) VALUES ('SOXS', ?, 60311.1, '2024-01-02', '2024-01-02T03:04:05', 60310,
            'SCIENCE', ?, ?, 60, 'M 77', ?, ?, 'SLIT1.0')""",
        (name, dprType, tech, f"./raw/{name}", arm),
    )


def _offset_run(
    arm: str,
    tech: str,
    dprTypes: Sequence[str],
    *,
    tplStart: str = TPL_START,
    prefix: str = "off",
    firstSecond: int = 0,
    exptime: float = 120.0,
) -> pd.DataFrame:
    """Return one offset template run: a frame per entry of ``dprTypes``, in template order."""
    template = raw_frame_table().iloc[0].to_dict()
    rows = []
    for index, dprType in enumerate(dprTypes):
        second = firstSecond + index
        name = f"{prefix}-{arm}-{second}.fits"
        rows.append(
            {
                **template,
                "eso seq arm": arm,
                "eso dpr catg": "SCIENCE",
                "eso dpr tech": tech,
                "eso dpr type": dprType,
                "eso obs id": "101",
                "eso obs name": "OB 101",
                "exptime": exptime,
                "slit": "SLIT1.0",
                "object": "M 77",
                "template": "SOXS_M77_offset",
                "eso tpl name": "SOXS_M77_offset",
                "eso tpl nexp": len(dprTypes),
                "eso tpl expno": index + 1,
                "eso tpl start": tplStart,
                "file": name,
                "filepath": f"./raw/2025-11-30/{name}",
                "date-obs": f"2025-12-01T02:44:{second:02d}.000",
                "mjd-obs": 60980.1 + second * 0.0001,
                "set_first_file": name,
            }
        )
    return pd.DataFrame(rows)


@pytest.fixture
def organiser(tmp_path: Path, log: Any) -> Iterator[data_organiser]:
    """A SOXS organiser on an empty workspace database, with default settings."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn, _ = organiser._get_or_create_db_connection()
    organiser.instrument = INSTRUMENT
    organiser.settings = {}
    yield organiser
    organiser.conn.close()


def _ignored(organiser: data_organiser) -> set[str]:
    """Return the names of the frames flagged to ignore."""
    rows = organiser.conn.execute("SELECT file FROM raw_frames WHERE ignore = 1").fetchall()
    return {row[0] for row in rows}


@pytest.mark.parametrize("tech", [FIXED_SKY_TECH, GENERIC_TECH])
@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_sky_frame_of_either_offset_template_is_not_ignored(organiser, tech, arm) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "offset-sky.fits", arm=arm, tech=tech)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


@pytest.mark.parametrize(
    ("arm", "tech"),
    [
        ("VIS", STARE_TECH),
        ("VIS", NOD_TECH),
        ("VIS", "IMAGE"),
        ("ACQ", "IMAGE,OFFSET,POINT"),
        ("ACQ", FIXED_SKY_TECH),
    ],
)
def test_sky_frame_that_is_not_an_offset_spectrum_is_still_ignored(organiser, arm, tech) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "other-sky.fits", arm=arm, tech=tech)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"other-sky.fits"}


def test_offset_sky_frame_is_still_ignored_for_xshooter(organiser) -> None:
    # ARRANGE
    organiser.instrument = "XSH"
    _insert_frame(organiser.conn, "xsh-sky.fits", tech=FIXED_SKY_TECH)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"xsh-sky.fits"}


def _flag_as_the_old_rule_did(organiser: data_organiser, name: str) -> None:
    """Set `ignore = 1` on one frame, as a workspace prepared before DY-2355 holds it."""
    organiser.conn.execute("UPDATE raw_frames SET ignore = 1 WHERE file = ?", (name,))
    organiser.conn.commit()


@pytest.mark.parametrize("tech", [FIXED_SKY_TECH, GENERIC_TECH])
def test_offset_sky_frame_flagged_by_the_old_rule_is_kept_again(organiser, tech) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "old-sky.fits", tech=tech)
    _flag_as_the_old_rule_did(organiser, "old-sky.fits")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


@pytest.mark.parametrize(
    ("arm", "tech", "slit"),
    [
        ("ACQ", FIXED_SKY_TECH, "SLIT1.0"),
        ("VIS", FIXED_SKY_TECH, "UNDEFINED"),
        ("VIS", STARE_TECH, "SLIT1.0"),
    ],
)
def test_sky_frame_that_another_rule_flags_stays_ignored_after_the_reset(organiser, arm, tech, slit) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "other-rule.fits", arm=arm, tech=tech)
    organiser.conn.execute("UPDATE raw_frames SET slit = ?", (slit,))
    _flag_as_the_old_rule_did(organiser, "other-rule.fits")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"other-rule.fits"}


def test_simulated_offset_sky_frame_stays_ignored_after_the_reset(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "simulated-sky.fits")
    organiser.conn.execute("UPDATE raw_frames SET simulation = 1")
    _flag_as_the_old_rule_did(organiser, "simulated-sky.fits")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"simulated-sky.fits"}


def test_ignored_frame_that_is_not_an_offset_sky_frame_is_not_reset(organiser) -> None:
    # ARRANGE: AN OBJECT FRAME FLAGGED BY HAND (A FLAG NO RULE SETS) IS NOT RESET: ONLY OFFSET SKY FRAMES ARE
    _insert_frame(organiser.conn, "hand-flagged.fits", dprType="OBJECT")
    _flag_as_the_old_rule_did(organiser, "hand-flagged.fits")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"hand-flagged.fits"}


def test_ignored_offset_sky_frame_stays_ignored_for_xshooter(organiser) -> None:
    # ARRANGE
    organiser.instrument = "XSH"
    _insert_frame(organiser.conn, "xsh-sky.fits")
    _flag_as_the_old_rule_did(organiser, "xsh-sky.fits")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"xsh-sky.fits"}


def _offset_groups(
    tmp_path: Path, log: Any, frames: pd.DataFrame, *, ttype: str = "OBJECT", filterName: str = "science_offset"
) -> pd.DataFrame:
    """Group ``frames`` with the technique list of the named SOXS offset entry of the SOF map."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser._select_instrument(inst=INSTRUMENT)
    organiser.conn = sqlite3.connect(":memory:")
    frames.to_sql("raw_frames_valid", organiser.conn, index=False)
    entry = organiser.sofMapLookup[filterName]
    _, rawGroups = organiser.get_raw_frames_and_groups(
        ttype=ttype,
        tech=entry["eso dpr tech"],
        recipe=entry["recipe"],
        recipeOrder=0,
        filterName=filterName,
        unprocessedOnly=True,
    )
    organiser.conn.close()
    return rawGroups


def _members(rawGroups: pd.DataFrame) -> list[str]:
    return sorted(path.split("/")[-1] for paths in rawGroups["filepaths"] for path in paths)


@pytest.mark.parametrize("tech", [FIXED_SKY_TECH, GENERIC_TECH])
def test_sky_frames_join_the_sof_of_the_object_frames_of_their_run(tmp_path, log, tech) -> None:
    # ARRANGE: TWO ON FRAMES AND TWO OFF FRAMES OF ONE TEMPLATE RUN
    frames = _offset_run("VIS", tech, ["OBJECT", "OBJECT", "SKY", "SKY"])

    # ACT
    rawGroups = _offset_groups(tmp_path, log, frames)

    # ASSERT
    assert rawGroups["counts"].tolist() == [4]
    assert rawGroups["eso dpr type"].tolist() == ["OBJECT"]
    assert _members(rawGroups) == ["off-VIS-0.fits", "off-VIS-1.fits", "off-VIS-2.fits", "off-VIS-3.fits"]


def test_sky_frames_of_a_run_keep_to_their_own_arm(tmp_path, log) -> None:
    # ARRANGE
    frames = pd.concat(
        [
            _offset_run("VIS", GENERIC_TECH, ["OBJECT", "SKY"]),
            _offset_run("NIR", GENERIC_TECH, ["OBJECT", "SKY"], firstSecond=10),
        ],
        ignore_index=True,
    )

    # ACT
    rawGroups = _offset_groups(tmp_path, log, frames)

    # ASSERT
    byArm = {row["eso seq arm"]: sorted(p.split("/")[-1] for p in row["filepaths"]) for _, row in rawGroups.iterrows()}
    assert byArm == {"NIR": ["off-NIR-10.fits", "off-NIR-11.fits"], "VIS": ["off-VIS-0.fits", "off-VIS-1.fits"]}


def test_sky_frame_of_a_run_with_no_object_frame_is_left_out(tmp_path, log) -> None:
    # ARRANGE: THE SKY FRAME'S RUN HAS NO OBJECT FRAME ON THIS ARM, SO NO SOF CAN USE IT
    frames = pd.concat(
        [
            _offset_run("VIS", GENERIC_TECH, ["OBJECT", "SKY"]),
            _offset_run("VIS", GENERIC_TECH, ["SKY"], tplStart="2025-12-01T03:43:49", prefix="lone", firstSecond=20),
        ],
        ignore_index=True,
    )

    # ACT
    rawGroups = _offset_groups(tmp_path, log, frames)

    # ASSERT
    assert _members(rawGroups) == ["off-VIS-0.fits", "off-VIS-1.fits"]


def test_sky_frames_of_a_flux_standard_run_join_the_standard_sof(tmp_path, log) -> None:
    # ARRANGE
    frames = _offset_run("VIS", FIXED_SKY_TECH, ["STD,FLUX", "STD,FLUX", "SKY", "SKY"])

    # ACT
    rawGroups = _offset_groups(tmp_path, log, frames, ttype="STD,FLUX", filterName="flux_standard_offset")

    # ASSERT
    assert rawGroups["counts"].tolist() == [4]
    assert rawGroups["eso dpr type"].tolist() == ["STD,FLUX"]


def test_sky_frames_of_every_run_of_an_ob_join_its_one_sof(tmp_path, log) -> None:
    # ARRANGE: AN SOF HOLDS ALL THE TEMPLATE RUNS OF ONE OB, AS IT DID BEFORE THE SKY FRAMES WERE KEPT
    frames = pd.concat(
        [
            _offset_run("VIS", GENERIC_TECH, ["OBJECT", "SKY"]),
            _offset_run("VIS", GENERIC_TECH, ["OBJECT", "SKY"], tplStart="2025-12-10T02:18:56", firstSecond=20),
        ],
        ignore_index=True,
    )

    # ACT
    rawGroups = _offset_groups(tmp_path, log, frames)

    # ASSERT
    assert rawGroups["counts"].tolist() == [4]


@pytest.mark.parametrize(
    ("tech", "filterName", "recipe"),
    [(STARE_TECH, "science_stare", "stare_obj"), (NOD_TECH, "science_nodding", "nod_obj")],
)
def test_stare_and_nod_sofs_take_no_sky_frames(tmp_path, log, tech, filterName, recipe) -> None:
    # ARRANGE: SKY FRAMES OF THE STARE AND NOD TEMPLATES SHARE THE RUN OF THEIR OBJECT FRAMES
    frames = _offset_run("VIS", tech, ["OBJECT", "OBJECT", "SKY"])

    # ACT
    rawGroups = _offset_groups(tmp_path, log, frames, filterName=filterName)

    # ASSERT
    assert rawGroups["counts"].tolist() == [2]
    assert _members(rawGroups) == ["off-VIS-0.fits", "off-VIS-1.fits"]


def test_the_stare_and_nod_entries_of_the_sof_map_keep_their_techniques(tmp_path, log) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)

    # ACT
    organiser._select_instrument(inst=INSTRUMENT)
    techs = {name: filters["eso dpr tech"] for name, filters in organiser.sofMapLookup.items()}

    # ASSERT
    assert techs["science_stare"] == [STARE_TECH]
    assert techs["science_nodding"] == [NOD_TECH]
    assert techs["flux_standard_stare"] == [STARE_TECH]
    assert techs["flux_standard_nodding"] == [NOD_TECH]


@pytest.mark.parametrize("name", ["flux_standard_offset", "telluric_standard_offset", "science_offset"])
def test_every_offset_entry_of_the_sof_map_lists_both_offset_techniques(tmp_path, log, name) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)

    # ACT
    organiser._select_instrument(inst=INSTRUMENT)

    # ASSERT
    assert organiser.sofMapLookup[name]["eso dpr tech"] == [FIXED_SKY_TECH, GENERIC_TECH]


def _workspace(tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch, frames: pd.DataFrame) -> data_organiser:
    """Return the organiser of a workspace of raw frames, before its first `prep`."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    insert_frames(connection, frames)
    organiser.session_create("science")
    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
    monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
    monkeypatch.setattr(organiser, "_sync_raw_frames", lambda: organiser._select_instrument(inst=INSTRUMENT))
    monkeypatch.setattr(organiser, "get_incomplete_sets_report", lambda: (None, None))
    return organiser


def _prepared(tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch, frames: pd.DataFrame) -> data_organiser:
    """Prepare a workspace of raw frames, as `soxspipe prep` does, and return its organiser."""
    organiser = _workspace(tmp_path, log, monkeypatch, frames)
    organiser.prepare(report=False)
    return organiser


def _set_sizes(organiser: data_organiser) -> dict[str, list[int]]:
    """Return the number of frames of each SOF, by recipe."""
    rows = organiser.conn.execute(
        "SELECT recipe, counts FROM raw_frame_sets WHERE `eso seq arm` = 'VIS' ORDER BY recipe, counts"
    ).fetchall()
    return {recipe: [count for name, count in rows if name == recipe] for recipe, _ in rows}


def test_prepare_puts_every_offset_frame_in_the_offset_sof_and_leaves_stare_and_nod_alone(
    tmp_path, log, monkeypatch
) -> None:
    # ARRANGE: A FIXEDSKYOFFSET RUN AND A GENERICOFFSET RUN OF ONE OB EACH, AND A STARE AND A NOD RUN WITH A SKY FRAME
    frames = pd.concat(
        [
            _offset_run("VIS", FIXED_SKY_TECH, ["OBJECT", "OBJECT", "SKY", "SKY"], prefix="fixed"),
            _offset_run(
                "VIS", GENERIC_TECH, ["OBJECT", "SKY"], tplStart="2025-12-10T02:18:56", prefix="generic", exptime=60.0
            ).assign(**{"eso obs id": "102", "eso obs name": "OB 102"}),
            _offset_run("VIS", STARE_TECH, ["OBJECT", "SKY"], tplStart="2025-12-11T02:18:56", prefix="stare").assign(
                **{"eso obs id": "103", "eso obs name": "OB 103"}
            ),
            _offset_run(
                "VIS", NOD_TECH, ["OBJECT", "OBJECT", "SKY"], tplStart="2025-12-12T02:18:56", prefix="nod"
            ).assign(**{"eso obs id": "104", "eso obs name": "OB 104"}),
        ],
        ignore_index=True,
    )

    # ACT
    organiser = _prepared(tmp_path, log, monkeypatch, frames)

    # ASSERT: EVERY OFFSET FRAME IS IN AN OFFSET SOF; THE SKY FRAMES OF THE OTHER TEMPLATES ARE STILL LEFT OUT
    sizes = _set_sizes(organiser)
    organiser.conn.close()
    assert sizes["offset_obj"] == [2, 4]
    assert sizes["stare_obj"] == [1]
    assert sizes["nod_obj"] == [2]


# THE UPSTREAM PRODUCTS A VIS OFFSET SOF NEEDS: RECIPE, PRODUCT CATEGORY AND SLIT (RESPONSE: SLIT5.0 STANDARD)
UPSTREAM = (
    ("mbias", "MASTER_BIAS_VIS", "SLIT1.0"),
    ("mflat", "MASTER_FLAT_VIS", "SLIT1.0"),
    ("mflat", "ORDER_TAB_VIS", "SLIT1.0"),
    ("spat_solution", "DISP_TAB_VIS", "SLIT1.0"),
    ("spat_solution", "DISP_IMAGE_VIS", "SLIT1.0"),
    ("nod_std", "RESP_TAB_VIS", "SLIT5.0"),
)


def _seed_upstream(organiser: data_organiser) -> None:
    """Add complete, passing upstream products of the VIS arm, as the recipes would have left them."""
    order = {filters["recipe"]: index for index, filters in enumerate(organiser.sofMapLookup.values())}
    columns = [row[1] for row in organiser.conn.execute("PRAGMA table_info(product_frames)")]
    values = organiser.conn.execute("SELECT * FROM product_frames WHERE recipe = 'offset_obj'").fetchone()
    template = dict(zip(columns, values, strict=True))
    rows = []
    for recipe, catg, slit in UPSTREAM:
        name = f"{recipe}-{catg}.fits"
        rows.append(
            {
                **template,
                "file": name,
                "sof": f"{recipe}-upstream.sof",
                "filepath": f"./reduced/2025-11-29/{name}",
                "eso pro catg": catg,
                "recipe": recipe,
                "recipe_order": order[recipe],
                "slit": slit,
                "complete": 1,
                "status": None,
                "mjd-obs": 60980.0,
            }
        )
    frame = pd.DataFrame(rows).drop(columns=["uuid"])
    frame.to_sql("product_frames", con=organiser.conn, index=False, if_exists="append")
    # A PRODUCT THAT IS NOT IN THE SESSION'S SOF MAP IS REQUEUED BY THE NEXT BUILD
    sofMap = frame[["sof", "file", "filepath", "eso pro catg"]].rename(columns={"eso pro catg": "tag"})
    sofMap.assign(complete=1).to_sql(
        f"sof_map_{organiser.sessionId}", con=organiser.conn, index=False, if_exists="append"
    )
    organiser.conn.commit()


def _offset_sofs(organiser: data_organiser) -> dict[str, set[str]]:
    """Return the raw frames of each complete VIS `offset_obj` SOF, by SOF name."""
    rows = organiser.conn.execute(
        "SELECT DISTINCT m.sof, m.file FROM sof_map_science m JOIN product_frames p ON p.sof = m.sof "
        "WHERE p.recipe = 'offset_obj' AND p.`eso seq arm` = 'VIS' AND m.filepath LIKE './raw/%'"
    ).fetchall()
    sofs: dict[str, set[str]] = {}
    for sof, name in rows:
        sofs.setdefault(sof, set()).add(name)
    return sofs


@pytest.mark.parametrize("tech", [FIXED_SKY_TECH, GENERIC_TECH])
def test_late_sky_frames_rejoin_the_offset_sof_of_their_run(tmp_path, log, monkeypatch, tech) -> None:
    # ARRANGE: THE ON FRAMES OF A RUN ARE GROUPED AND THEIR SOF IS COMPLETE; THE OFF FRAMES COME LATER
    onFrames = _offset_run("VIS", tech, ["OBJECT", "OBJECT"], exptime=60.0)
    offFrames = _offset_run("VIS", tech, ["SKY", "SKY"], exptime=60.0, firstSecond=2)
    organiser = _workspace(tmp_path, log, monkeypatch, onFrames)
    organiser.prepare(report=False)
    _seed_upstream(organiser)
    organiser.prepare(report=False)
    (earlySof,) = _offset_sofs(organiser)

    # ACT
    insert_frames(organiser.conn, offFrames)
    organiser.conn.commit()
    organiser.prepare(report=False)

    # ASSERT: ONE OFFSET SOF, WITH THE NAME IT HAD, HOLDS ALL 4 FRAMES; NO SKY-ONLY SOF EXISTS
    sofs = _offset_sofs(organiser)
    allSofs = {
        row[0] for row in organiser.conn.execute("SELECT DISTINCT sof FROM raw_frame_sets WHERE recipe = 'offset_obj'")
    }
    organiser.conn.close()
    assert sofs == {earlySof: {"off-VIS-0.fits", "off-VIS-1.fits", "off-VIS-2.fits", "off-VIS-3.fits"}}
    assert allSofs == {earlySof}
