"""The data organiser adds a matching NIR master dark to the NIR nod SOFs, and needs none to complete them (DY-2194)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from tests.factories import raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration

INSTRUMENT = "SOXS"
NOD_TECH = "ECHELLE,SLIT,NODDING"
# THE UPSTREAM RECIPES A NIR NOD SOF NEEDS, IN THE ORDER THE SOF MAP LISTS THEM
UPSTREAM = (
    ("mflat", "MASTER_FLAT_NIR"),
    ("mflat", "ORDER_TAB_NIR"),
    ("spat_solution", "DISP_TAB_NIR"),
    ("spat_solution", "DISP_IMAGE_NIR"),
    ("nod_std", "RESP_TAB_NIR"),
)


def _insert_frames(connection, frames):
    """Insert dataframe rows into `raw_frames`, keeping only the columns the table has."""
    columns = {row[1] for row in connection.execute("PRAGMA table_info(raw_frames);")}
    frames.loc[:, frames.columns.isin(columns)].to_sql("raw_frames", con=connection, index=False, if_exists="append")


def _nir_frames(kind, count, *, dprType, tech, tplName, exptime, firstSecond=0):
    """Return ``count`` NIR exposures of one template run."""
    template = raw_frame_table().iloc[0].to_dict()
    rows = []
    for index in range(count):
        second = firstSecond + index
        name = f"{kind}-NIR-{second}.fits"
        rows.append(
            {
                **template,
                "eso seq arm": "NIR",
                "eso dpr catg": "SCIENCE" if kind == "nod" else "CALIB",
                "eso dpr tech": tech,
                "eso dpr type": dprType,
                "eso obs id": "101",
                "eso obs name": "OB 101",
                "exptime": exptime,
                "slit": "SLIT5.0" if kind == "nod" else "--",
                "object": "SN2024abc" if kind == "nod" else "DARK",
                "template": tplName,
                "eso tpl name": tplName,
                "eso tpl nexp": count,
                "eso tpl expno": index + 1,
                "eso tpl start": "2024-01-02T03:00:00" if kind == "nod" else "2024-01-02T02:00:00",
                "file": name,
                "filepath": f"./raw/2024-01-01/{name}",
                "date-obs": f"2024-01-02T03:05:{second:02d}.000"
                if kind == "nod"
                else f"2024-01-02T02:05:{second:02d}.000",
                "mjd-obs": (60311.2118 if kind == "nod" else 60311.1701) + second * 0.0001,
                "set_first_file": name,
            }
        )
    return pd.DataFrame(rows)


def _workspace(tmp_path, log, monkeypatch, frames):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    _insert_frames(connection, frames)
    organiser.session_create("science")

    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
    monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
    monkeypatch.setattr(organiser, "_sync_raw_frames", lambda: organiser._select_instrument(inst=INSTRUMENT))
    monkeypatch.setattr(organiser, "get_incomplete_sets_report", lambda: (None, None))
    return organiser


def _nod_sof(organiser, recipe="nod_obj"):
    rows = organiser.conn.execute(
        "SELECT DISTINCT sof FROM product_frames WHERE recipe = ? AND \"eso seq arm\" = 'NIR' "
        "AND sof NOT LIKE '%-upstream.sof'",
        (recipe,),
    ).fetchall()
    assert len(rows) == 1
    return rows[0][0]


def _tags(organiser, sof):
    return sorted(row[0] for row in organiser.conn.execute("SELECT DISTINCT tag FROM sof_map WHERE sof = ?", (sof,)))


def _complete(organiser, sof):
    return {row[0] for row in organiser.conn.execute("SELECT complete FROM product_frames WHERE sof = ?", (sof,))}


def _seed_upstream(organiser, upstream):
    """Add complete, passing upstream products of the NIR arm, as the recipes would have left them."""
    order = {filters["recipe"]: index for index, filters in enumerate(data_organiser_sof_map(organiser).values())}
    template = dict(
        zip(
            [row[1] for row in organiser.conn.execute("PRAGMA table_info(product_frames)")],
            organiser.conn.execute(
                "SELECT * FROM product_frames WHERE recipe = 'nod_obj' OR recipe = 'nod_std'"
            ).fetchone(),
            strict=True,
        )
    )
    rows = []
    for recipe, catg in upstream:
        name = f"{recipe}-{catg}.fits"
        rows.append(
            {
                **template,
                "uuid": None,
                "file": name,
                "sof": f"{recipe}-upstream.sof",
                "filepath": f"./reduced/2024-01-01/{name}",
                "eso pro catg": catg,
                "recipe": recipe,
                "recipe_order": order[recipe],
                "slit": "SLIT5.0",
                "complete": 1,
                "status": None,
                "mjd-obs": 60311.1,
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


def data_organiser_sof_map(organiser):
    return organiser.sofMapLookup


def _prepared(tmp_path, log, monkeypatch, *, withDark, nodType="OBJECT"):
    """Prepare a workspace of NIR nod frames, with or without the raw darks of a master dark, and its upstream."""
    frames = [
        _nir_frames("nod", 4, dprType=nodType, tech=NOD_TECH, tplName="SOXS_slt_obs_Nod", exptime=300.0),
    ]
    if withDark:
        frames.insert(0, _nir_frames("dark", 3, dprType="DARK", tech="IMAGE", tplName="SOXS_cal_dark", exptime=300.0))
    organiser = _workspace(tmp_path, log, monkeypatch, pd.concat(frames, ignore_index=True))
    organiser.prepare(report=False)
    _seed_upstream(organiser, UPSTREAM)
    organiser.prepare(report=False)
    return organiser


@pytest.fixture
def organisers():
    built = []
    yield built
    for organiser in built:
        organiser.conn.close()


def test_a_nir_nod_sof_lists_the_matching_master_dark(tmp_path, log, monkeypatch, organisers) -> None:
    # ACT
    organiser = _prepared(tmp_path, log, monkeypatch, withDark=True)
    organisers.append(organiser)
    sof = _nod_sof(organiser)

    # ASSERT
    assert _complete(organiser, sof) == {1}
    assert _tags(organiser, sof) == [
        "DISP_IMAGE_NIR",
        "DISP_TAB_NIR",
        "MASTER_DARK_NIR",
        "MASTER_FLAT_NIR",
        "OBJECT_NIR",
        "ORDER_TAB_NIR",
        "RESP_TAB_NIR",
    ]


def test_a_nir_nod_standard_sof_lists_the_matching_master_dark(tmp_path, log, monkeypatch, organisers) -> None:
    # ACT
    organiser = _prepared(tmp_path, log, monkeypatch, withDark=True, nodType="STD,FLUX")
    organisers.append(organiser)
    sof = _nod_sof(organiser, "nod_std")

    # ASSERT
    assert _complete(organiser, sof) == {1}
    assert "MASTER_DARK_NIR" in _tags(organiser, sof)


def test_a_nir_nod_sof_with_no_master_dark_is_still_complete_without_one(
    tmp_path, log, monkeypatch, organisers
) -> None:
    # ACT
    organiser = _prepared(tmp_path, log, monkeypatch, withDark=False)
    organisers.append(organiser)
    sof = _nod_sof(organiser)

    # ASSERT
    assert _complete(organiser, sof) == {1}
    assert "MASTER_DARK_NIR" not in _tags(organiser, sof)
    assert {"MASTER_FLAT_NIR", "ORDER_TAB_NIR", "DISP_TAB_NIR", "DISP_IMAGE_NIR", "RESP_TAB_NIR"} <= set(
        _tags(organiser, sof)
    )


def test_the_predicted_product_of_a_nir_nod_sof_is_unchanged(tmp_path, log, monkeypatch, organisers) -> None:
    # ACT
    organiser = _prepared(tmp_path, log, monkeypatch, withDark=True)
    organisers.append(organiser)
    sof = _nod_sof(organiser)
    products = organiser.conn.execute(
        'SELECT "eso pro catg", "eso pro tech", file FROM product_frames WHERE sof = ?', (sof,)
    ).fetchall()

    # ASSERT: ONE NOD-MODE PRODUCT, NAMED AFTER THE SOF, WITH NO MARKER AND NO STARE TECHNIQUE
    assert products == [("OBJECT_TAB_NIR", NOD_TECH, sof.replace(".sof", "_EXTRACTED_MERGED.fits"))]


def _set_reduce_as_stare(sessionsDir, sessionId):
    """Turn `reduce-as-stare` on in the `soxs-nod` block of one session's settings file."""
    settingsPath = sessionsDir / sessionId / "soxspipe.yaml"
    text = settingsPath.read_text()
    blockStart = text.index("\nsoxs-nod:\n")
    flagStart = text.index("reduce-as-stare: False", blockStart)
    settingsPath.write_text(
        text[:flagStart] + "reduce-as-stare: True" + text[flagStart + len("reduce-as-stare: False") :]
    )


def _fresh_organiser(tmp_path, log):
    """Return a new organiser on the workspace of ``tmp_path``, as the command line builds one for a session command."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn, _ = organiser._get_or_create_db_connection()
    return organiser


def _nod_obj_rows(organiser):
    return organiser.conn.execute(
        "SELECT file, filepath FROM product_frames WHERE recipe = 'nod_obj' AND \"eso seq arm\" = 'NIR'"
    ).fetchall()


def _unloaded_settings_warnings(log):
    return [message for level, message in log.messages if level == "warning" and "settings are not loaded" in message]


@pytest.mark.parametrize("command", ["session_create", "session_switch"])
def test_a_session_command_on_a_fresh_organiser_records_the_stare_mode_row_when_reduce_as_stare_is_on(
    tmp_path, log, monkeypatch, organisers, command
) -> None:
    # ARRANGE
    organiser = _prepared(tmp_path, log, monkeypatch, withDark=True)
    organisers.append(organiser)
    sessionsDir = Path(organiser.sessionsDir)
    organiser.session_create("beta")
    organiser.session_switch("science")
    _set_reduce_as_stare(sessionsDir, "beta")
    fresh = _fresh_organiser(tmp_path, log)
    organisers.append(fresh)
    log.messages.clear()

    # ACT
    getattr(fresh, command)("beta")

    # ASSERT
    stemFile = _nod_sof(organiser).replace(".sof", "_STAREMODE_EXTRACTED_MERGED.fits")
    stareRows = [row for row in _nod_obj_rows(fresh) if row[0] == stemFile]
    assert stareRows
    assert {row[1].split("/")[3] for row in stareRows} == {"soxs-stare"}
    assert _unloaded_settings_warnings(log) == []
