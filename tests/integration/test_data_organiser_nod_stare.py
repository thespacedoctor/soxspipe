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


def _nir_frames(
    kind, count, *, dprType, tech, tplName, exptime, firstSecond=0, obsId="101", obsName=None, tplStart=None
):
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
                "eso obs id": obsId,
                "eso obs name": obsName or f"OB {obsId}",
                "exptime": exptime,
                "slit": "SLIT5.0" if kind == "nod" else "--",
                "object": "SN2024abc" if kind == "nod" else "DARK",
                "template": tplName,
                "eso tpl name": tplName,
                "eso tpl nexp": count,
                "eso tpl expno": index + 1,
                "eso tpl start": tplStart or ("2024-01-02T03:00:00" if kind == "nod" else "2024-01-02T02:00:00"),
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
    assert products == [("SCI_SLIT_FLUX_NIR", NOD_TECH, sof.replace(".sof", "_SPECTRUM.fits"))]


def _set_reduce_as_stare(sessionsDir, sessionId, block="soxs-nod"):
    """Turn `reduce-as-stare` on in one block (`soxs-nod` or `soxs-nod-std`) of one session's settings file."""
    settingsPath = sessionsDir / sessionId / "soxspipe.yaml"
    text = settingsPath.read_text()
    blockStart = text.index(f"\n{block}:\n")
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
    stemFile = _nod_sof(organiser).replace(".sof", "_STAREMODE_SPECTRUM.fits")
    stareRows = [row for row in _nod_obj_rows(fresh) if row[0] == stemFile]
    assert stareRows
    assert {row[1].split("/")[3] for row in stareRows} == {"soxs-stare"}
    assert _unloaded_settings_warnings(log) == []


def _nod_product_rows(organiser):
    """Return the NIR product rows of the nod object and nod standard SOFs, by recipe."""
    return {
        recipe: organiser.conn.execute(
            'SELECT "eso pro catg", file, filepath FROM product_frames WHERE recipe = ? AND "eso seq arm" = \'NIR\' '
            "AND sof NOT LIKE '%-upstream.sof' ORDER BY \"eso pro catg\"",
            (recipe,),
        ).fetchall()
        for recipe in ("nod_obj", "nod_std")
    }


def _object_and_standard_rows(tmp_path, log, monkeypatch, organisers, stareBlocks):
    """Return the NIR nod rows before and after a regroup with `reduce-as-stare` on in the ``stareBlocks``."""
    frames = [
        _nir_frames("dark", 3, dprType="DARK", tech="IMAGE", tplName="SOXS_cal_dark", exptime=300.0),
        _nir_frames("nod", 4, dprType="OBJECT", tech=NOD_TECH, tplName="SOXS_slt_obs_Nod", exptime=300.0),
        _nir_frames(
            "nod",
            4,
            dprType="STD,FLUX",
            tech=NOD_TECH,
            tplName="SOXS_slt_cal_Nod",
            exptime=60.0,
            firstSecond=20,
            obsId="102",
            # GD71 IS IN THE STATIC STANDARDS LIBRARY, SO `prepare` DOES NOT FLAG THE FRAMES TO IGNORE
            obsName="GD71",
            tplStart="2024-01-02T04:00:00",
        ),
    ]
    organiser = _workspace(tmp_path, log, monkeypatch, pd.concat(frames, ignore_index=True))
    organisers.append(organiser)
    organiser.prepare(report=False)
    _seed_upstream(organiser, UPSTREAM)
    organiser.prepare(report=False)
    before = _nod_product_rows(organiser)

    sessionsDir = Path(organiser.sessionsDir)
    organiser.session_create("beta")
    organiser.session_switch("science")
    for block in stareBlocks:
        _set_reduce_as_stare(sessionsDir, "beta", block)
    fresh = _fresh_organiser(tmp_path, log)
    organisers.append(fresh)
    fresh.session_switch("beta")
    return before, _nod_product_rows(fresh)


def _as_stare_rows(rows, folder):
    """Return the nod-mode ``rows`` as the stare recipe writes them: in ``folder``, with `_STAREMODE` in the name."""
    stareRows = []
    for catg, file, filepath in rows:
        stareFile = file.replace("_SPECTRUM.fits", "_STAREMODE_SPECTRUM.fits").replace(
            "_RESP.fits", "_STAREMODE_RESP.fits"
        )
        stareRows.append((catg, stareFile, f"{filepath.rsplit('/', 2)[0]}/{folder}/{stareFile}"))
    return stareRows


@pytest.mark.parametrize(
    ("stareBlocks", "isObjectStare", "isStandardStare"),
    [
        ((), False, False),
        (("soxs-nod-std",), False, True),
        (("soxs-nod",), True, False),
        (("soxs-nod", "soxs-nod-std"), True, True),
    ],
    ids=["neither", "standard-only", "object-only", "both"],
)
def test_each_nod_sof_records_the_stare_rows_only_when_its_own_reduce_as_stare_key_is_on(
    tmp_path, log, monkeypatch, organisers, stareBlocks, isObjectStare, isStandardStare
) -> None:
    # ACT
    before, after = _object_and_standard_rows(tmp_path, log, monkeypatch, organisers, stareBlocks)

    # ASSERT: A REGROUP KEEPS THE OLD NOD-MODE ROWS, SO A SOF REDUCED AS STARE HAS BOTH, AND ANY OTHER HAS ONLY THOSE
    assert {row[0] for row in before["nod_std"]} == {"SCI_SLIT_FLUX_NIR", "RESP_TAB_NIR"}
    assert {row[2].split("/")[3] for row in before["nod_obj"]} == {"soxs-nod"}
    assert {row[2].split("/")[3] for row in before["nod_std"]} == {"soxs-nod-std"}
    objectStare = _as_stare_rows(before["nod_obj"], "soxs-stare") if isObjectStare else []
    standardStare = _as_stare_rows(before["nod_std"], "soxs-stare-std") if isStandardStare else []
    assert sorted(after["nod_obj"]) == sorted(before["nod_obj"] + objectStare)
    assert sorted(after["nod_std"]) == sorted(before["nod_std"] + standardStare)
