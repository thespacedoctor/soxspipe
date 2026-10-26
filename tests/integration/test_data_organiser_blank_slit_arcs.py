"""Integration contracts for rejecting pinhole and multi-pinhole arcs taken with a blank slit (DY-1577)."""

from __future__ import annotations

import pytest

from tests.factories import workspace_organiser

pytestmark = pytest.mark.integration

MULTI_PINHOLE_TECH = "ECHELLE,MULTI-PINHOLE"


def _insert_frame(connection, name, *, arm, slit, tech=MULTI_PINHOLE_TECH, dprType="LAMP,WAVE") -> None:
    """Seed one raw-frame row with the columns the ignore rules read."""
    connection.execute(
        """INSERT INTO raw_frames (
            instrume, file, "mjd-obs", "mjd-date", "date-obs", "night start mjd",
            "eso dpr catg", "eso dpr type", "eso dpr tech", exptime, object, filepath,
            "eso seq arm", slit
        ) VALUES ('SOXS', ?, 60311.1, '2024-01-02', '2024-01-02T03:04:05', 60310,
            'CALIB', ?, ?, 10, 'LAMP', ?, ?, ?)""",
        (name, dprType, tech, f"./raw/{name}", arm, slit),
    )


@pytest.fixture
def organiser(tmp_path, log):
    """A SOXS organiser on an empty workspace database, with default settings."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn, _ = organiser._get_or_create_db_connection()
    organiser.instrument = "SOXS"
    organiser.settings = {}
    yield organiser
    organiser.conn.close()


def _ignored(organiser):
    """Return the names of the frames flagged to ignore."""
    rows = organiser.conn.execute("SELECT file FROM raw_frames WHERE ignore = 1").fetchall()
    return {row[0] for row in rows}


@pytest.mark.parametrize("arm", ["NIR", "VIS"])
def test_blank_slit_multi_pinhole_arc_is_ignored_on_each_arm(organiser, arm) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "blank-mph.fits", arm=arm, slit="BLANK")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"blank-mph.fits"}


def test_blank_slit_multi_pinhole_arc_is_left_out_of_the_valid_raw_frames(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "blank-mph.fits", arm="NIR", slit="BLANK")
    _insert_frame(organiser.conn, "good-mph.fits", arm="NIR", slit="MULTPIN")
    organiser.conn.execute('UPDATE raw_frames SET "eso tpl nexp" = 1, set_size = 1')

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    validFiles = {row[0] for row in organiser.conn.execute("SELECT file FROM raw_frames_valid")}
    assert validFiles == {"good-mph.fits"}


def test_multpin_multi_pinhole_arc_is_not_ignored(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "good-nir.fits", arm="NIR", slit="MULTPIN")
    _insert_frame(organiser.conn, "good-vis.fits", arm="VIS", slit="MULTPIN")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


def test_blank_slit_frames_of_other_techniques_are_not_rejected_by_the_arc_rule(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "blank-slit-arc.fits", arm="NIR", slit="BLANK", tech="ECHELLE,SLIT")
    _insert_frame(
        organiser.conn,
        "blank-stare.fits",
        arm="NIR",
        slit="BLANK",
        tech="ECHELLE,SLIT,STARE",
        dprType="OBJECT",
    )

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


@pytest.mark.parametrize("arm", ["NIR", "VIS"])
def test_blank_slit_single_pinhole_arc_is_ignored_on_each_arm(organiser, arm) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "blank-ph.fits", arm=arm, slit="BLANK", tech="ECHELLE,PINHOLE")
    _insert_frame(organiser.conn, "good-ph.fits", arm=arm, slit="PINHOLE", tech="ECHELLE,PINHOLE")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"blank-ph.fits"}


def test_rejection_logs_the_file_name_and_the_reason(organiser, log) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "blank-mph.fits", arm="NIR", slit="BLANK")
    _insert_frame(organiser.conn, "blank-ph.fits", arm="VIS", slit="BLANK", tech="ECHELLE,PINHOLE")
    _insert_frame(organiser.conn, "good-mph.fits", arm="NIR", slit="MULTPIN")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    rejections = [m for level, m in log.messages if level == "info" and "slit wheel reads BLANK" in m]
    assert len(rejections) == 2
    assert "ECHELLE,MULTI-PINHOLE frame `blank-mph.fits`" in rejections[0]
    assert "ECHELLE,PINHOLE frame `blank-ph.fits`" in rejections[1]
    assert not any("good-mph.fits" in m for _, m in log.messages if "slit wheel reads BLANK" in m)


def test_a_second_pass_keeps_the_frame_ignored_and_does_not_warn_again(organiser, log) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "blank-mph.fits", arm="VIS", slit="BLANK")
    organiser._flag_files_to_ignore()
    log.messages.clear()

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"blank-mph.fits"}
    assert [message for level, message in log.messages if level == "warning"] == []


def test_x_shooter_pinhole_frames_are_left_alone(organiser) -> None:
    # ARRANGE
    organiser.instrument = "XSH"
    _insert_frame(organiser.conn, "xsh-mph.fits", arm="VIS", slit="BLANK")
    _insert_frame(organiser.conn, "xsh-ph.fits", arm="VIS", slit="BLANK", tech="ECHELLE,PINHOLE")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()
