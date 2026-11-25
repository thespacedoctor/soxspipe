"""Integration contracts for the NIR warm-detector and slit-wheel ignore rules (DY-2394)."""

from __future__ import annotations

import pytest

from tests.factories import workspace_organiser

pytestmark = pytest.mark.integration

MAX_DETECTOR_TEMP_SETTING = "nir-max-detector-temp-k"


def _insert_frame(
    connection, name, *, arm="NIR", dprType="DARK", slit="SLIT1.0", simulation=0, nirTemp=48.0, tech="IMAGE"
) -> None:
    """Seed one raw-frame row with the columns the ignore rules read."""
    connection.execute(
        """INSERT INTO raw_frames (
            instrume, file, "mjd-obs", "mjd-date", "date-obs", "night start mjd", "night start date",
            "eso dpr catg", "eso dpr type", "eso dpr tech", exptime, object, filepath,
            "eso seq arm", slit, simulation, "nir temp k"
        ) VALUES ('SOXS', ?, 60311.1, '2024-01-02', '2024-01-02T03:04:05', 60310, '2024-01-01',
            'CALIB', ?, ?, 10, 'DARK', ?, ?, ?, ?, ?)""",
        (name, dprType, tech, f"./raw/{name}", arm, slit, simulation, nirTemp),
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


@pytest.mark.parametrize("dprType", ["DARK", "BIAS", "LAMP,FLAT", "OBJECT"])
def test_nir_frame_with_a_warm_detector_is_ignored(organiser, dprType) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "warm.fits", dprType=dprType, nirTemp=181.0)
    _insert_frame(organiser.conn, "cold.fits", dprType=dprType, nirTemp=48.0)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"warm.fits"}


def test_nir_frame_at_the_default_limit_is_kept(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "limit.fits", nirTemp=60.0)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


def test_nir_frame_without_a_detector_temperature_is_kept(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "unknown.fits", nirTemp=None)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


def test_warm_vis_frame_is_not_ignored_by_the_nir_temperature_rule(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "vis.fits", arm="VIS", slit="SLIT1.0", nirTemp=250.0)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


def test_the_temperature_limit_comes_from_the_settings(organiser) -> None:
    # ARRANGE
    organiser.settings = {MAX_DETECTOR_TEMP_SETTING: 40}
    _insert_frame(organiser.conn, "warmish.fits", nirTemp=48.0)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"warmish.fits"}


def test_the_warm_detector_rejection_logs_the_reason_once_per_night(organiser, log) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "warm-1.fits", nirTemp=181.0)
    _insert_frame(organiser.conn, "warm-2.fits", nirTemp=190.5)
    _insert_frame(organiser.conn, "cold.fits", nirTemp=48.0)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    rejections = [m for level, m in log.messages if level == "info" and "detector temperature" in m]
    assert len(rejections) == 1
    assert "2 NIR frames" in rejections[0]
    assert "2024-01-01" in rejections[0]
    assert "181.0 to 190.5 K" in rejections[0]
    assert "limit of 60 K" in rejections[0]


def test_a_second_pass_does_not_log_the_warm_frames_again(organiser, log) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "warm.fits", nirTemp=181.0)
    organiser._flag_files_to_ignore()
    log.messages.clear()

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"warm.fits"}
    assert [m for level, m in log.messages if "detector temperature" in m] == []


@pytest.mark.parametrize("dprType", ["DARK", "BIAS"])
def test_nir_bias_or_dark_with_an_undefined_slit_is_kept(organiser, dprType) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "dark.fits", dprType=dprType, slit="UNDEFINED")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == set()


@pytest.mark.parametrize("dprType", ["LAMP,FLAT", "LAMP,WAVE", "OBJECT"])
def test_nir_frame_that_needs_the_slit_is_ignored_when_the_slit_is_undefined(organiser, dprType) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "undefined.fits", dprType=dprType, slit="UNDEFINED")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"undefined.fits"}


@pytest.mark.parametrize("arm", ["VIS", "UVB"])
def test_vis_and_uvb_darks_with_an_undefined_slit_are_still_ignored(organiser, arm) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "dark.fits", arm=arm, slit="UNDEFINED")

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"dark.fits"}


def test_simulated_nir_dark_is_still_ignored_and_a_simulated_science_frame_is_kept(organiser) -> None:
    # ARRANGE
    _insert_frame(organiser.conn, "sim-dark.fits", simulation=1)
    _insert_frame(organiser.conn, "sim-object.fits", dprType="OBJECT", slit="SLIT1.0", simulation=1)

    # ACT
    organiser._flag_files_to_ignore()

    # ASSERT
    assert _ignored(organiser) == {"sim-dark.fits"}


def test_a_non_numeric_temperature_limit_stops_with_a_clear_error(organiser) -> None:
    # ARRANGE
    organiser.settings = {MAX_DETECTOR_TEMP_SETTING: "cold"}

    # ACT, ASSERT
    with pytest.raises(ValueError, match="nir-max-detector-temp-k.*number of kelvin"):
        organiser._flag_files_to_ignore()


def test_the_packaged_soxs_settings_set_the_limit_to_60_k() -> None:
    # ARRANGE
    from importlib.resources import files

    import yaml

    # ACT
    settings = yaml.safe_load(files("soxspipe").joinpath("soxs_default_settings.yaml").read_text(encoding="utf-8"))

    # ASSERT
    assert settings[MAX_DETECTOR_TEMP_SETTING] == 60
