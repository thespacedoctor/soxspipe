"""Fixtures shared by the integration tests of the data organiser."""

from __future__ import annotations

import pytest

from tests.factories import calibration_chain_frames, insert_frames, workspace_organiser

INSTRUMENT = "SOXS"


@pytest.fixture
def workspace_with_frames(tmp_path, log, monkeypatch):
    """Build a prepared two-session workspace: the calibration chain, the flat exposures given, and any extra frames.

    Call the returned function as `build(flatFrames, extraFrames=())`. Each frame table in ``extraFrames`` is
    inserted before the sessions are created and `prep` first runs.
    """
    organisers = []

    def build(flatFrames, extraFrames=()):
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        organiser = workspace_organiser(tmp_path, log=log)
        connection, _ = organiser._get_or_create_db_connection()
        organiser.conn = connection
        insert_frames(connection, calibration_chain_frames())
        insert_frames(connection, flatFrames)
        for frames in extraFrames:
            insert_frames(connection, frames)
        organiser.session_create("archive")
        organiser.session_create("science")

        monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
        monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
        monkeypatch.setattr(organiser, "_sync_raw_frames", lambda: organiser._select_instrument(inst=INSTRUMENT))
        monkeypatch.setattr(organiser, "get_incomplete_sets_report", lambda: (None, None))
        organiser.prepare(report=False)
        organisers.append(organiser)
        return organiser

    yield build
    for organiser in organisers:
        if getattr(organiser, "conn", None) is not None:
            organiser.conn.close()
