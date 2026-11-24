"""The on-target frame list and the frame-accounting cards it gives a merged spectrum (DY-2356).

An on-target frame is a frame whose light is in the product: a stare frame, a nod A or B frame, or an offset ON frame.
The list is built from the raw headers before the recipe edits them, and the cards are stamped on a copy of the
product header, never on a shared frame header.
"""

from __future__ import annotations

from typing import Any

import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.on_target_frames import (
    build_on_target_frames,
    frame_accounting_cards,
    frames_named,
    read_on_target_frame,
    with_frame_accounting,
    without_frames,
    without_untraced_positions,
)
from tests.factories.on_target_frames import raw_science_header

pytestmark = pytest.mark.unit

SECONDS_PER_DAY = 86400.0


def _kw(log: Any, instrument: str = "soxs") -> Any:
    return keyword_lookup(log=log, settings={"instrument": instrument}).get


def _names(entries: Any) -> list[str]:
    return [entry.name for entry in entries]


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def _entry(log: Any, name: str, arm: str = "VIS", **headerCards: Any) -> Any:
    header = raw_science_header(arm, arcfile=f"{name}.fits", **headerCards)
    return read_on_target_frame(header, _kw(log), f"{name}_raw.fits", arm, log)


def _cards(entries: Any) -> dict[str, Any]:
    return {keyword: value for keyword, value, _ in frame_accounting_cards(entries)}


# READING ONE FRAME


def test_a_vis_frame_spends_its_exptime_on_target(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS", exptime=1800.0)

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert entry.exposureSeconds == 1800.0
    assert not _warnings(log)


def test_a_nir_frame_spends_its_dit_times_ndit_on_target_not_its_exptime(log: Any) -> None:
    # ARRANGE: EXPTIME IS ONE DIT ON SOXS NIR
    header = raw_science_header("NIR", dit=900.0, ndit=2)

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "NIR", log)

    # ASSERT
    assert header["EXPTIME"] == 900.0
    assert entry.exposureSeconds == 1800.0


def test_an_x_shooter_nir_frame_reads_its_dit_from_the_x_shooter_card(log: Any) -> None:
    # ARRANGE
    header = fits.Header({"HIERARCH ESO DET DIT": 60.0, "HIERARCH ESO DET NDIT": 3, "MJD-OBS": 61000.0})

    # ACT
    entry = read_on_target_frame(header, _kw(log, "xsh"), "frame.fits", "NIR", log)

    # ASSERT
    assert entry.exposureSeconds == 180.0


def test_a_nir_frame_without_an_ndit_card_is_one_integration(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("NIR", dit=300.0, ndit=None)

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "NIR", log)

    # ASSERT
    assert entry.exposureSeconds == 300.0


def test_a_frame_without_a_time_card_has_no_time_on_target_and_a_warning_names_the_card(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS")
    del header["EXPTIME"]

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert entry.exposureSeconds is None
    assert any("frame.fits" in message and "EXPTIME" in message for message in _warnings(log))


def test_a_frame_records_its_identity_cards(log: Any) -> None:
    # ARRANGE
    header = raw_science_header(
        "VIS", mjd=61000.25, arcfile="SOXS.A.fits", obsId=7, progId="115.X.001", tplStart="2026-01-02T03:04:05"
    )

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert (entry.arcfile, entry.mjdObs, entry.obsId, entry.progId, entry.tplStart) == (
        "SOXS.A.fits",
        61000.25,
        7,
        "115.X.001",
        "2026-01-02T03:04:05",
    )


def test_a_frame_names_itself_by_its_arcfile(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS", arcfile="SOXS.A.fits", ORIGFILE="orig.fits")

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert entry.name == "SOXS.A.fits"
    assert entry.provenance == "SOXS.A"
    assert not _warnings(log)


@pytest.mark.parametrize(
    ("cards", "expectedName"),
    [
        ({"ORIGFILE": "orig.fits"}, "orig.fits"),
        ({}, "frame.fits"),
    ],
)
def test_a_frame_with_no_arcfile_is_named_by_its_origfile_then_its_sof_file_with_a_warning(
    log: Any, cards: dict[str, str], expectedName: str
) -> None:
    # ARRANGE
    header = raw_science_header("VIS", arcfile=None, **cards)

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert entry.arcfile is None
    assert entry.name == expectedName
    assert entry.provenance == expectedName.removesuffix(".fits")
    assert any("frame.fits" in message and "ARCFILE" in message for message in _warnings(log))


def test_a_frame_missing_an_identity_card_keeps_none_and_warns(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS", obsId=None, progId=None, tplStart=None)

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert (entry.obsId, entry.progId, entry.tplStart) == (None, None, None)
    message = " ".join(_warnings(log))
    assert "OBS ID" in message and "OBS PROG ID" in message and "TPL START" in message


def test_a_frame_ends_when_its_time_on_target_is_over(log: Any) -> None:
    # ARRANGE
    entry = _entry(log, "a", mjd=61000.0, exptime=864.0)

    # ACT / ASSERT
    assert entry.mjdEnd == pytest.approx(61000.01, abs=1e-9)


def test_a_frame_with_an_unknown_time_has_no_end(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS")
    del header["EXPTIME"]

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert entry.mjdEnd is None


# BUILDING AND SELECTING THE LIST


def test_the_list_is_sorted_by_observation_time(log: Any) -> None:
    # ARRANGE
    headers = [raw_science_header("VIS", mjd=mjd, arcfile=f"{mjd}.fits") for mjd in (61000.3, 61000.1, 61000.2)]

    # ACT
    entries = build_on_target_frames(headers, ["c.fits", "a.fits", "b.fits"], _kw(log), "VIS", log)

    # ASSERT
    assert [entry.mjdObs for entry in entries] == [61000.1, 61000.2, 61000.3]
    assert [entry.sourceName for entry in entries] == ["a.fits", "b.fits", "c.fits"]


def test_the_list_refuses_two_frames_with_one_file_name(log: Any) -> None:
    # ARRANGE
    headers = [raw_science_header("VIS"), raw_science_header("VIS", mjd=61000.2)]

    # ACT / ASSERT
    with pytest.raises(ValueError, match=r"same\.fits"):
        build_on_target_frames(headers, ["same.fits", "same.fits"], _kw(log), "VIS", log)


def test_the_list_refuses_headers_and_names_of_different_lengths(log: Any) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match="2 headers and 1 file names"):
        build_on_target_frames([raw_science_header("VIS")] * 2, ["a.fits"], _kw(log), "VIS", log)


def test_frames_named_selects_the_entries_of_those_files_in_list_order(log: Any) -> None:
    # ARRANGE
    entries = tuple(_entry(log, name, mjd=61000.0 + index) for index, name in enumerate("abc"))

    # ACT
    selected = frames_named(entries, ["c_raw.fits", "a_raw.fits"])

    # ASSERT
    assert [entry.name for entry in selected] == ["a.fits", "c.fits"]


def test_frames_named_refuses_a_file_with_no_entry(log: Any) -> None:
    # ARRANGE
    entries = (_entry(log, "a"),)

    # ACT / ASSERT
    with pytest.raises(ValueError, match=r"zzz\.fits"):
        frames_named(entries, ["a_raw.fits", "zzz.fits"])


def test_leaving_frames_out_gives_a_new_list_and_one_warning_that_names_them_and_the_reason(log: Any) -> None:
    # ARRANGE
    entries = tuple(_entry(log, name, mjd=61000.0 + index) for index, name in enumerate("abc"))

    # ACT
    kept = without_frames(entries, ["a_raw.fits", "c_raw.fits"], "no object trace was found in the B frames", log)

    # ASSERT
    assert [entry.name for entry in kept] == ["b.fits"]
    assert len(entries) == 3
    assert len(_warnings(log)) == 1
    assert "a.fits" in _warnings(log)[0] and "c.fits" in _warnings(log)[0]
    assert "no object trace was found in the B frames" in _warnings(log)[0]


def test_leaving_no_frames_out_gives_the_same_list_and_no_warning(log: Any) -> None:
    # ARRANGE
    entries = (_entry(log, "a"),)

    # ACT
    kept = without_frames(entries, [], "any reason", log)

    # ASSERT
    assert kept == entries
    assert not _warnings(log)


def test_only_the_positions_with_no_trace_leave_the_list_each_with_its_own_warning(log: Any) -> None:
    # ARRANGE
    entries = tuple(_entry(log, name, mjd=61000.0 + index) for index, name in enumerate("abcd"))
    positions = [
        (True, ["a_raw.fits"], "A frame"),
        (False, ["b_raw.fits"], "B frame"),
        (False, ["c_raw.fits"], "C frame"),
    ]

    # ACT
    kept = without_untraced_positions(entries, positions, log)

    # ASSERT
    assert _names(kept) == ["a.fits", "d.fits"]
    warnings = _warnings(log)
    assert len(warnings) == 2
    assert "no object trace was found in the B frame" in warnings[0] and "b.fits" in warnings[0]
    assert "no object trace was found in the C frame" in warnings[1] and "c.fits" in warnings[1]


def test_positions_that_all_have_a_trace_leave_the_list_alone(log: Any) -> None:
    # ARRANGE
    entries = (_entry(log, "a"),)

    # ACT / ASSERT
    assert without_untraced_positions(entries, [(True, ["a_raw.fits"], "A frame")], log) == entries
    assert not _warnings(log)


# THE CARDS


def test_the_exposure_and_the_total_exposure_are_the_sum_of_the_time_on_target(log: Any) -> None:
    # ARRANGE
    entries = tuple(_entry(log, name, mjd=61000.0 + index, exptime=600.0) for index, name in enumerate("abc"))

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert cards["EXPTIME"] == cards["TEXPTIME"] == 1800.0
    assert cards["NCOMBINE"] == 3


def test_the_observation_starts_with_the_first_frame_and_ends_when_the_last_frame_ends(log: Any) -> None:
    # ARRANGE
    entries = (
        _entry(log, "a", mjd=61000.0, exptime=600.0),
        _entry(log, "b", mjd=61000.5, exptime=900.0),
    )

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert cards["MJD-OBS"] == 61000.0
    assert cards["MJD-END"] == pytest.approx(61000.5 + 900.0 / SECONDS_PER_DAY, abs=1e-9)


def test_the_observation_ends_when_the_frame_that_ends_last_ends_even_if_it_did_not_start_last(log: Any) -> None:
    # ARRANGE: THE FIRST FRAME STARTS EARLIER AND LASTS LONGER, SO IT ENDS AFTER THE SECOND
    entries = (
        _entry(log, "a", mjd=61000.0, exptime=3600.0),
        _entry(log, "b", mjd=61000.01, exptime=300.0),
    )

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert cards["MJD-END"] == pytest.approx(61000.0 + 3600.0 / SECONDS_PER_DAY, abs=1e-9)


def test_the_total_exposure_is_not_longer_than_the_time_between_the_start_and_the_end(log: Any) -> None:
    # ARRANGE: THREE BACK-TO-BACK EXPOSURES WITH A GAP BETWEEN THEM
    entries = tuple(
        _entry(log, name, mjd=61000.0 + index * 1000.0 / SECONDS_PER_DAY, exptime=600.0)
        for index, name in enumerate("abc")
    )

    # ACT
    cards = _cards(entries)

    # ASSERT
    elapsedSeconds = (cards["MJD-END"] - cards["MJD-OBS"]) * SECONDS_PER_DAY
    assert cards["TEXPTIME"] <= elapsedSeconds + 1e-3


def test_a_nir_list_totals_dit_times_ndit_for_each_frame(log: Any) -> None:
    # ARRANGE
    entries = tuple(_entry(log, name, "NIR", mjd=61000.0 + index) for index, name in enumerate("ab"))

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert cards["EXPTIME"] == 3600.0


def test_the_observation_ids_are_listed_once_each_in_time_order(log: Any) -> None:
    # ARRANGE
    entries = (
        _entry(log, "a", mjd=61000.0, obsId=11),
        _entry(log, "b", mjd=61000.1, obsId=11),
        _entry(log, "c", mjd=61000.2, obsId=12),
    )

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert (cards["OBID1"], cards["OBID2"]) == (11, 12)
    assert "OBID3" not in cards


def test_one_programme_is_the_prog_id_and_has_no_progid_cards(log: Any) -> None:
    # ARRANGE
    entries = (_entry(log, "a", mjd=61000.0), _entry(log, "b", mjd=61000.1))

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert cards["PROG_ID"] == "115.ABCD.001"
    assert not [keyword for keyword in cards if keyword.startswith("PROGID")]


def test_several_programmes_make_the_prog_id_multi_and_are_listed_in_time_order(log: Any) -> None:
    # ARRANGE
    entries = (
        _entry(log, "a", mjd=61000.0, progId="115.X.001"),
        _entry(log, "b", mjd=61000.1, progId="115.Y.002"),
        _entry(log, "c", mjd=61000.2, progId="115.X.001"),
    )

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert cards["PROG_ID"] == "MULTI"
    assert (cards["PROGID1"], cards["PROGID2"]) == ("115.X.001", "115.Y.002")
    assert "PROGID3" not in cards


def test_each_frame_has_a_provenance_card_without_zero_padding_and_without_its_fits_suffix(log: Any) -> None:
    # ARRANGE
    entries = tuple(_entry(log, f"SOXS.{index}", mjd=61000.0 + index) for index in range(10))

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert cards["PROV1"] == "SOXS.0"
    assert cards["PROV9"] == "SOXS.8"
    assert cards["PROV10"] == "SOXS.9"
    assert "PROV01" not in cards and "PROV11" not in cards


def test_a_list_with_a_frame_of_unknown_time_gives_no_exposure_cards(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS", mjd=61000.1)
    del header["EXPTIME"]
    unknown = read_on_target_frame(header, _kw(log), "b_raw.fits", "VIS", log)
    entries = (_entry(log, "a", mjd=61000.0), unknown)

    # ACT
    cards = _cards(entries)

    # ASSERT
    assert not {"EXPTIME", "TEXPTIME", "MJD-END"} & set(cards)
    assert cards["NCOMBINE"] == 2


def test_an_empty_list_gives_no_cards() -> None:
    # ACT / ASSERT
    assert frame_accounting_cards(()) == []


# READING THE TARGET POSITION


def test_a_frame_reads_its_target_position_from_the_telescope_target_cards(log: Any) -> None:
    # ARRANGE: THE POINTING CARDS DIFFER FROM THE TARGET AND MUST NOT BE USED
    header = raw_science_header(
        "VIS",
        RA=7.58343,
        DEC=-33.25143,
        **{"HIERARCH ESO TEL TARG ALPHA": 3019.8, "HIERARCH ESO TEL TARG DELTA": -331501.08},
    )

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert entry.targetRa == pytest.approx(7.5825, abs=1e-9)
    assert entry.targetDec == pytest.approx(-33.2503, abs=1e-9)
    assert not _warnings(log)


def test_a_frame_with_no_target_cards_has_no_target_position_and_no_warning(log: Any) -> None:
    # ARRANGE: A CALIBRATION FRAME HAS NO TEL TARG CARDS, AND THE PRODUCT WARNS ONCE WHEN NO FRAME HAS THEM
    header = raw_science_header("VIS", RA=7.58343, DEC=-33.25143)

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert (entry.targetRa, entry.targetDec) == (None, None)
    assert not _warnings(log)


def test_a_target_card_that_does_not_decode_warns_with_the_frame_and_the_card(log: Any) -> None:
    # ARRANGE
    header = raw_science_header(
        "VIS", **{"HIERARCH ESO TEL TARG ALPHA": 3099.8, "HIERARCH ESO TEL TARG DELTA": -331501.08}
    )

    # ACT
    entry = read_on_target_frame(header, _kw(log), "frame.fits", "VIS", log)

    # ASSERT
    assert (entry.targetRa, entry.targetDec) == (None, pytest.approx(-33.2503, abs=1e-9))
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "frame.fits" in warnings[0] and "ESO TEL TARG ALPHA" in warnings[0]


# STAMPING A HEADER


def test_stamping_gives_a_copy_and_leaves_the_frame_header_unchanged(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS", exptime=600.0)
    before = header.copy()
    entries = (_entry(log, "a", mjd=61000.0, exptime=600.0), _entry(log, "b", mjd=61000.1, exptime=600.0))

    # ACT
    stamped = with_frame_accounting(header, entries)

    # ASSERT
    assert stamped is not header
    assert header == before
    assert header["EXPTIME"] == 600.0
    assert stamped["EXPTIME"] == 1200.0
    assert stamped["NCOMBINE"] == 2


def test_stamping_removes_the_accounting_cards_an_earlier_stamp_left(log: Any) -> None:
    # ARRANGE
    stale = raw_science_header(
        "VIS", PROV1="old1", PROV2="old2", PROV3="old3", OBID1=1, OBID2=2, PROGID1="old", PROG_ID="MULTI", NCOMBINE=3
    )
    stale["TEXPTIME"] = 99.0
    stale["MJD-END"] = 1.0
    entries = (_entry(log, "a", mjd=61000.0),)

    # ACT
    stamped = with_frame_accounting(stale, entries)

    # ASSERT
    assert stamped["PROV1"] == "a"
    assert not {"PROV2", "PROV3", "OBID2", "PROGID1"} & set(stamped)
    assert stamped["PROG_ID"] == "115.ABCD.001"
    assert stamped["NCOMBINE"] == 1


def test_stamping_a_list_of_unknown_time_removes_the_raw_exptime_so_it_cannot_pass_for_a_total(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS", exptime=600.0)
    noTime = raw_science_header("VIS")
    del noTime["EXPTIME"]
    entries = (read_on_target_frame(noTime, _kw(log), "a_raw.fits", "VIS", log),)

    # ACT
    stamped = with_frame_accounting(header, entries)

    # ASSERT
    assert "EXPTIME" not in stamped and "TEXPTIME" not in stamped


def test_stamping_replaces_the_pointing_with_the_target_position_and_sets_the_reference_frame(log: Any) -> None:
    # ARRANGE: A RAW HEADER WITH THE POINTING AND AN OLD REFERENCE FRAME
    header = raw_science_header("VIS", RA=7.58343, DEC=-33.25143, RADECSYS="FK4", EQUINOX=1950.0, OBJECT="SN2026ob")
    entries = (
        _entry(
            log, "a", mjd=61000.0, **{"HIERARCH ESO TEL TARG ALPHA": 3019.8, "HIERARCH ESO TEL TARG DELTA": -331501.08}
        ),
    )

    # ACT
    stamped = with_frame_accounting(header, entries, log)

    # ASSERT
    assert stamped["RA"] == pytest.approx(7.5825, abs=1e-9)
    assert stamped["DEC"] == pytest.approx(-33.2503, abs=1e-9)
    assert stamped.comments["RA"] == "[deg] Spectroscopic target position (J2000)"
    assert (stamped["RADESYS"], stamped["EQUINOX"], stamped["OBJECT"]) == ("FK5", 2000.0, "SN2026ob")
    assert "RADECSYS" not in stamped
    assert header["RA"] == 7.58343


def test_stamping_with_no_target_position_removes_the_pointing_and_warns_once(log: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS", RA=7.58343, DEC=-33.25143, RADECSYS="FK5", EQUINOX=2000.0)
    entries = (_entry(log, "a", mjd=61000.0), _entry(log, "b", mjd=61000.1))

    # ACT
    stamped = with_frame_accounting(header, entries, log)

    # ASSERT
    assert not {"RA", "DEC", "RADESYS", "RADECSYS", "EQUINOX"} & set(stamped)
    assert len([message for message in _warnings(log) if "ESO TEL TARG" in message]) == 1


def test_stamping_without_a_logger_does_not_lose_the_warning(log: Any, caplog: pytest.LogCaptureFixture) -> None:
    # ARRANGE
    header = raw_science_header("VIS", RA=7.58343, DEC=-33.25143)
    entry = read_on_target_frame(header, _kw(log), "a_raw.fits", "VIS", log)

    # ACT
    with caplog.at_level("WARNING"):
        stamped = with_frame_accounting(header, (entry,))

    # ASSERT
    assert "RA" not in stamped
    assert any("ESO TEL TARG" in record.getMessage() for record in caplog.records)


@pytest.mark.parametrize("entries", [None, ()])
def test_stamping_with_no_list_returns_the_header_itself(entries: Any) -> None:
    # ARRANGE
    header = raw_science_header("VIS")

    # ACT / ASSERT
    assert with_frame_accounting(header, entries) is header
