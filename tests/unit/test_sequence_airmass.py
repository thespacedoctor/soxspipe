"""Unit tests for the airmass of a frame sequence that the extinction correction uses (DY-1695)."""

from __future__ import annotations

from typing import Any

import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.sequence_airmass import mean_airmass

pytestmark = pytest.mark.unit


def _header(**cards: object) -> fits.Header:
    header = fits.Header()
    names = {
        "start": "HIERARCH ESO TEL AIRM START",
        "end": "HIERARCH ESO TEL AIRM END",
        "exptime": "EXPTIME",
        "seq": "HIERARCH ESO DET SEQ1 EXPTIME",
    }
    for name, value in cards.items():
        header[names[name]] = value
    return header


@pytest.fixture
def soxs_kw(log: Any) -> Any:
    return keyword_lookup(log=log, settings={"instrument": "soxs"}).get


def test_each_frame_contributes_the_midpoint_of_its_start_and_end_airmass(log: Any, soxs_kw: Any) -> None:
    # ARRANGE: FRAME MIDPOINTS 1.2 AND 1.6
    headers = [
        _header(start=1.1, end=1.3, exptime=300.0),
        _header(start=1.5, end=1.7, exptime=300.0),
    ]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.4)


def test_frames_are_weighted_by_their_exposure_time(log: Any, soxs_kw: Any) -> None:
    # ARRANGE: 300 S AT 1.2 AND 100 S AT 1.6 GIVE (360 + 160) / 400 = 1.3
    headers = [
        _header(start=1.2, end=1.2, exptime=300.0),
        _header(start=1.6, end=1.6, exptime=100.0),
    ]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.3)


def test_a_nir_frame_is_weighted_by_its_sequence_duration_not_one_dit(log: Any, soxs_kw: Any) -> None:
    # ARRANGE: EQUAL DITS, BUT THE FIRST FRAME RAN 3 TIMES LONGER, SO (3 * 1.2 + 1.6) / 4 = 1.3
    headers = [
        _header(start=1.2, end=1.2, exptime=600.0, seq=2400.0),
        _header(start=1.6, end=1.6, exptime=600.0, seq=800.0),
    ]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.3)


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_a_frame_without_a_start_airmass_falls_back_to_its_end_airmass_and_says_so(log: Any, soxs_kw: Any) -> None:
    # ARRANGE
    headers = [
        _header(end=1.3, exptime=300.0),
        _header(start=1.5, end=1.7, exptime=300.0),
    ]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT: (1.3 + 1.6) / 2
    assert airmass == pytest.approx(1.45)
    assert any("ESO TEL AIRM START" in message for message in _warnings(log))


@pytest.mark.parametrize("badValue", ["N/A", ""])
def test_an_unreadable_start_airmass_counts_as_missing(log: Any, soxs_kw: Any, badValue: object) -> None:
    # ARRANGE
    headers = [_header(start=badValue, end=1.3, exptime=300.0)]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.3)
    assert _warnings(log)


def test_a_frame_without_an_end_airmass_falls_back_to_its_start_airmass(log: Any, soxs_kw: Any) -> None:
    # ARRANGE
    headers = [_header(start=1.1, exptime=300.0)]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.1)
    assert any("ESO TEL AIRM END" in message for message in _warnings(log))


def test_a_frame_with_no_airmass_is_left_out_of_the_mean(log: Any, soxs_kw: Any) -> None:
    # ARRANGE
    headers = [_header(exptime=300.0), _header(start=1.5, end=1.7, exptime=300.0)]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.6)
    assert _warnings(log)


def test_no_readable_airmass_in_any_frame_gives_none(log: Any, soxs_kw: Any) -> None:
    # ARRANGE
    headers = [_header(exptime=300.0), _header(exptime=300.0)]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass is None
    assert _warnings(log)


def test_an_unreadable_exposure_time_gives_every_frame_equal_weight(log: Any, soxs_kw: Any) -> None:
    # ARRANGE: WITHOUT A DURATION ON THE SECOND FRAME, THE WEIGHTS ARE DROPPED FOR ALL FRAMES
    headers = [
        _header(start=1.2, end=1.2, exptime=300.0),
        _header(start=1.6, end=1.6),
    ]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.4)
    assert any("equal weight" in message for message in _warnings(log))


def test_xshooter_frames_without_a_sequence_duration_alias_are_weighted_by_exptime(log: Any) -> None:
    # ARRANGE
    xshKw = keyword_lookup(log=log, settings={"instrument": "xsh"}).get
    headers = [
        _header(start=1.2, end=1.2, exptime=300.0),
        _header(start=1.6, end=1.6, exptime=100.0),
    ]

    # ACT
    airmass = mean_airmass(headers, xshKw, log)

    # ASSERT
    assert airmass == pytest.approx(1.3)


@pytest.mark.parametrize("sentinel", [-99.99, 0.0, 0.5])
def test_an_airmass_below_one_counts_as_missing(log: Any, soxs_kw: Any, sentinel: float) -> None:
    # ARRANGE: AN AIRMASS IS AT LEAST 1, AND THE DATA ORGANISER FILLS MISSING VALUES WITH -99.99
    headers = [_header(start=sentinel, end=1.3, exptime=300.0)]

    # ACT
    airmass = mean_airmass(headers, soxs_kw, log)

    # ASSERT
    assert airmass == pytest.approx(1.3)
    assert _warnings(log)
