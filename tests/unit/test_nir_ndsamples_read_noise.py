"""The SOXS NIR read noise as a function of the up-the-ramp sample count, ESO DET NDSAMPLES (DY-2341).

The table is measured once from raw dark pairs and stored in the SOXS detector parameters. A frame's read noise is
looked up from its own sample count, interpolated linearly in log(NDSAMPLES) between the measured points.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.read_noise import read_noise_from_ndsamples, table_read_noise

pytestmark = pytest.mark.unit

MEASURED_TABLE = {2: 17.50, 4: 16.94, 7: 14.62, 10: 13.14, 39: 8.95, 77: 7.62, 128: 7.55}
FALLBACK_RON = 7.55
FRAME_LABEL = "SOXS.2026-01-31T19:16:38.123.fits"


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_the_soxs_keyword_map_names_the_ndsamples_card(log: Any) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get

    # ACT
    card = kw("DET_NDSAMPLES")

    # ASSERT
    assert card == "ESO DET NDSAMPLES"


def test_the_soxs_nir_detector_parameters_hold_the_measured_table_and_the_128_sample_fallback(log: Any) -> None:
    # ARRANGE
    detector = detector_lookup(log=log, settings={"instrument": "soxs"})

    # ACT
    nir = detector.get("NIR")

    # ASSERT
    assert nir["read-noise-by-ndsamples"] == MEASURED_TABLE
    assert nir["ron"] == FALLBACK_RON


@pytest.mark.parametrize(("ndsamples", "expected"), sorted(MEASURED_TABLE.items()))
def test_a_measured_sample_count_gives_its_table_value(ndsamples: int, expected: float) -> None:
    # ARRANGE / ACT
    ron = table_read_noise(ndsamples, MEASURED_TABLE)

    # ASSERT
    assert ron == pytest.approx(expected)


def test_an_unmeasured_sample_count_is_interpolated_linearly_in_log_ndsamples() -> None:
    # ARRANGE: 97 LIES BETWEEN THE 77 AND 128 SAMPLE POINTS
    fraction = math.log(97 / 77) / math.log(128 / 77)
    expected = 7.62 + fraction * (7.55 - 7.62)

    # ACT
    ron = table_read_noise(97, MEASURED_TABLE)

    # ASSERT
    assert ron == pytest.approx(expected)
    assert 7.55 < ron < 7.62


@pytest.mark.parametrize("ndsamples", [None, 1, 0, -5, 129, 1000, True, 10.5, float("nan"), "ten"])
def test_a_missing_or_out_of_range_sample_count_has_no_table_value(ndsamples: object) -> None:
    # ARRANGE / ACT
    ron = table_read_noise(ndsamples, MEASURED_TABLE)

    # ASSERT
    assert ron is None


def test_a_whole_number_float_sample_count_is_accepted() -> None:
    # ARRANGE / ACT
    ron = table_read_noise(10.0, MEASURED_TABLE)

    # ASSERT
    assert ron == pytest.approx(13.14)


def test_a_usable_sample_count_gives_the_table_value_without_a_warning(log: Any) -> None:
    # ARRANGE / ACT
    ron = read_noise_from_ndsamples(10, MEASURED_TABLE, FALLBACK_RON, FRAME_LABEL, log)

    # ASSERT
    assert ron == pytest.approx(13.14)
    assert _warnings(log) == []


@pytest.mark.parametrize(
    ("ndsamples", "reason"),
    [
        (None, "no ESO DET NDSAMPLES"),
        (1, "outside the measured range of 2 to 128"),
        (129, "outside the measured range of 2 to 128"),
        ("ten", "not a whole number"),
    ],
    ids=["missing card", "below the range", "above the range", "not a number"],
)
def test_an_unusable_sample_count_falls_back_with_a_warning_naming_the_frame_and_the_reason(
    log: Any, ndsamples: object, reason: str
) -> None:
    # ARRANGE / ACT
    ron = read_noise_from_ndsamples(ndsamples, MEASURED_TABLE, FALLBACK_RON, FRAME_LABEL, log)

    # ASSERT
    assert ron == pytest.approx(FALLBACK_RON)
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert FRAME_LABEL in warnings[0]
    assert reason in warnings[0]
    assert "7.55" in warnings[0]
