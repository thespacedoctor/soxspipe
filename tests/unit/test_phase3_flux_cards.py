"""The flux-calibration cards of a merged spectrum, and a header scrub that keeps the target position (DY-2357)."""

from __future__ import annotations

from typing import Any

import pytest
from astropy.io import fits

from soxspipe.commonutils.phase3 import basic_header_scrubbing, flux_calibration_cards

pytestmark = pytest.mark.unit


def _values(cards: list[tuple[str, Any, str]]) -> dict[str, Any]:
    return {keyword: value for keyword, value, _ in cards}


def test_a_calibrated_spectrum_is_absolute_with_a_flux_error_of_minus_two() -> None:
    # ACT
    values = _values(flux_calibration_cards(isCalibrated=True))

    # ASSERT
    assert values == {"FLUXCAL": "ABSOLUTE", "TOT_FLUX": False, "FLUXERR": -2.0}
    assert isinstance(values["FLUXERR"], float)


def test_a_spectrum_in_counts_is_uncalibrated_with_a_flux_error_of_minus_one() -> None:
    # ACT
    values = _values(flux_calibration_cards(isCalibrated=False))

    # ASSERT
    assert values == {"FLUXCAL": "UNCALIBRATED", "TOT_FLUX": False, "FLUXERR": -1.0}
    assert isinstance(values["FLUXERR"], float)


@pytest.mark.parametrize("isCalibrated", [True, False])
def test_no_flux_calibration_card_lacks_a_comment(isCalibrated: bool) -> None:
    # ACT / ASSERT
    assert all(comment for _, _, comment in flux_calibration_cards(isCalibrated))


def test_the_scrub_keeps_an_equinox_card_that_is_already_there(log: Any) -> None:
    # ARRANGE: THE TARGET POSITION SET EQUINOX, AND THE RAW TEL TARG EQUINOX IS STILL IN THE HEADER
    header = fits.Header()
    header["EQUINOX"] = (2000.0, "[yr] Equinox of RA and DEC")
    header["HIERARCH ESO TEL TARG EQUINOX"] = 1950.0
    header["RADESYS"] = ("FK5", "Reference frame of RA and DEC")
    header["RADECSYS"] = "ICRS"

    # ACT
    scrubbed = basic_header_scrubbing(log, {}, header)

    # ASSERT
    assert scrubbed["EQUINOX"] == 2000.0
    assert scrubbed.comments["EQUINOX"] == "[yr] Equinox of RA and DEC"
    assert scrubbed["RADESYS"] == "FK5"
    assert "ESO TEL TARG EQUINOX" not in scrubbed
    assert "RADECSYS" not in scrubbed


def test_the_scrub_still_renames_the_old_cards_when_the_new_ones_are_absent(log: Any) -> None:
    # ARRANGE
    header = fits.Header()
    header["HIERARCH ESO TEL TARG EQUINOX"] = (2000.0, "equinox of the pointing")
    header["RADECSYS"] = "ICRS"

    # ACT
    scrubbed = basic_header_scrubbing(log, {}, header)

    # ASSERT
    assert (scrubbed["EQUINOX"], scrubbed["RADESYS"]) == (2000.0, "ICRS")
    assert scrubbed.comments["EQUINOX"] == "equinox of the pointing"
    assert "RADECSYS" not in scrubbed and "ESO TEL TARG EQUINOX" not in scrubbed
