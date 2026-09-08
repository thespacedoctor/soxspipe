"""Characterization tests for instrument metadata lookups."""

from __future__ import annotations

import pytest

from soxspipe.commonutils import detector_lookup, keyword_lookup

pytestmark = pytest.mark.unit


def test_keyword_lookup_resolves_scalar_list_and_index(log: object) -> None:
    lookup = keyword_lookup(log=log, instrument="xsh").get

    assert lookup("DET_NDITSKIP") == "ESO DET NDITSKIP"
    assert lookup(["PROV", "DET_NDITSKIP"]) == ["PROV", "ESO DET NDITSKIP"]
    assert lookup("PROV", 1) == "PROV01"
    assert lookup("PROV", 14) == "PROV14"


def test_keyword_lookup_defaults_to_soxs(log: object) -> None:
    lookup = keyword_lookup(log=log).get

    assert lookup("DET_READ_SPEED") == "ESO DET READ CURID"


def test_keyword_lookup_rejects_unknown_alias(log: object) -> None:
    lookup = keyword_lookup(log=log, instrument="soxs").get

    with pytest.raises(
        LookupError,
        match="UNKNOWN is not in the list of known FITS Header keyword aliases",
    ):
        lookup("UNKNOWN")


def test_detector_lookup_is_case_insensitive_and_returns_expected_arm(
    log: object,
) -> None:
    lookup = detector_lookup(log=log, settings={"instrument": "soxs"}).get

    uppercase = lookup("NIR")
    lowercase = lookup("nir")

    assert lowercase == uppercase
    assert uppercase["binning"] == [1, 1]
    assert uppercase["dispersion-axis"] == "y"


def test_detector_lookup_rejects_unknown_arm(log: object) -> None:
    lookup = detector_lookup(log=log, settings={"instrument": "soxs"}).get

    with pytest.raises(
        LookupError,
        match="the detector 'RUBBISH' cannot be found",
    ):
        lookup("rubbish")
