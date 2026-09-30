"""Characterization tests for instrument metadata lookups."""

from __future__ import annotations

import decimal

import numpy
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


# CHARACTERIZATION TESTS PINNING CURRENT `get()` INDEX-FORMATTING BEHAVIOUR
# (DY-63). THESE PIN THE EXISTING BEHAVIOUR EXACTLY, INCLUDING QUIRKS, SO A
# LATER REFACTOR CAN BE VERIFIED AGAINST THEM. DO NOT "FIX" A SURPRISING ROW.


@pytest.mark.parametrize(
    "index",
    [0, False, None, ""],
    ids=["zero", "false", "none", "empty_string"],
)
def test_keyword_lookup_get_omits_suffix_for_falsy_index(
    log: object, index: object
) -> None:
    lookup = keyword_lookup(log=log, instrument="xsh").get

    result = lookup("PROV", index)

    assert result == "PROV"


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        (10, "PROV10"),
        (True, "PROV01"),
        (100, "PROV100"),
        (-3, "PROV-03"),
    ],
    ids=["two_digit", "true", "three_digit", "negative_int"],
)
def test_keyword_lookup_get_zero_pads_integer_index(
    log: object, index: object, expected: str
) -> None:
    lookup = keyword_lookup(log=log, instrument="xsh").get

    result = lookup("PROV", index)

    assert result == expected


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        (7.9, "PROV07"),
        (-7.9, "PROV-07"),
    ],
    ids=["positive_float_truncates", "negative_float_truncates_toward_zero"],
)
def test_keyword_lookup_get_truncates_float_index_toward_zero(
    log: object, index: float, expected: str
) -> None:
    lookup = keyword_lookup(log=log, instrument="xsh").get

    result = lookup("PROV", index)

    assert result == expected


@pytest.mark.parametrize(
    ("index", "expected"),
    [
        (decimal.Decimal("7.9"), "PROV07"),
        (numpy.int64(5), "PROV05"),
        (numpy.float32(2.5), "PROV02"),
    ],
    ids=["decimal", "numpy_int64", "numpy_float32"],
)
def test_keyword_lookup_get_accepts_decimal_and_numpy_numeric_index(
    log: object, index: object, expected: str
) -> None:
    lookup = keyword_lookup(log=log, instrument="xsh").get

    result = lookup("PROV", index)

    assert result == expected


@pytest.mark.parametrize(
    ("index", "exceptionType"),
    [
        ("3", TypeError),
        (b"3", TypeError),
        (bytearray(b"3"), TypeError),
        (memoryview(b"3"), TypeError),
        ([1], TypeError),
        (1j, TypeError),
        (float("nan"), ValueError),
        (float("inf"), OverflowError),
    ],
    ids=["numeric_string", "numeric_bytes", "numeric_bytearray", "numeric_memoryview", "list", "complex", "nan", "inf"],
)
def test_keyword_lookup_get_raises_on_unformattable_index(
    log: object, index: object, exceptionType: type[Exception]
) -> None:
    lookup = keyword_lookup(log=log, instrument="xsh").get

    with pytest.raises(exceptionType):
        lookup("PROV", index)


def test_keyword_lookup_get_applies_index_to_every_tag_in_a_list(
    log: object,
) -> None:
    lookup = keyword_lookup(log=log, instrument="xsh").get

    result = lookup(["PROV", "PROV"], 3)

    assert result == ["PROV03", "PROV03"]


def test_keyword_lookup_uses_instrument_from_settings(log: object) -> None:
    lookup = keyword_lookup(log=log, settings={"instrument": "xsh"})

    assert lookup.instrument == "xsh"
    assert lookup.get("DET_NDITSKIP") == "ESO DET NDITSKIP"


def test_keyword_lookup_explicit_instrument_wins_over_settings(
    log: object,
) -> None:
    lookup = keyword_lookup(
        log=log, instrument="soxs", settings={"instrument": "xsh"}
    )

    assert lookup.instrument == "soxs"


def test_keyword_lookup_falls_back_to_soxs_when_unspecified(log: object) -> None:
    lookup = keyword_lookup(log=log)

    assert lookup.instrument == "soxs"


def test_keyword_lookup_lowercases_instrument_name_for_dictionary_selection(
    log: object,
) -> None:
    lookup = keyword_lookup(log=log, instrument="XSH").get

    assert lookup("DET_NDITSKIP") == "ESO DET NDITSKIP"
