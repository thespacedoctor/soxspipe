"""The read noise a UVB/VIS reduction uses, resolved from the master bias it names (DY-1458).

The master bias carries `ESO QC RAW RON`, the single-frame read noise measured from a pair of raw
frames, and `ESO PRO DATANCOM`, the number of frames it combined. A master bias that lacks either
is an old product, and the resolver falls back to the nominal header value, then the detector default.
"""

from __future__ import annotations

from typing import Any

import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.read_noise import ReadNoise, resolve_read_noise

pytestmark = pytest.mark.unit

NOMINAL_RON = 3.3
DEFAULT_RON = 3.8
MEASURED_RON = 4.0


def _kw(log: Any):
    return keyword_lookup(log=log, settings={"instrument": "soxs"}).get


def _header(ron: object = MEASURED_RON, nCombined: object = 4, **cards: object) -> fits.Header:
    """Return a master bias header; pass None for a card the product lacks."""
    header = fits.Header()
    if ron is not None:
        header["HIERARCH ESO QC RAW RON"] = ron
    if nCombined is not None:
        header["HIERARCH ESO PRO DATANCOM"] = nCombined
    for key, value in cards.items():
        header[key] = value
    return header


def _resolve(log: Any, header: fits.Header | None, nominal: float | None = NOMINAL_RON) -> ReadNoise:
    return resolve_read_noise(header, nominal, DEFAULT_RON, _kw(log), log)


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_the_measured_read_noise_replaces_the_nominal_one(log: Any) -> None:
    # ACT
    result = _resolve(log, _header(ron=4.0, nCombined=4))

    # ASSERT
    assert result.singleFrame == pytest.approx(4.0)
    assert result.nCombined == 4
    assert result.source == "master bias"


@pytest.mark.parametrize(
    ("ron", "reason"),
    [
        (None, "missing"),
        ("nan", "not finite"),
        ("inf", "not finite"),
        (0.0, "not positive"),
        (-1.0, "not positive"),
        ("garbage", "not a number"),
        (True, "a boolean"),
    ],
)
def test_an_unusable_measured_value_falls_back_to_the_nominal_value_with_a_warning(
    log: Any, ron: object, reason: str
) -> None:
    # ACT
    result = _resolve(log, _header(ron=ron))

    # ASSERT
    assert result.singleFrame == pytest.approx(NOMINAL_RON), reason
    assert result.source == "raw-frame header"
    assert len(_warnings(log)) == 1


def test_a_master_bias_without_datancom_is_an_old_product_and_falls_back_to_nominal(log: Any) -> None:
    # ARRANGE
    header = _header(ron=4.0, nCombined=None)

    # ACT
    result = _resolve(log, header)

    # ASSERT
    assert result.singleFrame == pytest.approx(NOMINAL_RON)
    assert result.source == "raw-frame header"
    assert any("DATANCOM" in message for message in _warnings(log))


def test_the_old_detector_ron_keyword_is_never_read(log: Any) -> None:
    # ARRANGE
    header = _header(ron=None, nCombined=None)
    header["HIERARCH ESO QC RON DETECTOR"] = 9.9

    # ACT
    result = _resolve(log, header)

    # ASSERT
    assert result.singleFrame == pytest.approx(NOMINAL_RON)


def test_no_master_bias_header_falls_back_to_nominal_with_a_warning_and_no_master_term(log: Any) -> None:
    # ACT
    result = _resolve(log, None)

    # ASSERT
    assert result.singleFrame == pytest.approx(NOMINAL_RON)
    assert result.skyFit == pytest.approx(NOMINAL_RON)
    assert result.nCombined is None
    assert len(_warnings(log)) == 1


def test_a_missing_nominal_value_falls_back_to_the_detector_default(log: Any) -> None:
    # ACT
    result = _resolve(log, _header(ron=None), nominal=None)

    # ASSERT
    assert result.singleFrame == pytest.approx(DEFAULT_RON)
    assert result.source == "detector default"
    assert len(_warnings(log)) == 1


def test_a_measured_value_over_fifteen_percent_from_nominal_warns_but_is_used(log: Any) -> None:
    # ACT: 4.0 IS 21% ABOVE 3.3
    result = _resolve(log, _header(ron=4.0))

    # ASSERT
    assert result.singleFrame == pytest.approx(4.0)
    assert len(_warnings(log)) == 1
    assert "nominal" in _warnings(log)[0]


def test_a_measured_value_within_fifteen_percent_of_nominal_does_not_warn(log: Any) -> None:
    # ACT: 3.6 IS 9% ABOVE 3.3
    result = _resolve(log, _header(ron=3.6))

    # ASSERT
    assert result.singleFrame == pytest.approx(3.6)
    assert _warnings(log) == []


def test_the_combined_frame_count_falls_back_to_the_raw_frame_cards(log: Any) -> None:
    # ARRANGE: AN OLD PRODUCT, SO THE FALLBACK IS USED, BUT N IS STILL KNOWN FROM THE RAW CARDS
    header = _header(ron=4.0, nCombined=None)
    for index in (1, 2, 3):
        header[f"HIERARCH ESO PRO REC1 RAW{index} NAME"] = f"raw_{index}.fits"

    # ACT
    result = _resolve(log, header)

    # ASSERT
    assert result.nCombined == 3
    assert result.skyFit == pytest.approx(3.3 * 1.154701, rel=1e-5)


def test_the_combined_frame_count_defaults_to_one_with_a_warning(log: Any) -> None:
    # ACT
    result = _resolve(log, _header(ron=4.0, nCombined=None))

    # ASSERT
    assert result.nCombined == 1
    assert result.skyFit == pytest.approx(3.3 * 1.414214, rel=1e-5)
    assert any("combined" in message for message in _warnings(log))


@pytest.mark.parametrize(
    ("nCombined", "expectedSkyFit"),
    [(1, 5.656854), (4, 4.472136), (9, 4.216370)],
)
def test_the_sky_fit_noise_adds_the_master_bias_noise_in_quadrature(
    log: Any, nCombined: int, expectedSkyFit: float
) -> None:
    # ACT: WORKED BY HAND, 4.0 * SQRT(1 + 1/N)
    result = _resolve(log, _header(ron=4.0, nCombined=nCombined))

    # ASSERT
    assert result.singleFrame == pytest.approx(4.0)
    assert result.skyFit == pytest.approx(expectedSkyFit, rel=1e-6)
