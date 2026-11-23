"""Inferring the host redshift from the host-line windows found in a stare frame, and the companion windows (DY-2242).

Two or more windows that agree on one redshift for distinct lines of a short rest-frame list give the redshift. The
other lines of the list that fall in an order are then masked too, because a weak line such as [NII] 6548 is below
the detection limit on its own. Several pairs of lines in the list have almost the same wavelength ratio (for example
[NII] 6548 and H-alpha against the two [SII] lines), so two windows can fit two redshifts. The inference then gives
none.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import host_line_mask
from soxspipe.commonutils.host_line_mask import companion_windows, infer_host_redshift, load_host_lines
from soxspipe.commonutils.nod_host_contamination import host_feature

pytestmark = pytest.mark.unit

Z = 0.0166
REST_NM = dict(load_host_lines())
HALPHA_NM = REST_NM["Halpha"]
NII_6583_NM = REST_NM["[NII] 6583"]
NII_6548_NM = REST_NM["[NII] 6548"]
SII_6716_NM = REST_NM["[SII] 6716"]
SII_6731_NM = REST_NM["[SII] 6731"]
HALF_WIDTH_NM = 0.25
STEP_NM = 0.04
# THE CENTRE OF A WINDOW IS GOOD TO ABOUT 0.1 NM, WHICH IS ABOUT 1.5E-4 IN REDSHIFT AT 670 NM
JITTER = 1.5e-4


def _window(restNm: float, z: float = Z, order: int = 1, halfWidthNm: float = HALF_WIDTH_NM) -> host_feature:
    centre = restNm * (1 + z)
    return host_feature(order, centre - halfWidthNm, centre + halfWidthNm, 9.0)


def _order(firstNm: float, lastNm: float, orderNumber: int = 1) -> pd.DataFrame:
    """One synthetic order, 20 slit pixels wide, with every column that the usable-pixel rule reads."""
    columns = np.arange(int((lastNm - firstNm) / STEP_NM))
    frame = pd.DataFrame(
        {
            "x": np.repeat(np.arange(20), len(columns)),
            "y": np.tile(columns, 20),
            "order": float(orderNumber),
            "wavelength": np.tile(firstNm + STEP_NM * columns, 20),
            "slit_position": np.repeat(np.linspace(-2.0, 2.0, 20), len(columns)),
            "flux": 20.0,
            "error": 4.0,
        }
    )
    for flag in ("flagged_object_clipped", "flagged_bad_pixel_clipped", "flagged_edge_clipped"):
        frame[flag] = False
    return frame


# ----- THE REST-FRAME LINE LIST -----


def test_the_rest_line_list_holds_the_strong_host_lines_in_air_wavelengths() -> None:
    assert REST_NM["Halpha"] == pytest.approx(656.28)
    assert {"[NII] 6548", "[NII] 6583", "[SIII] 9069", "[SIII] 9531", "Paalpha"} <= set(REST_NM)
    assert all(wavelength > 0 for wavelength in REST_NM.values())


def test_a_missing_line_list_raises_a_value_error_naming_the_file(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    monkeypatch.setattr(host_line_mask, "HOST_LINES_FILE", tmp_path / "missing.yaml")

    with pytest.raises(ValueError, match="missing.yaml"):
        load_host_lines()


@pytest.mark.parametrize(
    "text",
    [
        "lines: 5\n",
        "other: []\n",
        "lines:\n  - {name: Halpha}\n",
        "lines:\n  - {name: Halpha, wavelength-nm: red}\n",
        "lines: [unclosed\n",
        # A YAML INTEGER TOO LARGE FOR A FLOAT; A FLOAT LITERAL SUCH AS 1E400 BECOMES INF AND DOES NOT OVERFLOW
        "lines:\n  - {name: Halpha, wavelength-nm: 1" + "0" * 400 + "}\n",
        "lines:\n  - {name: Halpha, wavelength-nm: 1e400}\n",
        "lines:\n  - {name: Halpha, wavelength-nm: .nan}\n",
        b"lines:\n  - {name: H\xe9alpha, wavelength-nm: 656.28}\n",
    ],
    ids=[
        "not-a-list",
        "no-lines-key",
        "no-wavelength",
        "wavelength-not-a-number",
        "invalid-yaml",
        "wavelength-overflows-a-float",
        "wavelength-infinite",
        "wavelength-nan",
        "not-utf-8",
    ],
)
def test_a_malformed_line_list_raises_a_value_error_naming_the_file(
    monkeypatch: pytest.MonkeyPatch, tmp_path, text: str | bytes
) -> None:
    path = tmp_path / "host_emission_lines.yaml"
    if isinstance(text, bytes):
        path.write_bytes(text)
    else:
        path.write_text(text)
    monkeypatch.setattr(host_line_mask, "HOST_LINES_FILE", path)

    with pytest.raises(ValueError, match="host_emission_lines.yaml"):
        load_host_lines()


# ----- THE REDSHIFT -----


def test_two_windows_that_agree_on_a_redshift_give_it() -> None:
    redshift = infer_host_redshift([_window(HALPHA_NM), _window(NII_6583_NM)], load_host_lines())

    assert redshift is not None
    assert redshift.z == pytest.approx(Z, abs=1e-5)
    assert redshift.lineNames == ("Halpha", "[NII] 6583")
    assert redshift.nWindows == 2


def test_the_redshift_does_not_depend_on_the_order_of_the_windows() -> None:
    windows = [_window(HALPHA_NM), _window(NII_6583_NM), _window(953.1)]

    forward = infer_host_redshift(windows, load_host_lines())
    backward = infer_host_redshift(windows[::-1], load_host_lines())

    assert forward == backward


def test_one_window_gives_no_redshift() -> None:
    assert infer_host_redshift([_window(HALPHA_NM)], load_host_lines()) is None


def test_no_window_gives_no_redshift() -> None:
    assert infer_host_redshift([], load_host_lines()) is None


def test_no_rest_line_gives_no_redshift() -> None:
    assert infer_host_redshift([_window(HALPHA_NM), _window(NII_6583_NM)], ()) is None


def test_two_windows_that_do_not_agree_on_a_redshift_give_none() -> None:
    windows = [host_feature(1, 600.0, 600.5, 9.0), host_feature(2, 701.3, 701.8, 9.0)]

    assert infer_host_redshift(windows, load_host_lines()) is None


@pytest.mark.parametrize("sign", [1, -1])
def test_two_windows_that_differ_by_less_than_the_tolerance_give_the_median_redshift(sign: int) -> None:
    offset = sign * 4e-4

    redshift = infer_host_redshift([_window(HALPHA_NM), _window(953.1, z=Z + offset)], load_host_lines())

    assert redshift is not None
    assert redshift.z == pytest.approx(Z + offset / 2, abs=1e-6)


@pytest.mark.parametrize("sign", [1, -1])
def test_two_windows_that_differ_by_more_than_the_tolerance_give_none(sign: int) -> None:
    offset = sign * 6e-4

    assert infer_host_redshift([_window(HALPHA_NM), _window(953.1, z=Z + offset)], load_host_lines()) is None


def test_two_windows_on_the_same_line_are_not_two_lines() -> None:
    windows = [_window(HALPHA_NM, order=1), _window(HALPHA_NM, order=2)]

    assert infer_host_redshift(windows, load_host_lines()) is None


def test_the_redshift_with_the_most_windows_wins() -> None:
    windows = [_window(HALPHA_NM), _window(NII_6583_NM), _window(953.1), host_feature(1, 860.0, 860.4, 9.0)]

    redshift = infer_host_redshift(windows, load_host_lines())

    assert redshift is not None
    assert redshift.lineNames == ("Halpha", "[NII] 6583", "[SIII] 9531")


def test_the_half_width_fraction_is_the_median_of_the_matched_windows() -> None:
    windows = [
        _window(HALPHA_NM, halfWidthNm=0.2),
        _window(NII_6583_NM, halfWidthNm=0.4),
        _window(953.1, halfWidthNm=0.3),
    ]

    redshift = infer_host_redshift(windows, load_host_lines())

    assert redshift is not None
    assert redshift.halfWidthFraction == pytest.approx(0.3 / (953.1 * (1 + Z)), rel=1e-3)


# ----- PAIRS OF LINES THAT HAVE ALMOST THE SAME WAVELENGTH RATIO -----

# EACH PAIR OF LINES HAS A SECOND PAIR IN THE LIST WHOSE WAVELENGTH RATIO DIFFERS BY LESS THAN THE TOLERANCE
ALIAS_PAIRS = [
    pytest.param(("[NII] 6548", "Halpha"), ("[SII] 6716", "[SII] 6731"), id="NII6548-Halpha-vs-SII"),
    pytest.param(("[SII] 6716", "[SII] 6731"), ("[NII] 6548", "Halpha"), id="SII-vs-NII6548-Halpha"),
    pytest.param(("Hbeta", "[OIII] 4959"), ("[NII] 6583", "[SII] 6716"), id="Hbeta-OIII-vs-NII6583-SII"),
    pytest.param(("Hbeta", "Halpha"), ("[SII] 6716", "[SIII] 9069"), id="Hbeta-Halpha-vs-SII-SIII9069"),
]


def _alias_redshift(pair: tuple[str, str], alias: tuple[str, str]) -> float:
    """A redshift of the pair at which the alias pair gives a redshift inside the limits too."""
    for z in (Z, 0.05, 0.1, 0.2, 0.4):
        aliasZ = (1 + z) * REST_NM[pair[0]] / REST_NM[alias[0]] - 1
        if host_line_mask.MIN_REDSHIFT <= aliasZ <= host_line_mask.MAX_REDSHIFT:
            return z
    raise AssertionError(f"no redshift of {pair} puts {alias} inside the limits")


@pytest.mark.parametrize(("pair", "alias"), ALIAS_PAIRS)
@pytest.mark.parametrize("jitter", [(JITTER, -JITTER), (-JITTER, JITTER), (0.0, 0.0)], ids=["+-", "-+", "none"])
def test_two_windows_that_fit_two_redshifts_equally_give_none(
    pair: tuple[str, str], alias: tuple[str, str], jitter: tuple[float, float]
) -> None:
    z = _alias_redshift(pair, alias)
    windows = [_window(REST_NM[name], z=z + shift) for name, shift in zip(pair, jitter, strict=True)]
    # THE PAIR REALLY DOES FIT THE ALIAS REDSHIFT TOO, OR THE TEST CHECKS NOTHING
    centres = [(window.wavelengthMinNm + window.wavelengthMaxNm) / 2 for window in windows]
    aliasZ = [centre / REST_NM[name] - 1 for centre, name in zip(centres, alias, strict=True)]
    assert abs(aliasZ[0] - aliasZ[1]) < host_line_mask.REDSHIFT_TOLERANCE
    assert all(host_line_mask.MIN_REDSHIFT <= value <= host_line_mask.MAX_REDSHIFT for value in aliasZ)

    assert infer_host_redshift(windows, load_host_lines()) is None


@pytest.mark.parametrize(("pair", "alias"), ALIAS_PAIRS)
def test_a_third_window_on_the_true_redshift_separates_the_aliases(
    pair: tuple[str, str], alias: tuple[str, str]
) -> None:
    z = _alias_redshift(pair, alias)
    windows = [_window(REST_NM[name], z=z + shift) for name, shift in zip(pair, (JITTER, -JITTER), strict=True)]
    windows.append(_window(REST_NM["[SIII] 9531"], z=z))

    redshift = infer_host_redshift(windows, load_host_lines())

    assert redshift is not None
    assert redshift.z == pytest.approx(z, abs=3 * JITTER)
    assert set(pair) <= set(redshift.lineNames)
    assert redshift.nWindows == 3


def test_the_alias_that_is_outside_the_redshift_limits_does_not_make_a_pair_ambiguous() -> None:
    # [NII] 6548 AND H-ALPHA AT Z = 0.0166 FIT THE [SII] PAIR ONLY AT Z = -0.009, WHICH IS BELOW THE LOWER LIMIT
    windows = [_window(NII_6548_NM), _window(HALPHA_NM)]

    redshift = infer_host_redshift(windows, load_host_lines())

    assert redshift is not None
    assert redshift.lineNames == ("[NII] 6548", "Halpha")


# ----- THE LIMITS OF THE REDSHIFT -----


def test_a_redshift_inside_the_limits_is_returned() -> None:
    z = host_line_mask.MAX_REDSHIFT - 1e-3

    redshift = infer_host_redshift([_window(HALPHA_NM, z=z), _window(953.1, z=z)], load_host_lines())

    assert redshift is not None
    assert redshift.z == pytest.approx(z, abs=1e-6)


def test_a_rival_redshift_outside_the_limits_does_not_make_the_inside_one_ambiguous() -> None:
    # THE RIVAL IS SUPPORTED BY TWO WINDOWS AND ITS SEED IS INSIDE THE LIMIT, BUT ITS MEDIAN IS NOT
    z = host_line_mask.MAX_REDSHIFT
    rival = [_window(REST_NM["Hbeta"], z=z - 1e-4), _window(REST_NM["[OIII] 5007"], z=z + 3e-4)]
    windows = [_window(HALPHA_NM), _window(953.1), *rival]

    redshift = infer_host_redshift(windows, load_host_lines())

    assert redshift is not None
    assert redshift.z == pytest.approx(Z, abs=1e-5)
    assert redshift.lineNames == ("Halpha", "[SIII] 9531")


def test_a_redshift_above_the_upper_limit_gives_none_even_when_its_seed_is_inside() -> None:
    # THE SEED OF THE FIRST WINDOW IS INSIDE THE LIMIT, BUT THE MEDIAN OF THE TWO WINDOWS IS NOT
    z = host_line_mask.MAX_REDSHIFT
    windows = [_window(HALPHA_NM, z=z - 1e-4), _window(953.1, z=z + 3e-4)]

    assert infer_host_redshift(windows, load_host_lines()) is None


def test_a_redshift_below_the_lower_limit_gives_none_even_when_its_seed_is_inside() -> None:
    z = host_line_mask.MIN_REDSHIFT
    windows = [_window(HALPHA_NM, z=z + 1e-4), _window(953.1, z=z - 3e-4)]

    assert infer_host_redshift(windows, load_host_lines()) is None


# ----- THE COMPANION WINDOWS -----


def _companions(order: pd.DataFrame, detected: list[host_feature], skylines=()) -> tuple[host_feature, ...]:
    redshift = infer_host_redshift(detected, load_host_lines())
    assert redshift is not None
    return companion_windows(order, detected, redshift, load_host_lines(), skylines)


def _companion_half_width(detected: list[host_feature], centreNm: float) -> float:
    """The half width the rule gives: the detected relative width and the redshift uncertainty, in quadrature."""
    redshift = infer_host_redshift(detected, load_host_lines())
    assert redshift is not None
    return math.hypot(redshift.halfWidthFraction, host_line_mask.REDSHIFT_UNCERTAINTY) * centreNm


def _bounds(windows: tuple[host_feature, ...]) -> list[tuple[float, float]]:
    return sorted((round(window.wavelengthMinNm, 1), round(window.wavelengthMaxNm, 1)) for window in windows)


def test_a_companion_line_gets_a_window_that_allows_for_the_redshift_uncertainty() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]

    (companion,) = _companions(_order(660.0, 668.0), detected)

    centre = NII_6548_NM * (1 + Z)
    assert companion.wavelengthMinNm < centre < companion.wavelengthMaxNm
    halfWidth = (companion.wavelengthMaxNm - companion.wavelengthMinNm) / 2
    assert halfWidth == pytest.approx(_companion_half_width(detected, centre), rel=1e-3)
    assert halfWidth > HALF_WIDTH_NM
    assert companion.order == 1


def test_a_detected_companion_gets_no_second_window_but_another_companion_does() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM), _window(NII_6548_NM)]

    windows = _companions(_order(655.0, 690.0), detected)

    centres = sorted(round((window.wavelengthMinNm + window.wavelengthMaxNm) / 2, 1) for window in windows)
    assert centres == sorted(round(rest * (1 + Z), 1) for rest in (SII_6716_NM, SII_6731_NM))


def test_a_companion_whose_window_only_partly_overlaps_a_detected_window_gets_no_window() -> None:
    centre = NII_6548_NM * (1 + Z)
    # THE DETECTED WINDOW REACHES THE EDGE OF THE COMPANION WINDOW BUT NOT ITS CENTRE
    edge = centre + 0.15
    detected = [_window(HALPHA_NM), _window(NII_6583_NM), host_feature(1, edge, edge + 0.5, 9.0)]

    assert _companions(_order(660.0, 668.0), detected) == ()


def test_a_companion_line_outside_the_order_gets_no_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]

    assert _companions(_order(669.0, 680.0), detected) == ()


def test_a_companion_line_at_the_edge_of_the_order_gets_no_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]
    centre = NII_6548_NM * (1 + Z)

    assert _companions(_order(centre - 0.1, centre + 5.0), detected) == ()


def test_the_searched_range_comes_from_the_usable_pixels_only() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]
    order = _order(655.0, 690.0)
    # THE PIXELS OF THE RED HALF ARE BAD, SO THE ORDER ENDS AT 668 NM AS FAR AS THE SKY FIT IS CONCERNED
    order.loc[order["wavelength"] > 668.0, "flagged_bad_pixel_clipped"] = True

    windows = _companions(order, detected)

    centres = [round((window.wavelengthMinNm + window.wavelengthMaxNm) / 2, 1) for window in windows]
    assert centres == [round(NII_6548_NM * (1 + Z), 1)]


def test_pixels_with_a_wavelength_that_is_not_a_number_do_not_set_the_searched_range() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]
    order = _order(660.0, 668.0)
    order.loc[order.index[::50], "wavelength"] = np.nan

    assert len(_companions(order, detected)) == 1


def test_an_order_with_no_usable_pixel_gets_no_window_and_no_warning() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]
    order = _order(655.0, 690.0).assign(wavelength=np.nan)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        windows = _companions(order, detected)

    assert windows == ()


def test_a_companion_line_next_to_a_catalogue_skyline_gets_no_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]
    centre = NII_6548_NM * (1 + Z)

    assert _companions(_order(660.0, 668.0), detected, skylines=[centre + 0.3]) == ()


def test_a_companion_line_with_a_catalogue_skyline_far_away_still_gets_its_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]
    centre = NII_6548_NM * (1 + Z)

    assert len(_companions(_order(660.0, 668.0), detected, skylines=[centre + 3.0])) == 1


def test_every_companion_line_in_a_wide_order_gets_a_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]

    windows = _companions(_order(655.0, 690.0), detected)

    centres = sorted(round((window.wavelengthMinNm + window.wavelengthMaxNm) / 2, 1) for window in windows)
    assert centres == sorted(round(rest * (1 + Z), 1) for rest in (NII_6548_NM, SII_6716_NM, SII_6731_NM))


# ----- WINDOWS THAT ARE NEW -----


def test_a_window_inside_a_detected_window_is_not_new() -> None:
    detected = [host_feature(1, 10.0, 11.0, 9.0)]

    assert host_line_mask.new_windows([host_feature(2, 10.2, 10.8, 0.0)], detected) == ()


def test_a_window_that_only_touches_a_detected_window_is_not_new() -> None:
    detected = [host_feature(1, 10.0, 11.0, 9.0)]

    assert host_line_mask.new_windows([host_feature(2, 11.0, 12.0, 0.0)], detected) == ()


def test_a_window_clear_of_the_detected_windows_is_new() -> None:
    detected = [host_feature(1, 10.0, 11.0, 9.0)]
    candidate = host_feature(2, 12.0, 13.0, 0.0)

    assert host_line_mask.new_windows([candidate], detected) == (candidate,)
