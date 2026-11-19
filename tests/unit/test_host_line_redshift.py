"""Inferring the host redshift from the host-line windows found in a stare frame, and the companion windows (DY-2242).

Two or more windows that agree on one redshift for distinct lines of a short rest-frame list give the redshift. The
other lines of the list that fall in an order are then masked too, because a weak line such as [NII] 6548 is below
the detection limit on its own.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.host_line_mask import (
    companion_windows,
    infer_host_redshift,
    load_host_lines,
)
from soxspipe.commonutils.nod_host_contamination import host_feature

pytestmark = pytest.mark.unit

Z = 0.0166
HALPHA_NM = 656.28
NII_6583_NM = 658.35
NII_6548_NM = 654.80
HALF_WIDTH_NM = 0.25
STEP_NM = 0.04


def _window(restNm: float, z: float = Z, order: int = 1, halfWidthNm: float = HALF_WIDTH_NM) -> host_feature:
    centre = restNm * (1 + z)
    return host_feature(order, centre - halfWidthNm, centre + halfWidthNm, 9.0)


def _order(firstNm: float, lastNm: float) -> pd.DataFrame:
    """The wavelengths of one synthetic order, 20 slit pixels wide, with the detector columns the step is read from."""
    columns = np.arange(int((lastNm - firstNm) / STEP_NM))
    return pd.DataFrame(
        {
            "x": np.repeat(np.arange(20), len(columns)),
            "y": np.tile(columns, 20),
            "order": 1.0,
            "wavelength": np.tile(firstNm + STEP_NM * columns, 20),
        }
    )


def test_the_rest_line_list_holds_the_strong_host_lines_in_air_wavelengths() -> None:
    lines = dict(load_host_lines())

    assert lines["Halpha"] == pytest.approx(HALPHA_NM)
    assert {"[NII] 6548", "[NII] 6583", "[SIII] 9069", "[SIII] 9531", "Paalpha"} <= set(lines)
    assert all(wavelength > 0 for wavelength in lines.values())


def test_two_windows_that_agree_on_a_redshift_give_it() -> None:
    redshift = infer_host_redshift([_window(HALPHA_NM), _window(NII_6583_NM)], load_host_lines())

    assert redshift is not None
    assert redshift.z == pytest.approx(Z, abs=1e-5)
    assert redshift.lineNames == ("Halpha", "[NII] 6583")


def test_one_window_gives_no_redshift() -> None:
    assert infer_host_redshift([_window(HALPHA_NM)], load_host_lines()) is None


def test_no_window_gives_no_redshift() -> None:
    assert infer_host_redshift([], load_host_lines()) is None


def test_two_windows_that_do_not_agree_on_a_redshift_give_none() -> None:
    windows = [host_feature(1, 600.0, 600.5, 9.0), host_feature(2, 701.3, 701.8, 9.0)]

    assert infer_host_redshift(windows, load_host_lines()) is None


def test_two_windows_for_one_redshift_that_differ_by_more_than_the_tolerance_give_none() -> None:
    windows = [_window(HALPHA_NM, z=0.0166), _window(953.1, z=0.0266)]

    assert infer_host_redshift(windows, load_host_lines()) is None


def test_two_windows_on_the_same_line_are_not_two_lines() -> None:
    windows = [_window(HALPHA_NM, order=1), _window(HALPHA_NM, order=2)]

    assert infer_host_redshift(windows, load_host_lines()) is None


def test_the_redshift_with_the_most_windows_wins() -> None:
    windows = [_window(HALPHA_NM), _window(NII_6583_NM), _window(953.1), host_feature(1, 860.0, 860.4, 9.0)]

    redshift = infer_host_redshift(windows, load_host_lines())

    assert redshift is not None
    assert redshift.lineNames == ("Halpha", "[NII] 6583", "[SIII] 9531")


def test_the_half_width_fraction_is_the_median_of_the_detected_windows() -> None:
    windows = [
        _window(HALPHA_NM, halfWidthNm=0.2),
        _window(NII_6583_NM, halfWidthNm=0.4),
        _window(953.1, halfWidthNm=0.3),
    ]

    redshift = infer_host_redshift(windows, load_host_lines())

    assert redshift is not None
    centre = 953.1 * (1 + Z)
    assert redshift.halfWidthFraction == pytest.approx(0.3 / centre, rel=1e-3)


def _companions(order: pd.DataFrame, detected: list[host_feature], skylines=(), **kwargs) -> tuple[host_feature, ...]:
    redshift = infer_host_redshift(detected, load_host_lines())
    assert redshift is not None
    return companion_windows(order, detected, redshift, load_host_lines(), skylines, **kwargs)


def test_a_companion_line_in_the_order_gets_a_window_of_the_width_of_the_detected_ones() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]

    (companion,) = _companions(_order(660.0, 668.0), detected)

    centre = NII_6548_NM * (1 + Z)
    assert companion.wavelengthMinNm < centre < companion.wavelengthMaxNm
    expected = redshift_half_width(detected, centre)
    assert (companion.wavelengthMaxNm - companion.wavelengthMinNm) / 2 == pytest.approx(expected, rel=1e-3)
    assert companion.order == 1


def redshift_half_width(detected: list[host_feature], centreNm: float) -> float:
    redshift = infer_host_redshift(detected, load_host_lines())
    assert redshift is not None
    return redshift.halfWidthFraction * centreNm


def test_a_companion_line_that_was_detected_gets_no_second_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM), _window(NII_6548_NM)]

    assert _companions(_order(660.0, 672.0), detected) == ()


def test_a_companion_line_outside_the_order_gets_no_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]

    assert _companions(_order(669.0, 680.0), detected) == ()


def test_a_companion_line_at_the_edge_of_the_order_gets_no_window() -> None:
    detected = [_window(HALPHA_NM), _window(NII_6583_NM)]
    centre = NII_6548_NM * (1 + Z)

    assert _companions(_order(centre - 0.1, centre + 5.0), detected) == ()


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

    names = sorted(round((window.wavelengthMinNm + window.wavelengthMaxNm) / 2, 1) for window in windows)
    expected = sorted(round(rest * (1 + Z), 1) for rest in (NII_6548_NM, 671.64, 673.08))
    assert names == expected
