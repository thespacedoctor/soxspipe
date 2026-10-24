"""Faint-order detection, anchoring and the binned trace fit (DY-1941)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.faint_order_profile import collapse_bin_profile, fit_fixed_width_gaussian
from soxspipe.commonutils.faint_order_trace import (
    VIS_BAND_GROUPS,
    anchor_neighbours,
    bright_edge_anchor,
    faint_trace_settings,
    find_faint_orders,
    make_pixels_at,
    pixels_per_arcsec,
)
from soxspipe.commonutils.toolkit import SOXS_VIS_ORDERS_BLUE_TO_RED
from tests.unit._continuum_vis_scene import vis_dispersion_map, vis_slit_scale

pytestmark = pytest.mark.unit

SETTINGS = faint_trace_settings.from_block({"faint-order-trace": True})


def _order_rows(
    order: float,
    count: int = 100,
    detectedFraction: float = 1.0,
    scatterPx: float = 0.0,
    curvaturePx: float = 0.0,
    seed: int = 1941,
) -> pd.DataFrame:
    """Return the sampled rows of one order: `detectedFraction` of them found, with Gaussian centre-shift noise."""
    rng = np.random.default_rng(seed + int(order))
    wavelength = np.linspace(400.0, 500.0, count)
    unit = (wavelength - 450.0) / 50.0
    shift = (
        2.0 + curvaturePx * unit**2 + rng.normal(0.0, scatterPx, count) if scatterPx else 2.0 + curvaturePx * unit**2
    )
    isDetected = np.arange(count) < round(detectedFraction * count)
    fitX = np.full(count, 100.0 * order)
    return pd.DataFrame(
        {
            "order": order,
            "wavelength": wavelength,
            "fit_x": fitX,
            "fit_y": np.linspace(10.0, 10.0 + count - 1, count),
            "cont_x": np.where(isDetected, fitX - shift, np.nan),
            "cont_y": np.where(isDetected, np.linspace(10.0, 10.0 + count - 1, count), np.nan),
            "centre_shift": np.where(isDetected, shift, np.nan),
            "gauss_stddev": np.where(isDetected, 1.0 + 2.0 * (wavelength - 400.0) / 100.0, np.nan),
            "pre-clipped": ~isDetected,
        }
    )


def _table(*orderRows: pd.DataFrame) -> pd.DataFrame:
    return pd.concat(orderRows, ignore_index=True)


def test_find_faint_orders_flags_an_order_below_the_detection_fraction() -> None:
    # ARRANGE: ORDER 2 IS DETECTED IN 40% OF ITS SAMPLES, ORDER 3 IN ALL
    table = _table(_order_rows(2.0, detectedFraction=0.4), _order_rows(3.0))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert list(faint) == [2.0]
    assert faint[2.0].detectionFraction == pytest.approx(0.4)


def test_find_faint_orders_flags_an_order_whose_centre_shift_scatters_more_than_the_limit() -> None:
    # ARRANGE: ORDER 2 IS DETECTED EVERYWHERE, BUT ITS CENTRES SCATTER BY 3 PX
    table = _table(_order_rows(2.0, scatterPx=3.0), _order_rows(3.0, scatterPx=0.2))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert list(faint) == [2.0]
    assert faint[2.0].scatterPx == pytest.approx(3.0, rel=0.3)


def test_find_faint_orders_leaves_a_curved_but_clean_trace_alone() -> None:
    # ARRANGE: A 4 PX PARABOLA IS THE FIT'S OWN SHAPE, NOT SCATTER
    table = _table(_order_rows(3.0, curvaturePx=4.0, scatterPx=0.2))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert faint == {}


def test_find_faint_orders_flags_an_order_with_no_detection() -> None:
    # ARRANGE
    table = _table(_order_rows(2.0, detectedFraction=0.0), _order_rows(3.0))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert faint[2.0].detectionFraction == 0.0
    assert np.isnan(faint[2.0].scatterPx)


def test_find_faint_orders_flags_an_order_with_too_few_accepted_samples_for_a_fit() -> None:
    # ARRANGE: THREE ACCEPTED SAMPLES IN 6 ARE DETECTED, AND A DEGREE-2 FIT NEEDS 4
    table = _table(_order_rows(2.0, count=6, detectedFraction=0.5), _order_rows(3.0))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert 2.0 in faint


def test_find_faint_orders_ignores_pre_clipped_samples_when_measuring_scatter() -> None:
    # ARRANGE: A CLEAN ORDER WHOSE FIVE OUTLIERS WERE PRE-CLIPPED BY THE SAMPLER
    rows = _order_rows(3.0, scatterPx=0.2)
    rows.loc[:4, "centre_shift"] = [30.0, -30.0, 25.0, -25.0, 20.0]
    rows.loc[:4, "pre-clipped"] = True

    # ACT
    faint = find_faint_orders(rows, "x", SETTINGS)

    # ASSERT
    assert faint == {}


def test_find_faint_orders_uses_the_limits_in_the_settings() -> None:
    # ARRANGE: 70% DETECTED IS FAINT UNDER A 0.8 LIMIT
    table = _table(_order_rows(2.0, detectedFraction=0.7), _order_rows(3.0))
    strict = faint_trace_settings.from_block({"faint-order-trace": True, "faint-order-min-detection-fraction": 0.8})

    # ACT
    faintDefault = find_faint_orders(table, "x", SETTINGS)
    faintStrict = find_faint_orders(table, "x", strict)

    # ASSERT
    assert faintDefault == {}
    assert list(faintStrict) == [2.0]


def _orders_table(*orders: float) -> pd.DataFrame:
    return pd.DataFrame({"order": [float(order) for order in orders for _ in range(3)]})


def test_band_groups_are_the_adjacent_pairs_of_the_blue_to_red_order_list() -> None:
    # ARRANGE
    blueToRed = SOXS_VIS_ORDERS_BLUE_TO_RED

    # ACT
    pairs = [tuple(blueToRed[start : start + 2]) for start in (0, 2)]

    # ASSERT
    assert tuple(pairs) == VIS_BAND_GROUPS


@pytest.mark.parametrize(
    ("faintOrder", "expected"),
    [(2.0, 3.0), (4.0, 1.0), (3.0, 2.0), (1.0, 4.0)],
    ids=["u-from-g", "i-from-r", "g-from-u", "r-from-i"],
)
def test_anchor_neighbours_takes_the_other_order_of_the_band_group(faintOrder: float, expected: float) -> None:
    # ARRANGE
    table = _orders_table(1, 2, 3, 4)

    # ACT
    anchors = anchor_neighbours(table, {faintOrder: None})

    # ASSERT
    assert anchors == {faintOrder: expected}


def test_anchor_neighbours_gives_no_anchor_when_both_orders_of_a_group_are_faint() -> None:
    # ARRANGE: U AND G ARE BOTH FAINT, AND R IS BRIGHT BUT ACROSS THE BAND-GROUP BOUNDARY FROM G
    table = _orders_table(1, 2, 3, 4)

    # ACT
    anchors = anchor_neighbours(table, {2.0: None, 3.0: None, 4.0: None})

    # ASSERT: ONLY I GETS AN ANCHOR
    assert anchors == {4.0: 1.0}


def test_anchor_neighbours_skips_an_order_whose_partner_is_not_in_the_table() -> None:
    # ARRANGE
    table = _orders_table(2, 1, 4)

    # ACT
    anchors = anchor_neighbours(table, {2.0: None})

    # ASSERT
    assert anchors == {}


def test_anchor_neighbours_ignores_orders_outside_the_vis_band_groups() -> None:
    # ARRANGE
    table = _orders_table(11, 12)

    # ACT
    anchors = anchor_neighbours(table, {11.0: None})

    # ASSERT
    assert anchors == {}


def _stub_pixels_at(pixelsPerArcsec: float):
    """Return a `pixelsAt` for a map with order `n` at 100 n px and a slit of `pixelsPerArcsec` px per arcsec."""

    def pixels_at(order: float, wavelengths: np.ndarray, slitArcsec: float | np.ndarray) -> tuple:
        wavelengths = np.asarray(wavelengths, dtype=float)
        slit = np.broadcast_to(np.asarray(slitArcsec, dtype=float), wavelengths.shape)
        return 100.0 * order + pixelsPerArcsec * slit, wavelengths - 390.0

    return pixels_at


@pytest.mark.parametrize(("pixelsPerArcsec", "expectedSlit"), [(-4.0, 0.5), (4.0, -0.5)], ids=["flipped", "unflipped"])
def test_bright_edge_anchor_converts_the_centre_shift_to_slit_arcsec_with_the_sign_of_the_map(
    pixelsPerArcsec: float, expectedSlit: float
) -> None:
    # ARRANGE: THE TRACE SITS 2 PX FROM THE MAP'S SLIT-ZERO TRACK AT EVERY WAVELENGTH
    rows = _order_rows(3.0)

    # ACT
    slitArcsec, _ = bright_edge_anchor(rows, 400.0, _stub_pixels_at(pixelsPerArcsec))

    # ASSERT: CENTRE = FIT - SHIFT, SO THE TRACE IS AT SLIT = -SHIFT / SCALE
    assert slitArcsec == pytest.approx(expectedSlit)


def test_bright_edge_anchor_evaluates_the_smooth_centre_shift_at_the_requested_edge() -> None:
    # ARRANGE: THE SHIFT IS 2 + 4 u^2 PX WITH u RUNNING -1 TO +1 ACROSS THE ORDER, SO IT IS 6 PX AT EACH END
    rows = _order_rows(3.0, curvaturePx=4.0, scatterPx=0.05)

    # ACT
    blueSlit, _ = bright_edge_anchor(rows, 400.0, _stub_pixels_at(-4.0))
    middleSlit, _ = bright_edge_anchor(rows, 450.0, _stub_pixels_at(-4.0))

    # ASSERT
    assert blueSlit == pytest.approx(1.5, abs=0.05)
    assert middleSlit == pytest.approx(0.5, abs=0.05)


def test_bright_edge_anchor_is_not_pulled_by_clipped_outliers() -> None:
    # ARRANGE
    rows = _order_rows(3.0, scatterPx=0.1)
    rows.loc[::10, "centre_shift"] = 40.0

    # ACT
    slitArcsec, _ = bright_edge_anchor(rows, 400.0, _stub_pixels_at(-4.0))

    # ASSERT
    assert slitArcsec == pytest.approx(0.5, abs=0.05)


def test_bright_edge_anchor_takes_the_width_from_the_samples_nearest_the_edge_in_arcsec() -> None:
    # ARRANGE: THE STDDEV GROWS FROM 1 PX AT 400 NM TO 3 PX AT 500 NM
    rows = _order_rows(3.0)

    # ACT
    _, blueSigma = bright_edge_anchor(rows, 400.0, _stub_pixels_at(-4.0))
    _, redSigma = bright_edge_anchor(rows, 500.0, _stub_pixels_at(-4.0))

    # ASSERT: THE MEDIAN OF THE NEAREST QUARTER OF THE ORDER, DIVIDED BY 4 PX/ARCSEC
    assert blueSigma == pytest.approx(1.25 / 4.0, abs=0.02)
    assert redSigma == pytest.approx(2.75 / 4.0, abs=0.02)


def test_make_pixels_at_reads_the_dispersion_map_in_frame_pixels(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE FRAME IS BINNED 2 IN y RELATIVE TO THE MAP, AND THE SLIT RUNS ALONG x
    mapPath = vis_dispersion_map(tmp_path)
    pixelsAt = make_pixels_at(log, mapPath, {"x": 1.0, "y": 2.0}, "x", "y")

    # ACT
    axisA, axisB = pixelsAt(3.0, np.array([410.0, 430.0]), 0.5)

    # ASSERT: x = -30 + 90 * 3 + 0.5 * SCALE(3) AND y = (WAVELENGTH - 390) / 2
    np.testing.assert_allclose(axisA, 240.0 + 0.5 * vis_slit_scale(3.0), atol=1e-3)
    np.testing.assert_allclose(axisB, [10.0, 20.0], atol=1e-3)


def test_make_pixels_at_keeps_samples_that_fall_off_the_detector(log: object, tmp_path: Path) -> None:
    # ARRANGE: A SLIT POSITION 100 ARCSEC OUT PUTS x BELOW ZERO
    pixelsAt = make_pixels_at(log, vis_dispersion_map(tmp_path), {"x": 1.0, "y": 1.0}, "x", "y")

    # ACT
    axisA, axisB = pixelsAt(1.0, np.array([410.0, 420.0]), 100.0)

    # ASSERT
    assert len(axisA) == len(axisB) == 2
    assert (axisA < 0).all()


def test_make_pixels_at_follows_the_axis_roles_it_is_given(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE SAME MAP WITH THE ROLES OF THE AXES SWAPPED
    pixelsAt = make_pixels_at(log, vis_dispersion_map(tmp_path), {"x": 1.0, "y": 1.0}, "y", "x")

    # ACT
    axisA, axisB = pixelsAt(2.0, np.array([410.0]), 0.0)

    # ASSERT
    np.testing.assert_allclose(axisA, [20.0], atol=1e-3)
    np.testing.assert_allclose(axisB, [150.0], atol=1e-3)


def test_pixels_per_arcsec_is_the_signed_slit_scale_of_the_map(log: object, tmp_path: Path) -> None:
    # ARRANGE
    pixelsAt = make_pixels_at(log, vis_dispersion_map(tmp_path), {"x": 1.0, "y": 1.0}, "x", "y")

    # ACT
    scale = pixels_per_arcsec(pixelsAt, 4.0, np.array([400.0, 420.0]))

    # ASSERT
    np.testing.assert_allclose(scale, [vis_slit_scale(4.0)] * 2, atol=1e-3)


def _gaussian_image(centre: float, sigma: float = 1.5, shape: tuple = (60, 80), amplitude: float = 50.0) -> np.ndarray:
    """Return an image with a vertical stripe: a Gaussian across x centred at `centre` in every row."""
    columns = np.arange(shape[1], dtype=float)
    return np.tile(amplitude * np.exp(-0.5 * ((columns - centre) / sigma) ** 2), (shape[0], 1))


def _centroid(offsets: np.ndarray, profile: np.ndarray) -> float:
    weights = np.clip(profile, 0.0, None)
    return float((offsets * weights).sum() / weights.sum())


def test_collapse_bin_profile_centres_the_profile_on_the_track_with_sub_pixel_accuracy() -> None:
    # ARRANGE: THE STRIPE SITS 1.3 PX ABOVE A TRACK THAT RUNS ALONG x = 40.3 FOR 20 ROWS
    image = _gaussian_image(41.6)
    trackA, trackB = np.full(20, 40.3), np.arange(10.0, 30.0)

    # ACT
    offsets, profile = collapse_bin_profile(image, trackA, trackB, 8, "x")

    # ASSERT
    np.testing.assert_array_equal(offsets, np.arange(-8.0, 9.0))
    assert _centroid(offsets, profile) == pytest.approx(1.3, abs=0.05)


def test_collapse_bin_profile_cuts_along_the_axis_it_is_told() -> None:
    # ARRANGE: THE SAME STRIPE ON A TRANSPOSED IMAGE, WITH THE PROFILE CUT ALONG y
    trackA, trackB = np.full(20, 40.3), np.arange(10.0, 30.0)
    _, expected = collapse_bin_profile(_gaussian_image(41.6), trackA, trackB, 8, "x")

    # ACT
    _, profile = collapse_bin_profile(_gaussian_image(41.6).T.copy(), trackA, trackB, 8, "y")

    # ASSERT
    np.testing.assert_allclose(profile, expected)


def test_collapse_bin_profile_ignores_a_cosmic_ray_row() -> None:
    # ARRANGE: ONE ROW HOLDS A 10000-COUNT HIT NEXT TO THE TRACE
    image = _gaussian_image(40.0)
    image[15, 38:42] += 10000.0
    trackA, trackB = np.full(20, 40.0), np.arange(10.0, 30.0)

    # ACT
    offsets, profile = collapse_bin_profile(image, trackA, trackB, 8, "x")

    # ASSERT
    assert profile.max() == pytest.approx(50.0, rel=0.02)
    assert _centroid(offsets, profile) == pytest.approx(0.0, abs=0.02)


def test_collapse_bin_profile_skips_masked_pixels_and_samples_off_the_image() -> None:
    # ARRANGE: NaN MARKS MASKED PIXELS, AND THE LAST TRACK POINTS LIE ABOVE THE IMAGE
    image = _gaussian_image(40.0)
    image[12:14, :] = np.nan
    trackA, trackB = np.full(24, 40.0), np.arange(10.0, 34.0) + np.r_[np.zeros(20), np.full(4, 100.0)]

    # ACT
    offsets, profile = collapse_bin_profile(image, trackA, trackB, 8, "x")

    # ASSERT
    assert np.isfinite(profile).all()
    assert _centroid(offsets, profile) == pytest.approx(0.0, abs=0.02)


def test_collapse_bin_profile_is_nan_where_no_row_has_a_value() -> None:
    # ARRANGE: THE TRACK RUNS ALONG THE LEFT EDGE, SO THE OFFSETS BELOW -2 LIE OFF THE IMAGE
    image = _gaussian_image(2.0)
    trackA, trackB = np.full(10, 2.0), np.arange(10.0, 20.0)

    # ACT
    _, profile = collapse_bin_profile(image, trackA, trackB, 6, "x")

    # ASSERT
    assert np.isnan(profile[:4]).all()
    assert np.isfinite(profile[4:]).all()


PROFILE_OFFSETS = np.arange(-16.0, 17.0)
PROFILE_SIGMA = 1.5
PROFILE_WINDOW = 6.5


def _profile(
    sources: tuple = ((0.0, 50.0),), noise: float = 0.0, seed: int = 1941, baseline: float = 0.0
) -> np.ndarray:
    """Return a collapsed profile: Gaussian sources given as (centre, amplitude) pairs, plus Gaussian noise."""
    values = np.full(PROFILE_OFFSETS.shape, baseline)
    for centre, amplitude in sources:
        values = values + amplitude * np.exp(-0.5 * ((PROFILE_OFFSETS - centre) / PROFILE_SIGMA) ** 2)
    return values + np.random.default_rng(seed).normal(0.0, noise, values.shape) if noise else values


def _fit(profile: np.ndarray, minSnr: float = 3.0) -> float | None:
    return fit_fixed_width_gaussian(profile, PROFILE_OFFSETS, PROFILE_SIGMA, PROFILE_WINDOW, minSnr)


@pytest.mark.parametrize("truth", [-3.2, 0.0, 2.7])
def test_fit_fixed_width_gaussian_recovers_the_centre_of_a_clean_profile(truth: float) -> None:
    # ARRANGE
    profile = _profile(((truth, 50.0),), baseline=3.0)

    # ACT
    centre = _fit(profile)

    # ASSERT
    assert centre == pytest.approx(truth, abs=0.05)


def test_fit_fixed_width_gaussian_is_unbiased_for_a_peak_six_times_the_noise() -> None:
    # ARRANGE: THE LIMIT ON ONE PROFILE'S CENTRE IS 0.22 PX AT A PEAK SNR OF 6, SO TEST THE MEAN OF FORTY
    centres = [_fit(_profile(((1.5, 12.0),), noise=2.0, seed=seed)) for seed in range(40)]

    # ACT
    errors = np.array([centre - 1.5 for centre in centres if centre is not None])

    # ASSERT: ALMOST NO PROFILE IS REJECTED, AND THEIR MEAN LANDS ON THE TRUTH
    assert len(errors) >= 38
    assert abs(errors.mean()) < 0.1
    assert np.median(np.abs(errors)) < 0.25


def test_fit_fixed_width_gaussian_returns_none_for_pure_noise() -> None:
    # ACT
    centres = [_fit(_profile(((0.0, 0.0),), noise=1.0, seed=seed)) for seed in range(40)]

    # ASSERT
    assert centres == [None] * 40


def test_fit_fixed_width_gaussian_returns_none_below_the_snr_limit_and_keeps_it_above() -> None:
    # ARRANGE: A PEAK OF 4 OVER NOISE OF 2 HAS A PEAK SNR OF ABOUT 2
    profile = _profile(((0.0, 4.0),), noise=2.0, seed=3)

    # ACT
    strict = _fit(profile, minSnr=3.0)
    lenient = _fit(profile, minSnr=1.0)

    # ASSERT
    assert strict is None
    assert lenient is not None


def test_fit_fixed_width_gaussian_does_not_land_on_a_brighter_source_beyond_the_window() -> None:
    # ARRANGE: THE TARGET IS AT +1 PX, AND A BRIGHTER SOURCE SITS 3.7 ARCSEC (13.3 PX) AWAY
    profile = _profile(((1.0, 5.0), (13.3, 12.0)), noise=0.5)

    # ACT
    centre = _fit(profile)

    # ASSERT
    assert centre == pytest.approx(1.0, abs=0.3)


def test_fit_fixed_width_gaussian_returns_none_when_only_a_source_beyond_the_window_is_present() -> None:
    # ARRANGE
    profile = _profile(((13.3, 12.0),), noise=0.5)

    # ACT
    centre = _fit(profile)

    # ASSERT: THE BEST CENTRE INSIDE THE WINDOW IS ITS EDGE, WHICH IS NOT A TRACE
    assert centre is None


def test_fit_fixed_width_gaussian_never_picks_a_negative_lobe() -> None:
    # ARRANGE: THE PARTNER NOD'S NEGATIVE TRACE FALLS INSIDE THE WINDOW AT -4 PX, WITH NO POSITIVE TRACE AT ALL
    profile = _profile(((-4.0, -30.0),), noise=0.5)

    # ACT
    centre = _fit(profile)

    # ASSERT
    assert centre is None or abs(centre + 4.0) > 2.0


def test_fit_fixed_width_gaussian_keeps_the_positive_trace_beside_a_negative_lobe() -> None:
    # ARRANGE: THE NEGATIVE LOBE SITS AT -9 PX, AT THE EDGE OF THE FITTED RANGE
    profile = _profile(((1.0, 20.0), (-9.0, -30.0)), noise=0.5)

    # ACT
    centre = _fit(profile)

    # ASSERT
    assert centre == pytest.approx(1.0, abs=0.3)


def test_fit_fixed_width_gaussian_skips_nan_samples() -> None:
    # ARRANGE
    profile = _profile(((2.0, 40.0),), noise=0.5)
    profile[::5] = np.nan

    # ACT
    centre = _fit(profile)

    # ASSERT
    assert centre == pytest.approx(2.0, abs=0.1)


def test_fit_fixed_width_gaussian_returns_none_when_too_few_samples_are_valid() -> None:
    # ARRANGE
    profile = _profile(((0.0, 50.0),))
    profile[3:] = np.nan

    # ACT
    centre = _fit(profile)

    # ASSERT
    assert centre is None


def test_fit_fixed_width_gaussian_is_deterministic() -> None:
    # ARRANGE
    profile = _profile(((1.7, 9.0),), noise=2.0, seed=5)

    # ACT
    first, second = _fit(profile), _fit(profile.copy())

    # ASSERT
    assert first == second
