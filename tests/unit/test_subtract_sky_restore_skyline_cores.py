"""Restoring skyline cores outside the object, and holding the sky spline at its boundary (DY-1972)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate as ip

from soxspipe.commonutils.subtract_sky import (
    SKYLINE_DILATE_SAMPLES,
    _evaluate_sky_spline,
    _skyline_pixel_mask,
    subtract_sky,
)
from tests.unit.test_subtract_sky_bspline_characterization import _subtractor as _fit_subtractor

pytestmark = pytest.mark.unit

RESTORE_KEY = "restore_skyline_cores_outside_object"
# MARKS A SETTINGS FILE THAT DOES NOT CARRY THE KEY AT ALL (NONE AND FALSE ARE VALUES A FILE CAN HOLD)
_MISSING = object()
PIXELS_PER_ORDER = 4000
OBJECT_LOW, OBJECT_HIGH = 1.0, 1.6
# THE DETECTED RANGE COMES FROM A 100-BIN HISTOGRAM (BINS OF 0.11 ARCSEC), SO ITS EDGES SIT A FEW BINS FROM THE STRIPE
RANGE_SLACK = 0.25
RON = 3.3
# THE SKYLINE OF THE CLIP-LEVEL ORDERS: A WIDE LINE WHOSE FLAGGED CORE IS WELL INSIDE THE LINE MASK
LINE_SIGMA_SAMPLES = 8.0
CLIP_LINE_CENTRE = 2000
CLIP_LINE_AMPLITUDE = 300.0
CLIP_LINE_FLAGGED_HALF_WIDTH = 30
CLIP_LINE_CORE_HALF_WIDTH = 10
CLIP_LINE_FAR_SAMPLES = 150

# THE SYNTHETIC SKY ORDER OF THE END-TO-END TESTS
LINE_COLUMNS = 1500
SLIT_ROWS = 30
LINE_COUNT = 12
LINE_AMPLITUDE = 3000.0
LINE_WIDTH_NM = 0.035
CONTINUUM = 80.0
OBJECT_AMPLITUDE = 400.0
OBJECT_CENTRE, OBJECT_WIDTH = 1.3, 0.3
LINE_PEAK_WINDOW_NM = 0.01
# A DARK ORDER ON A FAINT TARGET, WITH A POSITIVE NOISE TAIL OF HOT PIXELS THE ONE-SIDED PRE-FIT CLIP FLAGS. RESTORING
# EVERY CLIPPED PIXEL LIFTED ITS FITTED CONTINUUM BY ABOUT 1.3 E- (26%); SKYLINE-ONLY RESTORE LEAVES IT UNCHANGED
DARK_CONTINUUM = 5.0
DARK_OBJECT_AMPLITUDE = 4.0
HOT_PIXEL_FRACTION = 0.05
HOT_PIXEL_EXCESS = 25.0
DARK_CONTINUUM_TOLERANCE = 0.1
# THE FRACTION OF EVERY LINE PEAK THE FIT MUST RECOVER WITH THE CORES RESTORED. WITHOUT, THE LINES WHOSE CORES THE
# PRE-FIT CLIP FLAGS FALL BELOW UNDER_FIT_PEAK_FRACTION (5 OF THE 12; THE OTHERS KEEP ENOUGH CORE PIXELS)
RECOVERED_PEAK_FRACTION = 0.95
UNDER_FIT_PEAK_FRACTION = 0.8
MIN_UNDER_FIT_LINES = 4


def _clip_subtractor(log: Any, arm: str, restoreSetting: Any = _MISSING) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = arm
    subtractor.ron = RON
    skySettings = {"percentile_rolling_window_size": 15, "noise_rolling_window_size": 45}
    if restoreSetting is not _MISSING:
        skySettings[RESTORE_KEY] = restoreSetting
    subtractor.recipeSettings = {"sky-subtraction": skySettings}
    return subtractor


def _line_profile(index: np.ndarray, centre: int, amplitude: float = CLIP_LINE_AMPLITUDE) -> np.ndarray:
    return amplitude * np.exp(-0.5 * ((index - centre) / LINE_SIGMA_SAMPLES) ** 2)


def _flagged_orders(hasObject: bool) -> list[pd.DataFrame]:
    """Two orders: a flagged skyline core, 3% scattered flags, 2% bad pixels, and (optionally) an object stripe."""
    orders = []
    for seed in (21, 22):
        rng = np.random.default_rng(seed)
        index = np.arange(PIXELS_PER_ORDER)
        slitPosition = rng.uniform(-5.5, 5.5, PIXELS_PER_ORDER)
        scattered = rng.uniform(0, 1, PIXELS_PER_ORDER) < 0.03
        stripe = (
            (slitPosition > OBJECT_LOW) & (slitPosition < OBJECT_HIGH) & (rng.uniform(0, 1, PIXELS_PER_ORDER) < 0.7)
        )
        isLineCore = np.abs(index - CLIP_LINE_CENTRE) <= CLIP_LINE_FLAGGED_HALF_WIDTH
        isBad = rng.uniform(0, 1, PIXELS_PER_ORDER) < 0.02
        orders.append(
            pd.DataFrame(
                {
                    "order": 12,
                    "wavelength": index * 0.01,
                    "flux": 100.0 + rng.normal(0.0, 3.0, PIXELS_PER_ORDER) + _line_profile(index, CLIP_LINE_CENTRE),
                    "error": np.full(PIXELS_PER_ORDER, 3.0),
                    "slit_position": slitPosition,
                    "flagged_edge_clipped": False,
                    "flagged_bad_pixel_clipped": isBad,
                    "flagged_object_clipped": scattered | isLineCore | (stripe if hasObject else False),
                    "flagged_all_clipped": isBad,
                    "flux_minus_smoothed_residual": 3.0 * np.sin(index / 5.0),
                    "flux_percentile_smoothed": 100.0 + np.cos(index / 50.0),
                }
            )
        )
    return orders


def _outside_object(order: pd.DataFrame) -> pd.Series:
    return (order["slit_position"] < OBJECT_LOW - RANGE_SLACK) | (order["slit_position"] > OBJECT_HIGH + RANGE_SLACK)


def _on_line(order: pd.DataFrame) -> pd.Series:
    """The pixels that are certainly on the skyline: its core, well inside the line mask."""
    return pd.Series(np.abs(order.index - CLIP_LINE_CENTRE) <= CLIP_LINE_CORE_HALF_WIDTH, index=order.index)


def _off_line(order: pd.DataFrame) -> pd.Series:
    """The pixels that are certainly not on the skyline: beyond the line mask and its dilation."""
    return pd.Series(np.abs(order.index - CLIP_LINE_CENTRE) > CLIP_LINE_FAR_SAMPLES, index=order.index)


def test_nir_clears_skyline_flags_outside_the_detected_object_range_by_default(log: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=True)
    flaggedOnLine = [int((o["flagged_object_clipped"] & _outside_object(o) & _on_line(o)).sum()) for o in orders]

    # ACT
    result = _clip_subtractor(log, "NIR").clip_object_slit_positions(orders, aggressive=True)

    # ASSERT
    assert min(flaggedOnLine) > 0
    for order in result:
        assert not (order["flagged_object_clipped"] & _outside_object(order) & _on_line(order)).any()
        assert order["flagged_object_clipped"].sum() > 0


def test_flags_away_from_any_skyline_stay_when_cores_are_restored(log: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=True)
    keptBefore = [(o["flagged_object_clipped"] & _outside_object(o) & _off_line(o)).to_numpy().copy() for o in orders]

    # ACT
    result = _clip_subtractor(log, "NIR").clip_object_slit_positions(orders, aggressive=True)

    # ASSERT
    for order, before in zip(result, keptBefore, strict=True):
        assert before.sum() > 0
        assert order.loc[before, "flagged_object_clipped"].all()
        assert order.loc[before, "flagged_all_clipped"].all()


def test_object_stripe_pixels_stay_flagged_when_cores_are_restored(log: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=True)
    stripeBefore = [
        (o["flagged_object_clipped"] & o["slit_position"].between(OBJECT_LOW, OBJECT_HIGH)).to_numpy().copy()
        for o in orders
    ]

    # ACT
    result = _clip_subtractor(log, "NIR").clip_object_slit_positions(orders, aggressive=True)

    # ASSERT
    for order, before in zip(result, stripeBefore, strict=True):
        assert order.loc[before, "flagged_object_clipped"].all()
        assert order.loc[before, "flagged_all_clipped"].all()


def test_restored_pixels_keep_their_bad_pixel_clip(log: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=True)

    # ACT
    result = _clip_subtractor(log, "NIR").clip_object_slit_positions(orders, aggressive=True)

    # ASSERT
    for order in result:
        assert order.loc[order["flagged_bad_pixel_clipped"], "flagged_all_clipped"].all()
        isRestorable = _outside_object(order) & _on_line(order) & ~order["flagged_bad_pixel_clipped"]
        assert isRestorable.sum() > 0
        assert not order.loc[isRestorable, "flagged_all_clipped"].any()


@pytest.mark.parametrize(
    ("arm", "restoreSetting"),
    [("NIR", False), ("VIS", _MISSING), ("UVB", _MISSING), ("VIS", False)],
)
def test_object_flags_outside_the_range_stay_when_restore_is_off(log: Any, arm: str, restoreSetting: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=True)
    flaggedOutside = [int((o["flagged_object_clipped"] & _outside_object(o)).sum()) for o in orders]

    # ACT
    result = _clip_subtractor(log, arm, restoreSetting).clip_object_slit_positions(orders, aggressive=True)

    # ASSERT
    assert [int((o["flagged_object_clipped"] & _outside_object(o)).sum()) for o in result] == flaggedOutside


def test_the_setting_turns_restore_on_for_an_arm_that_defaults_to_off(log: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=True)

    # ACT
    result = _clip_subtractor(log, "VIS", True).clip_object_slit_positions(orders, aggressive=True)

    # ASSERT
    assert not any((o["flagged_object_clipped"] & _outside_object(o) & _on_line(o)).any() for o in result)


def test_flags_are_unchanged_when_no_object_range_is_detected(log: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=False)
    flaggedBefore = [o["flagged_object_clipped"].to_numpy().copy() for o in orders]

    # ACT
    result = _clip_subtractor(log, "NIR").clip_object_slit_positions(orders, aggressive=True)

    # ASSERT
    for order, before in zip(result, flaggedBefore, strict=True):
        assert order["flagged_object_clipped"].to_numpy().tolist() == before.tolist()
        assert order["flagged_object_clipped"].sum() > 0


def test_the_non_aggressive_path_does_not_restore_cores(log: Any) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=True)
    flaggedBefore = [o["flagged_object_clipped"].to_numpy().copy() for o in orders]

    # ACT
    result = _clip_subtractor(log, "NIR").clip_object_slit_positions(orders, aggressive=False)

    # ASSERT
    for order, before in zip(result, flaggedBefore, strict=True):
        assert order["flagged_object_clipped"].to_numpy().tolist() == before.tolist()


@pytest.mark.parametrize("badValue", ["False", "True", 0, 1, None])
def test_a_setting_that_is_not_a_bool_raises_a_value_error_naming_the_key(log: Any, badValue: Any) -> None:
    # ARRANGE
    subtractor = _clip_subtractor(log, "NIR", badValue)

    # ACT / ASSERT
    with pytest.raises(ValueError, match=RESTORE_KEY):
        subtractor._restore_skyline_cores_enabled()


# THE SKYLINE-MASK ORDER: A WIDE LINE ON A FLAT NOISY CONTINUUM, SAMPLED IN WAVELENGTH ORDER
MASK_SAMPLES = 3000
MASK_LINE_CENTRE = 1500
MASK_LINE_AMPLITUDE = 300.0
MASK_NOISE = 3.0
# THE LINE IS CERTAINLY ON THE MASK WITHIN MASK_LINE_CORE_SAMPLES OF ITS CENTRE, AND NOWHERE ON IT BEYOND
# MASK_LINE_WING_SAMPLES (ITS HALF-WIDTH ABOVE THE THRESHOLD IS ABOUT 20 SAMPLES)
MASK_LINE_CORE_SAMPLES = 5
MASK_LINE_WING_SAMPLES = 40


def _mask_inputs(
    withLine: bool, lineCentre: int = MASK_LINE_CENTRE, samples: int = MASK_SAMPLES
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(5)
    index = np.arange(samples)
    flux = 100.0 + rng.normal(0.0, MASK_NOISE, samples)
    if withLine:
        flux = flux + _line_profile(index, lineCentre, MASK_LINE_AMPLITUDE)
    return index * 0.01, flux, np.full(samples, MASK_NOISE)


def test_a_narrow_bright_line_is_marked_including_its_dilated_wings() -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=True)
    wingReach = SKYLINE_DILATE_SAMPLES + MASK_LINE_CORE_SAMPLES
    beyondReach = SKYLINE_DILATE_SAMPLES + MASK_LINE_WING_SAMPLES

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert mask[MASK_LINE_CENTRE - wingReach : MASK_LINE_CENTRE + wingReach + 1].all()
    assert not mask[: MASK_LINE_CENTRE - beyondReach].any()
    assert not mask[MASK_LINE_CENTRE + beyondReach :].any()


def test_flat_noisy_continuum_with_isolated_hot_pixels_is_not_marked() -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=False)
    flux[::97] += 150.0

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert not mask.any()


def test_the_mask_follows_the_row_order_of_the_input_not_the_wavelength_order() -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=True)
    expected = _skyline_pixel_mask(wavelength, flux, error)
    shuffle = np.random.default_rng(3).permutation(MASK_SAMPLES)

    # ACT
    shuffled = _skyline_pixel_mask(wavelength[shuffle], flux[shuffle], error[shuffle])

    # ASSERT
    assert expected.sum() > 0
    assert shuffled.tolist() == expected[shuffle].tolist()


def test_nan_flux_and_error_do_not_create_skyline_pixels() -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=False)
    flux[::53] = np.nan
    error[::71] = np.nan

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert mask.dtype == bool
    assert not mask.any()


def test_an_order_with_no_finite_error_has_no_skyline_pixels() -> None:
    # ARRANGE
    wavelength, flux, _ = _mask_inputs(withLine=True)

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, np.full(MASK_SAMPLES, np.nan))

    # ASSERT
    assert not mask.any()


@pytest.mark.parametrize("samples", [0, 10, 60, 100])
def test_an_order_shorter_than_the_dilation_kernel_gets_a_mask_of_its_own_length(samples: int) -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=True, lineCentre=samples // 2, samples=samples)

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert mask.shape == (samples,)
    assert mask.dtype == bool


def test_a_line_in_an_order_shorter_than_the_dilation_kernel_is_marked() -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=True, lineCentre=50, samples=100)

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert mask[50]


def test_rows_with_a_nan_flux_or_error_inside_a_line_are_not_marked() -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=True)
    nanFlux, nanError = MASK_LINE_CENTRE - 3, MASK_LINE_CENTRE + 4
    flux[nanFlux] = np.nan
    error[nanError] = np.nan

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert not mask[nanFlux]
    assert not mask[nanError]
    assert mask[MASK_LINE_CENTRE]
    assert mask[nanFlux - 1]
    assert mask[nanError + 1]


@pytest.mark.parametrize("badError", [np.inf, -np.inf])
def test_rows_with_an_infinite_error_are_not_marked_and_do_not_set_the_threshold(badError: float) -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=True)
    isInfinite = np.arange(MASK_SAMPLES) % 5 != 0
    error[isInfinite] = badError

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert not mask[isInfinite].any()
    assert mask[MASK_LINE_CENTRE]
    assert mask[~isInfinite].sum() > 0


@pytest.mark.parametrize("lineCentre", [5, MASK_SAMPLES - 6])
def test_a_line_in_the_first_or_last_samples_of_an_order_is_marked_with_its_dilation_truncated(
    lineCentre: int,
) -> None:
    # ARRANGE
    wavelength, flux, error = _mask_inputs(withLine=True, lineCentre=lineCentre)
    reach = SKYLINE_DILATE_SAMPLES + MASK_LINE_CORE_SAMPLES
    beyondReach = SKYLINE_DILATE_SAMPLES + MASK_LINE_WING_SAMPLES

    # ACT
    mask = _skyline_pixel_mask(wavelength, flux, error)

    # ASSERT
    assert mask.shape == (MASK_SAMPLES,)
    assert mask[max(lineCentre - reach, 0) : lineCentre + reach + 1].all()
    assert not mask[: max(lineCentre - beyondReach, 0)].any()
    assert not mask[lineCentre + beyondReach :].any()


# THE DIRECT RESTORE TEST: ONE ORDER WITH A SKYLINE, TWO OBJECT RANGES AND CASES PLACED AT AND AROUND THEIR EDGES
DIRECT_SAMPLES = 4000
DIRECT_LINE_CENTRE = 3000
DIRECT_RANGES = [[-3.0, -2.0], [1.0, 2.0]]
EPS = 1e-6
DEFAULT_SLIT = 4.0
DIRECT_HOT_PIXEL_EXCESS = 150.0
# THE FIRST SAMPLE OF THE SKYLINE CASES, AND OF THE DARK-CONTINUUM CASES (A NOISY STRETCH FAR FROM THE LINE)
LINE_CASE_START = 2994
DARK_CASE_START = 600
DARK_NOISE = 6.0
DARK_STRETCH = (500, 900)

AT_EDGES = [-3.0, -2.0, 1.0, 2.0]
JUST_INSIDE = [-3.0 + EPS, -2.0 - EPS, 1.0 + EPS, 2.0 - EPS]
JUST_OUTSIDE = [-3.0 - EPS, -2.0 + EPS, 1.0 - EPS, 2.0 + EPS]


def _direct_order() -> pd.DataFrame:
    rng = np.random.default_rng(11)
    index = np.arange(DIRECT_SAMPLES)
    noise = np.where((index >= DARK_STRETCH[0]) & (index < DARK_STRETCH[1]), DARK_NOISE, 3.0)
    flux = 100.0 + rng.normal(0.0, 1.0, DIRECT_SAMPLES) * noise + _line_profile(index, DIRECT_LINE_CENTRE)
    return pd.DataFrame(
        {
            "wavelength": index * 0.01,
            "flux": flux,
            "error": np.full(DIRECT_SAMPLES, 3.0),
            "slit_position": np.full(DIRECT_SAMPLES, DEFAULT_SLIT),
            "flagged_edge_clipped": False,
            "flagged_bad_pixel_clipped": False,
            "flagged_object_clipped": True,
            "flagged_all_clipped": True,
        }
    )


def _restored_direct_order(log: Any) -> tuple[pd.DataFrame, dict[str, list[int]]]:
    """Run the restore on the direct order and return it with the sample numbers of every case group."""
    order = _direct_order()
    groups: dict[str, list[int]] = {}
    nextLine = iter(range(LINE_CASE_START, LINE_CASE_START + 60))
    nextDark = iter(range(DARK_CASE_START, DARK_CASE_START + 60))

    def place(name: str, slits: list[float], nextSample: Any) -> None:
        samples = [next(nextSample) for _ in slits]
        order.loc[samples, "slit_position"] = slits
        groups[name] = samples

    place("line_at_edges", AT_EDGES, nextLine)
    place("line_just_inside", JUST_INSIDE, nextLine)
    place("line_just_outside", JUST_OUTSIDE + [-1.0, -5.0, 5.0], nextLine)
    place("line_no_slit", [np.nan, np.nan], nextLine)
    place("line_edge_clipped", [-1.0, -5.0], nextLine)
    place("line_bad_pixel", [-1.0, 5.0], nextLine)
    place("dark_just_outside", JUST_OUTSIDE + [-1.0, DEFAULT_SLIT], nextDark)
    order.loc[groups["line_edge_clipped"], "flagged_edge_clipped"] = True
    order.loc[groups["line_bad_pixel"], "flagged_bad_pixel_clipped"] = True
    order.loc[groups["dark_just_outside"], "flux"] += DIRECT_HOT_PIXEL_EXCESS
    _clip_subtractor(log, "NIR")._restore_skyline_cores_outside_object(order, DIRECT_RANGES)
    return order, groups


@pytest.mark.parametrize("group", ["line_at_edges", "line_just_inside"])
def test_skyline_pixels_at_and_inside_the_range_ends_stay_flagged(log: Any, group: str) -> None:
    # ARRANGE / ACT
    order, groups = _restored_direct_order(log)

    # ASSERT
    assert order.loc[groups[group], "flagged_object_clipped"].all()
    assert order.loc[groups[group], "flagged_all_clipped"].all()


def test_skyline_pixels_just_outside_the_range_ends_are_released(log: Any) -> None:
    # ARRANGE / ACT
    order, groups = _restored_direct_order(log)

    # ASSERT
    assert not order.loc[groups["line_just_outside"], "flagged_object_clipped"].any()
    assert not order.loc[groups["line_just_outside"], "flagged_all_clipped"].any()


def test_hot_dark_continuum_pixels_just_outside_the_range_ends_stay_flagged(log: Any) -> None:
    # ARRANGE / ACT
    order, groups = _restored_direct_order(log)

    # ASSERT
    assert order.loc[groups["dark_just_outside"], "flagged_object_clipped"].all()
    assert order.loc[groups["dark_just_outside"], "flagged_all_clipped"].all()


def test_skyline_pixels_with_no_slit_position_stay_flagged(log: Any) -> None:
    # ARRANGE / ACT
    order, groups = _restored_direct_order(log)

    # ASSERT
    assert order.loc[groups["line_no_slit"], "flagged_object_clipped"].all()
    assert order.loc[groups["line_no_slit"], "flagged_all_clipped"].all()


def test_released_skyline_pixels_keep_their_edge_and_bad_pixel_clips(log: Any) -> None:
    # ARRANGE / ACT
    order, groups = _restored_direct_order(log)

    # ASSERT
    for group in ("line_edge_clipped", "line_bad_pixel"):
        assert not order.loc[groups[group], "flagged_object_clipped"].any()
        assert order.loc[groups[group], "flagged_all_clipped"].all()


@pytest.mark.parametrize(
    ("hasObject", "aggressive"),
    [(False, True), (True, False), (True, True)],
    ids=["no-object-range", "not-aggressive", "restoring"],
)
def test_a_setting_that_is_not_a_bool_raises_on_every_path_through_the_clip(
    log: Any, hasObject: bool, aggressive: bool
) -> None:
    # ARRANGE
    orders = _flagged_orders(hasObject=hasObject)
    subtractor = _clip_subtractor(log, "NIR", "False")

    # ACT / ASSERT
    with pytest.raises(ValueError, match=RESTORE_KEY):
        subtractor.clip_object_slit_positions(orders, aggressive=aggressive)


def _sky_order_with_object(
    seed: int = 1,
    continuum: float = CONTINUUM,
    objectAmplitude: float = OBJECT_AMPLITUDE,
    hotFraction: float = 0.0,
) -> tuple[pd.DataFrame, list[float]]:
    """A tilted NIR-like order: narrow bright skylines on every slit row plus a faint object stripe.

    ``hotFraction`` of the pixels, anywhere on the slit, carry ``HOT_PIXEL_EXCESS`` extra counts (the positive tail).
    """
    rng = np.random.default_rng(seed)
    column = np.repeat(np.arange(LINE_COLUMNS), SLIT_ROWS)
    # CONTINUOUS SLIT POSITIONS, SO THE OBJECT STRIPE FILLS THE SLIT-POSITION HISTOGRAM THE WAY REAL PIXELS DO
    slit = rng.uniform(-5.5, 5.5, LINE_COLUMNS * SLIT_ROWS)
    wavelength = 1500.0 + column * 0.02 + slit * 0.004
    centres = [1500.0 + (k + 0.5) * 0.02 * LINE_COLUMNS / LINE_COUNT for k in range(LINE_COUNT)]
    sky = np.full(wavelength.shape, continuum)
    for centre in centres:
        sky = sky + LINE_AMPLITUDE * np.exp(-0.5 * ((wavelength - centre) / LINE_WIDTH_NM) ** 2)
    objectFlux = objectAmplitude * np.exp(-0.5 * ((slit - OBJECT_CENTRE) / OBJECT_WIDTH) ** 2)
    flux = rng.poisson(sky + objectFlux) + rng.normal(0.0, RON, len(wavelength))
    if hotFraction:
        flux = flux + HOT_PIXEL_EXCESS * (rng.uniform(0, 1, len(wavelength)) < hotFraction)
    order = pd.DataFrame(
        {
            "order": 12,
            "wavelength": wavelength,
            "slit_position": slit,
            "flux": flux.astype(float),
            "error": np.sqrt(sky + objectFlux + RON**2),
            "mask": False,
            "truth": sky,
        }
    )
    return order, centres


def _fit_after_masking(log: Any, order: pd.DataFrame, restoreSetting: Any) -> pd.DataFrame:
    """Fit the sky of one order after the real pre-fit clip and aggressive object masking."""
    subtractor = _fit_subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    skySettings = subtractor.recipeSettings["sky-subtraction"]
    skySettings.update(
        {
            "percentile_clipping_sigma": 3,
            "percentile_clipping_iterations": 15,
            "percentile_rolling_window_size": 35,
            "noise_rolling_window_size": 250,
            "min_points_per_knot": 15,
            "residual_floor_percentile": 80,
            RESTORE_KEY: restoreSetting,
        }
    )
    clipped = subtractor.get_over_sampled_sky_from_order(order)
    masked = subtractor.clip_object_slit_positions([clipped], aggressive=True)[0]
    return subtractor.fit_bspline_curve_to_sky(masked)[0]


def _fitted_line_peaks(log: Any, restoreSetting: Any) -> np.ndarray:
    """The fitted line peak over its true peak, for every line."""
    order, centres = _sky_order_with_object()
    fitted = _fit_after_masking(log, order, restoreSetting)
    ratios = []
    for centre in centres:
        isPeak = (fitted["wavelength"] - centre).abs() < LINE_PEAK_WINDOW_NM
        ratios.append(
            (fitted.loc[isPeak, "sky_model_wl"].mean() - CONTINUUM) / (fitted.loc[isPeak, "truth"].mean() - CONTINUUM)
        )
    return np.array(ratios)


def test_restoring_skyline_cores_lets_the_fit_recover_the_line_peaks(log: Any) -> None:
    # ARRANGE / ACT
    peaks = _fitted_line_peaks(log, restoreSetting=True)

    # ASSERT
    assert peaks.min() > RECOVERED_PEAK_FRACTION
    assert peaks.max() < 2 - RECOVERED_PEAK_FRACTION


def test_without_restoring_skyline_cores_the_fit_under_estimates_the_line_peaks(log: Any) -> None:
    # ARRANGE / ACT
    peaks = _fitted_line_peaks(log, restoreSetting=False)

    # ASSERT
    assert (peaks < UNDER_FIT_PEAK_FRACTION).sum() >= MIN_UNDER_FIT_LINES


def test_restoring_skyline_cores_does_not_lift_the_continuum_of_a_dark_order_with_hot_pixels(log: Any) -> None:
    # ARRANGE
    order, _ = _sky_order_with_object(
        continuum=DARK_CONTINUUM, objectAmplitude=DARK_OBJECT_AMPLITUDE, hotFraction=HOT_PIXEL_FRACTION
    )

    # ACT
    restored = _fit_after_masking(log, order, restoreSetting=True)["sky_model_wl"]
    notRestored = _fit_after_masking(log, order, restoreSetting=False)["sky_model_wl"]

    # ASSERT
    assert abs(restored.median() - notRestored.median()) < DARK_CONTINUUM_TOLERANCE


# THE FIRST OR LAST FRACTION OF THE ORDER'S PIXELS IS CLIPPED, SO THE SPLINE ENDS INSIDE THE ORDER
CLIPPED_END_FRACTION = 0.05
END_PIXEL_COUNT = 20000


def _order_with_clipped_end(side: str) -> pd.DataFrame:
    rng = np.random.default_rng(7)
    wavelength = np.sort(rng.uniform(500.0, 510.0, END_PIXEL_COUNT))
    # A CURVED SKY SO THE END CUBIC PIECE CARRIES A SLOPE THAT EXTRAPOLATION WOULD FOLLOW
    sky = CONTINUUM + 40.0 * ((wavelength - 500.0) / 10.0) ** 2 + 30.0 * np.sin((wavelength - 500.0) / 1.7)
    noise = np.sqrt(sky + RON**2)
    flux = rng.poisson(sky) + rng.normal(0.0, RON, END_PIXEL_COUNT)
    clippedCount = int(END_PIXEL_COUNT * CLIPPED_END_FRACTION)
    isClipped = (
        np.arange(END_PIXEL_COUNT) < clippedCount
        if side == "head"
        else np.arange(END_PIXEL_COUNT) >= END_PIXEL_COUNT - clippedCount
    )
    return pd.DataFrame(
        {
            "order": 12,
            "wavelength": wavelength,
            "flux": flux.astype(float),
            "error": noise,
            "residual_windowed_std": noise.copy(),
            "flux_percentile_smoothed": sky,
            "residual_windowed_long_median": np.zeros(END_PIXEL_COUNT),
            "flux_windowed_long_median": np.full(END_PIXEL_COUNT, CONTINUUM),
            "flagged_all_clipped": isClipped,
            "flagged_bspline_clipped": False,
        }
    )


def _fit_with_clipped_end(log: Any, side: str) -> pd.DataFrame:
    order = _order_with_clipped_end(side)
    subtractor = _fit_subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update({"min_points_per_knot": 155, "residual_floor_percentile": 90})
    return subtractor.fit_bspline_curve_to_sky(order)[0]


def _beyond_fitted_range(fitted: pd.DataFrame, side: str) -> tuple[np.ndarray, int]:
    """The positions of the pixels beyond the fitted samples at one end, and the position of the nearest fitted one."""
    kept = np.flatnonzero(~fitted["flagged_all_clipped"].to_numpy(dtype=bool))
    if side == "head":
        return np.arange(0, kept[0]), int(kept[0])
    return np.arange(kept[-1] + 1, len(fitted)), int(kept[-1])


@pytest.mark.parametrize("side", ["head", "tail"])
def test_pixels_beyond_the_fitted_samples_get_the_boundary_model_value(log: Any, side: str) -> None:
    # ARRANGE / ACT
    fitted = _fit_with_clipped_end(log, side)

    # ASSERT
    beyond, nearest = _beyond_fitted_range(fitted, side)
    model = fitted["sky_model_wl"].to_numpy()
    assert len(beyond) > 0
    assert np.all(np.isfinite(model[beyond]))
    assert np.allclose(model[beyond], model[beyond[0]], rtol=0, atol=1e-9)
    assert abs(model[beyond[0]] - model[nearest]) < 0.05 * model[beyond[0]]


@pytest.mark.parametrize("side", ["head", "tail"])
def test_the_model_derivative_is_zero_beyond_the_fitted_samples(log: Any, side: str) -> None:
    # ARRANGE / ACT
    fitted = _fit_with_clipped_end(log, side)

    # ASSERT
    beyond, _ = _beyond_fitted_range(fitted, side)
    assert len(beyond) > 0
    assert np.all(fitted["sky_model_wl_derivative"].to_numpy()[beyond] == 0.0)


def _spline_over_one_to_ten() -> tuple[tuple, float]:
    x = np.linspace(1.0, 10.0, 40)
    tck = ip.splrep(x, 5.0 + x**2, k=3)
    return tck, float(ip.splev(1.0, tck))


def test_the_spline_is_constant_below_its_first_fitted_sample() -> None:
    # ARRANGE
    tck, lowerValue = _spline_over_one_to_ten()
    below = np.array([0.99, 0.5, -3.0])

    # ACT
    values = _evaluate_sky_spline(below, tck)
    slopes = _evaluate_sky_spline(below, tck, derivative=True)

    # ASSERT
    assert np.allclose(values, lowerValue, rtol=0, atol=1e-9)
    assert np.all(slopes == 0.0)


def test_the_spline_slope_inside_the_fitted_range_is_the_true_derivative() -> None:
    # ARRANGE
    tck, _ = _spline_over_one_to_ten()

    # ACT
    slope = _evaluate_sky_spline(np.array([5.0]), tck, derivative=True)

    # ASSERT
    assert slope[0] == pytest.approx(10.0, rel=1e-3)
