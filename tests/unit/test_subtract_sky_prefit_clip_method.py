"""The `prefit_clip_method` setting: the offset-bin pre-filter beside the rolling-percentile pre-fit clip (DY-2055)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from tests.unit.test_subtract_sky_bspline_characterization import _subtractor

pytestmark = pytest.mark.unit

METHOD_KEY = "prefit_clip_method"
COLUMNS = 400
SLIT_ROWS = 30
# NIR-LIKE: 0.02 NM PER COLUMN, SO THE DEFAULT 5-PIXEL BIN IS 0.1 NM WIDE, AND A SLIT TILT OF 0.004 NM PER SLIT UNIT
# THAT INTERLEAVES THE WAVELENGTHS OF NEIGHBOURING COLUMNS
NM_PER_COLUMN = 0.02
SLIT_TILT_NM = 0.004
SKY_LEVEL = 80.0
RON = 3.3
LINE_AMPLITUDE = 3000.0
# THREE LINES OF FWHM ABOUT 4 COLUMNS, AT COLUMNS THAT FALL ON, NEAR AND AWAY FROM A 5-COLUMN BIN EDGE
LINE_COLUMNS = (100.0, 202.5, 305.0)
LINE_SIGMA_COLUMNS = 1.75
OBJECT_AMPLITUDE = 300.0
OBJECT_CENTRE, OBJECT_WIDTH = 1.0, 0.4
COSMIC_AMPLITUDE = 4000.0
COSMIC_ROW = 250 * SLIT_ROWS + 3
BAD_PIXEL_ROW = 50 * SLIT_ROWS + 7


def _order(seed: int = 2, wavelengthScale: float = 1.0) -> pd.DataFrame:
    """A NIR-like order: a flat sky with three narrow lines, an object stripe, a cosmic ray and a bad pixel."""
    rng = np.random.default_rng(seed)
    column = np.repeat(np.arange(COLUMNS), SLIT_ROWS)
    row = np.tile(np.arange(SLIT_ROWS), COLUMNS)
    slit = rng.uniform(-5.5, 5.5, column.size)
    wavelength = 1500.0 + NM_PER_COLUMN * column + SLIT_TILT_NM * slit
    sky = np.full(column.shape, SKY_LEVEL)
    for centre in LINE_COLUMNS:
        sky = sky + LINE_AMPLITUDE * np.exp(
            -0.5 * ((wavelength - 1500.0 - NM_PER_COLUMN * centre) / (NM_PER_COLUMN * LINE_SIGMA_COLUMNS)) ** 2
        )
    objectFlux = OBJECT_AMPLITUDE * np.exp(-0.5 * ((slit - OBJECT_CENTRE) / OBJECT_WIDTH) ** 2)
    flux = rng.poisson(sky + objectFlux) + rng.normal(0.0, RON, column.size)
    flux = flux.astype(float)
    flux[COSMIC_ROW] += COSMIC_AMPLITUDE
    return pd.DataFrame(
        {
            "order": 12,
            "x": column,
            "y": row,
            "wavelength": wavelengthScale * wavelength,
            "slit_position": slit,
            "flux": flux,
            "error": np.sqrt(sky + objectFlux + RON**2),
            "mask": np.arange(column.size) == BAD_PIXEL_ROW,
        }
    )


def _clip(log: Any, order: pd.DataFrame, **skySettings: Any) -> pd.DataFrame:
    subtractor = _subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_clipping_sigma": 3,
            "percentile_clipping_iterations": 15,
            "percentile_rolling_window_size": 35,
            "noise_rolling_window_size": 250,
        }
    )
    subtractor.recipeSettings["sky-subtraction"].update(skySettings)
    return subtractor.get_over_sampled_sky_from_order(order.copy(), clipBPs=True, clipSlitEdge=False)


def _in_line(order: pd.DataFrame, centre: float, halfWidth: float = 1.5) -> np.ndarray:
    return (np.abs(order["x"].to_numpy() - centre) <= halfWidth) & ~order["mask"].to_numpy()


def _direct_rolling_clip(log: Any, order: pd.DataFrame) -> pd.DataFrame:
    """The rolling-percentile clip called directly, on an identical copy with the flag columns the caller sets."""
    subtractor = _subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    frame = order.copy()
    for column in ("flagged_all_clipped", "flagged_object_clipped", "flagged_bspline_clipped", "flagged_edge_clipped"):
        frame[column] = False
    frame["flagged_bad_pixel_clipped"] = frame["mask"].to_numpy(bool)
    frame["flagged_all_clipped"] = frame["mask"].to_numpy(bool)
    return subtractor.rolling_window_clipping(
        imageMapOrderDF=frame, windowSize=35, sigma_clip_limit=3, max_iterations=15
    )


def test_a_run_without_the_setting_is_the_rolling_percentile_clip_called_directly(log: Any) -> None:
    # ARRANGE
    order = _order()
    onLine = _in_line(order, LINE_COLUMNS[0])

    # ACT
    withoutKey = _clip(log, order)
    direct = _direct_rolling_clip(log, order)

    # ASSERT: THE SKYLINE PIXELS ARE FLAGGED, WHICH ONLY THE ROLLING CLIP DOES, AND THE WHOLE FRAME MATCHES
    assert withoutKey.loc[order.index[onLine], "flagged_object_clipped"].any()
    pd.testing.assert_frame_equal(withoutKey, direct)


def test_the_rolling_percentile_setting_is_the_rolling_percentile_clip_called_directly(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    withKey = _clip(log, order, **{METHOD_KEY: "rolling_percentile"})

    # ASSERT
    pd.testing.assert_frame_equal(withKey, _direct_rolling_clip(log, order))


def test_the_offset_bin_method_flags_a_cosmic_ray_as_clipped(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    clipped = _clip(log, order, **{METHOD_KEY: "offset_bin"})

    # ASSERT
    assert clipped.loc[COSMIC_ROW, "flagged_object_clipped"]
    assert clipped.loc[COSMIC_ROW, "flagged_all_clipped"]


@pytest.mark.parametrize("centre", LINE_COLUMNS)
def test_the_offset_bin_method_keeps_the_pixels_of_a_narrow_skyline(log: Any, centre: float) -> None:
    # ARRANGE
    order = _order()
    onLine = _in_line(order, centre)

    # ACT
    clipped = _clip(log, order, **{METHOD_KEY: "offset_bin"}).loc[order.index[onLine]]

    # ASSERT
    assert not clipped["flagged_all_clipped"].any()


def test_the_rolling_percentile_clip_flags_skyline_pixels_that_the_offset_bin_method_keeps(log: Any) -> None:
    # ARRANGE: THE DY-1972 PROBLEM, AS THE CONTRAST FOR THE TEST ABOVE
    order = _order()
    onLine = _in_line(order, LINE_COLUMNS[0])

    # ACT
    rolling = _clip(log, order).loc[order.index[onLine]]

    # ASSERT
    assert rolling["flagged_object_clipped"].any()


def test_bad_pixels_stay_clipped_but_are_not_counted_as_object_flags(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    clipped = _clip(log, order, **{METHOD_KEY: "offset_bin"})

    # ASSERT
    assert clipped.loc[BAD_PIXEL_ROW, "flagged_bad_pixel_clipped"]
    assert clipped.loc[BAD_PIXEL_ROW, "flagged_all_clipped"]
    assert not clipped.loc[BAD_PIXEL_ROW, "flagged_object_clipped"]


def test_a_bad_pixel_does_not_set_the_bin_statistics(log: Any) -> None:
    # ARRANGE: A BAD PIXEL HOLDING A HUGE VALUE WOULD SKEW ITS BIN IF THE FILTER READ IT
    order = _order()
    order.loc[BAD_PIXEL_ROW, "flux"] = 1e9
    reference = _clip(log, _order(), **{METHOD_KEY: "offset_bin"})

    # ACT
    clipped = _clip(log, order, **{METHOD_KEY: "offset_bin"})

    # ASSERT
    assert np.array_equal(clipped["flagged_object_clipped"].to_numpy(), reference["flagged_object_clipped"].to_numpy())


def test_the_offset_bin_method_returns_the_columns_the_object_masking_and_the_fit_read(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    clipped = _clip(log, order, **{METHOD_KEY: "offset_bin"})

    # ASSERT
    isClipped = clipped["flagged_all_clipped"].to_numpy(dtype=bool)
    assert clipped["wavelength"].is_monotonic_increasing
    for column in ("flux_percentile_smoothed", "flux_minus_smoothed_residual", "residual_global_sigma"):
        assert column in clipped
        assert clipped.loc[~isClipped, column].notna().sum() > 0.95 * (~isClipped).sum()
        assert clipped.loc[isClipped, column].isna().all()


def test_the_bin_width_follows_the_dispersion_of_the_order(log: Any) -> None:
    # ARRANGE: THE SAME ORDER ON A WAVELENGTH SCALE TEN TIMES COARSER HAS THE SAME PIXELS, SO THE SAME BINS
    order = _order()
    coarse = _order(wavelengthScale=10.0)

    # ACT
    flags = _clip(log, order, **{METHOD_KEY: "offset_bin"})["flagged_object_clipped"]
    coarseFlags = _clip(log, coarse, **{METHOD_KEY: "offset_bin"})["flagged_object_clipped"]

    # ASSERT
    assert flags.any()
    assert np.array_equal(flags.to_numpy(), coarseFlags.to_numpy())


def test_the_sigma_setting_sets_the_rejection_limit(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    standard = _clip(log, order, **{METHOD_KEY: "offset_bin"})
    unreachable = _clip(log, order, **{METHOD_KEY: "offset_bin", "offset_bin_sigma": 1e6})

    # ASSERT
    assert standard["flagged_object_clipped"].sum() > 0
    assert unreachable["flagged_object_clipped"].sum() == 0


def test_the_bin_width_setting_sets_the_bin_width_in_dispersion_pixels(log: Any) -> None:
    # ARRANGE: A BIN OF ONE PIXEL HOLDS 30 PIXELS OF ONE COLUMN, SO THE WHOLE LINE COLUMN IS ONE POPULATION AND THE
    # LINE PIXELS ARE NOT OUTLIERS; A BIN OF 5 PIXELS CUTS LINE WINGS DIFFERENTLY
    order = _order()

    # ACT
    narrow = _clip(log, order, **{METHOD_KEY: "offset_bin", "offset_bin_width_pixels": 1})
    standard = _clip(log, order, **{METHOD_KEY: "offset_bin"})

    # ASSERT
    assert not narrow["flagged_object_clipped"].equals(standard["flagged_object_clipped"])


def test_an_unknown_method_raises_a_value_error_naming_the_setting(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT / ASSERT
    with pytest.raises(ValueError, match=METHOD_KEY):
        _clip(log, order, **{METHOD_KEY: "kelson"})


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("offset_bin_sigma", 0),
        ("offset_bin_sigma", float("nan")),
        ("offset_bin_sigma", "3"),
        ("offset_bin_width_pixels", -5),
        ("offset_bin_width_pixels", None),
        ("offset_bin_width_pixels", True),
    ],
)
def test_a_bad_offset_bin_number_raises_a_value_error_naming_the_setting(log: Any, key: str, value: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT / ASSERT
    with pytest.raises(ValueError, match=key):
        _clip(log, order, **{METHOD_KEY: "offset_bin", key: value})


def test_the_object_stripe_still_gives_the_object_masking_a_range_to_mask(log: Any) -> None:
    # ARRANGE
    order = _order()
    subtractor = _subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_clipping_sigma": 3,
            "percentile_clipping_iterations": 15,
            "percentile_rolling_window_size": 35,
            "noise_rolling_window_size": 250,
            METHOD_KEY: "offset_bin",
        }
    )
    clipped = subtractor.get_over_sampled_sky_from_order(order.copy())

    # ACT
    masked = subtractor.clip_object_slit_positions([clipped], aggressive=True)[0]

    # ASSERT: THE STRIPE CENTRE IS MASKED ON EVERY COLUMN AWAY FROM THE LINES
    atStripe = (masked["slit_position"] - OBJECT_CENTRE).abs() < 0.3 * OBJECT_WIDTH
    assert masked.loc[atStripe, "flagged_all_clipped"].mean() > 0.9
