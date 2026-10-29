"""Minimum knot spacing and chi-squared knot refinement of the B-spline sky fit (DY-2056)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate as ip

from soxspipe.commonutils.subtract_sky import (
    _chi_squared_refinement_knots,
    _dispersion_per_pixel,
    _enforce_minimum_knot_spacing,
)
from tests.unit.test_subtract_sky_bspline_characterization import (
    PIXEL_COUNT,
    _skyline_order,
    _subtractor,
)

pytestmark = pytest.mark.unit

SPACING_KEY = "bspline_min_knot_spacing_fraction"
CHI_SQUARED_KEY = "bspline_chi_squared_knot_refinement"
# THE FIXTURE ORDER RUNS FROM 500 TO 510 NM OVER PIXEL_COUNT ROWS
FIXTURE_DISPERSION = 10.0 / (PIXEL_COUNT - 1)
SPACING_FRACTION = 3.0


def _constant(spacing: float):
    return lambda wavelength: np.full(np.shape(wavelength), spacing)


def _dispersed_order(subtractor: Any) -> pd.DataFrame:
    """The skyline fixture with the detector coordinates the dispersion estimate reads."""
    pixels = _skyline_order()
    pixels["x"] = 0
    pixels["y"] = np.arange(PIXEL_COUNT)
    return pixels


def _knot_gaps(knots: np.ndarray) -> np.ndarray:
    return np.diff(np.sort(np.asarray(knots, dtype=float)))


# ----------------------------------------------------------------- A: minimum spacing


def test_enforce_minimum_knot_spacing_leaves_no_neighbours_closer_than_the_spacing() -> None:
    # ARRANGE
    knots = np.array([1.0, 1.05, 1.09, 2.0, 2.5, 2.52, 4.0])

    # ACT
    kept = _enforce_minimum_knot_spacing(knots, _constant(0.2))

    # ASSERT
    assert _knot_gaps(kept).min() >= 0.2
    assert list(kept) == [1.0, 2.0, 2.5, 4.0]


def test_enforce_minimum_knot_spacing_returns_the_knots_unchanged_when_all_are_far_enough() -> None:
    # ARRANGE
    knots = np.array([3.0, 1.0, 2.0])

    # ACT
    kept = _enforce_minimum_knot_spacing(knots, _constant(0.5))

    # ASSERT
    assert list(kept) == [1.0, 2.0, 3.0]


def test_enforce_minimum_knot_spacing_uses_the_spacing_at_each_knot() -> None:
    # ARRANGE: THE SPACING GROWS WITH WAVELENGTH, AS A COARSER DISPERSION WOULD
    knots = np.array([1.0, 1.3, 10.0, 10.3, 10.6, 12.5])

    # ACT
    kept = _enforce_minimum_knot_spacing(knots, lambda wavelength: 0.1 * wavelength)

    # ASSERT
    assert list(kept) == [1.0, 1.3, 10.0, 12.5]


def test_enforce_minimum_knot_spacing_returns_empty_for_no_knots() -> None:
    # ACT
    kept = _enforce_minimum_knot_spacing(np.array([]), _constant(0.5))

    # ASSERT
    assert kept.size == 0


def test_dispersion_per_pixel_is_the_wavelength_step_between_adjacent_rows() -> None:
    # ARRANGE: A TILTED SLIT SHIFTS THE WAVELENGTH OF EVERY ROW BY THE SLIT POSITION, NOT THE STEP BETWEEN ROWS
    rows = np.repeat(np.arange(400), 5)
    slit = np.tile(np.linspace(-1.0, 1.0, 5), 400)
    pixels = pd.DataFrame({"x": np.tile(np.arange(5), 400), "y": rows, "wavelength": 500.0 + 0.02 * rows + 0.5 * slit})

    # ACT
    spacingAt = _dispersion_per_pixel(pixels)

    # ASSERT
    assert spacingAt(np.array([501.0, 504.0, 507.0])) == pytest.approx(0.02, rel=1e-3)


def test_dispersion_per_pixel_follows_a_dispersion_that_changes_along_the_order() -> None:
    # ARRANGE
    rows = np.arange(1000)
    pixels = pd.DataFrame({"x": 0, "y": rows, "wavelength": 500.0 + 0.01 * rows + 2e-5 * rows**2})

    # ACT
    spacingAt = _dispersion_per_pixel(pixels)

    # ASSERT
    assert spacingAt(np.array([500.0]))[0] == pytest.approx(0.01, rel=0.05)
    assert spacingAt(np.array([pixels["wavelength"].max()]))[0] == pytest.approx(0.01 + 4e-5 * 999, rel=0.05)


def test_dispersion_per_pixel_finds_the_dispersion_axis_itself() -> None:
    # ARRANGE: THE WAVELENGTH RUNS ALONG x, ACROSS 400 COLUMNS, AND THE ORDER IS 20 ROWS (THE SLIT) THICK
    columns = np.repeat(np.arange(400), 20)
    pixels = pd.DataFrame({"x": columns, "y": np.tile(np.arange(20), 400), "wavelength": 800.0 - 0.05 * columns})

    # ACT
    spacingAt = _dispersion_per_pixel(pixels)

    # ASSERT
    assert spacingAt(np.array([790.0]))[0] == pytest.approx(0.05, rel=1e-3)


def test_fit_with_minimum_spacing_leaves_no_knots_closer_than_the_fraction_of_a_pixel(log: Any) -> None:
    # ARRANGE
    subtractor = _subtractor(log, iterationLimit=12)
    subtractor.recipeSettings["sky-subtraction"][SPACING_KEY] = SPACING_FRACTION
    pixels = _dispersed_order(subtractor)

    # ACT
    _, _, knots, _, _ = subtractor.fit_bspline_curve_to_sky(pixels)

    # ASSERT
    assert _knot_gaps(knots).min() >= SPACING_FRACTION * FIXTURE_DISPERSION * (1 - 1e-6)


def test_fit_without_minimum_spacing_places_knots_closer_than_the_fraction_of_a_pixel(log: Any) -> None:
    # ARRANGE: THE SAME ORDER AS ABOVE WITH THE SETTING ABSENT, SO THE TEST ABOVE BITES
    subtractor = _subtractor(log, iterationLimit=12)
    pixels = _dispersed_order(subtractor)

    # ACT
    _, _, knots, _, _ = subtractor.fit_bspline_curve_to_sky(pixels)

    # ASSERT
    assert _knot_gaps(knots).min() < SPACING_FRACTION * FIXTURE_DISPERSION


@pytest.mark.parametrize("settings", [{}, {SPACING_KEY: 0, CHI_SQUARED_KEY: False}])
def test_fit_with_both_prototypes_off_matches_a_fit_with_the_settings_absent(
    log: Any, settings: dict[str, Any]
) -> None:
    # ARRANGE
    reference = _subtractor(log, iterationLimit=12)
    candidate = _subtractor(log, iterationLimit=12)
    candidate.recipeSettings["sky-subtraction"].update(settings)

    # ACT
    modelRef, tckRef, knotsRef, _, _ = reference.fit_bspline_curve_to_sky(_dispersed_order(reference))
    modelOff, tckOff, knotsOff, _, _ = candidate.fit_bspline_curve_to_sky(_dispersed_order(candidate))

    # ASSERT
    np.testing.assert_array_equal(knotsOff, knotsRef)
    np.testing.assert_array_equal(tckOff[1], tckRef[1])
    np.testing.assert_array_equal(modelOff["sky_model"].to_numpy(), modelRef["sky_model"].to_numpy())


@pytest.mark.parametrize("value", [-1, True, "0.5", None, float("nan"), float("inf")])
def test_minimum_spacing_fraction_must_be_a_finite_number_of_zero_or_more(log: Any, value: Any) -> None:
    # ARRANGE
    subtractor = _subtractor(log, iterationLimit=12)
    subtractor.recipeSettings["sky-subtraction"][SPACING_KEY] = value
    pixels = _dispersed_order(subtractor)

    # ACT / ASSERT
    with pytest.raises(ValueError, match=SPACING_KEY):
        subtractor.fit_bspline_curve_to_sky(pixels)


@pytest.mark.parametrize("value", ["true", 1, None])
def test_chi_squared_refinement_setting_must_be_a_bool(log: Any, value: Any) -> None:
    # ARRANGE
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"][CHI_SQUARED_KEY] = value

    # ACT / ASSERT
    with pytest.raises(ValueError, match=CHI_SQUARED_KEY):
        subtractor._knot_placement_prototypes(_dispersed_order(subtractor))


def test_both_prototype_settings_absent_leave_the_knot_placement_alone_without_reading_the_map(log: Any) -> None:
    # ARRANGE: NO x OR y COLUMNS, SO ANY READ OF THE DISPERSION WOULD FAIL
    subtractor = _subtractor(log)

    # ACT
    placement = subtractor._knot_placement_prototypes(pd.DataFrame({"wavelength": [500.0, 501.0]}))

    # ASSERT
    assert placement == (None, None, False)


def test_an_order_with_no_dispersion_estimate_keeps_the_existing_knot_placement_and_warns(log: Any) -> None:
    # ARRANGE: A SINGLE PIXEL HAS NO DISPERSION, SO THE RECIPE MUST NOT ABORT
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"].update({SPACING_KEY: 0.5, CHI_SQUARED_KEY: True})

    # ACT
    placement = subtractor._knot_placement_prototypes(
        pd.DataFrame({"x": [3], "y": [3], "wavelength": [500.0], "order": 7})
    )

    # ASSERT
    assert placement == (None, None, False)
    assert any(level == "warning" and "order 7" in message for level, message in log.messages)


def test_chi_squared_refinement_without_a_spacing_setting_uses_a_floor_of_half_a_pixel(log: Any) -> None:
    # ARRANGE: 0.01 NM PER PIXEL, SO THE HIDDEN FLOOR IS 0.005 NM
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"][CHI_SQUARED_KEY] = True
    pixels = _dispersed_order(subtractor)
    pixels["wavelength"] = 500.0 + 0.01 * np.arange(PIXEL_COUNT)
    wavelength = np.linspace(500.0, 501.0, 1001)
    terms = np.full_like(wavelength, 9.0)

    # ACT
    minSpacingAt, refinementSpacingAt, isOn = subtractor._knot_placement_prototypes(pixels)
    narrow = _chi_squared_refinement_knots(
        np.array([500.4, 500.409]), wavelength, terms, refinementSpacingAt, minPoints=1
    )
    wide = _chi_squared_refinement_knots(
        np.array([500.4, 500.411]), wavelength, terms, refinementSpacingAt, minPoints=1
    )

    # ASSERT: THE 0.009 NM INTERVAL IS TOO NARROW TO HALVE TO 0.005 NM, THE 0.011 NM ONE IS NOT
    assert minSpacingAt is None
    assert isOn
    assert refinementSpacingAt(np.array([500.5]))[0] == pytest.approx(0.005)
    assert 500.4045 not in [round(knot, 4) for knot in narrow]
    assert 500.4055 in [round(knot, 4) for knot in wide]


def test_chi_squared_refinement_floor_follows_the_spacing_setting_when_it_is_on(log: Any) -> None:
    # ARRANGE
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"].update({CHI_SQUARED_KEY: True, SPACING_KEY: 2.0})
    pixels = _dispersed_order(subtractor)

    # ACT
    _, refinementSpacingAt, _ = subtractor._knot_placement_prototypes(pixels)

    # ASSERT
    assert refinementSpacingAt(np.array([505.0]))[0] == pytest.approx(2.0 * FIXTURE_DISPERSION, rel=1e-3)


# ----------------------------------------------------------------- B: chi-squared refinement


def _chi_terms(wavelength: np.ndarray, badInterval: tuple[float, float], badValue: float) -> np.ndarray:
    """Weighted squared residuals: 1 everywhere except the bad interval."""
    terms = np.ones_like(wavelength)
    terms[(wavelength > badInterval[0]) & (wavelength < badInterval[1])] = badValue
    return terms


def test_chi_squared_refinement_bisects_an_interval_with_reduced_chi_squared_above_two() -> None:
    # ARRANGE
    wavelength = np.linspace(0.0, 30.0, 3001)
    terms = _chi_terms(wavelength, (10.0, 20.0), 4.0)

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([10.0, 20.0]), wavelength, terms, _constant(0.1), minPoints=10)

    # ASSERT
    assert list(newKnots) == [15.0]


def test_chi_squared_refinement_leaves_an_interval_with_reduced_chi_squared_below_two() -> None:
    # ARRANGE
    wavelength = np.linspace(0.0, 30.0, 3001)
    terms = _chi_terms(wavelength, (10.0, 20.0), 1.5)

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([10.0, 20.0]), wavelength, terms, _constant(0.1), minPoints=10)

    # ASSERT
    assert newKnots.size == 0


def test_chi_squared_refinement_never_goes_below_the_minimum_spacing() -> None:
    # ARRANGE: THE BAD INTERVAL IS 10 WIDE, SO ITS HALVES ARE 5 WIDE
    wavelength = np.linspace(0.0, 30.0, 3001)
    terms = _chi_terms(wavelength, (10.0, 20.0), 4.0)
    knots = np.array([10.0, 20.0])

    # ACT
    allowed = _chi_squared_refinement_knots(knots, wavelength, terms, _constant(5.0), minPoints=10)
    refused = _chi_squared_refinement_knots(knots, wavelength, terms, _constant(5.1), minPoints=10)

    # ASSERT
    assert list(allowed) == [15.0]
    assert refused.size == 0


def test_chi_squared_refinement_bisects_every_bad_interval_and_only_those() -> None:
    # ARRANGE: THREE INTERVALS BETWEEN FOUR KNOTS; THE FIRST AND THIRD ARE BAD
    wavelength = np.linspace(0.0, 40.0, 4001)
    terms = _chi_terms(wavelength, (10.0, 20.0), 3.0)
    terms[(wavelength > 30.0) & (wavelength < 40.0)] = 3.0

    # ACT
    newKnots = _chi_squared_refinement_knots(
        np.array([10.0, 20.0, 30.0]), wavelength, terms, _constant(0.1), minPoints=10
    )

    # ASSERT
    assert list(newKnots) == [15.0, 35.0]


def test_chi_squared_refinement_leaves_an_interval_with_too_few_samples() -> None:
    # ARRANGE: 10 SAMPLES IN THE BAD INTERVAL AGAINST A FLOOR OF 50
    wavelength = np.linspace(0.0, 30.0, 31)
    terms = _chi_terms(wavelength, (10.0, 20.0), 4.0)

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([10.0, 20.0]), wavelength, terms, _constant(0.1), minPoints=50)

    # ASSERT
    assert newKnots.size == 0


def test_chi_squared_refinement_refines_the_end_intervals_out_to_the_data_range() -> None:
    # ARRANGE: THE BLUE END INTERVAL, FROM THE FIRST SAMPLE TO THE FIRST KNOT, IS BAD
    wavelength = np.linspace(0.0, 30.0, 3001)
    terms = _chi_terms(wavelength, (-1.0, 10.0), 5.0)

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([10.0, 20.0]), wavelength, terms, _constant(0.1), minPoints=10)

    # ASSERT
    assert list(newKnots) == [5.0]


def test_fit_with_chi_squared_refinement_adds_knots_at_skylines_without_the_residual_sum_rule(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE RESIDUAL-SUM RULE IS THE ONLY CALLER OF np.ma.masked_where IN THE FIT; SPY ON IT
    calls = {"legacy": 0}
    original = np.ma.masked_where

    def spy(*args: Any, **kwargs: Any) -> Any:
        calls["legacy"] += 1
        return original(*args, **kwargs)

    reference = _subtractor(log, iterationLimit=12)
    refined = _subtractor(log, iterationLimit=12)
    refined.recipeSettings["sky-subtraction"][CHI_SQUARED_KEY] = True
    monkeypatch.setattr(np.ma, "masked_where", spy)

    # ACT
    _, _, knotsRef, _, _ = reference.fit_bspline_curve_to_sky(_dispersed_order(reference))
    callsWithRuleOn = calls["legacy"]
    _, tck, knots, _, _ = refined.fit_bspline_curve_to_sky(_dispersed_order(refined))

    # ASSERT: THE SPY SEES THE LEGACY RULE WITH THE SETTING OFF AND NEVER WITH IT ON
    assert callsWithRuleOn > 0
    assert calls["legacy"] == callsWithRuleOn
    assert len(knots) != len(knotsRef)
    for centre in (502.0, 504.5, 507.3):
        assert np.any(np.abs(np.asarray(knots) - centre) < 0.1)
    assert np.all(np.isfinite(tck[1]))


def test_fit_with_chi_squared_refinement_keeps_the_knots_the_minimum_spacing_allows(log: Any) -> None:
    # ARRANGE
    subtractor = _subtractor(log, iterationLimit=12)
    subtractor.recipeSettings["sky-subtraction"].update({CHI_SQUARED_KEY: True, SPACING_KEY: SPACING_FRACTION})
    pixels = _dispersed_order(subtractor)

    # ACT
    _, _, knots, _, _ = subtractor.fit_bspline_curve_to_sky(pixels)

    # ASSERT
    assert _knot_gaps(knots).min() >= SPACING_FRACTION * FIXTURE_DISPERSION * (1 - 1e-6)


def _weighted_residual_order(weight: float) -> tuple[pd.DataFrame, tuple]:
    """A flat unit residual against a zero model, weighted by ``weight``, with a noisy region in the bad interval."""
    wavelength = np.linspace(0.0, 30.0, 3001)
    pixels = pd.DataFrame(
        {
            "wavelength": wavelength,
            "flux": 1.0,
            "slit_normalisation_ratio": 1.0,
            "weights": weight,
            "flagged_noisy_region": False,
        }
    )
    tck = ip.splrep(wavelength, np.zeros_like(wavelength), k=3)
    return pixels, tck


@pytest.mark.parametrize(("weight", "expectedKnots"), [(1.0, []), (2.0, [5.0, 15.0, 25.0])])
def test_chi_squared_uses_the_weights_of_the_fit(log: Any, weight: float, expectedKnots: list[float]) -> None:
    # ARRANGE: A UNIT RESIDUAL IS CHI-SQUARED 1 AT WEIGHT 1 AND 4 AT WEIGHT 2, IN ALL THREE INTERVALS
    subtractor = _subtractor(log)
    pixels, tck = _weighted_residual_order(weight)
    notClipped = pd.Series(False, index=pixels.index)

    # ACT
    newKnots = subtractor._chi_squared_new_knots(pixels, notClipped, np.array([10.0, 20.0]), tck, _constant(0.1), 10)

    # ASSERT
    assert list(newKnots) == [pytest.approx(knot) for knot in expectedKnots]


def test_chi_squared_ignores_clipped_and_noisy_pixels(log: Any) -> None:
    # ARRANGE: ONLY THE PIXELS BLUER THAN THE FIRST KNOT ARE NEITHER CLIPPED NOR NOISY
    subtractor = _subtractor(log)
    pixels, tck = _weighted_residual_order(2.0)
    isClipped = pd.Series((pixels["wavelength"] >= 10.0) & (pixels["wavelength"] < 20.0), index=pixels.index)
    pixels["flagged_noisy_region"] = pixels["wavelength"] >= 20.0

    # ACT
    newKnots = subtractor._chi_squared_new_knots(pixels, isClipped, np.array([10.0, 20.0]), tck, _constant(0.1), 10)

    # ASSERT
    assert list(newKnots) == [5.0]


def test_chi_squared_gives_no_knots_when_every_pixel_is_excluded(log: Any) -> None:
    # ARRANGE
    subtractor = _subtractor(log)
    pixels, tck = _weighted_residual_order(2.0)
    allClipped = pd.Series(True, index=pixels.index)

    # ACT
    newKnots = subtractor._chi_squared_new_knots(pixels, allClipped, np.array([10.0, 20.0]), tck, _constant(0.1), 10)

    # ASSERT
    assert newKnots.size == 0


def test_chi_squared_refinement_gives_no_knots_when_there_are_no_knots_to_bisect() -> None:
    # ARRANGE
    wavelength = np.linspace(0.0, 30.0, 3001)

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([]), wavelength, np.full(3001, 9.0), _constant(0.1), minPoints=10)

    # ASSERT
    assert newKnots.size == 0


def test_dispersion_per_pixel_raises_when_the_order_has_a_single_pixel() -> None:
    # ARRANGE
    pixels = pd.DataFrame({"x": [3], "y": [3], "wavelength": [500.0]})

    # ACT / ASSERT
    with pytest.raises(ValueError, match="at least two detector rows"):
        _dispersion_per_pixel(pixels)


def test_chi_squared_refinement_leaves_an_interval_whose_half_would_have_too_few_samples() -> None:
    # ARRANGE: THE PARENT (10, 20) HOLDS 203 SAMPLES, BUT ALL BUT THREE OF THEM LIE BLUER THAN ITS MIDPOINT
    wavelength = np.concatenate(
        (np.linspace(0.0, 10.0, 100), np.linspace(10.0, 14.9, 200), [15.5, 17.0, 19.0], np.linspace(20.0, 30.0, 100))
    )
    terms = _chi_terms(wavelength, (10.0, 20.0), 4.0)

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([10.0, 20.0]), wavelength, terms, _constant(0.1), minPoints=10)

    # ASSERT: AN EMPTY-ISH HALF WOULD LET THE KNOT DROPPER REMOVE AN ORIGINAL KNOT, SO THE INTERVAL IS NOT BISECTED
    assert newKnots.size == 0


def test_chi_squared_refinement_ignores_samples_that_are_not_finite() -> None:
    # ARRANGE
    wavelength = np.linspace(0.0, 30.0, 3001)
    terms = _chi_terms(wavelength, (10.0, 20.0), 4.0)
    wavelength[5] = np.nan
    terms[1500] = np.nan

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([10.0, 20.0]), wavelength, terms, _constant(0.1), minPoints=10)

    # ASSERT
    assert list(newKnots) == [15.0]


@pytest.mark.parametrize(("chiSquaredSum", "expectedKnots"), [(20.0, []), (20.5, [15.0])])
def test_chi_squared_refinement_divides_by_the_samples_less_one_and_needs_strictly_more_than_two(
    chiSquaredSum: float, expectedKnots: list[float]
) -> None:
    # ARRANGE: 11 SAMPLES IN THE INTERVAL, SO THE REDUCED CHI-SQUARED IS THE SUM OVER 10; 20 IS EXACTLY THE LIMIT
    wavelength = np.concatenate((np.arange(10.0, 20.0), [19.5]))
    terms = np.concatenate((np.full(10, 2.0), [0.0]))
    terms[-1] = chiSquaredSum - 20.0

    # ACT
    newKnots = _chi_squared_refinement_knots(np.array([10.0, 20.0]), wavelength, terms, _constant(0.1), minPoints=1)

    # ASSERT
    assert list(newKnots) == expectedKnots


def test_dispersion_per_pixel_never_falls_to_zero_where_the_cubic_folds() -> None:
    # ARRANGE: THE WAVELENGTH TURNS BACK AT ROW 500 (A FOLD), SO THE FITTED DERIVATIVE IS 0 THERE
    rows = np.arange(1001)
    pixels = pd.DataFrame({"x": 0, "y": rows, "wavelength": 500.0 + 1e-4 * (rows - 500.0) ** 2})

    # ACT
    spacingAt = _dispersion_per_pixel(pixels)

    # ASSERT: THE TYPICAL DISPERSION IS ABOUT 0.05 NM PER PIXEL; THE FLOOR KEEPS THE FOLD AT A TENTH OF IT
    values = spacingAt(np.linspace(500.0, 525.0, 50))
    assert np.all(np.isfinite(values))
    assert values.min() >= 0.1 * 0.05 * 0.9
