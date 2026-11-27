"""Contracts for the NIR sky correction as one velocity shift across the orders (DY-2477)."""

from __future__ import annotations

import copy
import pickle

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.nir_sky_refit import (
    MAX_GLOBAL_RMS_PX,
    MAX_TIGHT_TOLERANCE_NM,
    MIN_GLOBAL_LINES,
    fit_with_iterated_matching,
    nir_order_sky_lines,
    order_geometry,
    sky_tuning_settings,
)
from soxspipe.commonutils.nir_sky_velocity import (
    detector_px_to_velocity,
    fit_nir_sky_corrections,
    fit_velocity_model,
    sky_velocity,
    velocity_correction,
    velocity_to_detector_px,
)

pytestmark = pytest.mark.unit

VELOCITY = sky_velocity(centreNm=1400.0, halfRangeNm=500.0, coefficients=(1.0e-4, 2.0e-5))


# PIXEL-TO-VELOCITY CONVERSION


def test_a_detector_motion_is_the_velocity_that_moves_the_wavelength_by_the_same_nm():
    # ARRANGE: -0.06 nm/px TIMES 2 px IS A 0.12 nm CORRECTION, WHICH IS 1e-4 OF 1200 nm
    # ACT
    velocity = detector_px_to_velocity(2.0, -0.06, 1200.0)

    # ASSERT
    assert velocity == pytest.approx(0.06 * 2.0 / 1200.0)
    assert velocity == pytest.approx(1.0e-4)


def test_the_same_detector_motion_is_a_smaller_velocity_at_a_longer_wavelength():
    # ACT
    blue = detector_px_to_velocity(1.0, -0.06, 1000.0)
    red = detector_px_to_velocity(1.0, -0.06, 2000.0)

    # ASSERT
    assert blue == pytest.approx(2 * red)


@pytest.mark.parametrize("dispersion", [-0.09, 0.09])
def test_velocity_converts_back_to_the_detector_motion_it_came_from(dispersion):
    # ACT
    velocity = detector_px_to_velocity(-1.7, dispersion, 1450.0)

    # ASSERT
    assert velocity_to_detector_px(velocity, dispersion, 1450.0) == pytest.approx(-1.7)


@pytest.mark.parametrize("dispersion", [0.0, float("nan")])
def test_a_conversion_without_a_usable_dispersion_is_rejected(dispersion):
    # ACT / ASSERT
    with pytest.raises(ValueError, match="dispersion"):
        detector_px_to_velocity(1.0, dispersion, 1200.0)


# THE VELOCITY MODEL


def test_the_velocity_is_a_polynomial_in_the_normalised_wavelength():
    # ACT
    fraction = VELOCITY.fraction_at(np.array([1400.0, 1650.0, 900.0]))

    # ASSERT
    assert fraction == pytest.approx([1.0e-4, 1.0e-4 + 2.0e-5 * 0.5, 1.0e-4 - 2.0e-5])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"centreNm": 1400.0, "halfRangeNm": 0.0, "coefficients": (1.0e-4,)},
        {"centreNm": 1400.0, "halfRangeNm": -5.0, "coefficients": (1.0e-4,)},
        {"centreNm": 1400.0, "halfRangeNm": 500.0, "coefficients": ()},
        {"centreNm": 1400.0, "halfRangeNm": float("nan"), "coefficients": (1.0e-4,)},
    ],
)
def test_a_velocity_without_a_wavelength_range_or_coefficients_is_rejected(kwargs):
    # ACT / ASSERT
    with pytest.raises(ValueError, match="velocity"):
        sky_velocity(**kwargs)


# THE CORRECTION OF A VELOCITY ORDER


def _order(scaleNm: float, dispersion: float | None = None, centreNm: float = 1200.0) -> order_geometry:
    return order_geometry(
        centreNm=centreNm,
        halfRangeNm=15.0,
        dispersionNmPerPx=-scaleNm if dispersion is None else dispersion,
        medianPixelScaleNm=scaleNm,
    )


def test_a_velocity_order_is_corrected_by_the_wavelength_times_the_velocity():
    # ARRANGE
    correction = velocity_correction(12, VELOCITY, _order(0.08), nLines=9)
    wavelengths = np.array([1190.0, 1200.0, 1215.0])

    # ACT
    shiftNm = correction.correction_nm(wavelengths)

    # ASSERT
    assert shiftNm == pytest.approx(wavelengths * VELOCITY.fraction_at(wavelengths))


def test_a_velocity_order_ignores_the_pixel_scale_of_its_rows():
    # ARRANGE
    correction = velocity_correction(12, VELOCITY, _order(0.08), nLines=9)
    wavelengths = np.array([1190.0, 1200.0, 1215.0])

    # ACT
    withoutScales = correction.correction_nm(wavelengths)
    withScales = correction.correction_nm(wavelengths, np.array([0.05, 0.08, 0.2]))

    # ASSERT
    assert withScales == pytest.approx(withoutScales, abs=0.0)


def test_two_orders_with_different_pixel_scales_get_the_same_correction_at_the_same_wavelength():
    # ARRANGE: AN OVERLAP OF TWO ORDERS WHOSE PIXEL SCALES DIFFER BY 25 %
    narrow = velocity_correction(17, VELOCITY, _order(0.08, centreNm=1130.0), nLines=9)
    wide = velocity_correction(18, VELOCITY, _order(0.10, centreNm=1100.0), nLines=9)
    overlap = np.array([1110.0, 1115.0, 1120.0])

    # ACT
    narrowNm = narrow.correction_nm(overlap, np.full(3, 0.08))
    wideNm = wide.correction_nm(overlap, np.full(3, 0.10))

    # ASSERT
    assert narrowNm == pytest.approx(wideNm, abs=0.0)


def test_the_shift_at_the_centre_of_a_velocity_order_is_the_recorded_nm_and_px_shift():
    # ARRANGE
    geometry = _order(0.08)
    correction = velocity_correction(12, VELOCITY, geometry, nLines=9)

    # ACT / ASSERT
    expectedNm = geometry.centreNm * VELOCITY.fraction_at(geometry.centreNm)
    assert correction.shift_nm == pytest.approx(expectedNm)
    assert correction.shiftPx == pytest.approx(expectedNm / 0.08)
    assert velocity_to_detector_px(
        VELOCITY.fraction_at(geometry.centreNm), geometry.dispersionNmPerPx, geometry.centreNm
    ) == pytest.approx(correction.shiftPx)


def test_the_slope_of_a_velocity_order_is_the_px_change_per_unit_position_along_it():
    # ARRANGE
    geometry = _order(0.08)
    correction = velocity_correction(12, VELOCITY, geometry, nLines=9)
    ends = np.array([geometry.centreNm - geometry.halfRangeNm, geometry.centreNm + geometry.halfRangeNm])
    endsNm = correction.correction_nm(ends)

    # ACT
    slopePx = correction.slope_px

    # ASSERT
    assert slopePx == pytest.approx((endsNm[1] - endsNm[0]) / 2 / 0.08, rel=1e-6)


def test_a_velocity_correction_records_its_lines_and_is_a_global_one():
    # ACT
    correction = velocity_correction(12, VELOCITY, _order(0.08), nLines=9, globalSeedLines=7, ownSeedLines=3)

    # ASSERT
    assert (correction.order, correction.nLines, correction.source) == (12, 9, "global")
    assert (correction.globalSeedLines, correction.ownSeedLines) == (7, 3)
    assert correction.velocity == VELOCITY


def test_a_velocity_correction_survives_pickling_and_deep_copy():
    # ARRANGE
    correction = velocity_correction(12, VELOCITY, _order(0.08), nLines=9)

    # ACT
    restored = pickle.loads(pickle.dumps(correction))  # noqa: S301 - ROUND-TRIPS AN OBJECT THIS TEST BUILT
    copied = copy.deepcopy(correction)

    # ASSERT
    assert restored == correction
    assert copied == correction


# FITTING THE VELOCITY

ORDERS = list(range(12, 20))
RESOLVING_POWER = 12000.0
SETTINGS = sky_tuning_settings(ownCorrectionOrders=())
TRUTH = sky_velocity(centreNm=1400.0, halfRangeNm=500.0, coefficients=(1.2e-4, 3.0e-5))


def _scene_geometry() -> dict[int, order_geometry]:
    """*orders every 100 nm from 1100 nm whose pixel scale grows with wavelength, as the NIR orders do*"""
    return {
        o: order_geometry(
            centreNm=1100.0 + 100.0 * (o - 12),
            halfRangeNm=40.0,
            dispersionNmPerPx=-(1100.0 + 100.0 * (o - 12)) / RESOLVING_POWER,
            medianPixelScaleNm=(1100.0 + 100.0 * (o - 12)) / RESOLVING_POWER,
        )
        for o in ORDERS
    }


def _scene_lines(orders=ORDERS, perOrder=10, noisePx=0.01, seed=3) -> pd.DataFrame:
    """*matched lines whose catalogue wavelengths lie ``TRUTH`` away from the observed ones, with noise in px*"""
    rng = np.random.default_rng(seed)
    geometry = _scene_geometry()
    frames = []
    for o in orders:
        scale = geometry[o].medianPixelScaleNm
        observed = geometry[o].centreNm + np.linspace(-0.9, 0.9, perOrder) * geometry[o].halfRangeNm
        catalogue = observed * (1.0 + TRUTH.fraction_at(observed)) + rng.normal(0.0, noisePx * scale, perOrder)
        frames.append(
            pd.DataFrame({"order": o, "catalogueNm": catalogue, "observedNm": observed, "pixelScaleNm": scale})
        )
    return pd.concat(frames, ignore_index=True)


def test_the_fit_recovers_the_velocity_that_moved_every_line():
    # ARRANGE
    grid = np.linspace(1000.0, 1900.0, 7)

    # ACT
    result = fit_velocity_model(_scene_lines(), _scene_geometry(), SETTINGS)

    # ASSERT
    assert result.fallbackReason is None
    assert result.velocity.fraction_at(grid) == pytest.approx(TRUTH.fraction_at(grid), abs=3e-6)
    assert result.coefficients == result.velocity.coefficients
    assert result.nLines == 10 * len(ORDERS)
    assert result.rmsPx < 0.05


def test_every_order_gets_the_same_velocity_whatever_its_pixel_scale():
    # ARRANGE
    result = fit_velocity_model(_scene_lines(), _scene_geometry(), SETTINGS)
    wavelengths = np.array([1150.0, 1180.0])

    # ACT
    corrections = [result.corrections[o].correction_nm(wavelengths) for o in ORDERS]

    # ASSERT
    for correction in corrections:
        assert correction == pytest.approx(corrections[0], abs=0.0)
    assert corrections[0] == pytest.approx(wavelengths * TRUTH.fraction_at(wavelengths), abs=5e-4)


def test_an_order_without_lines_takes_the_velocity_at_its_wavelengths():
    # ARRANGE
    lines = _scene_lines(orders=[o for o in ORDERS if o != 15])

    # ACT
    result = fit_velocity_model(lines, _scene_geometry(), SETTINGS)

    # ASSERT
    assert result.corrections[15].nLines == 0
    assert result.corrections[15].shift_nm == pytest.approx(1400.0 * TRUTH.fraction_at(1400.0), abs=5e-4)


def test_a_line_far_from_the_velocity_is_clipped_and_does_not_move_the_fit():
    # ARRANGE
    lines = _scene_lines()
    lines.loc[7, "catalogueNm"] += 40 * lines.loc[7, "pixelScaleNm"]

    # ACT
    result = fit_velocity_model(lines, _scene_geometry(), SETTINGS)

    # ASSERT
    assert result.clipFraction == pytest.approx(1 / len(lines))
    assert result.lines["clipped"].sum() == 1
    assert bool(result.lines.loc[7, "clipped"])
    assert result.velocity.fraction_at(1400.0) == pytest.approx(TRUTH.fraction_at(1400.0), abs=3e-6)


def test_the_result_lines_carry_the_offset_before_and_the_model_after_in_px():
    # ACT
    result = fit_velocity_model(_scene_lines(), _scene_geometry(), SETTINGS)

    # ASSERT
    scale = result.lines["pixelScaleNm"]
    beforePx = (result.lines["catalogueNm"] - result.lines["observedNm"]) / scale
    assert result.lines["residualPx"].to_numpy() == pytest.approx(beforePx.to_numpy())
    assert (result.lines["residualPx"] - result.lines["modelPx"]).abs().max() < 0.05


def test_the_fit_judges_a_line_in_px_whatever_the_sign_of_its_orders_dispersion():
    # ARRANGE
    geometry = _scene_geometry()
    flipped = {
        o: order_geometry(g.centreNm, g.halfRangeNm, -g.dispersionNmPerPx, g.medianPixelScaleNm)
        for o, g in geometry.items()
    }

    # ACT
    result = fit_velocity_model(_scene_lines(), flipped, SETTINGS)

    # ASSERT
    assert result.fallbackReason is None
    assert result.velocity.fraction_at(1400.0) == pytest.approx(TRUTH.fraction_at(1400.0), abs=3e-6)
    assert result.lines["residualPx"].to_numpy() == pytest.approx(
        -((result.lines["catalogueNm"] - result.lines["observedNm"]) / result.lines["pixelScaleNm"]).to_numpy()
    )


def test_a_higher_degree_follows_a_curved_velocity():
    # ARRANGE
    curved = sky_velocity(1400.0, 500.0, (1.2e-4, 3.0e-5, 4.0e-5))
    lines = _scene_lines()
    observed = lines["observedNm"].to_numpy()
    lines["catalogueNm"] = observed * (1.0 + curved.fraction_at(observed))

    # ACT
    result = fit_velocity_model(lines, _scene_geometry(), sky_tuning_settings(ownCorrectionOrders=(), velocityDeg=2))

    # ASSERT
    assert result.fallbackReason is None
    assert result.velocity.fraction_at(np.array([1000.0, 1900.0])) == pytest.approx(
        curved.fraction_at(np.array([1000.0, 1900.0])), abs=1e-7
    )


def test_too_few_lines_falls_back_with_a_reason_and_no_corrections():
    # ARRANGE
    lines = _scene_lines(perOrder=2).head(MIN_GLOBAL_LINES - 1)

    # ACT
    result = fit_velocity_model(lines, _scene_geometry(), SETTINGS)

    # ASSERT
    assert result.corrections == {}
    assert result.velocity is None
    assert f"{MIN_GLOBAL_LINES - 1} sky lines matched" in result.fallbackReason


def test_lines_at_too_few_distinct_wavelengths_for_the_degree_fall_back_as_rank_deficient():
    # ARRANGE: TWO DISTINCT WAVELENGTHS CANNOT FIX A QUADRATIC
    lines = pd.concat([_scene_lines(perOrder=1).iloc[:2]] * MIN_GLOBAL_LINES, ignore_index=True)

    # ACT
    result = fit_velocity_model(lines, _scene_geometry(), sky_tuning_settings(ownCorrectionOrders=(), velocityDeg=2))

    # ASSERT
    assert result.corrections == {}
    assert "rank-deficient" in result.fallbackReason


def test_scatter_above_the_limit_falls_back_with_a_reason():
    # ARRANGE: LINES THAT ALTERNATE ABOUT THE VELOCITY BY MORE THAN THE LIMIT, SO NO CLIP REMOVES THEM
    lines = _scene_lines(noisePx=0.0)
    alternating = np.where(np.arange(len(lines)) % 2 == 0, 1.0, -1.0) * (MAX_GLOBAL_RMS_PX + 0.5)
    lines["catalogueNm"] += alternating * lines["pixelScaleNm"]

    # ACT
    result = fit_velocity_model(lines, _scene_geometry(), SETTINGS)

    # ASSERT
    assert result.corrections == {}
    assert "fit rms" in result.fallbackReason


def test_a_line_from_an_order_without_geometry_is_rejected():
    # ARRANGE
    geometry = _scene_geometry()
    del geometry[15]

    # ACT / ASSERT
    with pytest.raises(ValueError, match="order 15 has sky lines but no geometry"):
        fit_velocity_model(_scene_lines(), geometry, SETTINGS)


def test_the_fit_leaves_its_input_lines_unchanged():
    # ARRANGE
    lines = _scene_lines()
    before = lines.copy()

    # ACT
    fit_velocity_model(lines, _scene_geometry(), SETTINGS)

    # ASSERT
    pd.testing.assert_frame_equal(lines, before)


# THE WHOLE CORRECTION: ONE VELOCITY, AND THE PIXEL MODEL FOR THE ORDERS THAT KEEP THEIR OWN

ENTRY_ORDERS = list(range(10, 18))
OWN_ORDERS = (10, 11)
ENTRY_SETTINGS = sky_tuning_settings(ownCorrectionOrders=OWN_ORDERS)
# THE RIGID SEED (px) THAT PUTS EVERY ORDER'S LINES CLOSE TO WHERE THE VELOCITY MOVES THEM
SEED_PX = 1.7


def _entry_geometry() -> dict[int, order_geometry]:
    return {
        o: order_geometry(
            centreNm=1000.0 + 100.0 * (o - 10),
            halfRangeNm=40.0,
            dispersionNmPerPx=-(1000.0 + 100.0 * (o - 10)) / RESOLVING_POWER,
            medianPixelScaleNm=(1000.0 + 100.0 * (o - 10)) / RESOLVING_POWER,
        )
        for o in ENTRY_ORDERS
    }


def _entry_order_lines(order, perOrder=12, velocityOffsetPx=0.0, toleranceNm=MAX_TIGHT_TOLERANCE_NM):
    """*the lines of an order: where ``TRUTH`` puts them, or for an own-correction order a detector motion*"""
    geometry = _entry_geometry()[order]
    catalogue = geometry.centreNm + np.linspace(-0.9, 0.9, perOrder) * geometry.halfRangeNm
    if order in OWN_ORDERS:
        peaks = catalogue + geometry.dispersionNmPerPx * (SEED_PX + 0.3 * (order - 10))
    else:
        peaks = catalogue / (1.0 + TRUTH.fraction_at(catalogue)) + geometry.dispersionNmPerPx * velocityOffsetPx
    return nir_order_sky_lines(
        order=order,
        peakNm=peaks,
        catalogueNm=catalogue,
        pixelScaleNm=np.full(perOrder, geometry.medianPixelScaleNm),
        toleranceNm=toleranceNm,
    )


def _entry_scene(replacements=None) -> dict[int, nir_order_sky_lines]:
    """*the lines of every order, with the lines of the orders in ``replacements`` swapped in*"""
    return {o: (replacements or {}).get(o) or _entry_order_lines(o) for o in ENTRY_ORDERS}


def test_the_orders_that_keep_their_own_correction_get_exactly_the_pixel_models():
    # ARRANGE
    orderLines = _entry_scene()
    pixel = fit_with_iterated_matching(orderLines, _entry_geometry(), SEED_PX, ENTRY_SETTINGS)

    # ACT
    result = fit_nir_sky_corrections(orderLines, _entry_geometry(), SEED_PX, ENTRY_SETTINGS)

    # ASSERT
    assert result.fallbackReason is None
    for o in OWN_ORDERS:
        assert result.corrections[o] == pixel.corrections[o]
        assert result.corrections[o].velocity is None


def test_every_other_order_gets_the_one_velocity_that_moved_the_lines():
    # ACT
    result = fit_nir_sky_corrections(_entry_scene(), _entry_geometry(), SEED_PX, ENTRY_SETTINGS)

    # ASSERT
    others = [o for o in ENTRY_ORDERS if o not in OWN_ORDERS]
    assert result.velocity.fraction_at(np.array([1200.0, 1600.0])) == pytest.approx(
        TRUTH.fraction_at(np.array([1200.0, 1600.0])), abs=3e-6
    )
    for o in others:
        assert result.corrections[o].velocity == result.velocity
        assert result.corrections[o].nLines == 12


def test_neighbouring_velocity_orders_get_the_same_correction_where_they_overlap():
    # ARRANGE
    result = fit_nir_sky_corrections(_entry_scene(), _entry_geometry(), SEED_PX, ENTRY_SETTINGS)
    scales = {o: _entry_geometry()[o].medianPixelScaleNm for o in (14, 15)}
    wavelengths = np.array([1450.0, 1460.0])

    # ACT
    low = result.corrections[14].correction_nm(wavelengths, np.full(2, scales[14]))
    high = result.corrections[15].correction_nm(wavelengths, np.full(2, scales[15]))

    # ASSERT
    assert scales[14] != scales[15]
    assert low == pytest.approx(high, abs=0.0)


def test_the_result_counts_the_lines_and_scatter_of_both_fits():
    # ACT
    result = fit_nir_sky_corrections(_entry_scene(), _entry_geometry(), SEED_PX, ENTRY_SETTINGS)

    # ASSERT
    assert result.nLines == 12 * len(ENTRY_ORDERS)
    assert set(result.lines["order"]) == set(ENTRY_ORDERS)
    kept = ~result.lines["clipped"]
    expectedRms = np.sqrt(np.mean(((result.lines["residualPx"] - result.lines["modelPx"])[kept]) ** 2))
    assert result.rmsPx == pytest.approx(expectedRms)
    assert result.clipFraction == 0.0


def test_with_no_own_correction_orders_every_order_takes_the_velocity():
    # ARRANGE: THE LINES OF ORDERS 10 AND 11 SIT WHERE A DETECTOR MOTION PUT THEM, WHICH THE VELOCITY ONLY APPROXIMATES
    settings = sky_tuning_settings(ownCorrectionOrders=())

    # ACT
    result = fit_nir_sky_corrections(_entry_scene(), _entry_geometry(), SEED_PX, settings)

    # ASSERT
    assert result.fallbackReason is None
    assert all(result.corrections[o].velocity is not None for o in ENTRY_ORDERS)


def test_an_own_correction_order_that_is_not_in_the_frame_is_ignored():
    # ARRANGE
    settings = sky_tuning_settings(ownCorrectionOrders=(10, 11, 40))

    # ACT
    result = fit_nir_sky_corrections(_entry_scene(), _entry_geometry(), SEED_PX, settings)

    # ASSERT
    assert result.fallbackReason is None
    assert sorted(result.corrections) == ENTRY_ORDERS


def test_an_order_without_lines_takes_the_velocity_and_an_own_order_without_lines_the_pixel_model():
    # ARRANGE
    orderLines = {o: lines for o, lines in _entry_scene().items() if o not in (11, 15)}

    # ACT
    result = fit_nir_sky_corrections(orderLines, _entry_geometry(), SEED_PX, ENTRY_SETTINGS)

    # ASSERT
    assert result.corrections[15].velocity == result.velocity
    assert result.corrections[15].nLines == 0
    assert result.corrections[11].velocity is None


def test_a_velocity_order_with_few_lines_far_from_the_model_does_not_keep_its_own_shift():
    # ARRANGE: LIKE AN ORDER WITH THREE LINES THAT SIT 7 px OFF, WHICH THE PIXEL MODEL WOULD KEEP AS ITS OWN SEED
    orderLines = _entry_scene({14: _entry_order_lines(14, perOrder=3, velocityOffsetPx=7.0)})
    ownSeeds = {14: SEED_PX + 7.0}

    # ACT
    result = fit_nir_sky_corrections(orderLines, _entry_geometry(), SEED_PX, ENTRY_SETTINGS, ownSeedsPx=ownSeeds)

    # ASSERT
    assert result.corrections[14].velocity == result.velocity
    assert result.corrections[14].source == "global"
    assert result.velocity.fraction_at(1400.0) == pytest.approx(TRUTH.fraction_at(1400.0), abs=3e-6)


def test_the_seed_matches_of_a_velocity_order_are_recorded_on_its_correction():
    # ACT
    result = fit_nir_sky_corrections(
        _entry_scene(), _entry_geometry(), SEED_PX, ENTRY_SETTINGS, ownSeedsPx={14: SEED_PX}
    )

    # ASSERT
    assert result.corrections[14].globalSeedLines == 12
    assert result.corrections[14].ownSeedLines == 12
    assert result.corrections[15].ownSeedLines is None


def test_lines_beyond_the_tolerance_of_the_rigid_seed_are_matched_about_the_fitted_velocity():
    # ARRANGE: A SEED OF 1 px IS ABOUT 0.09 nm SHORT OF WHAT THE VELOCITY DOES AT ORDER 17, BEYOND ITS 0.04 nm TOLERANCE
    orderLines = _entry_scene({17: _entry_order_lines(17, toleranceNm=0.04)})

    # ACT
    firstRound = fit_nir_sky_corrections(orderLines, _entry_geometry(), 1.0, ENTRY_SETTINGS, iterations=1)
    result = fit_nir_sky_corrections(orderLines, _entry_geometry(), 1.0, ENTRY_SETTINGS)

    # ASSERT
    assert firstRound.corrections[17].nLines == 0
    assert result.corrections[17].nLines == 12


def test_too_few_lines_for_the_velocity_falls_back_with_a_reason():
    # ARRANGE: EVERY ORDER KEEPS ITS OWN CORRECTION, SO NO LINE IS LEFT FOR THE VELOCITY
    settings = sky_tuning_settings(ownCorrectionOrders=tuple(ENTRY_ORDERS))

    # ACT
    result = fit_nir_sky_corrections(_entry_scene(), _entry_geometry(), SEED_PX, settings)

    # ASSERT
    assert result.corrections == {}
    assert "sky lines matched" in result.fallbackReason


def test_a_pixel_fit_that_falls_back_makes_the_whole_correction_fall_back():
    # ACT
    result = fit_nir_sky_corrections(_entry_scene(), _entry_geometry(), None, ENTRY_SETTINGS)

    # ASSERT
    assert result.corrections == {}
    assert result.fallbackReason == "no order gave a seed shift"


def test_without_a_seed_the_velocity_alone_falls_back_too():
    # ACT
    result = fit_nir_sky_corrections(
        _entry_scene(), _entry_geometry(), None, sky_tuning_settings(ownCorrectionOrders=())
    )

    # ASSERT
    assert result.corrections == {}
    assert result.fallbackReason == "no order gave a seed shift"


def test_an_order_with_lines_but_no_geometry_is_rejected():
    # ARRANGE
    geometry = _entry_geometry()
    del geometry[14]

    # ACT / ASSERT
    with pytest.raises(ValueError, match="order 14 has sky lines but no geometry"):
        fit_nir_sky_corrections(_entry_scene(), geometry, SEED_PX, ENTRY_SETTINGS)


def test_an_order_with_lines_but_no_geometry_is_rejected_when_no_order_keeps_a_pixel_model():
    # ARRANGE
    geometry = _entry_geometry()
    del geometry[14]

    # ACT / ASSERT
    with pytest.raises(ValueError, match="order 14 has sky lines but no geometry"):
        fit_nir_sky_corrections(_entry_scene(), geometry, SEED_PX, sky_tuning_settings(ownCorrectionOrders=()))


def test_lines_without_a_pixel_scale_column_are_weighed_by_the_median_scale_of_their_order():
    # ARRANGE
    lines = _scene_lines().drop(columns="pixelScaleNm")

    # ACT
    result = fit_velocity_model(lines, _scene_geometry(), SETTINGS)

    # ASSERT
    assert result.fallbackReason is None
    assert result.velocity.fraction_at(1400.0) == pytest.approx(TRUTH.fraction_at(1400.0), abs=3e-6)
