"""Contracts for the SOXS NIR global sky-line tuning in the Horne extraction (DY-1968)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import nir_sky_refit
from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.nir_sky_refit import MAX_TIGHT_TOLERANCE_NM, MIN_GLOBAL_LINES
from tests.unit._horne_nir_scene import (
    BLEND_GAP_NM,
    ORDERS,
    applied_shift,
    build_catalogue_extractor,
    build_extractor,
    build_orders,
    catalogue_nm,
    order_start_nm,
    true_shift_nm,
    true_shift_px,
)
from tests.unit._horne_vis_scene import qc_by_name

pytestmark = pytest.mark.unit

# THE SYNTHETIC SKY LINES ARE CENTRED TO ABOUT A TENTH OF A SAMPLE
SHIFT_TOLERANCE_NM = 0.01
DECOY_ORDER = 15


def _tune(extractor: horne_extraction, scene: pd.DataFrame, arm: str = "NIR") -> pd.DataFrame:
    return extractor._tune_or_skip(scene.copy(), arm)


def _spy_global_fit(monkeypatch: pytest.MonkeyPatch) -> list[dict]:
    calls: list[dict] = []
    original = nir_sky_refit.fit_with_iterated_matching

    def spy(orderLines, geometry, seedPx, settings, **kwargs):
        result = original(orderLines, geometry, seedPx, settings, **kwargs)
        calls.append({"orderLines": orderLines, "geometry": geometry, "seedPx": seedPx, "result": result, **kwargs})
        return result

    monkeypatch.setattr(nir_sky_refit, "fit_with_iterated_matching", spy)
    return calls


# DISPATCH


def test_soxs_nir_is_tuned_by_one_global_fit(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    scene = build_orders()
    calls = _spy_global_fit(monkeypatch)

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert len(calls) == 1
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines-global"
    for order in ORDERS:
        assert applied_shift(scene, tuned, order) == pytest.approx(true_shift_nm(order), abs=SHIFT_TOLERANCE_NM)


def test_global_tuning_records_the_fit_statistics_and_keeps_the_per_order_shift_qcs(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert qc_by_name(extractor, "SKY FIT RMS PX")["qc_value"] < 0.1
    assert qc_by_name(extractor, "SKY FIT NLINES")["qc_value"] == len(ORDERS) * len(catalogue_nm(ORDERS[0]))
    order = ORDERS[2]
    assert qc_by_name(extractor, f"SKY SHIFT O{order}")["qc_value"] == pytest.approx(true_shift_nm(order), abs=0.01)
    assert qc_by_name(extractor, f"SKY SHIFT O{order} PX")["qc_value"] == pytest.approx(true_shift_px(order), abs=0.2)


def test_global_tuning_keeps_the_corrections_for_the_not_flattened_pass(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert set(extractor.skyCorrections) == set(ORDERS)
    assert extractor.skyShiftsNm[ORDERS[0]] == pytest.approx(true_shift_nm(ORDERS[0]), abs=SHIFT_TOLERANCE_NM)


def test_x_shooter_nir_keeps_the_per_order_tuning(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, inst="XSH")
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert calls == []
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines"
    assert extractor.skyCorrections is None


def test_soxs_vis_keeps_the_per_order_tuning(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, arm="VIS")
    calls = _spy_global_fit(monkeypatch)
    scene = build_orders()
    scene["order"] = scene["order"] - 10

    # ACT
    extractor._tune_or_skip(scene.loc[scene["order"] == 1].copy(), "VIS")

    # ASSERT
    assert calls == []
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines"


def test_switching_the_global_fit_off_restores_the_per_order_tuning(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, sky_tuning={"global-fit": False})
    calls = _spy_global_fit(monkeypatch)
    scene = build_orders()

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert calls == []
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines"
    assert extractor.skyCorrections is None
    assert applied_shift(scene, tuned, ORDERS[0]) == pytest.approx(true_shift_nm(ORDERS[0]), abs=SHIFT_TOLERANCE_NM)


def test_an_attached_arc_still_skips_the_tuning(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    extractor._arc_attachment = lambda: (3.0, 1.0)
    calls = _spy_global_fit(monkeypatch)
    scene = build_orders()

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert calls == []
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skipped: attached arc"
    pd.testing.assert_frame_equal(tuned, scene)


@pytest.mark.parametrize(
    ("inst", "arm", "block", "expectsFit"),
    [
        ("SOXS", "NIR", None, True),
        ("SOXS", "NIR", {"global-fit": False}, False),
        ("SOXS", "VIS", None, False),
        ("XSH", "NIR", None, False),
    ],
)
def test_only_soxs_nir_with_the_global_fit_on_uses_it(log, tmp_path, inst, arm, block, expectsFit):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, sky_tuning=block, inst=inst, arm=arm)

    # ACT / ASSERT
    assert extractor._uses_nir_global_fit(arm) is expectsFit


def test_a_horne_extraction_without_sky_tuning_settings_uses_the_per_order_path():
    # ARRANGE
    extractor = object.__new__(horne_extraction)

    # ACT / ASSERT
    assert extractor._uses_nir_global_fit("NIR") is False


# SETTINGS


@pytest.mark.parametrize(
    ("inst", "arm", "block", "expected"),
    [
        ("SOXS", "NIR", {"clip-sigma": 2.5}, 2.5),
        ("SOXS", "NIR", None, 3.0),
        ("SOXS", "VIS", {"clip-sigma": 2.5}, None),
        ("XSH", "NIR", {"clip-sigma": 2.5}, None),
    ],
)
def test_the_sky_tuning_block_is_read_for_soxs_nir_only(log, tmp_path, inst, arm, block, expected):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, inst=inst, arm=arm)
    extractor.recipeSettings = {} if block is None else {"sky-tuning": block}

    # ACT
    settings = extractor._read_nir_sky_tuning()

    # ASSERT
    assert (None if settings is None else settings.clipSigma) == expected


def test_an_invalid_sky_tuning_block_fails_fast(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    extractor.recipeSettings = {"sky-tuning": {"clip-sigma": float("nan")}}

    # ACT / ASSERT
    with pytest.raises(ValueError, match="sky-tuning"):
        extractor._read_nir_sky_tuning()


# MATCHING


def test_a_decoy_catalogue_line_half_a_nanometre_away_does_not_pull_an_order_in_the_global_fit(log, tmp_path):
    # ARRANGE
    scene = build_orders()
    (tmp_path / "off").mkdir()
    globalExtractor = build_extractor(log, tmp_path, decoyOrders=(DECOY_ORDER,))
    perOrderExtractor = build_extractor(
        log, tmp_path / "off", sky_tuning={"global-fit": False}, decoyOrders=(DECOY_ORDER,)
    )

    # ACT
    globalTuned = _tune(globalExtractor, scene)
    perOrderTuned = _tune(perOrderExtractor, scene)

    # ASSERT
    expected = true_shift_nm(DECOY_ORDER)
    assert applied_shift(scene, globalTuned, DECOY_ORDER) == pytest.approx(expected, abs=SHIFT_TOLERANCE_NM)
    assert abs(applied_shift(scene, perOrderTuned, DECOY_ORDER) - expected) > 0.05


def test_the_seed_is_the_shift_of_the_line_rich_orders_not_of_an_aliased_one(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, decoyOrders=(DECOY_ORDER,))
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    seed = calls[0]["seedPx"]
    assert seed == pytest.approx(true_shift_px(DECOY_ORDER), abs=0.5)
    assert min(true_shift_px(o) for o in ORDERS) - 0.2 <= seed <= max(true_shift_px(o) for o in ORDERS) + 0.2


def test_the_match_tolerance_scales_with_the_arc_rms(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, arcRmsPx=0.2)
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    tolerances = [line.toleranceNm for line in calls[0]["orderLines"].values()]
    assert tolerances == pytest.approx([2.5 * 0.2 * 0.05] * len(tolerances))


def test_the_match_tolerance_is_the_cap_without_an_arc_rms_card(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, arcRmsPx=None)
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert {line.toleranceNm for line in calls[0]["orderLines"].values()} == {MAX_TIGHT_TOLERANCE_NM}


def test_an_order_with_a_faint_sky_gives_no_lines_but_still_takes_the_model_shift(log, tmp_path, monkeypatch):
    # ARRANGE
    faint = ORDERS[-1]
    extractor = build_extractor(log, tmp_path)
    scene = build_orders(linelessOrders=(faint,))
    calls = _spy_global_fit(monkeypatch)

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert faint not in calls[0]["orderLines"]
    assert applied_shift(scene, tuned, faint) == pytest.approx(true_shift_nm(faint), abs=0.02)
    assert extractor.skyCorrections[faint].nLines == 0


def test_too_few_lines_fall_back_to_the_per_order_tuning_and_say_so(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    scene = build_orders(linesOnly=ORDERS[:1])

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert len(catalogue_nm(ORDERS[0])) < MIN_GLOBAL_LINES + 1
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines-per-order"
    assert extractor.skyCorrections is None
    assert applied_shift(scene, tuned, ORDERS[0]) == pytest.approx(true_shift_nm(ORDERS[0]), abs=SHIFT_TOLERANCE_NM)
    assert f"SKY SHIFT O{ORDERS[0]}" in set(extractor.qc["qc_name"])


# CATALOGUE


def _catalogue_extractor(log, tmp_path, sky_tuning=None):
    extractor = build_extractor(log, tmp_path, sky_tuning=sky_tuning)
    del extractor._get_local_skylines_for_order
    extractor.skylinesDF = pd.DataFrame(
        {
            "WAVELENGTH": [1010.0, 1010.03, 1015.0],
            "FLUX": [300.0, 300.0, 400.0],
            "ISOLATED": [True, True, True],
        }
    )
    return extractor


@pytest.mark.parametrize("sky_tuning", [None, {"global-fit": False}])
def test_nir_calibration_lines_are_the_isolated_catalogue_lines_with_or_without_the_global_fit(
    log, tmp_path, sky_tuning
):
    # ARRANGE
    extractor = _catalogue_extractor(log, tmp_path, sky_tuning=sky_tuning)

    # ACT
    _, calibration = extractor._get_local_skylines_for_order(1005.0, 1020.0, 11, "wavelength")

    # ASSERT
    assert sorted(calibration.tolist()) == pytest.approx([1010.0, 1010.03, 1015.0])


# DEBUG PLOT


def test_debug_tuning_plots_the_global_fit_residuals_once(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    extractor.debug = True
    extractor._plot_skyline_shift_diagnostic = lambda *args: None
    plotted: list[object] = []
    extractor._plot_nir_sky_fit = plotted.append

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert len(plotted) == 1
    assert plotted[0].fallbackReason is None


def test_tuning_without_debug_draws_no_fit_plot(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    plotted: list[object] = []
    extractor._plot_nir_sky_fit = plotted.append

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert plotted == []


def test_the_fit_plot_shows_the_offset_before_and_after_the_model_for_every_line(log, tmp_path, monkeypatch):
    # ARRANGE
    import matplotlib
    import matplotlib.pyplot as plt

    from tests.unit._plot_spies import quiet_show, spy_figures

    extractor = build_extractor(log, tmp_path)
    extractor.debug = True
    extractor._plot_skyline_shift_diagnostic = lambda *args: None
    fits: list[object] = []
    extractor._plot_nir_sky_fit = fits.append
    _tune(extractor, build_orders())
    result = fits[0]
    del extractor._plot_nir_sky_fit
    monkeypatch.setattr(matplotlib, "use", lambda *args, **kwargs: None)
    figures = spy_figures(monkeypatch)
    shown = quiet_show(monkeypatch)

    # ACT
    extractor._plot_nir_sky_fit(result)

    # ASSERT
    assert len(shown) == 1
    before, after = figures[0].axes[:2]
    assert len(before.collections[0].get_offsets()) == len(result.lines)
    assert len(after.collections[0].get_offsets()) == len(result.lines)
    assert plt.get_fignums() == []


# EDGE CASES


def test_a_fallback_is_noted_for_the_qc_plots(log, tmp_path):
    # ARRANGE
    from soxspipe.commonutils.extraction_notes import extraction_note_logger

    noting = extraction_note_logger(log)
    extractor = build_extractor(noting, tmp_path)

    # ACT
    _tune(extractor, build_orders(linesOnly=ORDERS[:1]))

    # ASSERT
    assert any("Global sky-line fit not used" in message for _, message in noting.notes)


def test_an_order_without_a_dispersion_is_left_unshifted_and_records_a_zero_shift(log, tmp_path):
    # ARRANGE
    unsigned = ORDERS[-1]
    extractor = build_extractor(log, tmp_path)
    realDispersions = extractor._signed_order_dispersions
    extractor._signed_order_dispersions = lambda df, orders: {
        o: d for o, d in realDispersions(df, orders).items() if o != unsigned
    }
    scene = build_orders()

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, tuned, unsigned) == 0.0
    assert qc_by_name(extractor, f"SKY SHIFT O{unsigned}")["qc_value"] == 0.0
    assert applied_shift(scene, tuned, ORDERS[0]) == pytest.approx(true_shift_nm(ORDERS[0]), abs=SHIFT_TOLERANCE_NM)


def test_an_order_with_a_single_wavelength_is_left_unshifted(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    scene = build_orders()
    flat = scene["order"] == ORDERS[-1]
    scene.loc[flat, "wavelengthMean"] = 1500.0

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, tuned, ORDERS[-1]) == 0.0
    assert ORDERS[-1] not in extractor.skyCorrections


def test_an_unreadable_arc_header_gives_the_capped_tolerance(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    extractor._dispersion_map_header = lambda: None
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert {line.toleranceNm for line in calls[0]["orderLines"].values()} == {MAX_TIGHT_TOLERANCE_NM}


def test_an_order_with_no_detectable_peaks_supplies_no_lines(log, tmp_path, monkeypatch):
    # ARRANGE
    quiet = ORDERS[-1]
    quietStartNm = order_start_nm(quiet)
    extractor = build_extractor(log, tmp_path)
    realCentres = extractor._measure_peak_centres

    def centres(waveVals, skyVals, peaks):
        return [] if waveVals.min() >= quietStartNm else realCentres(waveVals, skyVals, peaks)

    extractor._measure_peak_centres = centres
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert quiet not in calls[0]["orderLines"]


# PER-ORDER OFFSETS AND OWN SEEDS


OFFSET_ORDER = 12
GENUINE_OFFSETS = {OFFSET_ORDER: 8.0}


def test_an_order_with_a_genuine_arc_offset_is_tuned_from_its_own_measurement(log, tmp_path):
    # ARRANGE: ORDER 12 SITS 8 px FROM THE DETECTOR SHIFT OF ITS NEIGHBOURS, 4 TIMES THE TIGHT TOLERANCE
    extractor = build_extractor(log, tmp_path)
    scene = build_orders(offsets=GENUINE_OFFSETS)

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines-global"
    assert extractor.skyCorrections[OFFSET_ORDER].perOrderOffsetPx > 6.0
    for order in ORDERS:
        expected = true_shift_nm(order, GENUINE_OFFSETS)
        assert applied_shift(scene, tuned, order) == pytest.approx(
            expected, abs=0.02 if order == OFFSET_ORDER else 0.01
        )


def test_each_order_offers_its_own_measured_shift_as_a_seed(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders(offsets=GENUINE_OFFSETS))

    # ASSERT
    ownSeeds = calls[0]["ownSeedsPx"]
    assert set(ownSeeds) == set(ORDERS)
    assert ownSeeds[OFFSET_ORDER] == pytest.approx(true_shift_px(OFFSET_ORDER, GENUINE_OFFSETS), abs=0.5)


# UNUSABLE ORDERS


def test_an_order_with_a_dispersion_but_a_single_wavelength_is_skipped_and_left_unshifted(log, tmp_path):
    # ARRANGE: THE DISPERSION MAP GIVES THE ORDER A DISPERSION, BUT ITS WAVELENGTH COLUMN HAS NO RANGE
    flat = ORDERS[-1]
    extractor = build_extractor(log, tmp_path)
    realDispersions = extractor._signed_order_dispersions
    extractor._signed_order_dispersions = lambda df, orders: {**realDispersions(df, orders), flat: -0.05}
    scene = build_orders()
    scene.loc[scene["order"] == flat, "wavelengthMean"] = 1500.0

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, tuned, flat) == 0.0
    assert flat not in extractor.skyCorrections
    assert applied_shift(scene, tuned, ORDERS[0]) == pytest.approx(true_shift_nm(ORDERS[0]), abs=SHIFT_TOLERANCE_NM)


def test_an_order_with_an_all_nan_sky_supplies_no_lines(log, tmp_path, monkeypatch):
    # ARRANGE
    blind = ORDERS[-1]
    extractor = build_extractor(log, tmp_path)
    scene = build_orders()
    scene.loc[scene["order"] == blind, "skyFlux"] = np.nan
    calls = _spy_global_fit(monkeypatch)

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert blind not in calls[0]["orderLines"]
    assert applied_shift(scene, tuned, blind) == pytest.approx(true_shift_nm(blind), abs=0.02)


def test_an_order_with_missing_pixel_scales_does_not_stop_the_fit(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    scene = build_orders()
    patchy = (scene["order"] == ORDERS[2]) & (scene.index % 3 == 0)
    scene.loc[patchy, "pixelScaleNm"] = np.nan

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines-global"
    assert np.isfinite(tuned["wavelengthMean"]).all()


def test_an_order_with_only_missing_pixel_scales_is_left_unshifted(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    scene = build_orders()
    scene.loc[scene["order"] == ORDERS[2], "pixelScaleNm"] = np.nan

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, tuned, ORDERS[2]) == 0.0
    assert ORDERS[2] not in extractor.skyCorrections


def test_global_tuning_raises_no_numpy_warnings(log, tmp_path):
    # ARRANGE
    import warnings

    extractor = build_extractor(log, tmp_path)

    # ACT / ASSERT
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _tune(extractor, build_orders(linelessOrders=(ORDERS[-1],)))


def test_the_clipped_shift_of_an_order_with_no_matches_is_zero_without_a_warning(log, tmp_path):
    # ARRANGE
    import warnings

    extractor = build_extractor(log, tmp_path)

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _, _, medianShift, keptShifts = extractor._compute_clipped_shift([], [])

    # ASSERT
    assert medianShift == 0.0
    assert len(keptShifts) == 0


# ARC RESIDUAL ALONG THE DISPERSION


def test_the_arc_rms_is_read_along_the_dispersion_direction_of_the_map(log, tmp_path):
    # ARRANGE: IN THE DISPERSION MAP'S PIXEL COORDINATES NIR WAVELENGTH RUNS ALONG x (axisB), THE OPPOSITE OF THE
    # RAW-FRAME dispersion-axis PARAMETER (y), SO THE ARC RESIDUAL ALONG THE DISPERSION IS THE X ONE
    from astropy.io import fits

    extractor = build_extractor(log, tmp_path, arcRmsPx=0.082)
    fits.setval(extractor.dispersionMap, "ESO QC X RES MEDIAN", value=0.189)
    extractor.axisA, extractor.axisB = "y", "x"

    # ACT
    arcRms = extractor._arc_rms_px()

    # ASSERT
    assert arcRms == pytest.approx(0.189)


def test_a_large_arc_rms_is_still_capped(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_extractor(log, tmp_path, arcRmsPx=5.0)
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders())

    # ASSERT
    assert {line.toleranceNm for line in calls[0]["orderLines"].values()} == {MAX_TIGHT_TOLERANCE_NM}


# AN ORDER WITH TOO FEW LINES KEEPS ITS OWN SEED


def test_an_order_with_three_lines_and_a_genuine_arc_offset_keeps_its_own_measured_shift(log, tmp_path):
    # ARRANGE: LIKE THE REAL O10: 3 SKY LINES, 8 px FROM THE DETECTOR SHIFT
    extractor = build_extractor(log, tmp_path)
    scene = build_orders(offsets=GENUINE_OFFSETS, lineCounts={OFFSET_ORDER: 3})

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    assert extractor.skyCorrections[OFFSET_ORDER].source == "own-seed"
    for order in ORDERS:
        expected = true_shift_nm(order, GENUINE_OFFSETS)
        assert applied_shift(scene, tuned, order) == pytest.approx(
            expected, abs=0.02 if order == OFFSET_ORDER else 0.01
        )
    assert any("own measured shift" in message for _, message in log.messages)


def test_the_log_names_where_each_orders_shift_came_from(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)

    # ACT
    _tune(extractor, build_orders(offsets=GENUINE_OFFSETS))

    # ASSERT
    printed = [message for level, message in log.messages if level == "print" and "applied to order" in message]
    assert any("(own-offset" in message for message in printed)
    assert any("(global" in message for message in printed)


# AN ORDER WHOSE LINES MATCH ABOUT NEITHER SEED

# LINES SCATTERED +-1.5 AND +-3.5 px ABOUT THE ORDER'S OWN SHIFT, A CURVED ARC-SOLUTION ERROR LIKE O10 IN A 1 ARCSEC NOD
CURVED_LINE_OFFSETS = {OFFSET_ORDER: 8.0 + np.array([-3.5, -1.5, 1.5, 3.5])}


def test_an_order_with_curved_lines_that_match_about_neither_seed_keeps_its_own_measured_shift(
    log, tmp_path, monkeypatch
):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)
    scene = build_orders(offsets=CURVED_LINE_OFFSETS, lineCounts={OFFSET_ORDER: 4})
    calls = _spy_global_fit(monkeypatch)

    # ACT
    tuned = _tune(extractor, scene)

    # ASSERT
    correction = extractor.skyCorrections[OFFSET_ORDER]
    ownSeed = calls[0]["ownSeedsPx"][OFFSET_ORDER]
    assert correction.source == "own-seed"
    # WITHIN THE 0.1 nm (2 px) CAP THE OWN SEED MATCHES THE TWO LINES 1.5 px AWAY AND THE GLOBAL SEED NONE, BUT NONE IS
    # WITHIN THE 0.5 px TIGHT TOLERANCE, SO THE OWN SEED WINS AND KEEPS THE ORDER
    assert (correction.globalSeedLines, correction.ownSeedLines) == (0, 2)
    assert correction.shiftPx == pytest.approx(ownSeed, abs=1e-9)
    assert ownSeed != pytest.approx(calls[0]["seedPx"], abs=1.0)
    assert applied_shift(scene, tuned, ORDERS[0]) == pytest.approx(true_shift_nm(ORDERS[0]), abs=SHIFT_TOLERANCE_NM)


def test_the_log_gives_the_seed_match_counts_of_each_order(log, tmp_path):
    # ARRANGE
    extractor = build_extractor(log, tmp_path)

    # ACT
    _tune(extractor, build_orders(offsets=CURVED_LINE_OFFSETS, lineCounts={OFFSET_ORDER: 4}))

    # ASSERT
    printed = [message for level, message in log.messages if level == "print" and "applied to order" in message]
    assert any("(own-seed; seed matches within 0.1 nm: global 0, own 2)" in message for message in printed)
    assert any("(global; seed matches within 0.1 nm: global 10, own 10)" in message for message in printed)


# THE FIRST PASS IS DEVELOP'S: THE ISOLATED-LINE CATALOGUE


def _first_pass(extractor, scene):
    orders = np.sort(scene["order"].unique())
    dispersions = extractor._signed_order_dispersions(scene, orders)
    shiftsNm, orderStats = extractor._measure_orders(scene.copy(), orders, "NIR", False, dispersions)
    return shiftsNm, {o: stats["nLines"] for o, stats in orderStats.items()}


def test_the_global_fit_first_pass_measures_every_order_exactly_as_the_per_order_tuning_does(log, tmp_path):
    # ARRANGE: THE FIRST LINE OF EVERY ORDER IS A DOUBLET, TWO ISOLATED LINES BUT ONE MERGED CENTROID
    for name in ("on", "off"):
        (tmp_path / name).mkdir()
    scene = build_orders(blendedFirstLine=True)
    globalOn = build_catalogue_extractor(log, tmp_path / "on", blendedFirstLine=True)
    globalOff = build_catalogue_extractor(
        log, tmp_path / "off", sky_tuning={"global-fit": False}, blendedFirstLine=True
    )

    # ACT
    onShifts, onCounts = _first_pass(globalOn, scene)
    offShifts, offCounts = _first_pass(globalOff, scene)

    # ASSERT
    assert onCounts == offCounts
    assert set(onCounts.values()) == {len(catalogue_nm(ORDERS[0])) + 1}
    assert onShifts == pytest.approx(offShifts, abs=0)


def test_the_global_fit_tight_matching_receives_the_isolated_catalogue_lines(log, tmp_path, monkeypatch):
    # ARRANGE: THE FIRST LINE OF EVERY ORDER IS A DOUBLET, LISTED AS TWO ISOLATED LINES
    extractor = build_catalogue_extractor(log, tmp_path, blendedFirstLine=True)
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders(blendedFirstLine=True))

    # ASSERT
    for order, lines in calls[0]["orderLines"].items():
        assert len(lines.catalogueNm) == len(catalogue_nm(order)) + 1


def test_the_global_fit_snr_gate_runs_on_the_isolated_catalogue_lines(log, tmp_path, monkeypatch):
    # ARRANGE
    extractor = build_catalogue_extractor(log, tmp_path, blendedFirstLine=True)
    seen: list[tuple] = []
    realSnrs = extractor._order_sky_line_snrs

    def record(*args):
        seen.append(args[2:])
        return realSnrs(*args)

    extractor._order_sky_line_snrs = record

    # ACT
    _tune(extractor, build_orders(blendedFirstLine=True))

    # ASSERT: THE GATE IS ASKED FOR NO CATALOGUE CHOICE, SO IT TAKES THE DEFAULT, THE ISOLATED LINES
    assert seen
    assert set(seen) == {()}


def test_a_doublet_peak_is_matched_once_not_to_both_components(log, tmp_path, monkeypatch):
    # ARRANGE: THE FIRST SKY LINE OF EVERY ORDER IS A BLEND HALFWAY BETWEEN TWO ISOLATED CATALOGUE COMPONENTS
    extractor = build_catalogue_extractor(log, tmp_path, blendedFirstLine=True)
    calls = _spy_global_fit(monkeypatch)

    # ACT
    _tune(extractor, build_orders(blendedFirstLine=True))

    # ASSERT: 11 CATALOGUE LINES BUT 10 PEAKS, SO THE FIT GETS ONE ROW PER PEAK. ONLY ONE COMPONENT IS MATCHED
    lines = calls[0]["result"].lines
    assert len(lines) == len(ORDERS) * len(catalogue_nm(ORDERS[0]))
    for order, group in lines.groupby("order"):
        components = [catalogue_nm(order)[0], catalogue_nm(order)[0] + BLEND_GAP_NM]
        assert sum(np.isclose(group["catalogueNm"], c).any() for c in components) == 1
