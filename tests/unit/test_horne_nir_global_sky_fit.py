"""Contracts for the SOXS NIR global sky-line tuning in the Horne extraction (DY-1968)."""

from __future__ import annotations

import pandas as pd
import pytest

from soxspipe.commonutils import nir_sky_refit
from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.nir_sky_refit import MAX_TIGHT_TOLERANCE_NM, MIN_GLOBAL_LINES
from soxspipe.commonutils.toolkit import SKYLINE_MIN_LINE_FLUX
from tests.unit._horne_nir_scene import (
    ORDERS,
    applied_shift,
    build_extractor,
    build_orders,
    catalogue_nm,
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
        calls.append({"orderLines": orderLines, "geometry": geometry, "seedPx": seedPx})
        return original(orderLines, geometry, seedPx, settings, **kwargs)

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
    extractor._mergedSkylinesDF = pd.DataFrame({"WAVELENGTH": [1010.015, 1015.0], "FLUX": [600.0, 400.0]})
    return extractor


def test_global_nir_tuning_matches_blend_merged_calibration_lines(log, tmp_path):
    # ARRANGE
    extractor = _catalogue_extractor(log, tmp_path)

    # ACT
    _, calibration = extractor._get_local_skylines_for_order(1005.0, 1020.0, 11, "wavelength")

    # ASSERT
    assert sorted(calibration.tolist()) == pytest.approx([1010.015, 1015.0])


def test_per_order_nir_tuning_keeps_the_isolated_catalogue_lines(log, tmp_path):
    # ARRANGE
    extractor = _catalogue_extractor(log, tmp_path, sky_tuning={"global-fit": False})

    # ACT
    _, calibration = extractor._get_local_skylines_for_order(1005.0, 1020.0, 11, "wavelength")

    # ASSERT
    assert sorted(calibration.tolist()) == pytest.approx([1010.0, 1010.03, 1015.0])


@pytest.mark.parametrize(("arm", "keyword"), [("NIR", "minBrightnessNIR"), ("VIS", "minBrightnessVIS")])
def test_the_blend_merge_reads_the_catalogue_without_its_brightness_cut(log, tmp_path, monkeypatch, arm, keyword):
    # ARRANGE
    import soxspipe.commonutils.toolkit as toolkit

    extractor = build_extractor(log, tmp_path, arm=arm)
    extractor._mergedSkylinesDF = None
    seen: list[dict] = []

    def fake(log, settings, arm, **kwargs):
        seen.append({"arm": arm, **kwargs})
        return pd.DataFrame({"WAVELENGTH": [1010.0, 1010.03], "FLUX": [10.0, 10.0]})

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", fake)
    extractor.settings = {}

    # ACT
    merged = extractor._merged_skylines()

    # ASSERT
    assert seen == [{"arm": arm, keyword: SKYLINE_MIN_LINE_FLUX}]
    assert len(merged) == 1
    assert merged["WAVELENGTH"].iloc[0] == pytest.approx(1010.015)
