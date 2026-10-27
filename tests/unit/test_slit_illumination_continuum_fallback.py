"""Contracts for the continuum-edge fallback of a VIS order whose skyline fit fails (DY-1927)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import slit_illumination
from soxspipe.commonutils.slit_illumination import fit_order_offsets
from tests.factories import synthetic_illumination_scene
from tests.unit.test_slit_illumination_continuum import _faint_lines_scene, _skyline_inputs
from tests.unit.test_slit_illumination_offset import _measure, _offsets

pytestmark = pytest.mark.unit

# HALF THE LINE SPACING OF THE SCENE'S SKYLINE TABLE, SO A SHIFTED TABLE LINE FALLS ON THE CONTINUUM BETWEEN TWO LINES
HALF_LINE_SPACING_NM = 0.4


def _mixed_scene(tmp_path, *, deltas=None, slitName="SLIT1.0"):
    """*two orders: the table lines of order 11 are in the frame, those of order 12 fall between the frame's lines*"""
    scene = synthetic_illumination_scene(tmp_path, offsets=deltas or {11: 1.0, 12: 0.8}, slitName=slitName)
    lines = scene.skylinesDF.copy()
    isOrder12 = lines["WAVELENGTH"] > 550.0
    lines.loc[isOrder12, "WAVELENGTH"] += HALF_LINE_SPACING_NM
    return replace(scene, skylinesDF=lines)


def test_a_vis_order_whose_skyline_fit_fails_gets_the_continuum_offset_and_records_the_method(
    tmp_path, log, monkeypatch
):
    scene = _faint_lines_scene(tmp_path, slitName="SLIT1.0")

    result = _measure(scene, tmp_path, log, monkeypatch)

    row = _offsets(result).loc[11]
    assert row["status"] == "fitted"
    assert row["method"] == "continuum"
    assert row["offset"] == pytest.approx(0.8, abs=0.1)
    assert row["reason"] == ""
    assert result.isApplied


def test_a_vis_order_whose_skyline_fit_passes_never_builds_a_continuum_profile(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: -1.0}, slitName="SLIT1.0")

    def forbidden(*args, **kwargs):
        raise AssertionError("the continuum builder must not run when every skyline fit passes")

    monkeypatch.setattr(slit_illumination, "continuum_slit_profiles", forbidden)
    offsets = _offsets(_measure(scene, tmp_path, log, monkeypatch))

    assert list(offsets["method"]) == ["skyline", "skyline"]
    assert list(offsets["status"]) == ["fitted", "fitted"]


def test_the_fallback_leaves_the_offset_of_an_order_whose_skyline_fit_passes_exactly_as_it_was(
    tmp_path, log, monkeypatch
):
    scene = _mixed_scene(tmp_path)
    skylineOnly = fit_order_offsets(*_skyline_inputs(scene, scene.skylinesDF)).set_index("order")

    offsets = _offsets(_measure(scene, tmp_path, log, monkeypatch))

    assert skylineOnly.loc[12, "status"] == "failed"
    assert offsets.loc[11, "method"] == "skyline"
    assert offsets.loc[11, "offset"] == skylineOnly.loc[11, "offset"]
    assert offsets.loc[11, "offsetErr"] == skylineOnly.loc[11, "offsetErr"]
    assert offsets.loc[12, "method"] == "continuum"
    assert offsets.loc[12, "status"] == "fitted"
    assert offsets.loc[12, "offset"] == pytest.approx(0.8, abs=0.1)


def test_a_failed_continuum_fallback_leaves_the_order_unmeasured_with_both_reasons(tmp_path, log, monkeypatch):
    scene = _faint_lines_scene(tmp_path, slitName="SLIT1.0")

    result = _measure(
        scene, tmp_path, log, monkeypatch, illuminationSettings={"continuum-halves-max-diff-arcsec": 1e-9}
    )

    row = _offsets(result).loc[11]
    assert row["status"] == "failed"
    assert "halves" in row["reason"]
    assert "too faint" in row["reason"]
    assert not result.isApplied


def test_a_failed_continuum_fallback_does_not_borrow_from_a_partner_order(tmp_path, log, monkeypatch):
    scene = _mixed_scene(tmp_path, deltas={11: 1.0, 12: 0.8})

    result = _measure(
        scene,
        tmp_path,
        log,
        monkeypatch,
        illuminationSettings={"continuum-halves-max-diff-arcsec": 1e-9},
        partners={12: 11},
    )

    offsets = _offsets(result)
    assert offsets.loc[11, "status"] == "fitted"
    assert offsets.loc[12, "status"] == "failed"
    assert pd.isna(offsets.loc[12, "borrowedFrom"])


def test_a_nir_order_whose_skyline_fit_fails_never_builds_a_continuum_profile(tmp_path, log, monkeypatch):
    scene = _faint_lines_scene(tmp_path, slitName="SLIT1.0")

    def forbidden(*args, **kwargs):
        raise AssertionError("the skyline-failure fallback is for the VIS arm only")

    monkeypatch.setattr(slit_illumination, "continuum_slit_profiles", forbidden)
    result = _measure(scene, tmp_path, log, monkeypatch, arm="NIR")

    row = _offsets(result).loc[11]
    assert row["status"] == "failed"
    assert row["method"] == "skyline"
    assert "too faint" in row["reason"]
    assert not result.isApplied


@pytest.mark.parametrize("slitName", ["SLIT5.0", "SLIT1.5", None])
def test_a_vis_order_whose_skyline_fit_fails_is_not_measured_from_the_continuum_on_other_slits(
    tmp_path, log, monkeypatch, slitName
):
    scene = _faint_lines_scene(tmp_path, slitName=slitName)

    def forbidden(*args, **kwargs):
        raise AssertionError("the skyline-failure fallback must run on SLIT1.0 only")

    monkeypatch.setattr(slit_illumination, "continuum_slit_profiles", forbidden)
    row = _offsets(_measure(scene, tmp_path, log, monkeypatch)).loc[11]

    assert row["status"] == "failed"
    assert row["method"] == "skyline"


def test_the_continuum_order_list_adds_the_failed_skyline_orders_on_the_slit_1_0_only(log):
    notes = {10: "the skyline table has no lines between 500 and 512 nm", 13: "every skyline is too faint"}

    onSlit = slit_illumination._continuum_orders(log, notes, "SLIT1.0", failedSkylineOrders=[13, 14])
    onWideSlit = slit_illumination._continuum_orders(log, notes, "SLIT5.0", failedSkylineOrders=[13, 14])

    assert onSlit == [10, 13, 14]
    assert onWideSlit == []


def test_failed_skyline_orders_left_off_the_continuum_on_another_slit_are_logged(log):
    # ARRANGE
    notes = {13: "every skyline is too faint"}

    # ACT
    orders = slit_illumination._continuum_orders(log, notes, "SLIT5.0", failedSkylineOrders=[14, 13])

    # ASSERT
    infos = [message for kind, message in log.messages if kind == "info"]
    assert orders == []
    assert any("SLIT5.0" in message and "[13, 14]" in message for message in infos)


@pytest.mark.parametrize(
    ("slitName", "failedSkylineOrders"), [("SLIT1.0", [13]), ("SLIT5.0", [])], ids=["slit-1-0", "no-failed-orders"]
)
def test_no_failed_skyline_skip_is_logged_when_nothing_is_left_off(log, slitName, failedSkylineOrders):
    # ACT
    slit_illumination._continuum_orders(log, {}, slitName, failedSkylineOrders=failedSkylineOrders)

    # ASSERT
    assert not any("keep their failed skyline result" in message for _, message in log.messages)


def test_the_failed_skyline_orders_are_the_failed_rows_of_the_vis_arm_only():
    offsets = pd.DataFrame(
        {"order": [1, 2, 3], "status": ["fitted", "failed", "failed"], "reason": ["", "gate", "too faint"]}
    )

    assert slit_illumination._failed_skyline_orders(offsets, "VIS") == [2, 3]
    assert slit_illumination._failed_skyline_orders(offsets, "NIR") == []
    assert np.array_equal(offsets["order"], [1, 2, 3])


def test_a_fallback_order_keeps_only_its_continuum_profile_when_the_skyline_profile_was_built(
    tmp_path, log, monkeypatch
):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 0.8}, slitName="SLIT1.0")

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"max-offset-arcsec": 0.5})

    row = _offsets(result).loc[11]
    assert row["method"] == "continuum"
    assert "exceeds the maximum" in row["reason"]
    orderProfile = result.profiles.loc[result.profiles["order"] == 11]
    assert len(orderProfile) > 20
    assert orderProfile["slit"].is_unique
