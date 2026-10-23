"""Contracts for the continuum-edge offset of orders whose skyline fit has no lines to use (DY-1516)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import slit_illumination
from soxspipe.commonutils.slit_illumination import (
    OFFSET_COLUMNS,
    fit_continuum_offsets,
    fit_order_offsets,
    skyline_slit_profiles,
)
from soxspipe.commonutils.slit_illumination_profiles import (
    continuum_slit_profiles,
    has_no_skylines_note,
)
from tests.factories import synthetic_illumination_scene
from tests.unit.test_slit_illumination_offset import _flat_edges, _measure, _offsets, _profile

pytestmark = pytest.mark.unit

# A SKYLINE TABLE WHOSE ONLY LINE IS FAR REDWARD OF THE SCENE, AS THE VIS TABLE IS OF THE U AND G ORDERS
NO_LINES_TABLE = pd.DataFrame({"WAVELENGTH": [2000.0], "FLUX": [1.0]})
HALVES_LIMIT = 0.15


def _continuum_scene(tmp_path, *, delta=0.8, **kwargs):
    """A one-order scene with a continuum-only sky, a shifted top-hat illumination and no skyline in the table."""
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: delta}, lineAmplitude=0.0, **kwargs)
    return replace(scene, skylinesDF=NO_LINES_TABLE)


def _messages(log, level):
    """The messages the recording logger took at one level."""
    return [message for kind, message in log.messages if kind == level]


def _scene_profiles(scene, **kwargs):
    """Run ``continuum_slit_profiles`` on a scene with flat edges fitted from its master flat."""
    flatData = np.asarray(scene.masterFlat.data, dtype=float)
    flatEdges = slit_illumination.flat_slit_edges(flatData, scene.slit, scene.orderMap)
    usable = ~scene.frame.mask
    return continuum_slit_profiles(
        np.asarray(scene.frame.data, dtype=float),
        np.asarray(scene.frame.uncertainty.array, dtype=float),
        usable,
        flatData,
        scene.wavelength,
        scene.slit,
        scene.orderMap,
        flatEdges,
        [11],
        **kwargs,
    )


def test_an_order_without_skylines_gets_its_known_offset_from_the_continuum_edge(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8)

    result = _measure(scene, tmp_path, log, monkeypatch)

    row = _offsets(result).loc[11]
    assert row["status"] == "fitted"
    assert row["method"] == "continuum"
    assert row["offset"] == pytest.approx(0.8, abs=0.1)
    assert result.isApplied


def test_the_continuum_offset_has_the_right_sign_for_a_negative_shift(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=-0.8)

    row = _offsets(_measure(scene, tmp_path, log, monkeypatch)).loc[11]

    assert row["status"] == "fitted"
    assert row["offset"] == pytest.approx(-0.8, abs=0.1)


def test_a_bright_object_in_the_slit_is_masked_so_the_continuum_offset_still_holds(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8, objectSlit=1.0, objectCounts=400.0)

    row = _offsets(_measure(scene, tmp_path, log, monkeypatch)).loc[11]

    assert row["status"] == "fitted"
    assert row["offset"] == pytest.approx(0.8, abs=0.1)


def test_the_continuum_profile_leaves_out_the_bins_within_the_object_mask_half_width(tmp_path):
    scene = _continuum_scene(tmp_path, delta=0.8, objectSlit=1.0, objectCounts=400.0)

    profiles, _, notes = _scene_profiles(scene, objectMaskArcsec=2.5)

    assert notes == {}
    assert len(profiles) > 20
    assert (np.abs(profiles["slit"] - 1.0) > 2.5).all()
    assert list(profiles.columns) == slit_illumination.PROFILE_COLUMNS


def test_the_continuum_profile_has_at_least_the_minimum_chunks_in_every_bin(tmp_path):
    scene = _continuum_scene(tmp_path, delta=0.8)

    profiles, _, _ = _scene_profiles(scene, minChunksPerBin=5)

    assert (profiles["nLines"] >= 5).all()


def test_the_two_halves_are_built_from_the_blue_and_red_chunks_of_the_order(tmp_path):
    scene = _continuum_scene(tmp_path, delta=0.8)

    profiles, halves, _ = _scene_profiles(scene)

    blue, red = halves
    assert len(blue) > 20
    assert len(red) > 20
    assert (blue["nLines"] <= 10).all()
    assert (red["nLines"] <= 10).all()


def test_an_order_with_too_few_pixels_for_the_chunks_gets_a_note_and_no_profile(tmp_path):
    scene = _continuum_scene(tmp_path, delta=0.8)

    profiles, halves, notes = _scene_profiles(scene, nChunks=5000)

    assert profiles.empty
    assert halves[0].empty
    assert "too few usable pixels" in notes[11]


def test_a_rejected_continuum_order_falls_back_with_both_reasons_and_no_usable_offset(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8)

    result = _measure(
        scene, tmp_path, log, monkeypatch, illuminationSettings={"continuum-halves-max-diff-arcsec": 1e-9}
    )

    row = _offsets(result).loc[11]
    assert row["status"] == "failed"
    assert "halves" in row["reason"]
    assert "the skyline table has no lines between" in row["reason"]
    assert not result.isApplied


def test_the_continuum_settings_are_read_from_the_slit_illumination_block(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8)

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"continuum-min-chunks-per-bin": 100})

    row = _offsets(result).loc[11]
    assert row["status"] == "failed"
    assert row["method"] == "continuum"


def test_orders_with_skylines_are_measured_exactly_as_before_and_never_build_a_continuum_profile(
    tmp_path, log, monkeypatch
):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: -1.0})

    def forbidden(*args, **kwargs):
        raise AssertionError("the continuum builder must not run when every order has skylines")

    monkeypatch.setattr(slit_illumination, "continuum_slit_profiles", forbidden)
    offsets = _offsets(_measure(scene, tmp_path, log, monkeypatch))

    assert list(offsets["method"]) == ["skyline", "skyline"]
    assert offsets.loc[11, "offset"] == pytest.approx(1.0, abs=0.05)
    assert offsets.loc[12, "offset"] == pytest.approx(-1.0, abs=0.05)


def test_a_skyline_order_next_to_a_continuum_order_keeps_its_skyline_offset(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: 0.8})
    lines = scene.skylinesDF.iloc[:14]
    scene = replace(scene, skylinesDF=lines)
    skylineOnly = fit_order_offsets(*_skyline_inputs(scene, lines))

    offsets = _offsets(_measure(scene, tmp_path, log, monkeypatch))

    assert offsets.loc[11, "method"] == "skyline"
    assert offsets.loc[11, "offset"] == pytest.approx(skylineOnly.set_index("order").loc[11, "offset"])
    assert offsets.loc[12, "method"] == "continuum"


def _skyline_inputs(scene, lines):
    """The ``fit_order_offsets`` arguments the module builds for a scene's skylines."""
    flatData = np.asarray(scene.masterFlat.data, dtype=float)
    flatEdges = slit_illumination.flat_slit_edges(flatData, scene.slit, scene.orderMap)
    profiles, notes = skyline_slit_profiles(
        np.asarray(scene.frame.data, dtype=float),
        np.asarray(scene.frame.uncertainty.array, dtype=float),
        ~scene.frame.mask,
        flatData,
        scene.wavelength,
        scene.slit,
        scene.orderMap,
        lines["WAVELENGTH"].to_numpy(dtype=float),
    )
    return profiles, flatEdges


def test_a_skyline_offset_row_records_the_skyline_method():
    row = fit_order_offsets(_profile(11, 1.0), _flat_edges((11,))).iloc[0]

    assert row["method"] == "skyline"
    assert OFFSET_COLUMNS[-1] == "method"


def test_the_no_skylines_note_is_recognised_and_other_notes_are_not(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, hasNoise=False)
    flatData = np.asarray(scene.masterFlat.data, dtype=float)
    _, notes = skyline_slit_profiles(
        np.asarray(scene.frame.data, dtype=float),
        np.asarray(scene.frame.uncertainty.array, dtype=float),
        ~scene.frame.mask,
        flatData,
        scene.wavelength,
        scene.slit,
        scene.orderMap,
        np.array([2000.0]),
    )

    assert has_no_skylines_note(notes[11])
    assert not has_no_skylines_note("every skyline in the order is too faint (under 5 sigma per slit bin)")
    assert not has_no_skylines_note("")


def _halves(blueDelta, redDelta):
    """Blue-half and red-half profiles of order 11 with the given offsets."""
    return _profile(11, blueDelta), _profile(11, redDelta)


def test_continuum_offsets_are_accepted_when_the_two_halves_agree():
    rows = fit_continuum_offsets(
        _profile(11, 0.8), _halves(0.78, 0.84), _flat_edges((11,)), [11], 3.0, 0.3, HALVES_LIMIT
    )

    row = rows.iloc[0]
    assert row["status"] == "fitted"
    assert row["method"] == "continuum"
    assert row["offset"] == pytest.approx(0.8, abs=0.01)
    assert row["reason"] == ""


def test_continuum_offsets_are_rejected_with_the_halves_reason_when_the_halves_disagree():
    rows = fit_continuum_offsets(_profile(11, 0.9), _halves(0.4, 1.4), _flat_edges((11,)), [11], 3.0, 0.3, HALVES_LIMIT)

    row = rows.iloc[0]
    assert row["status"] == "failed"
    assert "halves" in row["reason"]
    assert "+0.40" in row["reason"]
    assert "+1.40" in row["reason"]
    assert "0.15" in row["reason"]


def test_continuum_offsets_within_the_halves_limit_are_accepted():
    rows = fit_continuum_offsets(
        _profile(11, 0.8), _halves(0.75, 0.85), _flat_edges((11,)), [11], 3.0, 0.3, HALVES_LIMIT
    )

    assert rows.iloc[0]["status"] == "fitted"


def test_continuum_offsets_are_rejected_when_a_half_cannot_be_fitted():
    emptyHalf = _profile(11, 0.8).iloc[:0]

    rows = fit_continuum_offsets(
        _profile(11, 0.8), (_profile(11, 0.8), emptyHalf), _flat_edges((11,)), [11], 3.0, 0.3, HALVES_LIMIT
    )

    row = rows.iloc[0]
    assert row["status"] == "failed"
    assert "halves" in row["reason"]
    assert "red half failed" in row["reason"]


def test_continuum_offsets_that_fail_the_production_gates_keep_the_gate_reason():
    rows = fit_continuum_offsets(
        _profile(11, 0.8, noise=0.5, err=0.5, seed=2), _halves(0.8, 0.8), _flat_edges((11,)), [11], 3.0, 0.001, 1.0
    )

    row = rows.iloc[0]
    assert row["status"] == "failed"
    assert "uncertainty" in row["reason"]


def test_continuum_offsets_cover_only_the_requested_orders():
    flatEdges = _flat_edges((11, 12))
    profiles = pd.concat([_profile(11, 0.8), _profile(12, 0.5)], ignore_index=True)
    blue = pd.concat([_profile(11, 0.8), _profile(12, 0.5)], ignore_index=True)

    rows = fit_continuum_offsets(profiles, (blue, blue), flatEdges, [12], 3.0, 0.3, HALVES_LIMIT)

    assert list(rows["order"]) == [12]


@pytest.mark.parametrize("slitName", ["SLIT5.0", "SLIT1.5", None])
def test_other_slits_never_build_a_continuum_profile_and_fall_back_as_before(tmp_path, log, monkeypatch, slitName):
    scene = _continuum_scene(tmp_path, delta=0.8, slitName=slitName)

    def forbidden(*args, **kwargs):
        raise AssertionError("the continuum builder must run on SLIT1.0 only")

    monkeypatch.setattr(slit_illumination, "continuum_slit_profiles", forbidden)
    result = _measure(scene, tmp_path, log, monkeypatch)

    row = _offsets(result).loc[11]
    assert row["status"] == "failed"
    assert row["method"] == "skyline"
    assert row["reason"].startswith("the skyline table has no lines between")
    assert not result.isApplied
    assert any("SLIT1.0" in message for message in _messages(log, "info"))


def test_the_slit_gate_message_names_the_frame_slit(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8, slitName="SLIT5.0")

    _measure(scene, tmp_path, log, monkeypatch)

    assert any("SLIT5.0" in message and "continuum" in message for message in _messages(log, "info"))


def test_a_continuum_offset_is_never_borrowed_by_a_partner_order():
    rows = [
        {**_profile_row(11, "fitted", "continuum", 0.8)},
        {**_profile_row(12, "failed", "skyline", np.nan)},
    ]

    result = slit_illumination._borrow_from_partners(rows, {12: 11})

    assert result[1]["status"] == "failed"
    assert pd.isna(result[1]["borrowedFrom"])


def test_a_skyline_offset_is_still_borrowed_by_a_partner_order():
    rows = [
        {**_profile_row(11, "fitted", "skyline", 0.8)},
        {**_profile_row(12, "failed", "skyline", np.nan)},
    ]

    result = slit_illumination._borrow_from_partners(rows, {12: 11})

    assert result[1]["status"] == "borrowed"
    assert result[1]["offset"] == 0.8


def _profile_row(order, status, method, offset):
    """An offsets-table row as a dictionary."""
    return {
        "order": order,
        "offset": offset,
        "offsetErr": 0.02,
        "scienceWidth": 0.2,
        "nBins": 100,
        "status": status,
        "reason": "" if status == "fitted" else "no lines",
        "borrowedFrom": np.nan,
        "method": method,
    }


@pytest.mark.parametrize("key", ["continuum-chunks", "continuum-min-chunks-per-bin"])
@pytest.mark.parametrize("value", [0.5, 1, 0, -3, 2.5, True, "20", None, float("nan")])
def test_unusable_continuum_count_settings_warn_and_use_the_default(tmp_path, log, monkeypatch, key, value):
    scene = _continuum_scene(tmp_path, delta=0.8)

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={key: value})

    assert any(key in message for message in _messages(log, "warning"))
    assert _offsets(result).loc[11, "status"] == "fitted"


def test_a_whole_number_float_is_a_usable_continuum_count(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8)

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"continuum-chunks": 20.0})

    assert not any("continuum-chunks" in message for message in _messages(log, "warning"))
    assert _offsets(result).loc[11, "status"] == "fitted"


def test_halves_that_cannot_reach_the_minimum_chunks_per_bin_are_warned_about(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8)

    result = _measure(
        scene,
        tmp_path,
        log,
        monkeypatch,
        illuminationSettings={"continuum-chunks": 8, "continuum-min-chunks-per-bin": 5},
    )

    assert any("halves" in message and "continuum-chunks" in message for message in _messages(log, "warning"))
    assert _offsets(result).loc[11, "status"] == "failed"


def test_the_default_continuum_settings_give_no_halves_warning(tmp_path, log, monkeypatch):
    scene = _continuum_scene(tmp_path, delta=0.8)

    _measure(scene, tmp_path, log, monkeypatch)

    assert not _messages(log, "warning")
