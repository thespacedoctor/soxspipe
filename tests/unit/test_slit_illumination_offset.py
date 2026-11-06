"""Contracts for the per-order slit offset: the fit, its fallbacks, and the measurement from skylines."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.nddata import CCDData, VarianceUncertainty

from soxspipe.commonutils import slit_illumination, toolkit
from soxspipe.commonutils.slit_illumination import (
    fit_order_offsets,
    measure_slit_illumination,
    skyline_slit_profiles,
)
from tests.factories import pipeline_settings, synthetic_illumination_scene
from tests.factories.slit_illumination import top_hat

pytestmark = pytest.mark.unit

FLAT_LOWER = -5.0
FLAT_UPPER = 5.0
FLAT_WIDTH = 0.15


def _flat_edges(orders=(11, 12), lower=FLAT_LOWER, upper=FLAT_UPPER, width=FLAT_WIDTH):
    """Build the flat-edge table that ``flat_slit_edges`` returns, with the same edges in every order."""
    return pd.DataFrame(
        {
            "order": list(orders),
            "lower": lower,
            "upper": upper,
            "width": width,
            "reason": "",
        }
    )


def _profile(order, delta, *, scienceWidth=FLAT_WIDTH, isLitOnly=True, noise=0.0, err=0.02, footprint=5.4, seed=0):
    """Build one order's science slit profile: the flat's top-hat shifted by ``delta``.

    With ``isLitOnly`` the bins outside the flat's lit range are dropped, as the flat-fielded frame is masked there.
    """
    centres = (np.arange(-round(footprint * 10), round(footprint * 10)) + 0.5) * 0.1
    if isLitOnly:
        centres = centres[(centres > FLAT_LOWER) & (centres < FLAT_UPPER)]
    values = top_hat(centres, 1.0, FLAT_LOWER + delta, FLAT_UPPER + delta, scienceWidth)
    values = values + np.random.default_rng(seed).normal(0.0, noise, centres.size)
    return pd.DataFrame(
        {"order": order, "slit": centres, "profile": values, "profileErr": err, "nLines": 14},
    )


@pytest.mark.parametrize("delta", [1.0, -1.0, 0.0, 0.4])
def test_recovers_the_offset_and_its_sign_from_a_profile_with_one_visible_edge(delta):
    offsets = fit_order_offsets(_profile(11, delta), _flat_edges((11,)))

    row = offsets.iloc[0]
    assert row["status"] == "fitted"
    assert row["offset"] == pytest.approx(delta, abs=0.01)
    assert row["scienceWidth"] == pytest.approx(FLAT_WIDTH, abs=0.01)
    assert row["nBins"] == len(_profile(11, delta))
    assert row["reason"] == ""


def test_recovers_opposite_offsets_in_two_orders_of_one_frame():
    profiles = pd.concat([_profile(11, -1.0), _profile(12, 1.0)], ignore_index=True)

    offsets = fit_order_offsets(profiles, _flat_edges()).set_index("order")

    assert offsets.loc[11, "offset"] == pytest.approx(-1.0, abs=0.01)
    assert offsets.loc[12, "offset"] == pytest.approx(1.0, abs=0.01)


def test_recovers_the_offset_when_a_science_edge_lies_outside_the_map_footprint():
    # THE LOWER SCIENCE EDGE IS AT -6.5, BEYOND THE -5.4 FOOTPRINT, SO ONLY THE UPPER EDGE IS IN THE DATA
    offsets = fit_order_offsets(_profile(11, -1.5, isLitOnly=False), _flat_edges((11,)))

    assert offsets.iloc[0]["status"] == "fitted"
    assert offsets.iloc[0]["offset"] == pytest.approx(-1.5, abs=0.01)


def test_recovers_the_offset_from_a_noisy_profile_and_reports_a_sensible_error():
    offsets = fit_order_offsets(_profile(11, 1.0, noise=0.03, err=0.03, seed=4), _flat_edges((11,)))

    row = offsets.iloc[0]
    assert row["offset"] == pytest.approx(1.0, abs=0.05)
    assert 0.0 < row["offsetErr"] < 0.05


def test_does_not_modify_its_inputs():
    profiles = _profile(11, 1.0)
    flatEdges = _flat_edges((11,))
    profilesBefore, flatEdgesBefore = profiles.copy(), flatEdges.copy()

    fit_order_offsets(profiles, flatEdges)

    pd.testing.assert_frame_equal(profiles, profilesBefore)
    pd.testing.assert_frame_equal(flatEdges, flatEdgesBefore)


def test_fails_an_order_with_too_few_profile_bins_and_says_so():
    profiles = _profile(11, 1.0).iloc[:10]

    row = fit_order_offsets(profiles, _flat_edges((11,))).iloc[0]

    assert row["status"] == "failed"
    assert np.isnan(row["offset"])
    assert "too few" in row["reason"]
    assert row["nBins"] == 10


def test_fails_an_order_that_has_no_profile_at_all():
    row = fit_order_offsets(_profile(11, 1.0).iloc[:0], _flat_edges((11,))).iloc[0]

    assert row["status"] == "failed"
    assert row["nBins"] == 0


def test_fails_an_order_without_a_profile_and_gives_the_profile_note_as_the_reason():
    empty = _profile(11, 1.0).iloc[:0]

    row = fit_order_offsets(empty, _flat_edges((11,)), profileNotes={11: "the skyline table has no lines here"}).iloc[0]

    assert row["status"] == "failed"
    assert "the skyline table has no lines here" in row["reason"]


def test_fails_an_order_without_flat_edges_and_carries_the_flat_reason():
    flatEdges = _flat_edges((11,))
    flatEdges.loc[0, ["lower", "upper", "width"]] = np.nan
    flatEdges.loc[0, "reason"] = "flat level is not positive"

    row = fit_order_offsets(_profile(11, 1.0), flatEdges).iloc[0]

    assert row["status"] == "failed"
    assert "no flat edges" in row["reason"]
    assert "flat level is not positive" in row["reason"]


def test_fails_an_order_whose_offset_exceeds_the_maximum():
    row = fit_order_offsets(_profile(11, 3.5), _flat_edges((11,))).iloc[0]

    assert row["status"] == "failed"
    assert "exceeds the maximum" in row["reason"]
    assert row["offset"] == pytest.approx(3.5, abs=0.05)


def test_accepts_a_larger_offset_when_the_maximum_is_raised():
    row = fit_order_offsets(_profile(11, 3.5), _flat_edges((11,)), maxOffsetArcsec=4.0).iloc[0]

    assert row["status"] == "fitted"


@pytest.mark.parametrize("scienceWidth", [2.0, 0.01])
def test_fails_an_order_whose_science_edge_width_is_degenerate(scienceWidth):
    row = fit_order_offsets(_profile(11, 1.0, scienceWidth=scienceWidth), _flat_edges((11,))).iloc[0]

    assert row["status"] == "failed"
    assert "edge width" in row["reason"]


def test_fails_an_order_whose_offset_error_is_too_large():
    profiles = _profile(11, 1.0, noise=0.2, err=0.2, seed=2)

    row = fit_order_offsets(profiles, _flat_edges((11,)), maxOffsetErrArcsec=0.001).iloc[0]

    assert row["status"] == "failed"
    assert "uncertainty" in row["reason"]
    assert row["offsetErr"] > 0.001


def test_fails_an_order_when_the_fit_does_not_converge(monkeypatch):
    def refuse(*args, **kwargs):
        raise RuntimeError("Optimal parameters not found")

    monkeypatch.setattr(slit_illumination, "curve_fit", refuse)

    row = fit_order_offsets(_profile(11, 1.0), _flat_edges((11,))).iloc[0]

    assert row["status"] == "failed"
    assert "did not converge" in row["reason"]


def test_fails_an_order_whose_offset_uncertainty_cannot_be_estimated(monkeypatch):
    singular = np.full((3, 3), np.inf)
    monkeypatch.setattr(slit_illumination, "curve_fit", lambda *args, **kwargs: (np.array([1.0, 0.5, 0.15]), singular))

    row = fit_order_offsets(_profile(11, 0.5), _flat_edges((11,))).iloc[0]

    assert row["status"] == "failed"
    assert "not finite" in row["reason"]


def test_borrows_the_partner_offset_for_a_failed_order_and_records_it():
    profiles = pd.concat([_profile(11, 1.0), _profile(12, 1.0).iloc[:5]], ignore_index=True)

    offsets = fit_order_offsets(profiles, _flat_edges(), partners={12: 11}).set_index("order")

    assert offsets.loc[11, "status"] == "fitted"
    assert offsets.loc[12, "status"] == "borrowed"
    assert offsets.loc[12, "offset"] == pytest.approx(offsets.loc[11, "offset"])
    assert offsets.loc[12, "offsetErr"] == pytest.approx(offsets.loc[11, "offsetErr"])
    assert offsets.loc[12, "borrowedFrom"] == 11
    assert "too few" in offsets.loc[12, "reason"]


def test_leaves_a_failed_order_failed_when_no_partner_mapping_is_given():
    profiles = pd.concat([_profile(11, 1.0), _profile(12, 1.0).iloc[:5]], ignore_index=True)

    offsets = fit_order_offsets(profiles, _flat_edges()).set_index("order")

    assert offsets.loc[12, "status"] == "failed"
    assert np.isnan(offsets.loc[12, "offset"])
    assert pd.isna(offsets.loc[12, "borrowedFrom"])


def test_leaves_a_failed_order_failed_when_its_partner_failed_too():
    profiles = pd.concat([_profile(11, 1.0).iloc[:5], _profile(12, 1.0).iloc[:5]], ignore_index=True)

    offsets = fit_order_offsets(profiles, _flat_edges(), partners={12: 11, 11: 12})

    assert (offsets["status"] == "failed").all()


def test_does_not_borrow_from_an_order_that_is_not_in_the_frame():
    offsets = fit_order_offsets(_profile(11, 1.0).iloc[:5], _flat_edges((11,)), partners={11: 99})

    assert offsets.iloc[0]["status"] == "failed"


def test_does_not_overwrite_a_fitted_order_with_its_partner():
    profiles = pd.concat([_profile(11, 1.0), _profile(12, -1.0)], ignore_index=True)

    offsets = fit_order_offsets(profiles, _flat_edges(), partners={12: 11}).set_index("order")

    assert offsets.loc[12, "status"] == "fitted"
    assert offsets.loc[12, "offset"] == pytest.approx(-1.0, abs=0.01)


def _scene_profiles_and_notes(scene, lineWavelengths=None, binArcsec=0.1, isUsable=None, flatData=None):
    """Run ``skyline_slit_profiles`` on a synthetic scene's frame, flat, map and skylines."""
    frame = scene.frame
    if isUsable is None:
        isUsable = ~frame.mask & np.isfinite(frame.data) & np.isfinite(frame.uncertainty.array)
    lines = scene.skylinesDF["WAVELENGTH"].to_numpy() if lineWavelengths is None else lineWavelengths
    return skyline_slit_profiles(
        frame.data,
        frame.uncertainty.array,
        isUsable,
        scene.masterFlat.data if flatData is None else flatData,
        scene.wavelength,
        scene.slit,
        scene.orderMap,
        lines,
        binArcsec=binArcsec,
    )


def _scene_profiles(scene, **kwargs):
    """The profile table alone."""
    return _scene_profiles_and_notes(scene, **kwargs)[0]


def _total_lines(profiles):
    """The number of line measurements behind a profile table, summed over bins and orders."""
    return int(profiles["nLines"].sum())


def _profile_at(profiles, order, slit):
    """The profile value of the bin nearest ``slit`` in one order."""
    orderRows = profiles.loc[profiles["order"] == order]
    return orderRows.iloc[(orderRows["slit"] - slit).abs().argmin()]["profile"]


def test_skyline_profile_is_the_flat_top_hat_shifted_by_the_offset(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: -1.0}, hasNoise=False)

    profiles = _scene_profiles(scene)

    # ORDER 11: LIT FROM -4 (FLAT LOWER EDGE + 1) UP TO THE END OF THE LIT FLAT AT +5
    assert _profile_at(profiles, 11, -4.6) == pytest.approx(0.0, abs=0.03)
    assert _profile_at(profiles, 11, -3.0) == pytest.approx(1.0, abs=0.03)
    assert _profile_at(profiles, 11, 4.0) == pytest.approx(1.0, abs=0.03)
    # ORDER 12: LIT UP TO +4, DARK BETWEEN +4 AND THE END OF THE LIT FLAT
    assert _profile_at(profiles, 12, 3.0) == pytest.approx(1.0, abs=0.03)
    assert _profile_at(profiles, 12, 4.6) == pytest.approx(0.0, abs=0.03)


def test_skyline_profile_reports_the_columns_and_the_number_of_lines_per_bin(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0}, hasNoise=False)

    profiles = _scene_profiles(scene)

    assert list(profiles.columns) == ["order", "slit", "profile", "profileErr", "nLines"]
    assert set(profiles["order"]) == {11, 12}
    assert (profiles["profileErr"] > 0).all()
    # A 0.1 ARCSEC BIN IS HALF A PIXEL, SO EACH LINE FILLS ONLY SOME BINS AND NO BIN HOLDS MORE THAN THE 14 LINES
    assert profiles["nLines"].between(3, 14).all()
    assert profiles["nLines"].max() >= 7


def test_skyline_profile_skips_lines_whose_side_bands_hold_a_neighbouring_line(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0}, hasNoise=False)
    lines = scene.skylinesDF["WAVELENGTH"].to_numpy()
    # A NEIGHBOUR 0.15 NM (SEVEN PIXELS) FROM THE FIRST LINE OF ORDER 11 SITS IN ITS SIDE BAND. IT HOLDS NO LIGHT
    crowded = np.append(lines, lines[0] + 0.15)

    allLines = _total_lines(_scene_profiles(scene))
    withoutFirstLine = _total_lines(_scene_profiles(scene, lineWavelengths=lines[1:]))
    withNeighbour = _total_lines(_scene_profiles(scene, lineWavelengths=crowded))

    assert withoutFirstLine < allLines
    assert withNeighbour == withoutFirstLine


def test_skyline_profile_ignores_lines_outside_every_order(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0}, hasNoise=False)
    lines = np.append(scene.skylinesDF["WAVELENGTH"].to_numpy(), [300.0, 1500.0])

    assert _total_lines(_scene_profiles(scene, lineWavelengths=lines)) == _total_lines(_scene_profiles(scene))


def test_skyline_profile_is_empty_when_the_wavelength_does_not_change_between_pixels(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, hasNoise=False)
    flatWavelength = np.where(np.isfinite(scene.wavelength), 600.0, np.nan)
    frame = scene.frame

    profiles, _ = skyline_slit_profiles(
        frame.data,
        frame.uncertainty.array,
        ~frame.mask,
        scene.masterFlat.data,
        flatWavelength,
        scene.slit,
        scene.orderMap,
        np.array([600.0]),
    )

    assert profiles.empty


def test_skyline_notes_say_when_the_skyline_table_has_no_lines_in_an_order(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, hasNoise=False)
    # THE TABLE COVERS ORDER 11 (500-512 NM) ONLY
    lines = scene.skylinesDF["WAVELENGTH"].to_numpy()[:14]

    profiles, notes = _scene_profiles_and_notes(scene, lineWavelengths=lines)

    assert set(profiles["order"]) == {11}
    assert list(notes) == [12]
    assert "no lines between 600 and 612 nm" in notes[12]


def test_skyline_notes_say_when_every_line_is_blended(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, hasNoise=False)
    lines = scene.skylinesDF["WAVELENGTH"].to_numpy()
    blended = np.concatenate([lines[:14], lines[14:], lines[:14] + 0.15])

    profiles, notes = _scene_profiles_and_notes(scene, lineWavelengths=blended)

    assert list(notes) == [11]
    assert "blended" in notes[11]
    assert set(profiles["order"]) == {12}


def test_skyline_notes_say_when_the_lines_are_too_faint(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), lineAmplitude=0.5)

    profiles, notes = _scene_profiles_and_notes(scene)

    assert profiles.empty
    assert "too faint" in notes[11]


def test_skyline_notes_say_when_too_few_lines_share_a_slit_bin(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), hasNoise=False)
    twoLines = scene.skylinesDF["WAVELENGTH"].to_numpy()[:2]

    profiles, notes = _scene_profiles_and_notes(scene, lineWavelengths=twoLines)

    assert profiles.empty
    assert "at least 3" in notes[11]


def test_skyline_notes_say_when_an_order_has_no_usable_pixels(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, hasNoise=False)
    isUsable = ~scene.frame.mask & (scene.orderMap != 12)

    profiles, notes = _scene_profiles_and_notes(scene, isUsable=isUsable)

    assert set(profiles["order"]) == {11}
    assert "no usable pixels" in notes[12]


def test_skyline_notes_carry_the_reason_the_flat_is_unusable(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, hasNoise=False)

    profiles, notes = _scene_profiles_and_notes(scene, flatData=np.zeros_like(scene.masterFlat.data))

    assert profiles.empty
    assert all("flat level is not positive" in note for note in notes.values())
    assert set(notes) == {11, 12}


def test_skyline_profile_is_empty_when_there_are_no_skylines(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, hasNoise=False)

    profiles = _scene_profiles(scene, lineWavelengths=np.array([]))

    assert profiles.empty
    assert list(profiles.columns) == ["order", "slit", "profile", "profileErr", "nLines"]


def test_skyline_profile_does_not_modify_its_inputs(tmp_path):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0})
    before = [array.copy() for array in (scene.frame.data, scene.frame.uncertainty.array, scene.slit, scene.wavelength)]

    profiles = _scene_profiles(scene)

    assert not profiles.empty

    after = (scene.frame.data, scene.frame.uncertainty.array, scene.slit, scene.wavelength)
    for expected, actual in zip(before, after, strict=True):
        np.testing.assert_array_equal(actual, expected)


def test_skyline_profile_is_unchanged_far_from_a_bright_continuum_object(tmp_path):
    clean = synthetic_illumination_scene(tmp_path / "clean", offsets={11: -1.0}, hasNoise=False)
    bright = synthetic_illumination_scene(
        tmp_path / "bright", offsets={11: -1.0}, hasNoise=False, objectSlit=3.5, objectCounts=3000.0
    )

    cleanProfile = _scene_profiles(clean).query("order == 11 and slit < 1.5").set_index("slit")["profile"]
    brightProfile = _scene_profiles(bright).query("order == 11 and slit < 1.5").set_index("slit")["profile"]

    assert len(cleanProfile) > 50
    assert (brightProfile - cleanProfile).abs().max() < 0.01


calls = []


def _measure(scene, tmp_path, log, monkeypatch, *, frame=None, masterFlat="scene", binx=1, biny=1, arm="VIS", **kwargs):
    """Run ``measure_slit_illumination`` on a scene with the skyline table patched to the scene's lines."""
    calls.clear()

    def skylines(log, settings, arm, **kwargs):
        calls.append(kwargs)
        return scene.skylinesDF

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", skylines)
    return measure_slit_illumination(
        log,
        scene.frame if frame is None else frame,
        scene.masterFlat if isinstance(masterFlat, str) else masterFlat,
        scene.twoDMapPath,
        pipeline_settings(tmp_path),
        arm,
        binx,
        biny,
        **kwargs,
    )


def _offsets(result):
    """The result's offsets table indexed by order."""
    return result.offsets.set_index("order")


@pytest.mark.parametrize("delta", [1.0, -1.0])
def test_measures_the_offset_and_pins_the_sign_convention(tmp_path, log, monkeypatch, delta):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: delta})

    result = _measure(scene, tmp_path, log, monkeypatch)

    assert result.isApplied
    assert result.reason == ""
    assert _offsets(result).loc[11, "offset"] == pytest.approx(delta, abs=0.05)
    assert _offsets(result).loc[11, "status"] == "fitted"


def test_measures_opposite_offsets_in_two_orders_of_one_frame(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: -1.0, 12: 1.0})

    offsets = _offsets(_measure(scene, tmp_path, log, monkeypatch))

    assert offsets.loc[11, "offset"] == pytest.approx(-1.0, abs=0.05)
    assert offsets.loc[12, "offset"] == pytest.approx(1.0, abs=0.05)


def test_measures_an_offset_whose_other_science_edge_lies_outside_the_map_footprint(tmp_path, log, monkeypatch):
    # THE LOWER SCIENCE EDGE IS AT -6.5 ARCSEC, BEYOND THE -5.4 ARCSEC FOOTPRINT
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: -1.5})

    offsets = _offsets(_measure(scene, tmp_path, log, monkeypatch))

    assert offsets.loc[11, "offset"] == pytest.approx(-1.5, abs=0.05)


def test_a_bright_continuum_object_near_the_lit_edge_biases_the_offset_by_less_than_0p05_arcsec(
    tmp_path, log, monkeypatch
):
    scene = synthetic_illumination_scene(
        tmp_path, orders=(11,), offsets={11: -1.0}, objectSlit=3.5, objectCounts=3000.0
    )

    offsets = _offsets(_measure(scene, tmp_path, log, monkeypatch))

    assert offsets.loc[11, "status"] == "fitted"
    assert offsets.loc[11, "offset"] == pytest.approx(-1.0, abs=0.05)


def test_asks_the_skyline_table_for_the_fainter_lines_that_the_blue_and_red_orders_need(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 1.0})

    _measure(scene, tmp_path, log, monkeypatch)

    assert calls == [{"minBrightnessVIS": slit_illumination.SKYLINE_MIN_FLUX_VIS}]


def test_says_why_an_order_without_skylines_in_the_table_fails(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: 1.0})
    scene = replace(scene, skylinesDF=scene.skylinesDF.iloc[:14])

    # THE HALVES LIMIT IS TINY, SO THE CONTINUUM-EDGE OFFSET (DY-1516) IS REJECTED AND THE ORDER FALLS BACK
    result = _measure(
        scene, tmp_path, log, monkeypatch, illuminationSettings={"continuum-halves-max-diff-arcsec": 1e-9}
    )

    assert result.isApplied
    assert _offsets(result).loc[12, "status"] == "failed"
    assert "no lines between 600 and 612 nm" in _offsets(result).loc[12, "reason"]


def test_does_not_modify_the_frame_flat_or_map(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: -1.0})
    frameData, frameMask = scene.frame.data.copy(), scene.frame.mask.copy()
    flatData = scene.masterFlat.data.copy()
    with fits.open(scene.twoDMapPath) as hdus:
        mapBefore = hdus["SLIT"].data.copy()

    _measure(scene, tmp_path, log, monkeypatch)

    np.testing.assert_array_equal(scene.frame.data, frameData)
    np.testing.assert_array_equal(scene.frame.mask, frameMask)
    np.testing.assert_array_equal(scene.masterFlat.data, flatData)
    with fits.open(scene.twoDMapPath) as hdus:
        np.testing.assert_array_equal(hdus["SLIT"].data, mapBefore)


def test_returns_the_flat_edges_and_profiles_it_measured(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: 1.0})

    result = _measure(scene, tmp_path, log, monkeypatch)

    assert list(result.flatEdges["order"]) == [11, 12]
    assert result.flatEdges["lower"].to_numpy() == pytest.approx([-5.0, -5.0], abs=0.03)
    assert set(result.profiles["order"]) == {11, 12}


def test_the_result_is_frozen(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 1.0})
    result = _measure(scene, tmp_path, log, monkeypatch)

    with pytest.raises(FrozenInstanceError):
        result.isApplied = False


def test_is_not_applied_without_a_master_flat_and_does_not_raise(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0})

    result = _measure(scene, tmp_path, log, monkeypatch, masterFlat=None)

    assert not result.isApplied
    assert "no master flat" in result.reason
    assert result.offsets.empty
    assert any(level == "info" and "no master flat" in message for level, message in log.messages)


def test_is_not_applied_when_the_flat_has_no_lit_level(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0})
    dark = CCDData(np.zeros_like(scene.masterFlat.data), unit=scene.masterFlat.unit)

    result = _measure(scene, tmp_path, log, monkeypatch, masterFlat=dark)

    assert not result.isApplied
    assert "flat level is not positive" in result.reason


def test_is_not_applied_when_the_frame_has_no_uncertainty(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0})
    bare = CCDData(scene.frame.data, unit=scene.frame.unit, mask=scene.frame.mask)

    result = _measure(scene, tmp_path, log, monkeypatch, frame=bare)

    assert not result.isApplied
    assert "uncertainty" in result.reason


def test_is_not_applied_when_the_frame_uncertainty_is_not_a_standard_deviation(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0})
    variance = VarianceUncertainty(scene.frame.uncertainty.array**2, unit=scene.frame.unit**2)
    odd = CCDData(scene.frame.data, unit=scene.frame.unit, mask=scene.frame.mask, uncertainty=variance)

    result = _measure(scene, tmp_path, log, monkeypatch, frame=odd)

    assert not result.isApplied
    assert "not std" in result.reason


def test_is_not_applied_when_the_frame_and_map_shapes_differ(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0})
    small = CCDData(
        scene.frame.data[:100], unit=scene.frame.unit, mask=scene.frame.mask[:100], uncertainty=scene.frame.uncertainty
    )

    result = _measure(scene, tmp_path, log, monkeypatch, frame=small)

    assert not result.isApplied
    assert "shape" in result.reason


def test_is_not_applied_when_every_order_fails_and_lists_the_reasons(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 3.5, 12: -3.5})

    result = _measure(scene, tmp_path, log, monkeypatch)

    assert not result.isApplied
    assert "11: " in result.reason and "12: " in result.reason
    assert "exceeds the maximum" in result.reason
    assert (result.offsets["status"] == "failed").all()
    assert any(level == "info" and "order 11" in message for level, message in log.messages)


def test_is_applied_when_at_least_one_order_has_an_offset(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: 3.5})

    result = _measure(scene, tmp_path, log, monkeypatch)

    assert result.isApplied
    assert _offsets(result).loc[12, "status"] == "failed"


def test_takes_the_maximum_offset_from_the_illumination_settings(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 3.5})

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"max-offset-arcsec": 4.0})

    assert result.isApplied
    assert _offsets(result).loc[11, "offset"] == pytest.approx(3.5, abs=0.05)


def test_takes_the_profile_bin_width_from_the_illumination_settings(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 1.0})

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"profile-bin-arcsec": 0.2})

    assert np.diff(result.profiles["slit"]).min() == pytest.approx(0.2)
    assert _offsets(result).loc[11, "offset"] == pytest.approx(1.0, abs=0.05)


def test_takes_the_largest_offset_error_from_the_illumination_settings(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 1.0})

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"max-order-scatter-arcsec": 1e-9})

    assert not result.isApplied
    assert "uncertainty" in result.reason


def test_carries_the_overlap_margin_for_the_mask_that_follows(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 1.0})

    default = _measure(scene, tmp_path, log, monkeypatch)
    custom = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"overlap-margin-arcsec": 0.5})

    assert default.overlapMarginArcsec == pytest.approx(0.45)
    assert custom.overlapMarginArcsec == pytest.approx(0.5)


@pytest.mark.parametrize("badValue", [-1.0, 0.0, "wide", None, float("nan")])
def test_falls_back_to_the_default_for_an_unusable_setting_and_logs_it(tmp_path, log, monkeypatch, badValue):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 3.5})

    result = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={"max-offset-arcsec": badValue})

    assert not result.isApplied
    assert any(level == "warning" and "max-offset-arcsec" in message for level, message in log.messages)


def test_a_missing_settings_block_means_the_module_constants(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 1.0})

    withNone = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings=None)
    withEmpty = _measure(scene, tmp_path, log, monkeypatch, illuminationSettings={})

    pd.testing.assert_frame_equal(withNone.offsets, withEmpty.offsets)
    assert withNone.overlapMarginArcsec == slit_illumination.OVERLAP_MARGIN_ARCSEC


def test_a_failed_order_borrows_its_partner_offset_when_a_mapping_is_given(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: 1.0})
    # NO USABLE PIXEL IN ORDER 12
    mask = scene.frame.mask | (scene.orderMap == 12)
    broken = CCDData(scene.frame.data, unit=scene.frame.unit, mask=mask, uncertainty=scene.frame.uncertainty)

    borrowed = _measure(scene, tmp_path, log, monkeypatch, frame=broken, partners={12: 11})
    alone = _measure(scene, tmp_path, log, monkeypatch, frame=broken)

    assert _offsets(borrowed).loc[12, "status"] == "borrowed"
    assert _offsets(borrowed).loc[12, "offset"] == pytest.approx(_offsets(borrowed).loc[11, "offset"])
    assert _offsets(borrowed).loc[12, "borrowedFrom"] == 11
    assert _offsets(alone).loc[12, "status"] == "failed"
    assert alone.isApplied


def test_measures_a_binned_frame_from_a_map_at_full_resolution(tmp_path, log, monkeypatch):
    scene = synthetic_illumination_scene(tmp_path, orders=(11,), offsets={11: 1.0})
    # THE MAP HOLDS EVERY ROW TWICE, SO ITS 1 X 2 BLOCK AVERAGE IS THE SCENE'S MAP
    from tests.factories.slit_illumination import write_map_images

    doubled = write_map_images(
        tmp_path / "doubled" / "map.fits",
        np.repeat(scene.wavelength, 2, axis=0),
        np.repeat(scene.slit, 2, axis=0),
        np.repeat(scene.orderMap, 2, axis=0),
    )
    binnedScene = replace(scene, twoDMapPath=doubled)

    result = _measure(binnedScene, tmp_path, log, monkeypatch, biny=2)

    assert _offsets(result).loc[11, "offset"] == pytest.approx(1.0, abs=0.05)
