"""Contracts for the offset error rule, malformed inputs, partial flat-edge rows and the result's frame copies."""

from __future__ import annotations

import warnings
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils import slit_illumination, slit_illumination_profiles
from soxspipe.commonutils.slit_illumination import fit_order_offsets, slit_illumination_result
from tests.factories import synthetic_illumination_scene
from tests.factories.slit_illumination import top_hat, write_map_images
from tests.unit.test_slit_illumination_offset import _flat_edges, _measure, _offsets, _profile

pytestmark = pytest.mark.unit

LOOSE_MAX_ERR = 10.0


def _propagated_error(profile, delta, width, err):
    """The 1-sigma error on delta from the Jacobian of the shifted top-hat at the solution, for fixed ``err``."""
    slit = profile["slit"].to_numpy()
    start = np.array([1.0, delta, width])

    def model(p):
        return top_hat(slit, p[0], -5.0 + p[1], 5.0 + p[1], p[2])

    step = 1e-6
    jacobian = np.column_stack([(model(start + step * e) - model(start - step * e)) / (2 * step) for e in np.eye(3)])
    return float(np.sqrt(np.linalg.inv(jacobian.T @ jacobian / err**2)[1, 1]))


def _fit(profile, **kwargs):
    """The single offsets row of one order fitted with a loose error limit."""
    kwargs.setdefault("maxOffsetErrArcsec", LOOSE_MAX_ERR)
    return fit_order_offsets(profile, _flat_edges((11,)), **kwargs).iloc[0]


def test_offset_error_of_a_perfect_fit_is_the_propagated_error_of_the_fixed_input_errors():
    profile = _profile(11, 1.0, err=0.02)

    row = _fit(profile)

    expected = _propagated_error(profile, row["offset"], row["scienceWidth"], 0.02)
    assert row["offsetErr"] == pytest.approx(expected, rel=0.02)
    assert row["offsetErr"] == pytest.approx(0.00465, rel=0.02)


def test_doubling_the_input_errors_doubles_the_offset_error_when_the_fit_is_good():
    small = _fit(_profile(11, 1.0, err=0.02))
    large = _fit(_profile(11, 1.0, err=0.04))

    assert large["offsetErr"] == pytest.approx(2.0 * small["offsetErr"], rel=0.01)


def test_offset_error_is_inflated_by_root_reduced_chi_squared_when_the_fit_is_poor():
    # THE SCATTER IS TEN TIMES THE QUOTED ERROR, SO CHI-SQUARED PER DEGREE OF FREEDOM IS ABOUT 100
    profile = _profile(11, 1.0, noise=0.2, err=0.02, seed=4)

    row = _fit(profile)

    floor = _propagated_error(profile, row["offset"], row["scienceWidth"], 0.02)
    assert 8.0 < row["offsetErr"] / floor < 12.0


def test_offset_error_does_not_depend_on_the_scale_of_the_input_errors_when_the_fit_is_poor():
    noisy = {"noise": 0.2, "seed": 4}

    first = _fit(_profile(11, 1.0, err=0.02, **noisy))
    second = _fit(_profile(11, 1.0, err=0.04, **noisy))

    assert second["offsetErr"] == pytest.approx(first["offsetErr"], rel=0.02)


def test_offset_error_is_never_below_the_propagated_error_when_the_scatter_is_small():
    # THE SCATTER IS A TENTH OF THE QUOTED ERROR, SO A SCALED ERROR WOULD FALL BELOW THE PROPAGATED FLOOR
    profile = _profile(11, 1.0, noise=0.002, err=0.02, seed=4)

    row = _fit(profile)

    floor = _propagated_error(profile, row["offset"], row["scienceWidth"], 0.02)
    assert row["offsetErr"] == pytest.approx(floor, rel=0.02)


def test_the_offset_error_limit_accepts_an_order_just_under_it_and_rejects_one_just_over_it():
    profile = _profile(11, 1.0, err=0.02)
    measured = _fit(profile)["offsetErr"]

    accepted = _fit(profile, maxOffsetErrArcsec=measured * 1.001)
    rejected = _fit(profile, maxOffsetErrArcsec=measured * 0.999)

    assert accepted["status"] == "fitted"
    assert rejected["status"] == "failed"
    assert "uncertainty" in rejected["reason"]


def test_the_default_error_limit_accepts_a_good_fit():
    row = fit_order_offsets(_profile(11, 1.0, noise=0.03, err=0.03, seed=4), _flat_edges((11,))).iloc[0]

    assert row["status"] == "fitted"
    assert row["offsetErr"] < slit_illumination.MAX_OFFSET_ERR_ARCSEC


@pytest.mark.parametrize(
    ("column", "value"),
    [("lower", np.nan), ("upper", np.nan), ("width", np.nan), ("width", 0.0), ("width", -0.1)],
)
def test_fails_an_order_whose_flat_edge_row_is_partial_or_not_positive(column, value):
    flatEdges = _flat_edges((11,))
    flatEdges.loc[0, column] = value

    row = fit_order_offsets(_profile(11, 1.0), flatEdges).iloc[0]

    assert row["status"] == "failed"
    assert row["reason"].startswith("no flat edges")


def test_fails_an_order_whose_flat_edges_are_not_ordered():
    flatEdges = _flat_edges((11,), lower=5.0, upper=-5.0)

    row = fit_order_offsets(_profile(11, 1.0), flatEdges).iloc[0]

    assert row["status"] == "failed"
    assert row["reason"].startswith("no flat edges")


def test_flat_edge_fit_does_not_emit_a_row_with_a_non_finite_width(monkeypatch):
    monkeypatch.setattr(slit_illumination, "curve_fit", lambda *a, **k: (np.array([1.0, -5.0, 5.0, np.nan]), None))
    flatData = np.ones((200, 60))
    slit = np.tile(np.linspace(-5.4, 5.4, 200)[:, None], (1, 60))
    orderMap = np.full(slit.shape, 11.0)

    edges = slit_illumination.flat_slit_edges(flatData, slit, orderMap)

    assert edges["width"].isna().all()
    assert edges["reason"].str.contains("not finite").all()


@pytest.fixture
def scene(tmp_path):
    """Two orders with offsets +1 and -1 and a valid 2D map file."""
    return synthetic_illumination_scene(tmp_path, offsets={11: 1.0, 12: -1.0})


def _scene_with_map(scene, path):
    """The scene pointing at another 2D map file."""
    return replace(scene, twoDMapPath=path)


def test_is_not_applied_and_does_not_raise_when_the_map_file_is_missing(scene, tmp_path, log, monkeypatch):
    broken = _scene_with_map(scene, tmp_path / "no" / "such" / "map.fits")

    result = _measure(broken, tmp_path, log, monkeypatch)

    assert not result.isApplied
    assert "2D map" in result.reason
    assert any(level == "warning" and "2D map" in message for level, message in log.messages)


def test_is_not_applied_and_does_not_raise_when_the_map_lacks_an_extension(scene, tmp_path, log, monkeypatch):
    path = tmp_path / "partial" / "map.fits"
    path.parent.mkdir()
    fits.HDUList(
        [
            fits.PrimaryHDU(scene.wavelength.astype("float32"), header=fits.Header({"EXTNAME": "WAVELENGTH"})),
            fits.ImageHDU(scene.slit.astype("float32"), name="SLIT"),
        ]
    ).writeto(path)

    result = _measure(_scene_with_map(scene, path), tmp_path, log, monkeypatch)

    assert not result.isApplied
    assert "2D map" in result.reason
    assert "ORDER" in result.reason


@pytest.mark.parametrize("badImage", ["wavelength", "slit", "order"])
def test_is_not_applied_when_a_map_image_has_a_different_shape_from_the_frame(
    scene, tmp_path, log, monkeypatch, badImage
):
    images = {"wavelength": scene.wavelength, "slit": scene.slit, "order": scene.orderMap}
    images[badImage] = images[badImage][:100]
    path = write_map_images(tmp_path / "short" / "map.fits", images["wavelength"], images["slit"], images["order"])

    result = _measure(_scene_with_map(scene, path), tmp_path, log, monkeypatch)

    assert not result.isApplied
    assert "shape" in result.reason
    assert badImage.upper() in result.reason.upper()


def test_checks_the_frame_and_flat_before_it_opens_the_map(scene, tmp_path, log, monkeypatch):
    broken = _scene_with_map(scene, tmp_path / "no" / "such" / "map.fits")

    result = _measure(broken, tmp_path, log, monkeypatch, masterFlat=None)

    assert result.reason == "no master flat"


def test_the_frame_and_flat_shape_check_does_not_need_the_map(scene, tmp_path, log, monkeypatch):
    broken = _scene_with_map(scene, tmp_path / "no" / "such" / "map.fits")
    shortFlat = replace_flat(scene, scene.masterFlat.data[:100])

    result = _measure(broken, tmp_path, log, monkeypatch, masterFlat=shortFlat)

    assert not result.isApplied
    assert "shape" in result.reason
    assert "2D map" not in result.reason


def replace_flat(scene, data):
    """A master flat CCDData holding ``data`` in the scene flat's unit."""
    from astropy.nddata import CCDData

    return CCDData(data, unit=scene.masterFlat.unit)


def _result(offsets, flatEdges, profiles):
    """A result holding the given frames."""
    return slit_illumination_result(offsets, flatEdges, profiles, True, "", 0.3)


def test_the_result_holds_copies_of_its_frames_so_a_later_edit_does_not_reach_it():
    offsets = pd.DataFrame({"order": [11], "offset": [1.0]})
    flatEdges = pd.DataFrame({"order": [11], "lower": [-5.0]})
    profiles = pd.DataFrame({"order": [11], "slit": [0.0]})
    result = _result(offsets, flatEdges, profiles)

    offsets.loc[0, "offset"] = 99.0
    flatEdges.loc[0, "lower"] = 99.0
    profiles.loc[0, "slit"] = 99.0

    assert result.offsets.loc[0, "offset"] == 1.0
    assert result.flatEdges.loc[0, "lower"] == -5.0
    assert result.profiles.loc[0, "slit"] == 0.0
    assert result.offsets is not offsets


def test_the_result_still_refuses_to_have_its_fields_rebound():
    result = _result(pd.DataFrame(), pd.DataFrame(), pd.DataFrame())

    with pytest.raises(AttributeError):
        result.offsets = pd.DataFrame()


def test_a_measured_result_does_not_alias_the_offsets_it_was_built_from(scene, tmp_path, log, monkeypatch):
    first = _measure(scene, tmp_path, log, monkeypatch)
    before = _offsets(first).loc[11, "offset"]

    first.offsets.loc[:, "offset"] = 0.0
    second = _measure(scene, tmp_path, log, monkeypatch)

    assert _offsets(second).loc[11, "offset"] == pytest.approx(before)
    assert before == pytest.approx(1.0, abs=0.05)


def test_fails_an_order_that_has_profile_bins_but_no_row_in_the_flat_edge_table():
    row = fit_order_offsets(_profile(11, 1.0), _flat_edges((12,))).set_index("order").loc[11]

    assert row["status"] == "failed"
    assert row["reason"] == "no flat edges: the order is missing"


def test_is_not_applied_and_does_not_raise_when_a_map_extension_holds_no_data(scene, tmp_path, log, monkeypatch):
    path = tmp_path / "empty" / "map.fits"
    path.parent.mkdir()
    fits.HDUList(
        [
            fits.PrimaryHDU(scene.wavelength.astype("float32"), header=fits.Header({"EXTNAME": "WAVELENGTH"})),
            fits.ImageHDU(scene.slit.astype("float32"), name="SLIT"),
            fits.ImageHDU(name="ORDER"),
        ]
    ).writeto(path)

    result = _measure(_scene_with_map(scene, path), tmp_path, log, monkeypatch)

    assert not result.isApplied
    assert "2D map could not be read" in result.reason


def test_a_result_compares_by_identity_and_is_hashable_because_it_holds_dataframes():
    first = _result(pd.DataFrame({"a": [1]}), pd.DataFrame(), pd.DataFrame())
    second = _result(pd.DataFrame({"a": [1]}), pd.DataFrame(), pd.DataFrame())

    assert first == first
    assert first != second
    assert isinstance(hash(first), int)


def test_the_mask_of_the_frame_must_match_the_frame_shape(scene, tmp_path, log, monkeypatch):
    # CCDDATA REFUSES A MISMATCHED MASK ITSELF, SO ONLY A FRAME-LIKE OBJECT CAN REACH THIS CHECK
    shortMask = SimpleNamespace(
        data=scene.frame.data, unit=scene.frame.unit, mask=scene.frame.mask[:100], uncertainty=scene.frame.uncertainty
    )

    result = _measure(scene, tmp_path, log, monkeypatch, frame=shortMask)

    assert not result.isApplied
    assert "frame mask" in result.reason
    assert "shapes differ" in result.reason


def test_a_scalar_frame_mask_is_accepted(scene, tmp_path, log, monkeypatch):
    from astropy.nddata import CCDData

    scalarMask = CCDData(scene.frame.data, unit=scene.frame.unit, mask=np.False_, uncertainty=scene.frame.uncertainty)

    result = _measure(scene, tmp_path, log, monkeypatch, frame=scalarMask)

    assert "mask" not in result.reason


def test_the_shared_helpers_and_constants_are_public():
    for name in ("in_order_pixels", "normalised_flat_profile"):
        assert callable(getattr(slit_illumination_profiles, name))
    for name in (
        "OFFSET_COLUMNS",
        "FLAT_EDGE_COLUMNS",
        "STATUS_FITTED",
        "STATUS_BORROWED",
        "STATUS_FAILED",
        "PROFILE_COLUMNS",
        "SKYLINE_MIN_FLUX_VIS",
    ):
        assert name in slit_illumination.__all__
        assert hasattr(slit_illumination, name)
    assert slit_illumination.FLAT_EDGE_COLUMNS == ["order", "lower", "upper", "width", "reason"]


def test_an_empty_flat_edge_table_has_the_shared_column_list(scene, tmp_path, log, monkeypatch):
    result = _measure(scene, tmp_path, log, monkeypatch, masterFlat=None)

    assert list(result.flatEdges.columns) == slit_illumination.FLAT_EDGE_COLUMNS


def test_skips_an_offset_fit_start_whose_chi_squared_is_not_finite(monkeypatch):
    realFit = slit_illumination.curve_fit
    calls = []

    def first_start_is_nan(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            return np.array([np.nan, 0.0, 0.15]), np.eye(3)
        return realFit(*args, **kwargs)

    monkeypatch.setattr(slit_illumination, "curve_fit", first_start_is_nan)

    row = _fit(_profile(11, 1.0))

    assert row["status"] == "fitted"
    assert row["offset"] == pytest.approx(1.0, abs=0.01)


def test_fails_an_order_when_every_offset_fit_start_has_a_non_finite_chi_squared(monkeypatch):
    monkeypatch.setattr(slit_illumination, "curve_fit", lambda *a, **k: (np.array([np.nan, 0.0, 0.15]), np.eye(3)))

    row = _fit(_profile(11, 1.0))

    assert row["status"] == "failed"
    assert "did not converge" in row["reason"]


@pytest.mark.parametrize("badReason", [np.nan, None, 5, ""])
def test_the_flat_edge_reason_falls_back_when_it_is_not_a_usable_string(badReason):
    flatEdges = _flat_edges((11,))
    flatEdges["reason"] = flatEdges["reason"].astype(object)
    flatEdges.loc[0, ["lower", "upper", "width"]] = np.nan
    flatEdges.loc[0, "reason"] = badReason

    row = fit_order_offsets(_profile(11, 1.0), flatEdges).iloc[0]

    assert row["status"] == "failed"
    assert row["reason"] == "no flat edges: the flat edges are not finite, ordered and positive"


def test_a_real_flat_edge_reason_is_still_reported():
    flatEdges = _flat_edges((11,))
    flatEdges.loc[0, ["lower", "upper", "width"]] = np.nan
    flatEdges.loc[0, "reason"] = "flat level is not positive"

    row = fit_order_offsets(_profile(11, 1.0), flatEdges).iloc[0]

    assert row["reason"] == "no flat edges: flat level is not positive"


def test_the_flat_edge_fit_does_not_leak_an_optimize_warning(monkeypatch):
    from scipy.optimize import OptimizeWarning

    def warns_then_fits(*args, **kwargs):
        warnings.warn("Covariance of the parameters could not be estimated", OptimizeWarning, stacklevel=2)
        return np.array([1.0, -5.0, 5.0, 0.15]), None

    monkeypatch.setattr(slit_illumination, "curve_fit", warns_then_fits)
    flatData = np.ones((200, 60))
    slit = np.tile(np.linspace(-5.4, 5.4, 200)[:, None], (1, 60))

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        edges = slit_illumination.flat_slit_edges(flatData, slit, np.full(slit.shape, 11.0))

    assert edges.loc[0, "lower"] == pytest.approx(-5.0)
