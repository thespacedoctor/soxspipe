"""The sky fit takes the along-slit illumination from the master flat when ``flat-slit-gradient`` is on (DY-2147).

A synthetic NIR order is flat-fielded with a flat that is 10 % brighter at the top of the slit than at the bottom,
while the sky lights the slit uniformly. The flat-fielded sky then falls by 10 % along the slit. The DY-1280 acceptance
test: with the setting on, the sky model follows that tilt and the residual has no slope along the slit.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate
import yaml

import soxspipe
from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

ORDER = 11
COLUMNS = 800
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
KEPT_RANGE = (-5.0, 5.0)
START_WAVELENGTH = 1700.0
DISPERSION = 0.37
LINE_SIGMA_PX = 1.06
LINE_SPACING_PX = 20
LINE_PEAK = 3000.0
CONTINUUM = 400.0
READ_NOISE = 8.6
# THE FLAT RISES BY THIS FRACTION OF ITS LEVEL FROM THE BOTTOM TO THE TOP OF THE 11 ARCSEC SLIT
FLAT_TILT = 0.10
SETTING_KEY = "flat-slit-gradient"

DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"


def _sky_settings(**overrides: Any) -> dict[str, Any]:
    """The shipped SOXS NIR stare sky-subtraction settings, with the QC plot turned off."""
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    skySettings = copy.deepcopy(settings["soxs-stare"]["nir"]["sky-subtraction"])
    skySettings["sky_model_qc_plot"] = False
    skySettings.pop(SETTING_KEY, None)
    skySettings.update(overrides)
    return skySettings


def _bare_subtractor(log: Any, skySettings: dict[str, Any], keptRanges: dict | None = None) -> subtract_sky:
    """A subtract_sky with only the attributes the per-order steps read."""
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "NIR"
    subtractor.debug = False
    subtractor.binx = 1
    subtractor.biny = 1
    subtractor.axisA = "y"
    subtractor.axisB = "x"
    subtractor.ron = READ_NOISE
    subtractor.recipeSettings = {"sky-subtraction": skySettings}
    subtractor.bspline_order = skySettings["bspline_order"]
    subtractor.qcPlotOrder = -1
    if keptRanges is not None:
        subtractor.flatFitRanges = keptRanges
    return subtractor


def _tilted_order(seed: int = 2147, tilt: float = FLAT_TILT) -> pd.DataFrame:
    """One flat-fielded sky-only order: a uniformly lit sky divided by a flat that tilts along the slit."""
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slitPixel = slitPixel.ravel()
    slitPosition = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)[slitPixel]
    wavelength = START_WAVELENGTH + DISPERSION * column
    centres = np.arange(LINE_SPACING_PX / 2, COLUMNS, LINE_SPACING_PX)
    lines = LINE_PEAK * np.exp(-0.5 * ((column[:, None] - centres[None, :]) / LINE_SIGMA_PX) ** 2).sum(axis=1)
    flatProfile = 1.0 + tilt * slitPosition / (2 * SLIT_HALF_LENGTH)
    flat = 1000.0 * (1.0 + 0.2 * np.sin(column / 150.0)) * flatProfile
    noiseless = (CONTINUUM + lines) / flatProfile
    error = np.sqrt(noiseless + READ_NOISE**2)
    frame = pd.DataFrame(
        {
            "order": float(ORDER),
            "x": slitPixel,
            "y": column,
            "wavelength": wavelength,
            "slit_position": slitPosition,
            "flux": noiseless + rng.standard_normal(noiseless.size) * error,
            "error": error,
            "flat": flat,
            "mask": False,
        }
    )
    return frame.sort_values("wavelength", kind="stable").reset_index(drop=True)


def _clipped(subtractor: subtract_sky, order: pd.DataFrame) -> pd.DataFrame:
    """Run the pre-fit clipping of ``subtract`` on one order."""
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = subtractor.get_over_sampled_sky_from_order(order.copy(), clipBPs=True, clipSlitEdge=False)
    return subtractor.clip_object_slit_positions([clipped], aggressive=settings["aggressive_object_masking"])[0]


def _fit_all(subtractor: subtract_sky, clipped: pd.DataFrame) -> tuple:
    """Run ``fit_bspline_curve_to_sky`` on a clipped order and return everything it returns."""
    return subtractor.fit_bspline_curve_to_sky(clipped.copy())


def _fit(subtractor: subtract_sky, order: pd.DataFrame) -> pd.DataFrame:
    """Run the pre-fit clipping of ``subtract`` and then ``fit_bspline_curve_to_sky`` on one order."""
    return _fit_all(subtractor, _clipped(subtractor, order))[0]


def _rows(modelled: pd.DataFrame) -> pd.DataFrame:
    """The mean flux, sky model and residual at each slit position, over the kept range."""
    inRange = modelled["slit_position"].between(*KEPT_RANGE)
    return modelled[inRange].groupby("slit_position")[["flux", "sky_model", "sky_subtracted_flux"]].mean()


def _residual_slope_percent_per_arcsec(modelled: pd.DataFrame) -> float:
    """The slope along the slit of the mean residual, as a percentage of the mean sky model per arcsec."""
    rows = _rows(modelled)
    slope = np.polyfit(rows.index.to_numpy(), rows["sky_subtracted_flux"].to_numpy(), 1)[0]
    return 100 * slope / rows["sky_model"].mean()


def test_sky_model_follows_the_flat_tilt_and_the_residual_has_no_slope_with_the_setting_on(log: Any) -> None:
    # ARRANGE
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {ORDER: KEPT_RANGE})

    # ACT
    modelled = _fit(subtractor, _tilted_order())

    # ASSERT: THE MODEL FALLS BY ABOUT 9 % FROM THE BOTTOM TO THE TOP OF THE KEPT RANGE, AS THE DATA DO
    rows = _rows(modelled)
    modelTilt = rows["sky_model"].iloc[-1] / rows["sky_model"].iloc[0]
    dataTilt = rows["flux"].iloc[-1] / rows["flux"].iloc[0]
    bottom, top = (1.0 + FLAT_TILT * rows.index[[0, -1]] / (2 * SLIT_HALF_LENGTH)).to_numpy()
    assert dataTilt == pytest.approx(bottom / top, rel=0.01)
    assert modelTilt == pytest.approx(dataTilt, rel=0.003)
    assert abs(_residual_slope_percent_per_arcsec(modelled)) < 0.03
    assert modelled["slit_normalisation_ratio"].iloc[0] != 1


def test_residual_keeps_the_slope_with_the_setting_off(log: Any) -> None:
    # ARRANGE
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: False}), {ORDER: KEPT_RANGE})

    # ACT
    modelled = _fit(subtractor, _tilted_order())

    # ASSERT: THE UNMODELLED 10 % TILT IS ABOUT 0.9 % PER ARCSEC OF THE SKY
    assert (modelled["slit_normalisation_ratio"] == 1).all()
    assert _residual_slope_percent_per_arcsec(modelled) < -0.5


def test_sky_fit_is_unchanged_when_the_key_is_missing(log: Any) -> None:
    # ARRANGE
    order = _tilted_order()
    withKeyOff = _fit(_bare_subtractor(log, _sky_settings(**{SETTING_KEY: False}), {ORDER: KEPT_RANGE}), order)

    # ACT
    withoutKey = _fit(_bare_subtractor(log, _sky_settings(), {ORDER: KEPT_RANGE}), order)

    # ASSERT
    pd.testing.assert_frame_equal(withoutKey, withKeyOff)


def test_ratio_is_one_and_the_reason_logged_when_the_order_has_no_kept_range(log: Any) -> None:
    # ARRANGE
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {})

    # ACT
    modelled = _fit(subtractor, _tilted_order())

    # ASSERT
    assert (modelled["slit_normalisation_ratio"] == 1).all()
    assert any("no slit range to fit" in message for _, message in log.messages)


def test_ratio_is_one_and_the_reason_logged_when_the_frame_has_no_master_flat(log: Any) -> None:
    # ARRANGE
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {ORDER: KEPT_RANGE})

    # ACT
    modelled = _fit(subtractor, _tilted_order().drop(columns="flat"))

    # ASSERT
    assert (modelled["slit_normalisation_ratio"] == 1).all()
    assert [message for _, message in log.messages].count(
        "order 11: the slit ratio from the master flat is 1, because there is no master flat"
    ) == 1


def test_ratio_is_one_and_the_reason_logged_when_the_flat_slope_is_beyond_the_guard(log: Any) -> None:
    # ARRANGE: THE FLAT TILTS BY 40 %, ABOUT 3.6 % PER ARCSEC
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {ORDER: KEPT_RANGE})
    order = _tilted_order()
    steepProfile = 1.0 + 4 * FLAT_TILT * order["slit_position"] / (2 * SLIT_HALF_LENGTH)
    order["flat"] = order["flat"] / (1.0 + FLAT_TILT * order["slit_position"] / (2 * SLIT_HALF_LENGTH)) * steepProfile

    # ACT
    modelled = _fit(subtractor, order)

    # ASSERT
    assert (modelled["slit_normalisation_ratio"] == 1).all()
    assert any("slope" in message for _, message in log.messages)


def test_a_missing_key_means_the_flat_slit_gradient_is_off(log: Any) -> None:
    subtractor = _bare_subtractor(log, _sky_settings())

    assert subtractor._flat_slit_gradient_enabled() is False


@pytest.mark.parametrize("isOn", [True, False])
def test_a_boolean_key_is_taken_as_it_is(log: Any, isOn: bool) -> None:
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: isOn}))

    assert subtractor._flat_slit_gradient_enabled() is isOn


@pytest.mark.parametrize("isOn", [np.True_, np.False_])
def test_a_numpy_boolean_key_is_taken_as_a_boolean(log: Any, isOn: np.bool_) -> None:
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: isOn}))

    assert subtractor._flat_slit_gradient_enabled() is bool(isOn)


@pytest.mark.parametrize("badValue", ["True", "false", 1, 0, None])
def test_a_key_that_is_not_a_bool_raises_a_value_error_naming_the_key(log: Any, badValue: Any) -> None:
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: badValue}), {ORDER: KEPT_RANGE})

    with pytest.raises(ValueError, match=SETTING_KEY):
        subtractor._flat_slit_gradient_enabled()

    with pytest.raises(ValueError, match=SETTING_KEY):
        _fit(subtractor, _tilted_order())


def test_the_residuals_of_the_fit_are_measured_in_the_space_the_spline_was_fitted_in(log: Any) -> None:
    # ARRANGE: A TILTED FLAT-FIELDED ORDER, FITTED WITH THE SETTING ON, AND THE SPLINE OF THE LAST FIT
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {ORDER: KEPT_RANGE})
    modelled, tck, *_ = _fit_all(subtractor, _clipped(subtractor, _tilted_order()))
    ratio = modelled["slit_normalisation_ratio"].to_numpy(dtype=float)
    assert (ratio != 1).any()

    # ACT
    floored, _ = subtractor.determine_residual_floor(modelled.copy(), tck, 3)

    # ASSERT: (FLUX / RATIO - SPLINE) / (ERROR / RATIO), AS THE FIT AND THE RESIDUAL CLIP USE
    isUnclipped = ~floored["flagged_all_clipped"].to_numpy(dtype=bool)
    spline = np.asarray(scipy.interpolate.splev(floored["wavelength"].to_numpy(), tck))
    expected = np.abs(floored["flux"].to_numpy() - spline * ratio) / floored["error"].to_numpy()
    np.testing.assert_allclose(floored["sky_residuals"].to_numpy()[isUnclipped], expected[isUnclipped], rtol=1e-9)


def _untilted_twin(clipped: pd.DataFrame) -> pd.DataFrame:
    """A copy of a clipped untilted order, with every flux-like column scaled by the inverse of its flat profile."""
    slit = clipped["slit_position"].to_numpy(dtype=float)
    profile = 1.0 + FLAT_TILT * slit / (2 * SLIT_HALF_LENGTH)
    isUnclipped = ~clipped["flagged_all_clipped"].to_numpy(dtype=bool)
    scale = (1.0 / profile) / np.mean((1.0 / profile)[isUnclipped])
    twin = clipped.copy()
    for column in ("flux", "error", "residual_windowed_std", "flux_percentile_smoothed", "flux_windowed_long_median"):
        twin[column] = clipped[column].to_numpy(dtype=float) * scale
    twin["flat"] = clipped["flat"].to_numpy(dtype=float) * profile
    return twin


def test_a_tilted_flat_with_the_setting_on_gives_the_residual_floor_noise_flags_and_knots_of_the_untilted_case(
    log: Any,
) -> None:
    # ARRANGE: THE SAME ORDER TWICE. UNTILTED, AND TILTED BY A FLAT PROFILE THAT THE FIT DIVIDES OUT AGAIN
    untiltedSubtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: False}), {ORDER: KEPT_RANGE})
    untilted = _clipped(untiltedSubtractor, _tilted_order(tilt=0.0))
    flatTilt = _untilted_twin(untilted)
    tiltedSubtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {ORDER: KEPT_RANGE})

    # ACT
    reference = _fit_all(untiltedSubtractor, untilted)
    modelled = _fit_all(tiltedSubtractor, flatTilt)

    # ASSERT
    assert not (modelled[0]["slit_normalisation_ratio"] == 1).all()
    assert len(modelled[2]) == pytest.approx(len(reference[2]), rel=0.05)
    assert modelled[4] == pytest.approx(reference[4], rel=0.05)
    for column in ("flagged_noisy_region", "flagged_all_clipped"):
        assert modelled[0][column].sum() == pytest.approx(reference[0][column].sum(), abs=0.01 * len(untilted))
    skyLines = lambda frame: int((frame["flagged_sky_line"].astype(str) == "line").sum())  # noqa: E731
    assert skyLines(modelled[0]) == pytest.approx(skyLines(reference[0]), rel=0.1)
