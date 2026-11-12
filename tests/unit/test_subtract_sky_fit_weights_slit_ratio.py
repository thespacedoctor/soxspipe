"""The sky B-spline fit weights match the data that is fitted, ``flux / slit_normalisation_ratio`` (DY-2152).

The spline is fitted to ``flux / r``, whose noise is ``sigma / r``, so each weight is the noise-model weight of the sky
in raw-flux units multiplied by ``r``. The sky that goes into the noise model is never the measured flux (DY-1282): it
is the percentile-smoothed sky for the first fit, and the previous spline multiplied by ``r`` for every later fit.
With ``r = 1`` the weights are those the fit had before DY-2152.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate

from soxspipe.commonutils.subtract_sky import ANCHOR_WEIGHT, subtract_sky
from tests.unit.test_subtract_sky_flat_slit_gradient import (
    KEPT_RANGE,
    ORDER,
    SETTING_KEY,
    _bare_subtractor,
    _clipped,
    _sky_settings,
    _tilted_order,
)

pytestmark = pytest.mark.unit


class _FitRecord:
    """What each spline fit of one ``fit_bspline_curve_to_sky`` call was given and returned, in call order."""

    def __init__(self) -> None:
        self.rowLabels: list[Any] = []
        self.weights: list[np.ndarray] = []
        self.tcks: list[Any] = []


def _fit_recording_weights(
    subtractor: subtract_sky, clipped: pd.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> tuple[pd.DataFrame, Any, _FitRecord]:
    """Run ``fit_bspline_curve_to_sky``, recording the rows and weights of every fit and the spline each returned."""
    record = _FitRecord()
    realSamples = subtractor._spline_fit_samples
    realSplrep = scipy.interpolate.splrep

    def samples(imageMapOrder: pd.DataFrame, mask: pd.Series, starterKnots: Any) -> tuple:
        goodWl, goodFlux, goodWeights = realSamples(imageMapOrder, mask, starterKnots)
        record.rowLabels.append(goodWl.index)
        record.weights.append(goodWeights.copy())
        return goodWl, goodFlux, goodWeights

    def splrep(*args: Any, **kwargs: Any) -> tuple:
        result = realSplrep(*args, **kwargs)
        record.tcks.append(result[0])
        return result

    monkeypatch.setattr(subtractor, "_spline_fit_samples", samples)
    monkeypatch.setattr(scipy.interpolate, "splrep", splrep)
    modelled, tck, *_ = subtractor.fit_bspline_curve_to_sky(clipped.copy())
    return modelled, tck, record


def _expected_weights(modelled: pd.DataFrame, record: _FitRecord, ron: float, index: int) -> np.ndarray:
    """The interior weights fit ``index`` should have: the noise-model weight of the raw-flux sky, times the ratio."""
    labels = record.rowLabels[index][1:-1]
    rows = modelled.loc[labels]
    ratio = rows["slit_normalisation_ratio"].to_numpy(dtype=float)
    if index == 0:
        filled = subtract_sky._fill_nans_by_interpolation(modelled["flux_percentile_smoothed"].to_numpy())
        rawSky = pd.Series(filled, index=modelled.index).loc[labels].to_numpy()
    else:
        rawSky = scipy.interpolate.splev(rows["wavelength"].to_numpy(), record.tcks[index - 1]) * ratio
    return subtract_sky._inverse_noise_weights(rawSky, ron, rows["flat"].to_numpy(dtype=float)) * ratio


def test_fit_weights_are_the_noise_model_weights_times_the_slit_ratio_for_a_non_uniform_ratio(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: A TILTED ORDER WITH THE SETTING ON, SO THE RATIO VARIES ALONG THE SLIT
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {ORDER: KEPT_RANGE})
    clipped = _clipped(subtractor, _tilted_order())

    # ACT
    modelled, _, record = _fit_recording_weights(subtractor, clipped, monkeypatch)

    # ASSERT: EVERY FIT, THE FIRST BY THE PERCENTILE-SMOOTHED SKY AND THE REST BY THE PREVIOUS SPLINE TIMES THE RATIO
    assert modelled["slit_normalisation_ratio"].nunique() > 1
    assert len(record.weights) > 2
    for index, weights in enumerate(record.weights):
        expected = _expected_weights(modelled, record, subtractor.ron, index)
        np.testing.assert_allclose(weights[1:-1], expected, rtol=1e-12, err_msg=f"fit {index}")
        assert weights[0] == ANCHOR_WEIGHT
        assert weights[-1] == ANCHOR_WEIGHT


@pytest.mark.parametrize("isOn", [False, True])
def test_fit_weights_and_model_are_unchanged_when_the_ratio_is_one(
    log: Any, monkeypatch: pytest.MonkeyPatch, isOn: bool
) -> None:
    # ARRANGE: THE SETTING OFF, AND THE SETTING ON FOR AN ORDER WITH NO TILT THAT THE FLAT-SLIT RATIO FALLS BACK FROM
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: isOn}), {})
    clipped = _clipped(subtractor, _tilted_order())

    # ACT
    modelled, tck, record = _fit_recording_weights(subtractor, clipped, monkeypatch)

    # ASSERT: THE WEIGHTS ARE EXACTLY THE NOISE-MODEL WEIGHTS OF THE RAW-FLUX SKY, NOT MULTIPLIED BY ANYTHING
    assert (modelled["slit_normalisation_ratio"] == 1).all()
    for index, weights in enumerate(record.weights):
        expected = _expected_weights(modelled, record, subtractor.ron, index)
        np.testing.assert_array_equal(weights[1:-1], expected, err_msg=f"fit {index}")
    wavelength = modelled["wavelength"].to_numpy()
    np.testing.assert_array_equal(modelled["sky_model"].to_numpy(), scipy.interpolate.splev(wavelength, tck))


def test_refit_weights_are_the_noise_model_weights_of_the_raw_flux_sky_times_the_slit_ratio(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: A FITTED TILTED ORDER WITH A FEW PIXELS FLAGGED SO THE CLIP-AND-REFIT STEP REFITS ONCE
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}), {ORDER: KEPT_RANGE})
    modelled, tck, allKnots, *_ = subtractor.fit_bspline_curve_to_sky(_clipped(subtractor, _tilted_order()).copy())
    flagged = modelled.index[:: len(modelled) // 20]
    modelled.loc[flagged, ["flagged_bspline_clipped", "flagged_all_clipped"]] = True
    seen: list[np.ndarray] = []
    realSamples = subtractor._spline_fit_samples

    def samples(imageMapOrder: pd.DataFrame, mask: pd.Series, starterKnots: Any) -> tuple:
        seen.append(imageMapOrder["weights"].to_numpy(dtype=float, copy=True))
        return realSamples(imageMapOrder, mask, starterKnots)

    monkeypatch.setattr(subtractor, "_spline_fit_samples", samples)

    # ACT: A CLIPPING LIMIT NO PIXEL EXCEEDS, SO THE STEP STOPS AFTER THE FIRST REFIT
    subtractor._clip_residual_outliers_and_refit(modelled, tck, allKnots, allKnots, 1e9)

    # ASSERT
    ratio = modelled["slit_normalisation_ratio"].to_numpy(dtype=float)
    assert (ratio != 1).any()
    rawSky = scipy.interpolate.splev(modelled["wavelength"].to_numpy(), tck) * ratio
    expected = (
        subtract_sky._inverse_noise_weights(rawSky, subtractor.ron, modelled["flat"].to_numpy(dtype=float)) * ratio
    )
    assert len(seen) == 1
    np.testing.assert_allclose(seen[0], expected, rtol=1e-12)
