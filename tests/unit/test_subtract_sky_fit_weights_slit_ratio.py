"""The sky B-spline fit weights match the data that is fitted, ``flux / slit_normalisation_ratio`` (DY-2152).

The spline is fitted to ``flux / r``, whose noise is ``sigma / r``, so each weight is the noise-model weight of the sky
in raw-flux units multiplied by ``r``. The sky that goes into the noise model is never the measured flux (DY-1282): it
is the percentile-smoothed sky for the first fit, and the previous spline multiplied by ``r`` for every later fit.
With ``r = 1`` the weights are those the fit had before DY-2152. The residual clip measures the noise of the same fitted
data, so its noise is the noise-model sigma of the raw-flux sky divided by ``r``.
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

# THE NUMBER OF EVENLY SPACED PIXELS FLAGGED TO MAKE THE CLIP-AND-REFIT STEP REFIT
REFIT_FLAGGED_PIXELS = 20


class _FitRecord:
    """What each spline fit of a fit call was given and returned, in call order."""

    def __init__(self) -> None:
        self.rowLabels: list[Any] = []
        self.weights: list[np.ndarray] = []
        self.weightColumns: list[np.ndarray] = []
        self.tcks: list[Any] = []


def _record_spline_fits(subtractor: subtract_sky, monkeypatch: pytest.MonkeyPatch) -> _FitRecord:
    """Spy on the sample selection and on ``splrep`` of ``subtractor``, recording every fit, and return the record."""
    record = _FitRecord()
    realSamples = subtractor._spline_fit_samples
    realSplrep = scipy.interpolate.splrep

    def samples(imageMapOrder: pd.DataFrame, mask: pd.Series, starterKnots: Any) -> tuple:
        record.weightColumns.append(imageMapOrder["weights"].to_numpy(dtype=float, copy=True))
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
    return record


def _fit_recording_weights(
    subtractor: subtract_sky, clipped: pd.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> tuple[pd.DataFrame, Any, _FitRecord]:
    """Run ``fit_bspline_curve_to_sky``, recording the rows and weights of every fit and the spline each returned."""
    record = _record_spline_fits(subtractor, monkeypatch)
    modelled, tck, *_ = subtractor.fit_bspline_curve_to_sky(clipped.copy())
    return modelled, tck, record


def _weights_of_raw_sky(frame: pd.DataFrame, rawSky: np.ndarray, ron: float) -> np.ndarray:
    """The weights the fit should have: the noise-model weight of the raw-flux sky, times the slit ratio."""
    ratio = frame["slit_normalisation_ratio"].to_numpy(dtype=float)
    return subtract_sky._inverse_noise_weights(rawSky, ron, frame["flat"].to_numpy(dtype=float)) * ratio


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
    return _weights_of_raw_sky(rows, rawSky, ron)


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
    flagged = modelled.index[:: len(modelled) // REFIT_FLAGGED_PIXELS]
    modelled.loc[flagged, ["flagged_bspline_clipped", "flagged_all_clipped"]] = True
    record = _record_spline_fits(subtractor, monkeypatch)

    # ACT: A CLIPPING LIMIT NO PIXEL EXCEEDS, SO THE STEP STOPS AFTER THE FIRST REFIT
    subtractor._clip_residual_outliers_and_refit(modelled, tck, allKnots, allKnots, 1e9)

    # ASSERT
    ratio = modelled["slit_normalisation_ratio"].to_numpy(dtype=float)
    assert (ratio != 1).any()
    rawSky = scipy.interpolate.splev(modelled["wavelength"].to_numpy(), tck) * ratio
    assert len(record.weightColumns) == 1
    np.testing.assert_allclose(
        record.weightColumns[0], _weights_of_raw_sky(modelled, rawSky, subtractor.ron), rtol=1e-12
    )


def test_fit_weights_match_hand_computed_values(log: Any) -> None:
    # ARRANGE: FOUR PIXELS WITH A KNOWN RAW-FLUX SKY, FLAT AND RATIO, AND A READ NOISE OF 10 E-
    #   PIXEL 0: SKY 0,    FLAT 1, RATIO 1.0 -> 1 / SQRT(0 + 100)           * 1.0 = 0.1
    #   PIXEL 1: SKY 300,  FLAT 1, RATIO 2.0 -> 1 / SQRT(300 + 100)         * 2.0 = 0.1
    #   PIXEL 2: SKY 300,  FLAT 2, RATIO 0.5 -> 2 / SQRT(300 * 2 + 100)     * 0.5 = 0.0377964473009227
    #   PIXEL 3: SKY -5,   FLAT 1, RATIO 1.5 -> 1 / SQRT(0 + 100)           * 1.5 = 0.15 (A NEGATIVE SKY COUNTS AS 0)
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}))
    subtractor.ron = 10.0
    frame = pd.DataFrame({"flat": [1.0, 1.0, 2.0, 1.0], "slit_normalisation_ratio": [1.0, 2.0, 0.5, 1.5]})

    # ACT
    weights = subtractor._fit_weights(frame, np.array([0.0, 300.0, 300.0, -5.0]))

    # ASSERT
    np.testing.assert_allclose(weights, [0.1, 0.1, 0.0377964473009227, 0.15], rtol=1e-12)


def _clip_noise_frame(outlierRatio: float, outlierOffset: float) -> tuple[pd.DataFrame, Any]:
    """41 pixels fitted by a constant fitted-space sky of 200 e-, one with a residual and a slit ratio of its own."""
    wavelength = np.linspace(0.0, 1.0, 41)
    ratio = np.ones(41)
    ratio[20] = outlierRatio
    residual = np.zeros(41)
    residual[20] = outlierOffset
    frame = pd.DataFrame(
        {
            "wavelength": wavelength,
            "flux": ratio * (200.0 + residual),
            "slit_normalisation_ratio": ratio,
            "flat": 1.0,
            "flagged_all_clipped": False,
            "flagged_bspline_clipped": False,
        }
    )
    return frame, (np.array([0.0] * 4 + [1.0] * 4), np.array([200.0] * 4 + [0.0] * 4), 3)


@pytest.mark.parametrize(
    ("ratio", "offset", "isFlagged"),
    [
        # RAW SKY 400 E-: NOISE SQRT(400 + 100) / 2 = 11.18, SO THE 3 SIGMA LIMIT IS 33.5 (UNSCALED: 52.0)
        (2.0, 40.0, True),
        # RAW SKY 100 E-: NOISE SQRT(100 + 100) / 0.5 = 28.28, SO THE 3 SIGMA LIMIT IS 84.9 (UNSCALED: 52.0)
        (0.5, 70.0, False),
    ],
)
def test_residual_clip_noise_is_the_raw_flux_noise_of_the_sky_divided_by_the_slit_ratio(
    log: Any, ratio: float, offset: float, isFlagged: bool
) -> None:
    # ARRANGE: A READ NOISE OF 10 E- AND A FLAT OF 1, SO THE NOISE OF THE RAW-FLUX SKY IS SQRT(SKY + 100)
    subtractor = _bare_subtractor(log, _sky_settings(**{SETTING_KEY: True}))
    subtractor.ron = 10.0
    frame, tck = _clip_noise_frame(ratio, offset)

    # ACT
    newlyClipped = subtractor._flag_residual_outliers(frame, tck, 3.0)

    # ASSERT: ONLY THE OFFSET PIXEL CAN BE FLAGGED, AND ONLY WHEN ITS RESIDUAL EXCEEDS 3 SIGMA OF THE FITTED-SPACE NOISE
    assert newlyClipped.sum() == int(isFlagged)
    assert newlyClipped[20] == isFlagged
