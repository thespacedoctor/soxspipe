"""Public QC contracts for base recipe utility methods."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.recipes import base_recipe
from tests.factories import qc_table, synthetic_ccd

pytestmark = pytest.mark.unit


def _recipe(log: Any) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.recipeName = "soxs-mbias"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.qc = qc_table().iloc[0:0].copy()
    return recipe


def test_flag_poor_data_adds_detector_metrics_and_marks_outside_range_failed(
    log: Any,
) -> None:
    """SOXS QC flagging records temperatures and rejects out-of-range metrics."""
    recipe = _recipe(log)
    recipe.inst = "SOXS"
    recipe.detectorTemp = 82.5
    recipe.cptemp = 11.2
    recipe.recipeSettings = {"qc-acceptable-ranges": {"trace-rms": [0.0, 2.0]}}
    recipe.qc = pd.DataFrame(
        [
            {
                "soxspipe_recipe": "soxs-mbias",
                "qc_name": "TRACE RMS",
                "qc_value": 2.5,
                "qc_unit": "pixel",
                "qc_order": None,
                "qc_comment": "Synthetic trace residual",
                "obs_date_utc": recipe.dateObs,
                "reduction_date_utc": recipe.dateObs,
                "to_header": True,
            },
            {
                "soxspipe_recipe": "soxs-mbias",
                "qc_name": "RON",
                "qc_value": 1.2,
                "qc_unit": "electron",
                "qc_order": None,
                "qc_comment": "Synthetic read noise",
                "obs_date_utc": recipe.dateObs,
                "reduction_date_utc": recipe.dateObs,
                "to_header": True,
            },
        ]
    )

    result = recipe.flag_poor_data()

    assert result is None
    flagged = recipe.qc.set_index("qc_name")
    assert flagged.loc["TRACE RMS", "qc_flag"] == "fail"
    assert flagged.loc["TRACE RMS", "qc_value_min"] == 0.0
    assert flagged.loc["TRACE RMS", "qc_value_max"] == 2.0
    assert flagged.loc["RON", "qc_flag"] == "pass"
    assert flagged.loc["DETECTOR TEMP", "qc_value"] == 82.5
    assert flagged.loc["CPATH TEMP", "qc_value"] == 11.2
    assert set(flagged["qc_order"]) == {"-1"}


def test_qc_median_flux_level_excludes_masked_pixels_and_records_metric(
    log: Any,
) -> None:
    """Median-flux QC ignores masked pixels before recording its result."""
    recipe = _recipe(log)
    recipe.kw = lambda keyword: keyword
    frame = synthetic_ccd(
        shape=(2, 2),
        prepared=True,
        headerOverrides={"EXPTIME": 2.0},
    )
    frame.data = np.array([[1.0, 3.0], [5.0, 1000.0]])
    frame.mask = np.array([[False, False], [False, True]])

    medianFlux = recipe.qc_median_flux_level(
        frame,
        frameType="MDARK",
        frameName="synthetic dark",
    )

    assert medianFlux == pytest.approx(3.0)
    metric = recipe.qc.iloc[-1]
    assert metric["qc_name"] == "MDARK MEDIAN"
    assert metric["qc_value"] == pytest.approx(3.0)
    assert metric["qc_comment"] == "[e-] Median flux level of synthetic dark"
    assert metric["qc_unit"] == "electrons"
    assert metric["to_header"] is True


def test_qc_ron_records_supplied_raw_and_master_noise_values(log: Any) -> None:
    """Precomputed noise values are retained without loading raw frame data."""
    recipe = _recipe(log)

    rawRon, masterRon = recipe.qc_ron(
        frameType="MBIAS",
        frameName="master bias",
        rawRon=2.5,
        masterRon=0.8,
    )

    assert rawRon == pytest.approx(2.5)
    assert masterRon == pytest.approx(0.8)
    metrics = recipe.qc.set_index("qc_name")
    assert metrics.loc["RAW RON", "qc_value"] == pytest.approx(2.5)
    assert metrics.loc["RAW RON", "qc_comment"] == "[e-] RON in single BIAS"
    assert metrics.loc["MASTER RON", "qc_value"] == pytest.approx(0.8)
    assert metrics.loc["MASTER RON", "qc_comment"] == "[e-] Combined RON in MBIAS"


def test_subtract_mean_flux_level_clips_outlier_and_preserves_frame_contract(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Remove a sigma-clipped raw-frame mean and retain the CCDData result."""
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    recipe.recipeSettings = {
        "frame-clipping-sigma": 2.0,
        "frame-clipping-iterations": 3,
    }
    frame = synthetic_ccd(shape=(3, 3), prepared=True)
    frame.data = np.full((3, 3), 10.0)
    frame.data[2, 2] = 1000.0
    calls: list[object] = []
    monkeypatch.setattr(toolkit, "frame_to_32", lambda value: calls.append(value))

    meanFlux, fluxStandardDeviation, noiseFrame = recipe.subtract_mean_flux_level(frame)

    assert meanFlux == pytest.approx(10.0)
    assert fluxStandardDeviation == pytest.approx(0.0)
    assert noiseFrame is frame
    assert noiseFrame.data[0, 0] == pytest.approx(0.0)
    assert noiseFrame.mask[2, 2]
    assert calls == [frame]


def test_xshooter_rotation_and_binned_trim_preserve_science_region(log: Any) -> None:
    """The X-shooter compatibility and trim helpers retain the requested pixels."""
    recipe = _recipe(log)
    recipe.kw = lambda keyword: keyword
    recipe.arm = "VIS"
    recipe.detectorParams = {
        "clockwise-rotation": 90.0,
        "binning": [2, 2],
        "science-pixels": {
            "rows": {"start": 2, "end": 6},
            "columns": {"start": 2, "end": 6},
        },
    }
    frame = CCDData(np.arange(36, dtype=float).reshape(6, 6), unit=u.adu)

    rotated = recipe.xsh2soxs(frame)
    trimmed = recipe._trim_frame(rotated)

    np.testing.assert_array_equal(rotated.data, np.rot90(np.arange(36).reshape(6, 6)))
    np.testing.assert_array_equal(trimmed.data, rotated.data[1:3, 1:3])
    assert trimmed.shape == (2, 2)
