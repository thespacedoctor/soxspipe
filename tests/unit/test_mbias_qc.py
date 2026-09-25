"""Characterization tests for master-bias quality-control calculations."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from astropy.nddata import CCDData

from soxspipe.recipes.soxs_mbias import soxs_mbias

pytestmark = pytest.mark.unit


class _Frames:
    """Expose the image-collection surface needed by periodic-noise QC."""

    def __init__(self, frames: list[CCDData]) -> None:
        self._frames = frames

    def ccds(self, *, ccd_kwargs: dict[str, str]) -> list[CCDData]:
        assert ccd_kwargs["hdu_uncertainty"] == "ERRS"
        assert ccd_kwargs["hdu_mask"] == "QUAL"
        return self._frames


def _recipe(log: object) -> soxs_mbias:
    """Return the minimal mutable state consumed by the QC methods."""
    recipe = soxs_mbias.__new__(soxs_mbias)
    recipe.log = log
    recipe.qc = pd.DataFrame()
    recipe.recipeName = "soxs-mbias"
    recipe.dateObs = "2024-01-02T03:04:05"
    return recipe


def test_qc_bias_structure_records_axis_slopes(log: object) -> None:
    """Collapsed linear bias gradients become named STRUCT QC measurements."""
    recipe = _recipe(log)
    yPixels, xPixels = np.indices((4, 5))
    frame = 2.0 * xPixels + 3.0 * yPixels

    structX, structY = recipe.qc_bias_structure(frame)

    expectedX = np.polyfit(
        np.linspace(0, frame.shape[1], frame.shape[1], dtype=int),
        np.nansum(frame, axis=0),
        deg=1,
    )[0]
    expectedY = np.polyfit(
        np.linspace(0, frame.shape[0], frame.shape[0], dtype=int),
        np.nansum(frame, axis=1),
        deg=1,
    )[0]
    assert structX == pytest.approx(expectedX)
    assert structY == pytest.approx(expectedY)
    assert recipe.qc["qc_name"].tolist() == ["STRUCTX", "STRUCTY"]
    assert recipe.qc["qc_value"].tolist() == pytest.approx([structX, structY])
    assert recipe.qc["to_header"].tolist() == [True, True]


def test_qc_periodic_pattern_noise_records_maximum_frame_ratio(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """FFT-derived periodic-noise QC is finite and retains its header contract."""
    recipe = _recipe(log)
    random = np.random.default_rng(42)
    first = CCDData(random.normal(size=(20, 20)), unit="electron")
    second = CCDData(random.normal(size=(20, 20)), unit="electron")
    first.mask = np.zeros(first.shape, dtype=bool)
    second.mask = np.zeros(second.shape, dtype=bool)
    quicklookCalls: list[object] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.quicklook_image",
        lambda **kwargs: quicklookCalls.append(kwargs["CCDObject"]),
    )

    periodicNoise = recipe.qc_periodic_pattern_noise(_Frames([first, second]))

    assert np.isfinite(periodicNoise)
    assert periodicNoise > 0
    assert len(quicklookCalls) == 4
    assert recipe.qc["qc_name"].tolist() == ["FPN FRACMAX"]
    assert recipe.qc.loc[0, "qc_value"] == pytest.approx(periodicNoise)
    # NUMPY BOOL, NOT PYTHON BOOL, NOW THE COLUMN CARRIES A REAL BOOL DTYPE
    assert bool(recipe.qc.loc[0, "to_header"]) is True
