"""The bias pedestal must not inflate the ERRS of a VIS nod or offset difference frame (DY-1390)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes import soxs_nod

pytestmark = pytest.mark.unit

PEDESTAL = 650.0
READ_NOISE = 3.5
SKY_A = 120.0
SKY_B = 80.0
SHAPE = (200, 200)


def _recipe(log: Any, *, arm: str = "VIS", recipeName: str = "soxs-nod") -> soxs_nod:
    """Return an unconstructed recipe carrying only what the A-B difference reads."""
    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.arm = arm
    recipe.recipeName = recipeName
    recipe.detectorParams = {"ron": READ_NOISE * u.electron}
    return recipe


def _prepared_frame(sky: float, seed: int) -> CCDData:
    """Return a raw sky frame on the pedestal, with the ERRS that preparation derives from its raw counts."""
    rng = np.random.default_rng(seed)
    data = PEDESTAL + rng.poisson(sky, SHAPE) + rng.normal(0.0, READ_NOISE, SHAPE)
    error = np.sqrt(np.clip(data, 0, None) + READ_NOISE**2)
    return CCDData(
        data.astype(np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(error.astype(np.float32), unit=u.electron),
        mask=np.zeros(SHAPE, dtype=bool),
    )


def _master_bias() -> CCDData:
    """Return a noiseless master bias at the pedestal level."""
    return CCDData(np.full(SHAPE, PEDESTAL, dtype=np.float32), unit=u.electron)


def test_a_vis_a_minus_b_error_matches_its_photon_and_read_noise(log: Any) -> None:
    """With the master bias, the A-B residuals over ERRS have unit scatter in sky pixels."""
    recipe = _recipe(log)
    recipe.masterBias = _master_bias()
    aFrame, bFrame = _prepared_frame(SKY_A, seed=1), _prepared_frame(SKY_B, seed=2)

    aMinusB, bMinusA = recipe._difference_cycle_frames(aFrame, bFrame)

    expectedError = np.sqrt(SKY_A + SKY_B + 2 * READ_NOISE**2)
    for difference, expectedMean in ((aMinusB, SKY_A - SKY_B), (bMinusA, SKY_B - SKY_A)):
        zScores = (difference.data - expectedMean) / difference.uncertainty.array
        assert np.std(zScores) == pytest.approx(1.0, abs=0.02)
        assert np.median(difference.uncertainty.array) == pytest.approx(expectedError, rel=0.01)


@pytest.mark.parametrize("recipeName", ["soxs-nod", "soxs-nod-std", "soxs-offset"])
def test_the_master_bias_leaves_the_difference_data_bit_identical(log: Any, recipeName: str) -> None:
    """Only the error changes: A-B (and, for nodding, B-A) data match a run without the master bias."""
    withBias = _recipe(log, recipeName=recipeName)
    withBias.masterBias = _master_bias()
    withoutBias = _recipe(log, recipeName=recipeName)
    aFrame, bFrame = _prepared_frame(SKY_A, seed=3), _prepared_frame(SKY_B, seed=4)

    correctedDifferences = withBias._difference_cycle_frames(aFrame, bFrame)
    plainDifferences = withoutBias._difference_cycle_frames(aFrame, bFrame)

    np.testing.assert_array_equal(correctedDifferences[0].data, plainDifferences[0].data)
    if "nod" in recipeName:
        np.testing.assert_array_equal(correctedDifferences[1].data, plainDifferences[1].data)
    else:
        assert correctedDifferences[1] is False


def test_without_a_master_bias_the_difference_error_is_the_raw_quadrature_sum(log: Any) -> None:
    """NIR, or a UVB/VIS SOF with no master bias, keeps today's error: the two raw-count errors in quadrature."""
    recipe = _recipe(log, arm="NIR")
    aFrame, bFrame = _prepared_frame(SKY_A, seed=5), _prepared_frame(SKY_B, seed=6)

    aMinusB, _ = recipe._difference_cycle_frames(aFrame, bFrame)

    expectedError = np.hypot(aFrame.uncertainty.array, bFrame.uncertainty.array)
    np.testing.assert_allclose(aMinusB.uncertainty.array, expectedError, rtol=1e-6)


def test_the_difference_leaves_the_input_frame_errors_untouched(log: Any) -> None:
    """The A and B frames are reused after the difference, so their own errors must not change."""
    recipe = _recipe(log)
    recipe.masterBias = _master_bias()
    aFrame, bFrame = _prepared_frame(SKY_A, seed=7), _prepared_frame(SKY_B, seed=8)
    aBefore, bBefore = aFrame.uncertainty.array.copy(), bFrame.uncertainty.array.copy()

    recipe._difference_cycle_frames(aFrame, bFrame)

    np.testing.assert_array_equal(aFrame.uncertainty.array, aBefore)
    np.testing.assert_array_equal(bFrame.uncertainty.array, bBefore)


class _InputFrames:
    """The one ``files_filtered`` lookup the master-bias reader makes, over a fixed catalogue of products."""

    def __init__(self, catalogue: dict[str, str]) -> None:
        self.catalogue = catalogue
        self.lookups: list[dict[str, object]] = []

    def files_filtered(self, include_path: bool = False, **filters: object) -> list[str]:
        self.lookups.append(filters)
        return [path for category, path in self.catalogue.items() if filters.get("PRO_CATG") == category]


def _bias_file(tmp_path: Path) -> str:
    """Write a master bias with an ERRS extension and return its path."""
    path = tmp_path / "MASTER_BIAS.fits"
    bias = _master_bias()
    bias.uncertainty = StdDevUncertainty(np.full(SHAPE, 0.4, dtype=np.float32), unit=u.electron)
    bias.write(path, hdu_uncertainty="ERRS", hdu_mask="QUAL", hdu_flags="FLAGS", key_uncertainty_type="UTYPE")
    return str(path)


@pytest.mark.parametrize("arm", ["VIS", "UVB"])
def test_a_uvb_vis_master_bias_named_in_the_sof_is_read(log: Any, tmp_path: Path, arm: str) -> None:
    """The master bias the nod and offset SOFs already name is read, data and error."""
    recipe = _recipe(log, arm=arm)
    recipe.inputFrames = _InputFrames({f"MASTER_BIAS_{arm}": _bias_file(tmp_path)})

    masterBias = recipe._read_master_bias(kw=lambda keyword: keyword, arm=arm)

    np.testing.assert_array_equal(masterBias.data, np.full(SHAPE, PEDESTAL))
    np.testing.assert_allclose(masterBias.uncertainty.array, np.full(SHAPE, 0.4), rtol=1e-6)
    assert not [message for level, message in log.messages if level == "warning"]


def test_a_vis_sof_without_a_master_bias_warns_and_keeps_todays_error(log: Any) -> None:
    """With no master bias the reader returns False, so the difference error is as before, and says why."""
    recipe = _recipe(log)
    recipe.inputFrames = _InputFrames({})

    masterBias = recipe._read_master_bias(kw=lambda keyword: keyword, arm="VIS")

    assert masterBias is False
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "MASTER_BIAS_VIS" in warnings[0]
    assert "bias pedestal" in warnings[0]


def test_a_nir_reduction_never_looks_for_a_master_bias(log: Any) -> None:
    """NIR frames carry no bias pedestal, so NIR is unchanged: no lookup, no warning."""
    recipe = _recipe(log, arm="NIR")
    recipe.inputFrames = _InputFrames({"MASTER_BIAS_NIR": "never-read.fits"})

    masterBias = recipe._read_master_bias(kw=lambda keyword: keyword, arm="NIR")

    assert masterBias is False
    assert recipe.inputFrames.lookups == []
    assert not [message for level, message in log.messages if level == "warning"]
