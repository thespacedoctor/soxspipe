"""The NIR raw read-out noise is a robust median over consecutive raw-frame pairs (DY-2322)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes import base_recipe
from tests.factories import StubInputFrames, instrument_header, qc_table

pytestmark = pytest.mark.unit

SHAPE = (128, 128)
TRUE_SIGMA = 7.0
# A 128 x 128 FRAME GIVES A MAD SIGMA PRECISE TO ABOUT 1%, SO 3% IS A SAFE BOUND FOR A SEEDED DRAW
SIGMA_TOLERANCE = 0.03
# A HEAVY-TAIL SPIKE OF 60 e- IS ABOUT 6 SIGMA OF A PAIR DIFFERENCE, SO THE 10 SIGMA CLIP OF THE BIAS STATISTIC KEEPS IT
TAIL_AMPLITUDE = 60.0


def _frame(data: np.ndarray, arm: str = "NIR") -> CCDData:
    return CCDData(
        data.astype(np.float32),
        unit=u.electron,
        meta=instrument_header(arm=arm),
        mask=np.zeros(data.shape, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(data.shape, 1.0), unit=u.electron),
    )


def _recipe(log: Any, frames: list[CCDData], arm: str | None = "NIR") -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.recipeName = "soxs-mdark"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.arm = arm
    recipe.inputFrames = StubInputFrames(frames)
    return recipe


def _gaussian_frames(
    sigmas: list[float],
    *,
    tailFraction: float = 0.0,
    seed: int = 3,
    arm: str = "NIR",
) -> list[CCDData]:
    """Return frames with a shared fixed pattern, white noise per frame and optional heavy-tail pixels."""
    rng = np.random.default_rng(seed)
    pattern = rng.normal(loc=0.0, scale=40.0, size=SHAPE)
    frames = []
    for sigma in sigmas:
        data = 10.0 + pattern + rng.normal(scale=sigma, size=SHAPE)
        if tailFraction:
            tail = rng.random(SHAPE) < tailFraction
            data = data + tail * rng.choice([-1.0, 1.0], size=SHAPE) * TAIL_AMPLITUDE
        frames.append(_frame(data, arm))
    return frames


def test_nir_raw_ron_skips_the_settling_first_frame(log: Any) -> None:
    """With three or more frames the noisy first frame does not enter the NIR raw noise."""
    # ARRANGE
    # THE FIRST FRAME IS STILL SETTLING, SO PAIR 1-2 WOULD READ ABOUT 11.9 e-
    recipe = _recipe(log, _gaussian_frames([15.0] + [TRUE_SIGMA] * 5))

    # ACT
    rawRon, _masterRon = recipe.qc_ron(frameType="MDARK", frameName="master dark")

    # ASSERT
    assert rawRon == pytest.approx(TRUE_SIGMA, rel=SIGMA_TOLERANCE)
    rawRow = recipe.qc.set_index("qc_name").loc["RAW RON"]
    assert rawRow["qc_value"] == pytest.approx(rawRon)
    assert rawRow["qc_comment"] == "[e-] RON in single DARK"


def test_nir_raw_ron_ignores_heavy_tail_pixels(log: Any) -> None:
    """Pixels far out in the tail of the pair difference do not inflate the NIR raw noise."""
    # ARRANGE
    # 1% OF THE PIXELS CARRY A SPIKE THE 10 SIGMA CLIP KEEPS. A CLIPPED STANDARD DEVIATION WOULD READ ABOUT 9.2 e-
    recipe = _recipe(log, _gaussian_frames([TRUE_SIGMA] * 6, tailFraction=0.01))

    # ACT
    rawRon, _masterRon = recipe.qc_ron(frameType="MDARK", frameName="master dark")

    # ASSERT
    assert rawRon == pytest.approx(TRUE_SIGMA, rel=SIGMA_TOLERANCE)


def test_nir_raw_ron_with_two_frames_uses_that_one_pair(log: Any) -> None:
    """Two frames give one pair, so the first frame is not skipped."""
    # ARRANGE
    recipe = _recipe(log, _gaussian_frames([TRUE_SIGMA] * 2, tailFraction=0.01))

    # ACT
    rawRon, _masterRon = recipe.qc_ron(frameType="MDARK", frameName="master dark")

    # ASSERT
    assert rawRon == pytest.approx(TRUE_SIGMA, rel=SIGMA_TOLERANCE)


def test_nir_raw_ron_is_the_median_over_pairs(log: Any) -> None:
    """One noisy pair among three leaves the NIR raw noise at the typical pair value."""
    # ARRANGE
    # PAIR 2-3 READS ABOUT 15 e-; PAIRS 3-4 AND 4-5 READ 7 e-. THE MEDIAN IS 7 e-, THE MEAN WOULD BE ABOUT 9.7 e-
    recipe = _recipe(log, _gaussian_frames([TRUE_SIGMA, 20.0, TRUE_SIGMA, TRUE_SIGMA, TRUE_SIGMA]))

    # ACT
    rawRon, _masterRon = recipe.qc_ron(frameType="MDARK", frameName="master dark")

    # ASSERT
    assert rawRon == pytest.approx(TRUE_SIGMA, rel=SIGMA_TOLERANCE)


def test_nir_master_ron_still_uses_the_mask_of_the_first_pair(log: Any) -> None:
    """The master noise is measured under the mask of the frame 1-2 difference, as before."""
    # ARRANGE
    frames = _gaussian_frames([TRUE_SIGMA] * 4)
    # A 5000 e- SPIKE IN THE FIRST FRAME ONLY IS CLIPPED OUT OF THE 1-2 DIFFERENCE, AND OF NO LATER PAIR
    frames[0].data[3, 3] += 5000.0
    masterFrame = _frame(np.full(SHAPE, 5.0))
    masterFrame.data[3, 3] = 1000.0
    masterFrame.data[4, 4] = 9.0
    recipe = _recipe(log, frames)

    # ACT
    _rawRon, masterRon = recipe.qc_ron(frameType="MDARK", frameName="master dark", masterFrame=masterFrame)

    # ASSERT
    # UNMASKED, THE 1000 e- PIXEL WOULD GIVE A MASTER NOISE OF ABOUT 7.7
    assert masterRon == pytest.approx(float(np.ma.std(np.ma.array(masterFrame.data, mask=_spike_mask()))))
    assert masterRon < 0.2


def _spike_mask() -> np.ndarray:
    mask = np.zeros(SHAPE, dtype=bool)
    mask[3, 3] = True
    return mask


def test_nir_raw_ron_rejects_frames_with_no_measurable_spread(log: Any) -> None:
    """Identical NIR frames raise rather than report zero read-out noise."""
    # ARRANGE
    flat = np.full(SHAPE, 12.0)
    recipe = _recipe(log, [_frame(flat), _frame(flat), _frame(flat)])

    # ACT / ASSERT
    with pytest.raises(ValueError, match="appear to be corrupted"):
        recipe.qc_ron(frameType="MDARK", frameName="master dark")


@pytest.mark.parametrize("arm", ["VIS", "UVB", None])
def test_non_nir_raw_ron_still_measures_only_the_first_pair(log: Any, arm: str | None) -> None:
    """The bias path keeps the clipped standard deviation of frames 1-2, whatever the later frames hold."""
    # ARRANGE
    # THE FIRST PAIR IS NOISE 7 e-; THE THIRD FRAME IS FAR NOISIER AND MUST NOT MOVE THE RESULT
    twoFrames = _gaussian_frames([TRUE_SIGMA] * 2, arm=arm or "VIS")
    threeFrames = _gaussian_frames([TRUE_SIGMA] * 2 + [60.0], arm=arm or "VIS")

    # ACT
    pairRon, _ = _recipe(log, twoFrames, arm).qc_ron(frameType="MBIAS", frameName="master bias")
    tripletRon, _ = _recipe(log, threeFrames, arm).qc_ron(frameType="MBIAS", frameName="master bias")

    # ASSERT
    assert tripletRon == pairRon
