"""The read noise of a UVB/VIS science frame, measured from the prescan of the raw frame (DY-2437).

The prescan is a strip of columns the detector reads without exposing it, so its scatter is the read noise of
that frame alone. It is measured on the raw frame in ADU, before the frame is trimmed, with a robust standard
deviation, and converted to e- with the gain.
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u

from soxspipe.commonutils.prescan_read_noise import (
    measure_prescan_read_noise,
    uses_prescan_read_noise,
)
from soxspipe.commonutils.read_noise import SOURCE_SCIENCE_PRESCAN, ReadNoise, science_read_noise

pytestmark = pytest.mark.unit

GAIN = 1.1 * u.electron / u.adu
PRESCAN_SIGMA_ADU = 4.6
SHAPE = (400, 120)
# THE PRESCAN IS COLUMNS 8 TO 46 OF EVERY ROW BUT THE 50 AT EACH END, AS IN THE VIS DETECTOR PARAMETERS
PRESCAN = {"rows": {"start": 50, "end": 350}, "columns": {"start": 8, "end": 46}}
TOLERANCE = 0.03


def _raw_frame(seed: int = 0, sigma: float = PRESCAN_SIGMA_ADU) -> np.ndarray:
    """A raw frame whose prescan has scatter `sigma` and whose science region is far noisier and brighter."""
    rng = np.random.default_rng(seed)
    frame = rng.normal(600.0, 30.0, SHAPE)
    frame[:, 8:46] = rng.normal(590.0, sigma, (SHAPE[0], 38))
    return frame


def test_the_read_noise_is_the_prescan_scatter_in_electrons() -> None:
    # ARRANGE
    frame = _raw_frame()

    # ACT
    ron = measure_prescan_read_noise(frame, PRESCAN, [1, 1], GAIN)

    # ASSERT
    assert ron == pytest.approx(PRESCAN_SIGMA_ADU * 1.1, rel=TOLERANCE)


def test_a_gain_given_as_a_plain_number_is_in_electrons_per_adu() -> None:
    # ACT
    ron = measure_prescan_read_noise(_raw_frame(), PRESCAN, [1, 1], 1.1)

    # ASSERT
    assert ron == pytest.approx(PRESCAN_SIGMA_ADU * 1.1, rel=TOLERANCE)


def test_cosmic_rays_in_the_prescan_do_not_bias_the_read_noise() -> None:
    # ARRANGE
    frame = _raw_frame()
    rng = np.random.default_rng(1)
    rows = rng.integers(60, 340, 450)
    columns = rng.integers(8, 46, 450)
    frame[rows, columns] += 5000.0

    # ACT
    ron = measure_prescan_read_noise(frame, PRESCAN, [1, 1], GAIN)

    # ASSERT
    assert ron == pytest.approx(PRESCAN_SIGMA_ADU * 1.1, rel=TOLERANCE)


def test_a_prescan_in_whole_adu_counts_is_not_rounded_to_the_median_absolute_deviation_steps() -> None:
    # ARRANGE: THE RAW FRAME IS INTEGER ADU, SO THE MEDIAN ABSOLUTE DEVIATION OF A PRESCAN WITH A 2.5 ADU SCATTER
    # CAN ONLY BE A WHOLE NUMBER OR A HALF, WHICH IS A 15 TO 20 % STEP IN THE READ NOISE
    frame = np.round(_raw_frame(sigma=2.5))

    # ACT
    ron = measure_prescan_read_noise(frame, PRESCAN, [1, 1], GAIN)

    # ASSERT: THE TRUE SCATTER IS 2.5 ADU, PLUS 1/12 ADU^2 OF ROUNDING NOISE
    assert ron == pytest.approx(np.sqrt(2.5**2 + 1.0 / 12.0) * 1.1, rel=TOLERANCE)


def test_the_prescan_window_follows_the_binning() -> None:
    # ARRANGE: 2 X 2 BINNING HALVES THE WINDOW ON BOTH AXES, AS THE TRIM OF THE SCIENCE PIXELS DOES
    frame = _raw_frame()[:200, :60].copy()
    frame[:, 4:23] = np.random.default_rng(2).normal(590.0, 2.0, (200, 19))

    # ACT
    ron = measure_prescan_read_noise(frame, PRESCAN, [2, 2], 1.0 * u.electron / u.adu)

    # ASSERT
    assert ron == pytest.approx(2.0, rel=TOLERANCE)


@pytest.mark.parametrize(
    ("window", "reason"),
    [
        ({"rows": {"start": 50, "end": 350}, "columns": {"start": 200, "end": 240}}, "outside the frame"),
        ({"rows": {"start": 50, "end": 52}, "columns": {"start": 8, "end": 10}}, "too few"),
    ],
)
def test_a_prescan_the_frame_does_not_hold_cannot_be_measured(window: dict, reason: str) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match=reason):
        measure_prescan_read_noise(_raw_frame(), window, [1, 1], GAIN)


def test_a_prescan_with_no_scatter_cannot_be_measured() -> None:
    # ARRANGE
    frame = _raw_frame()
    frame[:, 8:46] = 590.0

    # ACT / ASSERT
    with pytest.raises(ValueError, match="not a positive number"):
        measure_prescan_read_noise(frame, PRESCAN, [1, 1], GAIN)


def test_a_prescan_of_nan_cannot_be_measured() -> None:
    # ARRANGE
    frame = _raw_frame()
    frame[:, 8:46] = np.nan

    # ACT / ASSERT
    with pytest.raises(ValueError, match="too few"):
        measure_prescan_read_noise(frame, PRESCAN, [1, 1], GAIN)


@pytest.mark.parametrize(
    ("arm", "recipeName", "isExpected"),
    [
        ("VIS", "soxs-stare", True),
        ("VIS", "soxs-stare-std", True),
        ("VIS", "soxs-nod", True),
        ("VIS", "soxs-offset-std", True),
        ("UVB", "soxs-stare", True),
        ("NIR", "soxs-stare", False),
        ("VIS", "soxs-mflat", False),
        ("VIS", "soxs-mbias", False),
        ("VIS", "soxs-disp-solution", False),
    ],
)
def test_only_the_uvb_and_vis_science_recipes_use_the_prescan(arm: str, recipeName: str, isExpected: bool) -> None:
    # ARRANGE
    detectorParams = {"prescan-pixels": PRESCAN}

    # ACT / ASSERT
    assert uses_prescan_read_noise(arm, recipeName, detectorParams) is isExpected


def test_a_detector_with_no_prescan_window_does_not_use_the_prescan() -> None:
    # ACT / ASSERT
    assert uses_prescan_read_noise("VIS", "soxs-stare", {}) is False


def test_the_science_read_noise_is_the_root_mean_square_of_the_frames() -> None:
    # ARRANGE
    master = ReadNoise(singleFrame=3.83, nCombined=4, skyFit=3.83 * np.sqrt(1.25), source="master bias")

    # ACT
    result = science_read_noise([4.0, 6.0], master)

    # ASSERT
    expected = float(np.sqrt((16.0 + 36.0) / 2.0))
    assert result.singleFrame == pytest.approx(expected)
    assert result.source == SOURCE_SCIENCE_PRESCAN
    assert result.nCombined == 4
    assert result.skyFit == pytest.approx(expected * np.sqrt(1.25))


def test_the_science_read_noise_without_a_master_bias_is_also_the_sky_fit_noise() -> None:
    # ACT
    result = science_read_noise([5.1], None)

    # ASSERT
    assert result.singleFrame == pytest.approx(5.1)
    assert result.nCombined is None
    assert result.skyFit == pytest.approx(5.1)
