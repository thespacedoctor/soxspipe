"""Contracts for the robust set-level offset test on multi-pinhole arc-line sets (DY-1526, DY-2170)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.pinhole_set_offsets import (
    MIN_SETS_PER_ORDER,
    OFFSET_SET_CLIP_SIGMA,
    OFFSET_SET_MIN_SIGMA,
    flag_offset_pinhole_sets,
)
from tests.unit._pinhole_sets import NOISE_OFFSETS, PINHOLES

pytestmark = pytest.mark.unit

WITHIN_SET_SCATTER = 0.1


def _sets(order: int, offsetsY, offsetsX=None, firstWavelength: float = 500.0) -> pd.DataFrame:
    """One row per (set, pinhole); every pinhole of a set shares the set's offset plus a fixed within-set pattern."""
    offsetsY = np.asarray(offsetsY, dtype=float)
    offsetsX = np.zeros_like(offsetsY) if offsetsX is None else np.asarray(offsetsX, dtype=float)
    pattern = np.linspace(-WITHIN_SET_SCATTER, WITHIN_SET_SCATTER, PINHOLES)
    rows = []
    for setIndex, (offsetX, offsetY) in enumerate(zip(offsetsX, offsetsY, strict=True)):
        for pinhole in range(PINHOLES):
            rows.append(
                {
                    "wavelength": firstWavelength + setIndex,
                    "order": order,
                    "slit_index": pinhole,
                    "residuals_x": offsetX + pattern[pinhole],
                    "residuals_y": offsetY + pattern[::-1][pinhole],
                }
            )
    return pd.DataFrame(rows)


def _flagged_wavelengths(lines: pd.DataFrame, flagged: np.ndarray) -> list[float]:
    return sorted(set(lines.loc[flagged, "wavelength"]))


def test_drops_exactly_the_set_offset_in_every_pinhole() -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=[*NOISE_OFFSETS, 2.0])

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert _flagged_wavelengths(lines, flagged) == [511.0]
    assert flagged.sum() == PINHOLES


def test_drops_a_set_offset_in_x_only() -> None:
    # ARRANGE
    offsetsX = [*NOISE_OFFSETS, -2.0]
    lines = _sets(order=1, offsetsY=[*NOISE_OFFSETS, 0.0], offsetsX=offsetsX)

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert _flagged_wavelengths(lines, flagged) == [511.0]


def test_drops_nothing_when_set_offsets_are_consistent_with_noise() -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=NOISE_OFFSETS, offsetsX=NOISE_OFFSETS[::-1])

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert not flagged.any()


def test_skips_an_order_with_too_few_sets() -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=[0.0, 0.1, -0.1, 2.0])
    assert lines["wavelength"].nunique() < MIN_SETS_PER_ORDER

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert not flagged.any()


def test_a_smaller_sigma_limit_moves_the_cut() -> None:
    # ARRANGE: A +1 PX SET IS ABOUT 3.7 MAD-SIGMA OUT
    lines = _sets(order=1, offsetsY=[*NOISE_OFFSETS, 1.0])

    # ACT
    flaggedDefault = flag_offset_pinhole_sets(lines, sigmaLimit=OFFSET_SET_CLIP_SIGMA)
    flaggedTight = flag_offset_pinhole_sets(lines, sigmaLimit=3)

    # ASSERT
    assert not flaggedDefault.any()
    assert _flagged_wavelengths(lines, flaggedTight) == [511.0]


def test_each_order_is_judged_against_its_own_spread() -> None:
    # ARRANGE: THE SAME +1 PX SET IS AN OUTLIER IN A TIGHT ORDER BUT NOT IN A BROAD ONE
    tight = _sets(order=1, offsetsY=[*(NOISE_OFFSETS / 10), 1.0])
    broad = _sets(order=2, offsetsY=[*(NOISE_OFFSETS * 2), 1.0], firstWavelength=600.0)
    lines = pd.concat([tight, broad], ignore_index=True)

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert _flagged_wavelengths(lines, flagged) == [511.0]


def test_offsets_are_measured_from_the_order_median() -> None:
    # ARRANGE: THE WHOLE ORDER SITS +0.5 PX OFF THE MODEL, AND ONE SET A FURTHER +2 PX
    lines = _sets(order=1, offsetsY=[*(NOISE_OFFSETS + 0.5), 2.5])

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert _flagged_wavelengths(lines, flagged) == [511.0]


def test_single_pinhole_lines_are_never_flagged() -> None:
    # ARRANGE: ONE MEASUREMENT PER LINE, SO THERE ARE NO SETS
    lines = pd.DataFrame(
        {
            "wavelength": 500.0 + np.arange(12),
            "order": 1,
            "slit_index": 0,
            "residuals_x": 0.0,
            "residuals_y": [*NOISE_OFFSETS, 5.0],
        }
    )

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert not flagged.any()


def test_an_order_with_no_spread_is_skipped() -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=[0.0] * 11 + [2.0])

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert not flagged.any()


def test_the_result_is_positional_and_leaves_the_table_unchanged() -> None:
    # ARRANGE: A SHUFFLED, NON-DEFAULT INDEX
    lines = _sets(order=1, offsetsY=[*NOISE_OFFSETS, 2.0]).sample(frac=1, random_state=3)
    lines.index = lines.index + 1000
    before = lines.copy()

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert isinstance(flagged, np.ndarray)
    assert flagged.dtype == bool
    assert np.array_equal(flagged, (lines["wavelength"] == 511.0).to_numpy())
    pd.testing.assert_frame_equal(lines, before)


def test_an_empty_table_flags_nothing() -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=[]).reindex(columns=["wavelength", "order", "residuals_x", "residuals_y"])

    # ACT
    flagged = flag_offset_pinhole_sets(lines)

    # ASSERT
    assert flagged.shape == (0,)


# DY-2170: A MAD-SIGMA FLOOR MAKES A SET'S DECISION INDEPENDENT OF HOW TIGHTLY THE OTHER SETS IN ITS ORDER SCATTER
MIN_SIGMA = 0.1
# THE OTHER SETS SCATTER WITH A MAD-SIGMA OF ABOUT 0.08 PX, BELOW THE FLOOR
TIGHT_OFFSETS = NOISE_OFFSETS * 0.3


@pytest.mark.parametrize(("candidateOffset", "isClipped"), [(0.3, False), (0.8, True)])
@pytest.mark.parametrize("axis", ["x", "y"])
def test_a_set_decision_does_not_change_when_the_other_sets_scatter_halves(
    candidateOffset: float, isClipped: bool, axis: str
) -> None:
    # ARRANGE
    def lines_for(scatter: np.ndarray) -> pd.DataFrame:
        offsets = [*scatter, candidateOffset]
        zeros = [0.0] * len(offsets)
        if axis == "x":
            return _sets(order=1, offsetsY=zeros, offsetsX=offsets)
        return _sets(order=1, offsetsY=offsets)

    wide = lines_for(TIGHT_OFFSETS)
    halved = lines_for(TIGHT_OFFSETS / 2)

    # ACT
    flaggedWide = flag_offset_pinhole_sets(wide, minSigma=MIN_SIGMA)
    flaggedHalved = flag_offset_pinhole_sets(halved, minSigma=MIN_SIGMA)

    # ASSERT
    assert (511.0 in _flagged_wavelengths(wide, flaggedWide)) is isClipped
    assert (511.0 in _flagged_wavelengths(halved, flaggedHalved)) is isClipped


def test_without_a_floor_the_decision_follows_the_other_sets_scatter() -> None:
    # ARRANGE: THE DY-2170 SYMPTOM. +0.3 PX IS INSIDE 5 MAD-SIGMA OF THE WIDE ORDER, OUTSIDE IT ONCE THE SCATTER HALVES
    wide = _sets(order=1, offsetsY=[*TIGHT_OFFSETS, 0.3])
    halved = _sets(order=1, offsetsY=[*(TIGHT_OFFSETS / 2), 0.3])

    # ACT
    flaggedWide = flag_offset_pinhole_sets(wide)
    flaggedHalved = flag_offset_pinhole_sets(halved)

    # ASSERT
    assert not flaggedWide.any()
    assert _flagged_wavelengths(halved, flaggedHalved) == [511.0]


def test_a_set_offset_in_every_pinhole_is_still_clipped_with_a_floor() -> None:
    # ARRANGE: THE DY-1526 CASE, A SET OFF BY +1 PX IN EVERY PINHOLE OF A TIGHT ORDER
    lines = _sets(order=1, offsetsY=[*TIGHT_OFFSETS, 1.0])

    # ACT
    flagged = flag_offset_pinhole_sets(lines, minSigma=MIN_SIGMA)

    # ASSERT
    assert _flagged_wavelengths(lines, flagged) == [511.0]
    assert flagged.sum() == PINHOLES


def test_an_order_broader_than_the_floor_keeps_its_own_spread() -> None:
    # ARRANGE: THE MAD-SIGMA (ABOUT 0.27 PX) IS ABOVE THE FLOOR, SO THE CUT STAYS AT 5 MAD-SIGMA
    keptSet = _sets(order=1, offsetsY=[*NOISE_OFFSETS, 1.0])
    clippedSet = _sets(order=2, offsetsY=[*NOISE_OFFSETS, 2.0], firstWavelength=600.0)
    lines = pd.concat([keptSet, clippedSet], ignore_index=True)

    # ACT
    flagged = flag_offset_pinhole_sets(lines, minSigma=MIN_SIGMA)

    # ASSERT
    assert _flagged_wavelengths(lines, flagged) == [611.0]


def test_a_floor_gives_an_order_with_no_spread_a_scale() -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=[0.0] * 11 + [2.0])

    # ACT
    flagged = flag_offset_pinhole_sets(lines, minSigma=MIN_SIGMA)

    # ASSERT
    assert _flagged_wavelengths(lines, flagged) == [511.0]


def test_the_default_floor_keeps_the_mad_only_rule() -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=[*(TIGHT_OFFSETS / 2), 0.3])

    # ACT
    flaggedDefault = flag_offset_pinhole_sets(lines)
    flaggedExplicit = flag_offset_pinhole_sets(lines, minSigma=OFFSET_SET_MIN_SIGMA)

    # ASSERT
    assert OFFSET_SET_MIN_SIGMA == 0.0
    assert np.array_equal(flaggedDefault, flaggedExplicit)


@pytest.mark.parametrize("badFloor", [-0.1, float("inf"), float("nan")])
def test_an_invalid_floor_is_rejected(badFloor: float) -> None:
    # ARRANGE
    lines = _sets(order=1, offsetsY=[*NOISE_OFFSETS, 2.0])

    # ACT / ASSERT
    with pytest.raises(ValueError, match="minSigma"):
        flag_offset_pinhole_sets(lines, minSigma=badFloor)
