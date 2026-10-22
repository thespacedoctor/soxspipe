"""Contracts for the robust set-level offset test on multi-pinhole arc-line sets (DY-1526)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.pinhole_set_offsets import (
    MIN_SETS_PER_ORDER,
    OFFSET_SET_CLIP_SIGMA,
    flag_offset_pinhole_sets,
)

pytestmark = pytest.mark.unit

PINHOLES = 9
WITHIN_SET_SCATTER = 0.1
# ELEVEN SET OFFSETS SPREAD EVENLY OVER +-0.3 PX: A MAD-SIGMA OF ABOUT 0.27 PX
NOISE_OFFSETS = np.linspace(-0.3, 0.3, 11)


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
