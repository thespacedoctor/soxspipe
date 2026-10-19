"""Contracts for merging unresolved OH doublets in the sky-line catalogue into flux-weighted centroids (DY-1645)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.toolkit import merge_blended_skylines

pytestmark = pytest.mark.unit


def _catalogue(*lines: tuple[float, float]) -> pd.DataFrame:
    return pd.DataFrame({"WAVELENGTH": [w for w, _ in lines], "FLUX": [f for _, f in lines], "ISOLATED": False})


def test_a_doublet_merges_into_one_flux_weighted_centroid() -> None:
    # ARRANGE
    catalogue = _catalogue((785.32587, 8.526), (785.35240, 10.799))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    expected = (785.32587 * 8.526 + 785.35240 * 10.799) / (8.526 + 10.799)
    assert len(merged) == 1
    assert merged["WAVELENGTH"].iloc[0] == pytest.approx(expected)
    assert merged["FLUX"].iloc[0] == pytest.approx(8.526 + 10.799)
    assert list(merged.columns) == ["WAVELENGTH", "FLUX"]


def test_chained_lines_each_closer_than_the_gap_merge_into_one_group() -> None:
    # ARRANGE
    catalogue = _catalogue((600.00, 3.0), (600.08, 3.0), (600.16, 3.0))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert len(merged) == 1
    assert merged["WAVELENGTH"].iloc[0] == pytest.approx(600.08)
    assert merged["FLUX"].iloc[0] == pytest.approx(9.0)


def test_lines_further_apart_than_the_gap_stay_separate() -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, 20.0), (601.0, 30.0))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged["WAVELENGTH"].tolist() == pytest.approx([600.0, 601.0])


def test_a_group_with_a_comparable_neighbour_inside_the_isolation_window_is_dropped() -> None:
    # ARRANGE
    # EACH GROUP HAS MORE THAN 20 PERCENT OF THE OTHER'S FLUX 0.2 NM AWAY, SO NEITHER IS ISOLATED
    catalogue = _catalogue((600.0, 20.0), (600.2, 25.0), (610.0, 20.0))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged["WAVELENGTH"].tolist() == pytest.approx([610.0])


def test_a_faint_neighbour_below_the_flux_ratio_does_not_reject_a_bright_line() -> None:
    # ARRANGE
    # THE FAINT LINE HAS 11 PERCENT OF THE BRIGHT FLUX: IT CANNOT REJECT THE BRIGHT LINE, BUT IS REJECTED BY IT
    catalogue = _catalogue((600.0, 50.0), (600.2, 5.5))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged["WAVELENGTH"].tolist() == pytest.approx([600.0])


def test_a_group_whose_total_flux_is_too_faint_is_dropped() -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, 2.0), (600.05, 2.0), (610.0, 30.0))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged["WAVELENGTH"].tolist() == pytest.approx([610.0])


def test_faint_components_below_the_line_floor_neither_weight_nor_reject() -> None:
    # ARRANGE
    # A 0.4 FLUX LINE IS IGNORED, SO IT CANNOT PULL THE CENTROID OR ACT AS A NEIGHBOUR
    catalogue = _catalogue((600.0, 30.0), (600.05, 0.4), (600.2, 0.4))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged["WAVELENGTH"].tolist() == pytest.approx([600.0])
    assert merged["FLUX"].tolist() == pytest.approx([30.0])


def test_the_input_table_is_not_modified_and_an_empty_result_keeps_its_columns() -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, 1.0))
    before = catalogue.copy()

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    pd.testing.assert_frame_equal(catalogue, before)
    assert merged.empty
    assert list(merged.columns) == ["WAVELENGTH", "FLUX"]


def test_a_catalogue_with_no_line_above_the_floor_gives_an_empty_table() -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, 0.3), (601.0, 0.4))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged.empty
    assert list(merged.columns) == ["WAVELENGTH", "FLUX"]


def test_non_finite_rows_are_dropped_before_grouping() -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, 30.0), (600.05, np.nan), (np.nan, 40.0), (np.inf, 40.0), (610.0, np.inf))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged["WAVELENGTH"].tolist() == pytest.approx([600.0])
    assert merged["FLUX"].tolist() == pytest.approx([30.0])


# EACH THRESHOLD IS PROBED JUST BELOW, AT AND JUST ABOVE ITS VALUE. LINES USE EXACTLY REPRESENTABLE VALUES WHERE A
# COMPARISON IS "AT" THE THRESHOLD


@pytest.mark.parametrize(
    ("gapNm", "expectedRows"),
    [
        pytest.param(0.099, 1, id="just-below-merges"),
        pytest.param(0.125, 0, id="at-stays-separate"),
        pytest.param(0.101, 0, id="just-above-stays-separate"),
    ],
)
def test_merge_gap_threshold(gapNm: float, expectedRows: int) -> None:
    # ARRANGE
    # TWO FLUX-6 LINES: MERGED THEY HOLD FLUX 12 AND ARE KEPT, APART THEY REJECT EACH OTHER AS EQUAL NEIGHBOURS
    catalogue = _catalogue((600.0, 6.0), (600.0 + gapNm, 6.0))
    mergeNm = 0.125 if gapNm == 0.125 else 0.1

    # ACT
    merged = merge_blended_skylines(catalogue, mergeNm=mergeNm)

    # ASSERT
    assert len(merged) == expectedRows


@pytest.mark.parametrize(
    ("separationNm", "isolationNm", "expectedRows"),
    [
        pytest.param(0.349, 0.35, 0, id="just-inside-rejects"),
        pytest.param(0.5, 0.5, 2, id="at-is-isolated"),
        pytest.param(0.351, 0.35, 2, id="just-outside-is-isolated"),
    ],
)
def test_isolation_radius_threshold(separationNm: float, isolationNm: float, expectedRows: int) -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, 20.0), (600.0 + separationNm, 20.0))

    # ACT
    merged = merge_blended_skylines(catalogue, isolationNm=isolationNm)

    # ASSERT
    assert len(merged) == expectedRows


@pytest.mark.parametrize(
    ("neighbourFlux", "brightLineKept"),
    [
        pytest.param(9.9, True, id="just-below-the-ratio-keeps"),
        pytest.param(10.0, True, id="at-the-ratio-keeps"),
        pytest.param(10.1, False, id="just-above-the-ratio-rejects"),
    ],
)
def test_neighbour_flux_ratio_threshold(neighbourFlux: float, brightLineKept: bool) -> None:
    # ARRANGE
    # THE RATIO IS 0.2, SO A NEIGHBOUR NEEDS MORE THAN FLUX 10 TO REJECT THE FLUX-50 LINE
    catalogue = _catalogue((600.0, 50.0), (600.2, neighbourFlux))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert (50.0 in merged["FLUX"].tolist()) is brightLineKept


@pytest.mark.parametrize(
    ("companionFlux", "companionCounts"),
    [
        pytest.param(0.49, False, id="just-below-the-floor-is-ignored"),
        pytest.param(0.5, False, id="at-the-floor-is-ignored"),
        pytest.param(0.51, True, id="just-above-the-floor-takes-part"),
    ],
)
def test_line_flux_floor(companionFlux: float, companionCounts: bool) -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, 20.0), (600.05, companionFlux))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert len(merged) == 1
    assert bool(merged["FLUX"].iloc[0] > 20.0) == companionCounts


@pytest.mark.parametrize(
    ("groupFlux", "expectedRows"),
    [
        pytest.param(4.99, 0, id="just-below-is-dropped"),
        pytest.param(5.0, 0, id="at-is-dropped"),
        pytest.param(5.01, 1, id="just-above-is-kept"),
    ],
)
def test_group_flux_floor(groupFlux: float, expectedRows: int) -> None:
    # ARRANGE
    catalogue = _catalogue((600.0, groupFlux))

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert len(merged) == expectedRows


def test_missing_and_non_numeric_values_are_dropped_before_grouping() -> None:
    # ARRANGE
    catalogue = pd.DataFrame(
        {
            "WAVELENGTH": pd.Series([600.0, 600.05, pd.NA, "bad", 610.0], dtype=object),
            "FLUX": pd.Series([30.0, pd.NA, 40.0, 40.0, "x"], dtype=object),
        }
    )

    # ACT
    merged = merge_blended_skylines(catalogue)

    # ASSERT
    assert merged["WAVELENGTH"].tolist() == pytest.approx([600.0])
    assert merged["FLUX"].tolist() == pytest.approx([30.0])
