"""Slit position of an object trace in each order, read from the 2D map (DY-1612)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.slit_alignment import (
    MIN_TRACE_SAMPLES,
    build_two_d_map_images,
    order_slit_positions,
)
from tests.factories.slit_alignment import VIS_ORDERS, synthetic_alignment_map

pytestmark = pytest.mark.unit


def _trace_samples(scene, order: int, slitPosition: float, rows: np.ndarray) -> pd.DataFrame:
    """Return trace samples of ``order`` that sit at ``slitPosition`` arcsec, one per row of ``rows``."""
    columns = []
    for row in rows:
        onOrder = scene.orderMap[row] == order
        columnIndices = np.flatnonzero(onOrder)
        columns.append(np.interp(slitPosition, scene.slit[row][onOrder], columnIndices))
    return pd.DataFrame({"order": float(order), "cont_x": columns, "cont_y": rows.astype(float), "pre-clipped": False})


def test_an_order_position_is_the_median_slit_coordinate_of_its_trace_samples(tmp_path) -> None:
    # ARRANGE: THE TRACE OF ORDER 1 SITS AT 1.7 ARCSEC AND THE TRACE OF ORDER 3 AT -2.3 ARCSEC
    scene = synthetic_alignment_map(tmp_path)
    mapImages = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    rows = np.arange(40, 560, 13)
    samples = pd.concat([_trace_samples(scene, 1, 1.7, rows), _trace_samples(scene, 3, -2.3, rows)])

    # ACT
    positions = order_slit_positions(samples, mapImages).set_index("order")

    # ASSERT
    assert positions.loc[1, "position_arcsec"] == pytest.approx(1.7, abs=0.01)
    assert positions.loc[3, "position_arcsec"] == pytest.approx(-2.3, abs=0.01)
    assert positions.loc[1, "n_samples"] == len(rows)
    assert positions.index.tolist() == list(VIS_ORDERS)


def test_an_order_with_too_few_samples_has_no_position(tmp_path) -> None:
    # ARRANGE
    scene = synthetic_alignment_map(tmp_path)
    mapImages = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    samples = _trace_samples(scene, 2, 1.0, np.arange(40, 40 + MIN_TRACE_SAMPLES - 1))

    # ACT
    positions = order_slit_positions(samples, mapImages).set_index("order")

    # ASSERT
    assert np.isnan(positions.loc[2, "position_arcsec"])
    assert positions.loc[2, "n_samples"] == MIN_TRACE_SAMPLES - 1
    assert np.isnan(positions.loc[4, "position_arcsec"])
    assert positions.loc[4, "n_samples"] == 0


def test_an_order_whose_samples_scatter_widely_has_no_position(tmp_path) -> None:
    # ARRANGE: HALF THE SAMPLES SIT AT -3 ARCSEC AND HALF AT +3 ARCSEC, SO THERE IS NO TRACE TO FOLLOW
    scene = synthetic_alignment_map(tmp_path)
    mapImages = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    rows = np.arange(40, 560, 13)
    samples = pd.concat([_trace_samples(scene, 1, -3.0, rows), _trace_samples(scene, 1, 3.0, rows)])

    # ACT
    positions = order_slit_positions(samples, mapImages).set_index("order")

    # ASSERT
    assert np.isnan(positions.loc[1, "position_arcsec"])
    assert positions.loc[1, "scatter_arcsec"] > 1.0


def test_clipped_and_off_order_samples_do_not_count(tmp_path) -> None:
    # ARRANGE: THE TRACE IS AT 1.0 ARCSEC. A SECOND SET OF SAMPLES AT 4 ARCSEC IS CLIPPED, AND A THIRD SET SITS
    # IN THE GAP BETWEEN THE ORDERS
    scene = synthetic_alignment_map(tmp_path)
    mapImages = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    rows = np.arange(40, 560, 13)
    clipped = _trace_samples(scene, 1, 4.0, rows).assign(**{"pre-clipped": True})
    gap = _trace_samples(scene, 1, 1.0, rows).assign(cont_x=75.0)
    samples = pd.concat([_trace_samples(scene, 1, 1.0, rows), clipped, gap])

    # ACT
    positions = order_slit_positions(samples, mapImages).set_index("order")

    # ASSERT
    assert positions.loc[1, "position_arcsec"] == pytest.approx(1.0, abs=0.01)
    assert positions.loc[1, "n_samples"] == len(rows)
