"""Slit position of an object trace in each order, read from the 2D map (DY-1612)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.slit_alignment import (
    MIN_TRACE_SAMPLES,
    build_two_d_map_images,
    measure_trace_slit_positions,
    order_slit_positions,
)
from tests.factories import pipeline_settings
from tests.factories.slit_alignment import VIS_ORDERS, synthetic_alignment_map, synthetic_nod_frame
from tests.unit.test_slit_alignment_frames import DETECT_CONTINUUM_SETTINGS

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


@pytest.mark.parametrize("sampled", [(None, 0.0), (pd.DataFrame({"order": [1.0]}), 5.0)])
def test_a_frame_whose_trace_is_not_found_gives_every_order_no_position(log, tmp_path, monkeypatch, sampled) -> None:
    # ARRANGE: THE TRACE SAMPLING FINDS NOTHING, OR FINDS IT IN TOO FEW SAMPLES TO BE TRUSTED
    scene = synthetic_alignment_map(tmp_path)
    mapImages = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    monkeypatch.setattr(detect_continuum, "sample_trace", lambda self: sampled)
    frame = synthetic_nod_frame(scene, objectSlit=1.0, mjd=60000.0)

    # ACT
    positions = measure_trace_slit_positions(
        log=log,
        differenceFrame=frame,
        mapImages=mapImages,
        dispersionMap="unused.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"detect-continuum": DETECT_CONTINUUM_SETTINGS},
        recipeName="soxs-nod",
    )

    # ASSERT
    assert positions["order"].tolist() == list(scene.orders)
    assert positions["position_arcsec"].isna().all()
