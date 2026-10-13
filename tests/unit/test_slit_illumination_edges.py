"""Contracts for the master flat's slit edges, fitted per order by ``slit_illumination.flat_slit_edges``."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils import slit_illumination
from soxspipe.commonutils.slit_illumination import flat_slit_edges
from tests.factories import synthetic_illumination_scene

pytestmark = pytest.mark.unit

FLAT_LOWER = -5.0
FLAT_UPPER = 5.0
FLAT_WIDTH = 0.15


@pytest.fixture
def scene(tmp_path):
    """Two orders lit between -5 and +5 arcsec, with 0.15 arcsec edges."""
    return synthetic_illumination_scene(
        tmp_path,
        orders=(11, 12),
        flatEdges=(FLAT_LOWER, FLAT_UPPER, FLAT_WIDTH),
        hasNoise=False,
    )


def test_recovers_lower_upper_and_width_of_every_order(scene):
    edges = flat_slit_edges(scene.masterFlat.data, scene.slit, scene.orderMap)

    assert list(edges["order"]) == [11, 12]
    assert edges["lower"].to_numpy() == pytest.approx([FLAT_LOWER, FLAT_LOWER], abs=0.02)
    assert edges["upper"].to_numpy() == pytest.approx([FLAT_UPPER, FLAT_UPPER], abs=0.02)
    assert edges["width"].to_numpy() == pytest.approx([FLAT_WIDTH, FLAT_WIDTH], abs=0.03)
    assert (edges["reason"] == "").all()


def test_recovers_edges_that_differ_between_orders(tmp_path):
    flatData, slit, orderMap = _two_order_flat(
        {11: (-5.0, 4.8, 0.12), 12: (-4.6, 5.2, 0.2)},
    )

    edges = flat_slit_edges(flatData, slit, orderMap).set_index("order")

    assert edges.loc[11, "lower"] == pytest.approx(-5.0, abs=0.03)
    assert edges.loc[11, "upper"] == pytest.approx(4.8, abs=0.03)
    assert edges.loc[12, "lower"] == pytest.approx(-4.6, abs=0.03)
    assert edges.loc[12, "upper"] == pytest.approx(5.2, abs=0.03)
    assert edges.loc[12, "width"] == pytest.approx(0.2, abs=0.04)


def test_does_not_modify_its_inputs(scene):
    flatData = scene.masterFlat.data.copy()
    slit = scene.slit.copy()
    orderMap = scene.orderMap.copy()

    flat_slit_edges(flatData, slit, orderMap)

    np.testing.assert_array_equal(flatData, scene.masterFlat.data)
    np.testing.assert_array_equal(slit, scene.slit)
    np.testing.assert_array_equal(orderMap, scene.orderMap)


def test_marks_an_order_with_too_few_populated_bins_as_failed_with_a_reason(scene):
    orderMap = scene.orderMap.copy()
    sparse = np.zeros(orderMap.shape, dtype=bool)
    sparse[:3, :] = True
    orderMap[(orderMap == 12) & ~sparse] = np.nan

    edges = flat_slit_edges(scene.masterFlat.data, scene.slit, orderMap).set_index("order")

    assert edges.loc[11, "reason"] == ""
    assert np.isnan(edges.loc[12, "lower"])
    assert "too few" in edges.loc[12, "reason"]


def test_marks_a_dark_flat_as_failed_with_a_reason(scene):
    dark = np.zeros_like(scene.masterFlat.data)

    edges = flat_slit_edges(dark, scene.slit, scene.orderMap)

    assert edges["lower"].isna().all()
    assert edges["reason"].str.contains("flat level").all()


def test_ignores_pixels_outside_every_order(scene):
    edges = flat_slit_edges(scene.masterFlat.data, scene.slit, scene.orderMap)

    assert set(edges["order"]) == {11, 12}


def test_recovers_the_edges_with_a_coarser_bin(scene):
    edges = flat_slit_edges(scene.masterFlat.data, scene.slit, scene.orderMap, binArcsec=0.3)

    assert edges["lower"].to_numpy() == pytest.approx([FLAT_LOWER, FLAT_LOWER], abs=0.05)
    assert edges["upper"].to_numpy() == pytest.approx([FLAT_UPPER, FLAT_UPPER], abs=0.05)


def test_marks_an_order_as_failed_when_the_edge_fit_does_not_converge(scene, monkeypatch):
    def refuse(*args, **kwargs):
        raise RuntimeError("Optimal parameters not found")

    monkeypatch.setattr(slit_illumination, "curve_fit", refuse)

    edges = flat_slit_edges(scene.masterFlat.data, scene.slit, scene.orderMap)

    assert edges["lower"].isna().all()
    assert edges["reason"].str.contains("did not converge").all()


def test_marks_an_order_as_failed_when_the_fitted_edges_are_not_ordered(scene, monkeypatch):
    monkeypatch.setattr(slit_illumination, "curve_fit", lambda *args, **kwargs: (np.array([1.0, 2.0, -2.0, 0.1]), None))

    edges = flat_slit_edges(scene.masterFlat.data, scene.slit, scene.orderMap)

    assert edges["lower"].isna().all()
    assert edges["reason"].str.contains("not ordered").all()


def _two_order_flat(edgesByOrder):
    """Build a bare flat, slit and order image for two orders with their own edges."""
    from tests.factories.slit_illumination import top_hat

    rng = np.random.default_rng(5)
    slit = rng.uniform(-5.4, 5.4, size=(400, 2 * 40))
    orderMap = np.where(np.arange(slit.shape[1]) < 40, 11.0, 12.0)[None, :] * np.ones((400, 1))
    flat = np.zeros(slit.shape)
    for order, (lower, upper, width) in edgesByOrder.items():
        isOrder = orderMap == order
        flat[isOrder] = top_hat(slit[isOrder], 1.0, lower, upper, width)
    return flat, slit, orderMap
