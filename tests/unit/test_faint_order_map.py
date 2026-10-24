"""Reading the dispersion map as frame pixels for the faint-order trace (DY-1941)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from soxspipe.commonutils.faint_order_map import make_pixels_at
from tests.unit._continuum_vis_scene import vis_dispersion_map, vis_row, vis_slit_scale

pytestmark = pytest.mark.unit


def test_make_pixels_at_reads_the_dispersion_map_in_frame_pixels(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE FRAME IS BINNED 2 IN y RELATIVE TO THE MAP, AND THE SLIT RUNS ALONG x
    mapPath = vis_dispersion_map(tmp_path)
    pixelsAt = make_pixels_at(log, mapPath, {"x": 1.0, "y": 2.0}, "x", "y")

    # ACT
    axisA, axisB = pixelsAt(3.0, np.array([610.0, 630.0]), 0.5)

    # ASSERT: x = -30 + 90 * 3 + 0.5 * SCALE(3) AND y = ROW(WAVELENGTH) / 2
    np.testing.assert_allclose(axisA, 240.0 + 0.5 * vis_slit_scale(3.0), atol=1e-3)
    np.testing.assert_allclose(axisB, vis_row(3.0, np.array([610.0, 630.0])) / 2.0, atol=1e-3)


def test_make_pixels_at_keeps_samples_that_fall_off_the_detector(log: object, tmp_path: Path) -> None:
    # ARRANGE: A SLIT POSITION 100 ARCSEC OUT PUTS x BELOW ZERO
    pixelsAt = make_pixels_at(log, vis_dispersion_map(tmp_path), {"x": 1.0, "y": 1.0}, "x", "y")

    # ACT
    axisA, axisB = pixelsAt(1.0, np.array([810.0, 820.0]), 100.0)

    # ASSERT
    assert len(axisA) == len(axisB) == 2
    assert (axisA < 0).all()


def test_make_pixels_at_follows_the_axis_roles_it_is_given(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE SAME MAP WITH THE ROLES OF THE AXES SWAPPED
    pixelsAt = make_pixels_at(log, vis_dispersion_map(tmp_path), {"x": 1.0, "y": 1.0}, "y", "x")

    # ACT
    axisA, axisB = pixelsAt(2.0, np.array([410.0]), 0.0)

    # ASSERT
    np.testing.assert_allclose(axisA, vis_row(2.0, np.array([410.0])), atol=1e-3)
    np.testing.assert_allclose(axisB, [150.0], atol=1e-3)
