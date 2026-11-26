"""*tests for the line measurement that the SPEC_RES script of tools/ uses*

The measurement runs on synthetic frames: a wavelength map that is linear along the dispersion axis, a slit map in
arcsec across it, and Gaussian emission lines of a known FWHM that fill the slit. The resolving power it returns
must match `centre / FWHM` of the lines that were put in.

The maps have the layout of the real `DISP_IMAGE` files of each arm (DY-2455): the wavelength runs along numpy axis 0
for VIS and numpy axis 1 for NIR. The detector parameters name those arms `x` and `y`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tests.unit._arc_scene import LINE_COLUMNS, ROWS, Scene
from tests.unit._tools_loader import load_tool

line_fitting = load_tool("line_fitting")

pytestmark = pytest.mark.unit

DETECTOR_DISPERSION_AXIS = {"VIS": "x", "NIR": "y"}


def _measure(scene: Scene, **overrides: Any) -> list:
    """*measure a scene's lines with the default arguments, replacing the ones given*"""
    arguments = {
        "data": scene.data,
        "good": scene.good,
        "wavelength": scene.wavelength,
        "slit": scene.slit,
        "order": scene.order,
        "dispersionAxis": DETECTOR_DISPERSION_AXIS[scene.arm],
    }
    arguments.update(overrides)
    return line_fitting.measure_lines(**arguments)


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_measure_lines_recovers_the_resolving_power_of_isolated_lines(arm: str) -> None:
    # ARRANGE
    scene = Scene(arm=arm, fwhmPx=3.0)

    # ACT
    lines = _measure(scene)

    # ASSERT
    assert len(lines) == len(LINE_COLUMNS)
    for line in lines:
        assert line.order == 1
        assert line.fwhmPx == pytest.approx(3.0, rel=0.03)
        assert line.resolvingPower == pytest.approx(scene.expected_resolving_power(line.centreNm), rel=0.03)


def test_measure_lines_reports_a_wider_line_as_a_lower_resolving_power() -> None:
    # ARRANGE
    narrowScene, wideScene = Scene(fwhmPx=2.5), Scene(fwhmPx=5.0)

    # ACT
    narrow = np.median([line.resolvingPower for line in _measure(narrowScene)])
    wide = np.median([line.resolvingPower for line in _measure(wideScene)])

    # ASSERT
    assert narrow == pytest.approx(2 * wide, rel=0.06)


def test_measure_lines_leaves_out_a_line_with_a_bad_pixel_in_its_fit_window() -> None:
    # ARRANGE
    scene = Scene()
    scene.good[scene.pixel_at(scene.centres[3] + 0.02, row=ROWS // 2)] = False

    # ACT
    lines = _measure(scene)

    # ASSERT
    assert len(lines) == len(LINE_COLUMNS) - 1
    assert all(abs(line.centreNm - scene.centres[3]) > 0.05 for line in lines)


def test_measure_lines_leaves_out_two_lines_that_blend() -> None:
    # ARRANGE
    scene = Scene(lineColumns=[300, 306, 900, 1500])

    # ACT
    lines = _measure(scene)

    # ASSERT
    assert sorted(round(line.centreNm) for line in lines) == sorted(
        round(c) for c in (scene.centres[2], scene.centres[3])
    )


def test_measure_lines_ignores_light_outside_the_slit_range() -> None:
    # ARRANGE
    scene = Scene(minSlit=4.0)

    # ACT
    lines = _measure(scene)

    # ASSERT
    assert lines == []


def test_measure_lines_uses_the_slit_range_it_is_given() -> None:
    # ARRANGE: LINES ONLY BETWEEN 2.5 AND 4.5 ARCSEC FROM THE CENTRE
    scene = Scene(minSlit=2.5, maxSlit=4.5)

    # ACT
    lines = _measure(scene, slitRange=(2.0, 5.0))

    # ASSERT
    assert len(lines) == len(LINE_COLUMNS)


def test_measure_lines_returns_nothing_when_no_pixel_belongs_to_an_order() -> None:
    # ARRANGE
    scene = Scene()

    # ACT
    lines = _measure(scene, order=np.zeros_like(scene.order))

    # ASSERT
    assert lines == []


def test_median_resolving_power_is_the_median_over_lines() -> None:
    # ARRANGE
    lines = [line_fitting.LineFit(1, 500.0, 0.1, 3.0, power) for power in (4000.0, 5000.0, 9000.0)]

    # ACT
    median = line_fitting.median_resolving_power(lines)

    # ASSERT
    assert median == 5000.0
