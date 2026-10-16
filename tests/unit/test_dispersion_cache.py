"""Observable cache contracts for dispersion-map coefficient loading."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table
from numpy.testing import assert_allclose, assert_array_equal

from soxspipe.commonutils.dispersion_map_to_pixel_arrays import (
    _read_dispersion_map_axes,
    default_coefficients,
    dispersion_map_to_pixel_arrays,
)
from tests.factories import dispersion_map_fits

pytestmark = pytest.mark.unit


@pytest.fixture(autouse=True)
def clear_coefficient_cache() -> None:
    _read_dispersion_map_axes.cache_clear()
    yield
    _read_dispersion_map_axes.cache_clear()


def _coordinates() -> pd.DataFrame:
    return pd.DataFrame({"order": [1.0], "wavelength": [2.0], "slit_position": [0.0]})


def _convert(log: object, path: Path) -> pd.DataFrame:
    return dispersion_map_to_pixel_arrays(
        log,
        str(path),
        _coordinates(),
        removeOffDetectorLocation=False,
    )


def test_same_resolved_path_and_mtime_reuses_parsed_coefficients(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log: object,
) -> None:
    mapPath = dispersion_map_fits(tmp_path / "map.fits")
    aliasPath = tmp_path / "alias.fits"
    aliasPath.symlink_to(mapPath)
    originalRead = Table.read
    reads: list[Path] = []

    def counting_read(path: str, *args: object, **kwargs: object) -> Table:
        reads.append(Path(path))
        return originalRead(path, *args, **kwargs)

    monkeypatch.setattr(Table, "read", counting_read)

    first = _convert(log, mapPath)
    second = _convert(log, aliasPath)

    assert_allclose(first[["fit_x", "fit_y"]], second[["fit_x", "fit_y"]])
    assert reads == [mapPath.resolve()]
    assert _read_dispersion_map_axes.cache_info().hits == 1


def test_distinct_paths_changed_mtime_and_reset_each_reload_coefficients(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log: object,
) -> None:
    firstPath = dispersion_map_fits(tmp_path / "first.fits")
    secondPath = dispersion_map_fits(tmp_path / "second.fits")
    originalRead = Table.read
    reads: list[Path] = []

    def counting_read(path: str, *args: object, **kwargs: object) -> Table:
        reads.append(Path(path))
        return originalRead(path, *args, **kwargs)

    monkeypatch.setattr(Table, "read", counting_read)

    _convert(log, firstPath)
    _convert(log, secondPath)
    stat = firstPath.stat()
    os.utime(firstPath, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
    _convert(log, firstPath)
    _read_dispersion_map_axes.cache_clear()
    _convert(log, firstPath)

    assert reads == [
        firstPath.resolve(),
        secondPath.resolve(),
        firstPath.resolve(),
        firstPath.resolve(),
    ]


def test_fits_coefficients_preserve_axes_degrees_and_nan_padding(
    tmp_path: Path,
) -> None:
    coefficients = pd.DataFrame(
        {
            "axis": ["x", "y"],
            "order_deg": [0, 0],
            "wavelength_deg": [0, 0],
            "slit_deg": [0, 0],
            "c000": [3.5, 7.25],
            "unused": [np.nan, np.nan],
        }
    )
    mapPath = dispersion_map_fits(tmp_path / "constant.fits", coefficients=coefficients)

    axes = _read_dispersion_map_axes(str(mapPath.resolve()), mapPath.stat().st_mtime)

    assert axes == {
        "x": {"orderDeg": 0, "wavelengthDeg": 0, "slitDeg": 0, "coeff": (3.5,)},
        "y": {"orderDeg": 0, "wavelengthDeg": 0, "slitDeg": 0, "coeff": (7.25,)},
    }


def test_zero_degree_maps_use_one_coefficient_per_axis(tmp_path: Path, log: object) -> None:
    coefficients = pd.DataFrame(
        {
            "axis": ["x", "y"],
            "order_deg": [0, 0],
            "wavelength_deg": [0, 0],
            "slit_deg": [0, 0],
            "c000": [3.5, 7.25],
        }
    )
    mapPath = dispersion_map_fits(tmp_path / "constant.fits", coefficients=coefficients)

    result = _convert(log, mapPath)

    assert_allclose(result["fit_x"], [3.5])
    assert_allclose(result["fit_y"], [7.25])


def test_default_coefficients_are_unit_valued_with_axis_specific_degree_lengths() -> None:
    coefficients = default_coefficients(orderDeg=[1, 2], wavelengthDeg=[0, 1], slitDeg=[0, 0])

    assert_array_equal(coefficients["x"], np.ones(2))
    assert_array_equal(coefficients["y"], np.ones(6))


def test_default_coefficients_share_one_length_for_a_single_set_of_degrees() -> None:
    coefficients = default_coefficients(orderDeg=2, wavelengthDeg=3, slitDeg=1)

    assert_array_equal(coefficients["x"], np.ones(24))
    assert_array_equal(coefficients["y"], np.ones(24))
