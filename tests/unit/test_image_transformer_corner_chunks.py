"""The cell-corner conversion runs in row chunks, so its temporary arrays do not grow with the detector."""

from __future__ import annotations

import importlib
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.dispersion_map_to_pixel_arrays import (
    _read_dispersion_map_axes,
    dispersion_map_to_pixel_arrays,
)
from soxspipe.commonutils.image_transformer import image_transformer
from tests.factories import dispersion_map_fits

pytestmark = pytest.mark.unit

ORDER_NUMBER = 13
CORNER_COUNT = 25
SMALL_CHUNK_ROWS = 7


@pytest.fixture(autouse=True)
def clear_coefficient_cache() -> Iterator[None]:
    _read_dispersion_map_axes.cache_clear()
    yield
    _read_dispersion_map_axes.cache_clear()


def _curved_map(path: Path) -> Path:
    """Write a dispersion map with cross terms in order, wavelength and slit position."""
    rng = np.random.default_rng(11)
    coefficients = {
        "axis": ["x", "y"],
        "order_deg": [2, 2],
        "wavelength_deg": [2, 2],
        "slit_deg": [1, 1],
    }
    for index in range(18):
        coefficients[f"c{index:03d}"] = rng.normal(0.0, 1.0, 2) * 10.0 ** (1 - index % 3)
    return dispersion_map_fits(path, coefficients=pd.DataFrame(coefficients))


def _corners() -> tuple[np.ndarray, np.ndarray]:
    wavelengths = np.linspace(950.0, 1010.0, CORNER_COUNT)
    slitPositions = np.cos(np.linspace(0.0, 3.0, CORNER_COUNT)) * 4.5
    return wavelengths, slitPositions


def _transformer(log: object, mapPath: Path) -> image_transformer:
    transformer = image_transformer.__new__(image_transformer)
    transformer.log = log
    transformer.dispersionMap = str(mapPath)
    return transformer


def _single_call(log: object, mapPath: Path) -> pd.DataFrame:
    wavelengths, slitPositions = _corners()
    table = pd.DataFrame(
        {
            "order": np.full(CORNER_COUNT, ORDER_NUMBER, dtype=int),
            "wavelength": wavelengths,
            "slit_position": slitPositions,
        }
    )
    return dispersion_map_to_pixel_arrays(
        log=log,
        dispersionMapPath=str(mapPath),
        orderPixelTable=table,
        removeOffDetectorLocation=False,
        trimColumns=True,
    )


def test_chunked_corner_conversion_equals_the_single_call_conversion(log: object, tmp_path: Path) -> None:
    mapPath = _curved_map(tmp_path / "map.fits")
    wavelengths, slitPositions = _corners()
    expected = _single_call(log, mapPath)
    transformer = _transformer(log, mapPath)

    fitX, fitY = transformer._corner_pixel_positions(
        ORDER_NUMBER, wavelengths, slitPositions, chunkRows=SMALL_CHUNK_ROWS
    )

    np.testing.assert_array_equal(fitX, expected["fit_x"].to_numpy())
    np.testing.assert_array_equal(fitY, expected["fit_y"].to_numpy())
    assert fitX.dtype == np.float32
    assert fitY.dtype == np.float32


def test_corner_conversion_hands_the_converter_no_more_rows_than_the_chunk_size(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    mapPath = _curved_map(tmp_path / "map.fits")
    wavelengths, slitPositions = _corners()
    transformer = _transformer(log, mapPath)
    dispersionModule = importlib.import_module("soxspipe.commonutils.dispersion_map_to_pixel_arrays")
    rowCounts: list[int] = []

    def counting_conversion(**kwargs: object) -> pd.DataFrame:
        rowCounts.append(len(kwargs["orderPixelTable"]))
        return dispersion_map_to_pixel_arrays(**kwargs)

    monkeypatch.setattr(dispersionModule, "dispersion_map_to_pixel_arrays", counting_conversion)

    transformer._corner_pixel_positions(ORDER_NUMBER, wavelengths, slitPositions, chunkRows=SMALL_CHUNK_ROWS)

    wholeChunks, remainder = divmod(CORNER_COUNT, SMALL_CHUNK_ROWS)
    assert rowCounts == [SMALL_CHUNK_ROWS] * wholeChunks + [remainder]


def test_corner_conversion_applies_the_frame_binning_to_every_chunk(log: object, tmp_path: Path) -> None:
    mapPath = _curved_map(tmp_path / "map.fits")
    wavelengths, slitPositions = _corners()
    expected = _single_call(log, mapPath)
    transformer = _transformer(log, mapPath)
    transformer.dispersionBinRatio = (2.0, 1.0)

    fitX, fitY = transformer._corner_pixel_positions(
        ORDER_NUMBER, wavelengths, slitPositions, chunkRows=SMALL_CHUNK_ROWS
    )

    binnedX, binnedY = transformer._to_frame_binning(expected["fit_x"].to_numpy(), expected["fit_y"].to_numpy())
    np.testing.assert_array_equal(fitX, binnedX)
    np.testing.assert_array_equal(fitY, binnedY)


@pytest.mark.parametrize("chunkRows", [0, -3])
def test_corner_conversion_rejects_a_chunk_size_below_one(log: object, tmp_path: Path, chunkRows: int) -> None:
    mapPath = _curved_map(tmp_path / "map.fits")
    wavelengths, slitPositions = _corners()
    transformer = _transformer(log, mapPath)

    with pytest.raises(ValueError, match="chunkRows must be at least 1"):
        transformer._corner_pixel_positions(ORDER_NUMBER, wavelengths, slitPositions, chunkRows=chunkRows)
