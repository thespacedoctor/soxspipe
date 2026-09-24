"""Characterization tests for image rectification and overlap geometry."""

from __future__ import annotations

import importlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.base_util import base_util
from soxspipe.commonutils.image_transformer import (
    _clip_halfplane_nb,
    _clip_to_pixel_area_nb,
    _pixel_boundaries_grid,
    _resample_weights_kernel,
    image_transformer,
)

pytestmark = pytest.mark.unit


def _transformer(log: object, *, zoomFactor: int = 1) -> image_transformer:
    """Return a minimal transformer without invoking detector-sized setup."""
    transformer = image_transformer.__new__(image_transformer)
    transformer.log = log
    transformer.zoomFactor = zoomFactor
    transformer._cache_image_names = set()
    return transformer


@pytest.mark.parametrize(
    ("operation", "expected"),
    [
        ("sum", np.array([[10.0, 18.0], [42.0, 50.0]])),
        ("mean", np.array([[2.5, 4.5], [10.5, 12.5]])),
    ],
)
def test_unzoom_rebins_complete_blocks(
    log: object,
    operation: str,
    expected: np.ndarray,
) -> None:
    transformer = _transformer(log, zoomFactor=2)
    values = np.arange(16.0).reshape(4, 4)

    result = transformer._unzoom(values, operation=operation)

    np.testing.assert_array_equal(result, expected)


def test_unzoom_returns_small_array_and_rejects_unknown_operation(log: object) -> None:
    transformer = _transformer(log, zoomFactor=5)
    values = np.ones((2, 2))

    assert transformer._unzoom(values) is values
    transformer.zoomFactor = 1
    with pytest.raises(ValueError, match="Invalid operation"):
        transformer._unzoom(values, operation="median")


def test_pixel_boundary_grid_traces_each_cell_counter_clockwise() -> None:
    slitEdges = np.array([-1.0, 0.0, 1.0])
    wavelengthEdges = np.array([500.0, 502.0])

    slitBlock, wavelengthBlock = _pixel_boundaries_grid(
        slitEdges,
        wavelengthEdges,
        edge_samples=1,
    )

    assert slitBlock.shape == (2, 1, 4)
    np.testing.assert_array_equal(slitBlock[0, 0], [-1.0, 0.0, 0.0, -1.0])
    np.testing.assert_array_equal(
        wavelengthBlock[0, 0],
        [500.0, 500.0, 502.0, 502.0],
    )


def test_halfplane_clipping_interpolates_polygon_crossings() -> None:
    xCoordinates = np.array([-1.0, 1.0, 1.0, -1.0])
    yCoordinates = np.array([-1.0, -1.0, 1.0, 1.0])
    outputX = np.empty(8)
    outputY = np.empty(8)

    count = _clip_halfplane_nb.py_func(
        xCoordinates,
        yCoordinates,
        4,
        0,
        0.0,
        True,
        outputX,
        outputY,
    )

    assert count == 4
    assert np.all(outputX[:count] <= 0.0)
    np.testing.assert_allclose(np.sort(outputY[:count]), [-1.0, -1.0, 1.0, 1.0])


def test_polygon_clipping_returns_exact_pixel_overlap_areas() -> None:
    buffers = tuple(np.empty(12) for _ in range(4))
    squareX = np.array([-0.5, 0.5, 0.5, -0.5])
    squareY = np.array([-0.5, -0.5, 0.5, 0.5])

    fullArea = _clip_to_pixel_area_nb.py_func(
        squareX,
        squareY,
        4,
        0,
        0,
        *buffers,
    )
    outsideArea = _clip_to_pixel_area_nb.py_func(
        squareX + 3,
        squareY,
        4,
        0,
        0,
        *buffers,
    )

    assert fullArea == pytest.approx(1.0)
    assert outsideArea == 0.0


def test_resampling_kernel_records_positive_overlaps_and_skips_sentinels() -> None:
    xBoundaries = np.array([[[[-0.5, 0.5, 0.5, -0.5], [0.0, 0.0, 0.0, 0.0]]]]).reshape(
        1, 2, 4
    )
    yBoundaries = np.array([[[[-0.5, -0.5, 0.5, 0.5], [0.0, 0.0, 0.0, 0.0]]]]).reshape(
        1, 2, 4
    )
    pixelLow = np.array([[0, 0]], dtype=np.int32)
    pixelHigh = np.array([[0, -1]], dtype=np.int32)
    outputs = tuple(np.empty(2, dtype=np.int32) for _ in range(4))
    areas = np.empty(2)
    coverage = np.zeros((1, 2))

    count = _resample_weights_kernel.py_func(
        xBoundaries,
        yBoundaries,
        pixelLow,
        pixelHigh,
        pixelLow,
        pixelHigh,
        *outputs,
        areas,
        coverage,
    )

    assert count == 1
    assert areas[0] == pytest.approx(1.0)
    np.testing.assert_array_equal(coverage, [[1.0, 0.0]])


def test_cache_image_records_flux_mask_coverage_and_rectified_views(
    log: object,
) -> None:
    transformer = _transformer(log)
    transformer.uniqueOrders = [10]
    transformer.orderSlitEdges = [np.array([0.0, 1.0, 2.0])]
    transformer.orderWlEdges = [np.array([500.0, 501.0, 502.0])]
    transformer.orderSlices = [pd.DataFrame()]
    transformer._resamplingWeights = {
        10: {
            "py": np.array([0, 0, 1, 1]),
            "px": np.array([0, 1, 0, 1]),
            "area": np.ones(4),
            "flatIdx": np.arange(4),
            "coverage": np.ones((2, 2)),
        }
    }
    flux = np.array([[1.0, 2.0], [3.0, 4.0]])
    mask = np.array([[False, True], [False, False]])

    coverage = transformer.cache_image(
        "flux",
        flux,
        associatedMask=mask,
        returnCoverage=True,
    )
    rectified = transformer.get_order_rectified()

    np.testing.assert_array_equal(coverage, [np.ones((2, 2))])
    np.testing.assert_array_equal(rectified[0]["flux"], flux)
    np.testing.assert_array_equal(rectified[0]["bpMask"], mask)
    assert transformer.get_order_slices() is transformer.orderSlices


def test_cache_image_without_mask_or_coverage_preserves_optional_contract(
    log: object,
) -> None:
    transformer = _transformer(log)
    transformer.uniqueOrders = [1]
    transformer.orderSlitEdges = [np.array([0.0, 1.0])]
    transformer.orderWlEdges = [np.array([0.0, 1.0])]
    transformer.orderSlices = [pd.DataFrame()]
    transformer._resamplingWeights = {
        1: {
            "py": np.array([0]),
            "px": np.array([0]),
            "area": np.array([0.5]),
            "flatIdx": np.array([0]),
            "coverage": np.array([[0.5]]),
        }
    }

    result = transformer.cache_image("variance", np.array([[8.0]]))

    assert result is None
    assert "bpMask" not in transformer.orderSlices[0]
    np.testing.assert_array_equal(
        transformer.get_order_rectified()[0]["variance"],
        [[4.0]],
    )


def test_true_coordinate_cache_uses_bin_centres(log: object) -> None:
    transformer = _transformer(log)
    transformer.uniqueOrders = [12]
    transformer.orderSlitEdges = [np.array([-1.0, 0.0, 1.0])]
    transformer.orderWlEdges = [np.array([500.0, 502.0, 504.0])]
    transformer.orderSlices = [pd.DataFrame()]
    transformer.wlMinMax = [(500.0, 504.0)]

    transformer._cache_true_wavelength_slit_images()
    rectified = transformer.get_order_rectified()[0]

    np.testing.assert_array_equal(rectified["wavelength"], [[501.0, 503.0]] * 2)
    np.testing.assert_array_equal(rectified["slit"], [[-0.5, -0.5], [0.5, 0.5]])
    assert transformer.get_order_wavelength_ranges() == [(500.0, 504.0)]


@pytest.mark.parametrize("dispersionAxis", ["x", "y"])
def test_rectified_boundaries_follow_map_orientation_and_skip_absent_orders(
    log: object,
    dispersionAxis: str,
) -> None:
    transformer = _transformer(log)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = dispersionAxis
    transformer.pixelScale = 1.0
    transformer.slitLengthArcsec = 2.0
    transformer.orderPixelTable = pd.DataFrame(
        {"order": [10, 10], "xcoord_centre": [0.0, 1.0], "ycoord": [0, 0]}
    )
    transformer.mapDF = pd.DataFrame(
        {
            "x": [0, 1],
            "y": [0, 0],
            "slit_position": [-0.5, 0.5],
            "wavelength": [500.0, 502.0],
        }
    )
    transformer.orderNums = np.array([10, 11])
    transformer.amins = np.array([0.0, 0.0])
    transformer.amaxs = np.array([2.0, 2.0])
    transformer.waveLengthMin = np.array([500.0, 600.0])
    transformer.waveLengthMax = np.array([504.0, 604.0])
    transformer.uniqueOrders = np.array([10])

    slitEdges, wavelengthEdges = transformer._determine_rectified_image_boundaries()

    assert transformer.uniqueOrders == [10]
    assert len(slitEdges) == len(wavelengthEdges) == 1
    assert transformer.wlMinMax == [(500.0, 504.0)]
    assert len(transformer.orderSlices) == 1


def test_sigma_clip_bpm_preserves_existing_bad_pixels(log: object) -> None:
    transformer = _transformer(log)
    rawFlux = np.ones((2, 10))
    badPixels = np.zeros((2, 10), dtype=int)
    badPixels[0, 0] = 7

    result = transformer._sigma_clip_bpm(rawFlux, badPixels, order=10)

    assert result.shape == rawFlux.shape
    assert bool(result[0, 0])
    assert badPixels[0, 0] == 1


def test_constructor_prepares_geometry_weights_and_coordinate_cache(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[str] = []

    def fake_base_init(
        transformer: image_transformer,
        suppliedLog: object,
        settings: dict[str, object],
        *,
        associatedFrame: object,
        dispersionMap: str,
        twoDMapPath: str,
    ) -> None:
        transformer.log = suppliedLog
        transformer.twoDMap = {"WAVELENGTH": SimpleNamespace(data=np.zeros((3, 4)))}
        transformer.dispersionMap = dispersionMap

    def fake_boundaries(
        transformer: image_transformer,
    ) -> tuple[list[np.ndarray], list[np.ndarray]]:
        calls.append("boundaries")
        transformer.orderSlices = [pd.DataFrame()]
        return [np.array([0.0, 1.0])], [np.array([500.0, 501.0])]

    def fake_weights(transformer: image_transformer) -> dict[int, object]:
        calls.append("weights")
        return {10: {}}

    def fake_coordinate_cache(transformer: image_transformer) -> None:
        calls.append("coordinates")

    monkeypatch.setattr(base_util, "__init__", fake_base_init)
    monkeypatch.setattr(
        image_transformer,
        "_determine_rectified_image_boundaries",
        fake_boundaries,
    )
    monkeypatch.setattr(
        image_transformer,
        "_precompute_resampling_weights",
        fake_weights,
    )
    monkeypatch.setattr(
        image_transformer,
        "_cache_true_wavelength_slit_images",
        fake_coordinate_cache,
    )

    transformer = image_transformer(
        log=log,
        settings={},
        orderPixelTable=pd.DataFrame({"order": [10]}),
        twoDMapPath="map.fits",
        dispersionMap="coefficients.fits",
        associatedFrame=object(),
        slitHalfLength=5,
        edgeSamples=2,
    )

    assert calls == ["boundaries", "weights", "coordinates"]
    assert (transformer.ny, transformer.nx) == (3, 4)
    assert transformer.uniqueOrders.tolist() == [10]
    assert transformer.slitLengthArcsec == pytest.approx(2.8)
    assert transformer._resamplingWeights == {10: {}}


def test_precomputed_weights_convert_boundaries_once_and_preserve_area(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    transformer = _transformer(log)
    transformer.edgeSamples = 1
    transformer.uniqueOrders = [10]
    transformer.orderSlitEdges = [np.array([0.0, 1.0])]
    transformer.orderWlEdges = [np.array([0.0, 1.0])]
    transformer.dispersionMap = "coefficients.fits"
    transformer.nx = 3
    transformer.ny = 3
    received: list[pd.DataFrame] = []

    def fake_conversion(
        *,
        log: object,
        dispersionMapPath: str,
        orderPixelTable: pd.DataFrame,
        removeOffDetectorLocation: bool,
        trimColumns: bool,
    ) -> pd.DataFrame:
        received.append(orderPixelTable.copy(deep=True))
        return orderPixelTable.assign(
            fit_x=orderPixelTable["wavelength"],
            fit_y=orderPixelTable["slit_position"],
        )

    dispersionModule = importlib.import_module(
        "soxspipe.commonutils.dispersion_map_to_pixel_arrays"
    )
    monkeypatch.setattr(
        dispersionModule, "dispersion_map_to_pixel_arrays", fake_conversion
    )

    weights = transformer._precompute_resampling_weights()[10]

    assert len(received) == 1
    assert list(received[0].columns) == ["order", "wavelength", "slit_position"]
    assert len(weights["area"]) == 4
    assert weights["area"].sum() == pytest.approx(1.0)
    assert weights["coverage"][0, 0] == pytest.approx(1.0)
    np.testing.assert_array_equal(weights["flatIdx"], np.zeros(4, dtype=int))
