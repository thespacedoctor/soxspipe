"""Characterization tests for image rectification and overlap geometry."""

from __future__ import annotations

import importlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.base_util import base_util
from soxspipe.commonutils.image_transformer import (
    _arcsec_per_pixel,
    _clip_halfplane_nb,
    _clip_to_pixel_area_nb,
    _pixel_boundaries_grid,
    _pixel_uniform_wavelength_edges,
    _rebin_resampling_weights,
    _resample_weights_kernel,
    image_transformer,
)

pytestmark = pytest.mark.unit


def _transformer(
    log: object,
    *,
    zoomFactorSlit: int = 1,
    zoomFactorWavelength: int = 1,
) -> image_transformer:
    """Return a minimal transformer without invoking detector-sized setup."""
    transformer = image_transformer.__new__(image_transformer)
    transformer.log = log
    transformer.zoomFactorSlit = zoomFactorSlit
    transformer.zoomFactorWavelength = zoomFactorWavelength
    transformer._cache_image_names = set()
    transformer.slitHalfLength = 1
    transformer._measure_order_trace_geometry = _linear_trace_geometry(transformer)
    return transformer


def _linear_trace_geometry(transformer: image_transformer):
    """Stub the dispersion-map trace sampling with linear wavelength edges and 1 arcsec/pixel."""

    def measure(orderTraces: list[tuple]) -> list[tuple[np.ndarray, float]]:
        results = []
        for _order, wlmin, wlmax, pixelRange, _centreCoeffs in orderTraces:
            step = (wlmax - wlmin) / pixelRange / transformer.zoomFactorWavelength
            results.append((np.arange(wlmin, wlmax, step), 1.0))
        return results

    return measure


def _rebinned(
    i: list[int],
    j: list[int],
    px: list[int],
    py: list[int],
    area: list[float],
    shape: tuple[int, int],
) -> dict[str, np.ndarray]:
    """Build a rebinned-weights dict the way _precompute_resampling_weights stores it."""
    return _rebin_resampling_weights(
        np.array(i), np.array(j), np.array(px), np.array(py), np.array(area, dtype=float),
        nSp=shape[0], nWl=shape[1], zoomSlit=1, zoomWavelength=1, nx=8, ny=8,
    )


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
    transformer = _transformer(log, zoomFactorSlit=2, zoomFactorWavelength=2)
    values = np.arange(16.0).reshape(4, 4)

    result = transformer._unzoom(values, operation=operation)

    np.testing.assert_array_equal(result, expected)


def test_unzoom_returns_small_array_and_rejects_unknown_operation(log: object) -> None:
    transformer = _transformer(log, zoomFactorSlit=5, zoomFactorWavelength=5)
    values = np.ones((2, 2))

    assert transformer._unzoom(values) is values
    transformer.zoomFactorSlit = 1
    transformer.zoomFactorWavelength = 1
    with pytest.raises(ValueError, match="Invalid operation"):
        transformer._unzoom(values, operation="median")


@pytest.mark.parametrize(
    ("operation", "expected"),
    [
        ("sum", np.array([[24.0, 42.0], [96.0, 114.0]])),
        ("mean", np.array([[4.0, 7.0], [16.0, 19.0]])),
    ],
)
def test_unzoom_bins_slit_rows_and_wavelength_columns_independently(
    log: object,
    operation: str,
    expected: np.ndarray,
) -> None:
    transformer = _transformer(log, zoomFactorSlit=2, zoomFactorWavelength=3)
    values = np.arange(24.0).reshape(4, 6)

    result = transformer._unzoom(values, operation=operation)

    np.testing.assert_array_equal(result, expected)


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
            "shape": (2, 2),
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
            "shape": (1, 1),
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
    transformer.orderSlitCentreCoeffs = [np.array([0.0])]
    transformer.wlMinMax = [(500.0, 504.0)]

    transformer._cache_true_wavelength_slit_images()
    rectified = transformer.get_order_rectified()[0]

    np.testing.assert_array_equal(rectified["wavelength"], [[501.0, 503.0]] * 2)
    np.testing.assert_array_equal(rectified["slit"], [[-0.5, -0.5], [0.5, 0.5]])
    assert transformer.get_order_wavelength_ranges() == [(500.0, 504.0)]


def test_true_coordinate_cache_adds_slit_centre_per_column(log: object) -> None:
    transformer = _transformer(log)
    transformer.uniqueOrders = [12]
    transformer.orderSlitEdges = [np.array([-1.0, 0.0, 1.0])]
    transformer.orderWlEdges = [np.array([500.0, 502.0, 504.0])]
    transformer.orderSlices = [pd.DataFrame()]
    # SLIT CENTRE COEFFICIENTS FOR np.polyval: centre(wavelength) == wavelength
    transformer.orderSlitCentreCoeffs = [np.array([1.0, 0.0])]
    transformer.wlMinMax = [(500.0, 504.0)]

    transformer._cache_true_wavelength_slit_images()
    rectified = transformer.get_order_rectified()[0]

    # BIN-CENTRE SLIT OFFSETS (-0.5, 0.5) SHIFTED BY THE PER-COLUMN TRACE CENTRE (501.0, 503.0)
    np.testing.assert_array_equal(
        rectified["slit"], [[500.5, 502.5], [501.5, 503.5]]
    )


@pytest.mark.parametrize("dispersionAxis", ["x", "y"])
def test_rectified_boundaries_follow_map_orientation_and_skip_absent_orders(
    log: object,
    dispersionAxis: str,
) -> None:
    transformer = _transformer(log)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = dispersionAxis
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


def test_rectified_boundaries_use_separate_slit_and_wavelength_zoom(
    log: object,
) -> None:
    transformer = _transformer(log, zoomFactorSlit=4, zoomFactorWavelength=2)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
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
    transformer.orderNums = np.array([10])
    transformer.amins = np.array([0.0])
    transformer.amaxs = np.array([2.0])
    transformer.waveLengthMin = np.array([500.0])
    transformer.waveLengthMax = np.array([504.0])
    transformer.uniqueOrders = np.array([10])

    slitEdges, wavelengthEdges = transformer._determine_rectified_image_boundaries()

    np.testing.assert_allclose(np.diff(slitEdges[0]), 0.25)
    np.testing.assert_allclose(np.diff(wavelengthEdges[0]), 1.0)


def test_rectified_boundaries_centre_slit_on_trace_per_order(log: object) -> None:
    transformer = _transformer(log)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    transformer.orderPixelTable = pd.DataFrame(
        {
            "order": [10, 10, 11, 11],
            "xcoord_centre": [0.0, 1.0, 0.0, 1.0],
            "ycoord": [0, 0, 1, 1],
        }
    )
    # ORDER 10 SITS AT A CONSTANT SLIT POSITION OF -0.5, ORDER 11 AT A CONSTANT +1.5
    transformer.mapDF = pd.DataFrame(
        {
            "x": [0, 1, 0, 1],
            "y": [0, 0, 1, 1],
            "slit_position": [-0.5, -0.5, 1.5, 1.5],
            "wavelength": [500.0, 502.0, 600.0, 602.0],
        }
    )
    transformer.orderNums = np.array([10, 11])
    transformer.amins = np.array([0.0, 0.0])
    transformer.amaxs = np.array([2.0, 2.0])
    transformer.waveLengthMin = np.array([500.0, 600.0])
    transformer.waveLengthMax = np.array([504.0, 604.0])
    transformer.uniqueOrders = np.array([10, 11])

    transformer._determine_rectified_image_boundaries()

    assert len(transformer.orderSlitCentreCoeffs) == 2
    order10Centre = np.polyval(transformer.orderSlitCentreCoeffs[0], 501.0)
    order11Centre = np.polyval(transformer.orderSlitCentreCoeffs[1], 601.0)
    assert order10Centre == pytest.approx(-0.5)
    assert order11Centre == pytest.approx(1.5)


def test_rectified_boundaries_centre_follows_trace_along_wavelength(
    log: object,
) -> None:
    transformer = _transformer(log)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    transformer.orderPixelTable = pd.DataFrame(
        {"order": [10, 10], "xcoord_centre": [0.0, 1.0], "ycoord": [0, 0]}
    )
    # SLIT POSITION RISES LINEARLY WITH WAVELENGTH ALONG THE TRACE
    transformer.mapDF = pd.DataFrame(
        {
            "x": [0, 1],
            "y": [0, 0],
            "slit_position": [0.0, 2.0],
            "wavelength": [500.0, 502.0],
        }
    )
    transformer.orderNums = np.array([10])
    transformer.amins = np.array([0.0])
    transformer.amaxs = np.array([2.0])
    transformer.waveLengthMin = np.array([500.0])
    transformer.waveLengthMax = np.array([504.0])
    transformer.uniqueOrders = np.array([10])

    transformer._determine_rectified_image_boundaries()

    fittedCentre = np.polyval(transformer.orderSlitCentreCoeffs[0], 501.0)
    assert fittedCentre == pytest.approx(1.0)


def test_rectified_boundaries_fall_back_to_global_mean_for_order_without_valid_trace(
    log: object,
) -> None:
    transformer = _transformer(log)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    transformer.orderPixelTable = pd.DataFrame(
        {
            # ORDER 11's PIXEL COORDINATES ARE ABSENT FROM mapDF, SO ITS LOOKUP IS ALL-NaN
            "order": [10, 10, 11, 11],
            "xcoord_centre": [0.0, 1.0, 99.0, 98.0],
            "ycoord": [0, 0, 0, 0],
        }
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
    transformer.uniqueOrders = np.array([10, 11])

    transformer._determine_rectified_image_boundaries()

    # GLOBAL nanmean OF THE ONLY VALID (ORDER 10) SLIT POSITIONS IS 0.0
    fallbackCentre = np.polyval(transformer.orderSlitCentreCoeffs[1], 601.0)
    assert fallbackCentre == pytest.approx(0.0)
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("11" in message for message in warnings)


def test_rectified_boundaries_degrade_degree_on_rank_deficient_trace(
    log: object,
) -> None:
    transformer = _transformer(log)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    # TWO PIXELS SHARE THE SAME MAPPED WAVELENGTH (500.0) BUT DISAGREE ON SLIT POSITION —
    # A DEGREE >= 1 FIT ON THESE THREE POINTS IS RANK-DEFICIENT AT THE REQUESTED DEGREE
    transformer.orderPixelTable = pd.DataFrame(
        {"order": [10, 10, 10], "xcoord_centre": [0.0, 0.0, 1.0], "ycoord": [0, 1, 0]}
    )
    transformer.mapDF = pd.DataFrame(
        {
            "x": [0, 0, 1],
            "y": [0, 1, 0],
            "slit_position": [-0.5, 0.5, 0.0],
            "wavelength": [500.0, 500.0, 502.0],
        }
    )
    transformer.orderNums = np.array([10])
    transformer.amins = np.array([0.0])
    transformer.amaxs = np.array([2.0])
    transformer.waveLengthMin = np.array([500.0])
    transformer.waveLengthMax = np.array([504.0])
    transformer.uniqueOrders = np.array([10])

    transformer._determine_rectified_image_boundaries()

    # THE FIT MUST DEGRADE TO SOMETHING WELL-POSED RATHER THAN RETURNING GARBAGE/NaN COEFFICIENTS
    coeffs = transformer.orderSlitCentreCoeffs[0]
    assert np.all(np.isfinite(coeffs))


def test_rectified_boundaries_raise_when_every_order_trace_is_invalid(
    log: object,
) -> None:
    transformer = _transformer(log)
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    # NONE OF THESE PIXEL COORDINATES EXIST IN mapDF, SO EVERY LOOKUP IS NaN
    transformer.orderPixelTable = pd.DataFrame(
        {"order": [10, 10], "xcoord_centre": [99.0, 98.0], "ycoord": [0, 0]}
    )
    transformer.mapDF = pd.DataFrame(
        {
            "x": [0, 1],
            "y": [0, 0],
            "slit_position": [-0.5, 0.5],
            "wavelength": [500.0, 502.0],
        }
    )
    transformer.orderNums = np.array([10])
    transformer.amins = np.array([0.0])
    transformer.amaxs = np.array([2.0])
    transformer.waveLengthMin = np.array([500.0])
    transformer.waveLengthMax = np.array([504.0])
    transformer.uniqueOrders = np.array([10])

    with pytest.raises(ValueError, match="non-finite"):
        transformer._determine_rectified_image_boundaries()


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
    transformer.orderSlitCentreCoeffs = [np.array([0.0])]
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


def test_precomputed_weights_shift_corners_by_slit_centre(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    transformer = _transformer(log)
    transformer.edgeSamples = 1
    transformer.uniqueOrders = [10]
    transformer.orderSlitEdges = [np.array([0.0, 1.0])]
    transformer.orderWlEdges = [np.array([0.0, 1.0])]
    # CONSTANT SLIT CENTRE OFFSET OF 2.0 ARCSEC, INDEPENDENT OF WAVELENGTH
    transformer.orderSlitCentreCoeffs = [np.array([2.0])]
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

    transformer._precompute_resampling_weights()

    # RAW SLIT OFFSET CORNERS FOR A SINGLE CELL (BOTTOM, RIGHT, TOP, LEFT) ARE 0, 1, 1, 0 —
    # THE CONSTANT +2.0 CENTRE MUST BE ADDED BEFORE THE CORNERS ARE HANDED TO THE CONVERTER
    np.testing.assert_array_equal(
        received[0]["slit_position"].to_numpy(), [2.0, 3.0, 3.0, 2.0]
    )


def test_pixel_uniform_wavelength_edges_step_evenly_along_the_detector_trace() -> None:
    # NON-LINEAR DISPERSION: DETECTOR PATH LENGTH GROWS QUADRATICALLY WITH WAVELENGTH
    wavelengths = np.linspace(500.0, 510.0, 2001)
    pixelX = 100.0 * (wavelengths - 500.0) + 3.0 * (wavelengths - 500.0) ** 2
    pixelY = 0.1 * pixelX

    edges = _pixel_uniform_wavelength_edges(wavelengths, pixelX, pixelY, samplesPerPixel=2)

    pathAtEdges = np.interp(edges, wavelengths, np.hypot(pixelX, pixelY))
    assert edges[0] == pytest.approx(500.0)
    np.testing.assert_allclose(np.diff(pathAtEdges), 0.5, rtol=1e-3)
    # EDGES MUST BE NON-UNIFORM IN WAVELENGTH FOR A NON-LINEAR DISPERSION
    assert np.ptp(np.diff(edges)) > 1e-4


def test_pixel_uniform_wavelength_edges_use_longest_finite_run_without_bridging_gaps() -> None:
    wavelengths = np.array([500.0, 501.0, 502.0, 503.0, 504.0])
    pixelX = np.array([0.0, np.nan, 20.0, 30.0, 40.0])
    pixelY = np.zeros(5)

    edges = _pixel_uniform_wavelength_edges(wavelengths, pixelX, pixelY, samplesPerPixel=1)

    # THE ISOLATED SAMPLE AT 500 nm IS DISCARDED RATHER THAN JOINED ACROSS THE NaN GAP
    assert edges[0] == pytest.approx(502.0)
    assert np.all(np.isfinite(edges))
    assert len(edges) == 20


def test_pixel_uniform_wavelength_edges_raise_without_a_usable_trace() -> None:
    with pytest.raises(ValueError, match="trace"):
        _pixel_uniform_wavelength_edges(
            np.array([500.0, 501.0]), np.array([np.nan, 1.0]), np.zeros(2), samplesPerPixel=1
        )


def test_arcsec_per_pixel_divides_the_slit_probe_by_detector_distance() -> None:
    centreX = np.array([0.0, 10.0, 20.0])
    centreY = np.zeros(3)
    offsetX = centreX + np.array([0.0, 0.0, np.nan])
    offsetY = np.array([5.0, 5.0, 5.0])

    scale = _arcsec_per_pixel(1.0, centreX, centreY, offsetX, offsetY)

    assert scale == pytest.approx(0.2)


def test_arcsec_per_pixel_raises_when_no_probe_lands_on_the_detector() -> None:
    with pytest.raises(ValueError, match="arcsec"):
        _arcsec_per_pixel(1.0, np.zeros(2), np.zeros(2), np.full(2, np.nan), np.zeros(2))


def test_rebin_resampling_weights_merges_repeat_pixels_and_drops_trimmed_cells() -> None:
    # 2 x 5 ZOOMED GRID REBINNED 2 x 2 -> 1 x 2 OUTPUT; COLUMN j=4 IS TRIMMED BY _unzoom
    weights = _rebin_resampling_weights(
        i=np.array([0, 1, 0, 1, 0]),
        j=np.array([0, 1, 2, 3, 4]),
        px=np.array([0, 0, 1, 2, 3]),
        py=np.zeros(5, dtype=int),
        area=np.array([0.25, 0.25, 0.5, 0.5, 0.9]),
        nSp=2, nWl=5, zoomSlit=2, zoomWavelength=2, nx=8, ny=8,
    )

    assert weights["shape"] == (1, 2)
    order = np.lexsort((weights["px"], weights["flatIdx"]))
    np.testing.assert_array_equal(weights["flatIdx"][order], [0, 1, 1])
    np.testing.assert_array_equal(weights["px"][order], [0, 1, 2])
    np.testing.assert_allclose(weights["area"][order], [0.5, 0.5, 0.5])


def test_rebin_resampling_weights_pass_through_when_grid_is_smaller_than_zoom() -> None:
    weights = _rebin_resampling_weights(
        i=np.array([0]), j=np.array([0]), px=np.array([1]), py=np.array([2]),
        area=np.array([0.3]), nSp=1, nWl=1, zoomSlit=5, zoomWavelength=5, nx=8, ny=8,
    )

    assert weights["shape"] == (1, 1)
    np.testing.assert_array_equal(weights["py"], [2])
    np.testing.assert_allclose(weights["area"], [0.3])


def test_cache_variance_propagates_with_squared_merged_weights(log: object) -> None:
    transformer = _transformer(log)
    transformer.uniqueOrders = [3]
    transformer.orderSlitEdges = [np.array([0.0, 1.0])]
    transformer.orderWlEdges = [np.array([0.0, 1.0, 2.0])]
    transformer.orderSlices = [pd.DataFrame()]
    # CELL 0: HALF OF PIXEL (0,0) SPLIT OVER TWO SUB-CELLS (MERGED TO 0.5 BEFORE SQUARING)
    # CELL 1: HALF OF PIXEL (1,0) AND HALF OF PIXEL (2,0), INDEPENDENT
    rebinned = _rebin_resampling_weights(
        i=np.zeros(4, dtype=int), j=np.array([0, 0, 1, 1]),
        px=np.array([0, 0, 1, 2]), py=np.zeros(4, dtype=int),
        area=np.array([0.25, 0.25, 0.5, 0.5]),
        nSp=1, nWl=2, zoomSlit=1, zoomWavelength=1, nx=8, ny=8,
    )
    transformer._resamplingWeights = {3: {**rebinned, "coverage": np.ones((1, 2))}}
    variance = np.full((1, 8), 4.0)

    transformer.cache_variance("variance", variance)

    np.testing.assert_allclose(
        transformer.get_order_rectified()[0]["variance"], [[0.5**2 * 4.0, 2 * 0.5**2 * 4.0]]
    )


def test_rectified_boundaries_size_slit_bins_from_measured_arcsec_per_pixel(
    log: object,
) -> None:
    transformer = _transformer(log, zoomFactorSlit=4)
    transformer.slitHalfLength = 2
    measuredEdges = np.array([500.0, 500.5, 501.5, 503.0])
    transformer._measure_order_trace_geometry = lambda orderTraces: [(measuredEdges, 0.5)]
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    transformer.orderPixelTable = pd.DataFrame(
        {"order": [10, 10], "xcoord_centre": [0.0, 1.0], "ycoord": [0, 0]}
    )
    transformer.mapDF = pd.DataFrame(
        {"x": [0, 1], "y": [0, 0], "slit_position": [-0.5, 0.5], "wavelength": [500.0, 502.0]}
    )
    transformer.orderNums = np.array([10])
    transformer.amins = np.array([0.0])
    transformer.amaxs = np.array([2.0])
    transformer.waveLengthMin = np.array([500.0])
    transformer.waveLengthMax = np.array([504.0])
    transformer.uniqueOrders = np.array([10])

    slitEdges, wavelengthEdges = transformer._determine_rectified_image_boundaries()

    # 2 PIXELS EITHER SIDE OF THE TRACE AT 0.5 ARCSEC/PIXEL, 4 SUB-SAMPLES PER PIXEL
    assert slitEdges[0][0] == pytest.approx(-1.0)
    np.testing.assert_allclose(np.diff(slitEdges[0]), 0.125)
    np.testing.assert_array_equal(wavelengthEdges[0], measuredEdges)


def test_measure_order_trace_geometry_samples_trace_and_slit_probe_in_one_conversion(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    transformer = image_transformer.__new__(image_transformer)
    transformer.log = log
    transformer.zoomFactorWavelength = 1
    transformer.dispersionMap = "coefficients.fits"
    transformer.nx = 10000
    transformer.ny = 10000
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
        # 10 PIXELS PER nm ALONG X, 4 PIXELS PER ARCSEC ALONG Y
        return orderPixelTable.assign(
            fit_x=10.0 * orderPixelTable["wavelength"],
            fit_y=4.0 * orderPixelTable["slit_position"],
        )

    dispersionModule = importlib.import_module(
        "soxspipe.commonutils.dispersion_map_to_pixel_arrays"
    )
    monkeypatch.setattr(dispersionModule, "dispersion_map_to_pixel_arrays", fake_conversion)

    results = transformer._measure_order_trace_geometry(
        [(10, 500.0, 502.0, 20, np.array([0.0])), (11, 600.0, 601.0, 10, np.array([1.0]))]
    )

    assert len(received) == 1
    assert set(received[0]["order"]) == {10, 11}
    (edges10, scale10), (edges11, scale11) = results
    np.testing.assert_allclose(np.diff(edges10), 0.1)
    assert edges10[0] == pytest.approx(500.0)
    assert edges11[0] == pytest.approx(600.0)
    assert scale10 == pytest.approx(0.25)
    assert scale11 == pytest.approx(0.25)


def test_measure_order_trace_geometry_trims_trace_samples_that_fall_off_the_detector(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    transformer = image_transformer.__new__(image_transformer)
    transformer.log = log
    transformer.zoomFactorWavelength = 1
    transformer.dispersionMap = "coefficients.fits"
    transformer.nx = 50
    transformer.ny = 50

    def fake_conversion(**kwargs: object) -> pd.DataFrame:
        table = kwargs["orderPixelTable"]
        # 10 PIXELS PER nm FROM x = 0 AT 500 nm: EVERYTHING REDWARD OF ~505 nm IS OFF A 50-PIXEL DETECTOR
        return table.assign(
            fit_x=10.0 * (table["wavelength"] - 500.0),
            fit_y=10.0 + 4.0 * table["slit_position"],
        )

    dispersionModule = importlib.import_module(
        "soxspipe.commonutils.dispersion_map_to_pixel_arrays"
    )
    monkeypatch.setattr(dispersionModule, "dispersion_map_to_pixel_arrays", fake_conversion)

    ((edges, scale),) = transformer._measure_order_trace_geometry(
        [(10, 500.0, 510.0, 100, np.array([0.0]))]
    )

    assert edges[0] == pytest.approx(500.0)
    assert edges[-1] < 505.0
    np.testing.assert_allclose(np.diff(edges), 0.1)
    assert scale == pytest.approx(0.25)
