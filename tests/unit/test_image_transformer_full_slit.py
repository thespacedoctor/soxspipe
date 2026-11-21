"""The full-slit pass of the image transformer: fixed arcsec slit grid, squared-weight variance, skyline shifts."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.base_util import base_util
from soxspipe.commonutils.image_transformer import _column_pixel_scale_nm, image_transformer
from soxspipe.commonutils.nir_sky_refit import order_geometry, order_sky_correction

pytestmark = pytest.mark.unit

DETECTOR_ROWS = 36
DETECTOR_COLUMNS = 40
ORDER_COLUMNS = 30
ORDER_START_NM = {10: 500.0, 11: 600.0}
NM_PER_PIXEL = 0.1
ARCSEC_PER_PIXEL = {10: 0.5, 11: 0.6}
TRACE_ROW = 18.0
ZOOM = 5
SLIT_LENGTH_ARCSEC = 12.0


def _frame_arrays() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a deterministic flux image, variance image and bad-pixel mask on the synthetic detector."""
    flux = np.arange(DETECTOR_ROWS * DETECTOR_COLUMNS, dtype=float).reshape(DETECTOR_ROWS, DETECTOR_COLUMNS) * 0.5 + 3.0
    variance = 1.0 + (np.arange(DETECTOR_COLUMNS)[None, :] % 3) + 0.25 * np.arange(DETECTOR_ROWS)[:, None]
    mask = np.zeros((DETECTOR_ROWS, DETECTOR_COLUMNS), dtype=bool)
    mask[16, 5] = True
    return flux, variance, mask


@pytest.fixture
def synthetic_transformer(log: object, monkeypatch: pytest.MonkeyPatch):
    """Return a builder of real ``image_transformer`` instances on a linear two-order synthetic dispersion map.

    The map puts the wavelength axis along x at 0.1 nm per pixel and the slit axis along y, 12 rows being the
    trace, at 0.5 arcsec per pixel in order 10 and 0.6 in order 11. Cell edges fall on whole pixel columns.
    """

    def fake_base_init(transformer, suppliedLog, settings, *, associatedFrame, dispersionMap, twoDMapPath):
        transformer.log = suppliedLog
        transformer.settings = settings
        transformer.dispersionMap = dispersionMap
        transformer.twoDMap = {"WAVELENGTH": SimpleNamespace(data=np.zeros((DETECTOR_ROWS, DETECTOR_COLUMNS)))}
        transformer.detectorParams = {"slit_length": SLIT_LENGTH_ARCSEC}
        transformer.orderNums = np.array([10, 11])
        transformer.amins = np.array([0.0, 0.0])
        transformer.amaxs = np.array([float(ORDER_COLUMNS)] * 2)
        transformer.waveLengthMin = np.array([ORDER_START_NM[10], ORDER_START_NM[11]])
        transformer.waveLengthMax = np.array([ORDER_START_NM[10] + 3.0, ORDER_START_NM[11] + 3.0])

    def fake_geometry(transformer, orderTraces):
        edges = np.arange(ORDER_COLUMNS * ZOOM + 1) * NM_PER_PIXEL / ZOOM - NM_PER_PIXEL / 2
        return [(edges + ORDER_START_NM[trace[0]], ARCSEC_PER_PIXEL[trace[0]]) for trace in orderTraces]

    def fake_corners(transformer, order, wavelengths, slitPositions, chunkRows=None):
        x = (np.asarray(wavelengths) - ORDER_START_NM[order]) / NM_PER_PIXEL
        y = np.asarray(slitPositions) / ARCSEC_PER_PIXEL[order] + TRACE_ROW
        return x.astype(np.float32), y.astype(np.float32)

    monkeypatch.setattr(base_util, "__init__", fake_base_init)
    monkeypatch.setattr(image_transformer, "_read_dispersion_bin_ratio", lambda transformer: (1.0, 1.0))
    monkeypatch.setattr(image_transformer, "_measure_order_trace_geometry", fake_geometry)
    monkeypatch.setattr(image_transformer, "_corner_pixel_positions", fake_corners)

    def build(**options: object) -> image_transformer:
        arguments = {
            "log": log,
            "settings": {},
            "orderPixelTable": pd.DataFrame({"order": [10, 11]}),
            "twoDMapPath": "map.fits",
            "dispersionMap": "coefficients.fits",
            "associatedFrame": object(),
            "slitHalfLength": 3,
            "fixedSlitCentreArcsec": 0.0,
        }
        arguments.update(options)
        return image_transformer(**arguments)

    return build


def test_horne_configuration_rectifies_a_synthetic_scene_as_before(synthetic_transformer) -> None:
    # ARRANGE: THE TODAY'S HORNE CALL, 6 ROWS ABOUT THE TRACE. EVERY CELL EDGE CUTS THROUGH PIXEL CENTRES ALONG THE
    # SLIT, SO EACH CELL IS HALF OF TWO NEIGHBOURING DETECTOR PIXELS
    flux, variance, mask = _frame_arrays()
    transformer = synthetic_transformer()

    # ACT
    transformer.cache_image("fluxRaw", flux, associatedMask=mask)
    transformer.cache_image("variance", variance, associatedMask=mask)
    rectified = transformer.get_order_rectified()

    # ASSERT: THE GRID IS 6 PIXEL-SIZED ROWS OF 0.5 ARCSEC IN ORDER 10 AND 0.6 IN ORDER 11
    assert rectified[0]["fluxRaw"].shape == (6, ORDER_COLUMNS)
    assert transformer.orderSlitEdges[0][0] == -1.5
    assert transformer.orderSlitEdges[1][-1] == pytest.approx(1.8)
    # EACH CELL IS HALF OF TWO PIXELS: ROW r OF ORDER 10 HOLDS DETECTOR ROWS 15 + r AND 16 + r
    expectedFlux = 0.5 * (flux[15:21, :ORDER_COLUMNS] + flux[16:22, :ORDER_COLUMNS])
    expectedVariance = 0.5 * (variance[15:21, :ORDER_COLUMNS] + variance[16:22, :ORDER_COLUMNS])
    # THE FLAGGED PIXEL (16, 5) LEAVES THE CELLS IT FALLS IN (ROWS 0 AND 1, COLUMN 5): THE OTHER PIXEL TAKES THE CELL
    expectedFlux[0, 5] = flux[15, 5]
    expectedFlux[1, 5] = flux[17, 5]
    expectedVariance[0, 5] = variance[15, 5]
    expectedVariance[1, 5] = variance[17, 5]
    np.testing.assert_allclose(rectified[0]["fluxRaw"], expectedFlux, rtol=1e-12)
    np.testing.assert_allclose(rectified[0]["variance"], expectedVariance, rtol=1e-12)
    assert rectified[0]["bpMask"][:2, 5].all()
    assert rectified[0]["bpMask"].sum() == 2
    np.testing.assert_allclose(rectified[0]["wavelength"][0], 500.0 + NM_PER_PIXEL * np.arange(ORDER_COLUMNS))


def test_full_slit_grid_spans_the_whole_slit_in_fixed_arcsec_steps_for_every_order(synthetic_transformer) -> None:
    # ARRANGE / ACT: A 12 ARCSEC SLIT IN STEPS OF 0.5 ARCSEC, ALTHOUGH THE ORDERS SCALE AT 0.5 AND 0.6 ARCSEC PER PIXEL
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)

    # ASSERT: 24 ROWS, EACH 1/ZOOM OF A ROW WIDE IN THE ZOOMED GRID, IDENTICAL IN EVERY ORDER
    for slitEdges in transformer.orderSlitEdges:
        assert len(slitEdges) - 1 == 24 * ZOOM
        np.testing.assert_allclose(slitEdges, np.linspace(-6.0, 6.0, 24 * ZOOM + 1), rtol=0, atol=1e-12)
    np.testing.assert_array_equal(transformer.orderSlitEdges[0], transformer.orderSlitEdges[1])


def test_full_slit_step_defaults_to_the_median_arcsec_per_pixel_of_the_orders(synthetic_transformer) -> None:
    # ARRANGE / ACT: THE ORDERS SCALE AT 0.5 AND 0.6 ARCSEC PER PIXEL, SO THE MEDIAN IS 0.55 AND 12 / 0.55 NEEDS 22 ROWS
    transformer = synthetic_transformer(fullSlit=True)

    # ASSERT: THE GRID IS CENTRED AND COVERS THE SLIT, WITH THE STEP AND ROW COUNT RECORDED
    assert transformer.slitGridStepArcsec == pytest.approx(0.55)
    assert transformer.slitGridRows == 22
    for slitEdges in transformer.orderSlitEdges:
        assert slitEdges[0] == pytest.approx(-6.05)
        assert slitEdges[-1] == pytest.approx(6.05)


def test_full_slit_grid_is_the_same_for_frames_reduced_with_different_measured_plate_scales(
    synthetic_transformer, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: TWO DISPERSION MAPS THAT MEASURE DIFFERENT PLATE SCALES, ONE GIVEN STEP
    first = synthetic_transformer(fullSlit=True, slitStepArcsec=0.3)
    monkeypatch.setitem(ARCSEC_PER_PIXEL, 10, 0.52)
    monkeypatch.setitem(ARCSEC_PER_PIXEL, 11, 0.63)

    # ACT
    second = synthetic_transformer(fullSlit=True, slitStepArcsec=0.3)

    # ASSERT: 12 / 0.3 = 40 ROWS BOTH TIMES, TO THE LAST BIT
    assert first.slitGridRows == second.slitGridRows == 40
    for slitEdgesFirst, slitEdgesSecond in zip(first.orderSlitEdges, second.orderSlitEdges, strict=True):
        np.testing.assert_array_equal(slitEdgesFirst, slitEdgesSecond)


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"fullSlit": True, "fixedSlitCentreArcsec": None}, "fixedSlitCentreArcsec"),
        ({"slitHalfLength": None}, "slitHalfLength"),
        ({"slitStepArcsec": 0.5}, "only used with fullSlit"),
        ({"fullSlit": True, "slitStepArcsec": 0.0}, "positive"),
    ],
)
def test_slit_grid_options_that_cannot_give_a_grid_raise_before_any_setup(
    synthetic_transformer, options: dict[str, object], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        synthetic_transformer(**options)


def test_full_slit_grid_keeps_the_flux_of_every_covered_detector_pixel(synthetic_transformer) -> None:
    # ARRANGE: 0.5 ARCSEC ROWS ARE ONE PIXEL IN ORDER 10 AND 5/6 OF A PIXEL IN ORDER 11, 24 ROWS FROM -6 TO +6 ARCSEC
    flux, _, _ = _frame_arrays()
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)

    # ACT
    transformer.cache_image("fluxRaw", flux)
    rectified = transformer.get_order_rectified()

    # ASSERT: THE ROWS OF A COLUMN SUM TO THE DETECTOR PIXELS IT COVERS. ORDER 10 COVERS ROWS 6 TO 30 (THE ENDS HALF)
    # AND ORDER 11 ROWS 8 TO 28 (THE ENDS HALF), AS 12 ARCSEC IS 20 PIXELS AT 0.6 ARCSEC PER PIXEL
    columns = slice(0, ORDER_COLUMNS)
    expectedOrder10 = 0.5 * flux[6, columns] + flux[7:30, columns].sum(axis=0) + 0.5 * flux[30, columns]
    expectedOrder11 = 0.5 * flux[8, columns] + flux[9:28, columns].sum(axis=0) + 0.5 * flux[28, columns]
    np.testing.assert_allclose(rectified[0]["fluxRaw"].sum(axis=0), expectedOrder10, rtol=1e-12)
    np.testing.assert_allclose(rectified[1]["fluxRaw"].sum(axis=0), expectedOrder11, rtol=1e-12)
    assert rectified[0]["fluxRaw"].shape == (24, ORDER_COLUMNS)
    assert rectified[1]["fluxRaw"].shape == (24, ORDER_COLUMNS)


def _one_cell_transformer(log: object, areas: tuple[float, float]) -> image_transformer:
    """Return a transformer with one rectified cell built from detector pixels (0, 0) and (0, 1)."""
    transformer = image_transformer.__new__(image_transformer)
    transformer.log = log
    transformer._cache_image_names = set()
    transformer.uniqueOrders = [10]
    transformer.orderSlitEdges = [np.array([0.0, 1.0])]
    transformer.orderWlEdges = [np.array([0.0, 1.0])]
    transformer.orderSlices = [pd.DataFrame()]
    transformer._resamplingWeights = {
        10: {
            "flatIdx": np.array([0, 0]),
            "px": np.array([0, 1]),
            "py": np.array([0, 0]),
            "area": np.array(areas),
            "shape": (1, 1),
            "coverage": np.ones((1, 1)),
        }
    }
    return transformer


def test_squared_weights_give_the_variance_of_the_weighted_sum_exactly(log: object) -> None:
    # ARRANGE: A CELL OF 0.75 OF A PIXEL WITH VARIANCE 4 AND 0.25 OF A PIXEL WITH VARIANCE 16
    transformer = _one_cell_transformer(log, (0.75, 0.25))

    # ACT
    transformer.cache_image("variance", np.array([[4.0, 16.0]]), squaredWeights=True)

    # ASSERT: 0.75^2 * 4 + 0.25^2 * 16 = 2.25 + 1.0
    assert transformer.get_order_rectified()[0]["variance"][0, 0] == 3.25


def test_linear_weights_remain_the_default_for_variance(log: object) -> None:
    # ARRANGE
    transformer = _one_cell_transformer(log, (0.75, 0.25))

    # ACT
    transformer.cache_image("variance", np.array([[4.0, 16.0]]))

    # ASSERT: 0.75 * 4 + 0.25 * 16, THE LINEAR SUM THE HORNE LAYERS USE
    assert transformer.get_order_rectified()[0]["variance"][0, 0] == 7.0


def test_squared_weights_scale_a_masked_cell_by_the_square_of_its_renormalisation(log: object) -> None:
    # ARRANGE: THE FLAGGED PIXEL HOLDS 0.5 OF THE CELL, SO THE GOOD PIXEL'S WEIGHT GROWS FROM 0.5 TO 1.0
    transformer = _one_cell_transformer(log, (0.5, 0.5))

    # ACT
    transformer.cache_image(
        "variance", np.array([[9.0, 1.0e6]]), associatedMask=np.array([[False, True]]), squaredWeights=True
    )
    rectified = transformer.get_order_rectified()[0]

    # ASSERT: THE FLUX OF THE CELL IS THE GOOD PIXEL TIMES 1.0, SO ITS VARIANCE IS 1.0^2 * 9, WITH NO TRACE OF 1.0e6
    assert rectified["variance"][0, 0] == 9.0
    assert rectified["bpMask"][0, 0]


def test_squared_weights_give_a_lower_variance_than_linear_weights_for_independent_pixels(log: object) -> None:
    # ARRANGE: UNIT VARIANCE EVERYWHERE, SO Σ w² σ² IS Σ w²
    transformer = _one_cell_transformer(log, (0.6, 0.4))

    # ACT
    transformer.cache_image("variance", np.array([[1.0, 1.0]]), squaredWeights=True)

    # ASSERT: 0.36 + 0.16, BELOW THE 1.0 OF THE LINEAR WEIGHTS: INDEPENDENT PIXELS AVERAGE DOWN THEIR NOISE
    assert transformer.get_order_rectified()[0]["variance"][0, 0] == pytest.approx(0.52)


def test_full_slit_variance_is_the_squared_weight_sum_of_the_pixels_each_cell_covers(synthetic_transformer) -> None:
    # ARRANGE: EVERY CELL OF ORDER 10 IS HALF OF TWO DETECTOR PIXELS, ROW r HOLDING DETECTOR ROWS 6 + r AND 7 + r
    _, variance, _ = _frame_arrays()
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)

    # ACT
    transformer.cache_image("variance", variance, squaredWeights=True)
    transformer.cache_image("varianceLinear", variance)
    rectified = transformer.get_order_rectified()[0]

    # ASSERT: Σ w² σ² = 0.5² σa² + 0.5² σb², HALF THE Σ w σ² OF THE LINEAR WEIGHTS
    columns = slice(0, ORDER_COLUMNS)
    expected = 0.25 * (variance[6:30, columns] + variance[7:31, columns])
    np.testing.assert_allclose(rectified["variance"], expected, rtol=1e-12)
    np.testing.assert_allclose(rectified["varianceLinear"], 2 * expected, rtol=1e-12)


def _wavelength_layers(transformer: image_transformer) -> list[np.ndarray]:
    return [layers["wavelength"] for layers in transformer.get_order_rectified()]


def test_a_skyline_shift_moves_the_wavelength_layer_of_its_order_only(synthetic_transformer) -> None:
    # ARRANGE
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)
    flux, _, _ = _frame_arrays()
    transformer.cache_image("fluxRaw", flux)
    before = _wavelength_layers(transformer)
    fluxBefore = transformer.get_order_rectified()[0]["fluxRaw"].copy()

    # ACT: THE SHIFT IS ADDED TO THE WAVELENGTHS, THE SAME SIGN AS THE EXTRACTION'S SKYLINE TUNING
    transformer.apply_wavelength_shifts(shiftsNm={10: 0.125})
    after = _wavelength_layers(transformer)

    # ASSERT: ORDER 10 MOVES BY 0.125 NM IN EVERY ROW, ORDER 11 AND THE FLUX STAY WHERE THEY WERE
    np.testing.assert_allclose(after[0], before[0] + 0.125, rtol=0, atol=1e-12)
    np.testing.assert_array_equal(after[1], before[1])
    np.testing.assert_array_equal(transformer.get_order_rectified()[0]["fluxRaw"], fluxBefore)
    assert after[0].shape == (24, ORDER_COLUMNS)


def test_a_shift_for_all_orders_moves_every_wavelength_layer(synthetic_transformer) -> None:
    # ARRANGE
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)
    before = _wavelength_layers(transformer)

    # ACT: THE EXTRACTION KEYS A SHIFT OF EVERY ORDER BY "all", WHICH ADDS TO AN ORDER'S OWN SHIFT
    transformer.apply_wavelength_shifts(shiftsNm={"all": -0.25, 11: 0.5})
    after = _wavelength_layers(transformer)

    # ASSERT
    np.testing.assert_allclose(after[0], before[0] - 0.25, rtol=0, atol=1e-12)
    np.testing.assert_allclose(after[1], before[1] + 0.25, rtol=0, atol=1e-12)


def test_applying_wavelength_shifts_again_replaces_the_earlier_shifts(synthetic_transformer) -> None:
    # ARRANGE
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)
    before = _wavelength_layers(transformer)

    # ACT
    transformer.apply_wavelength_shifts(shiftsNm={10: 0.5})
    transformer.apply_wavelength_shifts(shiftsNm={11: 0.5})
    after = _wavelength_layers(transformer)

    # ASSERT: ORDER 10 IS BACK ON ITS ANALYTIC WAVELENGTHS, NOT 0.5 NM AWAY
    np.testing.assert_array_equal(after[0], before[0])
    np.testing.assert_allclose(after[1], before[1] + 0.5, rtol=0, atol=1e-12)


def test_a_row_by_row_sky_correction_uses_the_wavelengths_and_pixel_scale_of_the_layer(synthetic_transformer) -> None:
    # ARRANGE: AN ORDER 10 CORRECTION OF 0.5 PX AT ITS CENTRE, PLUS 1 PX PER UNIT t ALONG IT. THE ORDER'S MEDIAN PIXEL
    # SCALE IS 0.05 NM, HALF THE 0.1 NM OF THE LAYER'S COLUMNS, SO THE LOCAL DISPERSION DOUBLES (0.1 TO 0.2 NM/PX)
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)
    before = _wavelength_layers(transformer)
    geometry = order_geometry(centreNm=501.45, halfRangeNm=1.5, dispersionNmPerPx=0.1, medianPixelScaleNm=0.05)
    correction = order_sky_correction(
        order=10, shiftPx=0.5, alongPx=(1.0,), perOrderSlopePx=0.0, nLines=20, geometry=geometry
    )

    # ACT: THE SHIFTS FOR ORDER 11 ARE LEFT OUT, AS THE EXTRACTION LEAVES THEM OUT WHEN IT HAS CORRECTIONS
    transformer.apply_wavelength_shifts(shiftsNm={11: 9.0}, corrections={10: correction})
    after = _wavelength_layers(transformer)

    # ASSERT: -dispersion * (0.5 + t). THE FIRST AND LAST COLUMN HAVE NO CENTRAL DIFFERENCE, SO THEY TAKE THE MEDIAN
    t = (before[0][0] - 501.45) / 1.5
    expected = -0.2 * (0.5 + t)
    expected[[0, -1]] = -0.1 * (0.5 + t[[0, -1]])
    np.testing.assert_allclose(after[0][0] - before[0][0], expected, rtol=0, atol=1e-12)
    np.testing.assert_allclose(after[0] - before[0], np.broadcast_to(expected, after[0].shape), rtol=0, atol=1e-12)
    np.testing.assert_array_equal(after[1], before[1])


def test_empty_sky_corrections_fall_back_to_the_order_shifts(synthetic_transformer) -> None:
    # ARRANGE
    transformer = synthetic_transformer(fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None)
    before = _wavelength_layers(transformer)

    # ACT
    transformer.apply_wavelength_shifts(shiftsNm={10: 0.0625}, corrections={})

    # ASSERT
    np.testing.assert_allclose(_wavelength_layers(transformer)[0], before[0] + 0.0625, rtol=0, atol=1e-12)


@pytest.mark.parametrize("columns", [1, 2])
def test_column_pixel_scale_is_unknown_for_an_order_with_fewer_than_three_columns(columns: int) -> None:
    assert np.isnan(_column_pixel_scale_nm(np.linspace(500.0, 500.1, columns))).all()


def test_column_pixel_scale_is_half_the_difference_of_the_neighbouring_columns() -> None:
    scale = _column_pixel_scale_nm(np.array([500.0, 500.25, 500.75, 500.5]))

    np.testing.assert_array_equal(scale[1:-1], [0.375, 0.125])
    assert np.isnan(scale[[0, -1]]).all()


def test_full_slit_grid_is_centred_on_the_fixed_slit_centre(synthetic_transformer) -> None:
    # ARRANGE: A CENTRE 1 ARCSEC UP THE SLIT, 2 PIXELS IN ORDER 10, SO ITS 24 ROWS COVER DETECTOR ROWS 8 TO 32
    flux, _, _ = _frame_arrays()
    transformer = synthetic_transformer(
        fullSlit=True, slitStepArcsec=0.5, slitHalfLength=None, fixedSlitCentreArcsec=1.0
    )

    # ACT
    transformer.cache_image("fluxRaw", flux)
    rectified = transformer.get_order_rectified()[0]["fluxRaw"]

    # ASSERT
    columns = slice(0, ORDER_COLUMNS)
    expected = 0.5 * flux[8, columns] + flux[9:32, columns].sum(axis=0) + 0.5 * flux[32, columns]
    np.testing.assert_allclose(rectified.sum(axis=0), expected, rtol=1e-12)
