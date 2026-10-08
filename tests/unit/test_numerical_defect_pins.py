"""Strict-xfail pins for the numerical defects found by the DY-35 defect hunt.

Every test here pins a defect that has been filed but not fixed. Each test asserts the
CORRECT behaviour, so it fails today and is marked ``xfail(strict=True)``. When someone fixes
the defect the test XPASSes, strict mode turns that into a suite failure, and the fixer must
delete the ``xfail`` marker. If another test file pins the old (wrong) behaviour as a
characterization, the fixer must update or delete that pin in the same change.

The ``reason`` of each marker reads ``DY-<issue> (DY-35 <ID>): <defect>``, naming the Linear
issue that tracks the defect and its ID in the DY-35 findings table.

Each marker names the exception the defect produces today (``AssertionError`` for a wrong
value), so a setup, import or fixture error fails the suite instead of hiding as an xfail.
Some pins accept only one remedy (for example "refit on the kept points" for B1 and B2). If
the owner rules for a different remedy, rewrite the pin with the fix rather than deleting it.

Several builders are imported from neighbouring characterization test modules. That coupling
is accepted for now: if one of those helpers is renamed or moved, retarget the import here.
"""

from __future__ import annotations

from collections.abc import Callable
from numbers import Real
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import keyword_lookup, toolkit
from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.horne_extraction import fit_object_profile, horne_extraction
from soxspipe.commonutils.response_function import _fit_response_polynomial, _ResponseFitConvergenceError
from soxspipe.recipes import base_recipe
from soxspipe.recipes.soxs_mbias import soxs_mbias
from tests.factories import instrument_header, pipeline_settings
from tests.unit.test_create_dispersion_map_helpers import _mapper
from tests.unit.test_fit_global_polynomial import _design_matrix
from tests.unit.test_image_transformer_geometry import _transformer
from tests.unit.test_mbias_qc import _Frames
from tests.unit.test_toolkit_characterization import (
    _detector_lookup_stub,
    _identity_keyword_lookup,
    _order_frame,
    _unpack_order_table_stub,
)

pytestmark = pytest.mark.unit


class _SilentLog:
    """Accept any logger call, for code that only needs a logger-shaped object."""

    def __getattr__(self, name: str) -> Callable[..., None]:
        return lambda *args, **kwargs: None


# ---------------------------------------------------------------------------
# A1: MASTER RON
# ---------------------------------------------------------------------------


def _bias_frame(rng: np.random.Generator, shape: tuple[int, int], hotPixels: dict[tuple[int, int], float]) -> CCDData:
    data = rng.normal(500.0, 3.8, shape).astype(np.float32)
    for (row, column), value in hotPixels.items():
        # PERSISTENT HOT PIXEL, THE SAME IN EVERY FRAME
        data[row, column] = value + rng.normal(0, 3.8)
    header = instrument_header(overrides={"DPR_TYPE": "BIAS", "DPR_CATG": "CALIB", "DPR_TECH": "IMAGE"})
    return CCDData(
        data,
        unit=u.electron,
        meta=header,
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(shape, 3.8), unit=u.electron),
    )


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1249 (DY-35 A1): MASTER RON is the std of the unmasked array, so clipped and masked pixels count",
)
def test_master_ron_excludes_masked_pixels(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE: FIVE BIAS FRAMES SHARING 20 HOT PIXELS (0.05% OF THE DETECTOR)
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)
    rng = np.random.default_rng(7)
    shape = (200, 200)
    hotPixels = {(int(row), int(column)): 2000.0 for row, column in rng.integers(0, shape[0], (20, 2))}
    frames = [_bias_frame(rng, shape, hotPixels) for _ in range(5)]
    recipe = soxs_mbias.__new__(soxs_mbias)
    recipe.log = log
    recipe.settings = pipeline_settings(tmp_path)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.arm = "VIS"
    recipe.detectorParams = {}
    recipe.imageType = "BIAS"
    recipe.recipeSettings = {
        "stacked-clipping-sigma": 5.0,
        "stacked-clipping-iterations": 5,
        "frame-clipping-sigma": 5.0,
        "frame-clipping-iterations": 5,
    }
    recipe.verbose = False
    recipe.debug = False
    recipe.qc = pd.DataFrame()
    recipe.inputFrames = _Frames(frames)

    # ACT
    master, medianLevel, _rawRon, masterRon, _ = recipe._combine_bias_frames()

    # ASSERT: THE RON IS THE SPREAD OF THE PIXELS THAT SURVIVED CLIPPING
    maskedStd = float(np.ma.std(np.ma.array(master.data - medianLevel, mask=master.mask)))
    assert master.mask.any()
    assert masterRon == pytest.approx(maskedStd, rel=0.05)


# ---------------------------------------------------------------------------
# A2: clip_and_stack drops per-frame input masks
# ---------------------------------------------------------------------------

MASKED_PIXEL = (3, 3)
MASKED_FRAME_VALUES = [12.0, 9.0, 10.0, 11.0, 10.5]


def _stacker(tmp_path: Path, log: object) -> base_recipe:
    recipe = object.__new__(base_recipe)
    recipe.log = log
    recipe.settings = pipeline_settings(tmp_path)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.arm = "VIS"
    recipe.detectorParams = {}
    recipe.imageType = "BIAS"
    recipe.recipeSettings = {"stacked-clipping-sigma": 3.0, "stacked-clipping-iterations": 2}
    recipe.verbose = False
    recipe.debug = False
    return recipe


def _stack_with_one_masked_value(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, *, ignoreInputMasks: bool
) -> CCDData:
    """Stack five frames; frame 0 flags the pixel carrying 12.0 as bad but it is not a gross outlier."""
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)
    rng = np.random.default_rng(1)
    frames = []
    for frameIndex, value in enumerate(MASKED_FRAME_VALUES):
        data = 10.0 + rng.normal(0, 1.0, (8, 8))
        data[MASKED_PIXEL] = value
        mask = np.zeros((8, 8), dtype=bool)
        mask[MASKED_PIXEL] = frameIndex == 0
        frames.append(
            CCDData(
                data.astype(np.float32),
                unit=u.electron,
                meta=instrument_header(),
                mask=mask,
                uncertainty=StdDevUncertainty(np.full((8, 8), 3.0), unit=u.electron),
            )
        )
    return _stacker(tmp_path, log).clip_and_stack(
        frames, recipe="soxs_mbias", ignore_input_masks=ignoreInputMasks, post_stack_clipping=False
    )


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1254 (DY-35 A2): clip_and_stack discards pixels flagged in only some of the input frames",
)
def test_clip_and_stack_excludes_pixel_masked_in_one_input_frame(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT: FIVE FRAMES, THE 12.0 VALUE MASKED IN ITS OWN FRAME ONLY
    combined = _stack_with_one_masked_value(tmp_path, log, monkeypatch, ignoreInputMasks=False)

    # ASSERT: MEAN OF THE FOUR UNMASKED FRAMES (9.0, 10.0, 11.0, 10.5)
    assert combined.data[MASKED_PIXEL] == pytest.approx(10.125)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1254 (DY-35 A2): ignore_input_masks has no effect on the stacked value",
)
def test_clip_and_stack_ignore_input_masks_changes_the_stacked_value(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    honouringMasks = _stack_with_one_masked_value(tmp_path, log, monkeypatch, ignoreInputMasks=False)
    ignoringMasks = _stack_with_one_masked_value(tmp_path, log, monkeypatch, ignoreInputMasks=True)

    # ASSERT
    assert honouringMasks.data[MASKED_PIXEL] != ignoringMasks.data[MASKED_PIXEL]


# ---------------------------------------------------------------------------
# C1 (= A3): merged spectrum fabricates a zero-variance end bin
# ---------------------------------------------------------------------------


def _merge_one_order(log: object, firstWavelength: float) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.arm = "VIS"
    extractor.kw = lambda key: key
    wavelength = firstWavelength + 0.02 * np.arange(10)
    count = len(wavelength)
    extraction = pd.DataFrame(
        {
            "order": [10] * count,
            "wavelengthMean": wavelength,
            "pixelScaleNm": [0.02] * count,
            "extractedFluxOptimal": [100.0] * count,
            "extractedFluxBoxcar": [100.0] * count,
            "extractedFluxBoxcarRobust": [100.0] * count,
            "varianceSpectrum": [25.0] * count,
            "skyFlux": [1.0] * count,
        }
    )
    merged, _ = extractor.merge_extracted_orders(extraction)
    return wavelength, merged, extraction


def test_merged_spectrum_has_no_bin_outside_the_extracted_wavelength_range(log: object) -> None:
    # ARRANGE AND ACT: THE FIRST SAMPLE (500.004 NM) ROUNDS DOWN TO 500.00 ON THE 0.02 NM OUTPUT GRID
    wavelength, merged, _ = _merge_one_order(log, firstWavelength=500.004)

    # ASSERT
    mergedWavelength = merged["WAVE"].values.value
    mergedVariance = merged["VARIANCE"].values.value
    assert mergedWavelength.min() >= wavelength.min()
    assert mergedWavelength.max() <= wavelength.max()
    assert (mergedVariance > 0).all()


# ---------------------------------------------------------------------------
# A4 / B1: response-function fit
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1255 (DY-35 A4): response fit returns all-NaN coefficients from a non-finite response, no error",
)
def test_response_fit_returns_finite_coefficients_or_raises_for_a_non_finite_response() -> None:
    # ARRANGE: A ZERO-FLUX SAMPLE GIVES AN INFINITE RAW RESPONSE; THE SMOOTHING SPREADS IT OVER ABOUT 41 SAMPLES
    wavelength = np.linspace(1000.0, 2000.0, 400)
    rawResponse = np.full(400, 2.0)
    rawResponse[0] = np.inf

    # ACT: THE ROUTINE'S OWN FIT ERROR IS THE OTHER ACCEPTABLE OUTCOME; ANY OTHER ERROR PROPAGATES
    try:
        coefficients, _, _ = _fit_response_polynomial(
            wavelength, rawResponse, polynomialOrder=3, maxIterations=5, excludedRegions=[], smoothingSigma=5
        )
    except _ResponseFitConvergenceError:
        return

    # ASSERT
    assert np.isfinite(coefficients).all()


def test_response_fit_coefficients_match_a_refit_on_the_returned_points() -> None:
    # ARRANGE: A BLOCK OF HIGH POINTS PULLS THE FIRST FIT, SO THE SINGLE ALLOWED PASS STILL DELETES POINTS
    wavelength = np.linspace(500.0, 900.0, 1000)
    rawResponse = np.full(1000, 1.0)
    rawResponse[295:306] = 3.0

    # ACT: FAILING LOUDLY ON A CAP EXIT IS THE OTHER ACCEPTABLE OUTCOME
    try:
        coefficients, keptWavelength, keptResponse = _fit_response_polynomial(
            wavelength, rawResponse, polynomialOrder=1, maxIterations=1
        )
    except _ResponseFitConvergenceError:
        return

    # ASSERT: THE FIRST FIT OF A FRESH CALL ON THE KEPT POINTS IS THE FIT THOSE POINTS GIVE
    assert len(keptWavelength) < len(wavelength)
    refit, _, _ = _fit_response_polynomial(keptWavelength, keptResponse, polynomialOrder=1, maxIterations=1)
    np.testing.assert_allclose(coefficients, refit, atol=1e-9)


# ---------------------------------------------------------------------------
# A5: rectification leaks masked-pixel values
# ---------------------------------------------------------------------------


def test_rectified_unmasked_cell_does_not_carry_the_masked_pixel_value(log: object) -> None:
    # ARRANGE: ONE RECTIFIED CELL, 85% FROM PIXEL (0,0) AND 15% FROM THE FLAGGED HOT PIXEL (0,1)
    transformer = _transformer(log)
    transformer.uniqueOrders = [10]
    transformer.orderSlitEdges = [np.array([0.0, 1.0])]
    transformer.orderWlEdges = [np.array([0.0, 1.0])]
    transformer.orderSlices = [pd.DataFrame(index=[0])]
    transformer._resamplingWeights = {
        10: {
            "flatIdx": np.array([0, 0]),
            "px": np.array([0, 1]),
            "py": np.array([0, 0]),
            "area": np.array([0.85, 0.15]),
            "shape": (1, 1),
            "coverage": np.ones((1, 1)),
        }
    }
    flux = np.array([[100.0, 60000.0]])
    mask = np.array([[False, True]])

    # ACT
    transformer.cache_image("fluxRaw", flux, associatedMask=mask)
    rectified = transformer.get_order_rectified()[0]

    # ASSERT: THE CELL IS EITHER MASKED OR FREE OF THE HOT PIXEL (EXPECT 85 TO 100, NOT 9085)
    assert rectified["bpMask"][0, 0] or rectified["fluxRaw"][0, 0] < 200


# ---------------------------------------------------------------------------
# B2: global polynomial fit returns stale coefficients on a cap exit
# ---------------------------------------------------------------------------


def test_global_polynomial_coefficients_match_a_refit_on_the_kept_rows(log: object) -> None:
    # ARRANGE: FIVE ORDERS WITH A BLOCK OF CORRUPTED MEASUREMENTS AND THE NIR MFLAT TWO-PASS CAP
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.arm = "NIR"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.orderDeg = 2
    detector.axisBDeg = 2
    detector.recipeName = "soxs-mflat"
    detector.recipeSettings = {"poly-fitting-residual-clipping-sigma": 3, "poly-clipping-iteration-limit": 2}
    rng = np.random.default_rng(1)
    frames = []
    for order in range(11, 16):
        yValues = np.linspace(10, 2000, 60)
        xValues = 100 + 50 * order + 0.01 * yValues + 1e-6 * yValues**2 + rng.normal(0, 0.3, yValues.size)
        frames.append(pd.DataFrame({"order": float(order), "cont_y": yValues, "cont_x": xValues}))
    pixels = pd.concat(frames, ignore_index=True)
    corrupted = rng.choice(len(pixels), 25, replace=False)
    pixels.loc[corrupted, "cont_x"] += rng.normal(0, 8, corrupted.size)

    # ACT
    coefficients, kept, clipped = detector.fit_global_polynomial(pixels, axisACol="cont_x", axisBCol="cont_y")

    # ASSERT: THE RETURNED FIT IS THE LEAST-SQUARES FIT TO THE ROWS IT KEEPS
    assert len(clipped) > 0
    design = _design_matrix(kept, 2, 2)
    refit, *_ = np.linalg.lstsq(design, kept["cont_x"].to_numpy(), rcond=None)
    np.testing.assert_allclose(design @ np.asarray(coefficients), design @ refit, atol=1e-6)


# ---------------------------------------------------------------------------
# B3: dispersion solution raises TypeError when lines < coefficients
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason="DY-1256 (DY-35 B3): fit_polynomials raises TypeError instead of returning xerror when lines < coefficients",
)
def test_fit_polynomials_returns_xerror_when_fewer_lines_than_coefficients(log: object) -> None:
    # ARRANGE: TEN LINES CANNOT CONSTRAIN THE 16 COEFFICIENTS OF A 3 x 3 x 0 POLYNOMIAL
    mapper = _mapper(log)
    mapper.detectorParams = {}
    mapper.recipeSettings = {"poly-fitting-residual-clipping-sigma": 5, "poly-clipping-iteration-limit": 5}
    mapper.settings = {"workspace-root-dir": "/nonexistent"}
    rng = np.random.default_rng(0)
    lineCount = 10
    pixelTable = pd.DataFrame(
        {
            "dropped": False,
            "order": rng.integers(1, 10, lineCount).astype(float),
            "wavelength": rng.uniform(500, 900, lineCount),
            "slit_position": 0.0,
            "observed_x": rng.uniform(0, 2000, lineCount),
            "observed_y": rng.uniform(0, 4000, lineCount),
        }
    )

    # ACT
    result = mapper.fit_polynomials(orderPixelTable=pixelTable, wavelengthDeg=3, orderDeg=3, slitDeg=0)

    # ASSERT: THE CALLER REDUCES THE DEGREES WHEN IT SEES THE FALLBACK MARKER
    assert result[0] == "xerror"


# ---------------------------------------------------------------------------
# B4: Horne profile clipping never detects convergence
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1257 (DY-35 B4): fit_object_profile keeps iterating after a pass that clips nothing",
)
def test_object_profile_clipping_stops_once_a_pass_clips_nothing(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE: PURE GAUSSIAN NOISE, SO CLIPPING SETTLES AFTER A FEW PASSES WELL BEFORE THE CAP OF TEN
    slitRows, wavelengthColumns = 7, 2000
    rng = np.random.default_rng(3)
    profile = np.exp(-0.5 * ((np.arange(slitRows) - 3) / 1.2) ** 2)
    flux = profile[:, None] * 1000.0 + rng.normal(0, 10, (slitRows, wavelengthColumns))
    images = {"fluxRaw": flux, "mask": np.zeros_like(flux, dtype=bool)}
    # THE SPY COUNTS np.polyfit CALLS; A FIX THAT SWAPS THE FITTING ROUTINE MUST RETARGET IT
    pointCounts: list[int] = []
    realPolyfit = np.polyfit

    def recording_polyfit(x: np.ndarray, y: np.ndarray, deg: int, *args: object, **kwargs: object) -> np.ndarray:
        pointCounts.append(len(x))
        return realPolyfit(x, y, deg, *args, **kwargs)

    monkeypatch.setattr(np, "polyfit", recording_polyfit)

    # ACT
    fit_object_profile(
        pd.DataFrame({"x": range(wavelengthColumns)}),
        images,
        slitHalfLength=3,
        clippingSigma=3.0,
        clippingIterationLimit=10,
        hornePolyOrder=7,
        axisB="y",
        order=1,
        debug=False,
        plt=None,
    )

    # ASSERT: A ROW STARTS WITH ALL POINTS; THE POINTS ENTERING EACH LATER PASS MUST KEEP SHRINKING,
    # BECAUSE A PASS THAT CLIPS NOTHING ENDS THE LOOP
    rows: list[list[int]] = []
    for count in pointCounts:
        if count == wavelengthColumns:
            rows.append([])
        rows[-1].append(count)
    assert len(rows) == slitRows
    for rowCounts in rows:
        assert np.all(np.diff(rowCounts) < 0), rowCounts


# ---------------------------------------------------------------------------
# C2: rectified slit grid is one row short and off-centre
# ---------------------------------------------------------------------------


def test_rectified_slit_grid_has_two_h_rows_centred_on_the_trace(log: object) -> None:
    # ARRANGE: SLIT HALF LENGTH 3 AT 1 ARCSEC PER PIXEL, ONE ORDER ON A CONSTANT SLIT POSITION
    slitHalfLength = 3
    transformer = _transformer(log)
    transformer.slitHalfLength = slitHalfLength
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    transformer.orderPixelTable = pd.DataFrame({"order": [10, 10], "xcoord_centre": [0.0, 1.0], "ycoord": [0, 0]})
    transformer.mapDF = pd.DataFrame(
        {"x": [0, 1], "y": [0, 0], "slit_position": [0.0, 0.0], "wavelength": [500.0, 502.0]}
    )
    transformer.orderNums = np.array([10])
    transformer.amins = np.array([0.0])
    transformer.amaxs = np.array([2.0])
    transformer.waveLengthMin = np.array([500.0])
    transformer.waveLengthMax = np.array([504.0])
    transformer.uniqueOrders = np.array([10])

    # ACT
    slitEdges, _ = transformer._determine_rectified_image_boundaries()
    cellCentres = (slitEdges[0][:-1] + slitEdges[0][1:]) / 2.0
    rows = transformer._unzoom(np.broadcast_to(cellCentres[:, None], (len(cellCentres), 4)).copy(), operation="mean")

    # ASSERT
    assert rows.shape[0] == 2 * slitHalfLength
    assert rows[:, 0].mean() == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# C4: cut_image_slice reports the centre of the collapsed rows
# ---------------------------------------------------------------------------


def test_cut_image_slice_reports_the_centre_of_the_collapsed_rows(log: object) -> None:
    # ARRANGE: PIXEL VALUE EQUALS ROW INDEX, SO THE MEDIAN OF THE COLLAPSED SLICE IS THEIR CENTRE
    rowIndices, _ = np.mgrid[0:60, 0:80]
    frame = np.ma.masked_array(rowIndices.astype(float), mask=np.zeros((60, 80), dtype=bool))

    # ACT
    sliceOut, _, widthCentre = toolkit.cut_image_slice(
        log=log, frame=frame, width=3, length=20, x=40.0, y=30.0, sliceAxis="x", median=True
    )

    # ASSERT
    assert widthCentre == pytest.approx(float(np.ma.median(sliceOut)))


# ---------------------------------------------------------------------------
# C5: bias-structure pixel axis is not 0..N-1
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1258 (DY-35 C5): qc_bias_structure fits against a pixel axis ending at N, not N-1",
)
def test_bias_structure_slopes_use_a_zero_to_n_minus_one_pixel_axis() -> None:
    # ARRANGE: A SMALL FRAME MAKES THE AXIS ERROR LARGE; THE RAMP SLOPES ARE KNOWN EXACTLY
    captured: dict[str, float] = {}
    recipe = soxs_mbias.__new__(soxs_mbias)
    recipe.log = _SilentLog()
    recipe.add_qc = lambda **kwargs: captured.__setitem__(kwargs["qcName"], kwargs["qcValue"])
    rowCount, columnCount = 6, 8
    rowIndices, columnIndices = np.indices((rowCount, columnCount))
    bias = 100.0 + 0.01 * columnIndices + 0.02 * rowIndices

    # ACT
    structX, structY = recipe.qc_bias_structure(bias)

    # ASSERT: SUMS OVER THE OTHER AXIS MULTIPLY EACH SLOPE BY THE LENGTH OF THAT AXIS
    assert structX == pytest.approx(0.01 * rowCount)
    assert structY == pytest.approx(0.02 * columnCount)


# ---------------------------------------------------------------------------
# C6: QC values stored as strings
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1259 (DY-35 C6): inner-order pixel QC values are stored as formatted strings, not numbers",
)
def test_inner_order_qc_values_are_numeric(monkeypatch: pytest.MonkeyPatch, log: object) -> None:
    # ARRANGE: THE SAME STUBS AS THE CHARACTERIZATION TEST THAT PINS THE STRING FORM
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub("x"))
    pixels = pd.DataFrame({"xcoord_edgeup": [3, 4], "xcoord_edgelow": [1, 0], "ycoord": [1, 2]})
    monkeypatch.setattr(toolkit, "unpack_order_table", _unpack_order_table_stub(pixels))
    frame = _order_frame(badPixel=(2, 0))

    # ACT
    result = toolkit.spectroscopic_image_quality_checks(
        log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
    )

    # ASSERT
    assert len(result) > 0
    assert all(isinstance(value, Real) for value in result["qc_value"])


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1259 (DY-35 C6): dispersion-solution residual QC values are stored as formatted strings, not numbers",
)
def test_dispersion_residual_qc_values_are_numeric(log: object) -> None:
    # ARRANGE: THE TWO-LINE TABLE USED BY THE EXISTING calculate_residuals TESTS
    mapper = _mapper(log)
    mapper.arcFrame = False
    table = pd.DataFrame(
        {
            "order": [1.0, 1.0],
            "wavelength": [10.0, 20.0],
            "slit_position": [0.0, 0.0],
            "order_pow_x_0": [1.0, 1.0],
            "wavelength_pow_x_0": [1.0, 1.0],
            "wavelength_pow_x_1": [10.0, 20.0],
            "slit_position_pow_x_0": [1.0, 1.0],
            "order_pow_y_0": [1.0, 1.0],
            "wavelength_pow_y_0": [1.0, 1.0],
            "wavelength_pow_y_1": [10.0, 20.0],
            "slit_position_pow_y_0": [1.0, 1.0],
            "observed_x": [17.0, 40.0],
            "observed_y": [1.0, 5.0],
            "detector_x": [20.0, 40.0],
            "detector_y": [5.0, 5.0],
            "fwhm_pin_px": [2.0, 2.0],
        }
    )

    # ACT
    mapper.calculate_residuals(
        table,
        xcoeff=[0.0, 2.0],
        ycoeff=[5.0, 0.0],
        orderDeg=0,
        wavelengthDeg=1,
        slitDeg=0,
        pixelRange=True,
        writeQCs=True,
    )

    # ASSERT
    assert len(mapper.qc) > 0
    assert all(isinstance(value, Real) for value in mapper.qc["qc_value"])


# ---------------------------------------------------------------------------
# C7: binned map + same-binned frame raises
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="DY-1260 (DY-35 C7): a 2x2-binned map with a 2x2-binned frame is block-reduced twice and raises ValueError",
)
def test_binned_map_with_same_binned_frame_gives_one_row_per_pixel(tmp_path: Path) -> None:
    # ARRANGE: A MAP ALREADY AT THE FRAME'S 2x2 BINNING, WAVELENGTH INCREASING WITH X, ONE ORDER
    rowCount, columnCount = 4, 6
    rowIndices, columnIndices = np.mgrid[0:rowCount, 0:columnCount]
    wavelength = 500.0 + columnIndices.astype(float) + 0.1 * rowIndices
    header = fits.Header()
    header["ESO DET BINX"] = 2
    header["ESO DET BINY"] = 2
    mapPath = tmp_path / "map_2x2.fits"
    fits.HDUList(
        [
            fits.PrimaryHDU(wavelength, header=header),
            fits.ImageHDU(wavelength, name="WAVELENGTH"),
            fits.ImageHDU(np.zeros_like(wavelength), name="SLIT"),
            fits.ImageHDU(np.full_like(wavelength, 10.0), name="ORDER"),
        ]
    ).writeto(mapPath)
    keywords = {"SEQ_ARM": "ESO SEQ ARM", "WIN_BINX": "ESO DET BINX", "WIN_BINY": "ESO DET BINY"}
    frame = CCDData(
        np.ones((rowCount, columnCount)),
        unit=u.electron,
        meta=fits.Header({"ESO SEQ ARM": "VIS", "ESO DET BINX": 2, "ESO DET BINY": 2}),
        mask=np.zeros((rowCount, columnCount), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((rowCount, columnCount))),
    )

    # ACT
    mapDF, _ = toolkit.twoD_disp_map_image_to_dataframe(
        log=_SilentLog(),
        slit_length=11,
        twoDMapPath=str(mapPath),
        kw=lambda key: keywords[key],
        associatedFrame=frame,
    )

    # ASSERT
    assert len(mapDF) == rowCount * columnCount
