"""Characterization of rolling-window object clipping in subtract_sky (DY-40)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

PIXEL_COUNT = 2000
# PINNING EVERY ROW AS A LITERAL IS IMPRACTICAL; THESE ROWS SAMPLE THE FIRST
# OBJECT FLANK, A CLIPPED ROW, THE THIRD OBJECT AND THE LAST ROW
SAMPLED_ROWS = [5, 400, 1000, 1999]
FLOOR_WARNING = (
    "ORDER 10: OBJECT CLIPPING IS AT OR BELOW THE MINIMUM SIGMA LIMIT OF 1.5 "
    "WITH UNDER 5% OF PIXELS CLIPPED - NOT RETRYING"
)


def _subtractor(log: Any, *, arm: str) -> tuple[subtract_sky, list[str]]:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = arm
    subtractor.debug = True
    iterations: list[str] = []
    # DEBUG MODE PLOTS EVERY ITERATION; RECORD THE ITERATION LINE OF EACH TITLE INSTEAD
    subtractor.plot_order_skymodel_fitting_quicklook = lambda frame, spline, title=None: iterations.append(
        title.split("\n")[1]
    )
    return subtractor, iterations


def _object_order(level: float, *, seed: int = 5) -> pd.DataFrame:
    """Seeded shot noise at `level` plus four gaussian objects 4 to 16 pixels wide."""
    rng = np.random.default_rng(seed)
    index = np.arange(PIXEL_COUNT)
    flux = level + rng.normal(0.0, 1.0, PIXEL_COUNT) * np.sqrt(level)
    for position, width in enumerate((4, 8, 12, 16), start=1):
        flux += 30.0 * np.sqrt(level) * np.exp(-0.5 * ((index - position * PIXEL_COUNT / 5) / width) ** 2)
    return pd.DataFrame(
        {
            "order": 10,
            "flux": flux,
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
            "flagged_edge_clipped": False,
            "flagged_bad_pixel_clipped": False,
        }
    )


def _with_edge_and_bad_pixels(pixels: pd.DataFrame, *, edge: list[int], bad: list[int]) -> pd.DataFrame:
    """Flag `edge` and `bad` rows the way `get_over_sampled_sky_from_order` does before clipping."""
    flagged = pixels.copy()
    flagged.loc[edge, "flagged_edge_clipped"] = True
    flagged.loc[bad, "flagged_bad_pixel_clipped"] = True
    flagged.loc[edge + bad, "flagged_all_clipped"] = True
    return flagged


def _with_random_edge_pixels(pixels: pd.DataFrame, *, fraction: float, seed: int) -> pd.DataFrame:
    """Flag a seeded `fraction` of the rows as edge-clipped, so object clipping stays under 5% of the order."""
    flagged = pixels.copy()
    isEdge = np.random.default_rng(seed).random(len(flagged)) < fraction
    flagged.loc[isEdge, "flagged_edge_clipped"] = True
    flagged.loc[isEdge, "flagged_all_clipped"] = True
    return flagged


def _final_sigma(clipped: pd.DataFrame) -> float:
    """Return the sigma limit of the last clipping pass, recovered from the stored upper limit and std."""
    ratio = clipped["flux_minus_smoothed_residual_upper_limit"] / clipped["flux_minus_smoothed_residual_std"]
    return float(ratio.dropna().iloc[0])


def _record_pass_limits(subtractor: subtract_sky) -> list[tuple[pd.Series, pd.Series]]:
    """Record the stored (upper limit, std) columns of every clipping pass, which expose the exact sigma used."""
    passes: list[tuple[pd.Series, pd.Series]] = []
    previousSpy = subtractor.plot_order_skymodel_fitting_quicklook

    def spy(frame: pd.DataFrame, spline: Any, title: str | None = None) -> None:
        isFinite = frame["flux_minus_smoothed_residual_std"].notna()
        passes.append(
            (
                frame.loc[isFinite, "flux_minus_smoothed_residual_upper_limit"].copy(),
                frame.loc[isFinite, "flux_minus_smoothed_residual_std"].copy(),
            )
        )
        previousSpy(frame, spline, title=title)

    subtractor.plot_order_skymodel_fitting_quicklook = spy
    return passes


def _pass_used_exactly(limits: tuple[pd.Series, pd.Series], sigma: float) -> bool:
    """True when every row's upper limit is exactly its std times `sigma`."""
    upper, std = limits
    return bool((upper == std * sigma).all())


def _pass_used_less_than(limits: tuple[pd.Series, pd.Series], sigma: float) -> bool:
    """True when any row's upper limit is below its std times `sigma`, so the pass used a smaller sigma."""
    upper, std = limits
    return bool((upper < std * sigma).any())


def _component_union(pixels: pd.DataFrame) -> pd.Series:
    return pixels["flagged_edge_clipped"] | pixels["flagged_bad_pixel_clipped"] | pixels["flagged_object_clipped"]


def _clip(subtractor: subtract_sky, pixels: pd.DataFrame, *, sigma: float = 3.0, iterations: int = 10) -> Any:
    return subtractor.rolling_window_clipping(pixels, windowSize=11, sigma_clip_limit=sigma, max_iterations=iterations)


def _messages(log: Any, *levels: str) -> list[tuple[str, str]]:
    return [(level, message) for level, message in log.messages if level in levels]


def _approx_or_nan(values: list[float]) -> list[Any]:
    return [pytest.approx(value, rel=1e-12, abs=0, nan_ok=True) for value in values]


def test_a_faint_vis_order_retries_clipping_and_releases_pixels_from_the_abandoned_attempts(log: Any) -> None:
    """Under 5% clipped at iteration 5 restarts clipping before it converges.

    Each restart clears the object flags and releases their pixels, so only the
    pixels flagged by the final attempt are excluded. Releasing them lets the
    later attempts see the full order, so two restarts are enough, not four.
    """
    subtractor, iterations = _subtractor(log, arm="VIS")

    clipped = _clip(subtractor, _object_order(100.0))

    assert iterations == [
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.2% clipped",
        "iteration 3 - 2.7% clipped",
        "iteration 4 - 3.0% clipped",
        "iteration 5 - 3.2% clipped",
        "iteration 1 - 2.6% clipped",
        "iteration 2 - 3.1% clipped",
        "iteration 3 - 3.7% clipped",
        "iteration 4 - 4.2% clipped",
        "iteration 5 - 4.8% clipped",
        "iteration 1 - 4.1% clipped",
        "iteration 2 - 5.0% clipped",
        "iteration 3 - 5.7% clipped",
        "iteration 4 - 6.2% clipped",
        "iteration 5 - 6.6% clipped",
        "iteration 6 - 7.0% clipped",
        "iteration 7 - 7.2% clipped",
        "iteration 8 - 7.3% clipped",
    ]
    assert int(clipped["flagged_object_clipped"].sum()) == 147
    assert int(clipped["flagged_all_clipped"].sum()) == 147
    assert (clipped["flagged_all_clipped"] == _component_union(clipped)).all()
    assert _messages(log, "print", "warning") == [("print", "\tORDER 10: 147 pixels clipped in total = 7.3%)")]
    globalSigma = clipped["residual_global_sigma"]
    assert float(globalSigma.sum()) == pytest.approx(910.9525366238936, rel=1e-12, abs=0)
    assert int(globalSigma.isna().sum()) == 156
    assert globalSigma.iloc[SAMPLED_ROWS].tolist() == _approx_or_nan(
        [0.5251774431202234, np.nan, 0.31877349578230957, np.nan]
    )
    assert float(clipped["flux_percentile_smoothed"].sum()) == pytest.approx(189721.99422240542, rel=1e-12, abs=0)
    assert int(clipped["flux_percentile_smoothed"].isna().sum()) == 91


def test_a_bright_vis_order_is_skipped_at_iteration_five_and_the_warning_names_the_order(log: Any) -> None:
    """A smoothed sky above 2500 with under 5% clipped returns None, so the caller can skip that order alone."""
    subtractor, iterations = _subtractor(log, arm="VIS")

    result = _clip(subtractor, _object_order(3000.0))

    assert result is None
    assert iterations == [
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.2% clipped",
        "iteration 3 - 2.7% clipped",
        "iteration 4 - 3.0% clipped",
        "iteration 5 - 3.2% clipped",
    ]
    assert _messages(log, "print", "warning") == [
        ("warning", "ORDER 10: OBJECT IS LIKELY VERY BRIGHT - SKIPPING SKY-SUBTRACTION FOR THIS ORDER")
    ]


def test_the_same_bright_order_in_uvb_iterates_to_convergence(log: Any) -> None:
    """The retry and bright-object stop are VIS-only; UVB keeps clipping until nothing changes."""
    subtractor, iterations = _subtractor(log, arm="UVB")

    clipped = _clip(subtractor, _object_order(3000.0))

    assert iterations == [
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.2% clipped",
        "iteration 3 - 2.7% clipped",
        "iteration 4 - 3.0% clipped",
        "iteration 5 - 3.2% clipped",
        "iteration 6 - 3.4% clipped",
        "iteration 7 - 3.5% clipped",
        "iteration 8 - 3.7% clipped",
    ]
    assert int(clipped["flagged_object_clipped"].sum()) == 74
    assert int(clipped["flagged_all_clipped"].sum()) == 74
    assert _messages(log, "print", "warning") == [("print", "\tORDER 10: 74 pixels clipped in total = 3.7%)")]
    globalSigma = clipped["residual_global_sigma"]
    assert float(globalSigma.sum()) == pytest.approx(820.1670639502325, rel=1e-12, abs=0)
    assert int(globalSigma.isna().sum()) == 83
    assert globalSigma.iloc[SAMPLED_ROWS].tolist() == _approx_or_nan(
        [0.5072931128660357, np.nan, 0.22461213396868954, np.nan]
    )
    assert float(clipped["flux_percentile_smoothed"].sum()) == pytest.approx(5924705.139971411, rel=1e-12, abs=0)


def test_more_than_85_percent_clipped_logs_zero_percent_and_releases_every_object_pixel(log: Any) -> None:
    """A negative sigma clips 99.1% of pixels; the reset releases all of them from the exclusions."""
    subtractor, iterations = _subtractor(log, arm="UVB")

    clipped = _clip(subtractor, _object_order(100.0), sigma=-3.0, iterations=3)

    assert iterations == ["iteration 1 - 99.1% clipped"]
    assert int(clipped["flagged_object_clipped"].sum()) == 0
    assert int(clipped["flagged_all_clipped"].sum()) == 0
    assert _messages(log, "print", "warning") == [
        ("print", "\tORDER 10: 1982 pixels clipped in total = 99.1%)"),
        ("warning", "ORDER 10: More than 85% of pixels flagged to be clipped (99.1%). Clipping 0% instead."),
    ]
    globalSigma = clipped["residual_global_sigma"]
    assert float(globalSigma.sum()) == pytest.approx(7.526984597702706, rel=1e-12, abs=0)
    assert int(globalSigma.isna().sum()) == 1991
    assert float(clipped["flux_percentile_smoothed"].sum()) == pytest.approx(860.4901487716307, rel=1e-12, abs=0)


def test_more_than_85_percent_clipped_keeps_only_the_edge_and_bad_pixels_excluded(log: Any) -> None:
    """After the reset the all-clipped count equals the edge-or-bad-pixel count, with one row in both."""
    subtractor, _ = _subtractor(log, arm="UVB")
    pixels = _with_edge_and_bad_pixels(_object_order(100.0), edge=[0, 1, 2, 3, 4], bad=[4, 700, 1500])

    clipped = _clip(subtractor, pixels, sigma=-3.0, iterations=3)

    assert int(clipped["flagged_object_clipped"].sum()) == 0
    assert int(clipped["flagged_all_clipped"].sum()) == 7
    assert (
        clipped["flagged_all_clipped"] == (clipped["flagged_edge_clipped"] | clipped["flagged_bad_pixel_clipped"])
    ).all()


@pytest.mark.parametrize(
    "missingColumn",
    [
        pytest.param("flagged_edge_clipped", id="no-edge-column"),
        pytest.param("flagged_bad_pixel_clipped", id="no-bad-pixel-column"),
    ],
)
def test_a_reset_on_a_frame_lacking_a_flag_column_raises_a_key_error(log: Any, missingColumn: str) -> None:
    """Upstream always creates both flags, so a reset without one is a contract violation (DY-595)."""
    subtractor, _ = _subtractor(log, arm="UVB")
    pixels = _object_order(100.0).drop(columns=[missingColumn])

    with pytest.raises(KeyError, match=missingColumn):
        _clip(subtractor, pixels, sigma=-3.0, iterations=3)

    # THE 99.1% PRINT PRECEDES THE RESET BRANCH, SO THE ERROR COMES FROM THE RESET AND NOT FROM EARLIER CODE
    assert ("print", "\tORDER 10: 1982 pixels clipped in total = 99.1%)") in _messages(log, "print")
    assert _messages(log, "warning") == []


def test_a_vis_retry_keeps_the_edge_and_bad_pixels_excluded_and_leaves_no_unexplained_exclusions(log: Any) -> None:
    """After a retry every excluded pixel carries a component flag and the edge and bad pixels are still excluded."""
    subtractor, iterations = _subtractor(log, arm="VIS")
    pixels = _with_edge_and_bad_pixels(_object_order(100.0), edge=[0, 1, 2, 3, 4], bad=[4, 700, 1500])

    clipped = _clip(subtractor, pixels)

    assert sum(entry.startswith("iteration 1 ") for entry in iterations) > 1
    assert (clipped["flagged_all_clipped"] == _component_union(clipped)).all()
    assert clipped.loc[[0, 1, 2, 3, 4, 700, 1500], "flagged_all_clipped"].all()
    assert int(clipped["flagged_all_clipped"].sum()) == int(clipped["flagged_object_clipped"].sum()) + 7


def test_normal_clipping_without_a_reset_excludes_the_union_of_edge_bad_pixel_and_object_flags(log: Any) -> None:
    """Under 85% clipped on a non-retrying arm the exclusions are unchanged: edge, bad pixels and objects."""
    subtractor, iterations = _subtractor(log, arm="UVB")
    pixels = _with_edge_and_bad_pixels(_object_order(100.0), edge=[0, 1, 2, 3, 4], bad=[4, 700, 1500])

    clipped = _clip(subtractor, pixels)

    assert iterations[-1] == "iteration 8 - 3.7% clipped"
    assert int(clipped["flagged_object_clipped"].sum()) == 74
    assert int(clipped["flagged_all_clipped"].sum()) == 81
    assert (clipped["flagged_all_clipped"] == _component_union(clipped)).all()


def test_a_vis_retry_whose_first_pass_matches_the_last_pass_before_it_carries_on_clipping(log: Any) -> None:
    """The retry resets the last-clipped count, so an equal first pass does not end clipping (DY-1285)."""
    subtractor, iterations = _subtractor(log, arm="VIS")

    clipped = _clip(subtractor, _object_order(100.0, seed=121), iterations=10)

    assert iterations[:6] == [
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.1% clipped",
        "iteration 3 - 2.5% clipped",
        "iteration 4 - 2.8% clipped",
        "iteration 5 - 2.9% clipped",
        "iteration 1 - 2.9% clipped",
    ]
    assert len(iterations) == 18
    assert int(clipped["flagged_object_clipped"].sum()) == 145
    assert _final_sigma(clipped) == pytest.approx(2.8)


def test_vis_retries_stop_lowering_sigma_at_the_floor_and_warn(log: Any) -> None:
    """Every pass uses a sigma of at least 1.5 exactly, and reaching the floor logs a warning naming the order."""
    subtractor, iterations = _subtractor(log, arm="VIS")
    passes = _record_pass_limits(subtractor)
    pixels = _with_random_edge_pixels(_object_order(100.0, seed=41), fraction=0.95, seed=1041)

    clipped = _clip(subtractor, pixels, sigma=2.0, iterations=10)

    assert not any(_pass_used_less_than(limits, 1.5) for limits in passes)
    assert _pass_used_exactly(passes[-1], 1.5)
    assert len(iterations) == 31
    assert int(clipped["flagged_object_clipped"].sum()) == 69
    assert _messages(log, "warning") == [("warning", FLOOR_WARNING)]


def test_a_start_one_retry_above_the_floor_retries_twice_and_stops_exactly_at_the_floor(log: Any) -> None:
    """Starting at 1.7 the sigma goes 1.7, 1.6, 1.5 and then stops, so the third attempt is the last (DY-1285)."""
    subtractor, iterations = _subtractor(log, arm="VIS")
    passes = _record_pass_limits(subtractor)
    pixels = _with_random_edge_pixels(_object_order(100.0, seed=2), fraction=0.95, seed=1002)

    _clip(subtractor, pixels, sigma=1.7, iterations=10)

    firstPasses = [index for index, entry in enumerate(iterations) if entry.startswith("iteration 1 ")]
    assert len(firstPasses) == 3
    assert all(
        _pass_used_exactly(passes[index], sigma) for index, sigma in zip(firstPasses, (1.7, 1.6, 1.5), strict=True)
    )
    assert not any(_pass_used_less_than(limits, 1.5) for limits in passes)
    assert _messages(log, "warning") == [("warning", FLOOR_WARNING)]


def test_a_start_below_the_floor_is_not_clamped_and_does_not_retry(log: Any) -> None:
    """An initial sigma under the floor is used as given, with no retry and the same warning (DY-1285)."""
    subtractor, iterations = _subtractor(log, arm="VIS")
    passes = _record_pass_limits(subtractor)
    pixels = _with_random_edge_pixels(_object_order(100.0, seed=2), fraction=0.95, seed=1002)

    _clip(subtractor, pixels, sigma=1.2, iterations=10)

    assert sum(entry.startswith("iteration 1 ") for entry in iterations) == 1
    assert all(_pass_used_exactly(limits, 1.2) for limits in passes)
    assert _messages(log, "warning") == [("warning", FLOOR_WARNING)]


@pytest.mark.parametrize("sigma", [float("inf"), float("-inf"), float("nan")])
def test_a_non_finite_sigma_limit_raises_a_value_error_naming_the_setting(log: Any, sigma: float) -> None:
    """A sigma of inf or NaN would clip nothing or everything silently, so it is refused up front."""
    subtractor, iterations = _subtractor(log, arm="VIS")

    with pytest.raises(ValueError, match="percentile_clipping_sigma"):
        _clip(subtractor, _object_order(100.0), sigma=sigma)

    assert iterations == []


@pytest.mark.parametrize("arm", ["VI", "IS", "V"])
def test_an_arm_name_inside_the_word_vis_does_not_take_the_vis_only_branches(log: Any, arm: str) -> None:
    """Only the arm VIS retries and skips bright orders; a piece of the word VIS iterates like UVB (DY-1285)."""
    subtractor, iterations = _subtractor(log, arm=arm)
    uvbSubtractor, uvbIterations = _subtractor(log, arm="UVB")

    result = _clip(subtractor, _object_order(3000.0))
    uvbResult = _clip(uvbSubtractor, _object_order(3000.0))

    assert result is not None
    assert iterations == uvbIterations
    assert int(result["flagged_object_clipped"].sum()) == int(uvbResult["flagged_object_clipped"].sum())


BRIGHT_ORDER = 10
FAINT_ORDER = 11
FITTED_SKY = 100.0
PIXEL_ERROR = 2.0


def _two_order_subtractor(
    log: Any, monkeypatch: pytest.MonkeyPatch, *, bright: tuple[int, ...], aggressive: bool = False
) -> subtract_sky:
    """A subtractor over one faint and one bright VIS order, with the real object clipping and a flat stub sky fit.

    The order ids that reach the object-slit-range step and the sky fit are recorded in `clippedOrderIds` and
    `fittedOrderIds`.
    """
    subtractor, _ = _subtractor(log, arm="VIS")
    subtractor.debug = False
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.detectorParams = {"dispersion-axis": "x"}
    subtractor.dateObs = "2024-01-02T03:04:05"
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "bspline_order": 3,
            "clip-slit-edge-fraction": 0.0,
            "aggressive_object_masking": aggressive,
            "sky_model_qc_plot": False,
            "percentile_rolling_window_size": 11,
            "noise_rolling_window_size": 30,
        }
    }
    subtractor.objectFrame = CCDData(
        np.zeros((2, PIXEL_COUNT)),
        unit=u.electron,
        mask=np.zeros((2, PIXEL_COUNT), dtype=bool),
        uncertainty=StdDevUncertainty(np.full((2, PIXEL_COUNT), PIXEL_ERROR), unit=u.electron),
    )
    orders = []
    for row, order in enumerate((BRIGHT_ORDER, FAINT_ORDER)):
        pixels = _object_order(3000.0 if order in bright else 100.0)
        pixels["order"] = order
        pixels["slit_position"] = np.linspace(-1.0, 1.0, PIXEL_COUNT)
        pixels["x"] = np.arange(PIXEL_COUNT)
        pixels["y"] = row
        pixels["error"] = PIXEL_ERROR
        pixels["mask"] = False
        orders.append(pixels)
    subtractor.mapDF = pd.concat(orders, ignore_index=True)
    subtractor.qc = pd.DataFrame()
    subtractor.products = pd.DataFrame()

    def clip_order(imageMapOrder: pd.DataFrame, **_: Any) -> Any:
        return subtractor.rolling_window_clipping(imageMapOrder, windowSize=11, sigma_clip_limit=3, max_iterations=10)

    subtractor.clippedOrderIds = []
    subtractor.fittedOrderIds = []
    realClipObjects = subtractor.clip_object_slit_positions

    def clip_objects(orders: list[pd.DataFrame], **kwargs: Any) -> Any:
        subtractor.clippedOrderIds.append([int(order["order"].iloc[0]) for order in orders])
        return realClipObjects(orders, **kwargs)

    def fit_flat_sky(imageMapOrder: pd.DataFrame) -> Any:
        subtractor.fittedOrderIds.append(int(imageMapOrder["order"].iloc[0]))
        fitted = imageMapOrder.copy()
        fitted["sky_model"] = FITTED_SKY
        fitted["sky_subtracted_flux"] = fitted["flux"] - FITTED_SKY
        return fitted, object(), [1, 2], np.array([1.0]), 1.0

    monkeypatch.setattr(subtractor, "get_over_sampled_sky_from_order", clip_order)
    monkeypatch.setattr(subtractor, "clip_object_slit_positions", clip_objects)
    monkeypatch.setattr(subtractor, "fit_bspline_curve_to_sky", fit_flat_sky)
    return subtractor


def test_one_bright_order_keeps_the_sky_model_of_the_other_orders(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """The faint order is modelled and subtracted even though the bright order skips sky subtraction (DY-1285)."""
    subtractor = _two_order_subtractor(log, monkeypatch, bright=(BRIGHT_ORDER,))

    model, subtracted, residuals, _, _ = subtractor.subtract()

    assert (model.data[1] == FITTED_SKY).all()
    assert np.allclose(
        subtracted.data[1], subtractor.mapDF.loc[subtractor.mapDF["order"] == FAINT_ORDER, "flux"] - FITTED_SKY
    )
    assert not model.mask[1].any()


@pytest.mark.parametrize("aggressive", [False, True])
def test_a_skipped_bright_order_is_left_out_of_object_slit_range_pooling_and_the_sky_fit(
    log: Any, monkeypatch: pytest.MonkeyPatch, aggressive: bool
) -> None:
    """Only the faint order reaches the pooled object-slit-range step and the spline fit (DY-1285)."""
    subtractor = _two_order_subtractor(log, monkeypatch, bright=(BRIGHT_ORDER,), aggressive=aggressive)

    subtractor.subtract()

    assert subtractor.clippedOrderIds == [[FAINT_ORDER]]
    assert subtractor.fittedOrderIds == [FAINT_ORDER]


def test_modelled_orders_are_pooled_and_fitted_in_order_when_none_is_bright(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With no skipped order both orders go through, in map order."""
    subtractor = _two_order_subtractor(log, monkeypatch, bright=())

    subtractor.subtract()

    assert subtractor.clippedOrderIds == [[BRIGHT_ORDER, FAINT_ORDER]]
    assert subtractor.fittedOrderIds == [BRIGHT_ORDER, FAINT_ORDER]


def test_a_skipped_bright_order_passes_its_data_through_with_a_masked_zero_model(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A skipped order keeps its measured flux and error, and only its model is zero and flagged in the QUAL mask."""
    subtractor = _two_order_subtractor(log, monkeypatch, bright=(BRIGHT_ORDER,))
    brightFlux = subtractor.mapDF.loc[subtractor.mapDF["order"] == BRIGHT_ORDER, "flux"].to_numpy()

    model, subtracted, residuals, _, _ = subtractor.subtract()

    assert (model.data[0] == 0).all()
    assert model.mask[0].all()
    assert np.allclose(subtracted.data[0], brightFlux)
    assert np.allclose(residuals.data[0], brightFlux / PIXEL_ERROR)
    assert not subtracted.mask[0].any()
    assert (subtracted.uncertainty.array[0] == PIXEL_ERROR).all()
    assert (model.uncertainty.array[0] == PIXEL_ERROR).all()
    assert ("warning", "ORDER 10: OBJECT IS LIKELY VERY BRIGHT - SKIPPING SKY-SUBTRACTION FOR THIS ORDER") in _messages(
        log, "warning"
    )


def test_a_frame_whose_every_order_is_bright_returns_no_model(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """With no order modelled there is no sky model, so the recipe turns sky subtraction off as before."""
    subtractor = _two_order_subtractor(log, monkeypatch, bright=(BRIGHT_ORDER, FAINT_ORDER))

    result = subtractor.subtract()

    assert result == (None, None, None, subtractor.qc, subtractor.products)


def test_a_skipped_order_still_flags_its_model_when_the_object_frame_has_no_mask(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A frame without a mask gets one, flagging only the skipped order's sky-model pixels."""
    subtractor = _two_order_subtractor(log, monkeypatch, bright=(BRIGHT_ORDER,))
    subtractor.objectFrame.mask = None

    model, _, _, _, _ = subtractor.subtract()

    assert model.mask[0].all()
    assert not model.mask[1].any()
