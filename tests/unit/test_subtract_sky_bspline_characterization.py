"""Characterization of the iterative B-spline sky fit in subtract_sky (DY-40)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

PIXEL_COUNT = 3000
# PINNING ALL 3000 MODEL VALUES AS LITERALS IS IMPRACTICAL; THESE ROWS SAMPLE BOTH
# ORDER ENDS (WHERE THE FIT IS ANCHORED), A SKYLINE FLANK AND THE ORDER CENTRE
SAMPLED_ROWS = [0, 750, 1500, 2999]
# THE ORDER-END ANCHORS (DY-593) ARE WEIGHTED 1e5 AGAINST ~0.1 FOR EVERY OTHER SAMPLE, WHICH MAKES THE
# FIT SENSITIVE TO LAST-BIT DIFFERENCES BETWEEN CPU KERNELS; THE PINS BELOW THAT MOVED USE THIS TOLERANCE
ANCHORED_FIT_REL = 1e-9
# THE SKY-SUBTRACTED FLUX AT THE ENDS OF A NOISELESS SLOPE IS ZERO TO WITHIN THIS ABSOLUTE FLUX TOLERANCE (DY-593)
END_RESIDUAL_ABS = 1e-6
# A BRIGHT LINE MAY PULL AN END ANCHOR BY NO MORE THAN THIS FLUX AFTER CLIPPING (DY-593)
END_ANCHOR_LINE_TOLERANCE = 0.5


def _subtractor(
    log: Any,
    *,
    arm: str = "VIS",
    iterationLimit: int = 8,
    pointsPerKnot: int = 100,
    noiseSigma: float = 0,
    debug: bool = False,
    binning: int = 1,
) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = arm
    subtractor.debug = debug
    subtractor.debugInfo = "info"
    subtractor.binx = binning
    subtractor.biny = binning
    subtractor.bspline_order = 3
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "bspline_fitting_residual_clipping_sigma": 3,
            "bspline_iteration_limit": iterationLimit,
            "min_points_per_knot": 4,
            "noise_sigma": noiseSigma,
            "residual_floor_percentile": 50,
            "starting_points_per_knot": pointsPerKnot,
        }
    }
    return subtractor


def _default_like_subtractor(log: Any) -> subtract_sky:
    # THE VIS DEFAULTS FROM soxs_default_settings.yaml FOR THE SETTINGS THAT DRIVE KNOT PLACEMENT
    subtractor = _subtractor(log, iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update({"min_points_per_knot": 155, "residual_floor_percentile": 90})
    return subtractor


def _skyline_order(*, orderNumber: int = 10) -> pd.DataFrame:
    """A sloped continuum with three gaussian skylines and seeded shot noise."""
    rng = np.random.default_rng(42)
    wavelength = np.linspace(500.0, 510.0, PIXEL_COUNT)
    sky = 200.0 + 5.0 * (wavelength - 505.0)
    for centre, amplitude in [(502.0, 2000.0), (504.5, 1500.0), (507.3, 2500.0)]:
        sky = sky + amplitude * np.exp(-0.5 * ((wavelength - centre) / 0.02) ** 2)
    error = np.sqrt(sky)
    return pd.DataFrame(
        {
            "order": orderNumber,
            "wavelength": wavelength,
            "flux": sky + rng.normal(0.0, 1.0, PIXEL_COUNT) * error,
            "error": error,
            "residual_windowed_std": error.copy(),
            "flagged_all_clipped": False,
            "residual_windowed_long_median": np.zeros(PIXEL_COUNT),
            "flux_windowed_long_median": np.full(PIXEL_COUNT, 200.0),
        }
    )


def _sloped_order(*, noisy: bool) -> pd.DataFrame:
    """A line-free sky rising linearly from 175 to 225 across the order."""
    pixels = _skyline_order()
    sky = 200.0 + 5.0 * (pixels["wavelength"] - 505.0)
    pixels["error"] = np.sqrt(sky)
    pixels["flux"] = sky
    if noisy:
        pixels["flux"] = sky + np.random.default_rng(1).normal(0.0, 1.0, PIXEL_COUNT) * pixels["error"]
    pixels["residual_windowed_std"] = pixels["error"].copy()
    return pixels


def _info_messages(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "info"]


def _approx_list(values: list[float], rel: float = 1e-12) -> list[Any]:
    return [pytest.approx(value, rel=rel, abs=0) for value in values]


@pytest.mark.parametrize(
    ("arm", "expected"),
    [
        (
            "VIS",
            {
                "knotCount": 304,
                "knotEnds": [500.3225806451613, 509.89081287486084],
                "coefficientKnots": 312,
                "ratioSum": -667.792329091481,
                "ratioEnds": [-5.895442715411719, 6.615463372189216],
                "modelSum": 690156.4271830036,
                "model": [174.74606062860119, 192.89184103037647, 197.07661758171722, 225.8693315814219],
                "subtractedSum": -1124.2479776855848,
                "subtracted": [4.284967437719104, -5.23086796251323, -5.077067868990042, -5.881142582126529],
                "skyLineCounts": {"line": 2838, "False": 162},
            },
        ),
        (
            "NIR",
            {
                "knotCount": 242,
                "knotEnds": [500.3225806451613, 509.8097783173579],
                "coefficientKnots": 250,
                "ratioSum": 976.3504193849684,
                "ratioEnds": [-0.01088411500870891, 4.824839506859096],
                "modelSum": 689030.2939417951,
                "model": [174.86477575820976, 191.87741473173514, 196.5475479190557, 226.43588573954196],
                "subtractedSum": 1.8852635228096801,
                "subtracted": [4.166252308110529, -4.216441663871905, -4.547998206328515, -6.447696740246585],
                "skyLineCounts": {"line": 2763, "False": 237},
            },
        ),
    ],
)
def test_skylines_drive_knot_insertion_until_an_iteration_adds_none(
    log: Any, arm: str, expected: dict[str, Any]
) -> None:
    """VIS uses flux-boosted weights and 1000-point starter knots; NIR uses 300-point starter knots."""
    subtractor = _subtractor(log, arm=arm)

    modelled, spline, knots, fluxErrorRatio, residualFloor = subtractor.fit_bspline_curve_to_sky(_skyline_order())

    # THE FLOOR IS THE MEASURED ONE, NOT A HARDCODED 5 (DY-1286)
    unclipped = modelled["flagged_all_clipped"] == False  # noqa: E712
    assert residualFloor == pytest.approx(modelled.loc[unclipped, "sky_residual_floor"].median())
    assert knots.size == expected["knotCount"]
    assert knots[[0, -1]].tolist() == _approx_list(expected["knotEnds"], ANCHORED_FIT_REL)
    assert spline[0].size == expected["coefficientKnots"]
    assert spline[2] == 3
    # MORE THAN 2000 UNCLIPPED PIXELS, SO 1000 ARE TRIMMED FROM EACH END OF THE RATIO
    assert fluxErrorRatio.shape == (1000,)
    assert float(fluxErrorRatio.sum()) == pytest.approx(expected["ratioSum"], rel=ANCHORED_FIT_REL, abs=0)
    assert fluxErrorRatio[[0, -1]].tolist() == _approx_list(expected["ratioEnds"], ANCHORED_FIT_REL)
    assert float(modelled["sky_model"].sum()) == pytest.approx(expected["modelSum"], rel=ANCHORED_FIT_REL, abs=0)
    assert int(modelled["sky_model"].isna().sum()) == 0
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(expected["model"], ANCHORED_FIT_REL)
    assert float(modelled["sky_subtracted_flux"].sum()) == pytest.approx(
        expected["subtractedSum"], rel=ANCHORED_FIT_REL, abs=0
    )
    assert modelled["sky_subtracted_flux"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        expected["subtracted"], ANCHORED_FIT_REL
    )
    assert modelled["flagged_sky_line"].astype(str).value_counts().to_dict() == expected["skyLineCounts"]
    assert int(modelled["flagged_noisy_region"].sum()) == 0
    assert (modelled["slit_normalisation_ratio"] == 1).all()
    assert _info_messages(log) == ["\t\tNo new knots added on iteration 8. Stopping iterations.\n"]


def test_debug_binned_order_twelve_prints_knot_budgets_and_plots_every_fit(
    log: Any, capsys: pytest.CaptureFixture[str]
) -> None:
    """2x2 binning divides both knot budgets by four; debug prints them and plots order 12."""
    subtractor = _subtractor(log, debug=True, binning=2)
    plotCalls: list[tuple[str, int]] = []
    subtractor.plot_order_skymodel_fitting_quicklook = lambda frame, spline, title=None, knots=False: plotCalls.append(
        (title, len(knots))
    )

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(_skyline_order(orderNumber=12))

    printed = capsys.readouterr().out.replace("\x1b[1A\x1b[2K", "").splitlines()
    assert "defaultPointsPerKnot: 25.0, order: 12" in printed
    assert "min_points_per_knot: 4, order: 12" in printed
    assert "starterPointsPerKnot: 250.0, order: 12" in printed
    assert "default knot count 120" in printed
    assert "starter knot count 12" in printed
    assert sum(line.startswith("weights2 description") for line in printed) == 3
    assert printed.count("base knot count 12 -2") == 1
    # ONE PLOT PER ITERATION FROM 1 TO 7
    assert len(plotCalls) == 7
    assert plotCalls[:2] == [
        ("Fitting the sky model for order 12\niteration 1. #knots: 12. \ninfo", 12),
        ("Fitting the sky model for order 12\niteration 2. #knots: 28. \ninfo", 28),
    ]
    # THREE PROPOSED KNOTS BOUND AN INTERVAL WITH NO SAMPLE AND ARE DROPPED BEFORE THE FIT (DY-697)
    assert knots.size == 423
    assert spline[0].size == 431
    assert float(fluxErrorRatio.sum()) == pytest.approx(7398.83659145686, rel=ANCHORED_FIT_REL, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(690112.8772945863, rel=ANCHORED_FIT_REL, abs=0)
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [174.3872543237898, 192.2996782789385, 197.28061867741914, 226.28340470007237], ANCHORED_FIT_REL
    )


def test_zero_points_per_knot_drops_the_default_knots_from_iteration_five(log: Any) -> None:
    """A zero knot budget leaves only knots added from residuals once the starters retire."""
    subtractor = _subtractor(log, pointsPerKnot=0)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert knots.size == 274
    assert spline[0].size == 282
    assert float(fluxErrorRatio.sum()) == pytest.approx(3416.0518439631187, rel=ANCHORED_FIT_REL, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(690174.4422460035, rel=ANCHORED_FIT_REL, abs=0)
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [174.7460667297997, 192.1263418789959, 197.89248454496473, 225.8693274403997], ANCHORED_FIT_REL
    )
    assert _info_messages(log) == ["\t\tNo new knots added on iteration 8. Stopping iterations.\n"]


def test_a_poor_fitpack_fit_reverts_the_spline_and_returns_the_knots_of_the_reverted_spline(log: Any) -> None:
    """One default knot per pixel makes FITPACK report ier=10 on iteration 5.

    The spline reverts to iteration 4 and the returned knots are the 75 interior
    knots of that spline, not the 3078 knots of the rejected fit (DY-594).
    """
    subtractor = _subtractor(log, pointsPerKnot=1)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert _info_messages(log) == ["\t\tpoor fit on iteration 5 for order 10. Reverting to last iteration.\n"]
    assert knots.size == 75
    # FITPACK PADS THE INTERIOR KNOTS WITH k + 1 BOUNDARY KNOTS AT EACH END
    assert np.array_equal(spline[0][4:-4], knots)
    assert spline[0].size == 83
    assert float(fluxErrorRatio.sum()) == pytest.approx(-198301.60632199858, rel=ANCHORED_FIT_REL, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(691775.5579021351, rel=ANCHORED_FIT_REL, abs=0)
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [174.74608222469203, 189.76483897934492, 196.2662456790578, 225.8693425945859], ANCHORED_FIT_REL
    )


def _stub_splrep_poor_from_call(monkeypatch: pytest.MonkeyPatch, *, poorCall: int, ier: int) -> list[np.ndarray]:
    """Make the real splrep report ``ier`` from call number ``poorCall``; return the knots of every call."""
    realSplrep = scipy.interpolate.splrep
    knotsPerCall: list[np.ndarray] = []

    def stubbed(*args: Any, **kwargs: Any) -> tuple:
        knotsPerCall.append(np.array(kwargs["t"]))
        spline, residual, _, message = realSplrep(*args, **kwargs)
        if len(knotsPerCall) >= poorCall:
            return spline, residual, ier, "stubbed poor fit"
        return spline, residual, 0, message

    monkeypatch.setattr(scipy.interpolate, "splrep", stubbed)
    return knotsPerCall


@pytest.mark.parametrize("ier", [10, 30, 50, 99])
def test_a_poor_fitpack_fit_on_the_first_iteration_raises_a_value_error_naming_the_fit(
    log: Any, monkeypatch: pytest.MonkeyPatch, ier: int
) -> None:
    """Any ier >= 10 before an accepted fit leaves nothing to revert to, so it is a named pipeline error."""
    _stub_splrep_poor_from_call(monkeypatch, poorCall=1, ier=ier)
    subtractor = _subtractor(log)

    with pytest.raises(ValueError) as raised:
        subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert str(raised.value) == (
        f"BSpline fit failed for order 10 on iteration -2. FITPACK reported ier={ier}: stubbed poor fit"
    )
    assert _info_messages(log) == []


@pytest.mark.parametrize("ier", [10, 30, 50])
def test_a_poor_fitpack_fit_on_a_later_iteration_reverts_to_the_previous_spline_and_knots(
    log: Any, monkeypatch: pytest.MonkeyPatch, ier: int
) -> None:
    """A poor fit on iteration 3 logs, keeps the iteration 2 spline and its knots, and returns normally."""
    knotsPerCall = _stub_splrep_poor_from_call(monkeypatch, poorCall=6, ier=ier)
    subtractor = _subtractor(log)

    modelled, spline, knots, _, _ = subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert len(knotsPerCall) == 6
    # THE REJECTED FIT USED A DIFFERENT KNOT SET, SO MATCHING THE PREVIOUS ONE IS NOT VACUOUS
    assert knotsPerCall[4].size != knotsPerCall[5].size
    assert _info_messages(log) == ["\t\tpoor fit on iteration 3 for order 10. Reverting to last iteration.\n"]
    assert np.array_equal(knots, knotsPerCall[4])
    assert np.array_equal(spline[0][4:-4], knotsPerCall[4])
    assert int(modelled["sky_model"].isna().sum()) == 0


def test_a_fitpack_fit_below_ier_ten_is_accepted(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """FITPACK warnings below 10 (here ier=2) do not stop the iterations."""
    realSplrep = scipy.interpolate.splrep

    def warning_only(*args: Any, **kwargs: Any) -> tuple:
        spline, residual, _, message = realSplrep(*args, **kwargs)
        return spline, residual, 2, message

    monkeypatch.setattr(scipy.interpolate, "splrep", warning_only)
    subtractor = _subtractor(log)

    subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert _info_messages(log) == ["\t\tNo new knots added on iteration 8. Stopping iterations.\n"]


def test_coincident_starter_knots_are_collapsed_so_the_first_fit_succeeds(log: Any) -> None:
    """Coincident starter knots once made the unmocked splrep return ier=30 on the first fit (DY-602).

    The duplicates are now dropped before the fit (DY-697), so the fit runs. The named error for a
    poor first fit is pinned with a stubbed splrep in
    ``test_a_poor_fitpack_fit_on_the_first_iteration_raises_a_value_error_naming_the_fit``.
    """
    worker = _subtractor(log, arm="NIR", iterationLimit=7, pointsPerKnot=25)
    worker.recipeSettings["sky-subtraction"]["min_points_per_knot"] = 5
    rng = np.random.default_rng(1)
    wavelength = np.concatenate([np.full(900, 1000.0), np.linspace(1000.1, 1010.0, 300)])
    imageMapOrder = pd.DataFrame(
        {
            "order": 12,
            "wavelength": wavelength,
            "flux": rng.normal(100.0, 5.0, wavelength.size),
            "error": np.full(wavelength.size, 5.0),
            "residual_windowed_std": np.full(wavelength.size, 5.0),
            "flagged_all_clipped": False,
            "flagged_noisy_region": False,
            "residual_windowed_long_median": np.zeros(wavelength.size),
            "flux_windowed_long_median": np.full(wavelength.size, 100.0),
        }
    )

    modelled, spline, knots, _, _ = worker.fit_bspline_curve_to_sky(imageMapOrder)

    # ALL FOUR STARTER KNOTS FALL ON THE 900 SAMPLES AT 1000.0 NM, BELOW EVERY OTHER SAMPLE, SO NONE SURVIVES
    assert not np.any(knots == 1000.0)
    assert np.all(np.diff(knots) > 0)
    assert np.all(knots > 1000.0)
    assert np.array_equal(spline[0][4:-4], knots)
    assert int(modelled["sky_model"].isna().sum()) == 0


@pytest.mark.parametrize("failure", [ValueError, TypeError, RuntimeError])
def test_a_raising_fitpack_call_becomes_a_value_error_naming_the_knot_budget(
    log: Any, monkeypatch: pytest.MonkeyPatch, failure: type[Exception]
) -> None:
    """The three caught exception types are re-raised as one descriptive ValueError."""

    def raise_failure(*args: Any, **kwargs: Any) -> tuple:
        raise failure("fitpack refused")

    monkeypatch.setattr(scipy.interpolate, "splrep", raise_failure)
    subtractor = _subtractor(log)

    with pytest.raises(ValueError) as raised:
        subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert str(raised.value) == (
        "BSpline fit failed for order 10 on iteration -2. Possibly too many knots (3) "
        "for the number of data points (3000)."
    )
    assert (
        "debug",
        "fit_bspline_curve_to_sky: `tck, fp, ier, msg = ip.splrep(...` failed, continuing: fitpack refused",
    ) in log.messages


def test_a_clean_sloped_sky_is_anchored_to_extrapolated_end_values_at_both_order_ends(log: Any) -> None:
    """The end anchors are line fits to the end windows, evaluated at the end samples (DY-593).

    On a noiseless sky rising from 175 to 225 each end window is a straight
    line, so its fit extrapolates exactly to the end sample and the sky-subtracted
    flux at both ends is zero. The order median would leave -25 and +25 and the
    window medians -6.24 and +6.24.
    """
    subtractor = _default_like_subtractor(log)

    modelled, _, _, _, _ = subtractor.fit_bspline_curve_to_sky(_sloped_order(noisy=False))

    endResiduals = modelled["sky_subtracted_flux"].iloc[[0, -1]].tolist()
    assert endResiduals == pytest.approx([0.0, 0.0], abs=END_RESIDUAL_ABS)
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [175.00000000000585, 187.504168056021, 200.00833611203913, 224.99999999999952], ANCHORED_FIT_REL
    )


def test_a_noisy_line_free_order_follows_the_fitted_sloped_sky(log: Any) -> None:
    """With no extra knots the fitted spline is still the sky model, not a constant median (DY-593)."""
    subtractor = _default_like_subtractor(log)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(_sloped_order(noisy=True))

    assert knots.tolist() == _approx_list([502.5, 505.0, 507.5])
    assert modelled["sky_model"].iloc[-1] - modelled["sky_model"].iloc[0] > 40
    assert modelled["sky_subtracted_flux_weighted"].nunique() > 1
    assert np.unique(fluxErrorRatio).size > 1
    np.testing.assert_allclose(modelled["sky_model"], modelled["sky_model_wl"])
    np.testing.assert_allclose(modelled["sky_model_wl"], scipy.interpolate.splev(modelled["wavelength"].values, spline))


STARTER_KNOTS = np.array([502.5, 505.0, 507.5])


def test_end_anchors_are_line_fits_to_the_windows_outside_the_outer_starter_knots(log: Any) -> None:
    """Each anchor is the straight line fitted to its end window, evaluated at the end sample (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.array([500.0, 501.0, 502.0, 503.0, 506.0, 508.0, 509.0, 510.0])
    flux = np.array([10.0, 30.0, 20.0, 99.0, 99.0, 60.0, 80.0, 70.0])

    blueAnchor, redAnchor = subtractor._end_anchor_values(wavelength, flux, STARTER_KNOTS)

    # BLUE WINDOW: SLOPE 5 THROUGH MEAN 20 AT 501, SO 15 AT 500
    # RED WINDOW: SLOPE 5 THROUGH MEAN 70 AT 509, SO 75 AT 510
    assert (blueAnchor, redAnchor) == pytest.approx((15.0, 75.0))


def test_a_bright_line_inside_an_end_window_does_not_pull_the_end_anchors(log: Any) -> None:
    """The window fits clip skylines, so each anchor stays on the underlying slope at the end sample (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.linspace(500.0, 510.0, 201)
    flux = 200.0 + 5.0 * (wavelength - 505.0)
    # A BRIGHT LINE IN EACH END WINDOW, WELL INSIDE THE FIRST AND LAST STARTER KNOTS
    for centre, amplitude in [(501.0, 2000.0), (509.0, 1500.0)]:
        flux = flux + amplitude * np.exp(-0.5 * ((wavelength - centre) / 0.1) ** 2)

    blueAnchor, redAnchor = subtractor._end_anchor_values(wavelength, flux, STARTER_KNOTS)

    # THE UNDERLYING SLOPE IS 175 AT 500.0 AND 225 AT 510.0
    assert blueAnchor == pytest.approx(175.0, abs=END_ANCHOR_LINE_TOLERANCE)
    assert redAnchor == pytest.approx(225.0, abs=END_ANCHOR_LINE_TOLERANCE)


def test_an_end_anchor_is_the_window_median_when_the_window_has_one_distinct_wavelength(log: Any) -> None:
    """No line fits samples at a single wavelength, so the anchor is their median (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.array([500.0, 500.0, 500.0, 503.0, 506.0, 509.0, 509.0])
    flux = np.array([10.0, 50.0, 30.0, 99.0, 99.0, 60.0, 80.0])

    blueAnchor, redAnchor = subtractor._end_anchor_values(wavelength, flux, STARTER_KNOTS)

    assert (blueAnchor, redAnchor) == (30.0, 70.0)


def test_an_end_anchor_is_the_sample_flux_when_the_window_holds_one_sample(log: Any) -> None:
    """A window with a single unclipped sample anchors the end at that sample (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.array([501.0, 503.0, 506.0, 509.0])
    flux = np.array([12.0, 99.0, 99.0, 88.0])

    blueAnchor, redAnchor = subtractor._end_anchor_values(wavelength, flux, STARTER_KNOTS)

    assert (blueAnchor, redAnchor) == (12.0, 88.0)


def test_an_end_anchor_falls_back_to_the_nearest_unclipped_sample_when_its_window_is_empty(log: Any) -> None:
    """No unclipped sample is bluer than the first knot, so the bluest sample sets the blue anchor (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.array([503.0, 504.0, 506.0, 508.0, 509.0, 510.0])
    flux = np.array([41.0, 99.0, 99.0, 60.0, 80.0, 70.0])

    blueAnchor, redAnchor = subtractor._end_anchor_values(wavelength, flux, STARTER_KNOTS)

    assert blueAnchor == 41.0
    assert redAnchor == pytest.approx(75.0)


def test_an_end_anchor_falls_back_to_the_end_sample_when_the_red_window_is_empty(log: Any) -> None:
    """No unclipped sample is redder than the last knot, so the reddest sample sets the red anchor (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.array([500.0, 501.0, 502.0, 504.0, 506.0, 507.0])
    flux = np.array([10.0, 30.0, 20.0, 99.0, 99.0, 55.0])

    blueAnchor, redAnchor = subtractor._end_anchor_values(wavelength, flux, STARTER_KNOTS)

    assert blueAnchor == pytest.approx(15.0)
    assert redAnchor == 55.0


def test_end_anchors_are_the_end_samples_when_there_are_no_starter_knots(log: Any) -> None:
    """Without starter knots no window is defined, so each end keeps its own sample (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.array([500.0, 501.0, 502.0])
    flux = np.array([10.0, 30.0, 20.0])

    anchors = subtractor._end_anchor_values(wavelength, flux, np.array([]))

    assert anchors == (10.0, 20.0)


def test_blue_end_noise_prunes_only_the_first_knot(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Noise bluer than the first knot removes that knot and keeps the reddest knot (DY-601)."""
    pixels = _skyline_order()
    pixels.loc[200:599, "residual_windowed_long_median"] = 50.0
    subtractor = _subtractor(log, noiseSigma=3)
    # THE PRUNED KNOT SETS ARE LOCAL; RECORD THE BINS EACH np.digitize CALL RECEIVES
    realDigitize = np.digitize
    knotSets: list[list[float]] = []

    def record_digitize(values: Any, bins: Any, *args: Any, **kwargs: Any) -> Any:
        knotSets.append(np.asarray(bins).tolist())
        return realDigitize(values, bins, *args, **kwargs)

    monkeypatch.setattr(np, "digitize", record_digitize)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(pixels)

    noisyWavelengths = modelled.loc[modelled["flagged_noisy_region"], "wavelength"]
    assert int(modelled["flagged_noisy_region"].sum()) == 400
    assert [noisyWavelengths.min(), noisyWavelengths.max()] == _approx_list([500.6668889629877, 501.99733244414807])
    # FIRST PRUNING PASS: THE STARTER KNOTS, THEN THE KNOTS LEFT AFTER REMOVAL
    assert knotSets[0] == [502.5, 505.0, 507.5]
    assert knotSets[1] == [505.0, 507.5]
    # THE STARTER KNOT AT 502.5 IS GONE, SO THE FIT IS ILL-CONDITIONED. ITERATION 8 ONCE PROPOSED A
    # DUPLICATE KNOT, FITPACK REJECTED IT (ier=30) AND THE FIT REVERTED; THAT KNOT IS NOW DROPPED BEFORE
    # THE FIT (DY-697), SO NO ITERATION IS REJECTED AND THE KNOTS ARE THOSE OF THE RETURNED SPLINE
    assert _info_messages(log) == []
    assert knots.size == 223
    assert np.array_equal(spline[0][4:-4], knots)
    assert knots[[0, -1]].tolist() == _approx_list([500.0662038828806, 509.6774193548387])
    # THE DIVERGING FIT IS ILL-CONDITIONED, SO THESE NUMBERS ARE PINNED AT 1e-6, NOT 1e-12
    # THE SKY MODEL INSIDE THE DOWN-WEIGHTED NOISY BLOCK STILL HUMPS FAR ABOVE THE TRUE ~200 SKY, BECAUSE THE
    # BLOCK SHARES A KNOT INTERVAL WITH THE BRIGHT 502.0 NM SKYLINE (DY-695). THE HUMP IS HIGHER THAN BEFORE
    # THE KNOT FIX ONLY BECAUSE THE ITERATION-8 FIT IS NO LONGER REJECTED. THE REAL VIS AND NIR STARE FRAMES
    # SHOW NO SUCH HUMP, SO THESE VALUES CHARACTERISE A SYNTHETIC-ONLY FIT AND MUST FLIP IF THE HUMP IS FIXED
    assert float(fluxErrorRatio.sum()) == pytest.approx(-1293.940961304414, rel=1e-6, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(1300539.7029378437, rel=1e-6, abs=0)
    assert modelled["sky_model"].iloc[[0, 400, 1500, 2999]].tolist() == [
        pytest.approx(value, rel=1e-6, abs=0)
        for value in [174.74607314171362, 2379.2515818826296, 195.99985563512843, 225.86936119293233]
    ]


@pytest.mark.parametrize(
    ("knots", "noisyWavelengths", "expected"),
    [
        pytest.param([502.5, 505.0, 507.5], [509.0], [502.5, 505.0], id="red-of-the-last-knot-removes-only-the-last"),
        pytest.param([502.5, 505.0, 507.5], [506.0], [502.5], id="red-end-interval-removes-both-bounding-knots"),
        pytest.param(
            [501.0, 503.0, 505.0, 507.0, 509.0],
            [504.0],
            [501.0, 507.0, 509.0],
            id="middle-noise-keeps-every-distant-knot",
        ),
        pytest.param(
            [502.5, 505.0, 507.5], [500.0], [505.0, 507.5], id="blue-of-the-first-knot-removes-only-the-first"
        ),
        pytest.param(
            [501.0, 503.0, 505.0, 507.0, 509.0],
            [500.0, 504.0, 510.0],
            [507.0],
            id="noise-at-both-ends-and-the-middle",
        ),
        pytest.param([502.5, 505.0, 507.5], [], [502.5, 505.0, 507.5], id="no-noise-keeps-every-knot"),
    ],
)
def test_pruning_removes_only_the_knots_bounding_noisy_pixels(
    log: Any, knots: list[float], noisyWavelengths: list[float], expected: list[float]
) -> None:
    """Each noisy pixel removes the knots of its knot interval; no index wraps to the other end (DY-601)."""
    subtractor = _subtractor(log)

    pruned = subtractor._prune_knots_in_noise(np.array(knots), np.array(noisyWavelengths), order=10)

    assert pruned.tolist() == expected
    assert log.messages == []


def test_pruning_never_mutates_the_knots_it_is_given(log: Any) -> None:
    """The caller's knot array is left intact; pruning returns a new array."""
    subtractor = _subtractor(log)
    knots = np.array([502.5, 505.0, 507.5])

    subtractor._prune_knots_in_noise(knots, np.array([503.0]), order=10)

    assert knots.tolist() == [502.5, 505.0, 507.5]


@pytest.mark.parametrize(
    ("knots", "noisyWavelengths"),
    [
        pytest.param([502.5], [503.0], id="one-knot-and-noise-redward-of-it"),
        pytest.param([502.5, 505.0], [501.0, 506.0], id="noise-either-side-of-every-knot"),
    ],
)
def test_pruning_that_would_remove_every_knot_leaves_the_knots_and_names_the_order(
    log: Any, knots: list[float], noisyWavelengths: list[float]
) -> None:
    """An empty knot vector would break splrep, so the knots stay and the order is logged (DY-601)."""
    subtractor = _subtractor(log)

    pruned = subtractor._prune_knots_in_noise(np.array(knots), np.array(noisyWavelengths), order=13)

    assert pruned.tolist() == knots
    assert log.messages == [
        ("warning", "\t\tNoisy-region pruning would remove every knot for order 13. Keeping the knots unchanged.\n")
    ]


@pytest.mark.parametrize(
    ("knots", "expected"),
    [
        pytest.param([502.0, 505.0, 508.0], [502.0, 505.0, 508.0], id="valid-knots-are-kept"),
        pytest.param([502.0, 505.0, 505.0, 508.0], [502.0, 505.0, 508.0], id="an-exact-duplicate-is-dropped"),
        pytest.param([505.0, 502.0, 508.0], [502.0, 505.0, 508.0], id="knots-are-returned-sorted"),
        pytest.param([499.9, 500.0, 505.0], [505.0], id="knots-at-or-below-the-first-sample-are-dropped"),
        pytest.param([505.0, 510.0, 510.1], [505.0], id="knots-at-or-above-the-last-sample-are-dropped"),
        pytest.param(
            [502.0, 502.05, 502.4, 508.0], [502.0, 502.4, 508.0], id="a-knot-with-no-sample-since-the-last-is-dropped"
        ),
        pytest.param([502.3, 503.0, 508.0], [502.3, 508.0], id="a-sample-on-a-knot-counts-for-neither-interval"),
        pytest.param([502.0, np.nan, 508.0], [502.0, 508.0], id="a-nan-knot-is-dropped"),
        pytest.param([], [], id="no-knots"),
    ],
)
def test_knots_without_samples_in_their_interval_are_dropped(
    log: Any, knots: list[float], expected: list[float]
) -> None:
    """Every kept knot has a sample between it and the previous kept knot, and a sample after it (DY-697)."""
    subtractor = _subtractor(log)
    # NO SAMPLE LIES BETWEEN 502.0 AND 502.05 NM, SO A KNOT AT 502.05 AFTER ONE AT 502.0 BOUNDS AN EMPTY INTERVAL
    wavelength = np.array([500.0, 501.0, 501.5, 502.1, 502.3, 503.0, 505.5, 507.0, 509.0, 510.0])

    kept = subtractor._drop_knots_without_samples(np.array(knots), wavelength, order=10)

    assert kept.tolist() == expected
    assert log.messages == []


def test_dropping_knots_never_mutates_the_knots_it_is_given(log: Any) -> None:
    """The caller's knot array is left intact; a new array is returned."""
    subtractor = _subtractor(log)
    knots = np.array([505.0, 505.0, 499.0])

    subtractor._drop_knots_without_samples(knots, np.linspace(500.0, 510.0, 11), order=10)

    assert knots.tolist() == [505.0, 505.0, 499.0]


def test_tied_sample_wavelengths_keep_only_one_knot_among_them(log: Any) -> None:
    """Samples sharing one wavelength support no interval between knots placed on and among them (DY-697)."""
    subtractor = _subtractor(log)
    wavelength = np.concatenate([np.full(5, 1000.0), [1001.0, 1002.0, 1003.0]])

    kept = subtractor._drop_knots_without_samples(np.array([1000.0, 1000.0, 1000.5, 1001.5]), wavelength, order=10)

    # 1000.0 HAS NO SAMPLE BLUEWARD OF IT; 1000.5 HAS THE TIED SAMPLES; 1001.5 HAS THE SAMPLE AT 1001.0
    assert kept.tolist() == [1000.5, 1001.5]


def test_dropping_every_knot_warns_and_names_the_order(log: Any) -> None:
    """With no knot left the spline is one cubic across the order, so the loss is logged (DY-697)."""
    subtractor = _subtractor(log)

    wavelength = np.linspace(500.0, 510.0, 11)

    kept = subtractor._drop_knots_without_samples(np.array([499.0, 510.0, 511.0]), wavelength, order=13)

    assert kept.tolist() == []
    assert log.messages == [
        (
            "warning",
            "\t\tEvery proposed b-spline knot for order 13 lacks samples in its interval. Fitting without knots.\n",
        )
    ]


def test_a_fit_left_with_no_knots_is_one_cubic_across_the_order(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """When every knot is dropped the fit still runs, as a single cubic with only boundary knots (DY-697)."""
    subtractor = _subtractor(log)
    monkeypatch.setattr(subtractor, "_drop_knots_without_samples", lambda knots, wavelength, order: np.array([]))

    modelled, spline, knots, _, _ = subtractor.fit_bspline_curve_to_sky(_sloped_order(noisy=False))

    assert knots.size == 0
    # FITPACK PADS AN EMPTY INTERIOR KNOT VECTOR WITH k + 1 BOUNDARY KNOTS AT EACH END
    assert spline[0].size == 8
    assert int(modelled["sky_model"].isna().sum()) == 0


def test_an_empty_sample_list_keeps_no_knot(log: Any) -> None:
    """With no samples no knot interval can hold one."""
    subtractor = _subtractor(log)

    kept = subtractor._drop_knots_without_samples(np.array([505.0]), np.array([]), order=10)

    assert kept.tolist() == []


@pytest.mark.parametrize(
    ("subtractorFactory", "stopIteration"),
    [
        pytest.param(lambda log: _subtractor(log, noiseSigma=3, iterationLimit=25), 11, id="test-settings"),
        pytest.param(_default_like_subtractor, 9, id="default-like-settings"),
    ],
)
def test_a_dropped_proposal_leaves_the_extra_knots_so_knot_addition_stops_sooner(
    log: Any, subtractorFactory: Any, stopIteration: int
) -> None:
    """A dropped proposal is removed from the extra knots, so re-proposing it is not counted as a new
    knot. Keeping it made the "no new knots" stop fire on iterations 15 and 12 instead (DY-697)."""
    pixels = _skyline_order()
    pixels.loc[200:599, "residual_windowed_long_median"] = 50.0
    subtractor = subtractorFactory(log)

    subtractor.fit_bspline_curve_to_sky(pixels)

    assert _info_messages(log) == [f"\t\tNo new knots added on iteration {stopIteration}. Stopping iterations.\n"]


@pytest.mark.parametrize(
    "subtractorFactory",
    [
        pytest.param(lambda log: _subtractor(log, noiseSigma=3), id="test-settings"),
        pytest.param(_default_like_subtractor, id="default-like-settings"),
    ],
)
def test_knot_addition_next_to_noise_never_hands_fitpack_a_knot_without_samples(
    log: Any, monkeypatch: pytest.MonkeyPatch, subtractorFactory: Any
) -> None:
    """Knot addition beside a noisy block once made duplicate knots and empty intervals, so FITPACK
    returned ier=30 and the fit reverted (DY-697). Every knot vector FITPACK sees now has a sample in
    every interval, and no fit is rejected."""
    pixels = _skyline_order()
    pixels.loc[200:599, "residual_windowed_long_median"] = 50.0
    subtractor = subtractorFactory(log)
    realSplrep = scipy.interpolate.splrep
    calls: list[tuple[np.ndarray, np.ndarray, int]] = []

    def recording_splrep(*args: Any, **kwargs: Any) -> tuple:
        result = realSplrep(*args, **kwargs)
        calls.append((np.asarray(args[0]), np.asarray(kwargs["t"]), result[2]))
        return result

    monkeypatch.setattr(scipy.interpolate, "splrep", recording_splrep)

    subtractor.fit_bspline_curve_to_sky(pixels)

    assert len(calls) > 3
    for wavelength, knots, _ in calls:
        edges = np.concatenate(([-np.inf], knots, [np.inf]))
        samplesPerInterval = np.histogram(wavelength, bins=edges)[0]
        assert np.all(np.diff(knots) > 0)
        assert knots.size == 0 or (knots[0] > wavelength[0] and knots[-1] < wavelength[-1])
        assert np.all(samplesPerInterval > 0)
    assert [ier for _, _, ier in calls if ier >= 10] == []
    assert not any("poor fit" in message for message in _info_messages(log))


def test_a_rejected_anchor_sample_stays_rejected_in_later_clip_passes(log: Any) -> None:
    """Rejections accumulate across passes, so a refit line cannot re-admit a rejected sample (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.arange(20.0)
    # A FLAT SKY OF +-1 NOISE WITH A BRIGHT OUTLIER AT 3 AND A FAINTER ONE AT 0
    flux = np.tile([1.0, -1.0], 10)
    flux[3] += 40.0
    flux[0] += 8.0

    value = subtractor._clipped_line_value(wavelength, flux, wavelength[0])

    # BOTH OUTLIERS REJECTED LEAVES THE FLAT SKY, SO THE LINE IS NEAR 0 (RE-ADMITTING THE ONE AT 0 GIVES 1.99)
    assert value == pytest.approx(0.0, abs=0.5)


def test_an_end_anchor_is_clamped_to_the_flux_range_of_its_surviving_window_samples(log: Any) -> None:
    """A noisy short window whose line extrapolates beyond its data cannot set a wilder anchor (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.array([500.0, 501.0, 502.0, 503.0, 506.0, 508.0, 509.0, 510.0])
    flux = np.array([1.0, 0.0, 6.0, 8.0, 99.0, 60.0, 80.0, 70.0])

    blueAnchor, redAnchor = subtractor._end_anchor_values(wavelength, flux, STARTER_KNOTS)

    # THE UNCLAMPED BLUE LINE IS -1/6 AT 500, BELOW EVERY SAMPLE IN THE WINDOW (MINIMUM 0)
    assert blueAnchor == 0.0
    assert redAnchor == pytest.approx(75.0)


def test_the_clamp_uses_the_flux_range_of_the_survivors_and_not_of_the_whole_window(log: Any) -> None:
    """A clipped low outlier must not widen the clamp, or the line could fall below every survivor (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.arange(10.0)
    # A RISING SKY WITH ONE CLIPPED OUTLIER AT -40, WELL BELOW EVERY SURVIVOR (MINIMUM 5)
    flux = np.array([5.0, 6.0, 5.5, 7.0, 6.5, -40.0, 8.0, 7.5, 9.0, 8.5])

    value = subtractor._clipped_line_value(wavelength, flux, -10.0)

    # THE UNCLAMPED LINE IS 1.2 AT -10: INSIDE THE WHOLE-WINDOW RANGE (-40 TO 9) BUT BELOW THE SURVIVORS' MINIMUM
    assert value == 5.0


def test_the_clamp_stops_a_line_that_rises_above_the_maximum_surviving_flux(log: Any) -> None:
    """The upper bound binds as well as the lower one (DY-593)."""
    subtractor = _subtractor(log)
    wavelength = np.arange(8.0)
    flux = np.array([0.0, 1.2, 1.9, 3.1, 4.0, 4.9, 6.1, 7.0])

    value = subtractor._clipped_line_value(wavelength, flux, 11.0)

    # NO SAMPLE IS CLIPPED AND THE UNCLAMPED LINE IS 10.97 AT 11, ABOVE THE MAXIMUM SAMPLE FLUX OF 7
    assert value == 7.0


@pytest.mark.parametrize(
    ("wavelength", "flux", "survivors", "evaluationWavelength"),
    [
        pytest.param(
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [1.0, 55.0, 3.0, 35.0, -2.0],
            [0, 2, 4],
            2.0,
            id="second-pass-would-leave-two-samples",
        ),
        pytest.param(
            [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
            [0.1, 0.0, -0.2, 2.0, -6.6, -0.1],
            [0, 1, 2, 3, 5],
            1.5,
            id="second-pass-would-leave-one-distinct-wavelength",
        ),
    ],
)
def test_a_refused_second_clip_pass_keeps_the_fit_of_the_first_applied_pass(
    log: Any,
    wavelength: list[float],
    flux: list[float],
    survivors: list[int],
    evaluationWavelength: float,
) -> None:
    """Pass one clips and is applied; pass two is refused as unfittable, so the pass-one line is the result (DY-593)."""
    subtractor = _subtractor(log)
    wavelengths = np.array(wavelength)
    fluxes = np.array(flux)
    expectedLine = np.polyfit(wavelengths[survivors], fluxes[survivors], 1)

    value = subtractor._clipped_line_value(wavelengths, fluxes, evaluationWavelength)

    assert value == pytest.approx(np.polyval(expectedLine, evaluationWavelength), abs=1e-9)
