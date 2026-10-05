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
    subtractor.recipeSettings["sky-subtraction"].update(
        {"min_points_per_knot": 155, "residual_floor_percentile": 90}
    )
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


def _approx_list(values: list[float]) -> list[Any]:
    return [pytest.approx(value, rel=1e-12, abs=0) for value in values]


@pytest.mark.parametrize(
    ("arm", "expected"),
    [
        (
            "VIS",
            {
                "knotCount": 310,
                "knotEnds": [500.0211085701805, 509.9788914298195],
                "coefficientKnots": 318,
                "ratioSum": -667.7923290917166,
                "ratioEnds": [-5.895442715411514, 6.615463372189558],
                "modelSum": 690148.371412821,
                "model": [200.49696400397752, 192.8918410303765, 197.07661758171733, 200.49698499613217],
                "subtractedSum": -1116.1922075029995,
                "subtracted": [-21.465935937657235, -5.230867962513258, -5.077067868990156, 19.491204003163205],
                "skyLineCounts": {"line": 2839, "False": 161},
            },
        ),
        (
            "NIR",
            {
                "knotCount": 244,
                "knotEnds": [500.03872331265336, 509.96127668734664],
                "coefficientKnots": 252,
                "ratioSum": 976.3504193848743,
                "ratioEnds": [-0.010884115008691206, 4.8248395068587575],
                "modelSum": 689028.8674091145,
                "model": [200.49697134115678, 191.87741473173512, 196.5475479190558, 200.49697134118642],
                "subtractedSum": 3.311796203490644,
                "subtracted": [-21.46594327483649, -4.216441663871876, -4.547998206328629, 19.49121765810895],
                "skyLineCounts": {"line": 2783, "False": 217},
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

    assert residualFloor == 5
    assert knots.size == expected["knotCount"]
    assert knots[[0, -1]].tolist() == _approx_list(expected["knotEnds"])
    assert spline[0].size == expected["coefficientKnots"]
    assert spline[2] == 3
    # MORE THAN 2000 UNCLIPPED PIXELS, SO 1000 ARE TRIMMED FROM EACH END OF THE RATIO
    assert fluxErrorRatio.shape == (1000,)
    assert float(fluxErrorRatio.sum()) == pytest.approx(expected["ratioSum"], rel=1e-12, abs=0)
    assert fluxErrorRatio[[0, -1]].tolist() == _approx_list(expected["ratioEnds"])
    assert float(modelled["sky_model"].sum()) == pytest.approx(expected["modelSum"], rel=1e-12, abs=0)
    assert int(modelled["sky_model"].isna().sum()) == 0
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(expected["model"])
    assert float(modelled["sky_subtracted_flux"].sum()) == pytest.approx(expected["subtractedSum"], rel=1e-12, abs=0)
    assert modelled["sky_subtracted_flux"].iloc[SAMPLED_ROWS].tolist() == _approx_list(expected["subtracted"])
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
    subtractor.plot_order_skymodel_fitting_quicklook = lambda frame, spline, title=None, knots=False: (
        plotCalls.append((title, len(knots)))
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
    assert knots.size == 432
    assert spline[0].size == 440
    assert float(fluxErrorRatio.sum()) == pytest.approx(7844.2200848747725, rel=1e-12, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(690107.1687027367, rel=1e-12, abs=0)
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [200.49695792830985, 193.79577636769793, 197.2814518546773, 200.4969892156255]
    )


def test_zero_points_per_knot_drops_the_default_knots_from_iteration_five(log: Any) -> None:
    """A zero knot budget leaves only knots added from residuals once the starters retire."""
    subtractor = _subtractor(log, pointsPerKnot=0)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert knots.size == 280
    assert spline[0].size == 288
    assert float(fluxErrorRatio.sum()) == pytest.approx(3416.0518439632997, rel=1e-12, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(690166.1062583625, rel=1e-12, abs=0)
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [200.49696380443177, 192.12634187899585, 197.89248454496467, 200.49698505401028]
    )
    assert _info_messages(log) == ["\t\tNo new knots added on iteration 8. Stopping iterations.\n"]


def test_a_poor_fitpack_fit_reverts_the_spline_but_returns_the_rejected_knots(log: Any) -> None:
    """One default knot per pixel makes FITPACK report ier=10 on iteration 5.

    The spline reverts to iteration 4, but the returned knot array is the
    3078-knot set from the rejected fit, not the knots of the returned spline.
    """
    # SUSPICIOUS: RETURNED KNOTS DO NOT MATCH THE REVERTED SPLINE, FILED AS DY-594
    subtractor = _subtractor(log, pointsPerKnot=1)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert _info_messages(log) == ["\t\tpoor fit on iteration 5 for order 10. Reverting to last iteration.\n"]
    assert knots.size == 3078
    assert spline[0].size == 89
    assert float(fluxErrorRatio.sum()) == pytest.approx(-198301.60632865055, rel=1e-12, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(691899.9408698891, rel=1e-12, abs=0)
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [200.49687793408, 189.76482499233046, 196.26624567950512, 200.4969850529576]
    )


def test_a_poor_fitpack_fit_on_the_first_iteration_has_no_spline_to_revert_to(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """ier=10 before any accepted fit raises UnboundLocalError, not a pipeline error."""
    # SUSPICIOUS: CRASHES ON AN UNSET tck_previous, FILED AS DY-594
    realSplrep = scipy.interpolate.splrep

    def invalid_input(*args: Any, **kwargs: Any) -> tuple:
        spline, residual, _, message = realSplrep(*args, **kwargs)
        return spline, residual, 10, message

    monkeypatch.setattr(scipy.interpolate, "splrep", invalid_input)
    subtractor = _subtractor(log)

    with pytest.raises(UnboundLocalError, match="tck_previous"):
        subtractor.fit_bspline_curve_to_sky(_skyline_order())

    assert _info_messages(log) == ["\t\tpoor fit on iteration -2 for order 10. Reverting to last iteration.\n"]


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


def test_a_clean_sloped_sky_is_pinned_to_the_median_flux_at_both_order_ends(log: Any) -> None:
    """The first and last fit samples are replaced by the median flux at weight 1e5.

    On a noiseless sky rising from 175 to 225 the model is 200 at both ends, so
    the sky-subtracted flux is -25 at the blue end and +25 at the red end.
    """
    # SUSPICIOUS (WRONG SCIENCE): ORDER-END SKY FORCED TO THE MEDIAN, FILED AS DY-593
    subtractor = _default_like_subtractor(log)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(_sloped_order(noisy=False))

    assert knots.size == 19
    assert spline[0].size == 27
    assert modelled["sky_model"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [199.99992780248533, 187.49891924920803, 200.00818688926358, 200.00005908596637]
    )
    assert modelled["sky_subtracted_flux"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [-24.99992780248533, 0.00524880681069817, 0.0001492227738708607, 24.99994091403363]
    )
    assert float(modelled["sky_model"].sum()) == pytest.approx(599999.9075715989, rel=1e-12, abs=0)
    # CPU-SENSITIVE: THE RATIO CANCELS ALMOST COMPLETELY, SO ONE-ULP MODEL DIFFERENCES REACH 1.3e-9
    # RELATIVE IN THE SUM. PINNED TO THE AVX2 VALUE; THE AVX512 VALUE IS -0.0005019052654203701
    assert float(fluxErrorRatio.sum()) == pytest.approx(-0.000501905266075315, rel=1e-12, abs=0)


def test_a_noisy_line_free_order_collapses_to_a_constant_median_sky(log: Any) -> None:
    """With default-like settings no extra knots are added, so the fitted spline is discarded.

    The model becomes the median flux everywhere and the weighted flux-error
    ratio becomes exactly 1 for every pixel.
    """
    # SUSPICIOUS (WRONG SCIENCE): FLAT MEDIAN SKY WHEN NO KNOTS ARE ADDED, FILED AS DY-593
    subtractor = _default_like_subtractor(log)
    pixels = _sloped_order(noisy=True)

    modelled, spline, knots, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(pixels)

    assert knots.tolist() == _approx_list([502.5, 505.0, 507.5])
    assert spline[0].size == 11
    assert modelled["sky_model"].nunique() == 1
    assert modelled["sky_model"].iloc[0] == pytest.approx(199.86817357662824, rel=1e-12, abs=0)
    assert modelled["sky_model"].iloc[0] == float(np.median(pixels["flux"]))
    assert modelled["sky_subtracted_flux"].iloc[SAMPLED_ROWS].tolist() == _approx_list(
        [-20.296524430435227, -28.6797163555143, 5.63487821153646, 14.574571746419537]
    )
    assert (modelled["sky_subtracted_flux_weighted"] == 1).all()
    assert fluxErrorRatio.tolist() == [1] * 1000
    assert modelled["flagged_sky_line"].astype(str).value_counts().to_dict() == {
        "False": 2733,
        "peak": 207,
        "line": 60,
    }


def test_blue_end_noise_prunes_the_reddest_knot(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Noise bluer than the first knot maps to index -1, so knot pruning removes the last knot."""
    # SUSPICIOUS (WRONG SCIENCE): PRUNING DROPS THE REDDEST KNOT, NOT THE ONE IN THE NOISE, FILED AS DY-601
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
    assert knotSets[1] == [502.5, 505.0]
    assert knots.size == 21
    assert knots[[0, -1]].tolist() == _approx_list([499.9479254917661, 509.4736600214474])
    assert float(fluxErrorRatio.sum()) == pytest.approx(-2644.833214674917, rel=1e-12, abs=0)
    assert float(modelled["sky_model"].sum()) == pytest.approx(693026.6083934802, rel=1e-12, abs=0)
    assert modelled["sky_model"].iloc[[0, 400, 1500, 2999]].tolist() == _approx_list(
        [200.49679154663622, 321.9571447708863, 200.46468124773114, 200.49817069183024]
    )
