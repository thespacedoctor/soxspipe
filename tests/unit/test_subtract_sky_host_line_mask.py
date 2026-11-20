"""Host emission lines that fill the slit stay out of the stare sky fit and the sky is interpolated across (DY-2242).

The synthetic order has a flat sky and, in some tests, a host line that fills the slit, is strongest at the slit centre
and shifts along the slit (rotation). The tests run the real ``_fit_sky_to_orders``, because the point is that the
spline is fitted without the line pixels. The skyline catalogue is the one external input and is replaced.
"""

from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml

import soxspipe
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.nod_host_contamination import host_feature
from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

ORDER = 1
COLUMNS = 1200
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
FIRST_NM = 600.0
STEP_NM = 0.04
SKY_LEVEL = 20.0
READ_NOISE = 3.3
HOST_NM = 620.0
HOST_SIGMA_NM = 0.05
HOST_PEAK = 150.0
HOST_SLIT_SIGMA_ARCSEC = 3.0
HOST_ROTATION_NM_PER_ARCSEC = 0.04
# THE SKY MODEL AT THE LINE CENTRE MAY SIT THIS FAR ABOVE THE TRUE SKY WHEN THE LINE IS LEFT OUT OF THE FIT (E-)
MAX_MODEL_EXCESS = 4.0
# AND MUST SIT AT LEAST THIS FAR ABOVE IT WHEN THE LINE IS IN THE FIT
MIN_UNMASKED_EXCESS = 15.0

DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"


def _sky_settings(**overrides: Any) -> dict[str, Any]:
    """The shipped SOXS VIS stare sky-subtraction settings, with the QC plot and the wing correction off."""
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    skySettings = copy.deepcopy(settings["soxs-stare"]["vis"]["sky-subtraction"])
    skySettings["sky_model_qc_plot"] = False
    skySettings["object_wing_correction"] = False
    skySettings.update(overrides)
    return skySettings


def _subtractor(log: Any, **overrides: Any) -> subtract_sky:
    """A subtract_sky with only the attributes the per-order steps read."""
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.settings = {"instrument": "soxs"}
    subtractor.debug = False
    subtractor.binx = 1
    subtractor.biny = 1
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.ron = READ_NOISE
    subtractor.recipeSettings = {"sky-subtraction": _sky_settings(**overrides)}
    subtractor.bspline_order = subtractor.recipeSettings["sky-subtraction"]["bspline_order"]
    subtractor.qcPlotOrder = -1
    subtractor.qc = pd.DataFrame()
    subtractor.recipeName = "soxs-stare"
    subtractor.dateObs = "2026-01-01T00:00:00"
    return subtractor


def _order(hasHost: bool, seed: int = 2242) -> pd.DataFrame:
    """One synthetic order on a flat sky, with an optional host line that fills the slit."""
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slit = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)[slitPixel.ravel()]
    wavelength = FIRST_NM + STEP_NM * column
    noiseless = np.full(column.shape, SKY_LEVEL)
    if hasHost:
        centre = HOST_NM + HOST_ROTATION_NM_PER_ARCSEC * slit
        lineProfile = np.exp(-0.5 * ((wavelength - centre) / HOST_SIGMA_NM) ** 2)
        noiseless = noiseless + HOST_PEAK * lineProfile * np.exp(-0.5 * (slit / HOST_SLIT_SIGMA_ARCSEC) ** 2)
    error = np.sqrt(noiseless + READ_NOISE**2)
    frame = pd.DataFrame(
        {
            "order": float(ORDER),
            "x": column,
            "y": slitPixel.ravel(),
            "wavelength": wavelength + 0.002 * slit,
            "slit_position": slit,
            "flux": noiseless + rng.standard_normal(noiseless.size) * error,
            "error": error,
            "mask": False,
        }
    )
    return frame.sort_values("wavelength", kind="stable").reset_index(drop=True)


@pytest.fixture
def no_catalogue_skylines(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        toolkit, "get_skylines_dataframe", lambda *args, **kwargs: pd.DataFrame({"WAVELENGTH": [640.0, 700.0]})
    )


def _fit(subtractor: subtract_sky, order: pd.DataFrame) -> pd.DataFrame:
    """Run the clipping and the sky fit of ``subtract`` on one order and return the fitted order."""
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = subtractor.get_over_sampled_sky_from_order(
        order.copy(), clipBPs=True, clipSlitEdge=settings["clip-slit-edge-fraction"]
    )
    clipped = subtractor.clip_object_slit_positions([clipped], aggressive=settings["aggressive_object_masking"])
    fittedOrders, _, _ = subtractor._fit_sky_to_orders(clipped)
    return fittedOrders[0]


def _model_at_line_centre(fitted: pd.DataFrame) -> float:
    isCentre = (fitted["wavelength"] - HOST_NM).abs() < HOST_SIGMA_NM
    return float(fitted.loc[isCentre, "sky_model"].mean())


def test_the_pixels_of_a_host_line_are_left_out_of_the_sky_fit(log: Any, no_catalogue_skylines: None) -> None:
    fitted = _fit(_subtractor(log), _order(hasHost=True))

    isInLine = (fitted["wavelength"] - HOST_NM).abs() < 2 * HOST_SIGMA_NM
    assert fitted.loc[isInLine, "flagged_all_clipped"].all()


def test_the_detection_does_not_depend_on_the_dispersion_axis_the_subtractor_reports(
    log: Any, no_catalogue_skylines: None
) -> None:
    # THE SOXS VIS DETECTOR REPORTS THE DISPERSION AXIS AS x, YET ITS ORDER TABLE LISTS THE WAVELENGTH ALONG y (DY-2242)
    subtractor = _subtractor(log)
    subtractor.axisA = "y"

    _fit(subtractor, _order(hasHost=True))

    assert len(subtractor.hostLineWindows) == 1


def test_the_sky_model_interpolates_across_a_host_line_instead_of_taking_its_light(
    log: Any, no_catalogue_skylines: None
) -> None:
    fitted = _fit(_subtractor(log), _order(hasHost=True))

    assert _model_at_line_centre(fitted) - SKY_LEVEL < MAX_MODEL_EXCESS


def test_with_the_mask_switched_off_the_sky_model_takes_the_light_of_a_host_line(
    log: Any, no_catalogue_skylines: None
) -> None:
    fitted = _fit(_subtractor(log, **{"host-line-mask": False}), _order(hasHost=True))

    assert _model_at_line_centre(fitted) - SKY_LEVEL > MIN_UNMASKED_EXCESS


def test_an_order_with_no_host_line_is_fitted_exactly_as_with_the_mask_switched_off(
    log: Any, no_catalogue_skylines: None
) -> None:
    order = _order(hasHost=False)

    masked = _fit(_subtractor(log), order)
    unmasked = _fit(_subtractor(log, **{"host-line-mask": False}), order)

    pd.testing.assert_frame_equal(masked, unmasked)


def test_the_wing_refit_inherits_the_host_line_mask(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log, object_wing_correction=True)
    refitInputs = []

    def record_refit(preFitOrder: pd.DataFrame, fitResult: Any, minLeakSigma: float, minSkySigma: float) -> Any:
        refitInputs.append(preFitOrder)
        return fitResult

    subtractor._refit_without_object_wings = record_refit

    _fit(subtractor, _order(hasHost=True))

    isInLine = (refitInputs[0]["wavelength"] - HOST_NM).abs() < 2 * HOST_SIGMA_NM
    assert refitInputs[0].loc[isInLine, "flagged_all_clipped"].all()


def test_a_catalogue_skyline_next_to_the_host_line_is_not_taken_for_the_host(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        toolkit, "get_skylines_dataframe", lambda *args, **kwargs: pd.DataFrame({"WAVELENGTH": [HOST_NM]})
    )
    subtractor = _subtractor(log)

    _fit(subtractor, _order(hasHost=True))

    assert subtractor.hostLineWindows == ()


def test_the_whole_catalogue_is_searched_for_skylines_not_only_the_strong_lines(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = []

    def record(*args: Any, **kwargs: Any) -> pd.DataFrame:
        calls.append(kwargs)
        return pd.DataFrame({"WAVELENGTH": [640.0]})

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", record)

    _fit(_subtractor(log), _order(hasHost=False))

    assert calls == [{"minBrightnessVIS": 0, "minBrightnessNIR": 0}]


def test_the_catalogue_is_not_read_when_the_mask_is_switched_off(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(*args: Any, **kwargs: Any) -> pd.DataFrame:
        raise AssertionError("the skyline catalogue was read")

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", fail)

    _fit(_subtractor(log, **{"host-line-mask": False}), _order(hasHost=False))


def test_the_merged_windows_of_the_frame_are_kept_for_the_qc(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)

    _fit(subtractor, _order(hasHost=True))

    assert len(subtractor.hostLineWindows) == 1
    assert subtractor.hostLineWindows[0].wavelengthMinNm < HOST_NM < subtractor.hostLineWindows[0].wavelengthMaxNm


def test_a_detection_that_fails_is_logged_and_the_sky_is_fitted_without_a_mask(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail(*args: Any, **kwargs: Any) -> pd.DataFrame:
        raise OSError("NO CATALOGUE")

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", fail)
    order = _order(hasHost=True)
    unmasked = _fit(_subtractor(log, **{"host-line-mask": False}), order)

    fitted = _fit(_subtractor(log), order)

    pd.testing.assert_frame_equal(fitted, unmasked)
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("NO CATALOGUE" in message for message in warnings)


def test_the_windows_are_logged_with_their_wavelengths(log: Any, no_catalogue_skylines: None) -> None:
    _fit(_subtractor(log), _order(hasHost=True))

    printed = [message for level, message in log.messages if level == "print"]
    assert any("HOST EMISSION LINE" in message and "619." in message for message in printed)


def test_no_line_is_logged_when_none_is_found(log: Any, no_catalogue_skylines: None) -> None:
    _fit(_subtractor(log), _order(hasHost=False))

    printed = [message for level, message in log.messages if level == "print"]
    assert not any("HOST EMISSION LINE" in message for message in printed)


# ----- THE QC AND THE SETTINGS -----


def test_the_qc_records_the_number_of_windows(log: Any) -> None:
    subtractor = _subtractor(log)
    subtractor.hostLineWindows = (host_feature(None, 619.5, 620.5, 9.0), host_feature(None, 650.0, 651.0, 9.0))

    subtractor._record_host_line_windows()

    row = subtractor.qc.iloc[0]
    assert (row["qc_name"], row["qc_value"]) == ("SKY HOST LINES MASKED", 2)
    assert bool(row["to_header"]) is True


def test_the_qc_table_is_untouched_when_no_window_was_found(log: Any) -> None:
    subtractor = _subtractor(log)
    subtractor.qc = pd.DataFrame({"qc_name": ["RON"], "qc_value": [3.3]})
    before = subtractor.qc.copy()

    subtractor._record_host_line_windows()

    pd.testing.assert_frame_equal(subtractor.qc, before)


def test_settings_without_the_host_line_keys_use_the_shipped_defaults(log: Any) -> None:
    subtractor = _subtractor(log)
    for key in ("host-line-mask", "host-line-sigma", "host-line-min-slit-contrast", "host-line-padding-pixels"):
        del subtractor.recipeSettings["sky-subtraction"][key]

    assert subtractor._host_line_settings() == (True, 5.5, 0.3, 3)


def test_the_host_line_keys_are_read(log: Any) -> None:
    subtractor = _subtractor(
        log, **{"host-line-sigma": 7, "host-line-min-slit-contrast": 0.4, "host-line-padding-pixels": 5}
    )

    assert subtractor._host_line_settings() == (True, 7.0, 0.4, 5)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("host-line-sigma", 0),
        ("host-line-sigma", float("nan")),
        ("host-line-min-slit-contrast", -0.1),
        ("host-line-min-slit-contrast", 1.5),
        ("host-line-padding-pixels", -1),
        ("host-line-padding-pixels", 2.5),
    ],
)
def test_an_invalid_host_line_setting_raises_naming_the_key(log: Any, key: str, value: float) -> None:
    subtractor = _subtractor(log, **{key: value})

    with pytest.raises(ValueError, match=key):
        subtractor._host_line_settings()


@pytest.mark.parametrize("key", ["host-line-sigma", "host-line-min-slit-contrast"])
@pytest.mark.parametrize("value", [True, False, np.True_, np.False_], ids=["True", "False", "np.True_", "np.False_"])
def test_a_boolean_is_not_a_number_for_a_host_line_limit(log: Any, key: str, value: Any) -> None:
    subtractor = _subtractor(log, **{key: value})

    with pytest.raises(ValueError, match=key):
        subtractor._host_line_settings()


def test_a_host_line_switch_that_is_not_true_or_false_raises(log: Any) -> None:
    subtractor = _subtractor(log, **{"host-line-mask": "yes"})

    with pytest.raises(ValueError, match="host-line-mask"):
        subtractor._host_line_settings()


def test_the_shipped_defaults_keep_the_mask_on_in_stare_and_stare_std_for_both_arms() -> None:
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)

    for recipe in ("soxs-stare", "soxs-stare-std"):
        for arm in ("vis", "nir"):
            sky = settings[recipe][arm]["sky-subtraction"]
            assert sky["host-line-mask"] is True
            assert math.isclose(sky["host-line-sigma"], 5.5)
            assert math.isclose(sky["host-line-min-slit-contrast"], 0.3)
            assert sky["host-line-padding-pixels"] == 3


class _TracebackLogger:
    """A logger that keeps the keyword arguments of each warning, which the shared recording logger drops."""

    def __init__(self) -> None:
        self.warnings: list[tuple[str, dict[str, Any]]] = []

    def warning(self, message: object, **kwargs: Any) -> None:
        self.warnings.append((str(message), kwargs))

    def __getattr__(self, name: str) -> Any:
        return lambda *args, **kwargs: None


def test_a_detection_that_fails_logs_a_warning_with_its_traceback(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(*args: Any, **kwargs: Any) -> pd.DataFrame:
        raise OSError("NO CATALOGUE")

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", fail)
    log = _TracebackLogger()

    _subtractor(log)._mask_host_lines([_order(hasHost=False)])

    assert [kwargs.get("exc_info") for message, kwargs in log.warnings if "NO CATALOGUE" in message] == [True]


# ----- NOTHING DETECTED LEAVES THE OUTPUT UNCHANGED -----

STAR_FLUX_PER_COLUMN = 3000.0
STAR_CENTRE = 0.4
STAR_SIGMA_ARCSEC = 0.5
UNCATALOGUED_SKY_LINES_NM = (610.0, 650.0)
UNCATALOGUED_SKY_LINE_PEAK = 500.0
UNCATALOGUED_SKY_LINE_SIGMA_NM = 0.05


def _realistic_order() -> pd.DataFrame:
    """A flat sky with a bright star and two sky lines that are not in the catalogue, and no host line."""
    order = _order(hasHost=False)
    star = STAR_FLUX_PER_COLUMN * np.exp(-0.5 * ((order["slit_position"] - STAR_CENTRE) / STAR_SIGMA_ARCSEC) ** 2)
    lines = sum(
        UNCATALOGUED_SKY_LINE_PEAK
        * np.exp(-0.5 * ((order["wavelength"] - centre) / UNCATALOGUED_SKY_LINE_SIGMA_NM) ** 2)
        for centre in UNCATALOGUED_SKY_LINES_NM
    )
    error = np.sqrt(order["error"] ** 2 + star.clip(lower=0) + lines)
    return order.assign(flux=order["flux"] + star + lines, error=error)


def _fit_and_record(subtractor: subtract_sky, order: pd.DataFrame) -> tuple[Any, ...]:
    """Clip the order, run ``_fit_sky_to_orders`` and ``_record_host_line_windows``; return everything they produce."""
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = subtractor.get_over_sampled_sky_from_order(
        order.copy(), clipBPs=True, clipSlitEdge=settings["clip-slit-edge-fraction"]
    )
    clipped = subtractor.clip_object_slit_positions([clipped], aggressive=settings["aggressive_object_masking"])
    fittedOrders, alltck, allKnots = subtractor._fit_sky_to_orders(clipped)
    subtractor._record_host_line_windows()
    return fittedOrders, alltck, allKnots, subtractor.qc


def test_when_nothing_is_detected_the_fit_and_the_qc_table_are_those_of_a_run_with_the_mask_off(
    log: Any, no_catalogue_skylines: None
) -> None:
    order = _realistic_order()
    masked = _subtractor(log)
    unmasked = _subtractor(log, **{"host-line-mask": False})

    maskedFrames, maskedTck, maskedKnots, maskedQc = _fit_and_record(masked, order)
    unmaskedFrames, unmaskedTck, unmaskedKnots, unmaskedQc = _fit_and_record(unmasked, order)

    assert masked.hostLineWindows == ()
    pd.testing.assert_frame_equal(maskedFrames[0], unmaskedFrames[0])
    for maskedValue, unmaskedValue in zip(maskedTck[0], unmaskedTck[0], strict=True):
        np.testing.assert_array_equal(maskedValue, unmaskedValue)
    np.testing.assert_array_equal(maskedKnots[0], unmaskedKnots[0])
    pd.testing.assert_frame_equal(maskedQc, unmaskedQc)
    assert maskedQc.empty


# ----- THE REDSHIFT FROM THE DETECTED WINDOWS AND THE COMPANION LINES -----

HOST_Z = 0.01
HALPHA_REST_NM = 656.28
NII_6583_REST_NM = 658.35
NII_6548_REST_NM = 654.80
SII_REST_NM = (671.64, 673.08)
ORDER_FIRST_NM = 640.0
ORDER_COLUMNS = 1250


def _hosted_order(strongRestNm: tuple[float, ...], seed: int = 2243, orderNumber: int = ORDER) -> pd.DataFrame:
    """One order from 640 to 690 nm on a flat sky, with a slit-filling host line at each strong rest wavelength."""
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(ORDER_COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slit = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)[slitPixel.ravel()]
    wavelength = ORDER_FIRST_NM + STEP_NM * column
    noiseless = np.full(column.shape, SKY_LEVEL)
    for restNm in strongRestNm:
        centre = restNm * (1 + HOST_Z) + HOST_ROTATION_NM_PER_ARCSEC * slit
        lineProfile = np.exp(-0.5 * ((wavelength - centre) / HOST_SIGMA_NM) ** 2)
        noiseless = noiseless + HOST_PEAK * lineProfile * np.exp(-0.5 * (slit / HOST_SLIT_SIGMA_ARCSEC) ** 2)
    error = np.sqrt(noiseless + READ_NOISE**2)
    frame = pd.DataFrame(
        {
            "order": float(orderNumber),
            "x": column,
            "y": slitPixel.ravel(),
            "wavelength": wavelength + 0.002 * slit,
            "slit_position": slit,
            "flux": noiseless + rng.standard_normal(noiseless.size) * error,
            "error": error,
            "mask": False,
        }
    )
    return frame.sort_values("wavelength", kind="stable").reset_index(drop=True)


TWO_STRONG_LINES = (HALPHA_REST_NM, NII_6583_REST_NM)


def test_two_detected_lines_give_the_host_redshift(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)

    _fit(subtractor, _hosted_order(TWO_STRONG_LINES))

    assert subtractor.hostRedshift is not None
    assert subtractor.hostRedshift.z == pytest.approx(HOST_Z, abs=2e-4)


def test_the_weak_companion_lines_are_left_out_of_the_fit_too(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)

    fitted = _fit(subtractor, _hosted_order(TWO_STRONG_LINES))

    for restNm in (NII_6548_REST_NM, *SII_REST_NM):
        isInLine = (fitted["wavelength"] - restNm * (1 + HOST_Z)).abs() < 0.1
        assert fitted.loc[isInLine, "flagged_all_clipped"].all(), restNm
    assert len(subtractor.hostInferredWindows) == 3


def test_a_companion_next_to_a_catalogue_skyline_is_not_added(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda *args, **kwargs: pd.DataFrame({"WAVELENGTH": [NII_6548_REST_NM * (1 + HOST_Z) + 0.2]}),
    )
    subtractor = _subtractor(log)

    _fit(subtractor, _hosted_order(TWO_STRONG_LINES))

    assert len(subtractor.hostInferredWindows) == 2


def test_one_detected_line_gives_no_redshift_and_no_companion(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)
    order = _hosted_order((HALPHA_REST_NM,))

    fitted = _fit(subtractor, order)
    unmasked = _fit(_subtractor(log, **{"host-line-mask": False}), order)

    assert subtractor.hostRedshift is None
    assert subtractor.hostInferredWindows == ()
    isInCompanion = (fitted["wavelength"] - NII_6548_REST_NM * (1 + HOST_Z)).abs() < 0.1
    np.testing.assert_array_equal(
        fitted.loc[isInCompanion, "flagged_all_clipped"], unmasked.loc[isInCompanion, "flagged_all_clipped"]
    )


def test_the_log_gives_the_redshift_and_the_added_windows(log: Any, no_catalogue_skylines: None) -> None:
    _fit(_subtractor(log), _hosted_order(TWO_STRONG_LINES))

    printed = [message for level, message in log.messages if level == "print"]
    assert any("REDSHIFT" in message and "0.010" in message and "ADDED" in message for message in printed)


def test_the_qc_records_the_redshift_and_the_number_of_added_windows(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)
    _fit(subtractor, _hosted_order(TWO_STRONG_LINES))

    subtractor._record_host_line_windows()

    rows = subtractor.qc.set_index("qc_name")["qc_value"]
    assert rows["SKY HOST LINES MASKED"] == 5
    assert rows["SKY HOST LINES INFERRED"] == 3
    assert rows["SKY HOST REDSHIFT"] == pytest.approx(HOST_Z, abs=2e-4)


def test_the_qc_has_no_redshift_rows_when_none_was_inferred(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)
    _fit(subtractor, _hosted_order((HALPHA_REST_NM,)))

    subtractor._record_host_line_windows()

    assert list(subtractor.qc["qc_name"]) == ["SKY HOST LINES MASKED"]


def test_with_the_mask_off_no_redshift_is_inferred(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log, **{"host-line-mask": False})

    _fit(subtractor, _hosted_order(TWO_STRONG_LINES))

    assert subtractor.hostRedshift is None
    assert subtractor.hostInferredWindows == ()


def _clipped_orders(subtractor: subtract_sky, orders: list[pd.DataFrame]) -> list[pd.DataFrame]:
    """The orders as they reach ``_fit_sky_to_orders``: clipped, with the object ranges flagged."""
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = [
        subtractor.get_over_sampled_sky_from_order(
            order.copy(), clipBPs=True, clipSlitEdge=settings["clip-slit-edge-fraction"]
        )
        for order in orders
    ]
    return subtractor.clip_object_slit_positions(clipped, aggressive=settings["aggressive_object_masking"])


def test_a_line_detected_in_the_neighbouring_order_is_masked_but_not_counted_as_inferred(
    log: Any, no_catalogue_skylines: None
) -> None:
    subtractor = _subtractor(log)
    withNii6548 = _hosted_order((NII_6548_REST_NM, *TWO_STRONG_LINES), orderNumber=1)
    withoutNii6548 = _hosted_order(TWO_STRONG_LINES, seed=2244, orderNumber=2)
    orders = _clipped_orders(subtractor, [withNii6548, withoutNii6548])

    subtractor._mask_host_lines(orders)

    isInNii6548 = (orders[1]["wavelength"] - NII_6548_REST_NM * (1 + HOST_Z)).abs() < 0.1
    assert orders[1].loc[isInNii6548, "flagged_all_clipped"].all()
    assert subtractor.hostRedshift is not None
    assert subtractor.hostRedshift.nWindows == 3
    assert len(subtractor.hostLineWindows) == 5
    assert len(subtractor.hostInferredWindows) == 2


def test_the_log_gives_the_number_of_windows_that_support_the_redshift(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)
    orders = _clipped_orders(
        subtractor,
        [
            _hosted_order((NII_6548_REST_NM, *TWO_STRONG_LINES), orderNumber=1),
            _hosted_order(TWO_STRONG_LINES, seed=7, orderNumber=2),
        ],
    )

    subtractor._mask_host_lines(orders)

    printed = [message for level, message in log.messages if level == "print"]
    assert any("FROM 3 DETECTED WINDOW(S)" in message for message in printed)


def test_the_settings_are_read_and_checked_once_for_the_search(log: Any, no_catalogue_skylines: None) -> None:
    subtractor = _subtractor(log)
    orders = _clipped_orders(subtractor, [_hosted_order(TWO_STRONG_LINES)])
    calls = []
    original = subtractor._host_line_settings
    subtractor._host_line_settings = lambda: calls.append(1) or original()

    subtractor._mask_host_lines(orders)

    assert len(calls) == 1


@pytest.mark.parametrize("failure", ["missing", "malformed", "invalid-yaml"])
def test_a_line_list_that_cannot_be_read_keeps_the_detected_windows_and_logs_the_traceback(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, no_catalogue_skylines: None, failure: str
) -> None:
    from soxspipe.commonutils import host_line_mask

    path = tmp_path / "host_emission_lines.yaml"
    if failure == "malformed":
        path.write_text("lines: 5\n")
    elif failure == "invalid-yaml":
        path.write_text("lines: [unclosed\n")
    monkeypatch.setattr(host_line_mask, "HOST_LINES_FILE", path)
    log = _TracebackLogger()
    subtractor = _subtractor(log)
    orders = _clipped_orders(subtractor, [_hosted_order(TWO_STRONG_LINES)])

    subtractor._mask_host_lines(orders)

    assert len(subtractor.hostLineWindows) == 2
    assert subtractor.hostRedshift is None
    assert subtractor.hostInferredWindows == ()
    assert [kwargs.get("exc_info") for message, kwargs in log.warnings if "line list" in message] == [True]
