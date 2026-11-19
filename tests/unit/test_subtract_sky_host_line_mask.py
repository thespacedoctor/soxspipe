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

    assert subtractor._host_line_settings() == (True, 5.0, 0.3, 3)


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
            assert math.isclose(sky["host-line-sigma"], 5.0)
            assert math.isclose(sky["host-line-min-slit-contrast"], 0.3)
            assert sky["host-line-padding-pixels"] == 3
