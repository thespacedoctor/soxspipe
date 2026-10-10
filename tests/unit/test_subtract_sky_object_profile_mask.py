"""Object light must stay out of the sky model: the object mask comes from each order's spatial profile (DY-1279).

Each test drives a synthetic order through the three per-order steps that ``subtract_sky.subtract`` runs:
``get_over_sampled_sky_from_order``, ``clip_object_slit_positions`` and ``fit_bspline_curve_to_sky``.
The sky is flat and the object has Moffat wings along the slit.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml

import soxspipe
from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

ORDER = 12
COLUMNS = 800
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
SKY_LEVEL = 20.0
READ_NOISE = 3.3
OBJECT_CENTRE = 0.4
SEEING_FWHM = 1.0
# THE ATMOSPHERIC-TURBULENCE MOFFAT INDEX (TRUJILLO ET AL. 2001)
TURBULENCE_BETA = 4.765
# THE FRACTION OF THE OBJECT FLUX THAT MAY GO INTO THE SKY MODEL (DY-1279 ACCEPTANCE)
MAX_OBJECT_FLUX_IN_MODEL = 0.001

DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"


def _sky_settings() -> dict[str, Any]:
    """The shipped SOXS VIS stare-std sky-subtraction settings, with the QC plot turned off."""
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    skySettings = copy.deepcopy(settings["soxs-stare-std"]["vis"]["sky-subtraction"])
    skySettings["sky_model_qc_plot"] = False
    return skySettings


def _subtractor(log: Any) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.debug = False
    subtractor.binx = 1
    subtractor.biny = 1
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.recipeSettings = {"sky-subtraction": _sky_settings()}
    subtractor.bspline_order = subtractor.recipeSettings["sky-subtraction"]["bspline_order"]
    subtractor.qcPlotOrder = -1
    return subtractor


def _moffat_slit_profile(slitPosition: np.ndarray, beta: float, centre: float) -> np.ndarray:
    """A Moffat profile along the slit, normalised to sum to 1 over one wavelength column."""
    alpha = SEEING_FWHM / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))
    profile = (1.0 + ((slitPosition - centre) / alpha) ** 2) ** -beta
    return profile / profile.sum()


def _order(
    objectFluxPerColumn: float, beta: float = TURBULENCE_BETA, seed: int = 1279, centre: float = OBJECT_CENTRE
) -> pd.DataFrame:
    """One synthetic order on a flat sky, sorted by wavelength as the 2D-map dataframe is.

    The slit is tilted slightly against the wavelength axis, so neighbouring columns interleave in wavelength.
    The `object_flux` column holds each pixel's noiseless object flux.
    """
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slitPixel = slitPixel.ravel()
    slitPosition = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)[slitPixel]
    profile = _moffat_slit_profile(np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS), beta, centre)
    objectFlux = objectFluxPerColumn * profile[slitPixel]
    noiseless = SKY_LEVEL + objectFlux
    error = np.sqrt(noiseless + READ_NOISE**2)
    order = pd.DataFrame(
        {
            "order": ORDER,
            "x": column,
            "y": slitPixel,
            "wavelength": 400.0 + 0.01 * column + 0.002 * slitPosition,
            "slit_position": slitPosition,
            # ONE SET OF NORMAL DRAWS PER SEED, SO ORDERS WITH AND WITHOUT AN OBJECT SHARE THEIR NOISE PATTERN
            "flux": noiseless + rng.standard_normal(noiseless.size) * error,
            "error": error,
            "mask": False,
            "object_flux": objectFlux,
        }
    )
    return order.sort_values("wavelength", kind="stable").reset_index(drop=True)


def _model_order(subtractor: subtract_sky, order: pd.DataFrame) -> pd.DataFrame:
    """Run the per-order steps of `subtract` on one order and return the modelled order."""
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = subtractor.get_over_sampled_sky_from_order(
        order.copy(), clipBPs=True, clipSlitEdge=settings["clip-slit-edge-fraction"]
    )
    assert clipped is not None
    clipped = subtractor.clip_object_slit_positions([clipped], aggressive=settings["aggressive_object_masking"])[0]
    modelled, *_ = subtractor.fit_bspline_curve_to_sky(clipped)
    return modelled


def _object_flux_fraction_in_model(modelled: pd.DataFrame) -> float:
    """The object flux that goes into the sky model, as a fraction of the object flux.

    The sky fit uses only unclipped pixels, and the object spectrum is flat, so the object light in the model is the
    mean object flux of the unclipped pixels at every pixel of the order. This counts object light only. Comparing
    sky models made with and without the object would also count the shift in the upper-only clipping bias that the
    object's extra noise causes, which is a separate defect (DY-1281).
    """
    isFitted = ~modelled["flagged_all_clipped"]
    objectLightInModel = modelled.loc[isFitted, "object_flux"].mean() * len(modelled)
    return float(objectLightInModel / modelled["object_flux"].sum())


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_a_bright_star_with_moffat_wings_puts_under_a_tenth_of_a_percent_of_its_flux_into_the_sky_model(
    log: Any,
) -> None:
    """The wings of a bright star stay out of the sky model, so the star is not subtracted from itself."""
    subtractor = _subtractor(log)

    modelled = _model_order(subtractor, _order(objectFluxPerColumn=10000.0))

    assert _object_flux_fraction_in_model(modelled) < MAX_OBJECT_FLUX_IN_MODEL
    assert _warnings(log) == []


def test_an_empty_sky_loses_no_pixels_to_the_object_profile_mask(log: Any) -> None:
    """With no object on the slit, aggressive masking clips exactly the pixels the plain clipping does."""
    subtractor = _subtractor(log)
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = subtractor.get_over_sampled_sky_from_order(
        _order(objectFluxPerColumn=0.0), clipBPs=True, clipSlitEdge=settings["clip-slit-edge-fraction"]
    )
    plain = subtractor.clip_object_slit_positions([clipped.copy()], aggressive=False)[0]

    aggressive = subtractor.clip_object_slit_positions([clipped.copy()], aggressive=True)[0]

    pd.testing.assert_series_equal(aggressive["flagged_all_clipped"], plain["flagged_all_clipped"])
    assert _warnings(log) == []


def test_each_order_is_masked_where_its_own_object_lies(log: Any) -> None:
    """Two orders with the object at different slit positions each get a mask around their own object."""
    subtractor = _subtractor(log)
    settings = subtractor.recipeSettings["sky-subtraction"]
    lowerCentre, upperCentre = -2.0, 2.0
    clippedOrders = [
        subtractor.get_over_sampled_sky_from_order(
            _order(objectFluxPerColumn=10000.0, seed=seed, centre=centre),
            clipBPs=True,
            clipSlitEdge=settings["clip-slit-edge-fraction"],
        )
        for seed, centre in ((1279, lowerCentre), (1280, upperCentre))
    ]

    lowerOrder, upperOrder = subtractor.clip_object_slit_positions(clippedOrders, aggressive=True)

    def masked_slit_range(order: pd.DataFrame) -> tuple[float, float]:
        maskedFraction = (
            (order["flagged_object_clipped"] & ~order["flagged_edge_clipped"]).groupby(order["slit_position"]).mean()
        )
        maskedPositions = maskedFraction[maskedFraction > 0.5].index
        return float(maskedPositions.min()), float(maskedPositions.max())

    lowerLow, lowerHigh = masked_slit_range(lowerOrder)
    upperLow, upperHigh = masked_slit_range(upperOrder)
    # EACH MASK COVERS ITS OWN OBJECT'S WINGS, WHICH REACH ABOUT 2.8 ARCSEC FROM THE CENTRE AT 0.1 SIGMA
    assert lowerLow < lowerCentre - 2.0 and lowerHigh > lowerCentre + 2.0
    assert upperLow < upperCentre - 2.0 and upperHigh > upperCentre + 2.0
    # AND NO MASK REACHES THE OTHER ORDER'S OBJECT CENTRE
    assert lowerHigh < upperCentre
    assert upperLow > lowerCentre


def test_an_order_left_with_few_sky_pixels_per_wavelength_warns_and_names_the_order(log: Any) -> None:
    """Clipping 45% of the slit at each end leaves 4 pixels per wavelength, and the rolling clip takes more."""
    subtractor = _subtractor(log)
    subtractor.recipeSettings["sky-subtraction"]["clip-slit-edge-fraction"] = 0.45

    _model_order(subtractor, _order(objectFluxPerColumn=0.0))

    assert _warnings(log) == [
        f"ORDER {ORDER}: ONLY 2.0 SKY PIXELS PER WAVELENGTH ARE LEFT AFTER OBJECT MASKING (MINIMUM 5) "
        "- THE SKY MODEL MAY CONTAIN OBJECT LIGHT"
    ]


def test_a_workspace_settings_file_without_the_mask_threshold_still_masks_the_object_wings(log: Any) -> None:
    """Workspace settings are copied once at setup, so an older workspace lacks the setting; 0.1 sigma applies."""
    subtractor = _subtractor(log)
    del subtractor.recipeSettings["sky-subtraction"]["object_profile_mask_sigma"]

    modelled = _model_order(subtractor, _order(objectFluxPerColumn=10000.0))

    assert _object_flux_fraction_in_model(modelled) < MAX_OBJECT_FLUX_IN_MODEL


def test_an_object_above_the_threshold_across_the_whole_slit_leaves_the_faintest_slit_positions_as_sky(
    log: Any,
) -> None:
    """When the profile is above the threshold everywhere, the faintest slit positions stay in the sky fit.

    Masking the whole slit would leave the sky fit with no pixels at all.
    """
    subtractor = _subtractor(log)

    modelled = _model_order(subtractor, _order(objectFluxPerColumn=1000000.0, beta=1.0))

    isFitted = ~modelled["flagged_all_clipped"]
    assert isFitted.any()
    # THE FITTED PIXELS LIE AT THE SLIT END FARTHEST FROM THE OBJECT
    assert modelled.loc[isFitted, "slit_position"].min() < -3.5
    assert (
        f"ORDER {ORDER}: OBJECT LIGHT IS ABOVE THE MASK THRESHOLD ACROSS THE WHOLE SLIT - FITTING THE SKY TO THE "
        "5 FAINTEST SLIT BINS; THE SKY MODEL MAY CONTAIN OBJECT LIGHT"
    ) in _warnings(log)
