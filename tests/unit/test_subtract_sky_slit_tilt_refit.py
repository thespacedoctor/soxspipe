"""The slit-tilt refit fits the sky again on wavelengths corrected for a map tilt along the slit (DY-2065).

The synthetic NIR order is sky only: Gaussian skylines on a faint continuum. Its map wavelengths are off by a twist
s * (c0 + c1 * u) like the one measured in NIR order 11. The tests drive the order through the per-order steps of
``subtract_sky.subtract``.
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

ORDER = 11
COLUMNS = 800
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
START_WAVELENGTH = 1700.0
# NM PER PIXEL, AS IN NIR ORDER 11
DISPERSION = 0.37
LINE_SIGMA_PX = 1.06
LINE_SPACING_PX = 20
LINE_PEAK = 3000.0
CONTINUUM = 40.0
READ_NOISE = 8.6
# A TWIST OF +-0.04 PX PER ARCSEC AT THE ORDER ENDS, IN NM PER ARCSEC: TEN TIMES THE O11 TWIST, SO THAT ONE SYNTHETIC
# ORDER SHOWS IT WELL ABOVE ITS NOISE
TWIST_C1 = -0.04 * DISPERSION
# A STANDARD STAR BRIGHT ENOUGH FOR THE OBJECT-WING CORRECTION TO REFIT ITS ORDER, WITH GAUSSIAN SEEING ALONG THE SLIT
STAR_FLUX_PER_COLUMN = 100000.0
STAR_SIGMA_ARCSEC = 0.6
# THE SHARE OF THE OFF CHI2 THAT THE REFIT MUST REMOVE ON THE TWISTED ORDER
MIN_CHI2_GAIN = 0.1

DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"


def _sky_settings(**overrides: Any) -> dict[str, Any]:
    """The shipped SOXS NIR stare sky-subtraction settings, with the QC plot turned off."""
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    skySettings = copy.deepcopy(settings["soxs-stare"]["nir"]["sky-subtraction"])
    skySettings["sky_model_qc_plot"] = False
    skySettings.update(overrides)
    return skySettings


def _bare_subtractor(log: Any, skySettings: dict[str, Any]) -> subtract_sky:
    """A subtract_sky with only the attributes the per-order steps read."""
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "NIR"
    subtractor.debug = False
    subtractor.binx = 1
    subtractor.biny = 1
    subtractor.axisA = "y"
    subtractor.axisB = "x"
    subtractor.ron = READ_NOISE
    subtractor.recipeSettings = {"sky-subtraction": skySettings}
    subtractor.bspline_order = skySettings["bspline_order"]
    subtractor.qcPlotOrder = -1
    return subtractor


def _order(c1: float, seed: int = 2065, starFluxPerColumn: float = 0.0) -> pd.DataFrame:
    """One order whose map wavelengths are off by the twist s * c1 * u, with an optional star at the slit centre."""
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slitPixel = slitPixel.ravel()
    slitPosition = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)[slitPixel]
    # NEIGHBOURING SLIT PIXELS SIT AT SLIGHTLY DIFFERENT MAP WAVELENGTHS, AS ON THE DETECTOR
    pixel = column + 0.3 * slitPixel / SLIT_PIXELS
    wavelength = START_WAVELENGTH + DISPERSION * pixel
    u = 2 * (wavelength - wavelength.min()) / np.ptp(wavelength) - 1
    truePixel = pixel + slitPosition * c1 * u / DISPERSION
    centres = np.arange(LINE_SPACING_PX / 2, COLUMNS, LINE_SPACING_PX)
    lines = LINE_PEAK * np.exp(-0.5 * ((truePixel[:, None] - centres[None, :]) / LINE_SIGMA_PX) ** 2).sum(axis=1)
    starProfile = np.exp(-0.5 * (slitPosition / STAR_SIGMA_ARCSEC) ** 2)
    star = starFluxPerColumn * starProfile / np.exp(-0.5 * (np.unique(slitPosition) / STAR_SIGMA_ARCSEC) ** 2).sum()
    noiseless = CONTINUUM + lines + star
    error = np.sqrt(noiseless + READ_NOISE**2)
    frame = pd.DataFrame(
        {
            "order": float(ORDER),
            "x": slitPixel,
            "y": column,
            "wavelength": wavelength,
            "slit_position": slitPosition,
            "flux": noiseless + rng.standard_normal(noiseless.size) * error,
            "error": error,
            "mask": False,
        }
    )
    return frame.sort_values("wavelength", kind="stable").reset_index(drop=True)


def _modelled(subtractor: subtract_sky, order: pd.DataFrame) -> pd.DataFrame:
    """Run the per-order steps of ``subtract`` on one order and return the modelled order."""
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = subtractor.get_over_sampled_sky_from_order(
        order.copy(), clipBPs=True, clipSlitEdge=settings["clip-slit-edge-fraction"]
    )
    clipped = subtractor.clip_object_slit_positions([clipped], aggressive=settings["aggressive_object_masking"])
    fittedOrders, _, _ = subtractor._fit_sky_to_orders(clipped)
    return fittedOrders[0].sort_index()


def _chi2(modelled: pd.DataFrame) -> float:
    """The chi2 of the sky-subtracted flux over every pixel of the order."""
    return float((((modelled["flux"] - modelled["sky_model"]) / modelled["error"]) ** 2).sum())


def test_the_refit_removes_a_twist_from_the_sky_residual(log: Any) -> None:
    # ARRANGE
    order = _order(TWIST_C1)
    off = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=False)), order)

    # ACT
    on = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=True)), order)

    # ASSERT
    assert _chi2(on) < (1 - MIN_CHI2_GAIN) * _chi2(off)


def test_the_refit_is_off_when_the_settings_do_not_carry_its_key(log: Any) -> None:
    # ARRANGE
    order = _order(TWIST_C1)
    skySettings = _sky_settings()
    skySettings.pop("slit_tilt_refit", None)
    explicitOff = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=False)), order)

    # ACT
    keyAbsent = _modelled(_bare_subtractor(log, skySettings), order)

    # ASSERT
    pd.testing.assert_frame_equal(keyAbsent, explicitOff)


def test_the_refit_never_raises_the_chi2_of_an_untwisted_order(log: Any) -> None:
    # ARRANGE
    order = _order(0.0)
    off = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=False)), order)

    # ACT
    on = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=True)), order)

    # ASSERT
    assert _chi2(on) <= _chi2(off)


def test_the_refit_returns_the_map_wavelengths(log: Any) -> None:
    # ARRANGE
    order = _order(TWIST_C1)

    # ACT
    on = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=True)), order)

    # ASSERT
    np.testing.assert_array_equal(on["wavelength"].to_numpy(), order.loc[on.index, "wavelength"].to_numpy())


def test_with_the_object_wing_correction_the_refit_never_raises_the_chi2_of_the_sky(log: Any) -> None:
    # ARRANGE
    order = _order(TWIST_C1, starFluxPerColumn=STAR_FLUX_PER_COLUMN)
    off = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=False, object_wing_min_leak_sigma=0.0)), order)
    isSky = ~off["flagged_object_clipped"]

    # ACT
    on = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=True, object_wing_min_leak_sigma=0.0)), order)

    # ASSERT
    assert "object_wing_model" in on
    assert _chi2(on[isSky]) <= _chi2(off[isSky])
