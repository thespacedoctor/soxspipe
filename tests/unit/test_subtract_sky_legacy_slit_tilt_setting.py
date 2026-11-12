"""A workspace ``soxspipe.yaml`` that still holds the removed ``slit_tilt_refit`` key changes nothing (DY-2083).

A workspace keeps its own copy of the settings, so copies made while the key existed still carry it. The key is not
read: the sky is fitted as if it were absent, even on an order whose map wavelengths are twisted along the slit, the
case the removed refit corrected.
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
from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.flat_slit_gradient import load_gradient_scales
from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

ORDER = 11
COLUMNS = 800
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
START_WAVELENGTH = 1700.0
DISPERSION = 0.37
LINE_SIGMA_PX = 1.06
LINE_SPACING_PX = 20
LINE_PEAK = 3000.0
CONTINUUM = 40.0
READ_NOISE = 8.6
# A TWIST OF +-0.04 PX PER ARCSEC AT THE ORDER ENDS, IN NM PER ARCSEC, WELL ABOVE THE NOISE
TWIST_C1 = -0.04 * DISPERSION

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
    subtractor.gradientScales = load_gradient_scales(
        detector_lookup(log=log, settings={"instrument": "soxs"}).get("NIR")
    )
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


def _twisted_order(seed: int = 2083) -> pd.DataFrame:
    """One sky-only order whose map wavelengths are off by the twist s * TWIST_C1 * u."""
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slitPixel = slitPixel.ravel()
    slitPosition = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)[slitPixel]
    pixel = column + 0.3 * slitPixel / SLIT_PIXELS
    wavelength = START_WAVELENGTH + DISPERSION * pixel
    u = 2 * (wavelength - wavelength.min()) / np.ptp(wavelength) - 1
    truePixel = pixel + slitPosition * TWIST_C1 * u / DISPERSION
    centres = np.arange(LINE_SPACING_PX / 2, COLUMNS, LINE_SPACING_PX)
    lines = LINE_PEAK * np.exp(-0.5 * ((truePixel[:, None] - centres[None, :]) / LINE_SIGMA_PX) ** 2).sum(axis=1)
    noiseless = CONTINUUM + lines
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


@pytest.mark.parametrize("legacyValue", [True, False])
def test_a_leftover_slit_tilt_refit_key_does_not_change_the_sky_fit(log: Any, legacyValue: bool) -> None:
    # ARRANGE
    order = _twisted_order()
    withoutKey = _modelled(_bare_subtractor(log, _sky_settings()), order)

    # ACT
    withKey = _modelled(_bare_subtractor(log, _sky_settings(slit_tilt_refit=legacyValue)), order)

    # ASSERT
    pd.testing.assert_frame_equal(withKey, withoutKey)
