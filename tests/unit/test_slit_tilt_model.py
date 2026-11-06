"""The slit-tilt model measures a wavelength-map tilt along the slit from an order's own skylines (DY-2065).

The synthetic order has Gaussian skylines on a flat continuum. The true wavelength of each pixel is the map wavelength
plus an injected tilt s * (c0 + c1 * u), where u is the map wavelength scaled to [-1, 1] over the order.
"""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.slit_tilt_model import corrected_wavelength, estimate_slit_tilt

pytestmark = pytest.mark.unit

COLUMNS = 1500
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
START_WAVELENGTH = 1700.0
DISPERSION = 0.07
LINE_SIGMA = 0.09
LINE_PEAK = 5000.0
CONTINUUM = 50.0
READ_NOISE = 8.0
LINE_SPACING = 4.0
# THE LARGEST MAP ERROR, IN WAVELENGTH UNITS, THAT AN UNTILTED ORDER MAY BE GIVEN: 0.005 PX, ABOUT 5 TIMES THE NOISE OF
# THE FIT AND 8 TIMES BELOW THE O11 TWIST
MAX_SPURIOUS_SHIFT = 0.005 * DISPERSION
# AN O11-SIZED TWIST: ABOUT +-0.04 PX PER 10 ARCSEC OF SLIT AT THE ORDER ENDS, PLUS A SMALL MEAN TILT
INJECTED_C0 = 0.0001
INJECTED_C1 = -0.0003


def _sky(wavelength: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The noiseless sky and its wavelength derivative."""
    centres = np.arange(START_WAVELENGTH + 1.0, START_WAVELENGTH + COLUMNS * DISPERSION - 1.0, LINE_SPACING)
    offsets = wavelength[:, None] - centres[None, :]
    gauss = LINE_PEAK * np.exp(-0.5 * (offsets / LINE_SIGMA) ** 2)
    return CONTINUUM + gauss.sum(axis=1), (-offsets / LINE_SIGMA**2 * gauss).sum(axis=1)


def _order(c0: float, c1: float, seed: int = 2065) -> dict[str, np.ndarray]:
    """One order whose map wavelengths are off by the injected tilt; the sky model is the untilted sky."""
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    slitPosition = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)[slitPixel.ravel()]
    # NEIGHBOURING SLIT PIXELS SIT AT SLIGHTLY DIFFERENT MAP WAVELENGTHS, AS ON THE DETECTOR
    wavelength = START_WAVELENGTH + DISPERSION * (column.ravel() + 0.3 * slitPixel.ravel() / SLIT_PIXELS)
    u = 2 * (wavelength - wavelength.min()) / (wavelength.max() - wavelength.min()) - 1
    trueWavelength = wavelength + slitPosition * (c0 + c1 * u)
    noiseless, _ = _sky(trueWavelength)
    error = np.sqrt(noiseless + READ_NOISE**2)
    skyModel, skyDerivative = _sky(wavelength)
    return {
        "wavelength": wavelength,
        "slitPosition": slitPosition,
        "flux": noiseless + rng.normal(0, error),
        "skyModel": skyModel,
        "skyDerivative": skyDerivative,
        "error": error,
        "isFitted": np.ones(wavelength.size, dtype=bool),
    }


def test_recovers_an_injected_tilt_that_changes_along_the_order():
    # ARRANGE
    order = _order(INJECTED_C0, INJECTED_C1)

    # ACT
    estimate = estimate_slit_tilt(**order)

    # ASSERT
    assert estimate.c0 == pytest.approx(INJECTED_C0, abs=0.1 * abs(INJECTED_C1))
    assert estimate.c1 == pytest.approx(INJECTED_C1, rel=0.1)


def test_finds_no_tilt_in_an_untilted_order():
    # ARRANGE
    order = _order(0.0, 0.0)

    # ACT
    estimate = estimate_slit_tilt(**order)

    # ASSERT
    assert abs(estimate.c0) * SLIT_HALF_LENGTH < MAX_SPURIOUS_SHIFT
    assert abs(estimate.c1) * SLIT_HALF_LENGTH < MAX_SPURIOUS_SHIFT


def test_corrected_wavelength_puts_each_pixel_at_its_true_wavelength():
    # ARRANGE
    order = _order(INJECTED_C0, INJECTED_C1)
    u = 2 * (order["wavelength"] - order["wavelength"].min()) / np.ptp(order["wavelength"]) - 1
    trueWavelength = order["wavelength"] + order["slitPosition"] * (INJECTED_C0 + INJECTED_C1 * u)
    estimate = estimate_slit_tilt(**order)

    # ACT
    corrected = corrected_wavelength(estimate, order["wavelength"], order["slitPosition"])

    # ASSERT
    injectedError = np.abs(order["wavelength"] - trueWavelength).max()
    assert np.abs(corrected - trueWavelength).max() < 0.1 * injectedError


def test_corrected_wavelength_leaves_an_untilted_order_at_its_map_wavelengths():
    # ARRANGE
    order = _order(0.0, 0.0)
    estimate = estimate_slit_tilt(**order)

    # ACT
    corrected = corrected_wavelength(estimate, order["wavelength"], order["slitPosition"])

    # ASSERT
    assert np.abs(corrected - order["wavelength"]).max() < MAX_SPURIOUS_SHIFT


@pytest.mark.parametrize(
    "spoil",
    ["no-fitted-pixels", "one-wavelength", "nan-slit-positions", "flat-sky"],
)
def test_an_order_that_cannot_constrain_the_tilt_raises_a_value_error(spoil: str):
    # ARRANGE
    order = _order(INJECTED_C0, INJECTED_C1)
    if spoil == "no-fitted-pixels":
        order["isFitted"] = np.zeros_like(order["isFitted"])
    elif spoil == "one-wavelength":
        order["wavelength"] = np.full_like(order["wavelength"], START_WAVELENGTH)
    elif spoil == "nan-slit-positions":
        order["slitPosition"] = np.full_like(order["slitPosition"], np.nan)
    else:
        order["skyDerivative"] = np.zeros_like(order["skyDerivative"])

    # ACT / ASSERT
    with pytest.raises(ValueError):
        estimate_slit_tilt(**order)
