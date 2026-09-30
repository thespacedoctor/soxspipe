"""Analytic contracts for one-dimensional spectral extraction."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import (
    _sigma_clip_and_mask,
    compute_extractions,
)

pytestmark = pytest.mark.unit


def test_compute_extractions_recovers_a_known_profile_and_sorts_wavelengths() -> None:
    wavelengths = np.array([503.0, 501.0, -1.0, 502.0])
    profile = np.tile(np.array([[0.25], [0.5], [0.25]]), (1, wavelengths.size))
    amplitudes = np.array([40.0, 80.0, 20.0, 60.0])
    mask = np.zeros_like(profile, dtype=bool)
    mask[0, 1] = True
    orderImages = {
        "fluxRaw": profile * amplitudes,
        "fluxSky": np.tile(np.array([1.0, 2.0, 3.0, 4.0]), (3, 1)),
        "objectProfile": profile,
        "variance": np.full_like(profile, 4.0),
        "wavelength": np.tile(wavelengths, (3, 1)),
        "mask": mask,
    }
    slices = pd.DataFrame(
        {
            "order": np.full(wavelengths.size, 12),
            "pixelScaleNm": np.full(wavelengths.size, 0.1),
        }
    )

    result = compute_extractions(slices, orderImages, order=12)

    np.testing.assert_array_equal(result["wavelengthMean"], [501.0, 502.0, 503.0])
    np.testing.assert_allclose(
        result["extractedFluxOptimal"], [80.0, 60.0, 40.0], rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(
        result["extractedFluxBoxcar"], [80.0, 60.0, 40.0], rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(
        result["extractedFluxBoxcarRobust"],
        [60.0, 60.0, 40.0],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        result["skyFlux"], [2.0, 4.0, 1.0], rtol=1e-12, atol=1e-12
    )
    expectedVariance = np.array([12.8, 32.0 / 3.0, 32.0 / 3.0])
    np.testing.assert_allclose(
        result["varianceSpectrum"], expectedVariance, rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(
        result["snr"],
        result["extractedFluxOptimal"] / np.sqrt(expectedVariance),
        rtol=1e-12,
        atol=1e-12,
    )


def test_sigma_clip_mask_keeps_bad_pixels_and_nan_values_masked() -> None:
    flux = np.ones((3, 64), dtype=float)
    flux[0, 0] = np.nan
    flux[1, -1] = 10_000.0
    badPixels = np.zeros_like(flux, dtype=np.uint8)
    badPixels[2, 3] = 2

    result = _sigma_clip_and_mask(flux, badPixels)

    assert result[0, 0]
    assert result[1, -1]
    assert result[2, 3]
    assert badPixels[2, 3] == 1
