#!/usr/bin/env python
"""
*the combined photon-noise signal-to-noise of the catalogue sky lines in one extracted order*

Each catalogue line's signal is the highest ``skyFlux`` sample within ``SKY_SNR_PEAK_HALF_WINDOW`` samples of its
arc-solution position, minus the median of the surrounding ``SKY_SNR_BASELINE_HALF_WINDOW`` samples. Its noise is
the photon noise of that baseline plus the read noise, ``sqrt(baseline + ron^2)``, in electrons per pixel. The lines
add in quadrature, as the precision of a shift measured from all of them does.

A frame is the mean of NDIT integrations, so its variance is that of one integration divided by NDIT and every
signal-to-noise is the NDIT 1 value times ``sqrt(NDIT)`` (DY-2323).

``skyFlux`` is a smooth sky model, so the scatter of its own samples says nothing about its noise (DY-1663 probe:
the line-free residual is about 0.01 e-). The photon noise of the pixels it was fitted to is used instead.

Author
: David Young

Date Created
: October 20, 2026
"""

import numpy as np

# HALF-WIDTH IN SAMPLES OF THE WINDOW SEARCHED FOR A LINE'S PEAK; COVERS SOXS VIS NIGHT-TO-ARC FLEXURE (<= 3.3 px)
SKY_SNR_PEAK_HALF_WINDOW = 4
# HALF-WIDTH IN SAMPLES OF THE WINDOW WHOSE MEDIAN IS A LINE'S BASELINE
SKY_SNR_BASELINE_HALF_WINDOW = 25


def sky_line_snrs(wavelength, skyFlux, lineWavelengths, readNoise, ndit=1):
    """*return the photon-noise signal-to-noise of each catalogue line in the order, brightest first*

    **Key Arguments:**

    - ``wavelength`` -- the order's arc-solution wavelengths (nm), in any order
    - ``skyFlux`` -- the sky flux at each wavelength (e- per pixel); NaN samples are ignored
    - ``lineWavelengths`` -- the catalogue sky-line wavelengths (nm); lines outside the order are ignored
    - ``readNoise`` -- the detector read noise (e-)
    - ``ndit`` -- the number of integrations the frame averages, a whole number of 1 or more. Default *1*

    **Return:**

    - ``snrs`` -- one value per catalogue line in the order, sorted descending; a line below its baseline has 0
    """
    if ndit < 1:
        raise ValueError(f"NDIT must be 1 or more, not {ndit!r}")
    wavelength = np.asarray(wavelength, dtype=float)
    skyFlux = np.asarray(skyFlux, dtype=float)
    finite = np.isfinite(wavelength) & np.isfinite(skyFlux)
    ascending = np.argsort(wavelength[finite], kind="stable")
    wavelength, skyFlux = wavelength[finite][ascending], skyFlux[finite][ascending]
    if wavelength.size == 0:
        return np.array([])

    lines = np.asarray(lineWavelengths, dtype=float)
    lines = lines[(lines >= wavelength[0]) & (lines <= wavelength[-1])]
    snrs = []
    for line in lines:
        centre = int(np.argmin(np.abs(wavelength - line)))
        peak = np.max(skyFlux[max(centre - SKY_SNR_PEAK_HALF_WINDOW, 0) : centre + SKY_SNR_PEAK_HALF_WINDOW + 1])
        baseline = np.median(
            skyFlux[max(centre - SKY_SNR_BASELINE_HALF_WINDOW, 0) : centre + SKY_SNR_BASELINE_HALF_WINDOW + 1]
        )
        snrs.append(max((peak - baseline) / np.sqrt(max(baseline, 0.0) + readNoise**2), 0.0))
    snrs = np.sort(np.asarray(snrs, dtype=float))[::-1]
    # THE FRAME IS THE MEAN OF NDIT INTEGRATIONS; NDIT 1 IS RETURNED UNTOUCHED
    return snrs if ndit == 1 else snrs * np.sqrt(ndit)


def combine_line_snrs(lineSnrs):
    """*return the combined signal-to-noise of an order's lines: their quadrature sum, 0 for no lines*

    **Key Arguments:**

    - ``lineSnrs`` -- the per-line signal-to-noise values (``sky_line_snrs``)
    """
    return float(np.sqrt(np.sum(np.asarray(lineSnrs, dtype=float) ** 2)))
