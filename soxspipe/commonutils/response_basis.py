#!/usr/bin/env python
"""
*The counts basis shared by the response function and the flux calibration*

The response function and the flux calibration both work on the flat-fielded ``FLUX_COUNTS`` of an extraction,
divided by the detector binning along the dispersion axis: electrons per unbinned pixel. The master flat is scaled by
one global level, so it keeps the pixel scale (Angstrom per pixel) and the blaze of each order, and dividing by it
already removes the pixel scale. Only the binning is left: a pixel binned twice along the dispersion axis collects
twice the counts, and the flat-fielding does not remove that.

:Author:
    David Young

:Date Created:
    October 22, 2026
"""

from typing import Any

# THE PRIMARY-HEADER CARD OF A RESP PRODUCT THAT STATES THE BASIS OF ITS RESPONSE CURVE
RESPONSE_BASIS_KEYWORD = "ESO PRO RESP BASIS"
RESPONSE_BASIS = "PER_UNBINNED_PIXEL"
RESPONSE_BASIS_COMMENT = "flat-field e-/s/unbinned px"

# ARMS WHOSE DETECTOR IS NEVER BINNED
UNBINNED_ARMS = ("NIR",)


def dispersion_binning(header: Any, arm: str, kw: Any) -> int:
    """*return the detector binning along the dispersion axis of a frame*

    On the SOXS VIS detector the wavelength changes along the y-axis (about 0.03 nm per pixel, against about 0.005
    nm per pixel along x), so this is ``ESO DET BINY``. The NIR detector is never binned.

    **Key Arguments:**

    - ``header`` -- the FITS header of the frame, or of a product that keeps the raw-frame cards
    - ``arm`` -- the spectrograph arm
    - ``kw`` -- the keyword lookup function

    **Return:**

    - ``binning`` -- the binning along the dispersion axis, 1 or more
    """
    import math

    if arm.upper() in UNBINNED_ARMS:
        return 1
    keyword = kw("WIN_BINY")
    if header is None or keyword not in header:
        raise ValueError(
            f"The {arm} header has no '{keyword}' card, so the binning along the dispersion axis is unknown and the "
            "spectrum cannot be flux calibrated."
        )
    value = header[keyword]
    try:
        binning = float(value)
    except (TypeError, ValueError) as e:
        raise ValueError(f"The {arm} header card '{keyword}' is {value!r}, which is not a binning.") from e
    if not math.isfinite(binning) or not binning.is_integer() or binning < 1:
        raise ValueError(
            f"The {arm} header card '{keyword}' is {value!r}; a binning must be a whole number of 1 or more."
        )
    return int(binning)


def counts_per_unbinned_pixel(spectrum: Any, dispersionBinning: int) -> Any:
    """*return the flat-fielded electrons per unbinned pixel of an extracted spectrum*

    **Key Arguments:**

    - ``spectrum`` -- a flat-fielded extracted spectrum (data frame or table) with a ``FLUX_COUNTS`` column
    - ``dispersionBinning`` -- the detector binning along the dispersion axis (see ``dispersion_binning``)

    **Return:**

    - ``counts`` -- the electrons per unbinned pixel at each sample, as a float array
    """
    import numpy as np

    return np.asarray(spectrum["FLUX_COUNTS"], dtype=float) / dispersionBinning


def require_response_basis(header: Any, responsePath: str) -> None:
    """*raise a ValueError unless a RESP product header states the per-unbinned-pixel basis*

    A RESP product written before the basis card existed scales the counts of the standard's own binning. It is
    refused, not converted.

    **Key Arguments:**

    - ``header`` -- the primary header of the RESP product
    - ``responsePath`` -- the path of the RESP product, for the error message
    """
    basis = header.get(RESPONSE_BASIS_KEYWORD)
    if basis != RESPONSE_BASIS:
        raise ValueError(
            f"The response product {responsePath} has no '{RESPONSE_BASIS_KEYWORD} = {RESPONSE_BASIS}' card "
            f"(found {basis!r}). It was made by an older soxspipe release that did not correct for the detector "
            "binning. Please rebuild the response from the standard star with this release, then flux calibrate again."
        )
