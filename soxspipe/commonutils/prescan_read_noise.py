"""Measure the read noise of a UVB/VIS science frame from the prescan of the raw frame (DY-2437).

The prescan is a strip of columns the detector reads without exposing it, so its scatter is the read noise of that
frame alone. The night frames of the SOXS VIS arm have more read noise than the morning master bias, which sets the
read noise of the reduction otherwise. The prescan is measured on the raw frame in ADU, before the frame is trimmed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

# THE RECIPES WHOSE FRAMES ARE SCIENCE FRAMES. THE NAME MAY END IN ``-std``
SCIENCE_RECIPE_PREFIXES = ("soxs-stare", "soxs-nod", "soxs-offset")
# THE FIRST WORD OF THE DPR TYPE OF A RAW SCIENCE FRAME: ``OBJECT,ASYNC``, ``SKY,...``, ``STD,FLUX``, ``STD,TELLURIC``.
# A BIAS, DARK OR LAMP FRAME IN THE SOF IS A CALIBRATION FRAME
SCIENCE_DPR_TYPES = ("OBJECT", "SKY", "STD")
# THE ARMS WITH A PRESCAN. THE NIR ARM HAS NONE AND TAKES ITS READ NOISE FROM ITS SAMPLE COUNT. UVB IS KEPT FOR
# X-SHOOTER, WHICH HAS NO PRESCAN WINDOW YET
PRESCAN_ARMS = ("UVB", "VIS")
# THE DETECTOR-PARAMETER KEY OF THE PRESCAN WINDOW, IN THE SAME FORM AS ``science-pixels``
PRESCAN_KEY = "prescan-pixels"
# THE FEWEST FINITE PIXELS THAT GIVE A READ NOISE TO TRUST
MIN_PRESCAN_PIXELS = 100
# THE CLIP OF THE PRESCAN PIXELS BEFORE THEIR STANDARD DEVIATION IS TAKEN
CLIP_SIGMA = 4
CLIP_ITERATIONS = 5


def is_science_dpr_type(dprType: str | None) -> bool:
    """*say whether a ``DPR TYPE`` is that of a raw science frame*

    **Key Arguments:**

    - ``dprType`` -- the ``ESO DPR TYPE`` of the frame, for example ``OBJECT,ASYNC``, or *None* when it has none

    **Return:**

    - ``isScience`` -- True when its first word is OBJECT, SKY or STD

    **Usage:**

    ```python
    from soxspipe.commonutils.prescan_read_noise import is_science_dpr_type
    isScience = is_science_dpr_type("STD,FLUX")
    ```
    """
    if not isinstance(dprType, str):
        return False
    return dprType.split(",")[0].strip() in SCIENCE_DPR_TYPES


def uses_prescan_read_noise(
    arm: str, recipeName: str, detectorParams: Mapping[str, Any], dprType: str | None = None
) -> bool:
    """*say whether a frame takes its read noise from its prescan*

    **Key Arguments:**

    - ``arm`` -- the arm of the frames
    - ``recipeName`` -- the name of the recipe, for example ``soxs-stare-std``
    - ``detectorParams`` -- the detector parameters of the arm
    - ``dprType`` -- the ``ESO DPR TYPE`` of the frame. A bias, dark or lamp frame in the SOF of a science recipe is
      not a science frame. Default *None*, which is not a science frame either

    **Return:**

    - ``isUsed`` -- True for a UVB or VIS science frame of a science recipe, on a detector that has a prescan window

    **Usage:**

    ```python
    from soxspipe.commonutils.prescan_read_noise import uses_prescan_read_noise
    isUsed = uses_prescan_read_noise("VIS", "soxs-stare", detectorParams, "OBJECT")
    ```
    """
    return (
        str(arm).upper() in PRESCAN_ARMS
        and str(recipeName).startswith(SCIENCE_RECIPE_PREFIXES)
        and PRESCAN_KEY in detectorParams
        and is_science_dpr_type(dprType)
    )


def binned_window(window: Mapping[str, Mapping[str, int]], binning: Sequence[int]) -> tuple[int, int, int, int]:
    """*the row and column limits of a detector window, for the binning of the frame*

    **Key Arguments:**

    - ``window`` -- the ``rows`` and ``columns`` limits (``start`` and ``end``) of the window on an unbinned frame
    - ``binning`` -- the binning of the frame as ``[rows, columns]``

    **Return:**

    - ``limits`` -- the row start, row end, column start and column end

    **Usage:**

    ```python
    from soxspipe.commonutils.prescan_read_noise import binned_window
    rs, re, cs, ce = binned_window(detectorParams["science-pixels"], [1, 2])
    ```
    """
    rs, re = window["rows"]["start"], window["rows"]["end"]
    cs, ce = window["columns"]["start"], window["columns"]["end"]
    if binning[0] > 1:
        rs, re = int(rs / binning[0]), int(re / binning[0])
    if binning[1] > 1:
        cs, ce = int(cs / binning[1]), int(ce / binning[1])
    return rs, re, cs, ce


def _gain_in_electrons_per_adu(gain: Any) -> float:
    """*the gain as a number in e-/ADU: a quantity is converted from its unit, a plain number is taken as e-/ADU*"""
    from astropy import units as u

    if isinstance(gain, u.Quantity):
        try:
            return float(gain.to_value(u.electron / u.adu))
        except u.UnitConversionError as e:
            raise ValueError(f"the gain ({gain}) is not in electrons per ADU") from e
    return float(gain)


def measure_prescan_read_noise(
    data: Any, prescanPixels: Mapping[str, Mapping[str, int]], binning: Sequence[int], gain: Any
) -> float:
    """*measure the read noise (e-) of a raw frame from its prescan*

    The noise is the standard deviation of the prescan pixels that survive a 4 sigma clip about the median, so cosmic
    rays in the prescan do not bias it. The clip level comes from the median absolute deviation. The deviation itself
    is not used as the noise: the raw counts are whole numbers, so for a scatter of a few ADU it moves in steps of
    15 to 20 % of the read noise.

    **Key Arguments:**

    - ``data`` -- the raw frame (ADU), before it is trimmed
    - ``prescanPixels`` -- the prescan window on an unbinned frame, as ``rows`` and ``columns`` limits
    - ``binning`` -- the binning of the frame as ``[rows, columns]``
    - ``gain`` -- the gain, a quantity in any unit that converts to e-/ADU, or a number in e-/ADU

    **Return:**

    - ``ron`` -- the read noise (e-)

    A ``ValueError`` says why the noise cannot be measured: the gain is not in electrons per ADU, the window is
    outside the frame, it holds too few finite pixels, or its scatter is not a positive number.

    **Usage:**

    ```python
    from soxspipe.commonutils.prescan_read_noise import measure_prescan_read_noise
    ron = measure_prescan_read_noise(rawData, detectorParams["prescan-pixels"], [1, 1], detectorParams["gain"])
    ```
    """
    import math

    import numpy as np
    from astropy.stats import sigma_clip

    gainValue = _gain_in_electrons_per_adu(gain)
    rs, re, cs, ce = binned_window(prescanPixels, binning)
    if min(rs, cs) < 0 or re > data.shape[0] or ce > data.shape[1] or re <= rs or ce <= cs:
        raise ValueError(f"the prescan window (rows {rs}:{re}, columns {cs}:{ce}) is outside the frame {data.shape}")

    pixels = np.asarray(data[rs:re, cs:ce], dtype=np.float64)
    pixels = pixels[np.isfinite(pixels)]
    if pixels.size < MIN_PRESCAN_PIXELS:
        raise ValueError(f"the prescan holds too few finite pixels ({pixels.size}; at least {MIN_PRESCAN_PIXELS})")

    survivors = sigma_clip(pixels, sigma=CLIP_SIGMA, maxiters=CLIP_ITERATIONS, cenfunc="median", stdfunc="mad_std")
    ron = float(np.std(survivors.compressed())) * gainValue
    if not math.isfinite(ron) or ron <= 0:
        raise ValueError(f"the prescan scatter ({ron} e-) is not a positive number")
    return ron
