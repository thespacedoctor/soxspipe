"""Measure the read noise of a UVB/VIS science frame from the prescan of the raw frame (DY-2437).

The prescan is a strip of columns the detector reads without exposing it, so its scatter is the read noise of that
frame alone. The night frames of the SOXS VIS arm have more read noise than the morning master bias, which sets the
read noise of the reduction otherwise. The prescan is measured on the raw frame in ADU, before the frame is trimmed.
"""

# THE RECIPES WHOSE FRAMES ARE SCIENCE FRAMES. THE NAME MAY END IN ``-std``
SCIENCE_RECIPE_PREFIXES = ("soxs-stare", "soxs-nod", "soxs-offset")
# THE ARMS WITH A PRESCAN. THE NIR ARM HAS NONE AND TAKES ITS READ NOISE FROM ITS SAMPLE COUNT
PRESCAN_ARMS = ("UVB", "VIS")
# THE DETECTOR-PARAMETER KEY OF THE PRESCAN WINDOW, IN THE SAME FORM AS ``science-pixels``
PRESCAN_KEY = "prescan-pixels"
# THE FEWEST FINITE PIXELS THAT GIVE A READ NOISE TO TRUST
MIN_PRESCAN_PIXELS = 100
# THE CLIP OF THE PRESCAN PIXELS BEFORE THEIR STANDARD DEVIATION IS TAKEN
CLIP_SIGMA = 4
CLIP_ITERATIONS = 5


def uses_prescan_read_noise(arm, recipeName, detectorParams):
    """*say whether a recipe takes the read noise of its frames from their prescan*

    **Key Arguments:**

    - ``arm`` -- the arm of the frames
    - ``recipeName`` -- the name of the recipe, for example ``soxs-stare-std``
    - ``detectorParams`` -- the detector parameters of the arm

    **Return:**

    - ``isUsed`` -- True for a UVB or VIS science recipe on a detector that has a prescan window

    **Usage:**

    ```python
    from soxspipe.commonutils.prescan_read_noise import uses_prescan_read_noise
    isUsed = uses_prescan_read_noise("VIS", "soxs-stare", detectorParams)
    ```
    """
    return (
        str(arm).upper() in PRESCAN_ARMS
        and str(recipeName).startswith(SCIENCE_RECIPE_PREFIXES)
        and PRESCAN_KEY in detectorParams
    )


def binned_window(window, binning):
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


def measure_prescan_read_noise(data, prescanPixels, binning, gain):
    """*measure the read noise (e-) of a raw frame from its prescan*

    The noise is the standard deviation of the prescan pixels that survive a 4 sigma clip about the median, so cosmic
    rays in the prescan do not bias it. The clip level comes from the median absolute deviation. The deviation itself
    is not used as the noise: the raw counts are whole numbers, so for a scatter of a few ADU it moves in steps of
    15 to 20 % of the read noise.

    **Key Arguments:**

    - ``data`` -- the raw frame (ADU), before it is trimmed
    - ``prescanPixels`` -- the prescan window on an unbinned frame, as ``rows`` and ``columns`` limits
    - ``binning`` -- the binning of the frame as ``[rows, columns]``
    - ``gain`` -- the gain in e-/ADU, a quantity or a number

    **Return:**

    - ``ron`` -- the read noise (e-)

    A ``ValueError`` says why the noise cannot be measured: the window is outside the frame, it holds too few
    finite pixels, or its scatter is not a positive number.

    **Usage:**

    ```python
    from soxspipe.commonutils.prescan_read_noise import measure_prescan_read_noise
    ron = measure_prescan_read_noise(rawData, detectorParams["prescan-pixels"], [1, 1], detectorParams["gain"])
    ```
    """
    import math

    import numpy as np
    from astropy.stats import sigma_clip

    rs, re, cs, ce = binned_window(prescanPixels, binning)
    if min(rs, cs) < 0 or re > data.shape[0] or ce > data.shape[1] or re <= rs or ce <= cs:
        raise ValueError(f"the prescan window (rows {rs}:{re}, columns {cs}:{ce}) is outside the frame {data.shape}")

    pixels = np.asarray(data[rs:re, cs:ce], dtype=np.float64)
    pixels = pixels[np.isfinite(pixels)]
    if pixels.size < MIN_PRESCAN_PIXELS:
        raise ValueError(f"the prescan holds too few finite pixels ({pixels.size}; at least {MIN_PRESCAN_PIXELS})")

    survivors = sigma_clip(pixels, sigma=CLIP_SIGMA, maxiters=CLIP_ITERATIONS, cenfunc="median", stdfunc="mad_std")
    ron = float(np.std(survivors.compressed())) * float(getattr(gain, "value", gain))
    if not math.isfinite(ron) or ron <= 0:
        raise ValueError(f"the prescan scatter ({ron} e-) is not a positive number")
    return ron
