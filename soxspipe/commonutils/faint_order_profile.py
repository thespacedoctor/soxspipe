"""The profile steps of the faint-order trace: collapse the rows of a wavelength bin, then fit one Gaussian.

A faint VIS order is too noisy to fit slice by slice. Its rows are collapsed along an anchor track into one median
profile per wavelength bin, and a Gaussian of fixed width is fitted to that profile inside a narrow window. Both steps
are deterministic, need numpy only, and run no iterative optimiser.
"""

import warnings

import numpy as np

# THE BINNED GAUSSIAN FIT: THE SAMPLES USED REACH THIS MANY WIDTHS BEYOND THE SEARCH WINDOW, THE FEWEST SAMPLES THAT CAN
# BE FITTED, THE SPACING OF THE GRID OF TRIAL CENTRES (PIXELS) AND THE FITTED PARAMETERS (CENTRE, AMPLITUDE, BASELINE)
FIT_MARGIN_SIGMAS = 2.0
MIN_PROFILE_SAMPLES = 8
CENTRE_GRID_STEP = 0.05
FIT_PARAMETER_COUNT = 3
# THE RMS RESIDUAL NEVER FALLS BELOW THIS, SO A NOISE-FREE PROFILE HAS A FINITE SIGNAL-TO-NOISE RATIO
MIN_RMS = 1e-9


def collapse_bin_profile(image, trackA, trackB, halfLength, sliceAxis):
    """*collapse the rows along an anchor track into one profile, by their median*

    Each track point gives one cross-dispersion cut, ``2 * halfLength + 1`` pixels long and centred on the point,
    read from ``image`` by linear interpolation. The profile is the median of the cuts, so a cosmic-ray hit in one
    row does not move it. A pixel that is NaN or off the image does not count, and an offset with no value in any cut
    is NaN.

    **Key Arguments:**

    - ``image`` -- the 2D frame, NaN where masked
    - ``trackA`` -- the cross-dispersion position of each track point, in pixels
    - ``trackB`` -- the dispersion position of each track point, in pixels
    - ``halfLength`` -- the half-length of the cut, in pixels
    - ``sliceAxis`` -- the image axis the cut runs along, ``"x"`` (columns) or ``"y"`` (rows)

    **Return:**

    - ``offsets`` -- the pixel offsets from the track, ``-halfLength`` to ``+halfLength``
    - ``profile`` -- the median profile at those offsets
    """
    data = image if sliceAxis == "x" else image.T
    offsets = np.arange(-halfLength, halfLength + 1, dtype=float)
    rowIndex = np.rint(np.asarray(trackB, dtype=float)).astype(int)
    position = np.asarray(trackA, dtype=float)[:, np.newaxis] + offsets[np.newaxis, :]
    lower = np.floor(position).astype(int)
    fraction = position - lower
    isRowInside = ((rowIndex >= 0) & (rowIndex < data.shape[0]))[:, np.newaxis]
    isValid = isRowInside & (lower >= 0) & (lower + 1 < data.shape[1])
    rows = np.clip(rowIndex, 0, data.shape[0] - 1)[:, np.newaxis]
    left = np.clip(lower, 0, data.shape[1] - 1)
    right = np.clip(lower + 1, 0, data.shape[1] - 1)
    cuts = np.where(isValid, (1.0 - fraction) * data[rows, left] + fraction * data[rows, right], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return offsets, np.nanmedian(cuts, axis=0)


def fit_fixed_width_gaussian(profile, offsets, sigmaPx, windowPx, minSnr):
    """*find the centre of a Gaussian of fixed width in a profile, inside a window around the anchor*

    The centre runs over a fixed grid across the window. At each centre the amplitude and a flat baseline come from
    a linear least-squares solve, and the amplitude must not be negative. The centre with the smallest residual wins,
    so the answer is the same on every run. There is no iterative optimiser.

    The fit uses the samples within ``FIT_MARGIN_SIGMAS`` widths beyond the window. The peak signal-to-noise ratio is
    the amplitude over the rms residual. A fit that fails the limit, or whose best centre sits on the window edge
    (the sign of a source beyond the window), returns *None*.

    **Key Arguments:**

    - ``profile`` -- the collapsed profile, NaN where unknown
    - ``offsets`` -- the pixel offsets of the profile from the anchor track
    - ``sigmaPx`` -- the fixed Gaussian standard deviation, in pixels
    - ``windowPx`` -- the half-width of the search window around the anchor, in pixels
    - ``minSnr`` -- the smallest accepted peak signal-to-noise ratio

    **Return:**

    - ``centre`` -- the fitted centre, in pixels from the anchor track, or *None*
    """
    offsets = np.asarray(offsets, dtype=float)
    profile = np.asarray(profile, dtype=float)
    isUsed = np.isfinite(profile) & (np.abs(offsets) <= windowPx + FIT_MARGIN_SIGMAS * sigmaPx)
    if isUsed.sum() < MIN_PROFILE_SAMPLES:
        return None
    x, y = offsets[isUsed], profile[isUsed]

    grid = np.arange(-windowPx, windowPx + CENTRE_GRID_STEP / 2, CENTRE_GRID_STEP)
    shapes = np.exp(-0.5 * ((x[np.newaxis, :] - grid[:, np.newaxis]) / sigmaPx) ** 2)
    shapes = shapes - shapes.mean(axis=1, keepdims=True)
    centred = y - y.mean()
    cross = shapes @ centred
    power = (shapes**2).sum(axis=1)
    amplitude = cross / power
    gain = np.where(amplitude > 0, cross**2 / power, -np.inf)
    best = int(np.argmax(gain))
    if not np.isfinite(gain[best]):
        return None

    residualSumSquares = max(float(centred @ centred - gain[best]), 0.0)
    rms = max(np.sqrt(residualSumSquares / (len(x) - FIT_PARAMETER_COUNT)), MIN_RMS)
    isOnEdge = abs(grid[best]) > windowPx - CENTRE_GRID_STEP
    if amplitude[best] / rms < minSnr or isOnEdge:
        return None
    return float(grid[best])


def frame_to_nan_array(frame):
    """*the pixels of a frame as a float array, NaN where the frame is masked*

    **Key Arguments:**

    - ``frame`` -- a masked array or a ``CCDData`` frame, whose ``mask`` may be missing

    **Return:**

    - ``image`` -- a new 2D float array
    """
    data = np.array(getattr(frame, "data", frame), dtype=float)
    mask = getattr(frame, "mask", None)
    if mask is not None:
        data[np.broadcast_to(np.asarray(mask, dtype=bool), data.shape)] = np.nan
    return data
