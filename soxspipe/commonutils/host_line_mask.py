"""Host emission lines that fill the slit of a stare frame, found so the sky fit can leave them out (DY-2242).

The stare sky model is one b-spline in wavelength, fitted to every slit pixel that is not object. A host emission line
that fills the slit is in those pixels, so the model takes the line from the object. The line has to be found without a
redshift. Three things tell it from the sky and from a point source:

- it is present along most of the slit, which a cosmic ray, a hot pixel and the wings of a star are not
- it changes along the slit, in strength and in wavelength (rotation), which a sky line does not
- it is narrow and rises above the continuum, which telluric absorption and the star's smooth wings do not

For each order the pixels are binned by true wavelength (2 detector lines) and by slit position (1 arcsec), and each
cell is the median flux. The excess of a cell is the cell less its running median along wavelength. The median of the
excess over the slit is the common excess. The largest value over the slit of the excess less the common excess is the
slit variation. A run of bins is a feature when both are significant and, summed over the run, the slit variation is
a large enough share of the peak excess (the slit contrast). The share is taken over the run, not each bin: a bright
unresolved sky line is sampled differently in each slit bin because the slit is tilted, and its flanks then vary along
the slit, but its core, which holds most of the light, does not. Bins near a catalogue skyline, and the ends of the
order, are not searched.

The window of a feature runs over the wavelengths of its bins, which already include the rotation of the host, plus a
padding of detector lines. ``subtract_sky`` leaves the pixels inside the windows out of the sky fit, so the spline
interpolates the sky across them. Nothing here changes any data.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from soxspipe.commonutils.nod_host_contamination import (
    MAD_TO_SIGMA,
    MEDIAN_VARIANCE_FACTOR,
    SKY_MASK_HALF_WIDTH_LINES,
    host_feature,
)

# THE DEFAULTS OF THE host-line-sigma, host-line-min-slit-contrast AND host-line-padding-pixels SETTINGS
DEFAULT_SIGMA_LIMIT = 5.5
DEFAULT_MIN_SLIT_CONTRAST = 0.3
DEFAULT_PADDING_LINES = 3
# A NEIGHBOUR OF A SEED BIN, WHICH TAKES IN THE WEAK WINGS OF A LINE, NEEDS THIS MANY SIGMA
NEIGHBOUR_SIGMA = 3.0
# THE WIDTH OF A WAVELENGTH BIN IN DETECTOR LINES, AND OF A SLIT BIN IN ARCSEC
LINES_PER_BIN = 2
SLIT_BIN_ARCSEC = 1.0
# A SLIT BIN IS THE OBJECT TRACE WHEN MORE THAN THIS SHARE OF ITS PIXELS ARE FLAGGED AS OBJECT. THE CLIP FLAGS A WHOLE
# SLIT RANGE FOR AN OBJECT, BUT ONLY THE CORES OF SKYLINES AND HOST LINES, AND COSMIC RAYS, OUTSIDE IT
OBJECT_TRACE_SHARE = 0.5
# A CELL IS THE MEDIAN OF AT LEAST THIS MANY PIXELS, SO A COSMIC RAY OR A HOT PIXEL DOES NOT MOVE IT
MIN_CELL_PIXELS = 3
# THE RUNNING MEDIAN THAT REMOVES THE CONTINUUM, AND THE RUNNING MAD THAT GIVES THE LOCAL NOISE, SPAN THIS MANY BINS
# (50 DETECTOR LINES, AS THE 51 OF THE NOD CHECK). THE MEDIAN IS ONE-SIDED AT THE ORDER ENDS, SO HALF OF IT IS NOT
# SEARCHED
TREND_WINDOW_BINS = 25
# A FEATURE SPANS AT LEAST THIS MANY BINS, AND FLAGGED BINS THIS FAR APART OR CLOSER BELONG TO ONE FEATURE
MIN_FEATURE_BINS = 2
MAX_FEATURE_GAP_BINS = 1


@dataclass(frozen=True)
class _slit_bins:
    """The per-bin statistics of one order: wavelength of each bin centre and the three significance tests."""

    wavelength: np.ndarray
    commonSigma: np.ndarray
    variationSigma: np.ndarray
    variationPeak: np.ndarray
    excessPeak: np.ndarray
    lineStepNm: float
    binStepNm: float


def find_host_line_windows(
    imageMapOrder: Any,
    skylineWavelengthsNm: Sequence[float] = (),
    sigmaLimit: float = DEFAULT_SIGMA_LIMIT,
    minSlitContrast: float = DEFAULT_MIN_SLIT_CONTRAST,
    paddingLines: int = DEFAULT_PADDING_LINES,
) -> tuple[host_feature, ...]:
    """*find the wavelength windows of host emission lines that fill the slit of one order of a stare frame*

    **Key Arguments:**

    - ``imageMapOrder`` -- the order table that reaches the sky fit, one row per pixel, with the columns ``order``,
      ``x``, ``y``, ``wavelength``, ``slit_position``, ``flux``, ``error``, ``flagged_object_clipped``,
      ``flagged_bad_pixel_clipped`` and ``flagged_edge_clipped``. Pixels flagged as object, bad pixel or edge are
      left out. The detector axis the wavelength runs along is the one with more distinct pixel values
    - ``skylineWavelengthsNm`` -- wavelengths (nm) of catalogue skylines. Bins within 8 detector lines of one are not
      searched. Default *()*
    - ``sigmaLimit`` -- the common excess and the slit variation of a seed bin must both exceed this many times the
      local noise. Default *5.5*
    - ``minSlitContrast`` -- the slit variation, summed over a feature, must be at least this share of its summed peak
      excess. A sky line, which is constant along the slit, has a share near 0. Default *0.3*
    - ``paddingLines`` -- detector lines added to each end of a window, for the wings of the line. Default *3*

    **Return:**

    - ``windows`` -- a tuple of ``host_feature``, one per window and sorted by wavelength, with the order and the
      significance (the smaller of the two tests, at its largest in the window). Empty when no line is found

    **Usage:**

    ```python
    from soxspipe.commonutils.host_line_mask import find_host_line_windows
    windows = find_host_line_windows(imageMapOrder, skylineWavelengthsNm=skylines)
    ```
    """
    bins = _slit_bin_statistics(imageMapOrder)
    if bins is None:
        return ()
    neighbourSigma = min(NEIGHBOUR_SIGMA, sigmaLimit)
    isSearched = _searched_bins(bins, skylineWavelengthsNm)
    significance = np.minimum(bins.commonSigma, bins.variationSigma)
    isSeed = (significance > sigmaLimit) & isSearched
    isNeighbour = (significance > neighbourSigma) & isSearched
    order = int(imageMapOrder["order"].iloc[0])
    halfBinNm = bins.binStepNm / 2
    paddingNm = paddingLines * bins.lineStepNm
    return tuple(
        host_feature(
            order=order,
            wavelengthMinNm=float(bins.wavelength[first] - halfBinNm - paddingNm),
            wavelengthMaxNm=float(bins.wavelength[last] + halfBinNm + paddingNm),
            significance=float(np.max(significance[first : last + 1])),
        )
        for first, last in _feature_runs(isNeighbour, isSeed)
        if _slit_contrast(bins, first, last) >= minSlitContrast
    )


def host_line_pixel_mask(wavelength: Any, windows: Sequence[host_feature]) -> np.ndarray:
    """*mark the pixels whose wavelength lies inside any host-line window*

    **Key Arguments:**

    - ``wavelength`` -- the wavelength (nm) of each pixel
    - ``windows`` -- the ``host_feature`` windows, from ``find_host_line_windows``

    **Return:**

    - ``isInWindow`` -- a boolean array, one value per pixel, True inside a window (ends included). All False when
      there are no windows

    **Usage:**

    ```python
    from soxspipe.commonutils.host_line_mask import host_line_pixel_mask
    isMasked = host_line_pixel_mask(imageMapOrder["wavelength"].to_numpy(), windows)
    ```
    """
    wavelength = np.asarray(wavelength, dtype=float)
    isInWindow = np.zeros(wavelength.shape, dtype=bool)
    for window in windows:
        isInWindow |= (wavelength >= window.wavelengthMinNm) & (wavelength <= window.wavelengthMaxNm)
    return isInWindow


def _slit_bin_statistics(imageMapOrder: Any) -> _slit_bins | None:
    """*the significance tests of every wavelength bin of an order, or None when it has too few usable pixels*"""
    cells = _slit_cells(imageMapOrder)
    if cells is None:
        return None
    flux, variance, wavelength, lineStepNm = cells
    excess = flux - _running(flux, "median")
    common = excess.median(axis=1)
    slitCount = flux.notna().sum(axis=1)
    commonNoise = np.maximum(
        np.sqrt(MEDIAN_VARIANCE_FACTOR * variance.sum(axis=1) / slitCount**2),
        MAD_TO_SIGMA * _running((common - _running(common, "median")).abs(), "median"),
    )
    variation = excess.sub(common, axis=0)
    variationNoise = np.maximum(np.sqrt(variance), MAD_TO_SIGMA * _running(variation.abs(), "median"))
    with np.errstate(divide="ignore", invalid="ignore"):
        return _slit_bins(
            wavelength=wavelength,
            commonSigma=np.nan_to_num((common / commonNoise).to_numpy(), nan=0.0, posinf=0.0, neginf=0.0),
            variationSigma=np.nan_to_num((variation / variationNoise).max(axis=1).to_numpy(), nan=0.0, posinf=0.0),
            variationPeak=np.nan_to_num(variation.max(axis=1).to_numpy(), nan=0.0),
            excessPeak=np.nan_to_num(excess.max(axis=1).to_numpy(), nan=0.0),
            lineStepNm=lineStepNm,
            binStepNm=LINES_PER_BIN * lineStepNm,
        )


def _slit_cells(imageMapOrder: Any) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, float] | None:
    """*the median flux and the variance of each wavelength-by-slit cell, the bin centres and the line step (nm)*

    The pixels are binned by true wavelength, so the tilt of the slit does not move a sky line between bins. Bins
    with no usable cell are kept, empty, so the running median sees the real wavelength spacing.
    """
    wavelength = imageMapOrder["wavelength"].to_numpy(dtype=float)
    flux = imageMapOrder["flux"].to_numpy(dtype=float)
    error = imageMapOrder["error"].to_numpy(dtype=float)
    slit = imageMapOrder["slit_position"].to_numpy(dtype=float)
    isUsable = (
        np.isfinite(wavelength)
        & np.isfinite(flux)
        & np.isfinite(error)
        & np.isfinite(slit)
        & ~imageMapOrder["flagged_bad_pixel_clipped"].to_numpy(dtype=bool)
        & ~imageMapOrder["flagged_edge_clipped"].to_numpy(dtype=bool)
    )
    lineStepNm = _line_step_nm(imageMapOrder)
    if not isUsable.any() or lineStepNm is None:
        return None
    firstNm = float(np.min(wavelength[isUsable]))
    slitBin = np.zeros(len(slit), dtype=int)
    slitBin[isUsable] = np.floor((slit[isUsable] - np.min(slit[isUsable])) / SLIT_BIN_ARCSEC).astype(int)
    # THE CORES OF SKYLINES ARE KEPT: LEFT OUT, THE FLANKS THAT REMAIN ARE SAMPLED DIFFERENTLY IN EACH SLIT BIN
    isObjectPixel = imageMapOrder["flagged_object_clipped"].to_numpy(dtype=bool)
    objectShare = pd.Series(isObjectPixel).groupby(slitBin).transform("mean").to_numpy()
    isUsable &= objectShare <= OBJECT_TRACE_SHARE
    if not isUsable.any():
        return None
    pixels = pd.DataFrame(
        {
            "bin": np.floor(np.where(isUsable, wavelength - firstNm, 0.0) / (LINES_PER_BIN * lineStepNm)).astype(int),
            "slit": slitBin,
            "flux": flux,
            "variance": error**2,
        }
    )[isUsable]
    grouped = pixels.groupby(["bin", "slit"])
    count = grouped["flux"].count()
    isCell = count >= MIN_CELL_PIXELS
    binIndex = np.arange(pixels["bin"].max() + 1)
    cellFlux = grouped["flux"].median()[isCell].unstack("slit").reindex(binIndex)
    cellVariance = (grouped["variance"].sum()[isCell] / count[isCell] ** 2).unstack("slit").reindex(binIndex)
    centres = firstNm + (binIndex + 0.5) * LINES_PER_BIN * lineStepNm
    return cellFlux, cellVariance, centres, lineStepNm


def _line_step_nm(imageMapOrder: Any) -> float | None:
    """*the median wavelength increment (nm) of one detector line along the dispersion, or None when there is none*

    An order is tens of pixels across the slit and thousands along the dispersion, whichever detector axis that is, so
    the axis with more distinct pixel values is the dispersion axis. The increment is taken between neighbouring lines
    of one slit pixel, so the tilt of the slit, which adds its own spread to the wavelength range of the order, does
    not bias it.
    """
    lineAxis, crossAxis = ("x", "y") if imageMapOrder["x"].nunique() >= imageMapOrder["y"].nunique() else ("y", "x")
    ordered = imageMapOrder.sort_values([crossAxis, lineAxis])
    line = ordered[lineAxis].to_numpy()
    cross = ordered[crossAxis].to_numpy()
    wavelength = ordered["wavelength"].to_numpy(dtype=float)
    isNeighbour = (np.diff(line) == 1) & (np.diff(cross) == 0)
    increments = np.abs(np.diff(wavelength))[isNeighbour]
    increments = increments[np.isfinite(increments) & (increments > 0)]
    return float(np.median(increments)) if len(increments) else None


def _running(values: Any, how: str) -> Any:
    """*the centred running median over ``TREND_WINDOW_BINS`` bins, one-sided at the ends*"""
    return getattr(values.rolling(TREND_WINDOW_BINS, center=True, min_periods=1), how)()


def _searched_bins(bins: _slit_bins, skylineWavelengthsNm: Sequence[float]) -> np.ndarray:
    """*a mask of the bins that are searched: not near a catalogue skyline and not at an end of the order*"""
    isSearched = np.ones(len(bins.wavelength), dtype=bool)
    edge = TREND_WINDOW_BINS // 2
    isSearched[:edge] = False
    isSearched[len(isSearched) - edge :] = False
    skylines = np.asarray(skylineWavelengthsNm, dtype=float)
    if len(skylines):
        distance = np.min(np.abs(bins.wavelength[:, None] - skylines[None, :]), axis=1)
        isSearched &= distance > SKY_MASK_HALF_WIDTH_LINES * bins.lineStepNm
    return isSearched


def _slit_contrast(bins: _slit_bins, first: int, last: int) -> float:
    """*the share of the peak excess of bins ``first`` to ``last`` that changes along the slit*"""
    variation = np.clip(bins.variationPeak[first : last + 1], 0.0, None).sum()
    excess = np.clip(bins.excessPeak[first : last + 1], 0.0, None).sum()
    return float(variation / excess) if excess > 0 else 0.0


def _feature_runs(isNeighbour: np.ndarray, isSeed: np.ndarray) -> list[tuple[int, int]]:
    """*the (first, last) bin of each run of neighbour bins that is long enough and holds a seed bin*"""
    indices = np.flatnonzero(isNeighbour)
    if not len(indices):
        return []
    breaks = np.flatnonzero(np.diff(indices) > MAX_FEATURE_GAP_BINS + 1)
    starts = np.concatenate(([0], breaks + 1))
    ends = np.concatenate((breaks, [len(indices) - 1]))
    runs = []
    for start, end in zip(starts, ends, strict=True):
        first, last = int(indices[start]), int(indices[end])
        if last - first + 1 >= MIN_FEATURE_BINS and isSeed[first : last + 1].any():
            runs.append((first, last))
    return runs
