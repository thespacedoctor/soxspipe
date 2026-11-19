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

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from soxspipe.commonutils.nod_host_contamination import (
    MAD_TO_SIGMA,
    MEDIAN_VARIANCE_FACTOR,
    SKY_MASK_HALF_WIDTH_LINES,
    host_feature,
)

# THE RESOURCE THAT LISTS THE REST-FRAME HOST EMISSION LINES
HOST_LINES_FILE = Path(__file__).resolve().parent.parent / "resources" / "host_emission_lines.yaml"
# TWO WINDOWS THAT AGREE ON ONE REDSHIFT, FOR DISTINCT LINES, WITHIN THIS TOLERANCE (150 KM/S), GIVE THE HOST REDSHIFT.
# THE CENTRE OF A WINDOW IS GOOD TO ABOUT 0.1 NM BECAUSE THE HOST ROTATES, WHICH IS 2E-4 AT 500 NM
REDSHIFT_TOLERANCE = 5e-4
# THE REDSHIFT IS KNOWN TO HALF THE TOLERANCE, WHICH IS ABOUT 0.17 NM AT 670 NM. IT IS ADDED IN QUADRATURE TO THE
# RELATIVE HALF WIDTH OF THE DETECTED WINDOWS WHEN A COMPANION WINDOW IS SET
REDSHIFT_UNCERTAINTY = REDSHIFT_TOLERANCE / 2
MIN_REDSHIFT_WINDOWS = 2
# THE HOST REDSHIFTS THAT ARE SEARCHED: A LITTLE BLUE FOR A NEARBY GALAXY, AND UP TO WHERE THE LINES LEAVE THE ARMS
MIN_REDSHIFT = -0.005
MAX_REDSHIFT = 0.5
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
    isUsable = _usable_pixels(imageMapOrder)
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


def _usable_pixels(imageMapOrder: Any) -> np.ndarray:
    """*a mask of the pixels with a finite wavelength, flux, error and slit position that are not bad or edge pixels*"""
    return (
        np.isfinite(imageMapOrder["wavelength"].to_numpy(dtype=float))
        & np.isfinite(imageMapOrder["flux"].to_numpy(dtype=float))
        & np.isfinite(imageMapOrder["error"].to_numpy(dtype=float))
        & np.isfinite(imageMapOrder["slit_position"].to_numpy(dtype=float))
        & ~imageMapOrder["flagged_bad_pixel_clipped"].to_numpy(dtype=bool)
        & ~imageMapOrder["flagged_edge_clipped"].to_numpy(dtype=bool)
    )


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


@dataclass(frozen=True)
class host_redshift:
    """The host redshift inferred from detected windows.

    ``z`` is the median redshift of the matched windows, ``lineNames`` the rest-frame lines they were matched to in
    wavelength order, ``nWindows`` the number of windows that support the redshift (one for each line), and
    ``halfWidthFraction`` the median half width of the matched windows as a fraction of their centre wavelength,
    which sets the width of the companion windows.
    """

    z: float
    lineNames: tuple[str, ...]
    nWindows: int
    halfWidthFraction: float


@dataclass(frozen=True)
class _hypothesis:
    """One redshift and the (window, line) pairs that agree with it."""

    z: float
    matches: tuple[tuple[int, int], ...]
    scatter: float


def load_host_lines() -> tuple[tuple[str, float], ...]:
    """*the rest-frame host emission lines of ``resources/host_emission_lines.yaml``*

    **Return:**

    - ``lines`` -- a tuple of ``(name, wavelength in nm in air)``, in the order of the file

    **Raises:**

    - ``FileNotFoundError`` -- if the file is missing
    - ``ValueError`` -- if the file is not a ``lines`` list of entries with a ``name`` and a numeric ``wavelength-nm``

    **Usage:**

    ```python
    from soxspipe.commonutils.host_line_mask import load_host_lines
    restLines = load_host_lines()
    ```
    """
    import yaml

    with HOST_LINES_FILE.open() as stream:
        content = yaml.safe_load(stream)
    try:
        return tuple((str(entry["name"]), float(entry["wavelength-nm"])) for entry in content["lines"])
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError(
            f"{HOST_LINES_FILE.name} must hold a list `lines` of entries with a name and wavelength-nm: {e!r}"
        ) from e


def infer_host_redshift(
    windows: Sequence[host_feature], restLines: Sequence[tuple[str, float]]
) -> host_redshift | None:
    """*infer the host redshift from the centres of detected windows*

    Each window centre is matched to each rest-frame line, which gives a redshift per pair. A redshift is accepted
    when at least ``MIN_REDSHIFT_WINDOWS`` windows agree on it within ``REDSHIFT_TOLERANCE`` for distinct lines, and it
    lies between ``MIN_REDSHIFT`` and ``MAX_REDSHIFT``. A pair of windows can fit two redshifts, because some pairs of
    lines have wavelength ratios within the tolerance. When the best redshift and the best one that differs from it by
    more than the tolerance have the same number of windows, none is returned.

    **Key Arguments:**

    - ``windows`` -- the detected ``host_feature`` windows, of any orders
    - ``restLines`` -- the ``(name, wavelength)`` pairs from ``load_host_lines``

    **Return:**

    - ``redshift`` -- a ``host_redshift``, or *None* when no redshift is supported by two windows or two redshifts
      are supported equally

    **Usage:**

    ```python
    from soxspipe.commonutils.host_line_mask import infer_host_redshift, load_host_lines
    redshift = infer_host_redshift(windows, load_host_lines())
    ```
    """
    if len(windows) < MIN_REDSHIFT_WINDOWS or not len(restLines):
        return None
    centres = np.array([(window.wavelengthMinNm + window.wavelengthMaxNm) / 2 for window in windows])
    rest = np.array([wavelength for _, wavelength in restLines])
    candidates = centres[:, None] / rest[None, :] - 1.0
    hypotheses = _redshift_hypotheses(candidates)
    if not hypotheses:
        return None
    best = max(hypotheses, key=lambda hypothesis: (len(hypothesis.matches), -hypothesis.scatter))
    rivals = [h for h in hypotheses if abs(h.z - best.z) > REDSHIFT_TOLERANCE]
    if any(len(rival.matches) == len(best.matches) for rival in rivals):
        return None
    if not MIN_REDSHIFT <= best.z <= MAX_REDSHIFT:
        return None
    return _redshift_of(best, windows, centres, restLines)


def _redshift_of(
    hypothesis: _hypothesis,
    windows: Sequence[host_feature],
    centres: np.ndarray,
    restLines: Sequence[tuple[str, float]],
) -> host_redshift:
    """*the ``host_redshift`` of an accepted hypothesis: its lines in wavelength order and the median relative width*"""
    fractions = [
        (windows[window].wavelengthMaxNm - windows[window].wavelengthMinNm) / 2 / centres[window]
        for window, _ in hypothesis.matches
    ]
    ordered = sorted(hypothesis.matches, key=lambda match: restLines[match[1]][1])
    return host_redshift(
        z=hypothesis.z,
        lineNames=tuple(restLines[line][0] for _, line in ordered),
        nWindows=len(hypothesis.matches),
        halfWidthFraction=float(np.median(fractions)),
    )


def _redshift_hypotheses(candidates: np.ndarray) -> list[_hypothesis]:
    """*one hypothesis for each seed redshift in range that at least ``MIN_REDSHIFT_WINDOWS`` windows support*

    ``candidates`` holds the redshift of each (window, line) pair. A hypothesis keeps its seed's matches and takes the
    median of their redshifts as its own.
    """
    hypotheses = []
    for seed in candidates[(candidates >= MIN_REDSHIFT) & (candidates <= MAX_REDSHIFT)]:
        matches = _matches_of_redshift(candidates, seed)
        if len(matches) < MIN_REDSHIFT_WINDOWS:
            continue
        redshifts = [candidates[window, line] for window, line in matches]
        scatter = sum(abs(redshift - seed) for redshift in redshifts)
        hypotheses.append(_hypothesis(z=float(np.median(redshifts)), matches=tuple(matches), scatter=float(scatter)))
    return hypotheses


def _matches_of_redshift(candidates: np.ndarray, seed: float) -> list[tuple[int, int]]:
    """*the (window, line) pairs that agree with ``seed``, at most one window per line and one line per window*"""
    distance = np.abs(candidates - seed)
    pairs = sorted(
        ((window, line) for window, line in zip(*np.nonzero(distance <= REDSHIFT_TOLERANCE), strict=True)),
        key=lambda pair: distance[pair],
    )
    usedWindows, usedLines, matches = set(), set(), []
    for window, line in pairs:
        if window in usedWindows or line in usedLines:
            continue
        usedWindows.add(window)
        usedLines.add(line)
        matches.append((int(window), int(line)))
    return matches


def companion_windows(
    imageMapOrder: Any,
    detected: Sequence[host_feature],
    redshift: host_redshift,
    restLines: Sequence[tuple[str, float]],
    skylineWavelengthsNm: Sequence[float] = (),
) -> tuple[host_feature, ...]:
    """*the windows of the rest-frame lines of the host that were not detected in one order*

    A line gets a window when it lies inside the searched part of the usable pixels of the order, touches no
    detected window of the order, and has no catalogue skyline within 8 detector lines. The window is centred on the
    line at the inferred redshift. Its half width is the detected relative half width and ``REDSHIFT_UNCERTAINTY``
    added in quadrature, times the wavelength.

    **Key Arguments:**

    - ``imageMapOrder`` -- the order table, as for ``find_host_line_windows``
    - ``detected`` -- the windows detected in this order
    - ``redshift`` -- the ``host_redshift`` from ``infer_host_redshift``
    - ``restLines`` -- the pairs from ``load_host_lines``
    - ``skylineWavelengthsNm`` -- wavelengths (nm) of catalogue skylines. Default *()*

    **Return:**

    - ``windows`` -- a tuple of ``host_feature`` of this order, by wavelength

    **Usage:**

    ```python
    from soxspipe.commonutils.host_line_mask import companion_windows
    extra = companion_windows(imageMapOrder, windows, redshift, restLines, skylineWavelengthsNm=skylines)
    ```
    """
    lineStepNm = _line_step_nm(imageMapOrder)
    searched = _searched_range_nm(imageMapOrder, lineStepNm)
    if searched is None:
        return ()
    skylines = np.asarray(skylineWavelengthsNm, dtype=float)
    order = int(imageMapOrder["order"].iloc[0])
    fraction = math.hypot(redshift.halfWidthFraction, REDSHIFT_UNCERTAINTY)
    candidates = (
        host_feature(
            order, restNm * (1.0 + redshift.z) * (1.0 - fraction), restNm * (1.0 + redshift.z) * (1.0 + fraction), 0.0
        )
        for _, restNm in sorted(restLines, key=lambda line: line[1])
    )
    margin = SKY_MASK_HALF_WIDTH_LINES * lineStepNm
    return tuple(c for c in candidates if _is_wanted(c, searched, detected, skylines, margin))


def _is_wanted(
    window: host_feature,
    searched: tuple[float, float],
    detected: Sequence[host_feature],
    skylines: np.ndarray,
    margin: float,
) -> bool:
    """*True when a companion window is in the searched range, touches no detected window and has no skyline near*"""
    first, last = window.wavelengthMinNm, window.wavelengthMaxNm
    if first < searched[0] or last > searched[1]:
        return False
    if any(first <= other.wavelengthMaxNm and last >= other.wavelengthMinNm for other in detected):
        return False
    return not _near_catalogue_skyline(skylines, first, last, margin)


def new_windows(candidates: Sequence[host_feature], detected: Sequence[host_feature]) -> tuple[host_feature, ...]:
    """*the candidate windows that touch none of the detected windows*

    **Key Arguments:**

    - ``candidates`` -- the windows to check, for example the merged companion windows
    - ``detected`` -- the detected windows of every order

    **Return:**

    - ``windows`` -- the candidates that share no wavelength with any detected window, in their order

    **Usage:**

    ```python
    from soxspipe.commonutils.host_line_mask import new_windows
    added = new_windows(merge_host_features(companions), detected)
    ```
    """
    return tuple(
        candidate
        for candidate in candidates
        if not any(
            candidate.wavelengthMinNm <= window.wavelengthMaxNm and candidate.wavelengthMaxNm >= window.wavelengthMinNm
            for window in detected
        )
    )


def _searched_range_nm(imageMapOrder: Any, lineStepNm: float | None) -> tuple[float, float] | None:
    """*the wavelength range of the usable pixels of an order that the detector searches, or None when it has none*"""
    isUsable = _usable_pixels(imageMapOrder)
    if lineStepNm is None or not isUsable.any():
        return None
    wavelength = imageMapOrder["wavelength"].to_numpy(dtype=float)[isUsable]
    edgeNm = (TREND_WINDOW_BINS // 2) * LINES_PER_BIN * lineStepNm
    return float(wavelength.min()) + edgeNm, float(wavelength.max()) - edgeNm


def _near_catalogue_skyline(skylines: np.ndarray, firstNm: float, lastNm: float, marginNm: float) -> bool:
    """*True when a catalogue skyline lies within ``marginNm`` of the window from ``firstNm`` to ``lastNm``*"""
    return bool(len(skylines) and np.any((skylines >= firstNm - marginNm) & (skylines <= lastNm + marginNm)))
