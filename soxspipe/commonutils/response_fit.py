#!/usr/bin/env python
"""
*Robust spline fit of the instrument response curve*

The raw response (catalogue flux over extracted flux) is binned, the log of the bin medians is fitted with a cubic
least-squares spline, and bins that deviate from the fit are clipped in log space.

:Author:
    David Young
"""

import os
from dataclasses import dataclass
from typing import Any

os.environ["TERM"] = "vt100"

SPLINE_DEGREE = 3
# FLOOR ON THE PER-BIN LOG ERROR SO NOISELESS BINS DO NOT DOMINATE THE WEIGHTED FIT
BIN_ERROR_FLOOR = 1e-3
# FLOOR ON THE ROBUST SIGMA OF ERROR-SCALED RESIDUALS, SO OVERESTIMATED BIN ERRORS DO NOT TIGHTEN THE CLIP
SCALED_SIGMA_FLOOR = 1.0
MAD_TO_SIGMA = 1.4826


class _ResponseFitConvergenceError(Exception):
    """Mark failures raised inside the response fitting loop."""

    def __init__(self, originalError: Exception) -> None:
        super().__init__(str(originalError))
        self.originalError = originalError


def mask_response_samples(
    wave: Any,
    ratio: Any,
    telluricRegions: Any,
    stellarLines: Any,
    lineHalfWidth: float,
    log: Any = None,
) -> Any:
    """*flag the raw response samples that are safe to fit*

    **Key Arguments:**

    - ``wave`` -- sample wavelengths (nm)
    - ``ratio`` -- raw response (catalogue flux over extracted flux) at each wavelength
    - ``telluricRegions`` -- list of ``[lower, upper]`` wavelength ranges (nm) to drop
    - ``stellarLines`` -- wavelengths (nm) of stellar lines to drop
    - ``lineHalfWidth`` -- half width (nm) removed around each stellar line
    - ``log`` -- optional logger; warns when non-finite samples are dropped

    **Return:**

    - ``mask`` -- boolean array, True for samples to fit. Non-finite samples and samples with a ratio of zero are
      always False. A negative ratio (negative extracted counts from noise) is kept, so low-signal bins are not
      biased by dropping their negative half.
    """
    import numpy as np

    wave = np.asarray(wave, dtype=float)
    ratio = np.asarray(ratio, dtype=float)
    isFinite = np.isfinite(wave) & np.isfinite(ratio)
    nonFiniteCount = int(np.count_nonzero(~isFinite))
    if nonFiniteCount > 0 and log is not None:
        log.warning(
            f"Dropped {nonFiniteCount} raw response samples with a non-finite wavelength or response "
            "(for example from zero extracted flux) before fitting the response function."
        )
    with np.errstate(invalid="ignore"):
        mask = isFinite & (ratio != 0)
        for lower, upper in telluricRegions:
            mask &= ~((wave >= lower) & (wave <= upper))
        for line in stellarLines:
            mask &= ~(np.abs(wave - line) <= lineHalfWidth)
    return mask


def bin_log_response(wave: Any, ratio: Any, mask: Any, binWidth: float, minSamples: int) -> tuple[Any, Any, Any]:
    """*bin the unmasked response samples into median log-response points*

    Each bin takes the median of the inverse ratio (extracted over catalogue flux), which is linear in the noisy
    counts and stays defined when they are negative. Taking the median of log(ratio) instead would drop the
    negative-count samples and bias low-signal bins towards a low response. Bins whose median inverse ratio is not
    positive are dropped.

    **Key Arguments:**

    - ``wave`` -- sample wavelengths (nm)
    - ``ratio`` -- raw response at each wavelength
    - ``mask`` -- boolean array, True for samples to use
    - ``binWidth`` -- bin width (nm)
    - ``minSamples`` -- bins with fewer samples than this are dropped

    **Return:**

    - ``binWave`` -- median wavelength of each bin
    - ``binLog`` -- median log(ratio) of each bin
    - ``binErr`` -- MAD-based error of the median, with a small floor
    """
    import numpy as np

    mask = np.asarray(mask, dtype=bool)
    usedWave = np.asarray(wave, dtype=float)[mask]
    usedInverse = 1.0 / np.asarray(ratio, dtype=float)[mask]
    if usedWave.size == 0:
        return np.array([]), np.array([]), np.array([])

    binIndex = np.floor((usedWave - usedWave.min()) / binWidth).astype(int)
    binWave, binLog, binErr = [], [], []
    for index in np.unique(binIndex):
        inBin = binIndex == index
        count = int(np.count_nonzero(inBin))
        if count < minSamples:
            continue
        values = usedInverse[inBin]
        medianValue = np.median(values)
        if medianValue <= 0:
            continue
        mad = MAD_TO_SIGMA * np.median(np.abs(values - medianValue))
        binWave.append(np.median(usedWave[inBin]))
        binLog.append(-np.log(medianValue))
        # FRACTIONAL ERROR OF THE MEDIAN IS THE ERROR IN LOG SPACE
        binErr.append(mad / medianValue / np.sqrt(count) + BIN_ERROR_FLOOR)
    return np.array(binWave), np.array(binLog), np.array(binErr)


def _interior_knots(x: Any, knotSpacing: float) -> Any:
    """Return a regular knot grid, keeping only knots with data on both sides within half a spacing."""
    import numpy as np

    grid = np.arange(x.min() + knotSpacing / 2.0, x.max() - knotSpacing / 2.0, knotSpacing)
    halfSpacing = knotSpacing / 2.0
    keep = [
        knot
        for knot in grid
        if np.any((x >= knot - halfSpacing) & (x < knot)) and np.any((x > knot) & (x <= knot + halfSpacing))
    ]
    return np.array(keep, dtype=float)


def _fit_spline(binWave: Any, binLog: Any, binErr: Any, keep: Any, knotSpacing: float) -> Any:
    """Fit the weighted least-squares cubic spline to the kept bins."""
    from scipy.interpolate import LSQUnivariateSpline

    x = binWave[keep]
    knots = _interior_knots(x, knotSpacing) if x.size > 0 else []
    if x.size < SPLINE_DEGREE + 1 + len(knots) or x.size == 0:
        raise _ResponseFitConvergenceError(
            ValueError(
                f"Only {x.size} response bins remain for a cubic spline with {len(knots)} interior knots; "
                "too few bins to fit the response function."
            )
        )
    try:
        return LSQUnivariateSpline(x, binLog[keep], knots, w=1.0 / binErr[keep], k=SPLINE_DEGREE)
    except Exception as error:
        raise _ResponseFitConvergenceError(error) from error


def _clip_mask(
    spline: Any, binWave: Any, binLog: Any, binErr: Any, keep: Any, clipSigmaLow: float, clipSigmaHigh: float
) -> Any:
    """Return the bins within the asymmetric robust clipping limits of the spline, in log space.

    Each residual is divided by its bin error before the robust sigma is taken, so a noisy bin (for example at the
    low-signal end of an arm) is judged against its own noise, not against the scatter of the bright bins.
    """
    import numpy as np

    # UNUSABLE BINS (NON-POSITIVE OR NAN ERRORS) GIVE NON-FINITE RESIDUALS THAT THE CALLER MASKS OUT
    with np.errstate(divide="ignore", invalid="ignore"):
        scaled = (binLog - spline(binWave)) / binErr
        keptScaled = scaled[keep]
        sigma = max(MAD_TO_SIGMA * np.median(np.abs(keptScaled - np.median(keptScaled))), SCALED_SIGMA_FLOOR)
        return (scaled <= clipSigmaHigh * sigma) & (scaled >= -clipSigmaLow * sigma)


def _usable_bins(binWave: Any, binLog: Any, binErr: Any) -> tuple[Any, Any, Any, Any]:
    """Return the bins as float arrays plus a mask of those with finite values and a positive error."""
    import numpy as np

    binWave = np.asarray(binWave, dtype=float)
    binLog = np.asarray(binLog, dtype=float)
    binErr = np.asarray(binErr, dtype=float)
    usable = np.isfinite(binWave) & np.isfinite(binLog) & np.isfinite(binErr) & (binErr > 0)
    if not usable.any():
        raise _ResponseFitConvergenceError(
            ValueError("The response has no finite samples to fit: every sample is non-finite, non-positive or masked.")
        )
    return binWave, binLog, binErr, usable


def _refit_after_iteration_cap(
    binWave: Any, binLog: Any, binErr: Any, keep: Any, knotSpacing: float, maxIterations: int, log: Any
) -> Any:
    """Refit the bins that survived the last clip, because the spline in hand predates it."""
    if log is not None:
        log.warning(
            f"The response function fit reached the iteration cap ({int(maxIterations)}) while still "
            f"rejecting bins; refitting the {int(keep.sum())} surviving bins."
        )
    return _fit_spline(binWave, binLog, binErr, keep, knotSpacing)


def fit_log_spline(
    binWave: Any,
    binLog: Any,
    binErr: Any,
    knotSpacing: float,
    clipSigmaLow: float,
    clipSigmaHigh: float,
    maxIterations: int,
    log: Any = None,
) -> tuple[Any, Any]:
    """*fit the log response with a clipped cubic least-squares spline*

    The returned spline always fits the returned kept bins. If the iteration cap is reached while the clipping is
    still changing, the spline is refitted on the surviving bins and a warning is logged.

    **Key Arguments:**

    - ``binWave`` -- bin wavelengths (nm), increasing
    - ``binLog`` -- median log response of each bin
    - ``binErr`` -- error of each bin, used as ``1/err`` weights
    - ``knotSpacing`` -- spacing (nm) of the regular interior knot grid
    - ``clipSigmaLow`` -- bins this many robust sigma below the fit are clipped
    - ``clipSigmaHigh`` -- bins this many robust sigma above the fit are clipped
    - ``maxIterations`` -- maximum number of clipping iterations
    - ``log`` -- optional logger

    **Return:**

    - ``spline`` -- the fitted log-response spline
    - ``keptMask`` -- boolean array over the input bins, True for bins used by the final fit
    """
    import numpy as np

    binWave, binLog, binErr, usable = _usable_bins(binWave, binLog, binErr)
    keep = usable.copy()
    safeWave = np.where(usable, binWave, 0.0)

    for _ in range(max(int(maxIterations), 1)):
        spline = _fit_spline(binWave, binLog, binErr, keep, knotSpacing)
        newKeep = usable & _clip_mask(spline, safeWave, binLog, binErr, keep, clipSigmaLow, clipSigmaHigh)
        if np.array_equal(newKeep, keep):
            return spline, keep
        keep = newKeep

    # THE CAP WAS REACHED MID-REJECTION: THE SPLINE PREDATES THE LAST CLIP, SO REFIT THE SURVIVORS
    spline = _refit_after_iteration_cap(binWave, binLog, binErr, keep, knotSpacing, maxIterations, log)
    return spline, keep


def evaluate_log_spline(spline: Any, wave: Any, fitRange: tuple[float, float]) -> Any:
    """*evaluate the response (linear units) from a log-response spline*

    Inside ``fitRange`` the response is ``exp(spline)``. Outside it, the log response continues linearly from the
    range edge using the spline value and first derivative there, so a cubic cannot blow up.

    **Key Arguments:**

    - ``spline`` -- the fitted log-response spline
    - ``wave`` -- wavelengths (nm) to evaluate
    - ``fitRange`` -- (lower, upper) wavelength of the kept bins

    **Return:**

    - ``response`` -- the response at each wavelength
    """
    import numpy as np

    wave = np.asarray(wave, dtype=float)
    lower, upper = fitRange
    clamped = np.clip(wave, lower, upper)
    logResponse = spline(clamped)
    below = wave < lower
    above = wave > upper
    logResponse = np.where(below, spline(lower) + spline(lower, nu=1) * (wave - lower), logResponse)
    logResponse = np.where(above, spline(upper) + spline(upper, nu=1) * (wave - upper), logResponse)
    return np.exp(logResponse)


def resolve_response_settings(recipeSettings: Any, arm: str, instrument: str) -> dict:
    """*return the response fit settings for one arm, with missing keys filled from the packaged defaults*

    A workspace ``soxspipe.yaml`` written by an older release has no spline keys. Its ``poly_order`` becomes
    ``legacy_poly_order`` and every other missing key comes from the packaged default settings.

    **Key Arguments:**

    - ``recipeSettings`` -- the ``soxs-response`` block of the workspace settings
    - ``arm`` -- the arm name, for example ``VIS``
    - ``instrument`` -- ``soxs`` or ``xsh``, which selects the packaged default settings file

    **Return:**

    - ``armSettings`` -- the complete settings dictionary for the arm
    """
    from pathlib import Path

    import yaml

    packageDir = Path(__file__).resolve().parent.parent
    defaultsPath = packageDir / f"{instrument.lower()}_default_settings.yaml"
    with open(defaultsPath) as stream:
        defaults = yaml.safe_load(stream)["soxs-response"].get(arm.lower(), {})

    workspace = dict((recipeSettings or {}).get(arm.lower(), {}))
    if "poly_order" in workspace:
        legacyOrder = workspace.pop("poly_order")
        workspace.setdefault("legacy_poly_order", legacyOrder)
    return {**defaults, **workspace}


@dataclass(frozen=True)
class ResponseFit:
    """The raw samples, the binned points and the fitted spline of one response curve."""

    wave: Any
    ratio: Any
    mask: Any
    binWave: Any
    binLog: Any
    keptMask: Any
    spline: Any
    fitRange: tuple[float, float]

    def __post_init__(self) -> None:
        """Replace each array field with a read-only copy, so the fit cannot be changed through any alias."""
        import numpy as np

        for field in ("wave", "ratio", "mask", "binWave", "binLog", "keptMask"):
            frozenArray = np.array(getattr(self, field))
            frozenArray.flags.writeable = False
            object.__setattr__(self, field, frozenArray)

    def evaluate(self, wave: Any) -> Any:
        """Return the fitted response (linear units) at the given wavelengths."""
        return evaluate_log_spline(self.spline, wave, self.fitRange)


def fit_response_curve(wave: Any, ratio: Any, fitSettings: Any, log: Any = None) -> ResponseFit:
    """*mask, bin and spline-fit a raw response curve*

    **Key Arguments:**

    - ``wave`` -- sample wavelengths (nm)
    - ``ratio`` -- raw response (catalogue flux over extracted flux) at each wavelength
    - ``fitSettings`` -- the settings for the arm, as returned by ``resolve_response_settings``
    - ``log`` -- optional logger

    **Return:**

    - ``fit`` -- a ``ResponseFit``. Raises ``_ResponseFitConvergenceError`` if the spline cannot be fitted.
    """
    import numpy as np

    wave = np.asarray(wave, dtype=float)
    ratio = np.asarray(ratio, dtype=float)
    mask = mask_response_samples(
        wave,
        ratio,
        fitSettings["telluric_regions"],
        fitSettings["stellar_lines"],
        fitSettings["line_half_width_nm"],
        log=log,
    )
    binWave, binLog, binErr = bin_log_response(
        wave, ratio, mask, fitSettings["bin_width_nm"], fitSettings["min_bin_samples"]
    )
    spline, keptMask = fit_log_spline(
        binWave,
        binLog,
        binErr,
        fitSettings["knot_spacing_nm"],
        fitSettings["clip_sigma_low"],
        fitSettings["clip_sigma_high"],
        fitSettings["max_iteration"],
        log=log,
    )
    fitRange = (float(binWave[keptMask].min()), float(binWave[keptMask].max()))
    return ResponseFit(wave, ratio, mask, binWave, binLog, keptMask, spline, fitRange)
