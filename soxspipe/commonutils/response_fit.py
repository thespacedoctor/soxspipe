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
# FEWEST REFERENCE BINS THAT DEFINE THE STRAIGHT-LINE CONTINUUM BESIDE A TELLURIC REGION
WING_MIN_REFERENCE_BINS = 3
# A WING BIN MUST SIT THIS MANY BIN ERRORS ABOVE THE CONTINUUM, WHATEVER THE FIXED THRESHOLD
WING_ERROR_MULTIPLE = 3.0
# THE WING WALK STOPS AT A STEP LARGER THAN THIS MANY MEDIAN BIN SPACINGS (A HOLE IN THE BINS)
WING_MAX_STEP_BINS = 2


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


def _in_any_region(binWave: Any, telluricRegions: Any) -> Any:
    """Return a mask of the bins that lie inside any telluric region."""
    import numpy as np

    inRegion = np.zeros(binWave.size, dtype=bool)
    for lower, upper in telluricRegions:
        inRegion |= (binWave >= lower) & (binWave <= upper)
    return inRegion


def _wing_side_mask(
    binWave: Any,
    binLog: Any,
    binErr: Any,
    isUsable: Any,
    distance: Any,
    wingLimits: tuple[float, float, float, float],
    sideLabel: str,
    log: Any,
) -> Any:
    """Return the wing bins on one side of a region edge, by walking outward from the edge over raised bins.

    ``distance`` is each bin's distance (nm) outward from the edge, so it is positive on the outer side. The
    reference window and the wing reach are measured in that distance. ``wingLimits`` holds the wing reach, the
    reference width, the rise threshold and the largest step (nm) the walk may take between neighbouring bins.
    """
    import numpy as np

    wingMaxNm, wingRefNm, wingThreshold, maxStepNm = wingLimits
    wingMask = np.zeros(binWave.size, dtype=bool)
    # THE REFERENCE WINDOW STARTS WHERE THE WING REACH ENDS
    isReference = isUsable & (distance > wingMaxNm) & (distance <= wingMaxNm + wingRefNm)
    if np.count_nonzero(isReference) < WING_MIN_REFERENCE_BINS:
        if log is not None:
            log.debug(
                f"Skipped telluric wing growth on the {sideLabel}: fewer than {WING_MIN_REFERENCE_BINS} reference bins."
            )
        return wingMask

    coefficients = np.polyfit(binWave[isReference], binLog[isReference], 1, w=1.0 / binErr[isReference])
    # WALK OVER EVERY BIN IN REACH, USABLE OR NOT, SO A HOLE OR AN UNUSABLE BIN ENDS THE CONTIGUOUS WING
    inReach = np.flatnonzero((distance > 0) & (distance <= wingMaxNm))
    previousDistance = 0.0
    for index in inReach[np.argsort(distance[inReach])]:
        if not isUsable[index] or distance[index] - previousDistance > maxStepNm:
            break
        excess = binLog[index] - np.polyval(coefficients, binWave[index])
        if not excess > max(wingThreshold, WING_ERROR_MULTIPLE * binErr[index]):
            break
        wingMask[index] = True
        previousDistance = distance[index]
    return wingMask


def grow_telluric_wings(
    binWave: Any,
    binLog: Any,
    binErr: Any,
    telluricRegions: Any,
    wingMaxNm: float,
    wingRefNm: float,
    wingThreshold: float,
    log: Any = None,
) -> Any:
    """*flag the bins in the raised wings just outside each telluric region*

    Absorption raises the response in the wings of a masked band. For each side of each region, a weighted straight
    line is fitted to the bins in a reference window beyond the reach of the wing. Walking outward from the region
    edge, a bin is flagged while it sits above that continuum by more than the larger of ``wingThreshold`` and
    three times its own error. The walk stops at the first bin that is not raised, at the first unusable bin, and at
    a hole in the bins wider than twice the median bin spacing. Only positive excursions count.

    **Key Arguments:**

    - ``binWave`` -- bin wavelengths (nm)
    - ``binLog`` -- median log response of each bin
    - ``binErr`` -- error of each bin
    - ``telluricRegions`` -- list of ``[lower, upper]`` wavelength ranges (nm) that were masked
    - ``wingMaxNm`` -- the furthest distance (nm) from a region edge a wing can reach. Zero or less disables growth
    - ``wingRefNm`` -- width (nm) of the reference window that starts where the wing reach ends
    - ``wingThreshold`` -- minimum rise (log units) above the continuum for a bin to count as wing
    - ``log`` -- optional logger; a side with fewer than three reference bins is skipped with a debug message

    **Return:**

    - ``wingMask`` -- boolean array over the bins, True for wing bins to drop from the fit
    """
    import numpy as np

    binWave = np.asarray(binWave, dtype=float)
    binLog = np.asarray(binLog, dtype=float)
    binErr = np.asarray(binErr, dtype=float)
    wingMask = np.zeros(binWave.size, dtype=bool)
    if wingMaxNm <= 0:
        return wingMask

    finiteWave = np.sort(binWave[np.isfinite(binWave)])
    maxStepNm = WING_MAX_STEP_BINS * float(np.median(np.diff(finiteWave))) if finiteWave.size > 1 else np.inf
    with np.errstate(invalid="ignore"):
        isUsable = (
            np.isfinite(binWave)
            & np.isfinite(binLog)
            & np.isfinite(binErr)
            & (binErr > 0)
            & ~_in_any_region(binWave, telluricRegions)
        )
        for lower, upper in telluricRegions:
            sides = (("blue side", lower - binWave), ("red side", binWave - upper))
            for sideName, distance in sides:
                wingMask |= _wing_side_mask(
                    binWave,
                    binLog,
                    binErr,
                    isUsable,
                    distance,
                    (wingMaxNm, wingRefNm, wingThreshold, maxStepNm),
                    f"{sideName} of the {lower:.0f}-{upper:.0f} nm region",
                    log,
                )
    return wingMask


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


def _fit_lsq_spline(binWave: Any, binLog: Any, binErr: Any, keep: Any, knotSpacing: float) -> Any:
    """Fit the weighted least-squares cubic spline to the kept bins, with interior knots only where data are."""
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


def _second_difference_matrix(knots: Any, spacing: float) -> Any:
    """Return the matrix of second divided differences of the B-spline coefficients, on the Greville abscissae.

    On a uniform grid this equals the ordinary second difference. The clamped boundary knots make the Greville
    abscissae uneven at the ends, and the divided differences keep every straight line exactly unpenalised there.
    Rows are scaled by ``spacing**2`` so the penalty strength does not depend on the knot spacing.
    """
    import numpy as np

    greville = np.array(
        [np.mean(knots[index + 1 : index + SPLINE_DEGREE + 1]) for index in range(len(knots) - SPLINE_DEGREE - 1)]
    )
    coefficientCount = greville.size
    difference = np.zeros((coefficientCount - 2, coefficientCount))
    for row in range(coefficientCount - 2):
        lowerStep = greville[row + 1] - greville[row]
        upperStep = greville[row + 2] - greville[row + 1]
        scale = 2.0 * spacing**2 / (greville[row + 2] - greville[row])
        difference[row, row] = scale / lowerStep
        difference[row, row + 1] = -scale * (1.0 / lowerStep + 1.0 / upperStep)
        difference[row, row + 2] = scale / upperStep
    return difference


def _fit_penalised_spline(binWave: Any, binLog: Any, binErr: Any, keep: Any, knotSpacing: float, penalty: float) -> Any:
    """Fit a cubic P-spline to the kept bins: a regular knot grid and a second-difference penalty on the coefficients.

    The penalty is ``penalty * median(w**2)``, so it is dimensionless. Where no bin constrains the curve, the penalty
    leaves the coefficients, and so the log curve, close to a straight line.
    """
    import numpy as np
    from scipy.interpolate import BSpline

    x = binWave[keep]
    if x.size < SPLINE_DEGREE + 1:
        raise _ResponseFitConvergenceError(
            ValueError(
                f"Only {x.size} response bins remain for a cubic P-spline; too few bins to fit the response function."
            )
        )
    if x.max() <= x.min():
        raise _ResponseFitConvergenceError(
            ValueError("All remaining response bins lie at one wavelength; cannot fit a cubic P-spline.")
        )
    try:
        segmentCount = max(int(np.ceil((x.max() - x.min()) / knotSpacing)), 1)
        edges = np.linspace(x.min(), x.max(), segmentCount + 1)
        knots = np.concatenate([np.full(SPLINE_DEGREE, edges[0]), edges, np.full(SPLINE_DEGREE, edges[-1])])
        design = BSpline.design_matrix(x, knots, SPLINE_DEGREE).toarray()
        weights = 1.0 / binErr[keep]
        difference = _second_difference_matrix(knots, edges[1] - edges[0])
        # AUGMENTED SYSTEM [W B; SQRT(LAMBDA) D] C = [W Y; 0] AVOIDS SQUARING THE CONDITION NUMBER
        penaltyRoot = np.sqrt(penalty * np.median(weights**2))
        augmentedMatrix = np.vstack([weights[:, None] * design, penaltyRoot * difference])
        augmentedTarget = np.concatenate([weights * binLog[keep], np.zeros(difference.shape[0])])
        coefficients = np.linalg.lstsq(augmentedMatrix, augmentedTarget, rcond=None)[0]
        return BSpline(knots, coefficients, SPLINE_DEGREE)
    except (np.linalg.LinAlgError, ValueError) as error:
        raise _ResponseFitConvergenceError(error) from error


def _fit_spline(binWave: Any, binLog: Any, binErr: Any, keep: Any, knotSpacing: float, penalty: float = 0.0) -> Any:
    """Fit the weighted cubic spline to the kept bins, penalised when ``penalty`` is above zero."""
    if penalty > 0:
        return _fit_penalised_spline(binWave, binLog, binErr, keep, knotSpacing, penalty)
    return _fit_lsq_spline(binWave, binLog, binErr, keep, knotSpacing)


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
    binWave: Any, binLog: Any, binErr: Any, keep: Any, knotSpacing: float, maxIterations: int, log: Any, penalty: float
) -> Any:
    """Refit the bins that survived the last clip, because the spline in hand predates it."""
    if log is not None:
        log.warning(
            f"The response function fit reached the iteration cap ({int(maxIterations)}) while still "
            f"rejecting bins; refitting the {int(keep.sum())} surviving bins."
        )
    return _fit_spline(binWave, binLog, binErr, keep, knotSpacing, penalty)


def fit_log_spline(
    binWave: Any,
    binLog: Any,
    binErr: Any,
    knotSpacing: float,
    clipSigmaLow: float,
    clipSigmaHigh: float,
    maxIterations: int,
    log: Any = None,
    excludeMask: Any = None,
    penalty: float = 0.0,
) -> tuple[Any, Any]:
    """*fit the log response with a clipped cubic spline*

    The returned spline always fits the returned kept bins. If the iteration cap is reached while the clipping is
    still changing, the spline is refitted on the surviving bins and a warning is logged.

    **Key Arguments:**

    - ``binWave`` -- bin wavelengths (nm), increasing
    - ``binLog`` -- median log response of each bin
    - ``binErr`` -- error of each bin, used as ``1/err`` weights
    - ``knotSpacing`` -- spacing (nm) of the regular interior knot grid. With ``penalty`` above zero the knots are a
      dense regular grid over the kept bins, with no filter on where data lie
    - ``clipSigmaLow`` -- bins this many robust sigma below the fit are clipped
    - ``clipSigmaHigh`` -- bins this many robust sigma above the fit are clipped
    - ``maxIterations`` -- maximum number of clipping iterations
    - ``log`` -- optional logger
    - ``excludeMask`` -- optional boolean array over the input bins, True for bins to leave out of the fit from the
      start (for example telluric wing bins)
    - ``penalty`` -- dimensionless weight of a second-difference penalty on the spline coefficients. Zero fits the
      plain least-squares spline; above zero fits a penalised spline that is close to a straight line in log where
      no bins constrain it

    **Return:**

    - ``spline`` -- the fitted log-response spline, a ``LSQUnivariateSpline`` or, with a penalty, a ``BSpline``
    - ``keptMask`` -- boolean array over the input bins, True for bins used by the final fit
    """
    import numpy as np

    binWave, binLog, binErr, usable = _usable_bins(binWave, binLog, binErr)
    if excludeMask is not None:
        usable = usable & ~np.asarray(excludeMask, dtype=bool)
    keep = usable.copy()
    safeWave = np.where(usable, binWave, 0.0)

    for _ in range(max(int(maxIterations), 1)):
        spline = _fit_spline(binWave, binLog, binErr, keep, knotSpacing, penalty)
        newKeep = usable & _clip_mask(spline, safeWave, binLog, binErr, keep, clipSigmaLow, clipSigmaHigh)
        if np.array_equal(newKeep, keep):
            return spline, keep
        keep = newKeep

    # THE CAP WAS REACHED MID-REJECTION: THE SPLINE PREDATES THE LAST CLIP, SO REFIT THE SURVIVORS
    spline = _refit_after_iteration_cap(binWave, binLog, binErr, keep, knotSpacing, maxIterations, log, penalty)
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
    wingMask: Any

    def __post_init__(self) -> None:
        """Replace each array field with a read-only copy, so the fit cannot be changed through any alias."""
        import numpy as np

        for field in ("wave", "ratio", "mask", "binWave", "binLog", "keptMask", "wingMask"):
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
    # SETTINGS WITHOUT THE WING KEYS (HAND-BUILT BLOCKS) FIT WITH WING GROWTH DISABLED
    wingMask = grow_telluric_wings(
        binWave,
        binLog,
        binErr,
        fitSettings["telluric_regions"],
        fitSettings.get("telluric_wing_max_nm", 0),
        fitSettings.get("telluric_wing_ref_nm", 0),
        fitSettings.get("telluric_wing_threshold", 0),
        log=log,
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
        excludeMask=wingMask,
        penalty=fitSettings.get("smoothing_penalty", 0.0),
    )
    fitRange = (float(binWave[keptMask].min()), float(binWave[keptMask].max()))
    return ResponseFit(wave, ratio, mask, binWave, binLog, keptMask, spline, fitRange, wingMask)
