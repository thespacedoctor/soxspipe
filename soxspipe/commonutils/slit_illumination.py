"""Slit illumination: how far each order's sky image sits along the slit from the master flat's.

Spatial flexure shifts the science slit image along the slit, relative to the master flat and the 2D dispersion map,
by an offset δ that depends on the rotator and differs from order to order. This module measures δ for every order from
the frame itself. It is pure numpy and scipy on 2D arrays, and it returns new objects without modifying its inputs.

The measurement has three steps:

1. The master flat's slit profile in each order gives the flat's lit edges, from an erf top-hat fit.
2. The skylines give the science slit profile in each order. For every line and slit bin the side-band median is
   subtracted from the core median, which removes the object and the sky continuum without an object mask.
3. The science profile is fitted as the flat's top-hat shifted by δ. A positive δ means the science image sits
   toward +SLIT.
"""

import math
import warnings
from dataclasses import dataclass
from numbers import Real

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.special import ndtr

from soxspipe.commonutils.slit_illumination_profiles import (
    PROFILE_BIN_ARCSEC,
    PROFILE_COLUMNS,
    in_order_pixels,
    normalised_flat_profile,
    skyline_slit_profiles,
)

__all__ = [
    "FLAT_EDGE_COLUMNS",
    "OFFSET_COLUMNS",
    "PROFILE_BIN_ARCSEC",
    "PROFILE_COLUMNS",
    "SKYLINE_MIN_FLUX_VIS",
    "STATUS_BORROWED",
    "STATUS_FAILED",
    "STATUS_FITTED",
    "fit_order_offsets",
    "flat_slit_edges",
    "measure_slit_illumination",
    "skyline_slit_profiles",
    "slit_illumination_result",
    "top_hat",
]

FLAT_LIT_FRACTION = 0.5
FLAT_START_WIDTH_ARCSEC = 0.15
MIN_EDGE_WIDTH_ARCSEC = 0.01
MAX_FLAT_EDGE_WIDTH_ARCSEC = 2.0

MAX_OFFSET_ARCSEC = 3.0
# THE OVERLAP MASK THAT FOLLOWS THIS MEASUREMENT TRIMS THIS FAR INSIDE THE BOTH-LIT RANGE
OVERLAP_MARGIN_ARCSEC = 0.3
# THE LARGEST UNCERTAINTY ON A FITTED OFFSET THAT STILL GIVES A USABLE MASK EDGE
MAX_OFFSET_ERR_ARCSEC = 0.3
MIN_PROFILE_BINS = 20
# THE OFFSET FIT FREES THE AMPLITUDE, THE OFFSET AND THE EDGE WIDTH
N_FIT_PARAMETERS = 3
OFFSET_START_ARCSEC = (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
# THE FIT MAY WANDER THIS MANY TIMES THE MAXIMUM OFFSET, SO THAT AN OVER-LARGE OFFSET IS REPORTED, NOT CLIPPED
OFFSET_BOUND_FACTOR = 2.0
# THE STARTING OFFSETS ARE KEPT THIS FRACTION OF THE WAY INTO THE BOUNDS, AS THE OPTIMISER NEEDS A START INSIDE THEM
START_INSIDE_BOUND_FRACTION = 0.99
MAX_SCIENCE_EDGE_WIDTH_ARCSEC = 2.0
# THE SCIENCE EDGE WIDTH MUST STAY WITHIN THESE FACTORS OF THE FLAT'S EDGE WIDTH
MIN_WIDTH_RATIO = 0.25
MAX_WIDTH_RATIO = 4.0

# THE VIS SKYLINE TABLE HAS ONLY 5 LINES ABOVE 5 IN ITS R ORDER BUT 33 ABOVE 1, SO THE FAINTER LINES ARE NEEDED. THE
# NIR THRESHOLD IS LEFT AT THE TABLE LOADER'S DEFAULT
SKYLINE_MIN_FLUX_VIS = 1.0

STATUS_FITTED = "fitted"
STATUS_BORROWED = "borrowed"
STATUS_FAILED = "failed"
FLAT_EDGE_COLUMNS = ["order", "lower", "upper", "width", "reason"]
OFFSET_COLUMNS = ["order", "offset", "offsetErr", "scienceWidth", "nBins", "status", "reason", "borrowedFrom"]


def top_hat(slit, amplitude, lower, upper, width):
    """*an erf top-hat: ``amplitude`` between ``lower`` and ``upper``, with Gaussian-blurred edges of sigma ``width``*

    **Key Arguments:**

    - ``slit`` -- slit positions in arcsec
    - ``amplitude`` -- the level between the edges
    - ``lower`` -- the 50% position of the lower edge in arcsec
    - ``upper`` -- the 50% position of the upper edge in arcsec
    - ``width`` -- the Gaussian sigma of both edges in arcsec

    **Return:**

    - ``profile`` -- the top-hat evaluated at ``slit``

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination import top_hat
    profile = top_hat(slit, 1.0, -5.0, 5.0, 0.15)
    ```
    """
    return amplitude * (ndtr((slit - lower) / width) - ndtr((slit - upper) / width))


def _quiet_curve_fit(*args, **kwargs):
    """*``curve_fit`` with the ``OptimizeWarning`` suppressed for this call only*"""
    # A NON-FINITE PCOV BECOMES AN INF OFFSET ERROR, WHICH THE ERROR GATE REJECTS, SO THE WARNING IS REDUNDANT
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)
        return curve_fit(*args, **kwargs)


def _fit_flat_edges(centres, profile, binArcsec):
    """*fit an erf top-hat to a normalised flat profile; return ``(lower, upper, width, reason)``*"""
    isFinite = np.isfinite(profile)
    # THE PROFILE IS NORMALISED TO ITS POSITIVE INTERIOR MEDIAN, SO SOME BINS ALWAYS REACH THE LIT FRACTION
    litCentres = centres[isFinite & (profile >= FLAT_LIT_FRACTION)]
    start = (1.0, litCentres[0], litCentres[-1], FLAT_START_WIDTH_ARCSEC)
    low = (0.0, centres[0] - binArcsec, centres[0] - binArcsec, MIN_EDGE_WIDTH_ARCSEC)
    high = (np.inf, centres[-1] + binArcsec, centres[-1] + binArcsec, MAX_FLAT_EDGE_WIDTH_ARCSEC)
    try:
        popt, _ = _quiet_curve_fit(top_hat, centres[isFinite], profile[isFinite], p0=start, bounds=(low, high))
    except (RuntimeError, ValueError):
        return np.nan, np.nan, np.nan, "flat edge fit did not converge"
    _, lower, upper, width = popt
    if not np.isfinite(popt).all() or not width > 0:
        return np.nan, np.nan, np.nan, "flat edge fit gave values that are not finite"
    if not lower < upper:
        return np.nan, np.nan, np.nan, "flat edges are not ordered"
    return float(lower), float(upper), float(width), ""


def flat_slit_edges(flatData, slit, orderMap, binArcsec=PROFILE_BIN_ARCSEC):
    """*fit the lit slit edges of the master flat in every order*

    The flat is binned along the slit, normalised to its interior median, and fitted with an erf top-hat
    A·[Φ((s−l)/w) − Φ((s−u)/w)]. The flat's own quality mask is not applied, as it flags the low-level edge pixels
    this fit needs. Pixels that are not finite are skipped.

    The flat pixels used here need only a whole-numbered order and a finite slit position. The flat profile that
    ``skyline_slit_profiles`` uses also requires a positive finite wavelength, so the two pixel sets can differ.
    This is deliberate, as changing it could move the measured offsets.

    **Key Arguments:**

    - ``flatData`` -- the master flat's flux image
    - ``slit`` -- the 2D map's SLIT image in arcsec, the same shape as ``flatData``
    - ``orderMap`` -- the 2D map's ORDER image, the same shape as ``flatData``
    - ``binArcsec`` -- the slit bin width in arcsec. Default *PROFILE_BIN_ARCSEC*

    **Return:**

    - ``flatEdges`` -- dataframe with the columns ``order``, ``lower``, ``upper``, ``width`` (arcsec) and ``reason``.
      An order whose fit failed has NaN edges and a non-empty reason

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination import flat_slit_edges
    flatEdges = flat_slit_edges(masterFlat.data, slit, orderMap)
    ```
    """
    isInOrder = in_order_pixels(orderMap, slit)
    rows = []
    for order in np.unique(orderMap[isInOrder]):
        isOrder = isInOrder & (orderMap == order)
        centres, profile, reason = normalised_flat_profile(flatData, slit, isOrder, binArcsec)
        lower = upper = width = np.nan
        if not reason:
            lower, upper, width, reason = _fit_flat_edges(centres, profile, binArcsec)
        rows.append({"order": int(order), "lower": lower, "upper": upper, "width": width, "reason": reason})
    return pd.DataFrame(rows, columns=FLAT_EDGE_COLUMNS)


def _fit_shifted_top_hat(slit, profile, err, flatEdge, maxOffsetArcsec):
    """*fit A·T(s−δ; l, u, w) to a science profile from several δ starting points; return the best fit or None*

    The errors are used as absolute 1-sigma errors (``absolute_sigma=True``), so the returned covariance is the
    propagated error and is not rescaled by the scatter of the fit. ``_offset_error`` applies the scatter.

    **Return:**

    - ``best`` -- ``(chi2, popt, pcov)`` of the lowest chi-squared start, or *None* if no start converged
    """
    lower, upper, flatWidth = flatEdge
    limit = OFFSET_BOUND_FACTOR * maxOffsetArcsec

    def model(s, amplitude, delta, width):
        return top_hat(s, amplitude, lower + delta, upper + delta, width)

    bounds = ((0.0, -limit, MIN_EDGE_WIDTH_ARCSEC), (np.inf, limit, MAX_SCIENCE_EDGE_WIDTH_ARCSEC))
    startWidth = float(np.clip(flatWidth, MIN_EDGE_WIDTH_ARCSEC, MAX_SCIENCE_EDGE_WIDTH_ARCSEC))
    best = None
    for start in OFFSET_START_ARCSEC:
        p0 = (
            1.0,
            float(np.clip(start, -START_INSIDE_BOUND_FRACTION * limit, START_INSIDE_BOUND_FRACTION * limit)),
            startWidth,
        )
        try:
            popt, pcov = _quiet_curve_fit(model, slit, profile, p0=p0, sigma=err, absolute_sigma=True, bounds=bounds)
        except (RuntimeError, ValueError):
            continue
        chi2 = float(np.sum(((profile - model(slit, *popt)) / err) ** 2))
        if not np.isfinite(chi2):
            continue
        if best is None or chi2 < best[0]:
            best = (chi2, popt, pcov)
    return best


def _offset_error(pcov, chi2, nBins):
    """*the 1-sigma error on the offset: the propagated error, inflated by the fit's scatter when the fit is poor*

    The propagated error is ``sqrt(pcov[1, 1])`` for the quoted input errors. It is multiplied by
    ``sqrt(max(1, chi2 / (nBins - 3)))``, so it is a floor: a fit that is better than its quoted errors keeps the
    propagated error, and a fit that is worse (the top-hat does not describe the profile) has the error inflated.
    The error is *inf* if the covariance is not finite or negative.
    """
    variance = pcov[1, 1]
    if not (np.isfinite(variance) and variance >= 0):
        return np.inf
    reducedChi2 = chi2 / max(nBins - N_FIT_PARAMETERS, 1)
    return float(np.sqrt(variance) * np.sqrt(max(1.0, reducedChi2)))


def _flat_edge_problem(flatRow):
    """*why a flat-edge row cannot place the science edges, or an empty string*

    A usable row has finite ``lower``, ``upper`` and ``width``, ``lower < upper`` and ``width > 0``.
    """
    if flatRow is None:
        return "the order is missing"
    edges = np.array([flatRow["lower"], flatRow["upper"], flatRow["width"]], dtype=float)
    if np.isfinite(edges).all() and edges[0] < edges[1] and edges[2] > 0:
        return ""
    reason = flatRow["reason"]
    if isinstance(reason, str) and reason:
        return reason
    return "the flat edges are not finite, ordered and positive"


def _offset_failure_reason(offset, offsetErr, width, flatWidth, maxOffsetArcsec, maxOffsetErrArcsec):
    """*why a converged fit cannot be used, or an empty string*"""
    ratio = width / flatWidth
    if not MIN_WIDTH_RATIO <= ratio <= MAX_WIDTH_RATIO:
        return f"science edge width {width:.2f} arcsec is degenerate (flat edge width {flatWidth:.2f} arcsec)"
    if not np.isfinite(offsetErr):
        return "offset uncertainty is not finite"
    if offsetErr > maxOffsetErrArcsec:
        return f"offset uncertainty {offsetErr:.2f} arcsec exceeds the maximum {maxOffsetErrArcsec:.2f} arcsec"
    if abs(offset) > maxOffsetArcsec:
        return f"offset {offset:+.2f} arcsec exceeds the maximum {maxOffsetArcsec:.2f} arcsec"
    return ""


def _fit_one_order(orderProfile, flatRow, maxOffsetArcsec, maxOffsetErrArcsec, note=""):
    """*the offset-table row for one order*"""
    slit = orderProfile["slit"].to_numpy(dtype=float)
    profile = orderProfile["profile"].to_numpy(dtype=float)
    err = orderProfile["profileErr"].to_numpy(dtype=float)
    isValid = np.isfinite(slit) & np.isfinite(profile) & np.isfinite(err) & (err > 0)
    row = {
        "offset": np.nan,
        "offsetErr": np.nan,
        "scienceWidth": np.nan,
        "nBins": int(isValid.sum()),
        "status": STATUS_FAILED,
        "reason": "",
        "borrowedFrom": np.nan,
    }
    flatProblem = _flat_edge_problem(flatRow)
    if flatProblem:
        return {**row, "reason": f"no flat edges: {flatProblem}"}
    if row["nBins"] < MIN_PROFILE_BINS:
        return {**row, "reason": note or f"too few profile bins ({row['nBins']} < {MIN_PROFILE_BINS})"}
    flatEdge = (float(flatRow["lower"]), float(flatRow["upper"]), float(flatRow["width"]))
    best = _fit_shifted_top_hat(slit[isValid], profile[isValid], err[isValid], flatEdge, maxOffsetArcsec)
    if best is None:
        return {**row, "reason": "offset fit did not converge"}
    chi2, (_, offset, width), pcov = best
    offsetErr = _offset_error(pcov, chi2, row["nBins"])
    reason = _offset_failure_reason(offset, offsetErr, width, flatEdge[2], maxOffsetArcsec, maxOffsetErrArcsec)
    fitted = {**row, "offset": float(offset), "offsetErr": offsetErr, "scienceWidth": float(width), "reason": reason}
    return fitted if reason else {**fitted, "status": STATUS_FITTED}


def _borrow_from_partners(rows, partners):
    """*give each failed order the offset of its partner order if that partner was fitted on its own data*"""
    byOrder = {row["order"]: row for row in rows}
    borrowed = []
    for row in rows:
        partner = byOrder.get((partners or {}).get(row["order"]))
        if row["status"] != STATUS_FAILED or partner is None or partner["status"] != STATUS_FITTED:
            borrowed.append(row)
            continue
        reason = f"borrowed from order {partner['order']}; own fit failed: {row['reason']}"
        borrowed.append(
            {
                **row,
                "offset": partner["offset"],
                "offsetErr": partner["offsetErr"],
                "status": STATUS_BORROWED,
                "reason": reason,
                "borrowedFrom": partner["order"],
            }
        )
    return borrowed


def fit_order_offsets(
    profiles,
    flatEdges,
    maxOffsetArcsec=MAX_OFFSET_ARCSEC,
    maxOffsetErrArcsec=MAX_OFFSET_ERR_ARCSEC,
    partners=None,
    profileNotes=None,
):
    """*fit the science slit profile of every order as the flat's top-hat shifted by δ*

    The model is P(s) = A·T(s − δ; l_f, u_f, w), with the flat edges l_f and u_f fixed and A, δ and w free. Only the
    bins in ``profiles`` are used, so a science edge that lies outside the map footprint, or in a masked part of the
    frame, is simply not in the data. A positive δ means the science image sits toward +SLIT.

    An order fails, with a reason, when it has no flat edges, too few bins, no converging fit, a degenerate science
    edge width, too large an uncertainty, or an offset beyond ``maxOffsetArcsec``. A failed order takes the offset of
    its partner in ``partners`` if that partner was fitted on its own data, and is then marked as borrowed.

    **Key Arguments:**

    - ``profiles`` -- dataframe from ``skyline_slit_profiles`` with ``order``, ``slit``, ``profile`` and ``profileErr``
    - ``flatEdges`` -- dataframe from ``flat_slit_edges``
    - ``maxOffsetArcsec`` -- the largest accepted |δ| in arcsec. Default *MAX_OFFSET_ARCSEC*
    - ``maxOffsetErrArcsec`` -- the largest accepted uncertainty on δ in arcsec. Default *MAX_OFFSET_ERR_ARCSEC*
    - ``partners`` -- optional ``{order: partnerOrder}``, from which a failed order may borrow. Default *None*
    - ``profileNotes`` -- optional ``{order: reason}`` from ``skyline_slit_profiles``, used as the failure reason of
      an order that has too few bins. Default *None*

    **Return:**

    - ``offsets`` -- dataframe sorted by order with the columns ``order``, ``offset``, ``offsetErr``,
      ``scienceWidth``, ``nBins``, ``status`` (``fitted``, ``borrowed`` or ``failed``), ``reason`` and
      ``borrowedFrom``. A failed order has a NaN offset

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination import fit_order_offsets
    offsets = fit_order_offsets(profiles, flatEdges, partners={2: 3, 3: 2})
    ```
    """
    flatByOrder = {int(row["order"]): row for _, row in flatEdges.iterrows()}
    orders = sorted(set(flatByOrder) | {int(order) for order in profiles["order"].unique()})
    rows = []
    for order in orders:
        orderProfile = profiles.loc[profiles["order"] == order]
        note = (profileNotes or {}).get(order, "")
        fitted = _fit_one_order(orderProfile, flatByOrder.get(order), maxOffsetArcsec, maxOffsetErrArcsec, note)
        rows.append({"order": order, **fitted})
    return pd.DataFrame(_borrow_from_partners(rows, partners), columns=OFFSET_COLUMNS)


@dataclass(frozen=True, eq=False)
class slit_illumination_result:
    """The outcome of ``measure_slit_illumination``.

    - ``offsets`` -- per-order dataframe from ``fit_order_offsets``, empty if nothing could be measured
    - ``flatEdges`` -- per-order dataframe from ``flat_slit_edges``, empty if there was no flat
    - ``profiles`` -- the per-order, per-bin skyline slit profiles that were fitted
    - ``isApplied`` -- True only if at least one order has a usable offset (fitted or borrowed)
    - ``reason`` -- why ``isApplied`` is False, or an empty string
    - ``overlapMarginArcsec`` -- the margin, from the settings, for the overlap mask that uses these offsets

    The dataclass is frozen, so its fields cannot be rebound. The three dataframes are copies made when the result is
    built, so an in-place edit of the dataframes that were passed in does not reach the result. The dataframes held by
    the result are still mutable: the result is not deeply immutable.
    """

    offsets: pd.DataFrame
    flatEdges: pd.DataFrame
    profiles: pd.DataFrame
    isApplied: bool
    reason: str
    overlapMarginArcsec: float

    def __post_init__(self):
        """*replace the three dataframes with copies, so the result does not alias the caller's frames*"""
        for name in ("offsets", "flatEdges", "profiles"):
            object.__setattr__(self, name, getattr(self, name).copy())


# SETTINGS KEY -> (CONFIG NAME, DEFAULT, MAY BE ZERO)
ILLUMINATION_SETTINGS = {
    "profile-bin-arcsec": ("binArcsec", PROFILE_BIN_ARCSEC, False),
    "max-offset-arcsec": ("maxOffsetArcsec", MAX_OFFSET_ARCSEC, False),
    "max-order-scatter-arcsec": ("maxOffsetErrArcsec", MAX_OFFSET_ERR_ARCSEC, False),
    "overlap-margin-arcsec": ("overlapMarginArcsec", OVERLAP_MARGIN_ARCSEC, True),
}


def _resolve_settings(log, illuminationSettings):
    """*the measurement settings: each ``slit-illumination`` key if it holds a usable number, else its default*"""
    config = {name: default for name, default, _ in ILLUMINATION_SETTINGS.values()}
    for key, (name, default, mayBeZero) in ILLUMINATION_SETTINGS.items():
        if key not in (illuminationSettings or {}):
            continue
        value = illuminationSettings[key]
        isNumber = isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)
        if isNumber and (value > 0 or (mayBeZero and value == 0)):
            config[name] = float(value)
        else:
            log.warning(f"slit-illumination setting {key} = {value!r} is not usable, so the default {default} is used")
    return config


def _read_disp_map(twoDMapPath, binx, biny, kw):
    """*the WAVELENGTH, SLIT and ORDER images of a 2D map, block-reduced from its own binning to the frame's*

    ``kw`` is the fits keyword lookup used to read the map's own binning.
    """
    from soxspipe.commonutils.toolkit import _open_disp_map_hdus

    hdul, _, _ = _open_disp_map_hdus(twoDMapPath, binx, biny, kw=kw)
    try:
        return tuple(np.array(hdul[name].data, dtype=float) for name in ("WAVELENGTH", "SLIT", "ORDER"))
    finally:
        hdul.close()


def _input_problem(frame, masterFlat, mapImages=None):
    """*why the frame, flat and (if given) map images cannot be measured, or an empty string*

    The frame and flat are checked first, so that they can be checked before the map file is opened. ``mapImages``
    is the ``(wavelength, slit, orderMap)`` tuple, whose shapes must match the frame and the flat.
    """
    if masterFlat is None or masterFlat.data is None:
        return "no master flat"
    if frame.uncertainty is None:
        return "the frame has no uncertainty"
    if frame.uncertainty.uncertainty_type != "std":
        return f"the frame uncertainty is {frame.uncertainty.uncertainty_type}, not std"
    if frame.data.shape != masterFlat.data.shape:
        return f"image shapes differ: frame {frame.data.shape}, flat {masterFlat.data.shape}"
    maskShape = np.shape(frame.mask)
    if maskShape not in ((), frame.data.shape):
        return f"image shapes differ: frame {frame.data.shape}, frame mask {maskShape}"
    if mapImages is None:
        return ""
    for name, image in zip(("WAVELENGTH", "SLIT", "ORDER"), mapImages, strict=True):
        if image.shape != frame.data.shape:
            return f"image shapes differ: frame {frame.data.shape}, 2D map {name} {image.shape}"
    return ""


def _usable_pixels(frame):
    """*True where the frame has finite data and a positive finite uncertainty, and is not masked*"""
    errs = np.asarray(frame.uncertainty.array, dtype=float)
    isMasked = np.zeros(frame.data.shape, dtype=bool) if frame.mask is None else np.asarray(frame.mask, dtype=bool)
    return ~isMasked & np.isfinite(frame.data) & np.isfinite(errs) & (errs > 0)


def _unusable_reason(offsets):
    """*one line naming why each order has no usable offset*"""
    reasons = [f"{row['order']}: {row['reason']}" for _, row in offsets.iterrows()]
    return "no order has a usable offset (" + "; ".join(reasons) + ")" if reasons else "the 2D map holds no orders"


def _not_applied(log, reason, config, flatEdges=None, profiles=None, offsets=None):
    """*the result for a measurement that cannot be applied; the reason is logged*"""
    log.warning(f"slit illumination offset is not applied: {reason}")
    return slit_illumination_result(
        offsets=pd.DataFrame(columns=OFFSET_COLUMNS) if offsets is None else offsets,
        flatEdges=pd.DataFrame(columns=FLAT_EDGE_COLUMNS) if flatEdges is None else flatEdges,
        profiles=pd.DataFrame(columns=PROFILE_COLUMNS) if profiles is None else profiles,
        isApplied=False,
        reason=reason,
        overlapMarginArcsec=config["overlapMarginArcsec"],
    )


def _measure_orders(log, frame, masterFlat, mapImages, settings, arm, config, partners):
    """*the flat edges, skyline profiles and offsets of every order, for inputs that passed the checks*"""
    from soxspipe.commonutils.toolkit import get_skylines_dataframe

    wavelength, slit, orderMap = mapImages
    flatData = np.asarray(masterFlat.data, dtype=float)
    flatEdges = flat_slit_edges(flatData, slit, orderMap, binArcsec=config["binArcsec"])
    skylinesDF = get_skylines_dataframe(log, settings, arm, minBrightnessVIS=SKYLINE_MIN_FLUX_VIS)
    profiles, profileNotes = skyline_slit_profiles(
        np.asarray(frame.data, dtype=float),
        np.asarray(frame.uncertainty.array, dtype=float),
        _usable_pixels(frame),
        flatData,
        wavelength,
        slit,
        orderMap,
        skylinesDF["WAVELENGTH"].to_numpy(dtype=float),
        binArcsec=config["binArcsec"],
    )
    offsets = fit_order_offsets(
        profiles,
        flatEdges,
        config["maxOffsetArcsec"],
        config["maxOffsetErrArcsec"],
        partners=partners,
        profileNotes=profileNotes,
    )
    return flatEdges, profiles, offsets


def measure_slit_illumination(
    log,
    frame,
    masterFlat,
    twoDMapPath,
    settings,
    arm,
    binx,
    biny,
    illuminationSettings=None,
    partners=None,
):
    """*measure how far each order's sky image sits along the slit from the master flat's*

    The offsets are measured order by order and are never combined, as they differ in size and sign between orders.
    A positive offset means the science image sits toward +SLIT. Data-quality failures never raise: they are logged,
    and the result then has ``isApplied`` False with a ``reason``, so the reduction can carry on as before.

    The frame and flat are checked first (a flat is present, the frame has a StdDev uncertainty, their shapes match).
    Then the 2D map is read. A map that cannot be read (a missing file, a missing extension or an unreadable image,
    caught as ``OSError``, ``KeyError``, ``ValueError`` or ``AttributeError``) gives a not-applied result, as does a
    WAVELENGTH, SLIT or ORDER image whose shape differs from the frame's. Other exceptions are not caught, as they
    would be programming errors.

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- the detrended, flat-fielded science CCDData, with a mask and a StdDev uncertainty
    - ``masterFlat`` -- the master flat CCDData, or *None* if there is none
    - ``twoDMapPath`` -- path to the 2D dispersion map
    - ``settings`` -- the soxspipe settings dictionary, used to find the skyline table and the FITS keywords
    - ``arm`` -- the spectrograph arm
    - ``binx`` -- binning in the x-axis
    - ``biny`` -- binning in the y-axis
    - ``illuminationSettings`` -- the ``slit-illumination`` settings block, with the keys ``profile-bin-arcsec``,
      ``max-offset-arcsec``, ``max-order-scatter-arcsec`` (the largest accepted uncertainty on an offset) and
      ``overlap-margin-arcsec``. A missing key means the module constant. Default *None*
    - ``partners`` -- optional ``{order: partnerOrder}`` from which a failed order may borrow an offset. Default *None*

    **Return:**

    - ``result`` -- a frozen ``slit_illumination_result`` holding copies of its dataframes

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination import measure_slit_illumination
    result = measure_slit_illumination(log, frame, masterFlat, twoDMapPath, settings, arm, binx, biny)
    if result.isApplied:
        offsets = result.offsets
    ```
    """
    log.debug("starting the ``measure_slit_illumination`` function")
    config = _resolve_settings(log, illuminationSettings)
    problem = _input_problem(frame, masterFlat)
    if problem:
        return _not_applied(log, problem, config)
    from soxspipe.commonutils.keyword_lookup import keyword_lookup

    kw = keyword_lookup(log=log, settings=settings).get
    try:
        mapImages = _read_disp_map(twoDMapPath, binx, biny, kw)
    # ATTRIBUTEERROR COVERS A MAP EXTENSION THAT IS PRESENT BUT HOLDS NO DATA
    except (OSError, KeyError, ValueError, AttributeError) as error:
        return _not_applied(log, f"the 2D map could not be read ({type(error).__name__}: {error})", config)
    problem = _input_problem(frame, masterFlat, mapImages)
    if problem:
        return _not_applied(log, problem, config)
    flatEdges, profiles, offsets = _measure_orders(log, frame, masterFlat, mapImages, settings, arm, config, partners)
    for _, row in offsets.loc[offsets["status"] != STATUS_FITTED].iterrows():
        log.info(f"slit illumination offset of order {row['order']} is {row['status']}: {row['reason']}")
    isApplied = bool((offsets["status"] != STATUS_FAILED).any())
    if not isApplied:
        return _not_applied(log, _unusable_reason(offsets), config, flatEdges, profiles, offsets)
    log.debug("completed the ``measure_slit_illumination`` function")
    return slit_illumination_result(offsets, flatEdges, profiles, True, "", config["overlapMarginArcsec"])
