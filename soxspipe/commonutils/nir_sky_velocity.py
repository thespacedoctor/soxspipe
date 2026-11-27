"""The NIR sky correction as one velocity shift across the orders (DY-2477).

The sky-line fit of ``nir_sky_refit`` models the offset of a frame as a detector motion in px, which each row turns into
nm with its own pixel scale. Two orders that overlap in wavelength have different pixel scales there, so the same
wavelength got a different correction in each. A frame offset is a velocity (Dave, 2026-11-27): the correction of a
row is ``wavelength * v(wavelength)``, and ``v`` is one smooth polynomial in wavelength for all the orders that take
part. The orders named by the ``own-correction-orders`` setting are the exception and keep the pixel model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from soxspipe.commonutils import nir_sky_refit
from soxspipe.commonutils.nir_sky_refit import (
    MATCH_ITERATIONS,
    MAX_GLOBAL_RMS_PX,
    MIN_GLOBAL_LINES,
    _empty_lines,
    _fallback,
    _is_usable_geometry,
    _match_order_lines,
    _robust_fit,
    _usable_lines,
    order_sky_correction,
    seed_matches,
    sky_refit_result,
)
from soxspipe.commonutils.ug_wavelength_transfer import detector_px_to_nm_shift, nm_shift_to_detector_px


def detector_px_to_velocity(shiftPx, dispersionNmPerPx, wavelengthNm):
    """*the velocity (``Δλ/λ``, unitless) that moves a wavelength as a detector motion does*

    A detector motion of ``shiftPx`` corrects the wavelength by ``-dispersion * shiftPx`` nm
    (``detector_px_to_nm_shift``).
    Over the wavelength that is a velocity: positive when the correction adds to the wavelength.

    **Key Arguments:**

    - ``shiftPx`` -- the detector motion along the dispersion axis (px)
    - ``dispersionNmPerPx`` -- the signed dispersion d(wavelength)/d(pixel) of the order (nm/px)
    - ``wavelengthNm`` -- the wavelength (nm) at which the motion is taken

    **Return:**

    - ``velocity`` -- the fractional wavelength shift ``Δλ/λ``

    **Raises:**

    - ``ValueError`` -- if the dispersion is zero or not finite
    """
    return detector_px_to_nm_shift(shiftPx, dispersionNmPerPx) / wavelengthNm


def velocity_to_detector_px(velocity, dispersionNmPerPx, wavelengthNm):
    """*the detector motion (px) that moves a wavelength by a velocity; the inverse of ``detector_px_to_velocity``*

    **Key Arguments:**

    - ``velocity`` -- the fractional wavelength shift ``Δλ/λ``
    - ``dispersionNmPerPx`` -- the signed dispersion d(wavelength)/d(pixel) of the order (nm/px)
    - ``wavelengthNm`` -- the wavelength (nm) at which the velocity is taken

    **Return:**

    - ``shiftPx`` -- the detector motion along the dispersion axis (px)

    **Raises:**

    - ``ValueError`` -- if the dispersion is zero or not finite
    """
    return nm_shift_to_detector_px(velocity * wavelengthNm, dispersionNmPerPx)


@dataclass(frozen=True)
class sky_velocity:
    """*one velocity (``Δλ/λ``) that is a polynomial in wavelength: ``sum_k c_k x^k`` with ``x = (λ - centre) / half``*

    **Key Arguments:**

    - ``centreNm`` -- the wavelength (nm) where ``x`` is 0
    - ``halfRangeNm`` -- the wavelength (nm) where ``x`` is 1, away from the centre
    - ``coefficients`` -- ``c_0, c_1, ...`` (unitless)

    **Raises:**

    - ``ValueError`` -- if the half-range is not a finite number above 0, or there are no coefficients
    """

    centreNm: float
    halfRangeNm: float
    coefficients: tuple[float, ...]

    def __post_init__(self) -> None:
        if not (math.isfinite(self.halfRangeNm) and self.halfRangeNm > 0) or len(self.coefficients) == 0:
            raise ValueError(
                "a sky velocity needs a finite half-range above 0 nm and at least one coefficient, "
                f"got {self.halfRangeNm!r} and {self.coefficients!r}"
            )

    def fraction_at(self, wavelengthNm) -> np.ndarray:
        """*the velocity ``Δλ/λ`` at each wavelength*

        **Key Arguments:**

        - ``wavelengthNm`` -- the wavelengths (nm)

        **Return:**

        - ``velocity`` -- the velocity at each wavelength, unitless
        """
        x = (np.asarray(wavelengthNm, dtype=float) - self.centreNm) / self.halfRangeNm
        return np.polynomial.polynomial.polyval(x, self.coefficients)


def velocity_correction(
    order, velocity, geometry, nLines, globalSeedLines=None, ownSeedLines=None
) -> order_sky_correction:
    """*the ``order_sky_correction`` of an order that takes the velocity ``velocity``*

    ``shiftPx`` and ``alongPx`` are the detector motion at the order's centre and its linear change along the order that
    the velocity amounts to, at the order's median dispersion. They feed the ``SKY SHIFT O<n> PX`` QC and the logs, and
    do not take part in the correction.

    **Key Arguments:**

    - ``order`` -- the order
    - ``velocity`` -- the ``sky_velocity``
    - ``geometry`` -- the order's ``order_geometry``
    - ``nLines`` -- the lines of this order kept by the fit
    - ``globalSeedLines`` -- see ``order_sky_correction``. Default *None*
    - ``ownSeedLines`` -- see ``order_sky_correction``. Default *None*

    **Return:**

    - ``correction`` -- the ``order_sky_correction``, whose ``velocity`` is set
    """
    ends = geometry.centreNm + np.array([-1.0, 0.0, 1.0]) * geometry.halfRangeNm
    shiftsNm = ends * velocity.fraction_at(ends)
    shiftsPx = nm_shift_to_detector_px(shiftsNm, geometry.dispersionNmPerPx)
    return order_sky_correction(
        order=int(order),
        shiftPx=float(shiftsPx[1]),
        alongPx=(float((shiftsPx[2] - shiftsPx[0]) / 2),),
        perOrderSlopePx=0.0,
        nLines=int(nLines),
        geometry=geometry,
        globalSeedLines=globalSeedLines,
        ownSeedLines=ownSeedLines,
        velocity=velocity,
    )


def _velocity_range(geometry) -> tuple[float, float]:
    """*the centre and half-range (nm) of the wavelengths that the orders of ``geometry`` cover*

    **Key Arguments:**

    - ``geometry`` -- ``{order: order_geometry}``

    **Return:**

    - ``(centreNm, halfRangeNm)`` -- the middle and half the width of the range
    """
    low = min(g.centreNm - g.halfRangeNm for g in geometry.values())
    high = max(g.centreNm + g.halfRangeNm for g in geometry.values())
    return (high + low) / 2, (high - low) / 2


def _px_weights(lines, geometry) -> np.ndarray:
    """*the factor that turns a line's offset in nm into detector motion (px) on its own order's signed pixel scale*

    **Key Arguments:**

    - ``lines`` -- the matched lines, with ``order`` and optionally ``pixelScaleNm``
    - ``geometry`` -- ``{order: order_geometry}``

    **Return:**

    - ``weights`` -- one factor (px/nm) per line
    """
    orders = lines["order"].to_numpy()
    if "pixelScaleNm" in lines.columns:
        scales = lines["pixelScaleNm"].to_numpy(dtype=float)
    else:
        scales = np.array([geometry[o].medianPixelScaleNm for o in orders])
    signs = np.array([np.sign(geometry[o].dispersionNmPerPx) for o in orders])
    return -signs / scales


def fit_velocity_model(lines, geometry, settings) -> sky_refit_result:
    """*fit one velocity polynomial in wavelength to the matched sky lines of the orders in ``geometry``*

    A line at ``observedNm`` is ``catalogueNm - observedNm`` from where it belongs, which the model takes as
    ``observedNm * v(observedNm)``. The fit runs in px, each line's offset divided by its own pixel scale, so the
    clipping limit and the scatter limit are the same as for the pixel model, and a line weighs as its centroid
    precision does. Every order with a usable geometry gets a correction, an order with no lines taking the velocity
    at its wavelengths.

    **Key Arguments:**

    - ``lines`` -- one row per matched line: ``order``, ``catalogueNm``, ``observedNm`` (nm) and optionally
      ``pixelScaleNm`` at the line. Rows with a non-finite value or a non-positive pixel scale are dropped, and so are
      the rows of an order whose geometry is not usable
    - ``geometry`` -- ``{order: order_geometry}`` of every order to correct
    - ``settings`` -- the ``sky_tuning_settings``: ``velocityDeg`` and ``clipSigma`` are used

    **Return:**

    - ``result`` -- the ``sky_refit_result``, with ``velocity`` and ``coefficients`` of the velocity polynomial. Its
      ``fallbackReason`` is set, and its corrections are empty, when there are too few lines, the model is
      rank-deficient or the lines scatter more than ``MAX_GLOBAL_RMS_PX``

    **Raises:**

    - ``ValueError`` -- if a line belongs to an order that has no geometry
    """
    lines = lines.reset_index(drop=True)
    unknown = sorted(set(lines["order"]) - set(geometry))
    if unknown:
        raise ValueError(f"order {unknown[0]} has sky lines but no geometry")
    geometry = {o: g for o, g in geometry.items() if _is_usable_geometry(g)}
    lines = _usable_lines(lines, geometry)
    if len(lines) < MIN_GLOBAL_LINES:
        return _fallback(lines, f"{len(lines)} sky lines matched (min {MIN_GLOBAL_LINES})")

    centreNm, halfRangeNm = _velocity_range(geometry)
    observed = lines["observedNm"].to_numpy(dtype=float)
    x = (observed - centreNm) / halfRangeNm
    weights = _px_weights(lines, geometry)
    design = np.column_stack([weights * observed * x**k for k in range(settings.velocityDeg + 1)])
    residualPx = weights * (lines["catalogueNm"].to_numpy(dtype=float) - observed)

    fit = _robust_fit(design, residualPx, settings.clipSigma)
    if fit is None:
        return _fallback(lines, "the velocity model is rank-deficient (too few distinct wavelengths)")
    coefficients, kept = fit
    modelPx = design @ coefficients
    rmsPx = float(np.sqrt(np.mean((residualPx - modelPx)[kept] ** 2)))
    clipFraction = float(1.0 - kept.mean())
    if rmsPx > MAX_GLOBAL_RMS_PX:
        return _fallback(lines, f"fit rms {rmsPx:0.2f} px (max {MAX_GLOBAL_RMS_PX})", rmsPx, clipFraction)

    velocity = sky_velocity(centreNm, halfRangeNm, tuple(float(c) for c in coefficients))
    lines = lines.assign(residualPx=residualPx, modelPx=modelPx, clipped=~kept)
    counts = lines.loc[~lines["clipped"], "order"].value_counts()
    corrections = {
        int(o): velocity_correction(o, velocity, g, nLines=int(counts.get(o, 0))) for o, g in geometry.items()
    }
    return sky_refit_result(
        corrections, velocity.coefficients, int(kept.sum()), rmsPx, clipFraction, lines, None, velocity
    )


def _fit_velocity_with_iterated_matching(orderLines, geometry, seedPx, settings, iterations, ownSeedsPx):
    """*match the lines of the velocity orders and fit their velocity, matching again about each fit*

    The first match moves the peaks of each order by the rigid seed, or by the order's own seed when that matches more
    of its lines (``seed_matches``). Every order then goes into the fit, with no exception for an order that matched
    few lines about its own seed: a line that the velocity does not follow is clipped. Each later match moves the peaks
    by the previous fit, which brings in the lines of an order that the rigid seed missed.

    **Key Arguments:**

    - ``orderLines`` -- ``{order: nir_order_sky_lines}`` of the velocity orders
    - ``geometry`` -- ``{order: order_geometry}`` of the velocity orders
    - ``seedPx`` -- the rigid detector shift (px) of the first match
    - ``settings`` -- the ``sky_tuning_settings``
    - ``iterations`` -- the number of match-and-fit rounds
    - ``ownSeedsPx`` -- ``{order: shift (px)}`` that each order measured for itself

    **Return:**

    - ``result`` -- the ``sky_refit_result`` of the last round, or of the first one that fell back
    """
    orderLines = {o: lines for o, lines in orderLines.items() if _is_usable_geometry(geometry[o])}
    matches = seed_matches(orderLines, geometry, seedPx, ownSeedsPx)
    frames = [match.matched for match in matches.values()]
    result = None
    for _ in range(iterations):
        result = fit_velocity_model(
            pd.concat(frames, ignore_index=True) if frames else _empty_lines(), geometry, settings
        )
        if result.fallbackReason is not None:
            return result
        frames = [_match_order_lines(lines, result.corrections[o]) for o, lines in orderLines.items()]
    corrections = {
        o: replace(c, globalSeedLines=matches[o].globalLines, ownSeedLines=matches[o].ownLines) if o in matches else c
        for o, c in result.corrections.items()
    }
    return replace(result, corrections=corrections)


def _joined_result(velocity, pixel, ownOrders) -> sky_refit_result:
    """*the velocity result with the pixel corrections of ``ownOrders`` added, counting the lines of both fits*

    **Key Arguments:**

    - ``velocity`` -- the ``sky_refit_result`` of the velocity fit
    - ``pixel`` -- the ``sky_refit_result`` of the pixel fit
    - ``ownOrders`` -- the orders that keep their pixel correction

    **Return:**

    - ``result`` -- the joined ``sky_refit_result``
    """
    ownCorrections = {o: c for o, c in pixel.corrections.items() if o in ownOrders}
    ownLines = pixel.lines[pixel.lines["order"].isin(ownOrders)]
    lines = pd.concat([velocity.lines, ownLines], ignore_index=True)
    kept = ~lines["clipped"].to_numpy(dtype=bool)
    scatterPx = (lines["residualPx"] - lines["modelPx"]).to_numpy(dtype=float)
    return replace(
        velocity,
        corrections={**velocity.corrections, **ownCorrections},
        nLines=int(kept.sum()),
        rmsPx=float(np.sqrt(np.mean(scatterPx[kept] ** 2))),
        clipFraction=float(1.0 - kept.mean()),
        lines=lines,
    )


def fit_nir_sky_corrections(
    orderLines, geometry, seedPx, settings, iterations=MATCH_ITERATIONS, ownSeedsPx=None
) -> sky_refit_result:
    """*match the sky lines and fit the NIR correction: one velocity, and pixel models for the own-correction orders*

    The orders named by ``settings.ownCorrectionOrders`` keep the correction of ``fit_with_iterated_matching``: the
    pixel model fitted to the lines of every order, with the order's own offset and slope, unchanged. Every other order
    gets the velocity polynomial that ``fit_velocity_model`` fits to the lines of those orders alone, so the same
    wavelength has the same correction in two overlapping orders. A fallback of either fit is the fallback of the
    whole, and the pixel fit only runs when an own-correction order is in the frame.

    **Key Arguments:**

    - ``orderLines`` -- ``{order: nir_order_sky_lines}`` of the orders that may contribute lines
    - ``geometry`` -- ``{order: order_geometry}`` of every order to correct
    - ``seedPx`` -- the rigid detector shift (px) of the first match, or *None* when no order could give one
    - ``settings`` -- the ``sky_tuning_settings``
    - ``iterations`` -- the number of match-and-fit rounds. Default ``MATCH_ITERATIONS``
    - ``ownSeedsPx`` -- ``{order: shift (px)}`` that each order measured for itself (``own_seeds_px``). Default *None*

    **Return:**

    - ``result`` -- the ``sky_refit_result``. Its ``velocity`` and ``coefficients`` are those of the velocity fit, its
      ``nLines``, ``rmsPx``, ``clipFraction`` and ``lines`` count the lines of both fits, and its corrections hold every
      order with a geometry

    **Raises:**

    - ``ValueError`` -- if an order with lines has no geometry
    """
    if seedPx is None:
        return _fallback(_empty_lines(), "no order gave a seed shift")
    ownOrders = set(settings.ownCorrectionOrders) & set(geometry)
    ownSeedsPx = {} if ownSeedsPx is None else ownSeedsPx
    pixel = None
    if ownOrders:
        pixel = nir_sky_refit.fit_with_iterated_matching(
            orderLines, geometry, seedPx, settings, iterations=iterations, ownSeedsPx=ownSeedsPx
        )
        if pixel.fallbackReason is not None:
            return pixel
    missing = sorted(set(orderLines) - set(geometry))
    if missing:
        raise ValueError(f"order {missing[0]} has sky lines but no geometry")
    velocityGeometry = {o: g for o, g in geometry.items() if o not in ownOrders}
    velocityLines = {o: lines for o, lines in orderLines.items() if o in velocityGeometry}
    velocity = _fit_velocity_with_iterated_matching(
        velocityLines, velocityGeometry, seedPx, settings, iterations, ownSeedsPx
    )
    if velocity.fallbackReason is not None or pixel is None:
        return velocity
    return _joined_result(velocity, pixel, ownOrders)
