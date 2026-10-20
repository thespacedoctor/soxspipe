"""Predict the SOXS VIS u and g order wavelength shifts, which have no usable sky lines, from the r and i sky shifts.

The r order has only 3-4 usable sky lines, so its shift is also predicted from the i shift and the rotator angle change
(DY-1645):

    dy_r = i * di + cos * (cos(th_sci) - cos(th_cal)) + sin * (sin(th_sci) - sin(th_cal))

The u (order 2) and g (order 3) orders move against a daytime arc by up to 5-6 px at night. That motion is predicted
from the mean sky-measured detector shift of the r (order 1) and i (order 4) orders plus the change in the derotator
angle (``ESO ADA ABSROT``) between the science frame and the calibration frame (DY-1608):

    dy_order = m * (dr + di) / 2 + cos * (cos(th_sci) - cos(th_cal)) + sin * (sin(th_sci) - sin(th_cal))

All shifts are in detector px along the dispersion axis, science minus calibration, positive when the features move
to larger pixel coordinates. The coefficients live in the ``ug-wavelength-transfer`` block of the VIS detector
parameters.
"""

import math

from soxspipe.commonutils.slit_illumination_qc import (
    CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES,
    read_rotator_angle,
    rotator_angle_spread,
)

VIS_R_ORDER = 1
VIS_U_ORDER = 2
VIS_G_ORDER = 3
VIS_I_ORDER = 4
PREDICTED_ORDERS = (VIS_U_ORDER, VIS_G_ORDER)
MEASURED_ORDERS = (VIS_R_ORDER, VIS_I_ORDER)
# MINIMUM SKY LINES KEPT AFTER CLIPPING FOR AN r/i SHIFT TO DRIVE THE MODEL. THE HUB RAW-PIXEL STUDY USED 4 FOR r, BUT
# THE PIPELINE MATCHES ONLY ISOLATED CATALOGUE LINES AND KEEPS 2-3 IN r ON EVERY WELL-ATTACHED 2026 JUNE/JULY FRAME;
# THE ONE r SHIFT FROM A SINGLE LINE IN THAT SET WAS A +7 px OUTLIER (DY-1608 REAL-DATA CHECK)
MIN_MATCHED_LINES = {VIS_R_ORDER: 2, VIS_I_ORDER: 10}
# MINIMUM COMBINED PHOTON-NOISE SIGNAL-TO-NOISE OF AN ORDER'S CATALOGUE SKY LINES (``sky_line_snr``) FOR ITS OWN SKY
# SHIFT TO BE APPLIED OR TO DRIVE A MODEL (DY-1663). ON 25 JUNE/JULY 2026 VIS STARE STANDARDS THE ONE i SHIFT THAT
# TELLURICS SHOW TO BE NOISE (20260626T044348) HAS 4.7; THE LOWEST VERIFIED-GOOD i SHIFTS ARE 6.2 AND 6.3 (5-ARCSEC
# SLIT, r MEASURED ON ITS OWN AGREES) AND 8.8 (1-ARCSEC SLIT, ATTACHED ARC). 5.4 SITS ABOUT 15 PERCENT FROM EACH SIDE
MIN_SKY_LINE_SNR = 5.4
# MINIMUM SIGNAL-TO-NOISE OF THE SECOND-BRIGHTEST CATALOGUE LINE FOR A VIS r SHIFT MEASURED WITHOUT THE i PREDICTION
# (DY-1663). THE 557.7 nm LINE ALONE CARRIES THE COMBINED r SNR PAST ``MIN_SKY_LINE_SNR``, SO TWO KEPT LINES CAN BE ONE
# REAL LINE AND ONE NOISE PEAK. 20260626T044348 HAS 1.8 THERE AND A 1.2 px WRONG r SHIFT; THE OTHER r-ALONE FRAMES OF
# THE SET HAVE 3.4 OR MORE. 2.5 IS THEIR GEOMETRIC MIDDLE (SQRT(1.8 * 3.4) = 2.47)
MIN_R_SECOND_LINE_SNR = 2.5
# MINIMUM SKY LINES KEPT AFTER CLIPPING FOR ANY VIS ORDER SHIFT TO BE APPLIED (DY-1645)
MIN_SHIFT_LINES = 2
# FLOOR ON THE MEASURED r SHIFT UNCERTAINTY, SO A CHANCE-IDENTICAL SET OF LINE SHIFTS CANNOT GIVE IT INFINITE WEIGHT
MIN_MEASURED_SIGMA_PX = 0.05
# AN ARC IS ATTACHED TO A SCIENCE FRAME WHEN IT IS TAKEN WITHIN THIS MANY MINUTES (START TO START, EITHER SIDE) AND WITH
# A ROTATOR ANGLE WITHIN THIS MANY DEGREES (WRAP-AWARE). THESE ARE THE WELL-ATTACHED TEST OF THE HUB STUDY
# ``output/arc_transfer/report.md`` (``MAX_MINUTES = 7`` WITH ``<=``, ``MAX_DROT_DEG = 5`` WITH ``<``); DY-1665
MAX_ATTACHED_ARC_MINUTES = 7.0
MAX_ATTACHED_ARC_ROTATOR_DEGREES = 5.0
MINUTES_PER_DAY = 1440.0
# MJD-OBS IS ~61000 DAYS, SO ITS FLOAT RESOLUTION IS ~1E-6 MINUTES; ROUNDING KEEPS A GAP OF EXACTLY THE LIMIT ON IT
MINUTES_DECIMALS = 6
COEFFICIENTS_KEY = "ug-wavelength-transfer"
MODEL_TERMS = ("m", "cos", "sin")
R_COEFFICIENTS_KEY = "r-from-i-transfer"
R_MODEL_TERMS = ("i", "cos", "sin", "rms_px")


def load_ug_transfer_coefficients(detectorParams):
    """*read the u/g transfer coefficients from the arm's detector parameters*

    **Key Arguments:**

    - ``detectorParams`` -- the arm's detector parameters (the ``VIS`` block of ``soxs_detector_parameters.yaml``)

    **Return:**

    - ``coeffs`` -- ``{order: {"m": ..., "cos": ..., "sin": ...}}`` for orders 2 and 3, or *None* if the arm has no
      ``ug-wavelength-transfer`` block

    A block that is present but lacks an order or a term, or holds a non-numeric or non-finite value, raises
    ``ValueError``.

    **Usage:**

    ```python
    from soxspipe.commonutils.ug_wavelength_transfer import load_ug_transfer_coefficients
    coeffs = load_ug_transfer_coefficients(self.detectorParams)
    ```
    """
    block = detectorParams.get(COEFFICIENTS_KEY)
    if block is None:
        return None
    orders = block.get("orders") or {}
    coeffs = {}
    for order in PREDICTED_ORDERS:
        try:
            coeffs[order] = {term: float(orders[order][term]) for term in MODEL_TERMS}
        except (KeyError, TypeError, ValueError) as e:
            raise ValueError(
                f"the {COEFFICIENTS_KEY} detector parameters need numeric {MODEL_TERMS} terms for order {order} ({e})"
            ) from e
        if not all(math.isfinite(value) for value in coeffs[order].values()):
            raise ValueError(f"the {COEFFICIENTS_KEY} detector parameters hold a non-finite term for order {order}")
    return coeffs


def predict_ug_detector_shifts(deltaRpx, deltaIpx, scienceAngleDeg, calibrationAngleDeg, coeffs):
    """*predict the u and g detector shifts from the r and i shifts and the two rotator angles*

    **Key Arguments:**

    - ``deltaRpx`` -- the measured r-order detector shift (px, science minus calibration)
    - ``deltaIpx`` -- the measured i-order detector shift (px, science minus calibration)
    - ``scienceAngleDeg`` -- the science-frame rotator angle (degrees)
    - ``calibrationAngleDeg`` -- the calibration-frame rotator angle (degrees)
    - ``coeffs`` -- per-order coefficients, ``{order: {"m": ..., "cos": ..., "sin": ...}}`` for orders 2 and 3

    **Return:**

    - ``shifts`` -- ``{2: dy_u, 3: dy_g}``, the predicted detector shifts in px

    **Usage:**

    ```python
    from soxspipe.commonutils.ug_wavelength_transfer import predict_ug_detector_shifts
    shifts = predict_ug_detector_shifts(0.3, 0.5, -19.8, 0.0, coeffs)
    ```
    """
    meanShift = (deltaRpx + deltaIpx) / 2
    deltaCos, deltaSin = _rotator_terms(scienceAngleDeg, calibrationAngleDeg)
    return {
        order: coeffs[order]["m"] * meanShift + coeffs[order]["cos"] * deltaCos + coeffs[order]["sin"] * deltaSin
        for order in PREDICTED_ORDERS
    }


def load_r_transfer_coefficients(detectorParams):
    """*read the r-from-i transfer coefficients from the arm's detector parameters*

    **Key Arguments:**

    - ``detectorParams`` -- the arm's detector parameters (the ``VIS`` block of ``soxs_detector_parameters.yaml``)

    **Return:**

    - ``coeffs`` -- ``{"i": ..., "cos": ..., "sin": ..., "rms_px": ...}``, or *None* if the arm has no
      ``r-from-i-transfer`` block

    A block that is present but lacks a term, or holds a non-numeric, non-finite or (for ``rms_px``) non-positive
    value, raises ``ValueError``.

    **Usage:**

    ```python
    from soxspipe.commonutils.ug_wavelength_transfer import load_r_transfer_coefficients
    coeffs = load_r_transfer_coefficients(self.detectorParams)
    ```
    """
    block = detectorParams.get(R_COEFFICIENTS_KEY)
    if block is None:
        return None
    try:
        coeffs = {term: float(block[term]) for term in R_MODEL_TERMS}
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError(
            f"the {R_COEFFICIENTS_KEY} detector parameters need numeric {R_MODEL_TERMS} terms ({e})"
        ) from e
    if not all(math.isfinite(value) for value in coeffs.values()):
        raise ValueError(f"the {R_COEFFICIENTS_KEY} detector parameters hold a non-finite term")
    if coeffs["rms_px"] <= 0:
        raise ValueError(f"the {R_COEFFICIENTS_KEY} detector parameters need a positive rms_px")
    return coeffs


def _rotator_terms(scienceAngleDeg, calibrationAngleDeg):
    """*the change in cos and sin of the rotator angle between the science and calibration frames*"""
    scienceAngle = math.radians(scienceAngleDeg)
    calibrationAngle = math.radians(calibrationAngleDeg)
    return (
        math.cos(scienceAngle) - math.cos(calibrationAngle),
        math.sin(scienceAngle) - math.sin(calibrationAngle),
    )


def predict_r_detector_shift(deltaIpx, scienceAngleDeg, calibrationAngleDeg, coeffs):
    """*predict the r detector shift from the i shift and the two rotator angles*

    **Key Arguments:**

    - ``deltaIpx`` -- the measured i-order detector shift (px, science minus calibration)
    - ``scienceAngleDeg`` -- the science-frame rotator angle (degrees)
    - ``calibrationAngleDeg`` -- the calibration-frame rotator angle (degrees)
    - ``coeffs`` -- the ``r-from-i-transfer`` coefficients, ``{"i": ..., "cos": ..., "sin": ..., "rms_px": ...}``

    **Return:**

    - ``shiftPx`` -- the predicted r detector shift in px

    **Usage:**

    ```python
    from soxspipe.commonutils.ug_wavelength_transfer import predict_r_detector_shift
    shiftPx = predict_r_detector_shift(0.5, -19.8, 0.0, coeffs)
    ```
    """
    deltaCos, deltaSin = _rotator_terms(scienceAngleDeg, calibrationAngleDeg)
    return coeffs["i"] * deltaIpx + coeffs["cos"] * deltaCos + coeffs["sin"] * deltaSin


def combine_r_shift_px(measuredPx, measuredSigmaPx, predictedPx, predictionSigmaPx):
    """*combine the measured and predicted r shifts by inverse variance*

    **Key Arguments:**

    - ``measuredPx`` -- the r shift measured from the sky lines (px)
    - ``measuredSigmaPx`` -- its uncertainty (px), floored at ``MIN_MEASURED_SIGMA_PX``
    - ``predictedPx`` -- the r shift predicted from i and the rotator angle (px)
    - ``predictionSigmaPx`` -- the model's rms (px)

    **Return:**

    - ``shiftPx`` -- the inverse-variance weighted mean (px)
    """
    weightMeasured = 1 / max(measuredSigmaPx, MIN_MEASURED_SIGMA_PX) ** 2
    weightPredicted = 1 / predictionSigmaPx**2
    return (weightMeasured * measuredPx + weightPredicted * predictedPx) / (weightMeasured + weightPredicted)


def usable_rotator_angle(header, kw):
    """*the frame's rotator angle in degrees, or None if it cannot be used in the u/g model*

    The angle is unusable when ``ABSROT`` is missing, boolean or not a finite number. It is also unusable for the
    placeholder that some real night arcs carry (2025-11-06): ``ABSROT`` and ``ESO TEL PARANG START`` both exactly 0.
    Those arcs are 3-4 px outliers in the u/g fit.

    **Key Arguments:**

    - ``header`` -- the FITS header to read
    - ``kw`` -- the keyword lookup, mapping an alias such as ``ABSROT`` to the header keyword

    **Return:**

    - ``angle`` -- the rotator angle in degrees, or *None*

    **Usage:**

    ```python
    from soxspipe.commonutils.ug_wavelength_transfer import usable_rotator_angle
    angle = usable_rotator_angle(frame.header, self.kw)
    ```
    """
    angle = read_rotator_angle(header, kw("ABSROT"))
    if angle is None:
        return None
    parallacticAngle = read_rotator_angle(header, kw("PAR_ANG_START"))
    if angle == 0.0 and parallacticAngle == 0.0:
        return None
    return angle


def attached_arc_offsets(scienceHeader, arcHeader, kw):
    """*how far an arc is from a science frame in time and rotator angle, or None if either cannot be measured*

    The time offset is the absolute difference of the two ``MJD-OBS`` start times. The rotator offset is the
    wrap-aware difference of the two usable rotator angles (``usable_rotator_angle``), so 179 and -179 degrees are 2
    degrees apart.

    **Key Arguments:**

    - ``scienceHeader`` -- the science frame header
    - ``arcHeader`` -- the arc header (the dispersion-map primary header)
    - ``kw`` -- the keyword lookup, mapping an alias such as ``ABSROT`` to the header keyword

    **Return:**

    - ``offsets`` -- ``(minutes, rotatorDegrees)``, or *None* if a frame has no finite ``MJD-OBS`` or usable angle

    **Usage:**

    ```python
    from soxspipe.commonutils.ug_wavelength_transfer import attached_arc_offsets
    offsets = attached_arc_offsets(frame.header, arcHeader, self.kw)
    ```
    """
    # ``read_rotator_angle`` READS ANY FINITE NUMBER FROM A HEADER, NOT ONLY AN ANGLE
    scienceMjd = read_rotator_angle(scienceHeader, kw("MJDOBS"))
    arcMjd = read_rotator_angle(arcHeader, kw("MJDOBS"))
    scienceAngle = usable_rotator_angle(scienceHeader, kw)
    arcAngle = usable_rotator_angle(arcHeader, kw)
    if None in (scienceMjd, arcMjd, scienceAngle, arcAngle):
        return None
    minutes = round(abs(scienceMjd - arcMjd) * MINUTES_PER_DAY, MINUTES_DECIMALS)
    return minutes, rotator_angle_spread([scienceAngle, arcAngle])


def is_attached_arc(scienceHeader, arcHeader, kw):
    """*whether an arc is attached to a science frame: close in time and in rotator angle*

    An arc is attached when it is at most ``MAX_ATTACHED_ARC_MINUTES`` from the science frame and less than
    ``MAX_ATTACHED_ARC_ROTATOR_DEGREES`` from its rotator angle. A frame whose offsets cannot be measured
    (``attached_arc_offsets``) is never attached.

    **Key Arguments:**

    - ``scienceHeader`` -- the science frame header
    - ``arcHeader`` -- the arc header (the dispersion-map primary header)
    - ``kw`` -- the keyword lookup, mapping an alias such as ``ABSROT`` to the header keyword

    **Return:**

    - ``isAttached`` -- *True* if the arc is attached

    **Usage:**

    ```python
    from soxspipe.commonutils.ug_wavelength_transfer import is_attached_arc
    isAttached = is_attached_arc(frame.header, arcHeader, self.kw)
    ```
    """
    return within_attached_arc_limits(attached_arc_offsets(scienceHeader, arcHeader, kw))


def within_attached_arc_limits(offsets):
    """*whether offsets from ``attached_arc_offsets`` are within the attached-arc limits; None never is*"""
    if offsets is None:
        return False
    minutes, rotatorDegrees = offsets
    return minutes <= MAX_ATTACHED_ARC_MINUTES and rotatorDegrees < MAX_ATTACHED_ARC_ROTATOR_DEGREES


def calibration_rotator_angle(header, kw):
    """*the calibration frame's rotator angle in degrees, falling back to the park angle*

    Daytime calibrations carry no usable ``ABSROT`` and are taken with the rotator parked at
    ``CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES`` (DY-1531), so that angle is used whenever ``usable_rotator_angle``
    returns *None*.

    **Key Arguments:**

    - ``header`` -- the calibration frame header (the dispersion-map primary header)
    - ``kw`` -- the keyword lookup, mapping an alias such as ``ABSROT`` to the header keyword

    **Return:**

    - ``angle`` -- the rotator angle in degrees
    """
    angle = usable_rotator_angle(header, kw)
    return CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES if angle is None else angle


def _check_dispersion(dispersionNmPerPx):
    """*raise ``ValueError`` unless the signed dispersion is finite and non-zero*"""
    if not math.isfinite(dispersionNmPerPx) or dispersionNmPerPx == 0:
        raise ValueError(f"the signed dispersion must be finite and non-zero, not {dispersionNmPerPx}")


def nm_shift_to_detector_px(shiftNm, dispersionNmPerPx):
    """*the detector motion (px) implied by a sky-tuning wavelength shift (nm)*

    The tuning adds ``shiftNm`` to the wavelength solution. A feature that moved by ``dy`` px now reads
    ``dispersion * dy`` nm away from its catalogue wavelength, so the tuning shift is ``-dispersion * dy`` and
    ``dy = -shiftNm / dispersion``.

    **Key Arguments:**

    - ``shiftNm`` -- the wavelength shift added to the solution (nm)
    - ``dispersionNmPerPx`` -- the signed dispersion d(wavelength)/d(pixel) of the order (nm/px)

    **Return:**

    - ``shiftPx`` -- the detector motion along the dispersion axis (px), positive towards larger pixel coordinates
    """
    _check_dispersion(dispersionNmPerPx)
    return -shiftNm / dispersionNmPerPx


def detector_px_to_nm_shift(shiftPx, dispersionNmPerPx):
    """*the wavelength shift (nm) that corrects a detector motion (px); the inverse of ``nm_shift_to_detector_px``*

    **Key Arguments:**

    - ``shiftPx`` -- the detector motion along the dispersion axis (px), positive towards larger pixel coordinates
    - ``dispersionNmPerPx`` -- the signed dispersion d(wavelength)/d(pixel) of the order (nm/px)

    **Return:**

    - ``shiftNm`` -- the wavelength shift to add to the solution (nm)
    """
    _check_dispersion(dispersionNmPerPx)
    return -dispersionNmPerPx * shiftPx


def _order_guard_reason(order, stats):
    """*why an r/i order's sky shift cannot drive a model, or None if its guards pass*"""
    if stats is None:
        return f"order {order}: no sky shift"
    if stats["nLines"] < MIN_MATCHED_LINES[order]:
        return f"order {order}: {stats['nLines']} sky lines (min {MIN_MATCHED_LINES[order]})"
    # AN ORDER TUNED OUTSIDE THE SOXS VIS PATH (X-SHOOTER) CARRIES NO SKY-LINE SNR AND IS NOT GATED ON IT
    if stats.get("lineSnr") is None:
        return None
    return sky_line_snr_reason(order, stats["lineSnr"])


def sky_line_snr_reason(order, lineSnr):
    """*why an order's sky is too faint for its own sky shift to be trusted, or None if it is bright enough*

    **Key Arguments:**

    - ``order`` -- the order
    - ``lineSnr`` -- the order's combined sky-line signal-to-noise (``sky_line_snr.combine_line_snrs``)

    **Return:**

    - ``reason`` -- a short reason (it fits a FITS header comment), or *None*
    """
    if not math.isfinite(lineSnr) or lineSnr < MIN_SKY_LINE_SNR:
        return f"order {order}: sky-line SNR {lineSnr:.1f} (min {MIN_SKY_LINE_SNR})"
    return None


def r_second_line_snr_reason(secondLineSnr):
    """*why a VIS r shift measured on its own rests on a noise-level second line, or None if that line is bright enough*

    **Key Arguments:**

    - ``secondLineSnr`` -- the signal-to-noise of the r order's second-brightest catalogue line

    **Return:**

    - ``reason`` -- a short reason (it fits a FITS header comment), or *None*
    """
    if not math.isfinite(secondLineSnr) or secondLineSnr < MIN_R_SECOND_LINE_SNR:
        return f"order {VIS_R_ORDER}: 2nd sky line SNR {secondLineSnr:.1f} (min {MIN_R_SECOND_LINE_SNR})"
    return None


def ug_skip_reason(orderStats, scienceAngleDeg):
    """*why the u/g model must not be applied, or None if every guard passes*

    A guard failure never fails the recipe: u and g are then left unshifted and the reason is logged. An r order whose
    shift came from the i order (``"transferred": True``) has no sky-line guards of its own.

    **Key Arguments:**

    - ``orderStats`` -- per-order sky-tuning statistics, ``{order: {"nLines": int, "lineSnr": float, ...}}``;
      ``nLines`` is the number of sky lines kept after clipping and ``lineSnr`` the order's combined sky-line
      signal-to-noise
    - ``scienceAngleDeg`` -- the usable science rotator angle in degrees, or *None*

    **Return:**

    - ``reason`` -- a short reason (it fits a FITS header comment), or *None*
    """
    for order in MEASURED_ORDERS:
        stats = orderStats.get(order)
        if stats is not None and stats.get("transferred"):
            continue
        reason = _order_guard_reason(order, stats)
        if reason is not None:
            return reason
    if scienceAngleDeg is None:
        return "no usable science ABSROT"
    return None


def r_prediction_skip_reason(iStats, scienceAngleDeg):
    """*why the r shift must not be predicted from i, or None if every guard passes*

    **Key Arguments:**

    - ``iStats`` -- the i-order sky-tuning statistics, ``{"nLines": int, "lineSnr": float}``, or *None*
    - ``scienceAngleDeg`` -- the usable science rotator angle in degrees, or *None*

    **Return:**

    - ``reason`` -- a short reason (it fits a FITS header comment), or *None*
    """
    reason = _order_guard_reason(VIS_I_ORDER, iStats)
    if reason is not None:
        return reason
    if scienceAngleDeg is None:
        return "no usable science ABSROT"
    return None
