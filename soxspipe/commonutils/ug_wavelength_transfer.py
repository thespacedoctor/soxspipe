"""Predict the SOXS VIS u and g order wavelength shifts, which have no usable sky lines, from the r and i sky shifts.

The u (order 2) and g (order 3) orders move against a daytime arc by up to 5-6 px at night. That motion is predicted
from the mean sky-measured detector shift of the r (order 1) and i (order 4) orders plus the change in the derotator
angle (``ESO ADA ABSROT``) between the science frame and the calibration frame (DY-1608):

    dy_order = m * (dr + di) / 2 + cos * (cos(th_sci) - cos(th_cal)) + sin * (sin(th_sci) - sin(th_cal))

All shifts are in detector px along the dispersion axis, science minus calibration, positive when the features move
to larger pixel coordinates. The coefficients live in the ``ug-wavelength-transfer`` block of the VIS detector
parameters.
"""

import math

from soxspipe.commonutils.slit_illumination_qc import CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES, read_rotator_angle

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
# MAXIMUM STANDARD DEVIATION OF THE KEPT PER-LINE SHIFTS IN AN r/i ORDER
MAX_SHIFT_SCATTER_PX = 1.5
COEFFICIENTS_KEY = "ug-wavelength-transfer"
MODEL_TERMS = ("m", "cos", "sin")


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
    scienceAngle = math.radians(scienceAngleDeg)
    calibrationAngle = math.radians(calibrationAngleDeg)
    deltaCos = math.cos(scienceAngle) - math.cos(calibrationAngle)
    deltaSin = math.sin(scienceAngle) - math.sin(calibrationAngle)
    return {
        order: coeffs[order]["m"] * meanShift + coeffs[order]["cos"] * deltaCos + coeffs[order]["sin"] * deltaSin
        for order in PREDICTED_ORDERS
    }


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


def ug_skip_reason(orderStats, scienceAngleDeg):
    """*why the u/g model must not be applied, or None if every guard passes*

    A guard failure never fails the recipe: u and g are then left unshifted and the reason is logged.

    **Key Arguments:**

    - ``orderStats`` -- per-order sky-tuning statistics, ``{order: {"nLines": int, "scatterPx": float, ...}}``;
      ``nLines`` is the number of sky lines kept after clipping and ``scatterPx`` the standard deviation of their
      shifts (px)
    - ``scienceAngleDeg`` -- the usable science rotator angle in degrees, or *None*

    **Return:**

    - ``reason`` -- a short reason (it fits a FITS header comment), or *None*
    """
    for order in MEASURED_ORDERS:
        stats = orderStats.get(order)
        if stats is None:
            return f"order {order}: no sky shift"
        if stats["nLines"] < MIN_MATCHED_LINES[order]:
            return f"order {order}: {stats['nLines']} sky lines (min {MIN_MATCHED_LINES[order]})"
        scatter = stats["scatterPx"]
        if not math.isfinite(scatter) or scatter > MAX_SHIFT_SCATTER_PX:
            return f"order {order}: shift scatter {scatter:.2f} px (max {MAX_SHIFT_SCATTER_PX})"
    if scienceAngleDeg is None:
        return "no usable science ABSROT"
    return None
