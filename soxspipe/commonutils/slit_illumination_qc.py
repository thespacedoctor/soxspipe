"""QC for the slit-illumination measurement: the per-order QC plot and the rotator-angle spread.

The plot shows, for every order, the slit profile the offset was fitted to (from skylines, or from the sky continuum
for an order with no usable skylines) against the master flat's lit slit range, so a poor fit or a failed order can be
seen at a glance.
"""

import math

import numpy as np

from soxspipe.commonutils.slit_illumination import METHOD_SKYLINE, top_hat
from soxspipe.commonutils.slit_illumination_mask import kept_slit_range, usable_flat_edges, usable_offsets

__all__ = ["plot_slit_illumination_qc", "read_rotator_angle", "rotator_angle_spread"]

PANEL_COLUMNS = 3
PANEL_SIZE_INCHES = (5.5, 3.2)
# THE DRAWN TOP-HATS ARE SCALED TO THIS PERCENTILE OF THE MEASURED PROFILE, AS ITS AMPLITUDE IS ARBITRARY
PROFILE_LEVEL_PERCENTILE = 90
CURVE_SAMPLES = 400
# THE DRAWN SLIT RANGE EXTENDS THIS FAR BEYOND THE FLAT EDGES WHEN THERE IS NO PROFILE
CURVE_PADDING_ARCSEC = 2.0
KEPT_RANGE_ALPHA = 0.15
FULL_TURN_DEGREES = 360.0
# DAYTIME CALIBRATION FRAMES CARRY NO ABSROT KEYWORD AND ARE TAKEN WITH THE ROTATOR PARKED AT THIS ANGLE (DY-1531)
CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES = 0.0


def read_rotator_angle(header, key):
    """*the rotator angle in degrees from a header, or None if the keyword is missing or not a finite number*

    **Key Arguments:**

    - ``header`` -- the FITS header to read
    - ``key`` -- the rotator-angle keyword

    **Return:**

    - ``angle`` -- the angle in degrees, or *None* if the keyword is missing, boolean or not a finite number

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination_qc import read_rotator_angle
    angle = read_rotator_angle(header, "ESO ADA ABSROT END")
    ```
    """
    try:
        value = header[key]
        # A BOOLEAN IS NOT AN ANGLE, EVEN THOUGH float() ACCEPTS IT
        if isinstance(value, bool | np.bool_):
            return None
        angle = float(value)
    except (KeyError, TypeError, ValueError, OverflowError):
        return None
    return angle if math.isfinite(angle) else None


def rotator_angle_spread(angles):
    """*the spread of rotator angles: the smallest arc, in degrees, that holds every angle*

    The arc is found across the wrap at +/-180 degrees, so 179 and -179 degrees are 2 degrees apart. Angles that are
    not finite are ignored.

    **Key Arguments:**

    - ``angles`` -- the rotator angles in degrees

    **Return:**

    - ``spread`` -- the angular spread in degrees, 0 if fewer than two finite angles are given

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination_qc import rotator_angle_spread
    spread = rotator_angle_spread([179.0, -179.0])
    ```
    """
    finite = np.array([angle for angle in angles if math.isfinite(angle)], dtype=float)
    if finite.size < 2:
        return 0.0
    ordered = np.sort(np.mod(finite, FULL_TURN_DEGREES))
    gaps = np.diff(np.append(ordered, ordered[0] + FULL_TURN_DEGREES))
    return float(FULL_TURN_DEGREES - gaps.max())


def _panel_orders(result):
    """*the sorted order numbers that have profile rows or an offset row*"""
    orders = set(result.profiles["order"].astype(int)) | set(result.offsets["order"].astype(int))
    return sorted(orders)


def _slit_grid(profile, flatEdge):
    """*the slit positions at which the top-hats are drawn*"""
    if len(profile):
        low, high = profile["slit"].min(), profile["slit"].max()
    elif flatEdge is not None:
        low, high = flatEdge[0] - CURVE_PADDING_ARCSEC, flatEdge[1] + CURVE_PADDING_ARCSEC
    else:
        low, high = -1.0, 1.0
    return np.linspace(low, high, CURVE_SAMPLES)


def _order_method(offsetRow):
    """*the method that measured an order's profile (``skyline`` or ``continuum``), ``skyline`` if none is recorded*"""
    method = None if offsetRow is None else offsetRow.get("method")
    return method if isinstance(method, str) else METHOD_SKYLINE


def _draw_profile(axis, profile, method):
    """*draw the measured slit profile with its errors; return the level the top-hats are scaled to*"""
    if not len(profile):
        return 1.0
    axis.errorbar(
        profile["slit"],
        profile["profile"],
        yerr=profile["profileErr"],
        fmt=".",
        ms=3,
        color="black",
        lw=0.5,
        label=f"{method} profile",
    )
    values = profile["profile"].to_numpy(dtype=float)
    level = np.nanpercentile(values, PROFILE_LEVEL_PERCENTILE) if np.isfinite(values).any() else 1.0
    return level if np.isfinite(level) and level > 0 else 1.0


def _draw_hats(axis, grid, level, flatEdge, flatWidth, offset, marginArcsec):
    """*draw the flat's top-hat, the shifted top-hat and the kept range, for whichever of them can be drawn*"""
    if flatEdge is None:
        return
    lower, upper = flatEdge
    axis.plot(grid, top_hat(grid, level, lower, upper, flatWidth), color="tab:blue", lw=1, label="master flat")
    if offset is None:
        return
    axis.plot(
        grid,
        top_hat(grid, level, lower + offset, upper + offset, flatWidth),
        color="tab:red",
        lw=1,
        label=f"flat shifted by {offset:+.3f} arcsec",
    )
    keepLower, keepUpper = kept_slit_range(lower, upper, offset, marginArcsec)
    if keepLower < keepUpper:
        axis.axvspan(keepLower, keepUpper, color="tab:green", alpha=KEPT_RANGE_ALPHA, label="kept range")


def _panel_title(order, offsetRow):
    """*the panel title: the order and its status, with the reason when the offset was not fitted*"""
    if offsetRow is None:
        return f"order {order}: no offset measured"
    title = f"order {order}: {offsetRow['status']}"
    if isinstance(offsetRow.get("method"), str):
        title = f"{title} ({offsetRow['method']})"
    return f"{title} ({offsetRow['reason']})" if offsetRow["reason"] else title


def _draw_order(axis, order, result, usableOffsets, flatByOrder, flatWidths):
    """*draw one order's panel*"""
    profile = result.profiles.loc[result.profiles["order"] == order]
    offsetRows = result.offsets.loc[result.offsets["order"] == order]
    offsetRow = offsetRows.iloc[0] if len(offsetRows) else None
    flatEdge = flatByOrder.get(order)
    grid = _slit_grid(profile, flatEdge)
    level = _draw_profile(axis, profile, _order_method(offsetRow))
    _draw_hats(
        axis, grid, level, flatEdge, flatWidths.get(order, 0.15), usableOffsets.get(order), result.overlapMarginArcsec
    )
    isFailed = offsetRow is None or order not in usableOffsets
    axis.set_title(_panel_title(order, offsetRow), fontsize=7, color="red" if isFailed else "black")
    axis.set_xlabel("slit position (arcsec)", fontsize=7)
    axis.tick_params(labelsize=7)
    axis.legend(fontsize=6, loc="lower center")


def plot_slit_illumination_qc(result, title, filePath):
    """*plot, for every order, the slit profile against the master flat's lit range and the fitted offset*

    Each panel shows the measured skyline or continuum profile with its errors, the flat's top-hat, the top-hat shifted
    by the order's offset (when the offset is usable), and the slit range kept by the overlap mask. The panel title
    gives the order's status, the method that measured it and, for a failed or borrowed order, the reason.

    **Key Arguments:**

    - ``result`` -- the ``slit_illumination_result`` from ``measure_slit_illumination``
    - ``title`` -- the figure title
    - ``filePath`` -- the PDF path to write

    **Return:**

    - ``filePath`` -- the path the plot was written to

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_illumination_qc import plot_slit_illumination_qc
    plot_slit_illumination_qc(result, "title", "/path/to/plot.pdf")
    ```
    """
    import matplotlib.pyplot as plt

    orders = _panel_orders(result)
    nColumns = min(max(len(orders), 1), PANEL_COLUMNS)
    nRows = max(math.ceil(len(orders) / nColumns), 1)
    fig, axes = plt.subplots(
        nRows,
        nColumns,
        figsize=(PANEL_SIZE_INCHES[0] * nColumns, PANEL_SIZE_INCHES[1] * nRows),
        squeeze=False,
    )
    try:
        usableOffsets = usable_offsets(result.offsets)
        flatByOrder = usable_flat_edges(result.flatEdges)
        flatWidths = dict(
            zip(result.flatEdges["order"].astype(int), result.flatEdges["width"].astype(float), strict=True)
        )
        panels = axes.ravel()
        for axis, order in zip(panels, orders, strict=False):
            _draw_order(axis, order, result, usableOffsets, flatByOrder, flatWidths)
        for axis in panels[len(orders) :]:
            axis.set_visible(False)
        fig.suptitle(title, fontsize=10)
        fig.tight_layout()
        fig.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    return filePath
