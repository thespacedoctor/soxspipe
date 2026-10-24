"""How the global trace fit treats the synthetic rows of a transferred faint order (DY-1941, round 2)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum

pytestmark = pytest.mark.unit

ORDERS = (1.0, 2.0, 3.0)
SYNTHETIC_ORDER = 4.0
SMALL_SYNTHETIC_COUNT = 12
ROWS_PER_ORDER = 60
OUTLIER_ROWS = {5: 4.0, 17: -5.0, 33: 6.5, 48: -4.5}

# THE FIT OF `_measured_rows` AS IT WAS BEFORE THE faint_transfer COLUMN WAS HANDLED, COMPUTED FROM THE UNCHANGED FIT.
# THE COEFFICIENTS ARE COMPARED TO 1E-12 BECAUSE LAPACK ROUNDS THE LAST BIT DIFFERENTLY ON DIFFERENT CPUS
# fmt: off
BASELINE_X_COEFFICIENTS = [
    0.3629690089524446,
    0.009565695693405025,
    1.0647662906717527e-06,
    49.81470605474839,
    0.00023878790442548198,
    -3.8843147491783944e-08,
]
BASELINE_X_CLIPPED_ROWS = [
    5, 18, 21, 23, 55, 62, 65, 72, 77, 79, 84,
    88, 93, 99, 104, 108, 115, 121, 154, 165, 167, 178,
]
BASELINE_X_QCS = {"X RES MIN": -0.461, "X RES MAX": 0.47, "X RES SD": 0.239, "X RES MEDIAN": 0.203, "N ORDERS": 3.0}
BASELINE_STDDEV_COEFFICIENTS = [
    1.502138579123934,
    0.00010079249745908142,
    -1.4972212570940447e-10,
    -0.0023665286327763397,
    3.7836640439437196e-06,
    -1.397455145025835e-09,
]
BASELINE_STDDEV_CLIPPED_ROWS = [
    2, 9, 13, 23, 26, 42, 45, 47, 52, 70,
    84, 100, 120, 146, 153, 154, 170, 171, 174, 175, 178,
]
# fmt: on


def _detector(log: object) -> detect_continuum:
    detector = object.__new__(detect_continuum)
    detector.log, detector.arm, detector.axisA, detector.axisB = log, "VIS", "x", "y"
    detector.orderDeg, detector.axisBDeg = 1, 2
    detector.recipeSettings = {"poly-fitting-residual-clipping-sigma": 2.0, "poly-clipping-iteration-limit": 5}
    detector.recipeName = "soxs-stare"
    detector.dateObs = "2024-01-02T03:04:05"
    detector.qc = pd.DataFrame()
    return detector


def _measured_rows() -> pd.DataFrame:
    """Return noisy samples of three orders on a smooth trace, with a few outliers in the middle order."""
    rng = np.random.default_rng(1941)
    frames = []
    for order in ORDERS:
        y = np.linspace(0.0, 2000.0, ROWS_PER_ORDER)
        x = 50.0 * order + 0.01 * y + 1.0e-6 * y**2 + rng.normal(0.0, 0.3, ROWS_PER_ORDER)
        if order == 2.0:
            for row, jump in OUTLIER_ROWS.items():
                x[row] += jump
        frames.append(
            pd.DataFrame(
                {
                    "order": order,
                    "cont_x": x,
                    "cont_y": y,
                    "gauss_stddev": 1.5 + 0.0001 * y + rng.normal(0.0, 0.02, ROWS_PER_ORDER),
                    "pre-clipped": False,
                }
            )
        )
    table = pd.concat(frames, ignore_index=True)
    table["row_id"] = table.index
    return table


def _fit(detector: detect_continuum, table: pd.DataFrame) -> tuple:
    coefficients, fitted, clipped = detector.fit_global_polynomial(table, writeQCs=True)
    return coefficients, fitted, clipped, detector.qc.set_index("qc_name")["qc_value"].to_dict()


def test_a_table_without_the_faint_transfer_column_gives_the_pre_round_two_trace_fit(log: object) -> None:
    # ARRANGE
    table = _measured_rows()

    # ACT
    coefficients, fitted, clipped, qcs = _fit(_detector(log), table)

    # ASSERT: THE SAME COEFFICIENTS, THE SAME CLIPPED ROWS AND THE SAME QCS AS BEFORE THE CHANGE
    np.testing.assert_allclose(coefficients, BASELINE_X_COEFFICIENTS, rtol=1e-12, atol=0)
    assert sorted(clipped["row_id"]) == BASELINE_X_CLIPPED_ROWS
    assert len(fitted) == len(table) - len(BASELINE_X_CLIPPED_ROWS)
    assert qcs == BASELINE_X_QCS


def test_a_table_without_the_faint_transfer_column_gives_the_pre_round_two_stddev_fit(log: object) -> None:
    # ARRANGE
    table = _measured_rows()

    # ACT
    coefficients, fitted, clipped = _detector(log).fit_global_polynomial(table, axisACol="gauss_stddev")

    # ASSERT
    np.testing.assert_allclose(coefficients, BASELINE_STDDEV_COEFFICIENTS, rtol=1e-12, atol=0)
    assert sorted(clipped["row_id"]) == BASELINE_STDDEV_CLIPPED_ROWS
    assert len(fitted) == len(table) - len(BASELINE_STDDEV_CLIPPED_ROWS)


def _with_transferred_order(table: pd.DataFrame, offsetPx: float = 0.0, count: int = ROWS_PER_ORDER) -> pd.DataFrame:
    """Return `table` plus `count` rows of one more order on the noise-free trace, marked `faint_transfer`.

    `offsetPx` moves the synthetic rows off that trace.
    """
    y = np.linspace(0.0, 2000.0, count)
    synthetic = pd.DataFrame(
        {
            "order": SYNTHETIC_ORDER,
            "cont_x": 50.0 * SYNTHETIC_ORDER + 0.01 * y + 1.0e-6 * y**2 + offsetPx,
            "cont_y": y,
            "gauss_stddev": 1.5 + 0.0001 * y,
            "pre-clipped": False,
            "row_id": range(len(table), len(table) + count),
        }
    )
    return pd.concat([table.assign(faint_transfer=False), synthetic.assign(faint_transfer=True)], ignore_index=True)


def _measured_only(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.loc[~frame["faint_transfer"].astype(bool)]


def test_synthetic_rows_are_never_clipped_however_far_they_sit_from_the_fit(log: object) -> None:
    # ARRANGE: TWELVE SYNTHETIC ROWS 8 PX OFF THE TRACE, FAR BEYOND THE 2-SIGMA CLIP OF 180 MEASURED ROWS
    table = _with_transferred_order(_measured_rows(), offsetPx=8.0, count=SMALL_SYNTHETIC_COUNT)

    # ACT
    _, fitted, clipped, _ = _fit(_detector(log), table)

    # ASSERT
    assert not clipped["faint_transfer"].astype(bool).any()
    assert fitted["faint_transfer"].astype(bool).sum() == SMALL_SYNTHETIC_COUNT


def test_a_bright_order_keeps_its_samples_when_synthetic_rows_are_added(log: object) -> None:
    # ARRANGE: EXACT SYNTHETIC ROWS ON THE TRACE WOULD SHRINK THE CLIPPING SPREAD IF THEY COUNTED IN IT
    measured = _measured_rows()
    withSynthetic = _with_transferred_order(measured)

    # ACT
    _, fittedAlone, _, _ = _fit(_detector(log), measured.assign(faint_transfer=False))
    _, fittedWith, _, _ = _fit(_detector(log), withSynthetic)

    # ASSERT: THE KEPT COUNT OF EACH MEASURED ORDER DOES NOT DROP
    keptAlone = fittedAlone.groupby("order").size()
    keptWith = _measured_only(fittedWith).groupby("order").size()
    assert (keptWith >= keptAlone).all()


def test_synthetic_rows_are_left_out_of_the_x_res_qcs(log: object) -> None:
    # ARRANGE: SYNTHETIC ROWS 8 PX OFF THE TRACE WOULD DOMINATE THE RESIDUAL STATISTICS
    table = _with_transferred_order(_measured_rows(), offsetPx=8.0)

    # ACT
    _, fitted, _, qcs = _fit(_detector(log), table)

    # ASSERT: THE QCS ARE THE STATISTICS OF THE MEASURED RESIDUALS, AND THE SYNTHETIC ORDER STILL COUNTS AS AN ORDER
    residuals = _measured_only(fitted)["cont_x_fit_res"].to_numpy()
    assert qcs["X RES MIN"] == round(float(residuals.min()), 3)
    assert qcs["X RES MAX"] == round(float(residuals.max()), 3)
    assert qcs["X RES SD"] == round(float(residuals.std()), 3)
    assert qcs["X RES MEDIAN"] == round(float(np.abs(residuals).mean()), 3)
    assert qcs["N ORDERS"] == len(ORDERS) + 1


def test_synthetic_rows_still_count_toward_the_polynomial_fit(log: object) -> None:
    # ARRANGE: WITH A DEGREE-0 FIT IN BOTH AXES THE FIT IS THE MEAN, SO ROWS FAR FROM THE REST MOVE IT
    detector = _detector(log)
    detector.orderDeg, detector.axisBDeg = 0, 0
    measured = _measured_rows().assign(faint_transfer=False)
    shifted = _with_transferred_order(_measured_rows(), offsetPx=8.0)

    # ACT
    alone, _, _ = _detector(log).fit_global_polynomial(measured)
    withSynthetic, _, _ = detector.fit_global_polynomial(shifted)

    # ASSERT
    assert withSynthetic[0] != pytest.approx(alone[0], abs=1e-6)


def test_a_table_of_only_synthetic_rows_still_fits_and_reports_its_residuals(log: object) -> None:
    # ARRANGE: NO MEASURED ROW IS LEFT TO SET THE CLIPPING SPREAD OR THE QC STATISTICS
    table = _with_transferred_order(_measured_rows().iloc[0:0], offsetPx=0.0)
    detector = _detector(log)
    detector.orderDeg, detector.axisBDeg = 0, 2

    # ACT
    coefficients, fitted, clipped, qcs = _fit(detector, table)

    # ASSERT: ALL THE ROWS ARE KEPT, AND THE QCS FALL BACK TO THEIR RESIDUALS
    assert coefficients is not None
    assert len(fitted) == ROWS_PER_ORDER
    assert len(clipped) == 0
    assert abs(qcs["X RES MAX"]) < 1e-3
