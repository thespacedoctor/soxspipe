"""Slit shifts of a frame group, relative to its first frame in time (DY-1612)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.slit_alignment import slit_shifts

pytestmark = pytest.mark.unit


def _positions(rows: dict[str, tuple[float, list[float]]]) -> pd.DataFrame:
    """Return a position table from ``{frame: (mjd, [position of order 1, 2, ...])}``."""
    records = [
        {"frame": frame, "mjd": mjd, "order": order, "position_arcsec": position}
        for frame, (mjd, orderPositions) in rows.items()
        for order, position in enumerate(orderPositions, start=1)
    ]
    return pd.DataFrame(records)


def test_each_order_shift_is_its_position_minus_the_first_frame_in_time(log: object) -> None:
    # ARRANGE: THE FIRST FRAME IN TIME IS LISTED LAST
    positions = _positions(
        {
            "b.fits": (60000.2, [1.5, -0.5, 2.0, 0.0]),
            "c.fits": (60000.3, [0.0, 0.0, 0.0, 0.0]),
            "a.fits": (60000.1, [1.0, 1.0, 1.0, 1.0]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=True)

    # ASSERT
    byFrame = {frame: group.set_index("order")["shift_arcsec"] for frame, group in shifts.groupby("frame")}
    np.testing.assert_allclose(byFrame["a.fits"].to_numpy(), [0.0, 0.0, 0.0, 0.0])
    np.testing.assert_allclose(byFrame["b.fits"].to_numpy(), [0.5, -1.5, 1.0, -1.0])
    np.testing.assert_allclose(byFrame["c.fits"].to_numpy(), [-1.0, -1.0, -1.0, -1.0])
    assert shifts.loc[shifts["is_reference"], "frame"].unique().tolist() == ["a.fits"]
    assert fallbackCount == 0
    assert not shifts["is_fallback"].any()


def test_an_order_with_no_trace_takes_the_median_shift_of_the_other_orders_and_warns(log: object) -> None:
    # ARRANGE: ORDER 3 OF b.fits HAS NO RELIABLE TRACE. ITS OTHER ORDERS SHIFT BY 0.2, 0.4 AND 0.9
    positions = _positions(
        {
            "a.fits": (60000.1, [0.0, 0.0, 0.0, 0.0]),
            "b.fits": (60000.2, [0.2, 0.4, np.nan, 0.9]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=True)

    # ASSERT
    fallback = shifts.loc[shifts["is_fallback"]]
    assert fallback[["frame", "order"]].to_dict("records") == [{"frame": "b.fits", "order": 3}]
    assert fallback["shift_arcsec"].iloc[0] == pytest.approx(0.4)
    assert fallbackCount == 1
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "b.fits" in warnings[0]
    assert "order 3" in warnings[0]


def test_a_frame_with_no_trace_in_any_order_takes_the_median_of_the_other_frames_in_each_order(log: object) -> None:
    # ARRANGE: c.fits HAS NO TRACE AT ALL. THE OTHER FRAMES SHIFT ORDER 1 BY 0.0, 1.0 AND 3.0 AND ORDER 2 BY THE
    # NEGATIVE OF THAT. THE REFERENCE FRAME COUNTS AS A FRAME WITH A SHIFT OF 0
    positions = _positions(
        {
            "a.fits": (60000.1, [0.0, 0.0]),
            "b.fits": (60000.2, [1.0, -1.0]),
            "c.fits": (60000.3, [np.nan, np.nan]),
            "d.fits": (60000.4, [3.0, -3.0]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=True)

    # ASSERT
    unreliable = shifts.loc[shifts["frame"] == "c.fits"].set_index("order")
    assert unreliable["shift_arcsec"].to_dict() == {1: pytest.approx(1.0), 2: pytest.approx(-1.0)}
    assert unreliable["is_fallback"].all()
    assert fallbackCount == 2


def test_a_group_with_no_measured_shift_anywhere_is_left_unshifted(log: object) -> None:
    # ARRANGE
    positions = _positions(
        {
            "a.fits": (60000.1, [0.0, 0.0]),
            "b.fits": (60000.2, [np.nan, np.nan]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=True)

    # ASSERT
    assert shifts["shift_arcsec"].tolist() == [0.0, 0.0, 0.0, 0.0]
    assert fallbackCount == 2


def test_a_whole_frame_shift_is_the_median_over_orders_of_the_position_differences(log: object) -> None:
    # ARRANGE: THE ORDERS OF b.fits MOVE BY 1.0, 1.2 AND 5.0 (AN OUTLIER), SO THE FRAME MOVES BY 1.2
    positions = _positions(
        {
            "a.fits": (60000.1, [0.0, 0.0, 0.0]),
            "b.fits": (60000.2, [1.0, 1.2, 5.0]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=False)

    # ASSERT
    assert shifts["order"].tolist() == [-1, -1]
    assert shifts["frame"].tolist() == ["a.fits", "b.fits"]
    assert shifts["shift_arcsec"].tolist() == pytest.approx([0.0, 1.2])
    assert fallbackCount == 0


def test_a_whole_frame_with_no_trace_takes_the_median_shift_of_the_other_frames_and_warns(log: object) -> None:
    # ARRANGE: THE OTHER FRAMES SHIFT BY 0.0 (THE REFERENCE), 1.0 AND 3.0
    positions = _positions(
        {
            "a.fits": (60000.1, [0.0, 0.0]),
            "b.fits": (60000.2, [1.0, 1.0]),
            "c.fits": (60000.3, [np.nan, np.nan]),
            "d.fits": (60000.4, [3.0, 3.0]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=False)

    # ASSERT
    fallback = shifts.loc[shifts["is_fallback"]]
    assert fallback["frame"].tolist() == ["c.fits"]
    assert fallback["shift_arcsec"].tolist() == pytest.approx([1.0])
    assert fallbackCount == 1
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "c.fits" in warnings[0]


def test_a_reference_order_with_no_trace_is_warned_about_and_counted(log: object) -> None:
    # ARRANGE: THE FIRST FRAME IN TIME HAS NO TRACE IN ORDER 2, SO NO OTHER FRAME CAN BE MEASURED AGAINST IT THERE
    positions = _positions(
        {
            "a.fits": (60000.1, [0.0, np.nan]),
            "b.fits": (60000.2, [1.0, 5.0]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=True)

    # ASSERT: THE REFERENCE ROW AND THE B ROW OF ORDER 2 ARE BOTH FALLBACKS
    fallback = shifts.loc[shifts["is_fallback"]]
    assert fallback[["frame", "order"]].to_dict("records") == [
        {"frame": "a.fits", "order": 2},
        {"frame": "b.fits", "order": 2},
    ]
    assert fallbackCount == 2
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("reference" in message and "a.fits" in message and "order 2" in message for message in warnings)


def test_a_reference_frame_with_no_trace_is_warned_about_and_counted_in_whole_frame_mode(log: object) -> None:
    # ARRANGE
    positions = _positions(
        {
            "a.fits": (60000.1, [np.nan, np.nan]),
            "b.fits": (60000.2, [1.0, 1.0]),
        }
    )

    # ACT
    shifts, fallbackCount = slit_shifts(log, positions, perOrder=False)

    # ASSERT
    assert shifts.loc[shifts["is_reference"], "is_fallback"].tolist() == [True]
    assert fallbackCount == 2
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("reference" in message and "a.fits" in message for message in warnings)
