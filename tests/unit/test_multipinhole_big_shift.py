"""DY-1618: the first-iteration multi-pinhole shift picks the slit-centred pinhole identification (SOXS VIS).

When the pinhole at one end of the slit is missing from the multi-pinhole frame, its detection window latches onto
the neighbouring pinhole. The end-to-end shift test then sees a one-pinhole jump at that end only. That is also what
a real one-pinhole offset between the first-guess map and the frame looks like, so the shift is chosen between no
shift and one pinhole spacing either way by centring the 9-pinhole pattern between the order's slit edges.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.create_dispersion_map import PINHOLE_BEYOND_END_MIN_SETS, create_dispersion_map
from tests.factories import order_table_fits

pytestmark = pytest.mark.unit

SLIT_LOW_EDGE = 100.0
SLIT_HIGH_EDGE = 140.0
SLIT_CENTRE = 0.5 * (SLIT_LOW_EDGE + SLIT_HIGH_EDGE)
PINHOLE_SPACING = 4.0
LINE_ROWS = np.arange(10.0, 70.0, 2.0)
FRAME_ROWS = np.arange(6.0, 100.0, 6.0)
PINHOLE_FLUX = 800.0
REMAINDER = 0.8
FRAME_SHAPE = (110, 170)
SPOT_SIGMA = 1.2
# SLIT EDGES TOO NARROW FOR THE PROBE POSITION BEYOND THE PATTERN END TO LIE WITHIN THEIR MARGIN
NARROW_SLIT_LOW_EDGE = 114.0
NARROW_SLIT_HIGH_EDGE = 126.0
# WHERE THE OPPOSITE-END SEARCH LOOKS IN THE LATCHED-TOP TABLE: ONE SPACING BEYOND PINHOLE 0, SHIFTED BY THE REMAINDER
BEYOND_END_PROBE_X = SLIT_LOW_EDGE + REMAINDER


def _order_table(
    tmp_path: Path,
    orders: tuple[int, ...] = (2, 3),
    lowEdge: float = SLIT_LOW_EDGE,
    highEdge: float = SLIT_HIGH_EDGE,
) -> str:
    """Write an order table whose slit runs from ``lowEdge`` to ``highEdge`` at every row of every order."""
    polynomials = pd.DataFrame(
        [
            {
                "degorder_cent": 0,
                "degy_cent": 0,
                "cent_00": 0.5 * (lowEdge + highEdge),
                "degorder_edgeup": 0,
                "degy_edgeup": 0,
                "edgeup_00": highEdge,
                "degorder_edgelow": 0,
                "degy_edgelow": 0,
                "edgelow_00": lowEdge,
            }
        ]
    )
    metadata = pd.DataFrame({"order": list(orders), "ymin": 0.0, "ymax": 100.0})
    return str(order_table_fits(tmp_path / "orders.fits", polynomials=polynomials, metadata=metadata))


def _mapper(log: object, orderTable: str, pinholeFrame: np.ma.MaskedArray | None = None) -> create_dispersion_map:
    """Return a SOXS VIS multi-pinhole mapper that reads only what the shift step needs.

    With a ``pinholeFrame`` the mapper can also look for pinholes in the image, as the detection loop leaves it.
    """
    mapper = create_dispersion_map.__new__(create_dispersion_map)
    mapper.pinholeFrameMasked = pinholeFrame
    mapper.orderGeometryCache = None
    mapper.windowHalf = 15
    mapper.turnOffMP = True
    mapper.log = log
    mapper.debug = False
    mapper.arm = "VIS"
    mapper.inst = "SOXS"
    mapper.axisA = "x"
    mapper.axisB = "y"
    mapper.orderTable = orderTable
    mapper.kw = lambda name: name
    mapper.pinholeFrame = type("frame", (), {"header": {}})()
    return mapper


def _line_table(
    predictedCentre: float,
    observedSlitIndex: dict[int, int],
    order: int = 2,
    lineRows: np.ndarray = LINE_ROWS,
    observedOffset: float = 0.0,
) -> pd.DataFrame:
    """Return 9-pinhole line sets for ``order`` whose predictions are centred on ``predictedCentre``.

    ``observedSlitIndex`` maps each predicted slit index to the index of the real pinhole its window detected, and
    ``observedOffset`` moves every detection along the axis.
    """
    rows = []
    for wavelengthIndex, y in enumerate(lineRows):
        for slitIndex in range(9):
            predictedX = predictedCentre + (slitIndex - 4) * PINHOLE_SPACING
            observedX = SLIT_CENTRE + (observedSlitIndex[slitIndex] - 4) * PINHOLE_SPACING + observedOffset
            rows.append(
                {
                    "wavelength": 400.0 + wavelengthIndex,
                    "order": order,
                    "slit_index": slitIndex,
                    "detector_x_shifted": predictedX,
                    "detector_y_shifted": y,
                    "observed_x": observedX,
                    "observed_y": y,
                    "x_diff": predictedX - observedX,
                    "y_diff": 0.0,
                }
            )
    return pd.DataFrame(rows)


def _run(mapper: create_dispersion_map, table: pd.DataFrame) -> pd.DataFrame:
    """Run the first-iteration shift step on each order's shift group in turn, as the detection loop does."""
    result = table.copy()
    for order in result["order"].unique():
        result = mapper._handle_multipin_hole_big_shift(result, order, result["order"] == order, 1)
    return result


ALL_PINHOLES = {slitIndex: slitIndex for slitIndex in range(9)}


def test_missing_top_pinhole_does_not_shift_the_order(log: object, tmp_path: Path) -> None:
    """The top pinhole of order 2 is absent, so its window latches onto pinhole 7; neither order is shifted."""
    # ARRANGE
    mapper = _mapper(log, _order_table(tmp_path))
    observed = {slitIndex: slitIndex for slitIndex in range(8)} | {8: 7}
    table = pd.concat(
        [_line_table(SLIT_CENTRE, observed, order=2), _line_table(SLIT_CENTRE, ALL_PINHOLES, order=3)],
        ignore_index=True,
    )

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"])


def test_real_one_pinhole_offset_is_still_corrected(log: object, tmp_path: Path) -> None:
    """The first guess sits one pinhole low; the whole order moves up one spacing onto the slit centre."""
    # ARRANGE
    mapper = _mapper(log, _order_table(tmp_path))
    # PREDICTED INDEX k SITS ON REAL PINHOLE k-1; INDEX 0 HAS NO PINHOLE BELOW IT AND LATCHES ONTO PINHOLE 0
    observed = {0: 0} | {slitIndex: slitIndex - 1 for slitIndex in range(1, 9)}
    table = _line_table(SLIT_CENTRE - PINHOLE_SPACING, observed)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] + PINHOLE_SPACING)


def test_sub_spacing_offset_is_measured_when_the_end_shifts_have_opposite_signs(log: object, tmp_path: Path) -> None:
    """The prediction sits 1.5 px off and the top pinhole is missing: the shift is the 1.5 px, not 0 and not 2.5 px."""
    # ARRANGE
    mapper = _mapper(log, _order_table(tmp_path))
    offset = 1.5
    observed = {slitIndex: slitIndex for slitIndex in range(8)} | {8: 7}
    table = _line_table(SLIT_CENTRE - offset, observed)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] + offset)


def test_one_pinhole_offset_plus_a_sub_spacing_remainder_is_corrected_in_full(log: object, tmp_path: Path) -> None:
    """The first guess sits one spacing and 1 px low: the whole spacing and the 1 px both move the order."""
    # ARRANGE
    mapper = _mapper(log, _order_table(tmp_path))
    remainder = 1.0
    observed = {0: 0} | {slitIndex: slitIndex - 1 for slitIndex in range(1, 9)}
    table = _line_table(SLIT_CENTRE - PINHOLE_SPACING - remainder, observed)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] + PINHOLE_SPACING + remainder)


def test_shift_log_line_reports_the_whole_spacing_part_and_the_remainder(log: object, tmp_path: Path) -> None:
    """The log line gives the whole-spacing part (in spacings and px) and the sub-spacing remainder."""
    # ARRANGE
    mapper = _mapper(log, _order_table(tmp_path))
    observed = {0: 0} | {slitIndex: slitIndex - 1 for slitIndex in range(1, 9)}
    table = _line_table(SLIT_CENTRE - PINHOLE_SPACING - 1.0, observed)

    # ACT
    _run(mapper, table)

    # ASSERT
    messages = [message for level, message in log.messages if level == "info" and "multi-pinhole shift" in message]
    assert len(messages) == 1
    assert "whole-spacing part -1 x 4.00 = -4.00 px" in messages[0]
    assert "remainder -1.00 px" in messages[0]


def _add_spot(image: np.ndarray, x: float, y: float) -> None:
    """Add one Gaussian pinhole spot centred on ``x``, ``y`` to ``image`` in place."""
    yy, xx = np.mgrid[0 : image.shape[0], 0 : image.shape[1]]
    image += PINHOLE_FLUX * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * SPOT_SIGMA**2))


def _pinhole_frame(
    firstSpotIndex: int, lastSpotIndex: int, offset: float, rows: np.ndarray = FRAME_ROWS
) -> np.ma.MaskedArray:
    """Return a masked frame with a Gaussian spot for each real pinhole of every row of arc-line sets.

    Pinhole ``j`` sits at ``SLIT_CENTRE + (j - 4) * PINHOLE_SPACING + offset``, for ``j`` from ``firstSpotIndex`` to
    ``lastSpotIndex``.
    """
    rng = np.random.default_rng(1735)
    image = rng.normal(0.0, 4.0, size=FRAME_SHAPE)
    for y in rows:
        for j in range(firstSpotIndex, lastSpotIndex + 1):
            _add_spot(image, SLIT_CENTRE + (j - 4) * PINHOLE_SPACING + offset, y)
    return np.ma.array(image, mask=np.zeros(image.shape, dtype=bool))


def _latched_top_table() -> pd.DataFrame:
    """Return sets whose predictions are centred on the slit and whose top window latched onto pinhole 7."""
    observed = {slitIndex: slitIndex for slitIndex in range(8)} | {8: 7}
    return _line_table(SLIT_CENTRE, observed, lineRows=FRAME_ROWS, observedOffset=-REMAINDER)


def test_nine_visible_pinholes_one_spacing_off_the_prediction_move_the_order_a_whole_spacing(
    log: object, tmp_path: Path
) -> None:
    """All 9 real pinholes lie one spacing low of the slit-centred prediction, so the top window had no pinhole."""
    # ARRANGE
    frame = _pinhole_frame(0, 8, offset=-PINHOLE_SPACING - REMAINDER)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    table = _latched_top_table()

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(
        result["detector_x_shifted"], table["detector_x_shifted"] - PINHOLE_SPACING - REMAINDER, atol=0.3
    )
    assert any("pinhole visible one spacing beyond the opposite end" in message for _, message in log.messages)
    assert not [message for level, message in log.messages if level == "warning"]


def test_sets_without_a_detection_at_the_opposite_end_are_left_out_of_the_pinhole_search(
    log: object, tmp_path: Path
) -> None:
    """Two sets lack their pinhole 0 detection; the search looks at the other sets and still finds the shift."""
    # ARRANGE
    frame = _pinhole_frame(0, 8, offset=-PINHOLE_SPACING - REMAINDER)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    table = _latched_top_table()
    undetected = table["wavelength"].isin(table["wavelength"].unique()[:2]) & (table["slit_index"] == 0)
    table = table.loc[~undetected].reset_index(drop=True)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(
        result["detector_x_shifted"], table["detector_x_shifted"] - PINHOLE_SPACING - REMAINDER, atol=0.3
    )


def test_eight_visible_pinholes_keep_the_slit_centred_identification(log: object, tmp_path: Path) -> None:
    """The top pinhole is really absent: no pinhole shows beyond the opposite end, so no whole spacing is applied."""
    # ARRANGE
    frame = _pinhole_frame(0, 7, offset=-REMAINDER)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    table = _latched_top_table()

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] - REMAINDER)
    assert any("from the slit-centred pinhole identification" in message for _, message in log.messages)


def test_all_sets_hit_but_fewer_than_the_minimum_judged_do_not_move_the_order_a_whole_spacing(
    log: object, tmp_path: Path
) -> None:
    """A pinhole beyond the end in every one of too few sets is not enough evidence for a whole-spacing shift."""
    # ARRANGE
    rows = FRAME_ROWS[: PINHOLE_BEYOND_END_MIN_SETS - 1]
    assert len(rows) * 9 >= 50, "the shift step needs at least 50 detections"
    frame = _pinhole_frame(0, 8, offset=-PINHOLE_SPACING - REMAINDER, rows=rows)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    observed = {slitIndex: slitIndex for slitIndex in range(8)} | {8: 7}
    table = _line_table(SLIT_CENTRE, observed, lineRows=rows, observedOffset=-REMAINDER)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] - REMAINDER)
    assert any("from the slit-centred pinhole identification" in message for _, message in log.messages)


def test_a_pinhole_beyond_the_slit_edge_of_the_order_is_not_counted(log: object, tmp_path: Path) -> None:
    """A source one spacing beyond the end but 3.7 spacings outside the order's slit edges is not a pinhole of it."""
    # ARRANGE
    frame = _pinhole_frame(0, 8, offset=-PINHOLE_SPACING - REMAINDER)
    narrowSlit = _order_table(tmp_path, lowEdge=NARROW_SLIT_LOW_EDGE, highEdge=NARROW_SLIT_HIGH_EDGE)
    mapper = _mapper(log, narrowSlit, frame)
    table = _latched_top_table()

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] - REMAINDER)
    assert any("from the slit-centred pinhole identification" in message for _, message in log.messages)


def test_chance_sources_beyond_the_end_in_a_minority_of_sets_are_not_a_pinhole(log: object, tmp_path: Path) -> None:
    """Bright blobs at the probe position in 3 of 16 sets, as cosmic rays would leave, do not move the order."""
    # ARRANGE
    frame = _pinhole_frame(0, 7, offset=-REMAINDER)
    for y in FRAME_ROWS[:3]:
        _add_spot(frame.data, BEYOND_END_PROBE_X, y)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    table = _latched_top_table()

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] - REMAINDER)
    assert any("from the slit-centred pinhole identification" in message for _, message in log.messages)


TWO_SPACING_TEXT = "spacings from the interior pinholes"
FAR_END_TEXT = "that pinhole is missing or the order is still offset"


def _warnings(log: object) -> list[str]:
    """Return the warning messages logged so far."""
    return [message for level, message in log.messages if level == "warning"]


@pytest.mark.parametrize(
    ("predictedCentre", "observedOffset", "appliedText", "appliedShift"),
    [
        (
            SLIT_CENTRE,
            2 * PINHOLE_SPACING,
            "whole-spacing part applied is +0 x 4.00 = 0.00 px from the slit-centred",
            0.0,
        ),
        (
            SLIT_CENTRE - PINHOLE_SPACING,
            PINHOLE_SPACING,
            "whole-spacing part applied is -1 x 4.00 = -4.00 px from the slit-centred",
            -PINHOLE_SPACING,
        ),
    ],
    ids=["slit-centred-step-applies-nothing", "slit-centred-step-applies-one-spacing"],
)
def test_an_end_two_spacings_off_the_interior_says_which_shift_was_applied(
    log: object, tmp_path: Path, predictedCentre: float, observedOffset: float, appliedText: str, appliedShift: float
) -> None:
    """The low end window is 2 spacings out: the warning gives the whole-spacing part that was then applied."""
    # ARRANGE
    frame = _pinhole_frame(0, 8, offset=observedOffset)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    observed = {0: 0, 1: 0} | {slitIndex: slitIndex - 2 for slitIndex in range(2, 9)}
    table = _line_table(predictedCentre, observed, lineRows=FRAME_ROWS, observedOffset=observedOffset)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert TWO_SPACING_TEXT in warnings[0]
    assert FAR_END_TEXT not in warnings[0]
    assert appliedText in warnings[0]
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"] - appliedShift)


def test_a_far_end_without_a_pinhole_after_the_whole_spacing_shift_is_reported(log: object, tmp_path: Path) -> None:
    """Only 8 real pinholes: after the one-spacing shift the high end of the pattern has no pinhole, so warn."""
    # ARRANGE
    frame = _pinhole_frame(0, 7, offset=-PINHOLE_SPACING - REMAINDER)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    table = _latched_top_table()

    # ACT
    _run(mapper, table)

    # ASSERT
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert FAR_END_TEXT in warnings[0]
    assert TWO_SPACING_TEXT not in warnings[0]


def test_a_latched_end_with_too_few_sets_is_not_checked_after_the_whole_spacing_shift(
    log: object, tmp_path: Path
) -> None:
    """Only 9 sets have a detection at the latched end: no far-end warning, and a debug line says it was skipped."""
    # ARRANGE
    frame = _pinhole_frame(0, 7, offset=-PINHOLE_SPACING - REMAINDER)
    mapper = _mapper(log, _order_table(tmp_path), frame)
    table = _latched_top_table()
    firstSets = table["wavelength"].unique()[: len(FRAME_ROWS) - (PINHOLE_BEYOND_END_MIN_SETS - 1)]
    table = table.loc[~(table["wavelength"].isin(firstSets) & (table["slit_index"] == 8))].reset_index(drop=True)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    assert not _warnings(log)
    assert any(level == "debug" and "is skipped" in message for level, message in log.messages)
    np.testing.assert_allclose(
        result["detector_x_shifted"], table["detector_x_shifted"] - PINHOLE_SPACING - REMAINDER, atol=0.3
    )


def test_the_order_table_is_read_once_however_many_probes_look_at_the_image(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The opposite-end search and the far-end check both probe the image; the order table is read once."""
    # ARRANGE
    module = sys.modules[create_dispersion_map.__module__]
    reads: list[str] = []
    realUnpack = module.unpack_order_table

    def counting_unpack(**kwargs: object) -> object:
        reads.append(str(kwargs["orderTablePath"]))
        return realUnpack(**kwargs)

    monkeypatch.setattr(module, "unpack_order_table", counting_unpack)
    frame = _pinhole_frame(0, 8, offset=-PINHOLE_SPACING - REMAINDER)
    mapper = _mapper(log, _order_table(tmp_path), frame)

    # ACT
    _run(mapper, _latched_top_table())

    # ASSERT
    assert len(reads) == 1


def test_order_missing_from_the_order_table_is_not_shifted(log: object, tmp_path: Path) -> None:
    """Without slit edges for the order the shift cannot be decided, so none is applied and a warning says so."""
    # ARRANGE
    mapper = _mapper(log, _order_table(tmp_path, orders=(3,)))
    observed = {0: 0} | {slitIndex: slitIndex - 1 for slitIndex in range(1, 9)}
    table = _line_table(SLIT_CENTRE - PINHOLE_SPACING, observed)

    # ACT
    result = _run(mapper, table)

    # ASSERT
    np.testing.assert_allclose(result["detector_x_shifted"], table["detector_x_shifted"])
    assert any(level == "warning" and "No slit edges match" in message for level, message in log.messages)


def _bootstrapping_mapper(log: object, pinholeCount: int, pinholeCountRange: list[float]) -> create_dispersion_map:
    """Return a spatial-solution mapper whose fit finds ``pinholeCount`` pinholes for a single arc line."""
    mapper = create_dispersion_map.__new__(create_dispersion_map)
    mapper.log = log
    lines = pd.DataFrame(
        {
            "wavelength": [500.0] * pinholeCount,
            "order": [2] * pinholeCount,
            "slit_index": range(pinholeCount),
            "slit_position": np.linspace(-1.0, 1.0, pinholeCount),
            "observed_x": np.arange(100.0, 100.0 + pinholeCount),
        }
    )
    mapper.lineDetectionTable = lines.copy()
    mapper.recipeName = "soxs-spat-solution"
    mapper.recipeSettings = {
        "mph_line_set_min": 1,
        "qc-acceptable-ranges": {"pinhole-count-min": pinholeCountRange},
    }
    mapper.settings = {"tune-pipeline": False}
    mapper.firstGuessMap = "first-guess.fits"
    mapper.orderTable = "orders.fits"
    mapper.create2DMap = True
    mapper.qc = pd.DataFrame()
    mapper.products = pd.DataFrame()
    mapper.dateObs = "2025-01-02T03:04:05"
    mapper.inst = "SOXS"
    mapper._initialize_recipe_settings = lambda: (False, False, 1, 1, 1)
    mapper.get_predicted_line_list = lambda: lines.copy()
    mapper._prepare_pinhole_frame = lambda: "prepared-frame"
    mapper._generate_quicklook_and_debug_plots = lambda *args: None
    mapper._clip_on_measured_line_metrics = lambda table: table
    mapper.fit_polynomials = lambda **kwargs: (np.array([1.0]), np.array([2.0]), lines.copy(), lines.iloc[0:0].copy())
    mapper._get_output_filenames = lambda: ("good.fits", "missing.fits")
    mapper._write_qc_metrics = lambda *args: None
    mapper._write_line_list_qc = lambda *args: None
    mapper._prepare_line_list_columns = lambda *args: (lines.copy(), lines.copy())
    mapper._write_fitted_lines_file = lambda *args: None
    mapper._write_missing_lines_file = lambda *args: None
    mapper.write_map_to_file = lambda *args: "map.fits"
    mapper.map_to_image = lambda **kwargs: "map_IMAGE.fits"
    mapper._create_dispersion_map_qc_plot = lambda **kwargs: ["residuals.pdf"]
    return mapper


@pytest.mark.parametrize(
    ("pinholeCount", "pinholeCountRange", "expectedImage"),
    [
        (8, [7.5, 9.5], "map_IMAGE.fits"),
        (8, [8.5, 9.5], None),
        (9, [8.5, 9.5], "map_IMAGE.fits"),
        (10, [8.5, 9.5], None),
    ],
)
def test_2d_map_image_is_written_when_the_pinhole_count_passes_qc(
    log: object,
    pinholeCount: int,
    pinholeCountRange: list[float],
    expectedImage: str | None,
) -> None:
    """An 8-pinhole solution writes its 2D image where the QC range accepts 8 pinholes, and only there."""
    # ARRANGE
    mapper = _bootstrapping_mapper(log, pinholeCount, pinholeCountRange)

    # ACT
    result = mapper.get()

    # ASSERT
    assert mapper.minpin == pinholeCount
    assert result[1] == expectedImage


@pytest.mark.parametrize(
    ("arm", "pinholeCount", "expected"),
    [("vis", 8, True), ("vis", 7, False), ("nir", 8, False), ("nir", 9, True)],
)
def test_shipped_settings_accept_8_vis_pinholes_and_still_require_9_nir_pinholes(
    log: object, arm: str, pinholeCount: int, expected: bool
) -> None:
    """The shipped SOXS spatial-solution QC ranges let an 8-pinhole VIS solution through and nothing less."""
    # ARRANGE
    import yaml

    import soxspipe

    with open(Path(soxspipe.__file__).parent / "soxs_default_settings.yaml") as stream:
        settings = yaml.safe_load(stream)
    mapper = create_dispersion_map.__new__(create_dispersion_map)
    mapper.recipeSettings = settings["soxs-spat-solution"][arm]
    mapper.minpin = pinholeCount

    # ACT
    passes = mapper._pinhole_count_passes_qc()

    # ASSERT
    assert passes is expected
