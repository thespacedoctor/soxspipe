"""DY-1618: the first-iteration multi-pinhole shift picks the slit-centred pinhole identification (SOXS VIS).

When the pinhole at one end of the slit is missing from the multi-pinhole frame, its detection window latches onto
the neighbouring pinhole. The end-to-end shift test then sees a one-pinhole jump at that end only. That is also what
a real one-pinhole offset between the first-guess map and the frame looks like, so the shift is chosen between no
shift and one pinhole spacing either way by centring the 9-pinhole pattern between the order's slit edges.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.create_dispersion_map import create_dispersion_map
from tests.factories import order_table_fits

pytestmark = pytest.mark.unit

SLIT_LOW_EDGE = 100.0
SLIT_HIGH_EDGE = 140.0
SLIT_CENTRE = 0.5 * (SLIT_LOW_EDGE + SLIT_HIGH_EDGE)
PINHOLE_SPACING = 4.0
LINE_ROWS = np.arange(10.0, 70.0, 2.0)


def _order_table(tmp_path: Path, orders: tuple[int, ...] = (2, 3)) -> str:
    """Write an order table whose slit runs from x=100 to x=140 at every row of every order."""
    polynomials = pd.DataFrame(
        [
            {
                "degorder_cent": 0,
                "degy_cent": 0,
                "cent_00": SLIT_CENTRE,
                "degorder_edgeup": 0,
                "degy_edgeup": 0,
                "edgeup_00": SLIT_HIGH_EDGE,
                "degorder_edgelow": 0,
                "degy_edgelow": 0,
                "edgelow_00": SLIT_LOW_EDGE,
            }
        ]
    )
    metadata = pd.DataFrame({"order": list(orders), "ymin": 0.0, "ymax": 100.0})
    return str(order_table_fits(tmp_path / "orders.fits", polynomials=polynomials, metadata=metadata))


def _mapper(log: object, orderTable: str) -> create_dispersion_map:
    """Return a SOXS VIS multi-pinhole mapper that reads only what the shift step needs."""
    mapper = create_dispersion_map.__new__(create_dispersion_map)
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


def _line_table(predictedCentre: float, observedSlitIndex: dict[int, int], order: int = 2) -> pd.DataFrame:
    """Return 9-pinhole line sets for ``order`` whose predictions are centred on ``predictedCentre``.

    ``observedSlitIndex`` maps each predicted slit index to the index of the real pinhole its window detected.
    """
    rows = []
    for wavelengthIndex, y in enumerate(LINE_ROWS):
        for slitIndex in range(9):
            predictedX = predictedCentre + (slitIndex - 4) * PINHOLE_SPACING
            observedX = SLIT_CENTRE + (observedSlitIndex[slitIndex] - 4) * PINHOLE_SPACING
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
