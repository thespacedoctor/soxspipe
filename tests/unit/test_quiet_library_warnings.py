"""Expected library warnings stay out of the logs, and the two real ones are fixed (DY-2409).

`logging.captureWarnings` sends every Python warning to the console and the recipe logs, so a warning that is
expected and harmless buries the ones that matter. The filters here match an exact message and category: a
similar warning with another message still reaches the log.

pytest wraps every test in `warnings.catch_warnings`, and each test that installs a filter checks it from inside its
own `warnings.catch_warnings` block, so none of these filters leaks into another test.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy.io.fits.verify import VerifyWarning
from astropy.modeling import fitting, models
from astropy.utils.exceptions import AstropyUserWarning

from soxspipe.commonutils import toolkit
from soxspipe.recipes import base_recipe
from tests.unit.test_base_recipe_init_characterization import _isolate, _settings
from tests.unit.test_cli_dispatch import CliLogger, _run_cli
from tests.unit.test_continuum_detection import _detector, continuumModule
from tests.unit.test_create_dispersion_map_helpers import _mapper
from tests.unit.test_create_dispersion_map_public import _qc_plot_mapper, _run_qc_plot

pytestmark = pytest.mark.unit

CARD_TOO_LONG = "Card is too long, comment will be truncated."
SIGMA_CLIP_NAN = "Input data contains invalid values (NaNs or infs), which were automatically clipped."
FIT_UNSUCCESSFUL = "The fit may be unsuccessful; check fit_info['message'] for more information."
COLLAPSED_LAYOUT = "constrained_layout not applied because axes sizes collapsed to zero."


def _raised_by(
    call: Callable[[], object],
    category: type[Warning],
    message: str,
    *,
    setup: Callable[[], object] | None = None,
) -> list[warnings.WarningMessage]:
    """Return the warnings of ``category`` with ``message`` that ``call`` raised.

    Every warning is switched on first, then ``setup`` installs its filters on top, all inside a
    `warnings.catch_warnings` block that restores the process-wide filters on exit.
    """
    with warnings.catch_warnings(record=True) as raised:
        warnings.simplefilter("always")
        if setup is not None:
            setup()
        call()
    return [w for w in raised if issubclass(w.category, category) and message in str(w.message)]


def _failing_fit() -> object:
    """Run a Gaussian fit that LevMar cannot finish, so astropy warns that it may be unsuccessful."""
    x = np.arange(30.0)
    y = np.where(x == 2, 1e3, 0.0)
    return fitting.LevMarLSQFitter()(models.Gaussian1D(1.0, 10, 2.0), x, y)


# ---------------------------------------------------------------------------
# ignore_expected_library_warnings
# ---------------------------------------------------------------------------


def test_a_truncated_header_comment_raises_no_warning_once_the_expected_warnings_are_ignored() -> None:
    # ACT
    raised = _raised_by(
        lambda: warnings.warn(CARD_TOO_LONG, VerifyWarning, stacklevel=2),
        VerifyWarning,
        "Card is too long",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert raised == []


def test_a_sigma_clip_nan_notice_raises_no_warning_once_the_expected_warnings_are_ignored() -> None:
    # ARRANGE
    from astropy.stats import sigma_clip

    # ACT
    raised = _raised_by(
        lambda: sigma_clip(np.array([1.0, np.nan, 2.0, 3.0])),
        AstropyUserWarning,
        "Input data contains invalid",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert raised == []


def test_a_sigma_clip_nan_notice_reaches_the_log_before_the_expected_warnings_are_ignored() -> None:
    # ARRANGE
    from astropy.stats import sigma_clip

    # ACT
    raised = _raised_by(
        lambda: sigma_clip(np.array([1.0, np.nan, 2.0, 3.0])), AstropyUserWarning, "Input data contains invalid"
    )

    # ASSERT
    assert len(raised) == 1


def test_an_unlisted_verify_warning_still_reaches_the_log() -> None:
    # ACT
    raised = _raised_by(
        lambda: warnings.warn("Keyword name is too long", VerifyWarning, stacklevel=2),
        VerifyWarning,
        "too long",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


def test_an_unlisted_astropy_user_warning_still_reaches_the_log() -> None:
    # ACT
    raised = _raised_by(
        lambda: warnings.warn("Input data contains a surprise", AstropyUserWarning, stacklevel=2),
        AstropyUserWarning,
        "surprise",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


def test_the_card_message_is_not_hidden_under_a_different_category() -> None:
    # ACT
    raised = _raised_by(
        lambda: warnings.warn(CARD_TOO_LONG, UserWarning, stacklevel=2),
        UserWarning,
        "Card is too long",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


def test_ignoring_the_expected_warnings_twice_adds_each_filter_once() -> None:
    # ARRANGE
    with warnings.catch_warnings():
        toolkit.ignore_expected_library_warnings()
        filterCount = len(warnings.filters)

        # ACT
        toolkit.ignore_expected_library_warnings()

        # ASSERT
        assert len(warnings.filters) == filterCount


# ---------------------------------------------------------------------------
# where the filters are installed
# ---------------------------------------------------------------------------


def test_building_a_recipe_ignores_the_expected_library_warnings(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    # ARRANGE
    _isolate(monkeypatch, tmp_path, productPath=tmp_path / "reduced" / "synthetic.fits")

    def build_recipe() -> None:
        base_recipe(
            log=log,
            settings=_settings(tmp_path),
            inputFrames=str(tmp_path / "synthetic.sof"),
            recipeName="soxs-mbias",
            turnOffMP=True,
        )

    # ACT
    raised = _raised_by(
        lambda: warnings.warn(CARD_TOO_LONG, VerifyWarning, stacklevel=2),
        VerifyWarning,
        "Card is too long",
        setup=build_recipe,
    )

    # ASSERT
    assert raised == []


def test_the_command_line_ignores_the_expected_library_warnings(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    def run_command() -> None:
        _run_cli(monkeypatch, tmp_path, ["soxspipe", "session", "list"], logger=CliLogger())

    # ACT
    raised = _raised_by(
        lambda: warnings.warn(CARD_TOO_LONG, VerifyWarning, stacklevel=2),
        VerifyWarning,
        "Card is too long",
        setup=run_command,
    )

    # ASSERT
    assert raised == []


# ---------------------------------------------------------------------------
# quiet_least_squares_fit
# ---------------------------------------------------------------------------


def test_a_levmar_fit_that_may_be_unsuccessful_raises_no_warning_inside_the_quiet_fit() -> None:
    # ARRANGE
    def quiet_fit() -> None:
        with toolkit.quiet_least_squares_fit():
            _failing_fit()

    # ACT
    loud = _raised_by(_failing_fit, AstropyUserWarning, "The fit may be unsuccessful")
    quiet = _raised_by(quiet_fit, AstropyUserWarning, "The fit may be unsuccessful")

    # ASSERT
    assert len(loud) == 1
    assert quiet == []


def test_a_numpy_matmul_overflow_raises_no_warning_inside_the_quiet_fit() -> None:
    # ARRANGE
    huge = np.array([[1e200, 1e200], [0.0, 1e200]])

    def quiet_matmul() -> None:
        with toolkit.quiet_least_squares_fit():
            huge @ huge.T

    # ACT
    loud = _raised_by(lambda: huge @ huge.T, RuntimeWarning, "overflow encountered in matmul")
    quiet = _raised_by(quiet_matmul, RuntimeWarning, "overflow encountered in matmul")

    # ASSERT
    assert len(loud) == 1
    assert quiet == []


def test_another_runtime_warning_still_reaches_the_log_inside_the_quiet_fit() -> None:
    # ARRANGE
    def quiet_sqrt() -> None:
        with toolkit.quiet_least_squares_fit():
            warnings.warn("something else broke", RuntimeWarning, stacklevel=2)

    # ACT
    raised = _raised_by(quiet_sqrt, RuntimeWarning, "something else broke")

    # ASSERT
    assert len(raised) == 1


def test_the_quiet_fit_restores_the_warning_filters_and_the_numpy_error_state_on_exit() -> None:
    # ARRANGE
    filtersBefore = list(warnings.filters)
    errorStateBefore = np.geterr()

    # ACT
    with toolkit.quiet_least_squares_fit():
        pass

    # ASSERT
    assert warnings.filters == filtersBefore
    assert np.geterr() == errorStateBefore


# ---------------------------------------------------------------------------
# the fit sites
# ---------------------------------------------------------------------------


def test_the_continuum_slice_fit_raises_no_unsuccessful_fit_warning_for_a_slice_levmar_cannot_fit(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    detector = _detector(log)
    detector.traceFrame = object()
    detector.sliceWidth = 3
    detector.detectorParams = {"dispersion-axis": "x"}
    detector.debug = False
    detector.peakSigmaLimit = 3.0
    detector.sliceAxis = "x"
    detector.sliceAntiAxis = "y"
    positions = np.arange(21, dtype=float)
    unfittable = np.where(positions == 10, 1e3, 0.0)
    monkeypatch.setattr(continuumModule, "cut_image_slice", lambda **kwargs: (unfittable, 5, 7.0))
    pixels = pd.DataFrame({"fit_x": [15.0], "fit_y": [7.0]})

    # ACT
    raised = _raised_by(
        lambda: detector.fit_1d_gaussian_to_slices(pixels, sliceLength=21),
        AstropyUserWarning,
        "The fit may be unsuccessful",
    )

    # ASSERT
    assert raised == []


def test_the_arc_line_resolution_fit_raises_no_unsuccessful_fit_warning_for_a_line_levmar_cannot_fit(
    log: object,
) -> None:
    # ARRANGE
    mapper = _mapper(log)
    mapper.slitWidth = 1.0
    mapper.detectorParams = {"dispersion-axis": "x"}
    row = pd.Series({"observed_x": 20.0, "observed_y": 20.0, "wavelength": 500.0, "pixelScaleNm": 0.01})
    image = np.zeros((40, 40))
    image[20, 20] = 1e3
    mapper.arcFrame = SimpleNamespace(data=image)

    # ACT
    raised = _raised_by(
        lambda: mapper._calculate_resolution_on_slit(row),
        AstropyUserWarning,
        "The fit may be unsuccessful",
    )

    # ASSERT
    assert raised == []


# ---------------------------------------------------------------------------
# quiet_constrained_layout
# ---------------------------------------------------------------------------


def _collapsed_figure() -> plt.Figure:
    figure = plt.figure(figsize=(0.3, 0.3), layout="constrained")
    axes = figure.add_subplot(111)
    axes.set_title("a long title")
    axes.set_xlabel("x label")
    axes.set_ylabel("y label")
    return figure


def test_a_collapsed_constrained_layout_raises_no_warning_inside_the_quiet_layout() -> None:
    # ARRANGE
    figure = _collapsed_figure()

    def quiet_draw() -> None:
        with toolkit.quiet_constrained_layout():
            figure.draw_without_rendering()

    # ACT
    try:
        loud = _raised_by(figure.draw_without_rendering, UserWarning, "constrained_layout not applied")
        quiet = _raised_by(quiet_draw, UserWarning, "constrained_layout not applied")
    finally:
        plt.close(figure)

    # ASSERT
    assert len(loud) >= 1
    assert quiet == []


def test_measuring_the_lowest_axes_of_a_collapsed_figure_raises_no_layout_warning() -> None:
    # ARRANGE
    figure = _collapsed_figure()

    # ACT
    try:
        raised = _raised_by(
            lambda: toolkit._lowest_axes_bottom_in_figure_fraction(figure),
            UserWarning,
            "constrained_layout not applied",
        )
    finally:
        plt.close(figure)

    # ASSERT
    assert raised == []


# ---------------------------------------------------------------------------
# the matplotlib "w-" format string
# ---------------------------------------------------------------------------

FMT_CLASH = "is redundantly defined by the 'color' keyword argument and the fmt string"


def _grid_line_table() -> pd.DataFrame:
    return pd.DataFrame(
        {"line": [0, 0, 1, 1, 2], "fit_x": [1.0, 2.0, 1.0, 2.0, 3.0], "fit_y": [3.0, 4.0, 5.0, 6.0, 7.0]}
    )


@pytest.mark.parametrize("inst", ["SOXS", "XSHOOTER"])
def test_the_quicklook_grid_lines_are_drawn_without_a_colour_clash_between_the_format_string_and_the_keyword(
    inst: str,
) -> None:
    # ARRANGE
    figure, ax = plt.subplots()

    # ACT
    try:
        raised = _raised_by(
            lambda: toolkit._draw_dispersion_grid_lines(ax, _grid_line_table(), inst), UserWarning, FMT_CLASH
        )
        colours = {line.get_color() for line in ax.lines}
    finally:
        plt.close(figure)

    # ASSERT
    assert raised == []
    assert colours == {"black"}


def test_the_dispersion_map_qc_plot_draws_its_grid_lines_without_a_colour_clash(
    log: object,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    mapper = _qc_plot_mapper(log, tmp_path)
    mapper.pinholeFrame = None
    mapper.uniqueSlitPos = []
    gridLines = _grid_line_table().rename(columns={"fit_x": "fit_x", "fit_y": "fit_y"})
    monkeypatch.setattr(toolkit, "create_dispersion_solution_grid_lines_for_plot", lambda **_: (gridLines, None))

    # ACT
    try:
        raised = _raised_by(
            lambda: _run_qc_plot(mapper, dispMap="map.fits", dispMapImage="map_image.fits"), UserWarning, FMT_CLASH
        )
    finally:
        plt.close("all")

    # ASSERT
    assert raised == []
