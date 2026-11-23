"""Expected library warnings stay out of the logs, and the two real ones are fixed (DY-2409).

`logging.captureWarnings` sends every Python warning to the console and the recipe logs, so a warning that is
expected and harmless buries the ones that matter. Each filter here matches the full message and the category, and the
astropy and SciPy ones also the module that raises the warning: a similar warning with another message, or the same
message from another module, still reaches the log. The one deliberate prefix match is astropy's "The fit may be
unsuccessful", which carries a variable tail.

The process-wide filters are appended after any filter that is already set, so `-W error` and a user's own filters keep
precedence. Every test that checks a filter first clears the filter list inside its own `warnings.catch_warnings` block
and restores it on exit, so none of these filters leaks into another test.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from contextlib import contextmanager
from importlib import import_module
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.modeling import fitting, models
from astropy.stats import sigma_clip
from astropy.utils.exceptions import AstropyUserWarning
from scipy.optimize import leastsq

from soxspipe.commonutils import toolkit
from soxspipe.recipes import base_recipe
from tests.unit.test_base_recipe_init_characterization import _isolate, _settings
from tests.unit.test_cli_dispatch import CliLogger, _run_cli
from tests.unit.test_continuum_detection import _detector, continuumModule
from tests.unit.test_create_dispersion_map_helpers import _mapper
from tests.unit.test_create_dispersion_map_public import _qc_plot_mapper, _run_qc_plot

pytestmark = pytest.mark.unit

CARD_TOO_LONG = "Card is too long, comment will be truncated."
COLLAPSED_LAYOUT = (
    "constrained_layout not applied because axes sizes collapsed to zero.  "
    "Try making figure larger or Axes decorations smaller."
)
MATMUL_MODULE = "scipy.optimize._minpack_py"


def _raised_by(
    call: Callable[[], object],
    category: type[Warning],
    message: str,
    *,
    setup: Callable[[], object] | None = None,
) -> list[warnings.WarningMessage]:
    """Return the warnings of ``category`` with ``message`` that ``call`` raised.

    The filter list is cleared, ``setup`` installs its filters, and an appended "always" filter then records whatever
    the installed filters let through, all inside a `warnings.catch_warnings` block that restores the process-wide
    filters on exit.
    """
    with warnings.catch_warnings(record=True) as raised:
        warnings.resetwarnings()
        if setup is not None:
            setup()
        warnings.simplefilter("always", append=True)
        call()
    return [w for w in raised if issubclass(w.category, category) and message in str(w.message)]


def _failing_fit() -> object:
    """Run a Gaussian fit that LevMar cannot finish, so astropy warns that it may be unsuccessful."""
    x = np.arange(30.0)
    y = np.where(x == 2, 1e3, 0.0)
    return fitting.LevMarLSQFitter()(models.Gaussian1D(1.0, 10, 2.0), x, y)


def _overflowing_covariance_fit() -> object:
    """Run a SciPy least-squares fit whose covariance step overflows, as a near-singular LevMar Jacobian does."""
    x = np.arange(10.0)
    tiny = 1e-200
    return leastsq(
        lambda p: p[0] * x + tiny * p[1] * x - x,
        [0.0, 0.0],
        Dfun=lambda p: np.array([x, tiny * x]).T,
        full_output=True,
    )


def _truncated_header_card() -> object:
    """Write a header card whose comment is too long, so astropy cuts it and warns."""
    header = fits.Header()
    header["HIERARCH ESO QC TEST"] = (1, "x" * 100)
    return header.tostring()


def _sigma_clip_with_a_nan() -> object:
    """Sigma-clip data that holds a NaN, so astropy warns that it clipped invalid values."""
    return sigma_clip(np.array([1.0, np.nan, 2.0, 3.0]))


def _warn_from(message: str, category: type[Warning], module: str) -> None:
    """Raise a warning the way a library does, attributed to ``module``."""
    warnings.warn_explicit(message, category, filename=f"{module}.py", lineno=1, module=module, registry={})


# ---------------------------------------------------------------------------
# IGNORE_EXPECTED_LIBRARY_WARNINGS
# ---------------------------------------------------------------------------


def test_a_truncated_header_comment_raises_no_warning_once_the_expected_warnings_are_ignored() -> None:
    # ARRANGE
    loud = _raised_by(_truncated_header_card, VerifyWarning, "Card is too long")

    # ACT
    quiet = _raised_by(
        _truncated_header_card,
        VerifyWarning,
        "Card is too long",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(loud) >= 1
    assert quiet == []


def test_a_sigma_clip_nan_notice_raises_no_warning_once_the_expected_warnings_are_ignored() -> None:
    # ARRANGE
    loud = _raised_by(_sigma_clip_with_a_nan, AstropyUserWarning, "Input data contains invalid")

    # ACT
    quiet = _raised_by(
        _sigma_clip_with_a_nan,
        AstropyUserWarning,
        "Input data contains invalid",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(loud) == 1
    assert quiet == []


@pytest.mark.parametrize(
    ("message", "category", "module"),
    [
        (CARD_TOO_LONG + " Extra.", VerifyWarning, "astropy.io.fits.card"),
        (toolkit.SIGMA_CLIP_INVALID_MESSAGE + " Extra.", AstropyUserWarning, "astropy.stats.sigma_clipping"),
        ("Extra. " + toolkit.SIGMA_CLIP_INVALID_MESSAGE, AstropyUserWarning, "astropy.stats.sigma_clipping"),
    ],
)
def test_a_listed_message_with_an_unexpected_extra_text_still_reaches_the_log(
    message: str, category: type[Warning], module: str
) -> None:
    # ACT
    raised = _raised_by(
        lambda: _warn_from(message, category, module),
        category,
        "",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


@pytest.mark.parametrize(
    ("message", "category"),
    [
        (CARD_TOO_LONG, VerifyWarning),
        (toolkit.SIGMA_CLIP_INVALID_MESSAGE, AstropyUserWarning),
    ],
)
def test_a_listed_message_from_another_module_still_reaches_the_log(message: str, category: type[Warning]) -> None:
    # ACT
    raised = _raised_by(
        lambda: _warn_from(message, category, "some.other.module"),
        category,
        "",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


def test_an_unlisted_verify_warning_still_reaches_the_log() -> None:
    # ACT
    raised = _raised_by(
        lambda: _warn_from("Keyword name is too long", VerifyWarning, "astropy.io.fits.card"),
        VerifyWarning,
        "too long",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


def test_an_unlisted_astropy_user_warning_still_reaches_the_log() -> None:
    # ACT
    raised = _raised_by(
        lambda: _warn_from("Input data contains a surprise", AstropyUserWarning, "astropy.stats.sigma_clipping"),
        AstropyUserWarning,
        "surprise",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


def test_the_card_message_is_not_hidden_under_a_different_category() -> None:
    # ACT
    raised = _raised_by(
        lambda: _warn_from(CARD_TOO_LONG, UserWarning, "astropy.io.fits.card"),
        UserWarning,
        "Card is too long",
        setup=toolkit.ignore_expected_library_warnings,
    )

    # ASSERT
    assert len(raised) == 1


def test_a_filter_set_before_the_expected_warnings_are_ignored_keeps_precedence() -> None:
    # ARRANGE
    with warnings.catch_warnings():
        warnings.resetwarnings()
        warnings.simplefilter("error")

        # ACT
        toolkit.ignore_expected_library_warnings()

        # ASSERT
        with pytest.raises(VerifyWarning, match="Card is too long"):
            _truncated_header_card()


def test_ignoring_the_expected_warnings_twice_adds_each_filter_once() -> None:
    # ARRANGE
    with warnings.catch_warnings():
        toolkit.ignore_expected_library_warnings()
        filtersAfterFirstCall = list(warnings.filters)

        # ACT
        toolkit.ignore_expected_library_warnings()

        # ASSERT
        assert warnings.filters == filtersAfterFirstCall


def test_ignoring_the_expected_warnings_again_does_not_touch_the_filter_list() -> None:
    # ARRANGE
    # CPYTHON 3.12 DETAIL: `filterwarnings` CALLS THE MODULE-LEVEL `_filters_mutated`. IF A LATER PYTHON RENAMES IT,
    # THE PATCH BELOW WOULD DO NOTHING AND THE TEST WOULD PASS VACUOUSLY, SO FAIL LOUDLY INSTEAD
    assert callable(getattr(warnings, "_filters_mutated", None))
    mutations: list[str] = []
    realMutated = warnings._filters_mutated
    with warnings.catch_warnings():
        toolkit.ignore_expected_library_warnings()
        warnings._filters_mutated = lambda: mutations.append("mutated")
        try:
            # ACT
            toolkit.ignore_expected_library_warnings()
        finally:
            warnings._filters_mutated = realMutated

    # ASSERT
    assert mutations == []


# ---------------------------------------------------------------------------
# WHERE THE FILTERS ARE INSTALLED
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
        _truncated_header_card,
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
        _truncated_header_card,
        VerifyWarning,
        "Card is too long",
        setup=run_command,
    )

    # ASSERT
    assert raised == []


# ---------------------------------------------------------------------------
# QUIET_LEAST_SQUARES_FIT
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


def test_the_scipy_covariance_overflow_raises_no_warning_inside_the_quiet_fit() -> None:
    # ARRANGE
    def quiet_fit() -> None:
        with toolkit.quiet_least_squares_fit():
            _overflowing_covariance_fit()

    # ACT
    loud = _raised_by(_overflowing_covariance_fit, RuntimeWarning, "overflow encountered in matmul")
    quiet = _raised_by(quiet_fit, RuntimeWarning, "overflow encountered in matmul")

    # ASSERT
    assert len(loud) == 1
    assert quiet == []


@pytest.mark.parametrize("message", ["invalid value encountered in matmul", "overflow encountered in matmul"])
def test_each_listed_matmul_message_from_scipy_is_hidden_inside_the_quiet_fit(message: str) -> None:
    # ARRANGE
    def quiet_warn() -> None:
        with toolkit.quiet_least_squares_fit():
            _warn_from(message, RuntimeWarning, MATMUL_MODULE)

    # ACT
    raised = _raised_by(quiet_warn, RuntimeWarning, message)

    # ASSERT
    assert raised == []


@pytest.mark.parametrize(
    ("message", "module"),
    [
        ("overflow encountered in matmul", "some.other.module"),
        ("overflow encountered in matmul in the model", MATMUL_MODULE),
        ("overflow encountered in multiply", MATMUL_MODULE),
    ],
)
def test_a_similar_runtime_warning_still_reaches_the_log_inside_the_quiet_fit(message: str, module: str) -> None:
    # ARRANGE
    def quiet_warn() -> None:
        with toolkit.quiet_least_squares_fit():
            _warn_from(message, RuntimeWarning, module)

    # ACT
    raised = _raised_by(quiet_warn, RuntimeWarning, message)

    # ASSERT
    assert len(raised) == 1


@pytest.mark.parametrize(
    ("event", "message"),
    [
        (lambda: np.log(-1.0), "invalid value encountered in log"),
        (lambda: np.float64(1e308) * 10, "overflow encountered in scalar multiply"),
    ],
)
def test_a_real_numpy_event_from_the_model_still_reaches_the_log_inside_the_quiet_fit(
    event: Callable[[], object], message: str
) -> None:
    # ARRANGE
    def quiet_event() -> None:
        with toolkit.quiet_least_squares_fit():
            event()

    # ACT
    raised = _raised_by(quiet_event, RuntimeWarning, message)

    # ASSERT
    assert len(raised) == 1


def test_another_astropy_warning_still_reaches_the_log_inside_the_quiet_fit() -> None:
    # ARRANGE
    def quiet_warn() -> None:
        with toolkit.quiet_least_squares_fit():
            _warn_from("Something else about the fit", AstropyUserWarning, "astropy.modeling.fitting")

    # ACT
    raised = _raised_by(quiet_warn, AstropyUserWarning, "Something else")

    # ASSERT
    assert len(raised) == 1


def test_the_quiet_fit_restores_the_warning_filters_on_exit() -> None:
    # ARRANGE
    filtersBefore = list(warnings.filters)

    # ACT
    with toolkit.quiet_least_squares_fit():
        pass

    # ASSERT
    assert warnings.filters == filtersBefore


def test_the_quiet_fit_restores_the_warning_filters_when_the_block_raises() -> None:
    # ARRANGE
    filtersBefore = list(warnings.filters)

    # ACT
    with pytest.raises(ValueError, match="fit failed"), toolkit.quiet_least_squares_fit():
        raise ValueError("fit failed")

    # ASSERT
    assert warnings.filters == filtersBefore


# ---------------------------------------------------------------------------
# THE FIT SITES
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
    row = pd.DataFrame({"observed_x": [20.0], "observed_y": [20.0], "wavelength": [500.0], "pixelScaleNm": [0.01]})
    image = np.zeros((40, 40))
    image[20, 20] = 1e3
    mapper.arcFrame = SimpleNamespace(data=image)

    # ACT
    raised = _raised_by(
        lambda: mapper._measure_slit_resolutions(row),
        AstropyUserWarning,
        "The fit may be unsuccessful",
    )

    # ASSERT
    assert raised == []


def _continuum_detector(log: object, monkeypatch: pytest.MonkeyPatch, sliceCount: int) -> tuple[object, pd.DataFrame]:
    """Return a detector and a pixel table of ``sliceCount`` slices that LevMar cannot fit."""
    detector = _detector(log)
    detector.traceFrame = object()
    detector.sliceWidth = 3
    detector.detectorParams = {"dispersion-axis": "x"}
    detector.debug = False
    detector.peakSigmaLimit = 3.0
    detector.sliceAxis = "x"
    detector.sliceAntiAxis = "y"
    unfittable = np.where(np.arange(21) == 10, 1e3, 0.0)
    monkeypatch.setattr(continuumModule, "cut_image_slice", lambda **kwargs: (unfittable, 5, 7.0))
    pixels = pd.DataFrame({"fit_x": [15.0] * sliceCount, "fit_y": [7.0] * sliceCount})
    return detector, pixels


def _count_quiet_fit_entries(module: object, monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Replace ``module.quiet_least_squares_fit`` with a counting wrapper and return the entry log."""
    entries: list[str] = []
    real = module.quiet_least_squares_fit

    @contextmanager
    def counting() -> object:
        entries.append("entered")
        with real():
            yield

    monkeypatch.setattr(module, "quiet_least_squares_fit", counting)
    return entries


def test_the_continuum_slice_loop_enters_the_quiet_fit_once_for_all_its_slices(
    log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    detector, pixels = _continuum_detector(log, monkeypatch, sliceCount=3)
    entries = _count_quiet_fit_entries(continuumModule, monkeypatch)

    # ACT
    raised = _raised_by(
        lambda: detector.fit_1d_gaussian_to_slices(pixels, sliceLength=21),
        AstropyUserWarning,
        "The fit may be unsuccessful",
    )

    # ASSERT
    assert entries == ["entered"]
    assert raised == []


def test_the_slit_resolution_of_every_arc_line_is_measured_inside_one_quiet_fit(
    log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    mapper = _mapper(log)
    mapper.slitWidth = 1.0
    mapper.detectorParams = {"dispersion-axis": "x"}
    image = np.zeros((40, 40))
    image[20, 20] = 1e3
    image[20, 25] = 1e3
    mapper.arcFrame = SimpleNamespace(data=image)
    lines = pd.DataFrame(
        {
            "observed_x": [20.0, 25.0, 20.0],
            "observed_y": [20.0, 20.0, 20.0],
            "wavelength": [500.0] * 3,
            "pixelScaleNm": [0.01] * 3,
        }
    )
    entries = _count_quiet_fit_entries(import_module("soxspipe.commonutils.create_dispersion_map"), monkeypatch)

    # ACT
    raised = _raised_by(
        lambda: mapper._measure_slit_resolutions(lines),
        AstropyUserWarning,
        "The fit may be unsuccessful",
    )

    # ASSERT
    assert entries == ["entered"]
    assert raised == []


# ---------------------------------------------------------------------------
# QUIET_CONSTRAINED_LAYOUT
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
    assert str(loud[0].message) == COLLAPSED_LAYOUT
    assert quiet == []


@pytest.mark.parametrize("message", [COLLAPSED_LAYOUT + " Extra.", "Extra. " + COLLAPSED_LAYOUT])
def test_a_layout_notice_with_an_unexpected_extra_text_still_reaches_the_log_inside_the_quiet_layout(
    message: str,
) -> None:
    # ARRANGE
    def quiet_warn() -> None:
        with toolkit.quiet_constrained_layout():
            _warn_from(message, UserWarning, "matplotlib.figure")

    # ACT
    raised = _raised_by(quiet_warn, UserWarning, "constrained_layout not applied")

    # ASSERT
    assert len(raised) == 1


def test_the_quiet_layout_restores_the_warning_filters_when_the_block_raises() -> None:
    # ARRANGE
    filtersBefore = list(warnings.filters)

    # ACT
    with pytest.raises(ValueError, match="draw failed"), toolkit.quiet_constrained_layout():
        raise ValueError("draw failed")

    # ASSERT
    assert warnings.filters == filtersBefore


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
# THE MATPLOTLIB "w-" FORMAT STRING
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
