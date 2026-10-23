"""A QC plot function closes its figure when plotting raises (DY-1710).

Each case opens the figure for real, then forces a failure after it is open. Most cases make the
save or show call raise. Where a function cannot reach that point with a small input, the case
feeds it a table that fails once the figure is open. The test then checks that no figure is left
registered in pyplot.
"""

from __future__ import annotations

import importlib
import sys
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.create_dispersion_map import _plot_slit_index_comparisons, create_dispersion_map
from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.detect_order_edges import detect_order_edges
from soxspipe.commonutils.figure_guard import closes_its_figures
from soxspipe.commonutils.horne_extraction import horne_extraction, plot_rectified_images
from tests.factories import qc_table
from tests.unit._plot_spies import spy_figures
from tests.unit._qc_plot_builders import (
    call_plot_merged_spectrum_qc,
    continuum_detector,
    continuum_plot_inputs,
    order_strip,
    quicklook_ccd,
    quicklook_order_pixels,
    quicklook_skylines,
    quicklook_spline,
    quicklook_subtractor,
    strip_subtractor,
    stub_merged_skylines,
    synthetic_merged_orders,
)

pytestmark = pytest.mark.unit

CaseRunner = Callable[[Path, object, pytest.MonkeyPatch], None]


class OutputFailureError(Exception):
    """Raised by the patched save and show calls, after the figure is open."""


def _fail(*args: object, **kwargs: object) -> None:
    raise OutputFailureError("forced failure after the figure opened")


@pytest.fixture(autouse=True)
def isolated_figures():
    """Start and end every test with no figure registered in pyplot."""
    plt.close("all")
    yield
    plt.close("all")


@pytest.fixture
def failing_output(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make `plt.savefig` and `plt.show` raise, and `plt.pause` a no-op."""
    monkeypatch.setattr(plt, "savefig", _fail)
    monkeypatch.setattr(plt, "show", _fail)
    monkeypatch.setattr(plt, "pause", lambda *args, **kwargs: None)
    # THE DEBUG PLOTS SWITCH TO AN INTERACTIVE BACKEND, WHICH A HEADLESS RUN CANNOT USE
    monkeypatch.setattr(matplotlib, "use", lambda *args, **kwargs: None)


# ---------------------------------------------------------------------------
# toolkit
# ---------------------------------------------------------------------------


def _quicklook_image(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    toolkit.quicklook_image(log, quicklook_ccd(), show=True)


def _quicklook_image_surface(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    frame = np.arange(64, dtype=float).reshape(8, 8) + 20.0
    toolkit.quicklook_image(log, frame, show=True, inst="SOXS", surfacePlot=True)


def _plot_merged_spectrum_qc(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    stub_merged_skylines(monkeypatch)
    call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, synthetic_merged_orders())


# ---------------------------------------------------------------------------
# subtract_sky
# ---------------------------------------------------------------------------


def _plot_sky_sampling(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    # THE SURFACE QUICKLOOK OF THE OBJECT FRAME OPENS AND CLOSES ITS OWN FIGURE
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **kwargs: None)
    subtractor = strip_subtractor(log, tmp_path, dispersionAxis="x", rotate=False)
    subtractor.plot_sky_sampling(order=11, imageMapOrderDF=order_strip(), knotLocations=np.array([502.0, 508.0]))


def _plot_image_comparison(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    shape = (6, 8)
    subtractor = strip_subtractor(log, tmp_path, dispersionAxis="x", rotate=False, shape=shape)
    subtractor.mapDF = pd.DataFrame({"x": [1, 2], "y": [1, 4]})
    frame = subtractor.objectFrame
    skyModelFrame = CCDData(frame.data * 2.0, unit="electron", mask=np.zeros(shape, dtype=bool))
    skySubFrame = CCDData(frame.data * -3.0, unit="electron", mask=np.zeros(shape, dtype=bool))
    subtractor.plot_image_comparison(frame, skyModelFrame, skySubFrame)


def _plot_order_skymodel_fitting_quicklook(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(toolkit, "get_skylines_dataframe", lambda log, settings, arm: quicklook_skylines())
    subtractor = quicklook_subtractor(log)
    subtractor.plot_order_skymodel_fitting_quicklook(
        quicklook_order_pixels(), quicklook_spline(), knots=np.array([503.0, 507.0])
    )


# ---------------------------------------------------------------------------
# horne_extraction
# ---------------------------------------------------------------------------


def _plot_extracted_spectrum_qc(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.traceFallback = False
    # CONCATENATING NO EXTRACTIONS FAILS ONCE THE FIGURE IS OPEN
    extractor.plot_extracted_spectrum_qc([])


def _plot_skyline_shift_diagnostic(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.arm = "VIS"
    shiftArray = np.linspace(500.0, 510.0, 20)
    extractor._plot_skyline_shift_diagnostic(
        shiftArray=shiftArray,
        skyValsOriginal=np.ones(20),
        skyVals=np.ones(20),
        peaks=np.array([3, 9]),
        objectVals=np.ones(20),
        localSkylines=[503.0],
        calibrationSkylines=[504.0],
        matchedSkylinePixels=[503.0],
        matchedShifts=[0.1],
        clippedWave=[],
        clippedShifts=[],
        medianShift=0.1,
        order=10,
        wmin=500.0,
        wmax=510.0,
        pixelScaleMedian=0.5,
        shiftLabel="wavelength (nm)",
        iteration=0,
        shiftUnits="nm",
    )


def _plot_rectified_images(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    plot_rectified_images(orderRectifiedImages={"flux": np.arange(12.0).reshape(3, 4)}, order=10)


# ---------------------------------------------------------------------------
# detect_continuum and detect_order_edges
# ---------------------------------------------------------------------------


def _continuum_plot_results(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(toolkit, "qc_settings_plot_tables", lambda **_: None)
    detector = continuum_detector(log, tmp_path, recipeName="soxs-order-centres")
    pixels, polynomials, clipped = continuum_plot_inputs()
    detector.plot_results(orderPixelTable=pixels, orderPolyTable=polynomials, clippedData=clipped)


def _order_edges_plot_results(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    detector = object.__new__(detect_order_edges)
    detector.log = log
    detector.axisA = "x"
    detector.axisB = "y"
    detector.axisAbin = 1
    detector.axisBbin = 1
    detector.arm = "VIS"
    detector.orderDeg = 1
    detector.axisBDeg = 1
    detector.detectorParams = {"rotate-qc-plot": 0, "flip-qc-plot": 0, "dispersion-axis": "x"}
    detector.flatFrame = CCDData(np.full((32, 32), 20.0), unit=u.electron)
    detector.qc = qc_table()
    detector.recipeSettings = {}
    edges = pd.DataFrame(
        {
            "order": [10.0, 10.0],
            "xcoord_lower": [4.0, 5.0],
            "xcoord_upper": [20.0, 21.0],
            "xcoord_lower_fit_res": [0.1, -0.1],
            "xcoord_upper_fit_res": [0.2, -0.2],
            "ycoord": [3.0, 9.0],
        }
    )
    # AN EMPTY META TABLE FAILS THE PER-ORDER LOOKUP, WHICH RUNS AFTER THE FIGURE IS OPEN
    detector.plot_results(
        orderPixelTableUpper=edges,
        orderPixelTableLower=edges,
        orderPolyTable=pd.DataFrame(),
        orderMetaTable=pd.DataFrame(),
        clippedDataUpper=edges,
        clippedDataLower=edges,
    )


# ---------------------------------------------------------------------------
# create_dispersion_map
# ---------------------------------------------------------------------------


def _dispersion_map_qc_plot(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.kw = lambda keyword: keyword
    mapper.detectorParams = {"rotate-qc-plot": 0, "flip-qc-plot": 0}
    mapper.pinholeFrameMasked = CCDData(np.full((32, 32), 20.0), unit=u.electron, mask=np.zeros((32, 32), dtype=bool))
    mapper.qc = qc_table()
    mapper.recipeSettings = {}
    mapper.exptime = 60.0
    mapper.arcFrame = False
    lines = pd.DataFrame({"residuals_x": [0.1, -0.2], "residuals_y": [0.0, 0.3]})
    # THE FRAME FLUX STATISTICS ARE ABSENT, SO THE FIRST IMAGE PANEL FAILS ONCE THE FIGURE IS OPEN
    mapper._create_dispersion_map_qc_plot(
        xcoeff=None,
        ycoeff=None,
        orderDeg=1,
        wavelengthDeg=1,
        slitDeg=False,
        orderPixelTable=lines,
        missingLines=pd.DataFrame(),
        allClippedLines=pd.DataFrame(),
    )


def _clip_on_measured_line_metrics(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.firstGuessMap = False
    # A TABLE WITHOUT THE LINE METRIC COLUMNS FAILS ONCE THE FIGURE IS OPEN
    mapper._clip_on_measured_line_metrics(pd.DataFrame({"wavelength": [500.0]}))


def _plot_predicted_line_windows(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    mapper = object.__new__(create_dispersion_map)
    mapper.meanFrameFlux = 10.0
    mapper.stdFrameFlux = 2.0
    mapper.windowSize = 5
    lines = pd.DataFrame({"detector_x": [8.0, 20.0], "detector_y": [8.0, 20.0]})
    mapper._plot_predicted_line_windows(lines, np.full((32, 32), 12.0))


def _convert_and_fit_plot(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.axisA = "x"
    mapper.axisB = "y"
    mapper.dispersionMapPath = "synthetic.fits"
    mapper.map_to_image_displacement_threshold = 0.3
    slitMap = SimpleNamespace(data=np.full((3, 3), np.nan, dtype=float))
    wavelengthMap = SimpleNamespace(data=np.full((3, 3), np.nan, dtype=float))
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")

    def project(*, orderPixelTable: pd.DataFrame, **_: object) -> pd.DataFrame:
        return orderPixelTable.assign(fit_x=[1.1, 3.6], fit_y=[1.1, 0.2])

    monkeypatch.setattr(module, "dispersion_map_to_pixel_arrays", project)
    mapper.convert_and_fit(
        order=10,
        bigWlArray=np.array([500.0, 501.0]),
        bigSlitArray=np.array([0.0, 1.0]),
        slitMap=slitMap,
        wlMap=wavelengthMap,
        iteration=1,
        plots=True,
    )


def _slit_index_comparisons(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    differences = pd.DataFrame(
        {
            "order": [10, 10, 10, 10],
            "slit_index": [0, 0, 1, 1],
            "xy_diff": [0.1, 0.2, 0.3, 0.4],
            "x_diff": [0.05, 0.1, 0.15, 0.2],
            "y_diff": [0.05, 0.1, 0.15, 0.2],
        }
    )
    _plot_slit_index_comparisons(differences)


CASES: list[tuple[str, CaseRunner, type[Exception]]] = [
    ("toolkit.quicklook_image", _quicklook_image, OutputFailureError),
    ("toolkit.quicklook_image-surface", _quicklook_image_surface, OutputFailureError),
    ("toolkit.plot_merged_spectrum_qc", _plot_merged_spectrum_qc, OutputFailureError),
    ("subtract_sky.plot_sky_sampling", _plot_sky_sampling, OutputFailureError),
    ("subtract_sky.plot_image_comparison", _plot_image_comparison, OutputFailureError),
    ("subtract_sky.plot_order_skymodel_fitting_quicklook", _plot_order_skymodel_fitting_quicklook, OutputFailureError),
    ("horne_extraction.plot_extracted_spectrum_qc", _plot_extracted_spectrum_qc, ValueError),
    ("horne_extraction._plot_skyline_shift_diagnostic", _plot_skyline_shift_diagnostic, OutputFailureError),
    ("horne_extraction.plot_rectified_images", _plot_rectified_images, OutputFailureError),
    ("detect_continuum.plot_results", _continuum_plot_results, OutputFailureError),
    ("detect_order_edges.plot_results", _order_edges_plot_results, KeyError),
    ("create_dispersion_map._create_dispersion_map_qc_plot", _dispersion_map_qc_plot, AttributeError),
    ("create_dispersion_map._clip_on_measured_line_metrics", _clip_on_measured_line_metrics, KeyError),
    ("create_dispersion_map._plot_predicted_line_windows", _plot_predicted_line_windows, OutputFailureError),
    ("create_dispersion_map.convert_and_fit", _convert_and_fit_plot, OutputFailureError),
    ("create_dispersion_map._plot_slit_index_comparisons", _slit_index_comparisons, OutputFailureError),
]


@pytest.mark.parametrize(("runCase", "expectedError"), [case[1:] for case in CASES], ids=[case[0] for case in CASES])
def test_a_qc_plot_closes_its_figure_when_plotting_raises(
    runCase: CaseRunner,
    expectedError: type[Exception],
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    failing_output: None,
) -> None:
    # ARRANGE
    openedFigures = spy_figures(monkeypatch)
    figureNumbersBefore = plt.get_fignums()

    # ACT
    with pytest.raises(expectedError):
        runCase(tmp_path, log, monkeypatch)

    # ASSERT
    assert openedFigures, "the case must fail after the plot opened its figure"
    assert plt.get_fignums() == figureNumbersBefore


def test_a_qc_plot_leaves_unrelated_open_figures_alone_when_plotting_raises(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, failing_output: None
) -> None:
    # ARRANGE
    unrelatedFigure = plt.figure()
    figureNumbersBefore = plt.get_fignums()

    # ACT
    with pytest.raises(OutputFailureError):
        _plot_image_comparison(tmp_path, log, monkeypatch)

    # ASSERT
    assert plt.get_fignums() == figureNumbersBefore == [unrelatedFigure.number]


def test_a_qc_plot_that_closes_its_figure_early_closes_it_once(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, failing_output: None
) -> None:
    # ARRANGE: WITHOUT A SPLINE THE QUICKLOOK PAUSES, CLOSES EVERYTHING AND RETURNS EARLY
    monkeypatch.setattr(toolkit, "get_skylines_dataframe", lambda log, settings, arm: quicklook_skylines())
    closeCalls: list[object] = []
    realClose = plt.close

    def recording_close(*args: object, **kwargs: object) -> None:
        closeCalls.append(args)
        realClose(*args, **kwargs)

    monkeypatch.setattr(plt, "close", recording_close)
    figureNumbersBefore = plt.get_fignums()

    # ACT
    quicklook_subtractor(log).plot_order_skymodel_fitting_quicklook(quicklook_order_pixels(), None, title="clipping")

    # ASSERT
    assert closeCalls == [("all",)]
    assert plt.get_fignums() == figureNumbersBefore


def test_a_qc_plot_still_writes_its_pdf_when_plotting_succeeds(tmp_path: Path, log: object) -> None:
    # ARRANGE
    subtractor = strip_subtractor(log, tmp_path, dispersionAxis="x", rotate=False, shape=(6, 8))
    subtractor.mapDF = pd.DataFrame({"x": [1, 2], "y": [1, 4]})
    frame = subtractor.objectFrame
    skyModelFrame = CCDData(frame.data * 2.0, unit="electron", mask=np.zeros((6, 8), dtype=bool))
    skySubFrame = CCDData(frame.data * -3.0, unit="electron", mask=np.zeros((6, 8), dtype=bool))

    # ACT
    filePath = subtractor.plot_image_comparison(frame, skyModelFrame, skySubFrame)

    # ASSERT
    assert Path(filePath).read_bytes().startswith(b"%PDF")
    assert plt.get_fignums() == []


def test_the_guard_returns_the_result_and_closes_a_figure_the_function_left_open() -> None:
    # ARRANGE
    @closes_its_figures
    def plot_and_return_a_value() -> str:
        plt.figure()
        return "result"

    # ACT
    result = plot_and_return_a_value()

    # ASSERT
    assert result == "result"
    assert plt.get_fignums() == []


def test_the_guard_does_not_close_again_a_figure_the_function_closed_itself(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    closeCalls: list[object] = []
    realClose = plt.close

    def recording_close(figure: object) -> None:
        closeCalls.append(figure)
        realClose(figure)

    monkeypatch.setattr(plt, "close", recording_close)

    @closes_its_figures
    def plot_and_close() -> None:
        plt.close(plt.figure())

    # ACT
    plot_and_close()

    # ASSERT
    assert len(closeCalls) == 1
    assert plt.get_fignums() == []


def test_the_guard_closes_a_new_figure_that_reuses_the_number_of_a_figure_closed_inside_the_call() -> None:
    # ARRANGE: PYPLOT HANDS OUT THE NUMBER OF A CLOSED FIGURE AGAIN
    unrelatedFigure = plt.figure()

    @closes_its_figures
    def close_everything_then_plot_and_raise() -> None:
        plt.close("all")
        newFigure = plt.figure()
        assert newFigure.number == unrelatedFigure.number
        raise OutputFailureError

    # ACT
    with pytest.raises(OutputFailureError):
        close_everything_then_plot_and_raise()

    # ASSERT
    assert plt.get_fignums() == []


def test_a_guarded_function_inside_a_guarded_function_closes_the_figure_after_a_real_close_all() -> None:
    # ARRANGE
    unrelatedFigure = plt.figure()

    @closes_its_figures
    def inner() -> None:
        plt.close("all")
        plt.figure()
        raise OutputFailureError

    @closes_its_figures
    def outer() -> None:
        plt.figure()
        inner()

    # ACT
    with pytest.raises(OutputFailureError):
        outer()

    # ASSERT
    assert unrelatedFigure.number == 1
    assert plt.get_fignums() == []


def test_a_failure_while_closing_a_figure_does_not_mask_the_original_error(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    # ARRANGE
    @closes_its_figures
    def plot_and_raise() -> None:
        plt.figure()
        raise OutputFailureError("the plot failed")

    def failing_close(*args: object, **kwargs: object) -> None:
        raise RuntimeError("the close failed")

    # ACT: THE PATCH ENDS WITH THE BLOCK, BEFORE THE FIXTURE THAT CLOSES EVERY FIGURE RUNS
    with monkeypatch.context() as patch, caplog.at_level("WARNING"):
        patch.setattr(plt, "close", failing_close)
        with pytest.raises(OutputFailureError, match="the plot failed"):
            plot_and_raise()

    # ASSERT
    assert "could not close figure" in caplog.text


def test_a_failure_while_closing_one_figure_does_not_stop_the_others_closing(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    realClose = plt.close
    attempts: list[object] = []

    def close_first_fails(figure: object) -> None:
        attempts.append(figure)
        if len(attempts) == 1:
            raise RuntimeError("the close failed")
        realClose(figure)

    @closes_its_figures
    def plot_two_figures() -> None:
        plt.figure()
        plt.figure()

    # ACT: THE PATCH ENDS WITH THE BLOCK, BEFORE THE FIXTURE THAT CLOSES EVERY FIGURE RUNS
    with monkeypatch.context() as patch:
        patch.setattr(plt, "close", close_first_fails)
        plot_two_figures()

    # ASSERT
    assert len(attempts) == 2
    assert len(plt.get_fignums()) == 1


def test_the_guard_does_not_import_pyplot_when_the_function_never_plots(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE: HIDE PYPLOT, AS IT IS IN A RUN THAT HAS NOT PLOTTED YET
    monkeypatch.delitem(sys.modules, "matplotlib.pyplot")

    # ACT
    result = closes_its_figures(lambda: "no plot")()

    # ASSERT
    assert result == "no plot"
    assert "matplotlib.pyplot" not in sys.modules


def test_sample_trace_closes_a_figure_that_it_opens_when_it_raises(log: object) -> None:
    # ARRANGE: THE DEBUG PLOTS ARE DEEP IN A LONG BODY, SO AN EARLY CALL OPENS THE FIGURE AND FAILS
    detector = object.__new__(detect_continuum)
    detector.log = log

    def open_a_figure_and_fail() -> None:
        plt.figure()
        raise OutputFailureError

    detector.create_pixel_arrays = open_a_figure_and_fail

    # ACT
    with pytest.raises(OutputFailureError):
        detector.sample_trace()

    # ASSERT
    assert plt.get_fignums() == []
