"""Characterization tests freezing the current rendering behaviour of
`quicklook_image` and `plot_merged_spectrum_qc` in
`soxspipe/commonutils/toolkit.py`, ahead of the DY-79 part-2 split of those
two functions into smaller helpers.

These tests do not judge whether the current behaviour is good design -- they
pin it, including behaviour that looks like a defect (noted inline as
BUG-LIKE). Any change in the split that alters a pinned artist property,
axis limit, label, colour, or call argument must fail one of these tests.

Unlike `tests/unit/test_toolkit_plotting.py` and
`tests/unit/test_toolkit_characterization.py` (which pin the PDF output, the
`savefig`/`clf`/`close` call sequence, and the returned product-table rows),
this module pins the actual matplotlib *artists* produced -- figure sizes,
axes, line data, colours, text, and axis limits -- by spying on
`matplotlib.pyplot.figure` and inspecting the real (closed-but-inspectable)
`Figure` objects it returns.

IMPORTANT: `quicklook_image` calls `plt.clf()` after `save_qc_plot()` when
`saveToPath` is given, which strips every axis back off the figure. Any test
that needs to inspect the rendered artists must therefore call the function
with `show=True` (with `plt.show` monkeypatched to a no-op) and *without*
`saveToPath`, rather than the `saveToPath` pattern used by the PDF-existence
tests in the sibling modules.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.stats import sigma_clip, sigma_clipped_stats
from matplotlib.figure import Figure

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import toolkit
from tests.factories import product_table

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------


def _spy_figures(monkeypatch: pytest.MonkeyPatch) -> list[Figure]:
    """Patch `matplotlib.pyplot.figure`, delegating to the real function while
    recording every `Figure` it returns.

    **Key Arguments:**

    - ``monkeypatch`` -- the pytest monkeypatch fixture

    **Return:**

    - a list that is appended to (in call order) with each `Figure` produced
    """
    original = plt.figure
    figures: list[Figure] = []

    def _wrapper(*args: object, **kwargs: object) -> Figure:
        figure = original(*args, **kwargs)
        figures.append(figure)
        return figure

    monkeypatch.setattr(plt, "figure", _wrapper)
    return figures


def _quiet_show(monkeypatch: pytest.MonkeyPatch) -> list[tuple]:
    """Patch `matplotlib.pyplot.show` to a no-op recorder so `show=True` code
    paths can be exercised without a display, without clearing the figure."""
    calls: list[tuple] = []
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: calls.append((args, kwargs)))
    return calls


def _image_axis(figure: Figure):
    """Return the first axis on `figure` that holds an `AxesImage`."""
    return next(axis for axis in figure.axes if axis.images)


def _ccd(instrument: str = "SOXS", shape: tuple[int, int] = (8, 8), extraHeader: dict | None = None) -> CCDData:
    size = shape[0] * shape[1]
    data = np.arange(size, dtype=float).reshape(shape)
    mask = np.zeros(shape, dtype=bool)
    header = {"INSTRUME": instrument}
    if extraHeader:
        header.update(extraHeader)
    return CCDData(
        data,
        unit=u.electron,
        mask=mask,
        uncertainty=StdDevUncertainty(np.ones(shape)),
        meta=header,
    )


class _StubKeywordLookup:
    """Identity keyword lookup: `.get(key)` returns `key` unchanged."""

    def __init__(self, **kwargs: object) -> None:
        pass

    def get(self, key: str) -> str:
        return key


class _StubDetectorLookup:
    """Detector lookup stub returning a fixed `science-pixels` entry."""

    def __init__(self, **kwargs: object) -> None:
        pass

    def get(self, arm: str) -> dict[str, tuple[int, int, int, int]]:
        return {"science-pixels": (0, 8, 0, 8)}


def _stub_header_lookups(monkeypatch: pytest.MonkeyPatch) -> None:
    """Stub the *function-local* `from soxspipe.commonutils import
    detector_lookup, keyword_lookup` re-import inside `quicklook_image`.

    Both names are imported locally inside the function body, shadowing the
    module-level `toolkit.keyword_lookup`/`toolkit.detector_lookup` bindings
    used elsewhere in the module -- so the source attributes on
    `soxspipe.commonutils` must be patched, not the `toolkit` module ones.
    """
    monkeypatch.setattr(commonutils, "keyword_lookup", _StubKeywordLookup)
    monkeypatch.setattr(commonutils, "detector_lookup", _StubDetectorLookup)


def _stub_grid_lines(monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    """Stub `toolkit.create_dispersion_solution_grid_lines_for_plot`, recording
    the kwargs it is called with and returning a synthetic 3-line grid table
    plus a benign `interOrderMask`.

    **Return:**

    - a dict that is populated with the kwargs the stub was called with
    """
    captured: dict[str, object] = {}
    gridLinePixelTable = pd.DataFrame(
        {
            "line": [0, 0, 1, 1, 2, 2],
            "fit_x": [10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
            "fit_y": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        }
    )
    interOrderMask = np.zeros((8, 8), dtype=bool)
    interOrderMask[0, 0] = True

    def _stub(**kwargs: object) -> tuple[pd.DataFrame, np.ndarray]:
        captured.update(kwargs)
        return gridLinePixelTable, interOrderMask

    monkeypatch.setattr(toolkit, "create_dispersion_solution_grid_lines_for_plot", _stub)
    captured["interOrderMask"] = interOrderMask
    captured["gridLinePixelTable"] = gridLinePixelTable
    return captured


# ---------------------------------------------------------------------------
# quicklook_image
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("instrument", "expectedRotation"),
    [
        ("SOXS", lambda frame: frame),
        ("XSHOOTER", lambda frame: np.flipud(np.rot90(frame, 1))),
        ("OTHER", lambda frame: np.flipud(frame)),
    ],
)
def test_quicklook_image_orients_pixels_and_sets_clim_per_instrument(
    monkeypatch: pytest.MonkeyPatch, log: object, instrument: str, expectedRotation
) -> None:
    """The displayed image data and colour limits are pinned per instrument
    rotation convention and the sigma-clipped-stats-derived `vmin`/`vmax`."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    ccd = _ccd(instrument)

    toolkit.quicklook_image(log, ccd, show=True)

    axis = _image_axis(figures[-1])
    image = axis.images[0]
    np.testing.assert_array_equal(image.get_array().data, expectedRotation(ccd.data))

    mean, median, std = sigma_clipped_stats(
        ccd.data, sigma=50.0, stdfunc="mad_std", cenfunc="median", maxiters=3
    )
    vmin, vmax = image.get_clim()
    assert vmin == pytest.approx(median - 3 * 0.5 * std, rel=1e-12)
    assert vmax == pytest.approx(median + 3 * 0.5 * std, rel=1e-12)


def test_quicklook_image_clim_respects_custom_std_window(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """A non-default `stdWindow` scales the colour-limit half-width."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    ccd = _ccd("SOXS")

    toolkit.quicklook_image(log, ccd, show=True, stdWindow=7)

    image = _image_axis(figures[-1]).images[0]
    mean, median, std = sigma_clipped_stats(
        ccd.data, sigma=50.0, stdfunc="mad_std", cenfunc="median", maxiters=3
    )
    vmin, vmax = image.get_clim()
    assert vmin == pytest.approx(median - 7 * 0.5 * std, rel=1e-12)
    assert vmax == pytest.approx(median + 7 * 0.5 * std, rel=1e-12)


def test_quicklook_image_figure_size_and_box_aspect_flip_for_tall_images(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """A `rotatedImg` whose height exceeds its width by more than 1000px flips
    the figure to portrait and the image-axis box aspect to 2.0."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    normalFrame = np.zeros((8, 8))
    tallFrame = np.zeros((1002, 1))

    toolkit.quicklook_image(log, normalFrame, show=True, inst="OTHER")
    normalFigure = figures[-1]
    toolkit.quicklook_image(log, tallFrame, show=True, inst="OTHER")
    tallFigure = figures[-1]

    assert tuple(normalFigure.get_size_inches()) == (12.0, 5.0)
    assert _image_axis(normalFigure).get_box_aspect() == pytest.approx(0.5)
    assert tuple(tallFigure.get_size_inches()) == (5.0, 12.0)
    assert _image_axis(tallFigure).get_box_aspect() == pytest.approx(2.0)


@pytest.mark.parametrize(
    ("frameValues", "expectedFormatterType", "expectedFmt"),
    [
        (np.linspace(0.0, 5.0, 64), "ScalarFormatter", None),
        (np.linspace(20.0, 40.0, 64), "FormatStrFormatter", "%1.0f"),
    ],
)
def test_quicklook_image_colorbar_formatter_switches_on_mean_threshold(
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    frameValues: np.ndarray,
    expectedFormatterType: str,
    expectedFmt: str | None,
) -> None:
    """The colorbar uses a fixed `"%1.0f"` formatter when the sigma-clipped
    mean exceeds 10, and matplotlib's default `ScalarFormatter` otherwise."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    frame = frameValues.reshape(8, 8)

    toolkit.quicklook_image(log, frame, show=True, inst="OTHER")

    image = _image_axis(figures[-1]).images[0]
    formatter = image.colorbar.formatter
    assert type(formatter).__name__ == expectedFormatterType
    if expectedFmt is not None:
        assert formatter.fmt == expectedFmt


def test_quicklook_image_title_sets_suptitle_with_fixed_fontsize(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)

    toolkit.quicklook_image(log, np.zeros((8, 8)), show=True, inst="OTHER", title="My Title")

    suptitle = figures[-1]._suptitle
    assert suptitle.get_text() == "My Title"
    assert suptitle.get_fontsize() == 20


def test_quicklook_image_without_title_has_no_suptitle(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)

    toolkit.quicklook_image(log, np.zeros((8, 8)), show=True, inst="OTHER")

    assert figures[-1]._suptitle is None


@pytest.mark.parametrize(
    ("instrument", "expectedXLabel", "expectedYLabel", "expectedYlim"),
    [
        # `imshow` DEFAULTS TO `origin="upper"`, SO THE UN-INVERTED YLIM IS
        # ALREADY DESCENDING (7.5, -0.5); ONLY XSHOOTER'S EXPLICIT
        # `invert_yaxis()` CALL FLIPS IT BACK TO ASCENDING
        ("SOXS", "x-axis", "y-axis", (7.5, -0.5)),
        ("XSHOOTER", "y-axis", "x-axis", (-0.5, 7.5)),
        ("OTHER", "y-axis", "x-axis", (7.5, -0.5)),
    ],
)
def test_quicklook_image_axis_labels_and_yaxis_inversion_per_instrument(
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    instrument: str,
    expectedXLabel: str,
    expectedYLabel: str,
    expectedYlim: tuple[float, float],
) -> None:
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)

    toolkit.quicklook_image(log, np.zeros((8, 8)), show=True, inst=instrument)

    axis = _image_axis(figures[-1])
    assert axis.get_xlabel() == expectedXLabel
    assert axis.get_ylabel() == expectedYLabel
    assert axis.get_ylim() == pytest.approx(expectedYlim)


def test_quicklook_image_ext_none_plain_array_falls_back_to_xshooter_orientation(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """A plain ndarray with no header, and an `ext` that is not one of
    `"data"`/`"mask"`/`"uncertainty"`, is used as-is (the `else: frame =
    CCDObject` branch) and defaults to the `"XSHOOTER"` rotation convention
    when `inst` is not supplied."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    frame = np.arange(64, dtype=float).reshape(8, 8)

    toolkit.quicklook_image(log, frame, show=True, ext=None)

    image = _image_axis(figures[-1]).images[0]
    np.testing.assert_array_equal(image.get_array().data, np.flipud(np.rot90(frame, 1)))


def test_quicklook_image_inst_argument_overrides_header_instrume(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """An explicit `inst=` argument takes precedence over `header['INSTRUME']`."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    ccd = _ccd("SOXS")

    toolkit.quicklook_image(log, ccd, show=True, inst="XSHOOTER")

    image = _image_axis(figures[-1]).images[0]
    np.testing.assert_array_equal(image.get_array().data, np.flipud(np.rot90(ccd.data, 1)))


@pytest.mark.parametrize("instrument", ["SOXS", "OTHER"])
def test_quicklook_image_surface_plot_axes_and_orientation(
    monkeypatch: pytest.MonkeyPatch, log: object, instrument: str
) -> None:
    """`surfacePlot=True` adds a 3D axis (plus the 2D image axis and its
    colorbar axis) with azimuth/elevation/z-limits and axis labels pinned per
    instrument."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    frame = np.arange(64, dtype=float).reshape(8, 8) + 20.0

    toolkit.quicklook_image(log, frame, show=True, inst=instrument, surfacePlot=True)

    figure = figures[-1]
    assert tuple(figure.get_size_inches()) == (20.0, 8.0)
    # THREE AXES: THE 3D SURFACE, THE 2D IMAGE, AND THE COLORBAR
    assert len(figure.axes) == 3
    axis3d = figure.axes[0]
    assert axis3d.azim == (70 if instrument == "SOXS" else -120)
    assert axis3d.elev == 30

    mean, median, std = sigma_clipped_stats(
        frame, sigma=50.0, stdfunc="mad_std", cenfunc="median", maxiters=3
    )
    vmax = median + 3 * 0.5 * std
    vmin = median - 3 * 0.5 * std
    zlo, zhi = axis3d.get_zlim()
    assert zlo == pytest.approx(vmin, rel=1e-12)
    assert zhi == pytest.approx(min(np.nanmax(frame), vmax * 1.2), rel=1e-12)

    if instrument == "SOXS":
        assert axis3d.get_xlabel() == "x-axis"
        assert axis3d.get_ylabel() == "y-axis"
    else:
        assert axis3d.get_xlabel() == "y-axis"
        assert axis3d.get_ylabel() == "x-axis"

    imageAxis = _image_axis(figure)
    assert imageAxis.get_box_aspect() == pytest.approx(0.5)


def test_quicklook_image_surface_plot_restores_rcparams_after_call(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`surfacePlot=True` forces `axes.edgecolor` (etc.) to a fixed palette
    during rendering, but the original rcParams are restored afterwards."""
    _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    monkeypatch.setitem(matplotlib.rcParams, "axes.edgecolor", "orange")

    toolkit.quicklook_image(
        log, np.arange(64, dtype=float).reshape(8, 8) + 20.0, show=True, inst="SOXS", surfacePlot=True
    )

    assert matplotlib.rcParams["axes.edgecolor"] == "orange"


@pytest.mark.parametrize(
    ("instrument", "expectedXY"),
    [
        ("SOXS", [([10.0, 20.0], [1.0, 2.0]), ([30.0, 40.0], [3.0, 4.0])]),
        ("XSHOOTER", [([1.0, 2.0], [10.0, 20.0]), ([3.0, 4.0], [30.0, 40.0])]),
    ],
)
def test_quicklook_image_dispmap_branch_forwards_arguments_and_draws_grid_lines(
    monkeypatch: pytest.MonkeyPatch, log: object, instrument: str, expectedXY: list
) -> None:
    """The `dispMapImage` branch forwards `dispMap`/`dispMapImage`/
    `associatedFrame`/`kw`/`skylines` to the grid-lines helper, plots
    `fit_x` vs `fit_y` for SOXS (swapped otherwise), and -- because the
    helper's synthetic table tops out at `line == 2` -- `range(int(max))`
    drops the last line: only lines 0 and 1 are drawn."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    _stub_header_lookups(monkeypatch)
    captured = _stub_grid_lines(monkeypatch)
    ccd = _ccd(instrument, extraHeader={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-01"})
    settings = {"instrument": "soxs"}

    toolkit.quicklook_image(
        log, ccd, show=True, dispMap="disp.fits", dispMapImage="disp_image.fits", settings=settings
    )

    assert captured["dispMap"] == "disp.fits"
    assert captured["dispMapImage"] == "disp_image.fits"
    assert captured["associatedFrame"] is ccd
    assert captured["kw"]("SEQ_ARM") == "SEQ_ARM"
    assert captured["skylines"] is False

    axis = _image_axis(figures[-1])
    assert len(axis.lines) == 2
    for line, (expectedX, expectedY) in zip(axis.lines, expectedXY, strict=True):
        np.testing.assert_array_equal(line.get_xdata(), expectedX)
        np.testing.assert_array_equal(line.get_ydata(), expectedY)
        assert line.get_color() == "black"
        assert line.get_linewidth() == pytest.approx(0.5)
        assert line.get_alpha() == pytest.approx(0.8)


def test_quicklook_image_dispmap_branch_never_actually_masks_the_ndarray_frame(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """BUG-LIKE: the inter-order mask is combined into `mask = (frame.mask ==
    1) | (interOrderMask == 1)` and then assigned back with `frame.mask =
    mask`. Because `frame` is `CCDObject.data` -- a plain `numpy.ndarray`,
    which has no settable `.mask` attribute -- both the read and the write
    raise `AttributeError`, which is swallowed. The `interOrderMask` computed
    by the grid-lines helper is therefore never applied to the displayed
    image, no matter what it contains. Pinned as-is, not fixed."""
    figures = _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    _stub_header_lookups(monkeypatch)
    captured = _stub_grid_lines(monkeypatch)
    assert captured["interOrderMask"].any()  # THE STUB MASK HAS SOME TRUE PIXELS
    ccd = _ccd("SOXS", extraHeader={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-01"})

    toolkit.quicklook_image(
        log, ccd, show=True, dispMap="disp.fits", dispMapImage="disp_image.fits", settings={"instrument": "soxs"}
    )

    image = _image_axis(figures[-1]).images[0]
    array = image.get_array()
    assert isinstance(array, np.ma.MaskedArray)
    assert not np.ma.is_masked(array)


def test_quicklook_image_skylines_true_forwards_dataframe_to_grid_lines_helper(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`skylines=True` fetches the skylines dataframe via
    `get_skylines_dataframe(log, settings, arm)` and forwards the exact same
    object as the `skylines=` argument to the grid-lines helper."""
    _spy_figures(monkeypatch)
    _quiet_show(monkeypatch)
    _stub_header_lookups(monkeypatch)
    captured = _stub_grid_lines(monkeypatch)
    skylinesDF = pd.DataFrame({"WAVELENGTH": [500.0], "ISOLATED": [True]})
    skylineCalls: list[tuple] = []

    def _stub_skylines(log: object, settings: object, arm: str) -> pd.DataFrame:
        skylineCalls.append((settings, arm))
        return skylinesDF

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", _stub_skylines)
    ccd = _ccd("SOXS", extraHeader={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-01"})
    settings = {"instrument": "soxs"}

    toolkit.quicklook_image(
        log,
        ccd,
        show=True,
        dispMap="disp.fits",
        dispMapImage="disp_image.fits",
        settings=settings,
        skylines=True,
    )

    assert skylineCalls == [(settings, "VIS")]
    assert captured["skylines"] is skylinesDF


def test_quicklook_image_show_true_calls_plt_show(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    _spy_figures(monkeypatch)
    showCalls = _quiet_show(monkeypatch)

    toolkit.quicklook_image(log, np.zeros((8, 8)), show=True, inst="OTHER")

    assert len(showCalls) == 1


# ---------------------------------------------------------------------------
# plot_merged_spectrum_qc
# ---------------------------------------------------------------------------


def _synthetic_merged_orders(*, includeSky: bool = True) -> pd.DataFrame:
    """~50-row synthetic merged spectrum with a couple of outliers in FLUX and
    SNR, deterministic (no RNG) so sigma-clipped statistics can be
    recomputed independently in each test."""
    wave = np.linspace(500.0, 600.0, 50)
    flux = 10.0 + 2.0 * np.sin(np.linspace(0.0, 10.0, 50))
    flux[5] = 500.0
    flux[40] = -500.0
    snr = 5.0 + 0.1 * np.arange(50)
    snr[10] = 200.0
    columns = {"WAVE": wave, "FLUX_COUNTS": flux, "SNR": snr}
    if includeSky:
        columns["SKY_COUNTS"] = 10.0 + np.linspace(0.0, 5.0, 50)
    return pd.DataFrame(columns)


def _stub_skylines(monkeypatch: pytest.MonkeyPatch, extra: tuple = ()) -> None:
    wavelengths = [520.0, 540.0, *extra]
    isolated = [True, False] + [False] * len(extra)
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame({"WAVELENGTH": wavelengths, "ISOLATED": isolated}),
    )


def _expected_flux_ylims(flux: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    """Independently recompute the top-panel and middle-panel ylims from the
    same sigma-clip/sigma-clipped-stats calls the production code uses."""
    arrayMask = sigma_clip(flux, sigma_lower=3, sigma_upper=15.0, maxiters=1, cenfunc="mean", stdfunc="std")
    _, _, std = sigma_clipped_stats(flux, sigma=5.0, stdfunc="std", cenfunc="mean", maxiters=3)
    topYlim = (arrayMask.min() - 3 * std, arrayMask.max() + 3 * std)
    middleYlim = (max(arrayMask.min() * 0.5, 0), arrayMask.max() * 2)
    return topYlim, middleYlim


def _expected_snr_ylim(snr: np.ndarray) -> tuple[float, float]:
    mean, _, std = sigma_clipped_stats(snr, sigma=5.0, stdfunc="std", cenfunc="mean", maxiters=3)
    return (0.0, mean + 4 * std)


def _call_plot_merged_spectrum_qc(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: object,
    merged: pd.DataFrame,
    **overrides: object,
) -> tuple[pd.DataFrame, str]:
    kwargs: dict[str, object] = {
        "products": product_table().iloc[0:0].copy(),
        "log": log,
        "qcDir": str(tmp_path),
        "filenameTemplate": "synthetic.fits",
        "noddingSequence": "_AB",
        "dateObs": "2024-01-02T03:04:05",
        "arm": "VIS",
        "recipeName": "soxs-nod",
        "orderJoins": {"11-12": 550.0},
        "fluxCalibrated": False,
        "qcTable": False,
        "settings": {},
    }
    kwargs.update(overrides)
    return toolkit.plot_merged_spectrum_qc(merged, **kwargs)


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_figure_layout_size_dpi_and_axis_order(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    """The figure is 14x10in at 180dpi with exactly four axes, created in
    top/middle/bottom(SNR)/sky order."""
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    figure = figures[-1]
    assert tuple(figure.get_size_inches()) == (14.0, 10.0)
    assert figure.dpi == 180
    assert len(figure.axes) == 4
    top, middle, bottom, sky = figure.axes
    assert middle.get_yscale() == "log"
    assert sky.get_yscale() == "log"
    assert bottom.get_ylabel() == "SNR"
    assert sky.get_ylabel() == "sky flux ($e^{-}$)"


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_xlim_matches_wave_range_on_every_panel(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    for axis in figures[-1].axes:
        assert axis.get_xlim() == pytest.approx((500.0, 600.0), rel=1e-12)


def test_plot_merged_spectrum_qc_wave_as_quantity_column_raises_type_error(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """NOT PRACTICAL AS A `.value`-BRANCH CHARACTERIZATION TEST: a pandas
    column of per-row `astropy.units.Quantity` scalars (the shape the
    `try: ...min().value` fallback appears to defend against) cannot reach
    that fallback at all -- `ax.plot()` itself raises `TypeError` while
    building the line path, before any `set_xlim()` call runs. Pinned as a
    crash, not a `.value` codepath."""
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()
    merged["WAVE"] = pd.Series([value * u.nm for value in merged["WAVE"]])

    with pytest.raises(TypeError):
        _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged)


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_ylim_and_yscale_pin_sigma_clipped_statistics(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, fluxCalibrated: bool
) -> None:
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    top, middle, bottom, sky = figures[-1].axes
    expectedTopYlim, expectedMiddleYlim = _expected_flux_ylims(merged["FLUX_COUNTS"].to_numpy())
    assert top.get_ylim() == pytest.approx(expectedTopYlim, rel=1e-12)
    assert middle.get_ylim() == pytest.approx(expectedMiddleYlim, rel=1e-12)
    assert middle.get_yscale() == "log"

    expectedBottomYlim = _expected_snr_ylim(merged["SNR"].to_numpy())
    assert bottom.get_ylim() == pytest.approx(expectedBottomYlim, rel=1e-12)

    expectedSkyYlim = (10.0, merged["SKY_COUNTS"].max() * 1.1)
    assert sky.get_ylim() == pytest.approx(expectedSkyYlim, rel=1e-12)
    assert sky.get_yscale() == "log"


def test_plot_merged_spectrum_qc_sky_panel_has_no_line_and_pinned_default_ylim_when_sky_counts_absent(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """When `SKY_COUNTS` is absent, the sky flux line is never drawn and the
    axis keeps matplotlib's own autoscaled default for an empty log-scale
    axis (recorded here as a literal, since it is not derived from the
    input data at all)."""
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders(includeSky=False)

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

    sky = figures[-1].axes[3]
    # ONLY THE TWO SKYLINE AXVLINES ARE PRESENT; NO SKY-FLUX LINE IS DRAWN
    assert len(sky.lines) == 2
    assert sky.get_ylim() == pytest.approx((0.8912509381337456, 11.220184543019636), rel=1e-9)


@pytest.mark.parametrize(
    ("fluxCalibrated", "expectedLabelSuffix", "expectedTopColour", "expectedSkyColour"),
    [
        (False, "flux ($e^{-}$)", "#dc322f", "#859900"),
        (True, "flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)", "#2aa198", "#2aa198"),
    ],
)
def test_plot_merged_spectrum_qc_labels_title_and_line_colours(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    fluxCalibrated: bool,
    expectedLabelSuffix: str,
    expectedTopColour: str,
    expectedSkyColour: str,
) -> None:
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, fluxCalibrated=fluxCalibrated)

    top, middle, bottom, sky = figures[-1].axes
    assert top.get_ylabel() == expectedLabelSuffix
    assert top.get_title() == "Optimally Extracted Order-Merged Object Spectrum (VIS)\nsynthetic"
    assert bottom.get_xlabel() == "wavelength (nm)"

    topFluxLine = top.lines[0]
    assert topFluxLine.get_color() == expectedTopColour
    skyFluxLine = sky.lines[0]
    assert skyFluxLine.get_color() == expectedSkyColour
    snrLine = bottom.lines[0]
    assert snrLine.get_color() == "black"


def test_plot_merged_spectrum_qc_snr_text_sorts_global_vis_letters_then_numeric_orders(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The SNR text box sorts `GLOBAL` first, then VIS arm letters in
    `_vis_rank` order, then numeric orders ascending, then any order name
    that is neither a VIS letter nor numeric (the `_order_key` fallback
    `except` branch); duplicate `(qc_name, qc_order)` rows keep the *last*
    occurrence's value."""
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()
    qcTable = pd.DataFrame(
        {
            "qc_name": ["SNR MEDIAN"] * 8,
            "qc_order": [np.nan, "12", "11", "12", "u", "g", "r", "extra"],
            # THE SECOND "12" ROW (VALUE 9.0) MUST WIN OVER THE FIRST (5.0)
            "qc_value": [6.0, 5.0, 7.0, 9.0, 1.0, 2.0, 3.0, 4.0],
        }
    )

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=qcTable, orderJoins=False)

    bottom = figures[-1].axes[2]
    snrTexts = [text for text in bottom.texts if ":" in text.get_text()]
    assert len(snrTexts) == 1
    assert snrTexts[0].get_text() == "GLOBAL: 6\nu: 1\ng: 2\nr: 3\n11: 7\n12: 9\nextra: 4"


def test_plot_merged_spectrum_qc_snr_text_absent_when_qc_table_is_false(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=False, orderJoins=False)

    bottom = figures[-1].axes[2]
    assert len(bottom.texts) == 0


def test_plot_merged_spectrum_qc_snr_text_absent_when_no_snr_median_rows(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()
    qcTable = pd.DataFrame({"qc_name": ["OTHER"], "qc_order": ["GLOBAL"], "qc_value": [1.0]})

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, qcTable=qcTable, orderJoins=False)

    bottom = figures[-1].axes[2]
    assert len(bottom.texts) == 0


def test_plot_merged_spectrum_qc_order_joins_draw_axvline_and_text_on_three_panels_not_sky(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Each order join draws a dashed black `axvline` and an "ORDER JOIN"
    text on the top/middle/bottom panels -- but *not* the sky panel -- with
    the text placed at `(v + 5, 0.9 * ylim[1])` using each panel's own ylim
    at that point."""
    figures = _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins={"11-12": 550.0})

    top, middle, bottom, sky = figures[-1].axes
    expectedTopYlim, expectedMiddleYlim = _expected_flux_ylims(merged["FLUX_COUNTS"].to_numpy())
    expectedBottomYlim = _expected_snr_ylim(merged["SNR"].to_numpy())

    for panel, expectedYlimTop in (
        (top, expectedTopYlim[1]),
        (middle, expectedMiddleYlim[1]),
        (bottom, expectedBottomYlim[1]),
    ):
        joinLines = [line for line in panel.lines if line.get_linestyle() == "--"]
        assert len(joinLines) == 1
        assert joinLines[0].get_color() == "black"
        assert joinLines[0].get_xdata()[0] == pytest.approx(550.0)

        joinTexts = [text for text in panel.texts if text.get_text() == "ORDER JOIN"]
        assert len(joinTexts) == 1
        x, y = joinTexts[0].get_position()
        assert x == pytest.approx(555.0)
        assert y == pytest.approx(0.9 * expectedYlimTop, rel=1e-12)

    assert not any(line.get_linestyle() == "--" for line in sky.lines)
    assert not any(text.get_text() == "ORDER JOIN" for text in sky.texts)


def test_plot_merged_spectrum_qc_skylines_isolated_vs_other_colour_and_label(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Isolated skylines draw blue lines labelled "calibration skyline" (only
    on the top panel); non-isolated ones draw grey lines labelled
    "skyline" (also only on the top panel)."""
    figures = _spy_figures(monkeypatch)
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame(
            {"WAVELENGTH": [520.0, 540.0], "ISOLATED": [True, False]}
        ),
    )
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

    top = figures[-1].axes[0]
    isolatedLine = next(line for line in top.lines if line.get_xdata()[0] == pytest.approx(520.0))
    otherLine = next(line for line in top.lines if line.get_xdata()[0] == pytest.approx(540.0))
    assert isolatedLine.get_color() == "blue"
    assert isolatedLine.get_label() == "calibration skyline"
    assert otherLine.get_color() == "grey"
    assert otherLine.get_label() == "skyline"

    middle = figures[-1].axes[1]
    middleIsolatedLine = next(line for line in middle.lines if line.get_xdata()[0] == pytest.approx(520.0))
    # AN EMPTY-STRING LABEL IS TREATED BY MATPLOTLIB AS "NO LABEL", SO IT COMES
    # BACK AS AN AUTO-GENERATED "_childN" NAME RATHER THAN "" -- PIN THAT
    assert middleIsolatedLine.get_label().startswith("_child")


def test_plot_merged_spectrum_qc_skylines_drops_non_numeric_wavelengths(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`WAVELENGTH` values that fail `pd.to_numeric` are coerced to `NaN` and
    dropped, rather than raising or plotting at a bogus position."""
    figures = _spy_figures(monkeypatch)
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame(
            {"WAVELENGTH": [520.0, "not-a-number"], "ISOLATED": [True, True]}
        ),
    )
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, orderJoins=False)

    top = figures[-1].axes[0]
    skylineLines = [line for line in top.lines if line.get_color() == "blue"]
    assert len(skylineLines) == 1
    assert skylineLines[0].get_xdata()[0] == pytest.approx(520.0)


def test_plot_merged_spectrum_qc_debug_calls_show_before_savefig(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    order: list[str] = []
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: order.append("show"))
    realSavefig = plt.savefig

    def _savefig_spy(*args: object, **kwargs: object) -> object:
        order.append("savefig")
        return realSavefig(*args, **kwargs)

    monkeypatch.setattr(plt, "savefig", _savefig_spy)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _call_plot_merged_spectrum_qc(monkeypatch, tmp_path, log, merged, debug=True)

    assert order == ["show", "savefig"]


@pytest.mark.parametrize(
    ("noddingSequence", "expectedSuffix"),
    [(False, ""), ("_AB", "_AB")],
)
def test_plot_merged_spectrum_qc_filename_suffix_reflects_nodding_sequence(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, noddingSequence, expectedSuffix: str
) -> None:
    _spy_figures(monkeypatch)
    _stub_skylines(monkeypatch)
    merged = _synthetic_merged_orders()

    _, outputPath = _call_plot_merged_spectrum_qc(
        monkeypatch, tmp_path, log, merged, noddingSequence=noddingSequence
    )

    assert outputPath == str(tmp_path / f"synthetic_EXTRACTED_MERGED_QC_PLOT{expectedSuffix}.pdf")
