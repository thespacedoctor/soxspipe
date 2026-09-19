"""Characterization tests freezing the current rendering behaviour of
`quicklook_image` in `soxspipe/commonutils/toolkit.py`, ahead of the DY-79 split of
that function into smaller helpers.

These tests pin current behaviour, including behaviour that looks like a
defect (noted inline as BUG-LIKE). They pin the matplotlib artists produced
-- figure sizes, axes, line data, colours, text and axis limits -- by spying
on `matplotlib.pyplot.figure` and inspecting the closed-but-inspectable
`Figure` objects it returns. The PDF output and the savefig/clf/close call
sequence are pinned in `tests/unit/test_toolkit_characterization.py`.
"""

from __future__ import annotations

import matplotlib
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.stats import sigma_clipped_stats
from matplotlib.figure import Figure

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import toolkit
from tests.unit._plot_spies import quiet_show, spy_figures

pytestmark = pytest.mark.unit


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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)

    toolkit.quicklook_image(log, np.zeros((8, 8)), show=True, inst="OTHER", title="My Title")

    suptitle = figures[-1]._suptitle
    assert suptitle.get_text() == "My Title"
    assert suptitle.get_fontsize() == 20


def test_quicklook_image_without_title_has_no_suptitle(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)

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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)

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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
    frame = np.arange(64, dtype=float).reshape(8, 8)

    toolkit.quicklook_image(log, frame, show=True, ext=None)

    image = _image_axis(figures[-1]).images[0]
    np.testing.assert_array_equal(image.get_array().data, np.flipud(np.rot90(frame, 1)))


def test_quicklook_image_inst_argument_overrides_header_instrume(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """An explicit `inst=` argument takes precedence over `header['INSTRUME']`."""
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    figures = spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    spy_figures(monkeypatch)
    quiet_show(monkeypatch)
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
    spy_figures(monkeypatch)
    showCalls = quiet_show(monkeypatch)

    toolkit.quicklook_image(log, np.zeros((8, 8)), show=True, inst="OTHER")

    assert len(showCalls) == 1
