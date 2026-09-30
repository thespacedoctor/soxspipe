"""Characterization tests pinning the current behaviour of the lint-flagged
sites in `toolkit.py` (UP031, E712, SIM108, B905/E741, RET503, RET504, F841)
ahead of the DY-80 lint-fix commit.

These tests describe what the code does TODAY -- including behaviour that
looks like a defect (noted inline as BUG-LIKE) -- so that commit can be
verified as behaviour-preserving. They must all PASS unchanged against the
current implementation; they are not RED/GREEN tests.

Sites already pinned exactly by `test_toolkit.py`, `test_toolkit_characterization.py`,
`test_quicklook_image_characterization.py`, `test_toolkit_plotting.py`, or
`tests/integration/test_toolkit_fits.py` are not duplicated here -- see the
module docstring notes inline at each skipped site.
"""

from __future__ import annotations

import logging
import warnings
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import toolkit
from tests.factories import dispersion_map_fits, order_table_fits

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# SHARED HELPERS (LOCAL TO THIS FILE -- SEE `test_toolkit_characterization.py`
# FOR THE NEAR-IDENTICAL HELPERS USED TO PIN THE SURROUNDING QC-ROW CONTRACT)
# ---------------------------------------------------------------------------


def _identity_keyword_lookup(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(toolkit, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key))


def _detector_lookup_stub(**detectorParams: object) -> type:
    class _StubDetectorLookup:
        def __init__(self, **kwargs: object) -> None:
            pass

        def get(self, arm: str) -> dict[str, object]:
            return detectorParams

    return _StubDetectorLookup


def _unpack_order_table_stub(pixelsFrame: pd.DataFrame):
    def _stub(**kwargs: object) -> tuple[None, pd.DataFrame, None]:
        return None, pixelsFrame, None

    return _stub


def _spectroscopic_frame(*, shape: tuple[int, int] = (4, 4)) -> SimpleNamespace:
    return SimpleNamespace(
        data=np.arange(np.prod(shape), dtype=float).reshape(shape),
        mask=np.zeros(shape, dtype=bool),
        header={
            "SEQ_ARM": "VIS",
            "DATE_OBS": "2024-01-02T03:04:05",
            "INSTRUME": "SOXS",
            "WIN_BINX": 1,
            "WIN_BINY": 1,
        },
    )


# ---------------------------------------------------------------------------
# 1. UP031 -- `"%0.*f" % (3, value)` IN `spectroscopic_image_quality_checks`
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("label", "value"),
    [
        ("python_float", 2.5),
        ("np_float64", np.float64(2.5)),
        ("np_float32", np.float32(2.5)),
        ("negative", -3.14159),
        ("large", 1e20),
        ("round_half_up_2.0005", 2.0005),
        ("round_half_up_2.0015", 2.0015),
        ("nan", float("nan")),
        ("inf", float("inf")),
        ("neg_inf", float("-inf")),
    ],
)
def test_percent_format_and_fstring_format_agree_for_ordinary_floats(label: str, value: object) -> None:
    """`"%0.*f" % (3, value)` and `f"{value:.3f}"` are byte-for-byte identical
    for every ordinary numeric input reachable from `np.ma.mean`/`np.ma.sum`
    of a non-fully-masked frame -- both delegate to the same float formatter,
    including the half-up rounding of `2.0005`/`2.0015` to `2.001`/`2.002`."""
    # THE PERCENT FORM IS THE EXPRESSION UNDER TEST
    percentForm = "%0.*f" % (3, value)  # noqa: UP031
    fStringForm = f"{value:.3f}"

    assert percentForm == fStringForm


def test_percent_format_and_fstring_format_disagree_for_fully_masked_constant() -> None:
    """The reason toolkit.py keeps its `%` form under `noqa: UP031`: the two candidate rewrites do *not* agree on
    the one input this call site can actually receive from a fully-masked
    frame (`np.ma.mean`/`np.ma.sum` return `np.ma.masked`, a `MaskedConstant`).
    `"%0.*f" % (3, np.ma.masked)` silently coerces to NaN and renders `"nan"`
    (with a `UserWarning`); `f"{np.ma.masked:.3f}"` renders `"--"` (with a
    `FutureWarning`) and ignores the format spec entirely. A straight UP031
    rewrite would silently change the QC string for the fully-masked-frame
    edge case -- pinned here so the later commit can choose deliberately."""
    maskedScalar = np.ma.masked

    with pytest.warns(UserWarning, match="converting a masked element to nan"):
        # THE PERCENT FORM IS THE EXPRESSION UNDER TEST
        percentForm = "%0.*f" % (3, maskedScalar)  # noqa: UP031
    with pytest.warns(FutureWarning, match="Format strings passed to MaskedConstant"):
        fStringForm = f"{maskedScalar:.3f}"

    assert percentForm == "nan"
    assert fStringForm == "--"
    assert percentForm != fStringForm


@pytest.mark.parametrize(
    ("dispersionAxis", "pixelColumns", "badPixel", "expectedMean", "expectedSum"),
    [
        (
            "x",
            {"xcoord_edgeup": [3, 4], "xcoord_edgelow": [1, 0], "ycoord": [1, 2]},
            (2, 0),
            "8.200",
            "41.000",
        ),
        (
            "y",
            {"ycoord_edgeup": [3, 4], "ycoord_edgelow": [1, 0], "xcoord": [1, 2]},
            (0, 2),
            "8.800",
            "44.000",
        ),
    ],
)
def test_spectroscopic_image_quality_checks_qc_value_matches_percent_format(
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    dispersionAxis: str,
    pixelColumns: dict[str, list[float]],
    badPixel: tuple[int, int],
    expectedMean: str,
    expectedSum: str,
) -> None:
    """Restates the table pinned in `test_toolkit_characterization.py` as a
    direct cross-check against the standalone `"%0.*f"` table above, so the
    two live side by side and cannot silently drift apart."""
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub(**{"dispersion-axis": dispersionAxis}))
    monkeypatch.setattr(toolkit, "unpack_order_table", _unpack_order_table_stub(pd.DataFrame(pixelColumns)))
    frame = _spectroscopic_frame()
    frame.mask[badPixel] = True

    result = toolkit.spectroscopic_image_quality_checks(
        log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
    )

    assert result["qc_value"].tolist() == [expectedMean, expectedSum]


def test_spectroscopic_image_quality_checks_reports_nan_string_for_fully_masked_frame(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """When no pixel is ever unmasked (an empty order-pixel table), the
    combined mask excludes every pixel, so `np.ma.mean`/`np.ma.sum` return
    `np.ma.masked` and the QC value ends up as the literal string `"nan"`
    (via the `UserWarning`-emitting coercion pinned above) -- not a crash,
    not a `NaN` float, a `"nan"` string. Pinned as-is."""
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub(**{"dispersion-axis": "x"}))
    emptyPixels = pd.DataFrame({"xcoord_edgeup": [], "xcoord_edgelow": [], "ycoord": []})
    monkeypatch.setattr(toolkit, "unpack_order_table", _unpack_order_table_stub(emptyPixels))
    frame = _spectroscopic_frame()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = toolkit.spectroscopic_image_quality_checks(
            log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
        )

    assert result["qc_value"].tolist() == ["nan", "nan"]
    assert all(isinstance(value, str) for value in result["qc_value"])
    assert any("converting a masked element to nan" in str(w.message) for w in caught)


# ---------------------------------------------------------------------------
# 2. E712 -- `mapDF["mask"] == False` IN `twoD_disp_map_image_to_dataframe`
# ---------------------------------------------------------------------------

# NOTE: THE PRODUCTION BEHAVIOUR OF `removeMaskedPixels=True` (MASKED PIXELS
# REMOVED, UNMASKED PIXELS KEPT) IS ALREADY PINNED END-TO-END AGAINST A REAL
# FITS-BACKED FRAME BY
# `tests/integration/test_toolkit_fits.py::test_two_d_map_dataframe_preserves_associated_frame_arrays`
# -- NOT DUPLICATED HERE. WHAT IS PINNED BELOW IS THE GENERAL EQUIVALENCE
# QUESTION THE LATER E712 REWRITE NEEDS ANSWERED: DOES `Series == False`
# AGREE WITH `Series.eq(False)` ACROSS THE DTYPES `mapDF["mask"]` CAN
# PLAUSIBLY CARRY?


@pytest.mark.parametrize(
    ("label", "series"),
    [
        ("bool", pd.Series([True, False, True])),
        ("int_zero_one", pd.Series([1, 0, 1])),
        ("float_zero_one", pd.Series([1.0, 0.0, 1.0])),
        ("object_with_none", pd.Series([True, None, False], dtype=object)),
    ],
)
def test_series_equals_false_agrees_with_series_eq_false(label: str, series: pd.Series) -> None:
    """`Series == False` (E712) and `Series.eq(False)` produce identical
    boolean masks for bool, 0/1 int, 0.0/1.0 float, and object-with-`None`
    dtypes -- including `None`, which compares not-equal-to `False` under
    plain Python `==` and so is *not* selected by the `mask` filter."""
    equalityOperator = series == False  # noqa: E712 -- INTENTIONALLY PINNING THE FLAGGED FORM
    equalityMethod = series.eq(False)

    pd.testing.assert_series_equal(equalityOperator, equalityMethod, check_names=False)


# ---------------------------------------------------------------------------
# 3A. SIM108 -- `halfwidth` IN `cut_image_slice`
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("width", "expectedMedians", "expectedOffset", "expectedCentre"),
    [
        # ODD WIDTH: HALFWIDTH = (WIDTH - 1) / 2 = 1
        (3, [22.0, 23.0, 24.0, 25.0], 1, 3.5),
        # EVEN WIDTH: HALFWIDTH = WIDTH / 2 = 2.0 -- SAME RESULT HERE BECAUSE
        # THE SYNTHETIC FRAME IS A LINEAR RAMP, SYMMETRIC AROUND THE CENTRE
        (4, [22.0, 23.0, 24.0, 25.0], 1, 3.5),
        # FLOAT WIDTH (ODD BRANCH, SINCE 3.5 % 2 != 0): HALFWIDTH = 1.25 --
        # THE INT() TRUNCATION IN `slice_width_centre` THEN SHIFTS THE
        # REPORTED CENTRE TO 3.0 INSTEAD OF 3.5
        (3.5, [18.5, 19.5, 20.5, 21.5], 1, 3.0),
    ],
)
def test_cut_image_slice_halfwidth_across_odd_even_and_float_widths(
    log: object,
    width: float,
    expectedMedians: list[float],
    expectedOffset: int,
    expectedCentre: float,
) -> None:
    frame = np.arange(49, dtype=float).reshape(7, 7)

    medianSlice, offset, centre = toolkit.cut_image_slice(log, frame, width=width, length=4, x=3, y=3, median=True)

    np.testing.assert_array_equal(np.ma.filled(medianSlice, np.nan), expectedMedians)
    assert (offset, centre) == (expectedOffset, expectedCentre)


# ---------------------------------------------------------------------------
# 3B. SIM108 -- `shrink` IN `quicklook_image`
# ---------------------------------------------------------------------------


def test_quicklook_image_colorbar_shrink_is_smaller_for_surface_plots(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`shrink = 0.5` when `surfacePlot=True`, else `shrink = 1.0` -- pinned
    by spying on `Figure.colorbar` and inspecting the `shrink` kwarg it is
    called with."""
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    calls: list[dict[str, object]] = []
    originalColorbar = Figure.colorbar

    def _spy(self: Figure, *args: object, **kwargs: object) -> object:
        calls.append(kwargs)
        return originalColorbar(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "colorbar", _spy)
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: None)
    frame = np.zeros((8, 8))

    toolkit.quicklook_image(log, frame, show=True, inst="OTHER", surfacePlot=True)
    toolkit.quicklook_image(log, frame, show=True, inst="OTHER", surfacePlot=False)

    assert calls[0]["shrink"] == 0.5
    assert calls[1]["shrink"] == 1.0


# ---------------------------------------------------------------------------
# 3C / 4A. SIM108 `ratio` + B905 AXIS-COORDINATE ZIP/CLAMPS IN
# `unpack_order_table`
# ---------------------------------------------------------------------------


def _order_table_path(tmp_path, *, metadata: pd.DataFrame):
    polynomials = pd.DataFrame(
        [
            {
                "degorder_cent": 1,
                "degy_cent": 1,
                "cent_00": 1.0,
                "cent_01": 2.0,
                "cent_10": 3.0,
                "cent_11": 4.0,
                "std_00": 2.0,
                "std_01": 0.0,
                "std_10": 0.0,
                "std_11": 0.0,
            }
        ]
    )
    return order_table_fits(tmp_path / "orders.fits", polynomials=polynomials, metadata=metadata)


def test_unpack_order_table_prebinned_ratio_changes_candidate_pixel_range(tmp_path, log: object) -> None:
    """`ratio = axisBbin if prebinned else 1` (the SIM108 `ratio` site)
    changes the raw `ymin`/`ymax` bounds fed into the coordinate `arange`
    before the final `/= axisBbin` division at the end of the function, so
    `prebinned=True` and `prebinned=False` give different pixel ranges for
    the same metadata table and `biny`."""
    metadata = pd.DataFrame({"order": [10, 11], "ymin": [0.0, 2.0], "ymax": [8.0, 10.0]})
    orderPath = _order_table_path(tmp_path, metadata=metadata)

    _, notPrebinned, _ = toolkit.unpack_order_table(
        log=log, orderTablePath=str(orderPath), pixelDelta=1, biny=2, prebinned=False, order=11
    )
    _, prebinned, _ = toolkit.unpack_order_table(
        log=log, orderTablePath=str(orderPath), pixelDelta=1, biny=2, prebinned=True, order=11
    )

    assert list(notPrebinned["ycoord"]) == [1, 2, 3, 4]
    assert list(prebinned["ycoord"]) == [2, 3, 4, 5, 6, 7, 8, 9]


def test_unpack_order_table_clamps_axis_coordinates_to_0_and_4200(tmp_path, log: object) -> None:
    """The B905 zip site building `axisBcoords` clamps the lower bound to `0`
    (`math.floor(l) - int(r * extend) < 0`) and the upper bound to `4200`
    (`math.ceil(u) + int(r * extend) > 4200`) independently per order. A
    large `extend` pushes order 10's lower bound below zero and order 11's
    upper bound above 4200, exercising both clamps in one call."""
    metadata = pd.DataFrame({"order": [10, 11], "ymin": [1.0, 4190.0], "ymax": [3.0, 4199.0]})
    orderPath = _order_table_path(tmp_path, metadata=metadata)

    _, pixelTable, _ = toolkit.unpack_order_table(log=log, orderTablePath=str(orderPath), pixelDelta=1, extend=100.0)

    order10 = pixelTable.loc[pixelTable["order"] == 10, "ycoord"]
    order11 = pixelTable.loc[pixelTable["order"] == 11, "ycoord"]
    # ORDER 10: UNCLAMPED LOWER BOUND WOULD BE FLOOR(1) - INT(2*100) = -199 --
    # CLAMPED TO 0
    assert order10.min() == 0
    assert order10.max() == 202
    # ORDER 11: UNCLAMPED UPPER BOUND WOULD BE CEIL(4199) + INT(9*100) = 5099
    # -- CLAMPED TO 4200 (EXCLUSIVE), SO THE LAST POINT IS 4199
    assert order11.min() == 3290
    assert order11.max() == 4199


# ---------------------------------------------------------------------------
# 4B. B905/E741 -- MASK LOOP IN `spectroscopic_image_quality_checks`
# ---------------------------------------------------------------------------


def test_spectroscopic_image_quality_checks_clamps_and_skips_out_of_range_rows_axis_x(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """For `axisA == "x"`, the `zip(axisACoords_up, axisACoords_low,
    axisBCoords)` loop: (1) clamps `l` up to `0` and `u` down to
    `mask.shape[1]`; (2) silently skips any row whose `y` falls outside
    `mask.shape[0]` (no crash, no unmasking); (3) silently skips any row
    where the clamped `l >= u` (no unmasking). Only the row that survives
    all three checks unmasks its pixels."""
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub(**{"dispersion-axis": "x"}))
    pixelsFrame = pd.DataFrame(
        {
            # ROW 0: CLAMPS L FROM -5 TO 0 AND U FROM 10 TO 4, Y=1 VALID
            # ROW 1: Y=8 IS OUT OF BOUNDS (MASK.SHAPE[0] == 4) -- SKIPPED
            # ROW 2: L=3 >= U=1 AFTER NO CLAMPING NEEDED -- SKIPPED
            "xcoord_edgeup": [10, 3, 1],
            "xcoord_edgelow": [-5, 1, 3],
            "ycoord": [1, 8, 2],
        }
    )
    monkeypatch.setattr(toolkit, "unpack_order_table", _unpack_order_table_stub(pixelsFrame))
    frame = _spectroscopic_frame()

    result = toolkit.spectroscopic_image_quality_checks(
        log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
    )

    # ONLY ROW 1 OF THE 4X4 ARANGE(16) FRAME ([4, 5, 6, 7]) IS UNMASKED
    assert result["qc_value"].tolist() == ["5.500", "22.000"]


def test_spectroscopic_image_quality_checks_clamps_and_skips_out_of_range_rows_axis_y(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """Mirror of the `axisA == "x"` case above for `axisA == "y"`, where the
    loop clamps/unmasks columns instead of rows."""
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub(**{"dispersion-axis": "y"}))
    pixelsFrame = pd.DataFrame(
        {
            "ycoord_edgeup": [10, 3, 1],
            "ycoord_edgelow": [-5, 1, 3],
            "xcoord": [1, 8, 2],
        }
    )
    monkeypatch.setattr(toolkit, "unpack_order_table", _unpack_order_table_stub(pixelsFrame))
    frame = _spectroscopic_frame()

    result = toolkit.spectroscopic_image_quality_checks(
        log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
    )

    # ONLY COLUMN 1 OF THE 4X4 ARANGE(16) FRAME ([1, 5, 9, 13]) IS UNMASKED
    assert result["qc_value"].tolist() == ["7.000", "28.000"]


# ---------------------------------------------------------------------------
# 4C. B905 -- `read_spectral_format` WITH A `dispersionMap`
# ---------------------------------------------------------------------------


def test_read_spectral_format_with_dispersion_map_clamps_to_science_pixels(
    tmp_path, monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """The `zip(orderNums, waveLengthMin, waveLengthMax)` loop feeding the
    `dispersionMap` branch produces one `amin`/`amax` pair per order, clamped
    to `[0, science-pixels[rowCol]["end"]]`. Exercised end-to-end against a
    real spectral-format FITS table and a real dispersion-map FITS table
    (`tests/factories.dispersion_map_fits`)."""
    _identity_keyword_lookup(monkeypatch)
    detector = {
        "science-pixels": {"rows": {"end": 31}, "columns": {"end": 31}},
        "spectral format table": "format.fits",
        "dispersion-axis": "x",
    }
    monkeypatch.setattr(toolkit, "detector_lookup", lambda **kwargs: SimpleNamespace(get=lambda arm: detector))
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **kwargs: str(tmp_path))
    Table(
        {
            "ORDER": [10, 11],
            "WLMIN": [500.0, 600.0],
            "WLMAX": [502.0, 602.0],
            "WLMINFUL": [500.0, 600.0],
            "WLMAXFUL": [502.0, 602.0],
        }
    ).write(tmp_path / "format.fits")
    dispersionMapPath = dispersion_map_fits(tmp_path / "disp.fits")

    orderNums, waveMin, waveMax, amins, amaxs = toolkit.read_spectral_format(
        log, {}, "VIS", dispersionMap=str(dispersionMapPath), extended=False
    )

    assert list(orderNums) == [10, 11]
    assert list(waveMin) == [500.0, 600.0]
    assert list(waveMax) == [502.0, 602.0]
    # BOTH ORDERS' FIT_Y RANGE EXCEEDS THE 31-PIXEL "ROWS" SCIENCE-PIXEL
    # LIMIT, SO AMAX CLAMPS TO 31.0 FOR BOTH -- AMIN IS THE UNCLAMPED
    # POLYNOMIAL EVALUATION AT EACH ORDER'S MINIMUM WAVELENGTH
    assert amins == pytest.approx([112.2, 135.42], rel=1e-6)
    assert amaxs == [31.0, 31.0]


# ---------------------------------------------------------------------------
# 5. RET503 -- `MaxFilter.filter`
# ---------------------------------------------------------------------------


def test_max_filter_returns_none_not_false_for_records_above_the_limit() -> None:
    """`test_toolkit.py` already pins `True` below the limit and `is None` at
    the limit; this adds the `above`-the-limit case to complete the domain,
    and is explicit that the falsy return is `None`, not `False`."""
    filterObject = toolkit.MaxFilter(logging.WARNING)

    result = filterObject.filter(logging.LogRecord("x", logging.ERROR, "", 1, "", (), None))

    assert result is None
    assert result is not False


# ---------------------------------------------------------------------------
# 6. RET504 -- `extinction_correction_factor`
# ---------------------------------------------------------------------------
#
# ALREADY PINNED END-TO-END AGAINST A REAL FITS EXTINCTION TABLE BY
# `tests/integration/test_toolkit_fits.py::test_extinction_correction_reads_fits_table_and_uses_next_sample`
# -- NOT DUPLICATED HERE.


# ---------------------------------------------------------------------------
# 7. F841 -- DEAD-STORE ASSIGNMENTS WHOSE RIGHT-HAND SIDE CAN RAISE
# ---------------------------------------------------------------------------


def test_quicklook_image_raises_key_error_for_missing_date_obs_header(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`dateObs = CCDObject.header[kw("DATE_OBS")]` (L191) is otherwise
    unused, but a header missing `DATE_OBS` still raises `KeyError` today --
    deleting the dead assignment would silently swallow this."""
    monkeypatch.setattr(commonutils, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key))
    monkeypatch.setattr(commonutils, "detector_lookup", lambda **kwargs: SimpleNamespace(get=lambda arm: {}))
    frame = SimpleNamespace(header={"SEQ_ARM": "VIS"}, data=np.zeros((4, 4)), mask=np.zeros((4, 4), dtype=bool))

    with pytest.raises(KeyError, match="DATE_OBS"):
        toolkit.quicklook_image(log, frame, show=False, saveToPath="unused.pdf", settings={"instrument": "soxs"})


def test_quicklook_image_raises_key_error_for_missing_science_pixels(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`science_pixels = dp["science-pixels"]` (L198) is otherwise unused,
    but a detector-parameters dict missing that key still raises `KeyError`
    today."""
    monkeypatch.setattr(commonutils, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key))
    monkeypatch.setattr(commonutils, "detector_lookup", lambda **kwargs: SimpleNamespace(get=lambda arm: {}))
    frame = SimpleNamespace(
        header={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-01"},
        data=np.zeros((4, 4)),
        mask=np.zeros((4, 4), dtype=bool),
    )

    with pytest.raises(KeyError, match="science-pixels"):
        toolkit.quicklook_image(log, frame, show=False, saveToPath="unused.pdf", settings={"instrument": "soxs"})


def test_generic_quality_checks_raises_key_error_for_missing_seq_arm_header(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`arm = frame.header[kw("SEQ_ARM")]` (L682) is otherwise unused before
    `.lower()` is called on it, but a header missing `SEQ_ARM` still raises
    `KeyError` today."""
    _identity_keyword_lookup(monkeypatch)
    frame = SimpleNamespace(mask=np.zeros((2, 2), dtype=bool), header={})

    with pytest.raises(KeyError, match="SEQ_ARM"):
        toolkit.generic_quality_checks(log, frame, {}, "soxs-stare", pd.DataFrame())


def test_spectroscopic_image_quality_checks_raises_key_error_for_missing_instrume_header(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`inst = frame.header[kw("INSTRUME")]` (L801) is otherwise unused, but
    a header missing `INSTRUME` (while `SEQ_ARM`/`DATE_OBS`/`WIN_BINX`/
    `WIN_BINY` are present, isolating this specific site) still raises
    `KeyError` today."""
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub(**{"dispersion-axis": "x"}))
    frame = SimpleNamespace(
        data=np.zeros((4, 4)),
        mask=np.zeros((4, 4), dtype=bool),
        header={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-01", "WIN_BINX": 1, "WIN_BINY": 1},
    )

    with pytest.raises(KeyError, match="INSTRUME"):
        toolkit.spectroscopic_image_quality_checks(
            log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
        )


def test_read_spectral_format_raises_key_error_for_missing_science_pixels(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    """`science_pixels = dp["science-pixels"]` (L938) is otherwise unused,
    but a detector-parameters dict missing that key still raises `KeyError`
    today, before the spectral-format FITS file is ever read."""
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(
        toolkit,
        "detector_lookup",
        lambda **kwargs: SimpleNamespace(get=lambda arm: {"spectral format table": "x.fits"}),
    )

    with pytest.raises(KeyError, match="science-pixels"):
        toolkit.read_spectral_format(log, {}, "VIS")


def test_get_calibration_lamp_raises_key_error_for_missing_instrume_header() -> None:
    """`inst = frame.header["INSTRUME"]` (L1493) is otherwise unused, but a
    header missing `INSTRUME` still raises `KeyError` today."""
    with pytest.raises(KeyError, match="INSTRUME"):
        toolkit.get_calibration_lamp(logging.getLogger("x"), SimpleNamespace(header={}), lambda name: name)


# ---------------------------------------------------------------------------
# GENERAL B905 NOTE -- `strict=False` TRUNCATION SEMANTICS
# ---------------------------------------------------------------------------


def test_zip_without_strict_silently_truncates_to_the_shortest_iterable() -> None:
    """None of the B905-flagged `zip()` call sites above have production
    inputs of genuinely unequal length reachable today -- every zipped tuple
    of lists is built from the same source table/array, so lengths already
    match. This documents the semantics a `strict=False` (i.e. unchanged)
    rewrite must preserve, and what `strict=True` would newly reject, for
    whichever of those sites the later commit chooses not to make strict."""
    longer = [1, 2, 3, 4]
    shorter = ["a", "b"]

    # THE IMPLICIT DEFAULT IS THE BEHAVIOUR UNDER TEST
    truncated = list(zip(longer, shorter))  # noqa: B905

    assert truncated == [(1, "a"), (2, "b")]
    with pytest.raises(ValueError, match="zip"):
        list(zip(longer, shorter, strict=True))
