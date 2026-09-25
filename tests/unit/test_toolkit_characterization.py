"""Characterization tests pinning the current behaviour of `toolkit.py` ahead
of the DY-77 refactor.

These tests describe what the four naive `datetime.utcnow()`/`datetime.now(UTC)`
timestamp sites and the two hand-rolled `plt.savefig(...)` call sites in
`generic_quality_checks`, `spectroscopic_image_quality_checks`,
`add_snr_efficiency_qcs`, `plot_merged_spectrum_qc`, and `quicklook_image` do
today, before they are swapped for `toolkit.utcnow_string()` and
`toolkit.save_qc_plot()`. They must all pass unchanged against the current
implementation -- they are not RED/GREEN tests, they exist to freeze
behaviour (including behaviour that looks like a defect, noted inline) ahead
of the refactor.
"""

from __future__ import annotations

import datetime as datetime_module
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import toolkit
from tests.factories import product_table, qc_row

pytestmark = pytest.mark.unit

# THE FROZEN INSTANT CARRIES MICROSECONDS SO A CHANGE THAT STARTS RENDERING
# `%f` (E.G. SWITCHING TO `utcnow_string(microseconds=True)`) IS CAUGHT BY THE
# EXACT-STRING ASSERTIONS BELOW
_FROZEN_INSTANT = datetime_module.datetime(2025, 6, 7, 8, 9, 10, 123456)
_FROZEN_TIMESTAMP = "2025-06-07T08:09:10"
_TICK_START = datetime_module.datetime(2025, 1, 1, 0, 0, 0)
_TIMESTAMP_PATTERN = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d$")


def _build_controlled_datetime(
    start: datetime_module.datetime, *, tick: bool
) -> type[datetime_module.datetime]:
    """Return a `datetime.datetime` subclass whose `utcnow()`/`now()` are controlled.

    **Key Arguments:**

    - ``start`` -- the first instant the controlled clock returns
    - ``tick`` -- advance the returned instant by one second on every call
      (frozen when *False*)

    **Return:**

    - a `datetime.datetime` subclass suitable for monkeypatching both the
      stdlib `datetime` module attribute (so function-local
      ``from datetime import datetime`` imports pick it up) and
      `toolkit.datetime` (the module-level import)
    """
    realTimedelta = datetime_module.timedelta

    class _ControlledDatetime(datetime_module.datetime):
        _next_value = start

        @classmethod
        def _advance(cls) -> datetime_module.datetime:
            current = cls._next_value
            if tick:
                cls._next_value = current + realTimedelta(seconds=1)
            return current

        @classmethod
        def utcnow(cls) -> datetime_module.datetime:
            return cls._advance()

        @classmethod
        def now(cls, tz: object = None) -> datetime_module.datetime:
            value = cls._advance()
            if tz is not None:
                return value.replace(tzinfo=tz)
            return value

    return _ControlledDatetime


@pytest.fixture
def frozen_clock(
    monkeypatch: pytest.MonkeyPatch,
) -> type[datetime_module.datetime]:
    """Patch both the stdlib `datetime.datetime` and `toolkit.datetime` to a fixed instant."""
    controlled = _build_controlled_datetime(_FROZEN_INSTANT, tick=False)
    monkeypatch.setattr(datetime_module, "datetime", controlled)
    monkeypatch.setattr(toolkit, "datetime", controlled)
    # VERIFY THE PATCH TOOK EFFECT BEFORE ANY TEST RELIES ON IT
    assert controlled.utcnow() == _FROZEN_INSTANT
    assert controlled.now(tz=None) == _FROZEN_INSTANT
    return controlled


@pytest.fixture
def ticking_clock(
    monkeypatch: pytest.MonkeyPatch,
) -> type[datetime_module.datetime]:
    """Patch the clock so each call advances by one second, catching per-row minting."""
    controlled = _build_controlled_datetime(_TICK_START, tick=True)
    monkeypatch.setattr(datetime_module, "datetime", controlled)
    monkeypatch.setattr(toolkit, "datetime", controlled)
    # VERIFY THE PATCH TOOK EFFECT (AND ACTUALLY TICKS) BEFORE THE TEST RUNS,
    # THEN RESET THE COUNTER SO THE TEST OBSERVES A CLEAN SEQUENCE
    firstTick = controlled.utcnow()
    secondTick = controlled.utcnow()
    assert secondTick - firstTick == datetime_module.timedelta(seconds=1)
    controlled._next_value = _TICK_START
    return controlled


def _record_calls(
    monkeypatch: pytest.MonkeyPatch, attributeNames: list[str]
) -> dict[str, list[tuple[tuple, dict]]]:
    """Patch each named `matplotlib.pyplot` attribute, delegating to the real
    implementation while recording call order and arguments.

    **Key Arguments:**

    - ``monkeypatch`` -- the pytest monkeypatch fixture
    - ``attributeNames`` -- the `matplotlib.pyplot` attribute names to spy on

    **Return:**

    - a mapping of attribute name to its list of `(args, kwargs)` calls, plus
      a special ``"__order__"`` key holding the interleaved call order
    """
    import matplotlib.pyplot as plt

    order: list[str] = []
    calls: dict[str, list[tuple[tuple, dict]]] = {name: [] for name in attributeNames}

    for attributeName in attributeNames:
        original = getattr(plt, attributeName)

        def _wrapper(*args: object, __original=original, __name=attributeName, **kwargs: object) -> object:
            order.append(__name)
            calls[__name].append((args, kwargs))
            return __original(*args, **kwargs)

        monkeypatch.setattr(plt, attributeName, _wrapper)

    calls["__order__"] = order
    return calls


# ---------------------------------------------------------------------------
# generic_quality_checks
# ---------------------------------------------------------------------------


def _identity_keyword_lookup(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        toolkit, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key)
    )


def _bad_pixel_frame() -> SimpleNamespace:
    """Return a synthetic frame with 1 bad pixel out of 7 (a non-trivial mask fraction)."""
    return SimpleNamespace(
        mask=np.array([True, False, False, False, False, False, False]),
        header={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02T03:04:05"},
    )


def test_generic_quality_checks_pins_full_row_contract_under_frozen_clock(
    monkeypatch: pytest.MonkeyPatch, log: object, frozen_clock: object
) -> None:
    _identity_keyword_lookup(monkeypatch)
    frame = _bad_pixel_frame()

    result = toolkit.generic_quality_checks(log, frame, {}, "soxs-stare", pd.DataFrame())

    assert list(result.columns) == [
        "soxspipe_recipe",
        "qc_name",
        "qc_value",
        "qc_comment",
        "qc_unit",
        "obs_date_utc",
        "reduction_date_utc",
        "to_header",
    ]
    assert result["soxspipe_recipe"].tolist() == ["soxs-stare", "soxs-stare"]
    assert result["qc_name"].tolist() == ["BADPIX NUM", "BADPIX FRAC"]
    assert result["qc_value"].iloc[0] == pytest.approx(1.0, rel=1e-12)
    # THE FRACTION IS ROUNDED TO SIX DECIMALS BY `float(f"{percent:.6f}")`
    assert result["qc_value"].iloc[1] == pytest.approx(0.142857, rel=1e-12)
    # BUG-LIKE: THE CODE DOES `int(badCount)`, BUT CONCATENATING THAT ROW WITH
    # THE FLOAT FRACTION ROW UPCASTS THE WHOLE `qc_value` COLUMN TO float64, SO
    # THE COUNT IS NOT ACTUALLY STORED AS A PYTHON `int` -- PINNED AS-IS
    assert isinstance(result["qc_value"].iloc[0], float)
    assert result["qc_comment"].tolist() == [
        "Number of bad pixels",
        "Fraction of bad pixels",
    ]
    assert result["qc_unit"].tolist() == ["", ""]
    assert result["obs_date_utc"].tolist() == [
        "2024-01-02T03:04:05",
        "2024-01-02T03:04:05",
    ]
    assert result["reduction_date_utc"].tolist() == [_FROZEN_TIMESTAMP, _FROZEN_TIMESTAMP]
    assert result["to_header"].tolist() == [True, True]


def test_generic_quality_checks_shares_single_timestamp_across_rows_under_ticking_clock(
    monkeypatch: pytest.MonkeyPatch, log: object, ticking_clock: object
) -> None:
    _identity_keyword_lookup(monkeypatch)
    frame = _bad_pixel_frame()

    result = toolkit.generic_quality_checks(log, frame, {}, "soxs-stare", pd.DataFrame())

    timestamps = set(result["reduction_date_utc"])
    assert len(timestamps) == 1
    (timestamp,) = timestamps
    assert _TIMESTAMP_PATTERN.match(timestamp)


def test_generic_quality_checks_appends_to_existing_qc_table_without_mutating_it(
    monkeypatch: pytest.MonkeyPatch, log: object, frozen_clock: object
) -> None:
    _identity_keyword_lookup(monkeypatch)
    frame = _bad_pixel_frame()
    priorRow = qc_row(
        recipeName="soxs-mbias",
        name="RON",
        value=3.2,
        unit="electron",
        comment="Synthetic read noise",
    )
    priorRowSnapshot = priorRow.copy(deep=True)

    result = toolkit.generic_quality_checks(log, frame, {}, "soxs-stare", priorRow)

    pd.testing.assert_frame_equal(priorRow, priorRowSnapshot)
    assert len(result) == 3
    assert result.iloc[0]["qc_name"] == "RON"
    assert result.iloc[1]["qc_name"] == "BADPIX NUM"
    assert result.iloc[2]["qc_name"] == "BADPIX FRAC"


# ---------------------------------------------------------------------------
# spectroscopic_image_quality_checks
# ---------------------------------------------------------------------------


def _detector_lookup_stub(dispersionAxis: str) -> type:
    class _StubDetectorLookup:
        def __init__(self, **kwargs: object) -> None:
            pass

        def get(self, arm: str) -> dict[str, str]:
            return {"dispersion-axis": dispersionAxis}

    return _StubDetectorLookup


def _unpack_order_table_stub(pixelsFrame: pd.DataFrame):
    def _stub(**kwargs: object) -> tuple[None, pd.DataFrame, None]:
        return None, pixelsFrame, None

    return _stub


def _order_frame(
    *, badPixel: tuple[int, int], winBin: bool = True, arm: str = "VIS"
) -> SimpleNamespace:
    header: dict[str, object] = {
        "SEQ_ARM": arm,
        "DATE_OBS": "2024-01-02T03:04:05",
        "INSTRUME": "SOXS",
    }
    if winBin:
        header["WIN_BINX"] = 1
        header["WIN_BINY"] = 1
    data = np.arange(16, dtype=float).reshape(4, 4)
    maskArray = np.zeros((4, 4), dtype=bool)
    maskArray[badPixel] = True
    return SimpleNamespace(data=data, mask=maskArray, header=header)


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
def test_spectroscopic_image_quality_checks_pins_full_row_contract_per_dispersion_axis(
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    frozen_clock: object,
    dispersionAxis: str,
    pixelColumns: dict[str, list[float]],
    badPixel: tuple[int, int],
    expectedMean: str,
    expectedSum: str,
) -> None:
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub(dispersionAxis))
    monkeypatch.setattr(
        toolkit, "unpack_order_table", _unpack_order_table_stub(pd.DataFrame(pixelColumns))
    )
    frame = _order_frame(badPixel=badPixel)

    result = toolkit.spectroscopic_image_quality_checks(
        log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
    )

    assert list(result.columns) == [
        "soxspipe_recipe",
        "qc_name",
        "qc_value",
        "qc_comment",
        "qc_unit",
        "obs_date_utc",
        "reduction_date_utc",
        "to_header",
    ]
    assert result["soxspipe_recipe"].tolist() == ["soxs-stare", "soxs-stare"]
    assert result["qc_name"].tolist() == ["INNER ORDER PIX MEAN", "INNER ORDER PIX SUM"]
    assert result["qc_value"].tolist() == [expectedMean, expectedSum]
    # BUG-LIKE: THESE ARE NUMERIC MEASUREMENTS BUT ARE RECORDED AS FORMATTED
    # STRINGS (`"%0.*f" % (3, value)`), NOT FLOATS -- PINNED AS-IS
    assert isinstance(result["qc_value"].iloc[0], str)
    assert isinstance(result["qc_value"].iloc[1], str)
    assert result["qc_comment"].tolist() == [
        "[e-] Mean inner-order pixel value",
        "[e-] Sum of all inner-order pixel values",
    ]
    assert result["qc_unit"].tolist() == ["electrons", "electrons"]
    assert result["obs_date_utc"].tolist() == [
        "2024-01-02T03:04:05",
        "2024-01-02T03:04:05",
    ]
    assert result["reduction_date_utc"].tolist() == [_FROZEN_TIMESTAMP, _FROZEN_TIMESTAMP]
    assert result["to_header"].tolist() == [True, True]


def test_spectroscopic_image_quality_checks_falls_back_to_unbinned_pixels_for_nir_arm_missing_win_binx(
    monkeypatch: pytest.MonkeyPatch, log: object, frozen_clock: object
) -> None:
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub("y"))
    pixelsFrame = pd.DataFrame(
        {"ycoord_edgeup": [3, 4], "ycoord_edgelow": [1, 0], "xcoord": [1, 2]}
    )
    monkeypatch.setattr(toolkit, "unpack_order_table", _unpack_order_table_stub(pixelsFrame))
    frame = _order_frame(badPixel=(0, 2), winBin=False, arm="NIR")

    result = toolkit.spectroscopic_image_quality_checks(
        log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
    )

    assert result["qc_value"].tolist() == ["8.800", "44.000"]


def test_spectroscopic_image_quality_checks_raises_for_non_nir_arm_missing_win_binx(
    monkeypatch: pytest.MonkeyPatch, log: object, frozen_clock: object
) -> None:
    """BUG-LIKE: `binx`/`biny` are only defaulted to *1* in the `except KeyError`
    branch when `arm.lower() == "nir"`. Any other arm whose header omits
    `WIN_BINX`/`WIN_BINY` leaves both names unbound, so the function raises
    `UnboundLocalError` instead of falling back -- pinned here, not fixed."""
    _identity_keyword_lookup(monkeypatch)
    monkeypatch.setattr(commonutils, "detector_lookup", _detector_lookup_stub("x"))
    pixelsFrame = pd.DataFrame(
        {"xcoord_edgeup": [3, 4], "xcoord_edgelow": [1, 0], "ycoord": [1, 2]}
    )
    monkeypatch.setattr(toolkit, "unpack_order_table", _unpack_order_table_stub(pixelsFrame))
    frame = _order_frame(badPixel=(2, 0), winBin=False, arm="VIS")

    with pytest.raises(UnboundLocalError):
        toolkit.spectroscopic_image_quality_checks(
            log, frame, "unused-order-table.fits", {}, "soxs-stare", pd.DataFrame()
        )


# ---------------------------------------------------------------------------
# add_snr_efficiency_qcs
# ---------------------------------------------------------------------------

_EXPECTED_SNR_EFFICIENCY_QC_COLUMNS = [
    "soxspipe_recipe",
    "qc_name",
    "qc_value",
    "qc_comment",
    "qc_unit",
    "obs_date_utc",
    "reduction_date_utc",
    "to_header",
    "qc_order",
]


def _snr_efficiency_spectrum() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "WAVE": [400.0, 450.0, 550.0, 600.0],
            "SNR": [10.0, 20.0, 30.0, 40.0],
            "EFFICIENCY": [0.1, 0.2, 0.3, 0.4],
        }
    )


def test_add_snr_efficiency_qcs_pins_full_row_contract_under_frozen_clock(
    log: object, frozen_clock: object
) -> None:
    result = toolkit.add_snr_efficiency_qcs(
        log,
        _snr_efficiency_spectrum(),
        pd.DataFrame(),
        {"1011": 500.0},
        "soxs-stare",
        "2024-01-02T03:04:05",
    )

    assert list(result.columns) == _EXPECTED_SNR_EFFICIENCY_QC_COLUMNS
    assert result["soxspipe_recipe"].tolist() == ["soxs-stare"] * 6
    assert result["qc_name"].tolist() == [
        "EFF MEDIAN",
        "EFF MEDIAN",
        "EFF MEDIAN",
        "SNR MEDIAN",
        "SNR MEDIAN",
        "SNR MEDIAN",
    ]
    assert result["qc_value"].tolist() == [0.25, 0.35, 0.15, 25.0, 35.0, 15.0]
    assert result["qc_comment"].tolist() == [
        "Median efficiency across all orders",
        "Median efficiency in order 10.0",
        "Median efficiency in order 11.0",
        "Median SNR across all orders",
        "Median SNR in order 10.0",
        "Median SNR in order 11.0",
    ]
    assert result["qc_unit"].isna().all()
    assert result["obs_date_utc"].tolist() == ["2024-01-02T03:04:05"] * 6
    assert result["reduction_date_utc"].tolist() == [_FROZEN_TIMESTAMP] * 6
    assert result["to_header"].tolist() == [True] * 6
    assert result["qc_order"].iloc[[0, 3]].isna().all()
    assert result["qc_order"].iloc[[1, 2, 4, 5]].tolist() == [10.0, 11.0, 10.0, 11.0]


def test_add_snr_efficiency_qcs_shares_single_timestamp_across_all_rows_under_ticking_clock(
    log: object, ticking_clock: object
) -> None:
    result = toolkit.add_snr_efficiency_qcs(
        log,
        _snr_efficiency_spectrum(),
        pd.DataFrame(),
        {"1011": 500.0},
        "soxs-stare",
        "2024-01-02T03:04:05",
    )

    timestamps = set(result["reduction_date_utc"])
    assert len(timestamps) == 1
    (timestamp,) = timestamps
    assert _TIMESTAMP_PATTERN.match(timestamp)


# ---------------------------------------------------------------------------
# plot_merged_spectrum_qc
# ---------------------------------------------------------------------------


def _merged_orders() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "WAVE": [500.0, 501.0, 502.0, 503.0, 504.0],
            "FLUX_COUNTS": [10.0, 12.0, 11.0, 13.0, 12.0],
            "SNR": [5.0, 6.0, 5.5, 7.0, 6.5],
            "SKY_COUNTS": [11.0, 12.0, 13.0, 12.0, 11.0],
        }
    )


def _qc_table_for_plot() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "qc_name": ["SNR MEDIAN", "SNR MEDIAN", "SNR MEDIAN"],
            "qc_order": ["GLOBAL", "12", "11"],
            "qc_value": [6.0, 5.0, 7.0],
        }
    )


def _stub_skylines(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda *_: pd.DataFrame({"WAVELENGTH": [501.0, 503.0], "ISOLATED": [True, False]}),
    )


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_plot_merged_spectrum_qc_saves_pdf_with_pinned_savefig_arguments(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    fluxCalibrated: bool,
) -> None:
    _stub_skylines(monkeypatch)
    recorder = _record_calls(monkeypatch, ["savefig"])

    products, outputPath = toolkit.plot_merged_spectrum_qc(
        _merged_orders(),
        product_table().iloc[0:0].copy(),
        log,
        str(tmp_path),
        "synthetic.fits",
        "_AB",
        "2024-01-02T03:04:05",
        "VIS",
        "soxs-nod",
        orderJoins={"11-12": 502.0},
        fluxCalibrated=fluxCalibrated,
        qcTable=_qc_table_for_plot(),
        settings={},
    )

    assert len(recorder["savefig"]) == 1
    savefigArgs, savefigKwargs = recorder["savefig"][0]
    assert savefigArgs == (outputPath,)
    assert savefigKwargs == {"dpi": 120, "bbox_inches": "tight", "format": "pdf"}
    assert Path(outputPath).read_bytes().startswith(b"%PDF")
    assert products is not None


@pytest.mark.parametrize(
    ("fluxCalibrated", "expectedLabel"),
    [
        (False, "EXTRACTED_MERGED_QC_PLOT_AB"),
        (True, "EXTRACTED_MERGED_FLUXCALIBRATED_QC_PLOT_AB"),
    ],
)
def test_plot_merged_spectrum_qc_pins_product_row_under_frozen_clock(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    frozen_clock: object,
    fluxCalibrated: bool,
    expectedLabel: str,
) -> None:
    _stub_skylines(monkeypatch)

    products, outputPath = toolkit.plot_merged_spectrum_qc(
        _merged_orders(),
        product_table().iloc[0:0].copy(),
        log,
        str(tmp_path),
        "synthetic.fits",
        "_AB",
        "2024-01-02T03:04:05",
        "VIS",
        "soxs-nod",
        orderJoins={"11-12": 502.0},
        fluxCalibrated=fluxCalibrated,
        qcTable=_qc_table_for_plot(),
        settings={},
    )

    assert list(products.columns) == [
        "soxspipe_recipe",
        "product_label",
        "file_name",
        "file_type",
        "obs_date_utc",
        "reduction_date_utc",
        "product_desc",
        "file_path",
        "label",
    ]
    row = products.iloc[-1]
    expectedFileName = f"synthetic_{expectedLabel}.pdf"
    assert row["soxspipe_recipe"] == "soxs-nod"
    assert row["product_label"] == expectedLabel
    assert row["file_name"] == expectedFileName
    assert row["file_type"] == "PDF"
    assert row["obs_date_utc"] == "2024-01-02T03:04:05"
    assert row["reduction_date_utc"] == _FROZEN_TIMESTAMP
    assert row["product_desc"] == "QC plot of extracted order-merged source"
    assert row["file_path"] == outputPath
    assert row["label"] == "QC"


def test_plot_merged_spectrum_qc_closes_all_figures_after_saving(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda *_: pd.DataFrame({"WAVELENGTH": [501.0], "ISOLATED": [True]}),
    )
    recorder = _record_calls(monkeypatch, ["savefig", "close"])

    toolkit.plot_merged_spectrum_qc(
        _merged_orders(),
        product_table().iloc[0:0].copy(),
        log,
        str(tmp_path),
        "synthetic.fits",
        False,
        "2024-01-02T03:04:05",
        "VIS",
        "soxs-nod",
        settings={},
    )

    assert recorder["__order__"] == ["savefig", "close"]
    closeArgs, closeKwargs = recorder["close"][0]
    assert closeArgs == ("all",)
    assert closeKwargs == {}


# ---------------------------------------------------------------------------
# quicklook_image
# ---------------------------------------------------------------------------


def _quicklook_frame(instrument: str = "SOXS") -> CCDData:
    data = np.arange(64, dtype=float).reshape(8, 8)
    return CCDData(
        data,
        unit=u.electron,
        mask=data % 7 == 0,
        uncertainty=StdDevUncertainty(np.ones((8, 8))),
        meta={"INSTRUME": instrument},
    )


def test_quicklook_image_saves_pdf_with_pinned_savefig_arguments_then_clears_the_figure(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    outputPath = tmp_path / "quicklook.pdf"
    recorder = _record_calls(monkeypatch, ["savefig", "clf"])

    result = toolkit.quicklook_image(
        log,
        _quicklook_frame(),
        show=False,
        saveToPath=outputPath,
    )

    assert result is None
    assert recorder["__order__"] == ["savefig", "clf"]
    savefigArgs, savefigKwargs = recorder["savefig"][0]
    assert savefigArgs == (outputPath,)
    assert savefigKwargs == {"dpi": 120, "format": "pdf", "bbox_inches": "tight"}
    assert Path(outputPath).read_bytes().startswith(b"%PDF")
