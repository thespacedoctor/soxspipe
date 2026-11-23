"""Object-trace sampling counts only slices with a finite Gaussian width as usable (DY-2281).

A slice keeps a NaN width when it has no data, no peak, or a rejected fit. Counted as usable, those slices hide an order
with almost no real fits from the longer-slice retry. The tests drive ``sample_trace`` through the SOXS NIR path, where
all orders share one ``find_centre_points`` call, and stub the slice fitter.
"""

from __future__ import annotations

import importlib
import warnings

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.commonutils.detect_continuum import detect_continuum

pytestmark = pytest.mark.unit
continuumModule = importlib.import_module("soxspipe.commonutils.detect_continuum")

ROWS_PER_ORDER = 60
SOLUTION_CENTRE = 16.0
SLICE_LENGTH = 12
ORDERS = [10.0, 11.0]
MEDIAN_SHIFT = 1.9847457627118636
MEDIAN_STDDEV = 1.4847457627118643


def _shift_for(yValues: np.ndarray) -> np.ndarray:
    """The object sits 1.9 to 2.1 pixels below the dispersion-solution centre."""
    return 1.9 + 0.2 * (yValues % ROWS_PER_ORDER) / (ROWS_PER_ORDER - 1)


def _stddev_for(yValues: np.ndarray) -> np.ndarray:
    """The trace width varies smoothly between 1.4 and 1.6 pixels."""
    return 1.4 + 0.2 * ((7 * yValues) % ROWS_PER_ORDER) / (ROWS_PER_ORDER - 1)


def _trace_sampler(log: object, monkeypatch: pytest.MonkeyPatch) -> detect_continuum:
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.inst = "SOXS"
    detector.arm = "NIR"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.recipeName = "soxs-stare"
    detector.recipeSettings = {"slice-length": SLICE_LENGTH, "peak-sigma-limit": 3.0, "slice-width": 5}
    detector.traceFrame = CCDData(np.ones((32, 32)), unit=u.electron)
    detector.kw = lambda key: {"DPR_TYPE": "DPR TYPE", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    detector.debug = False
    detector.qc = pd.DataFrame()
    detector.dateObs = "2024-01-02T03:04:05"
    rowCount = ROWS_PER_ORDER * len(ORDERS)
    samples = pd.DataFrame(
        {
            "order": np.repeat(ORDERS, ROWS_PER_ORDER),
            "wavelength": np.linspace(1000.0, 2000.0, rowCount),
            "fit_x": np.full(rowCount, SOLUTION_CENTRE),
            "fit_y": np.arange(rowCount, dtype=float),
        }
    )
    monkeypatch.setattr(detector, "create_pixel_arrays", lambda: (samples.copy(), 1, 1))
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None)
    return detector


def _found(orderPixelTable: pd.DataFrame) -> pd.DataFrame:
    """The slice-fit result for a trace that is found at every sample, as ``fit_1d_gaussian_to_slices`` builds it."""
    yValues = orderPixelTable["fit_y"].to_numpy()
    return orderPixelTable.assign(
        cont_x=SOLUTION_CENTRE - _shift_for(yValues),
        cont_y=yValues,
        gauss_stddev=_stddev_for(yValues),
    )


def _lose(table: pd.DataFrame, rows: np.ndarray, width: float = np.nan) -> pd.DataFrame:
    """Give ``rows`` the result of a failed slice fit: no position, and ``width`` as the Gaussian width."""
    return table.assign(
        cont_x=table["cont_x"].where(~rows),
        cont_y=table["cont_y"].where(~rows),
        gauss_stddev=table["gauss_stddev"].where(~rows, width),
    )


def test_an_order_with_no_finite_widths_triggers_the_longer_slice_retry(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: ORDER 11 FINDS NOTHING WITH THE DEFAULT SLICE, BUT THE LONGER RETRY SLICE FINDS IT
    detector = _trace_sampler(log, monkeypatch)
    sliceLengths = []

    def fit_slices(orderPixelTable: pd.DataFrame, sliceLength: float, **_: object) -> pd.DataFrame:
        sliceLengths.append(sliceLength)
        found = _found(orderPixelTable)
        if sliceLength > SLICE_LENGTH:
            return found
        return _lose(found, (orderPixelTable["order"] == 11.0).to_numpy())

    monkeypatch.setattr(detector, "fit_1d_gaussian_to_slices", fit_slices)

    # ACT
    result, detectionPercentage = detector.sample_trace()

    # ASSERT: THE FIRST PROBE RETURNS A NaN SHIFT, SO THE SECOND PROBE USES THE 1.2 TIMES LONGER SLICE
    assert sliceLengths[0] == SLICE_LENGTH
    assert sliceLengths[1] == pytest.approx(SLICE_LENGTH * 1.2)


LOST_ROWS = {20: np.nan, 75: np.nan, 80: np.inf, 33: -np.inf}


def test_rows_with_non_finite_widths_come_back_pre_clipped_and_reach_no_sigma_clip(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: FOUR SLICES FAIL, TWO WITH A NaN WIDTH AND TWO WITH AN INFINITE ONE; EVERY ORDER KEEPS 5 PROBE FITS
    detector = _trace_sampler(log, monkeypatch)

    def fit_slices(orderPixelTable: pd.DataFrame, **_: object) -> pd.DataFrame:
        found = _found(orderPixelTable)
        for row, width in LOST_ROWS.items():
            found = _lose(found, (orderPixelTable["fit_y"] == row).to_numpy(), width)
        return found

    monkeypatch.setattr(detector, "fit_1d_gaussian_to_slices", fit_slices)

    # ACT: THE DY-2409 FILTER IS OVERRIDDEN, SO A sigma_clip "invalid values" WARNING FAILS THE TEST
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", category=DeprecationWarning, message=r"datetime\.datetime\.utcnow")
        result, detectionPercentage = detector.sample_trace()

    # ASSERT
    lost = result["fit_y"].isin(list(LOST_ROWS)).to_numpy()
    assert result.loc[lost, "pre-clipped"].eq(True).all()
    assert result.loc[~lost, "pre-clipped"].eq(False).all()
    assert ("debug", "find_centre_points: 4 slices have no finite Gaussian width") in log.messages


def test_a_table_with_only_finite_widths_keeps_its_mask_shift_and_width(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: NO WIDTH IS NaN; TWO NARROW (UNDER 0.3 PIXEL) FITS ARE PRE-CLIPPED BY THE WIDTH RULE
    detector = _trace_sampler(log, monkeypatch)
    calls = []

    def fit_slices(orderPixelTable: pd.DataFrame, medianStddev: object = False, **_: object) -> pd.DataFrame:
        calls.append((orderPixelTable["fit_x"].to_numpy().copy(), medianStddev))
        found = _found(orderPixelTable)
        narrow = orderPixelTable["fit_y"].isin([7, 64]).to_numpy()
        return found.assign(gauss_stddev=found["gauss_stddev"].where(~narrow, 0.2))

    monkeypatch.setattr(detector, "fit_1d_gaussian_to_slices", fit_slices)

    # ACT
    result, detectionPercentage = detector.sample_trace()

    # ASSERT: THE VALUES DEVELOP 2c169e33 GIVES FOR THIS TABLE
    fullFitX, fullFitWidth = calls[-1]
    assert len(calls) == 2
    assert SOLUTION_CENTRE - fullFitX[0] == pytest.approx(MEDIAN_SHIFT, abs=1e-6)
    assert fullFitWidth == pytest.approx(MEDIAN_STDDEV, abs=1e-6)
    assert sorted(result.loc[result["pre-clipped"], "fit_y"]) == [7.0, 64.0]
