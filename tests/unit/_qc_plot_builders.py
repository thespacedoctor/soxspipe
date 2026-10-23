"""Builders shared by the QC plot tests: synthetic frames, tables and plot-owning objects.

The sky-subtraction, merged-spectrum, quicklook and continuum plot tests each need a plot-owning object with just the
state its plot reads. They live here so that the characterization tests and the figure-cleanup test use one copy.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from scipy.interpolate import splrep

import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header, pipeline_settings, product_table, qc_table

FLAG_COLUMNS = (
    "flagged_all_clipped",
    "flagged_object_clipped",
    "flagged_bad_pixel_clipped",
    "flagged_edge_clipped",
    "flagged_bspline_clipped",
)

ORDER_PIXEL_COUNT = 200

SOF_NAME = "SOXS_VIS_FLAT_SOF"

QC_DIR_NAME = "qc"

PRODUCT_DIR_NAME = "products"


def quicklook_ccd(
    instrument: str = "SOXS", shape: tuple[int, int] = (8, 8), extraHeader: dict | None = None
) -> CCDData:
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


def synthetic_merged_orders(*, includeSky: bool = True) -> pd.DataFrame:
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


def stub_merged_skylines(monkeypatch: pytest.MonkeyPatch, extra: tuple = ()) -> None:
    wavelengths = [520.0, 540.0, *extra]
    isolated = [True, False] + [False] * len(extra)
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame({"WAVELENGTH": wavelengths, "ISOLATED": isolated}),
    )


def call_plot_merged_spectrum_qc(
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


def strip_subtractor(
    log: Any,
    outputPath: Path,
    *,
    dispersionAxis: str,
    rotate: int | bool,
    flip: int | bool = False,
    shape: tuple[int, int] = (16, 16),
) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    # THE CONSTRUCTOR SWAPS THE AXES FOR A Y-DISPERSION ARM
    subtractor.axisA, subtractor.axisB = ("x", "y") if dispersionAxis == "x" else ("y", "x")
    subtractor.filenameTemplate = "science.fits"
    subtractor.qcDir = str(outputPath)
    subtractor.detectorParams = {
        "dispersion-axis": dispersionAxis,
        "flip-qc-plot": flip,
        "rotate-qc-plot": rotate,
    }
    subtractor.objectFrame = CCDData(
        np.arange(np.prod(shape), dtype=float).reshape(shape),
        unit="electron",
        mask=np.zeros(shape, dtype=bool),
    )
    return subtractor


def order_strip(
    *, pixelCount: int = 64, columns: int = 4, flagged: bool = True, largeScatter: bool = False
) -> pd.DataFrame:
    """An order `columns` pixels wide running down the rows, with every ninth pixel flagged."""
    index = np.arange(pixelCount)
    pixels = pd.DataFrame(
        {
            "x": index % columns,
            "y": index // columns,
            "wavelength": np.linspace(500.0, 510.0, pixelCount),
            "flux": 12.0 + np.sin(index),
            "slit_position": np.linspace(-1.0, 1.0, pixelCount),
            "residual_windowed_std": np.full(pixelCount, 3000.0 if largeScatter else 1.0),
            "residual_global_sigma": np.cos(index),
            "flux_percentile_smoothed": np.full(pixelCount, 12.0),
            "sky_model": 100.0 + index,
            "sky_subtracted_flux": np.sin(index / 3.0),
            "sky_subtracted_flux_weighted": np.sin(index / 3.0),
        }
    )
    for column in FLAG_COLUMNS:
        pixels[column] = (index % 9 == 0) if flagged else False
    return pixels


def quicklook_subtractor(log: Any) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.settings = {}
    subtractor.arm = "VIS"
    return subtractor


def quicklook_skylines(*, empty: bool = False) -> pd.DataFrame:
    if empty:
        return pd.DataFrame({"WAVELENGTH": [], "FLUX": [], "ISOLATED": []})
    return pd.DataFrame(
        {
            "WAVELENGTH": [502.0, 505.0, 508.0],
            "FLUX": [10.0, 20.0, 30.0],
            "ISOLATED": [True, False, True],
        }
    )


def quicklook_order_pixels() -> pd.DataFrame:
    index = np.arange(ORDER_PIXEL_COUNT)
    flux = 100.0 + 20.0 * np.sin(index / 7.0)
    flux[::37] = -80.0
    skyLine = np.where(index % 11 == 0, "line", np.where(index % 13 == 0, "node", "")).astype(object)
    skyLine[skyLine == ""] = False
    return pd.DataFrame(
        {
            "wavelength": np.linspace(500.0, 510.0, ORDER_PIXEL_COUNT),
            "flux": flux,
            "flagged_object_clipped": index % 41 == 0,
            "flagged_all_clipped": index % 29 == 0,
            "flagged_bad_pixel_clipped": index % 53 == 0,
            "flagged_edge_clipped": index < 3,
            "flux_percentile_smoothed": 100.0 + np.cos(index / 5.0),
            "flux_minus_smoothed_residual": 4.0 * np.sin(index / 3.0),
            "flux_minus_smoothed_residual_upper_limit": np.full(ORDER_PIXEL_COUNT, 12.0),
            "flux_upper_limit": np.full(ORDER_PIXEL_COUNT, 112.0),
            "flagged_noisy_region": index % 17 == 0,
            "flagged_sky_line": skyLine,
            "flagged_bspline_clipped": index % 19 == 0,
            "sky_residuals": 3.0 * np.abs(np.cos(index / 4.0)),
            "sky_residual_floor": np.full(ORDER_PIXEL_COUNT, 1.5),
            "sky_residual_rolling_average": 2.0 * np.abs(np.sin(index / 9.0)),
        }
    )


def quicklook_spline() -> tuple:
    index = np.arange(ORDER_PIXEL_COUNT)
    wavelength = np.linspace(500.0, 510.0, ORDER_PIXEL_COUNT)
    return splrep(wavelength, 100.0 + 20.0 * np.sin(index / 7.0), k=3, s=ORDER_PIXEL_COUNT * 400)


def continuum_detector(
    log: object,
    tmp_path: Path,
    *,
    recipeName: str,
    noddingSequence: str = "",
    notFlattened: bool = False,
) -> detect_continuum:
    """Return a continuum detector holding the state the file-naming code reads."""
    settings = pipeline_settings(tmp_path)
    settings["tune-pipeline"] = False
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.settings = settings
    detector.kw = keyword_lookup(log=log, settings=settings).get
    detector.recipeName = recipeName
    detector.noddingSequence = noddingSequence
    detector.notFlattened = "_NOTFLAT" if notFlattened else ""
    detector.sofName = SOF_NAME
    detector.lampTag = False
    detector.inst = "SOXS"
    detector.arm = "VIS"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.orderDeg = 0
    detector.axisBDeg = 1
    detector.sliceLength = 10
    detector.exptime = 60.0
    detector.lamp = ""
    detector.debug = False
    detector.recipeSettings = {}
    detector.qc = qc_table()
    detector.dateObs = "2024-01-02T03:04:05.678"
    detector.products = pd.DataFrame()
    detector.qcDir = str(tmp_path / QC_DIR_NAME)
    detector.productDir = str(tmp_path / PRODUCT_DIR_NAME)
    Path(detector.qcDir).mkdir(exist_ok=True)
    detector.detectorParams = {"rotate-qc-plot": 0, "flip-qc-plot": 0, "dispersion-axis": "x"}
    detector.traceFrame = CCDData(
        np.full((32, 32), 20.0),
        unit=u.electron,
        meta=instrument_header(arm="VIS", overrides={"EXPTIME": 60.0}),
    )
    return detector


def continuum_plot_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return a one-order pixel table, polynomial table and clipped-sample table for ``plot_results``."""
    yValues = np.arange(3.0, 30.0, 3.0)
    pixels = pd.DataFrame(
        {
            "order": np.full(len(yValues), 10.0),
            "cont_x": np.full(len(yValues), 16.0),
            "cont_y": yValues,
            "cont_x_fit_res": np.zeros(len(yValues)),
            "wavelength": np.linspace(500.0, 510.0, len(yValues)),
            "gauss_stddev_fit": np.full(len(yValues), 2.0),
        }
    )
    polynomials = pd.DataFrame([{"cent_00": 16.0, "cent_01": 0.0, "std_00": 2.0, "std_01": 0.0}])
    clipped = pd.DataFrame({"cont_x": [17.0], "cont_y": [6.0], "fit_x": [16.0], "fit_y": [6.0]})
    return pixels, polynomials, clipped
