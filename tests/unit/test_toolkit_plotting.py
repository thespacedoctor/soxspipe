"""Behavioral contracts for toolkit plotting paths."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils import toolkit
from tests.factories import product_table

pytestmark = pytest.mark.unit


def _frame(instrument: str = "SOXS") -> CCDData:
    data = np.arange(64, dtype=float).reshape(8, 8)
    return CCDData(
        data,
        unit=u.electron,
        mask=data % 7 == 0,
        uncertainty=StdDevUncertainty(np.ones((8, 8))),
        meta={"INSTRUME": instrument},
    )


def test_quicklook_returns_before_importing_plotting_when_disabled(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    monkeypatch.setattr(plt, "figure", lambda *args, **kwargs: pytest.fail("plotted"))

    assert toolkit.quicklook_image(log, _frame(), show=False) is None


@pytest.mark.parametrize("extension", ["data", "mask", "uncertainty"])
def test_quicklook_renders_each_ccd_extension_without_showing(
    tmp_path: Path, log: object, extension: str
) -> None:
    outputPath = tmp_path / f"{extension}.pdf"

    result = toolkit.quicklook_image(
        log,
        _frame("XSHOOTER"),
        show=False,
        ext=extension,
        title=extension,
        saveToPath=outputPath,
    )

    assert result is None
    assert outputPath.read_bytes().startswith(b"%PDF")


@pytest.mark.parametrize("instrument", ["SOXS", "XSHOOTER", "OTHER"])
def test_quicklook_surface_plot_supports_each_orientation(
    tmp_path: Path, log: object, instrument: str
) -> None:
    outputPath = tmp_path / f"surface-{instrument}.pdf"

    toolkit.quicklook_image(
        log,
        np.arange(64, dtype=float).reshape(8, 8) + 20,
        show=False,
        inst=instrument,
        surfacePlot=True,
        saveToPath=outputPath,
    )

    assert outputPath.is_file()


def test_quicklook_uses_default_instrument_for_plain_array(
    tmp_path: Path, log: object
) -> None:
    outputPath = tmp_path / "array.pdf"

    toolkit.quicklook_image(
        log,
        np.zeros((8, 8)),
        show=False,
        saveToPath=outputPath,
    )

    assert outputPath.is_file()


def test_merged_spectrum_plot_skips_work_when_no_product_table_is_supplied(
    log: object,
) -> None:
    """QC plotting retains the no-product shortcut used by lightweight callers."""
    products, path = toolkit.plot_merged_spectrum_qc(
        pd.DataFrame(),
        False,
        log,
        "unused",
        "unused.fits",
        False,
        "2024-01-02T03:04:05",
        "VIS",
        "soxs-stare",
    )

    assert products is False
    assert path is None


@pytest.mark.parametrize("fluxCalibrated", [False, True])
def test_merged_spectrum_plot_writes_qc_pdf_and_product_record(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    fluxCalibrated: bool,
) -> None:
    """Merged spectra generate a local plot with sky lines, joins, and SNR labels."""
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda *_: pd.DataFrame(
            {"WAVELENGTH": [501.0, 503.0], "ISOLATED": [True, False]}
        ),
    )
    merged = pd.DataFrame(
        {
            "WAVE": [500.0, 501.0, 502.0, 503.0, 504.0],
            "FLUX_COUNTS": [10.0, 12.0, 11.0, 13.0, 12.0],
            "SNR": [5.0, 6.0, 5.5, 7.0, 6.5],
            "SKY_COUNTS": [11.0, 12.0, 13.0, 12.0, 11.0],
        }
    )
    qc = pd.DataFrame(
        {
            "qc_name": ["SNR MEDIAN", "SNR MEDIAN", "SNR MEDIAN"],
            "qc_order": ["GLOBAL", "12", "11"],
            "qc_value": [6.0, 5.0, 7.0],
        }
    )

    products, outputPath = toolkit.plot_merged_spectrum_qc(
        merged,
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
        qcTable=qc,
        settings={},
    )

    assert Path(outputPath).read_bytes().startswith(b"%PDF")
    product = products.iloc[-1]
    assert product["soxspipe_recipe"] == "soxs-nod"
    assert product["label"] == "QC"
    assert product["file_path"] == outputPath
    assert ("FLUXCALIBRATED" in product["product_label"]) is fluxCalibrated


@pytest.mark.parametrize("useSkylines", [False, True])
def test_dispersion_grid_builds_order_edges_and_cross_lines(
    monkeypatch: pytest.MonkeyPatch, log: object, useSkylines: bool
) -> None:
    mapFrame = pd.DataFrame(
        {
            "order": [10, 10, 10, 11, 11, 11],
            "wavelength": [500.0, 501.0, 502.0, 600.0, 601.0, 602.0],
            "slit_position": [-1.0, 0.0, 1.0, -1.0, 0.0, 1.0],
        }
    )
    interOrderMask = np.array([[False, True]])
    monkeypatch.setattr(
        toolkit,
        "twoD_disp_map_image_to_dataframe",
        lambda **kwargs: (mapFrame, interOrderMask),
    )
    monkeypatch.setattr(
        toolkit,
        "dispersion_map_to_pixel_arrays",
        lambda **kwargs: kwargs["orderPixelTable"].assign(
            fit_x=lambda table: table["wavelength"],
            fit_y=lambda table: table["slit_position"],
        ),
    )
    skylines = pd.DataFrame({"WAVELENGTH": [501.0, 601.0]}) if useSkylines else False

    result, mask = toolkit.create_dispersion_solution_grid_lines_for_plot(
        log,
        "solution.fits",
        "image.fits",
        _frame(),
        lambda name: name,
        skylines=skylines,
        slitPositions=(-0.5, 0.5),
    )

    assert mask is interOrderMask
    assert result["order"].unique().tolist() == [10.0, 11.0]
    assert {"line", "fit_x", "fit_y"}.issubset(result.columns)
    if useSkylines:
        assert {501.0, 601.0}.issubset(result["wavelength"])
