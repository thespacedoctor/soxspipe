"""A stare extraction with no sky model keeps NaN in SKY_COUNTS from the extraction to the ASCII file (DY-1361)."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.horne_extraction import extract_single_order
from soxspipe.commonutils.spectrum_ascii import write_ascii_spectrum_product
from tests.factories import pipeline_settings, product_table, response_product_fits, soxs_merge_extractor
from tests.factories.flux_calibration import bare_flux_calibration

pytestmark = pytest.mark.integration

KEYWORDS = {
    "DPR_CATG": "ESO DPR CATG",
    "DPR_TYPE": "ESO DPR TYPE",
    "DET_READ_SPEED": "ESO DET READ CURID",
    "CONAD": "ESO DET OUT1 CONAD",
    "GAIN": "ESO DET OUT1 GAIN",
    "RON": "ESO DET OUT1 RON",
    "SEQ_ARM": "ESO SEQ ARM",
    "PRO_TYPE": "ESO PRO TYPE",
    "PRO_CATG": "ESO PRO CATG",
    "PRO_AIRM_MEAN": "ESO PRO AIRM MEAN",
    "WIN_BINY": "ESO DET BINY",
}


def _orders_without_a_sky_model() -> pd.DataFrame:
    """Extract one real order from rectified images that hold no sky frame, so ``skyFlux`` is NaN."""
    wavelength = np.array([[504.0, 503.0, 502.0, 501.0, 500.0]] * 2)
    images = {
        "bpMask": np.zeros((2, 5), dtype=bool),
        "fluxRaw": np.array([[1.0, 2.0, 3.0, 4.0, 5.0]] * 2),
        "variance": np.ones((2, 5)),
        "wavelength": wavelength,
    }
    orders = extract_single_order(
        (pd.DataFrame({"order": [10] * 5}), images),
        logging.getLogger("dy1361"),
        slitHalfLength=1,
        clippingSigma=3.0,
        clippingIterationLimit=2,
        axisA="y",
        axisB="x",
        hornePolyOrder=0,
        boxcarOnly=True,
    )
    assert orders["skyFlux"].isna().all()
    return orders


def _merge(log: object, orders: pd.DataFrame) -> pd.DataFrame:
    """Order-merge ``orders`` with the real merge, as the stare extraction does."""
    merged, _ = soxs_merge_extractor(log, "NIR").merge_extracted_orders(orders)
    return merged


def _calibrate(tmp_path: Path, log: object, merged: pd.DataFrame) -> str:
    """Flux calibrate ``merged`` for real, writing ``<stem>_SPECTRUM.fits``, and return its path."""
    calibration = bare_flux_calibration(
        log=log,
        settings=pipeline_settings(tmp_path),
        responseFunction=str(response_product_fits(tmp_path / "RESP.fits")),
        extractedSpectrum=merged,
        exptime=10.0,
        airmass=1.0,
        extinctionPath="",
        arm="NIR",
        header=fits.Header({"DATE-OBS": "2024-01-02T03:04:05"}),
        recipeName="soxs-stare",
        sofName="SYNTHETIC",
        productDir=str(tmp_path),
        kw=KEYWORDS.__getitem__,
    )
    filePath, _ = calibration.calibrate()
    return filePath


def test_a_stare_extraction_with_no_sky_model_keeps_nan_sky_counts_in_the_spectrum_table_and_the_ascii_file(
    tmp_path: Path, log: object
) -> None:
    # ARRANGE: A REAL EXTRACTION AND MERGE OF AN ORDER WITH NO SKY FRAME
    merged = _merge(log, _orders_without_a_sky_model())

    # ACT
    fluxcalPath = _calibrate(tmp_path, log, merged)
    products = write_ascii_spectrum_product(fluxcalPath, True, product_table().iloc[0:0], "soxs-stare", "2024")

    # ASSERT: THE SKY IS NAN, NOT 0, IN THE FLUX-CALIBRATED TABLE AND IN THE ASCII FILE, AND THE FLUX IS NOT AFFECTED
    fluxcal = Table.read(fluxcalPath, hdu=1)
    assert len(fluxcal) > 0
    assert np.isnan(np.asarray(fluxcal["SKY_COUNTS"], dtype=float)).all()
    assert np.isfinite(np.asarray(fluxcal["FLUX_CALIBRATED"], dtype=float)).all()
    lines = Path(products.loc[0, "file_path"]).read_text().splitlines()
    skyColumn = lines[0].split().index("SKY_COUNTS")
    assert {line.split()[skyColumn].lower() for line in lines[1:]} == {"nan"}
