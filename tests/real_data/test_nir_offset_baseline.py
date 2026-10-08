"""Acceptance checks for the representative real NIR offset reduction.

The reduction is not bit-reproducible across hardware, so each value is asserted
within a band rather than exactly. The bands absorb the variation seen across CI
runner architectures and stay narrow enough that a real regression in extraction,
wavelength calibration or flux calibration still trips them. The centres were
re-recorded from bit-identical CI runs on 2026-10-07, after the astropy 7.2.2 upgrade
(DY-694). astropy 7 stopped running `sigma_clip` through bottleneck for float32 data,
which had lost precision, so every product moved.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from tests.real_data.reporting import report

pytestmark = pytest.mark.real_data


_PRODUCT_STEM = "20250514T040311_NIR_3_OFFSET_OBJ_SLIT1_0_300_0S_SOXS_CD-4412736"


def _offset_product_directory(workspace_path: Path) -> Path:
    product_paths = list(workspace_path.rglob(f"{_PRODUCT_STEM}_EXTRACTED_MERGED.fits"))
    assert len(product_paths) == 1
    return product_paths[0].parent


def test_nir_offset_reduction_matches_approved_baseline(reduced_workspace: Path) -> None:
    product_directory = _offset_product_directory(reduced_workspace)
    expected_products = {
        f"{_PRODUCT_STEM}.log",
        f"{_PRODUCT_STEM}_EXTRACTED_MERGED.fits",
        f"{_PRODUCT_STEM}_EXTRACTED_MERGED.txt",
        f"{_PRODUCT_STEM}_FLUXCAL.fits",
        f"{_PRODUCT_STEM}_OBJTRACE_1.fits",
        f"{_PRODUCT_STEM}_OBJTRACE_2.fits",
        f"{_PRODUCT_STEM}_ONOFF_1.fits",
        f"{_PRODUCT_STEM}_ONOFF_2.fits",
    }
    assert {path.name for path in product_directory.iterdir()} == expected_products

    merged_path = product_directory / f"{_PRODUCT_STEM}_EXTRACTED_MERGED.fits"
    fluxcal_path = product_directory / f"{_PRODUCT_STEM}_FLUXCAL.fits"
    with fits.open(merged_path) as merged_hdus:
        assert len(merged_hdus) == 2
        assert merged_hdus[0].header["INSTRUME"] == "SOXS"
        assert merged_hdus[0].header["ESO PRO REC1 ID"] == "soxs-offset"
        assert merged_hdus[0].header["DATE-OBS"] == "2025-05-14T04:03:11.8311"
        merged_table = merged_hdus[1].data
        # THE MERGED GRID RUNS FROM THE BLUEST TO THE REDDEST EXTRACTED SAMPLE IN
        # 0.06 NM STEPS, SO THE ROW COUNT FOLLOWS THE RED END. ±25 ROWS IS ±1.5 NM
        # OF SPECTRAL COVERAGE
        report("merged rows", len(merged_table))
        assert len(merged_table) == pytest.approx(20_608, abs=25)
        assert set(merged_table.names) == {
            "WAVE",
            "FLUX_COUNTS",
            "VARIANCE",
            "SKY_COUNTS",
            "SNR",
            "FLUX_DENSITY_COUNTS",
        }
        mergedWave = np.asarray(merged_table["WAVE"], dtype=float)
        # THE BLUE END IS PINNED BY THE ORDER LAYOUT, SO IT IS HELD TO ONE GRID STEP
        report("wave min", float(np.nanmin(mergedWave)))
        report("wave max", float(np.nanmax(mergedWave)))
        report("wave step median", float(np.nanmedian(np.diff(mergedWave))))
        report("snr median", float(np.nanmedian(merged_table["SNR"])))
        assert float(np.nanmin(mergedWave)) == pytest.approx(795.06, abs=0.12)
        assert float(np.nanmax(mergedWave)) == pytest.approx(2031.48, abs=1.0)
        assert float(np.nanmedian(np.diff(mergedWave))) == pytest.approx(0.06, abs=1e-6)
        assert float(np.nanmedian(merged_table["SNR"])) == pytest.approx(103.505, abs=2.0)

    with fits.open(fluxcal_path) as fluxcal_hdus:
        assert len(fluxcal_hdus) == 2
        fluxcal_table = fluxcal_hdus[1].data
        assert len(fluxcal_table) == pytest.approx(20_608, abs=25)
        assert len(fluxcal_table) == len(merged_table)
        assert set(fluxcal_table.names) == {"WAVE", "FLUX_CALIBRATED"}
        report("fluxcal rows", len(fluxcal_table))
        report("fluxcal median", float(np.nanmedian(fluxcal_table["FLUX_CALIBRATED"])))
        assert float(np.nanmedian(fluxcal_table["FLUX_CALIBRATED"])) == pytest.approx(8.425130566412588e-15, rel=0.05)

    with sqlite3.connect(reduced_workspace / "soxspipe.db") as connection:
        qc_values = dict(
            connection.execute(
                "SELECT qc_name, qc_value FROM quality_control "
                "WHERE soxspipe_recipe = ? AND obs_date_utc LIKE ? "
                "AND qc_name IN (?, ?)",
                ("soxs-offset", "2025-05-14%", "N ORDERS", "SAMPLES DET FRAC"),
            )
        )
    # EVERY ORDER MUST BE PRESENT: A MISSING ORDER IS A FAILURE, NOT DRIFT
    report("n orders", float(qc_values["N ORDERS"]))
    report("samples det frac", float(qc_values["SAMPLES DET FRAC"]))
    assert float(qc_values["N ORDERS"]) == 15
    assert float(qc_values["SAMPLES DET FRAC"]) == pytest.approx(0.975, abs=0.02)
