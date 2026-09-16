"""Acceptance checks for the representative real NIR offset reduction.

The reduction is not bit-reproducible across hardware, so each value is asserted
within a band rather than exactly. The bands absorb the variation seen across CI
runner architectures and stay narrow enough that a real regression in extraction,
wavelength calibration or flux calibration still trips them. The centres were
recorded from eight agreeing CI runs on 2026-09-16, after the deterministic-sort
and residual-quantisation fixes.
"""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

pytestmark = pytest.mark.real_data


_PRODUCT_STEM = "20250514T040311_NIR_3_OFFSET_OBJ_SLIT1_0_300_0S_SOXS_CD-4412736"


@pytest.fixture
def reduced_workspace() -> Path:
    """Return the workflow-created real-data workspace, or skip when absent."""
    workspace_value = os.environ.get("SOXSPIPE_REAL_DATA_DIR")
    if workspace_value is None:
        pytest.skip("SOXSPIPE_REAL_DATA_DIR is required for real-data tests")
    workspace_path = Path(workspace_value).resolve()
    if not workspace_path.is_dir():
        pytest.skip("SOXSPIPE_REAL_DATA_DIR does not name a workspace directory")
    return workspace_path


def _offset_product_directory(workspace_path: Path) -> Path:
    product_paths = list(workspace_path.rglob(f"{_PRODUCT_STEM}_EXTRACTED_MERGED.fits"))
    assert len(product_paths) == 1
    return product_paths[0].parent


def test_nir_offset_reduction_matches_approved_baseline(reduced_workspace: Path) -> None:
    product_directory = _offset_product_directory(reduced_workspace)
    expected_products = {
        f"{_PRODUCT_STEM}.fits",
        f"{_PRODUCT_STEM}.log",
        f"{_PRODUCT_STEM}_EXTRACTED_MERGED.fits",
        f"{_PRODUCT_STEM}_EXTRACTED_MERGED.txt",
        f"{_PRODUCT_STEM}_FLUXCAL.fits",
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
        assert len(merged_table) == pytest.approx(20_609, abs=25)
        assert set(merged_table.names) == {
            "WAVE", "FLUX_COUNTS", "VARIANCE", "SKY_COUNTS", "SNR", "FLUX_DENSITY_COUNTS"
        }
        mergedWave = np.asarray(merged_table["WAVE"], dtype=float)
        # THE BLUE END IS PINNED BY THE ORDER LAYOUT, SO IT IS HELD TO ONE GRID STEP
        assert float(np.nanmin(mergedWave)) == pytest.approx(794.76, abs=0.12)
        assert float(np.nanmax(mergedWave)) == pytest.approx(2031.24, abs=1.0)
        assert float(np.nanmedian(np.diff(mergedWave))) == pytest.approx(0.06, abs=1e-6)
        assert float(np.nanmedian(merged_table["SNR"])) == pytest.approx(104.56, abs=2.0)

    with fits.open(fluxcal_path) as fluxcal_hdus:
        assert len(fluxcal_hdus) == 2
        fluxcal_table = fluxcal_hdus[1].data
        assert len(fluxcal_table) == pytest.approx(20_609, abs=25)
        assert len(fluxcal_table) == len(merged_table)
        assert set(fluxcal_table.names) == {"WAVE", "FLUX_CALIBRATED"}
        assert float(np.nanmedian(fluxcal_table["FLUX_CALIBRATED"])) == pytest.approx(
            8.704858750481978e-15, rel=0.05
        )

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
    assert float(qc_values["N ORDERS"]) == 15
    assert float(qc_values["SAMPLES DET FRAC"]) == pytest.approx(0.97, abs=0.02)
