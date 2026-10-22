"""Acceptance checks for the representative real NIR offset reduction.

The reduction is not bit-reproducible across hardware, so each value is asserted
within a band rather than exactly. The bands absorb the variation seen across CI
runner architectures and stay narrow enough that a real regression in extraction,
wavelength calibration or flux calibration still trips them. The centres were
re-recorded on 2026-10-16 from two bit-identical CI runs (38915352201 and
38915356992), after DY-1612 made the set (two ON positions) align its ON frames
along the slit and stack them: one ONOFF and one OBJTRACE product instead of one
per cycle, and the SNR median rose from 104.0 to 149.8. The flux-calibrated median
was re-recorded on 2026-10-18 from six bit-identical CI runs (39165616944,
39165618635, 39165620359, 39165622153, 39165623841 and 39165625784), after the
response curve became a log-space spline (DY-1614). It fell 3.0%. It was
re-recorded on 2026-10-19 from six bit-identical CI runs (39290451388, 39290461991,
39290472137, 39290482158, 39290492461 and 39290503644), after the NIR response fit
change (DY-1647), and rose 0.6% (DY-1648). It was re-recorded on 2026-10-21 from six
bit-identical CI runs (39680454564, 39680459968, 39680464192, 39680468399, 39680472907
and 39680476407), after the nod shift-and-add (DY-1674) and the flux-conserving slit
shift (DY-1679), and rose 0.15%. The SNR median did not move (DY-1691).
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
        f"{_PRODUCT_STEM}_OBJTRACE_STACK.fits",
        f"{_PRODUCT_STEM}_ONOFF_STACK.fits",
    }
    productNames = {path.name for path in product_directory.iterdir()}

    merged_path = product_directory / f"{_PRODUCT_STEM}_EXTRACTED_MERGED.fits"
    fluxcal_path = product_directory / f"{_PRODUCT_STEM}_FLUXCAL.fits"
    with fits.open(merged_path) as merged_hdus:
        assert len(merged_hdus) == 2
        assert merged_hdus[0].header["INSTRUME"] == "SOXS"
        assert merged_hdus[0].header["ESO PRO REC1 ID"] == "soxs-offset"
        dateObs = merged_hdus[0].header["DATE-OBS"]
        merged_table = merged_hdus[1].data
        assert set(merged_table.names) == {
            "WAVE",
            "FLUX_COUNTS",
            "VARIANCE",
            "SKY_COUNTS",
            "SNR",
        }
        mergedRows = len(merged_table)
        mergedWave = np.asarray(merged_table["WAVE"], dtype=float)
        snrMedian = float(np.nanmedian(merged_table["SNR"]))
    waveMin = float(np.nanmin(mergedWave))
    waveMax = float(np.nanmax(mergedWave))
    waveStep = float(np.nanmedian(np.diff(mergedWave)))

    with fits.open(fluxcal_path) as fluxcal_hdus:
        assert len(fluxcal_hdus) == 2
        fluxcal_table = fluxcal_hdus[1].data
        assert set(fluxcal_table.names) == {"WAVE", "FLUX_CALIBRATED"}
        fluxcalRows = len(fluxcal_table)
        fluxcalMedian = float(np.nanmedian(fluxcal_table["FLUX_CALIBRATED"]))

    with sqlite3.connect(reduced_workspace / "soxspipe.db") as connection:
        qc_values = dict(
            connection.execute(
                "SELECT qc_name, qc_value FROM quality_control "
                "WHERE soxspipe_recipe = ? AND obs_date_utc LIKE ? "
                "AND qc_name IN (?, ?)",
                ("soxs-offset", "2025-05-14%", "N ORDERS", "SAMPLES DET FRAC"),
            )
        )
    nOrders = float(qc_values["N ORDERS"])
    detectedFraction = float(qc_values["SAMPLES DET FRAC"])

    # EVERY VALUE IS REPORTED BEFORE THE FIRST BAND ASSERTION, SO A FAILING RUN STILL
    # RECORDS ALL OF THEM FOR THE NEXT RE-RECORD
    report("merged rows", mergedRows)
    report("wave min", waveMin)
    report("wave max", waveMax)
    report("wave step median", waveStep)
    report("snr median", snrMedian)
    report("fluxcal rows", fluxcalRows)
    report("fluxcal median", fluxcalMedian)
    report("n orders", nOrders)
    report("samples det frac", detectedFraction)
    report("products", sorted(productNames))
    report("date-obs", dateObs)

    assert productNames == expected_products
    assert dateObs == "2025-05-14T04:03:11.8311"

    # THE MERGED GRID RUNS FROM THE BLUEST TO THE REDDEST EXTRACTED SAMPLE IN
    # 0.06 NM STEPS, SO THE ROW COUNT FOLLOWS THE RED END. ±25 ROWS IS ±1.5 NM
    # OF SPECTRAL COVERAGE
    assert mergedRows == pytest.approx(20_614, abs=25)
    # THE BLUE END IS PINNED BY THE ORDER LAYOUT, SO IT IS HELD TO TWO GRID STEPS
    assert waveMin == pytest.approx(795.06, abs=0.12)
    assert waveMax == pytest.approx(2031.84, abs=1.0)
    assert waveStep == pytest.approx(0.06, abs=1e-6)
    assert snrMedian == pytest.approx(149.8, abs=2.0)
    assert fluxcalRows == pytest.approx(20_614, abs=25)
    assert fluxcalRows == mergedRows
    assert fluxcalMedian == pytest.approx(8.243045025961036e-15, rel=0.05)
    # EVERY ORDER MUST BE PRESENT: A MISSING ORDER IS A FAILURE, NOT DRIFT
    assert nOrders == 15
    assert detectedFraction == pytest.approx(0.984, abs=0.02)
