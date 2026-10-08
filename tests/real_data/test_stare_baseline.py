"""Acceptance checks for the real VIS and NIR stare reductions.

The stare recipe is the only recipe that models the sky with the b-spline fit in
`subtract_sky`, so these checks are the end-to-end guard on that fit. As with the offset
check, the reduction is not bit-reproducible across hardware, so each value is asserted
within a band. The sky bands are narrow enough that a sky model of zero, or one doubled,
fails them. The centres were re-recorded from four bit-identical CI runs on
2026-10-07, after the astropy 7.2.2 upgrade (DY-694). A local arm64 reduction lands up
to 5% away, so expect local runs to sit near the band edges. The NIR detected fraction is
the exception: a local arm64 run gives about 0.92 against the CI value of 0.786, so it
fails locally (DY-994).
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from tests.real_data.reporting import report

pytestmark = pytest.mark.real_data

# A SKY OR FLUX LEVEL MAY DRIFT BY THIS FRACTION ACROSS RUNNER HARDWARE. A BROKEN SKY
# MODEL (ZERO, OR DOUBLED) MOVES THESE LEVELS BY 100%, SO IT CANNOT HIDE INSIDE THE BAND
LEVEL_REL = 0.05
# SAMPLES DET FRAC IS A FRACTION OF ~600 (VIS) OR ~3700 (NIR) TRACE SAMPLES
DETECTED_FRACTION_ABS = 0.02


@dataclass(frozen=True)
class StareBaseline:
    """The approved values for one stare reduction."""

    stem: str
    obsDate: str
    nOrders: int
    rows: int
    rowsAbs: int
    waveMin: float
    waveMinAbs: float
    waveMax: float
    waveStep: float
    snrMedian: float
    fluxcalMedian: float
    skyCountsMedian: float
    skyModelMedian: float
    skyModelP90: float
    detectedFraction: float
    # THE VIS OBJECT-TRACE PRODUCT IS NAMED FROM THE UPPER-CASED OBJECT NAME
    objtraceStem: str


BASELINES = [
    pytest.param(
        StareBaseline(
            stem="20251013T043511_VIS_1X1_1_STARE_OBJ_SLIT5_0_300_0S_SOXS_Feige110",
            obsDate="2025-10-13T04:35:11.440",
            nOrders=4,
            # 0.02 NM STEPS, SO ±50 ROWS IS ±1 NM OF SPECTRAL COVERAGE
            rows=25_105,
            rowsAbs=50,
            waveMin=347.98,
            waveMinAbs=0.04,
            waveMax=850.06,
            waveStep=0.02,
            snrMedian=139.01,
            fluxcalMedian=6.378164095002718e-14,
            skyCountsMedian=34.83055114746094,
            skyModelMedian=65.08501815795898,
            skyModelP90=874.5562133789062,
            detectedFraction=0.98,
            objtraceStem="20251013T043511_VIS_1X1_1_STARE_OBJ_SLIT5_0_300_0S_SOXS_FEIGE110",
        ),
        id="vis",
    ),
    pytest.param(
        StareBaseline(
            stem="20260111T080347_NIR_3_STARE_OBJ_SLIT5_0_30_0S_SOXS_CD-325613",
            obsDate="2026-01-11T08:03:47.6362",
            nOrders=15,
            # 0.06 NM STEPS, SO ±25 ROWS IS ±1.5 NM OF SPECTRAL COVERAGE
            rows=20_619,
            rowsAbs=25,
            waveMin=795.12,
            waveMinAbs=0.12,
            waveMax=2032.2,
            waveStep=0.06,
            snrMedian=15.46,
            fluxcalMedian=2.9754314300782318e-15,
            skyCountsMedian=20.080078125,
            skyModelMedian=23.640494346618652,
            skyModelP90=135.9283905029297,
            detectedFraction=0.786,
            objtraceStem="20260111T080347_NIR_3_STARE_OBJ_SLIT5_0_30_0S_SOXS_CD-325613",
        ),
        id="nir",
    ),
]


def _stare_product_directory(workspace_path: Path, stem: str) -> Path:
    product_paths = list(workspace_path.rglob(f"{stem}_EXTRACTED_MERGED.fits"))
    assert len(product_paths) == 1, f"expected one merged product for {stem}, found {product_paths}"
    return product_paths[0].parent


def _sky_model_levels(sky_model_path: Path) -> tuple[float, float]:
    """Return the median and 90th percentile of the sky model over good in-order pixels."""
    with fits.open(sky_model_path) as sky_hdus:
        sky = np.asarray(sky_hdus["FLUX"].data, dtype=float)
        quality = np.asarray(sky_hdus["QUAL"].data)
    # THE MODEL IS ZERO BETWEEN THE ORDERS, SO ZEROS WOULD SWAMP THE STATISTICS
    in_order = (quality == 0) & np.isfinite(sky) & (sky != 0)
    assert in_order.any(), f"the sky model in {sky_model_path.name} is empty or zero everywhere"
    return float(np.median(sky[in_order])), float(np.percentile(sky[in_order], 90))


@pytest.mark.parametrize("baseline", BASELINES)
def test_stare_reduction_matches_approved_baseline(reduced_workspace: Path, baseline: StareBaseline) -> None:
    stem = baseline.stem
    product_directory = _stare_product_directory(reduced_workspace, stem)
    expected_products = {
        f"{stem}_EXTRACTED_MERGED.fits",
        f"{stem}_EXTRACTED_MERGED.log",
        f"{stem}_EXTRACTED_MERGED.txt",
        f"{stem}_EXTRACTED_ORDERS.fits",
        f"{stem}_FLUXCAL.fits",
        f"{stem}_SKYMODEL.fits",
        f"{stem}_SKYSUB.fits",
        f"{stem}_SKYSUB_RESIDUALS.fits",
        f"{baseline.objtraceStem}_OBJTRACE.fits",
    }
    found_products = {path.name for path in product_directory.iterdir()}
    assert found_products == expected_products, f"unexpected or missing products: {found_products ^ expected_products}"

    with fits.open(product_directory / f"{stem}_EXTRACTED_MERGED.fits") as merged_hdus:
        assert len(merged_hdus) == 2
        assert merged_hdus[0].header["INSTRUME"] == "SOXS"
        assert merged_hdus[0].header["ESO PRO REC1 ID"] == "soxs-stare"
        assert merged_hdus[0].header["DATE-OBS"] == baseline.obsDate
        merged_table = merged_hdus[1].data
        assert set(merged_table.names) == {
            "WAVE",
            "FLUX_COUNTS",
            "VARIANCE",
            "SKY_COUNTS",
            "SNR",
            "FLUX_DENSITY_COUNTS",
        }
        merged_wave = np.asarray(merged_table["WAVE"], dtype=float)
        merged_rows = len(merged_table)
        snr_median = float(np.nanmedian(merged_table["SNR"]))
        sky_counts_median = float(np.nanmedian(merged_table["SKY_COUNTS"]))

    report("merged rows", merged_rows)
    report("wave min", float(np.nanmin(merged_wave)))
    report("wave max", float(np.nanmax(merged_wave)))
    report("wave step median", float(np.nanmedian(np.diff(merged_wave))))
    report("snr median", snr_median)
    report("sky counts median", sky_counts_median)
    assert merged_rows == pytest.approx(baseline.rows, abs=baseline.rowsAbs)
    # THE BLUE END IS PINNED BY THE ORDER LAYOUT, SO IT IS HELD TO TWO GRID STEPS
    assert float(np.nanmin(merged_wave)) == pytest.approx(baseline.waveMin, abs=baseline.waveMinAbs)
    assert float(np.nanmax(merged_wave)) == pytest.approx(baseline.waveMax, abs=1.0)
    assert float(np.nanmedian(np.diff(merged_wave))) == pytest.approx(baseline.waveStep, abs=1e-6)
    assert snr_median == pytest.approx(baseline.snrMedian, rel=LEVEL_REL)
    assert sky_counts_median == pytest.approx(baseline.skyCountsMedian, rel=LEVEL_REL)

    with fits.open(product_directory / f"{stem}_FLUXCAL.fits") as fluxcal_hdus:
        fluxcal_table = fluxcal_hdus[1].data
        assert len(fluxcal_table) == merged_rows
        assert set(fluxcal_table.names) == {"WAVE", "FLUX_CALIBRATED"}
        fluxcal_median = float(np.nanmedian(fluxcal_table["FLUX_CALIBRATED"]))
    report("fluxcal median", fluxcal_median)
    assert fluxcal_median == pytest.approx(baseline.fluxcalMedian, rel=LEVEL_REL)

    sky_median, sky_p90 = _sky_model_levels(product_directory / f"{stem}_SKYMODEL.fits")
    report("sky model median", sky_median)
    report("sky model p90", sky_p90)
    assert sky_median == pytest.approx(baseline.skyModelMedian, rel=LEVEL_REL)
    assert sky_p90 == pytest.approx(baseline.skyModelP90, rel=LEVEL_REL)

    with sqlite3.connect(reduced_workspace / "soxspipe.db") as connection:
        qc_rows = connection.execute(
            "SELECT qc_name, qc_value FROM quality_control "
            "WHERE soxspipe_recipe = ? AND sof_name = ? "
            "AND qc_name IN (?, ?)",
            ("soxs-stare", f"{stem}.sof", "N ORDERS", "SAMPLES DET FRAC"),
        ).fetchall()
    # ONE ROW PER QC NAME: A DUPLICATE OR A MISSING ROW FAILS HERE RATHER THAN BEING SILENTLY COLLAPSED
    assert sorted(name for name, _ in qc_rows) == ["N ORDERS", "SAMPLES DET FRAC"]
    qc_values = dict(qc_rows)
    report("n orders", float(qc_values["N ORDERS"]))
    report("samples det frac", float(qc_values["SAMPLES DET FRAC"]))
    # EVERY ORDER MUST BE PRESENT: A MISSING ORDER IS A FAILURE, NOT DRIFT
    assert float(qc_values["N ORDERS"]) == baseline.nOrders
    assert float(qc_values["SAMPLES DET FRAC"]) == pytest.approx(baseline.detectedFraction, abs=DETECTED_FRACTION_ABS)
