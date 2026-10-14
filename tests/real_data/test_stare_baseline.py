"""Acceptance checks for the real VIS and NIR stare reductions.

The stare recipe is the only recipe that models the sky with the b-spline fit in
`subtract_sky`, so these checks are the end-to-end guard on that fit. As with the offset
check, the reduction is not bit-reproducible across hardware, so each value is asserted
within a band. The sky bands are narrow enough that a sky model of zero, or one doubled,
fails them. The centres were re-recorded on 2026-10-14 from six bit-identical CI
runs (38566513736, 38566517851, 38566520911, 38566524861,
38566528464 and 38566531516), after the DY-1279 object-wing correction and the
DY-1518 background mask (DY-1521). DY-1279 removed object light from the VIS sky
model of this bright standard, so its sky levels fell by about 85%. A local arm64
reduction lands up to 5% away, so expect local runs to sit near the band edges. The
NIR detected fraction and blue end swing with small upstream changes (DY-994).
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
            rows=25_102,
            rowsAbs=50,
            waveMin=348.02,
            waveMinAbs=0.04,
            waveMax=850.04,
            waveStep=0.02,
            snrMedian=156.91,
            fluxcalMedian=5.925303507351052e-14,
            skyCountsMedian=6.8395094871521,
            skyModelMedian=9.14101791381836,
            skyModelP90=93.79542541503906,
            detectedFraction=0.931,
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
            rows=20_620,
            rowsAbs=25,
            waveMin=794.58,
            waveMinAbs=0.12,
            waveMax=2031.72,
            waveStep=0.06,
            snrMedian=16.91,
            fluxcalMedian=3.0643544845287745e-15,
            skyCountsMedian=17.535743713378906,
            skyModelMedian=20.37826919555664,
            skyModelP90=89.83494873046874,
            detectedFraction=0.938,
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


def _measure_stare_reduction(reduced_workspace: Path, product_directory: Path, baseline: StareBaseline) -> dict:
    """Read every baselined value of one stare reduction, checking product structure on the way."""
    stem = baseline.stem
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
        measured = {
            "merged rows": len(merged_table),
            "snr median": float(np.nanmedian(merged_table["SNR"])),
            "sky counts median": float(np.nanmedian(merged_table["SKY_COUNTS"])),
        }
    measured["wave min"] = float(np.nanmin(merged_wave))
    measured["wave max"] = float(np.nanmax(merged_wave))
    measured["wave step median"] = float(np.nanmedian(np.diff(merged_wave)))

    with fits.open(product_directory / f"{stem}_FLUXCAL.fits") as fluxcal_hdus:
        fluxcal_table = fluxcal_hdus[1].data
        assert set(fluxcal_table.names) == {"WAVE", "FLUX_CALIBRATED"}
        measured["fluxcal rows"] = len(fluxcal_table)
        measured["fluxcal median"] = float(np.nanmedian(fluxcal_table["FLUX_CALIBRATED"]))

    measured["sky model median"], measured["sky model p90"] = _sky_model_levels(
        product_directory / f"{stem}_SKYMODEL.fits"
    )

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
    measured["n orders"] = float(qc_values["N ORDERS"])
    measured["samples det frac"] = float(qc_values["SAMPLES DET FRAC"])
    return measured


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

    measured = _measure_stare_reduction(reduced_workspace, product_directory, baseline)
    # EVERY VALUE IS REPORTED BEFORE THE FIRST BAND ASSERTION, SO A FAILING RUN STILL
    # RECORDS ALL OF THEM FOR THE NEXT RE-RECORD
    for label, value in measured.items():
        report(label, value)

    assert measured["merged rows"] == pytest.approx(baseline.rows, abs=baseline.rowsAbs)
    # THE BLUE END IS PINNED BY THE ORDER LAYOUT, SO IT IS HELD TO TWO GRID STEPS
    assert measured["wave min"] == pytest.approx(baseline.waveMin, abs=baseline.waveMinAbs)
    assert measured["wave max"] == pytest.approx(baseline.waveMax, abs=1.0)
    assert measured["wave step median"] == pytest.approx(baseline.waveStep, abs=1e-6)
    assert measured["snr median"] == pytest.approx(baseline.snrMedian, rel=LEVEL_REL)
    assert measured["sky counts median"] == pytest.approx(baseline.skyCountsMedian, rel=LEVEL_REL)
    assert measured["fluxcal rows"] == measured["merged rows"]
    assert measured["fluxcal median"] == pytest.approx(baseline.fluxcalMedian, rel=LEVEL_REL)
    assert measured["sky model median"] == pytest.approx(baseline.skyModelMedian, rel=LEVEL_REL)
    assert measured["sky model p90"] == pytest.approx(baseline.skyModelP90, rel=LEVEL_REL)
    # EVERY ORDER MUST BE PRESENT: A MISSING ORDER IS A FAILURE, NOT DRIFT
    assert measured["n orders"] == baseline.nOrders
    assert measured["samples det frac"] == pytest.approx(baseline.detectedFraction, abs=DETECTED_FRACTION_ABS)
