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

Re-recorded again on 2026-10-16 from two bit-identical CI runs (38915352201 and
38915356992) in the DY-1612 pull request. These moves predate DY-1612: a gate run on
develop itself (38918109330, bc80f1d1) gave the same values. They arrived through
merges into develop, where the gate does not run. The VIS blue end moved from 348.02
to 348.14 nm and its sky counts median fell 7%; the NIR blue end moved from 794.58 to
795.06 nm and its SNR median fell from 16.91 to 14.09.

The two `fluxcalMedian` values were re-recorded on 2026-10-18 from six bit-identical
CI runs (39165616944, 39165618635, 39165620359, 39165622153, 39165623841 and
39165625784), after the response curve became a log-space spline (DY-1614). The VIS
median rose 2.9% and the NIR median fell 3.2%.

The NIR `fluxcalMedian` was re-recorded on 2026-10-19 from six bit-identical CI runs
(39290451388, 39290461991, 39290472137, 39290482158, 39290492461 and 39290503644),
after the NIR response fit change (DY-1647). It rose 0.6%. The VIS values moved by
0.01% or less, inside their bands, and were left as they are (DY-1648).

The extinction correction became linear instead of a staircase on 2026-10-21 (DY-1684).
In gate run 39616315251 the VIS `fluxcalMedian` fell 0.41%, inside its band, and was left
as it is. NIR does not use the extinction correction.

Both `fluxcalMedian` values were re-recorded on 2026-10-21 from six bit-identical CI
runs (39680454564, 39680459968, 39680464192, 39680468399, 39680472907 and 39680476407),
after DY-1674, DY-1675, DY-1679, DY-1680 and DY-1684 changed the standard-star
responses. The VIS median fell 0.42% (DY-1684) and the NIR median rose 0.12% (DY-1690,
DY-1691).
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
            rows=25_096,
            rowsAbs=50,
            waveMin=348.14,
            waveMinAbs=0.04,
            waveMax=850.04,
            waveStep=0.02,
            snrMedian=156.8,
            fluxcalMedian=6.078459105261144e-14,
            skyCountsMedian=6.347151279449463,
            skyModelMedian=9.057429313659668,
            skyModelP90=93.78917541503911,
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
            rows=20_617,
            rowsAbs=25,
            waveMin=795.06,
            waveMinAbs=0.12,
            waveMax=2032.02,
            waveStep=0.06,
            snrMedian=14.09,
            fluxcalMedian=3.003861084057503e-15,
            skyCountsMedian=17.42337417602539,
            skyModelMedian=20.28994369506836,
            skyModelP90=89.84078826904295,
            detectedFraction=0.939,
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


def _measure_stare_reduction(
    reduced_workspace: Path, product_directory: Path, baseline: StareBaseline
) -> dict[str, float]:
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
