"""Acceptance checks for the calibration and standard-star reductions of the real-data gate (DY-257).

The gate reduces three science SOFs, and each `reduce sof` builds the calibration chain
that SOF needs. The science baselines catch a calibration change only when it moves the
final spectrum. These checks hold every calibration SOF and both standard-star `nod`
SOFs to their own QC values, read from the workspace's `quality_control` table, and
check the response curve of each standard star. The gate reduces 23 calibration SOFs
and 2 standard-star SOFs.

Three QCs are checked only where the recipe writes them. The VIS order-centre runs write
no `Y RES SD`, and the VIS master flats write no `N LOW SENS`. The dispersion solution
writes no QC for its missing arc lines (`DISP_MAP_LINES_MISSING` is a product, the
`_MISSED_LINES.fits` table), so its detected-line count `DETLINES NUM` is held exactly
instead. `MDARK MEDIAN` sits within 2 e- of zero, so its 5% band is a few hundredths of
an electron wide; it holds only because the CI runs are bit-identical.

The response product holds polynomial coefficients, so the check evaluates the curve on
a fixed 1 nm grid inside each arm with the pipeline's own flux-calibration function, and
bands the median and three fixed wavelengths. The polynomial diverges outside the arm,
so the grids stop well short of the arm edges.

Bands follow the stare convention: levels, residuals and RMS values to 5% relative,
fractions to ±0.02 absolute, and counts exactly. The centres were recorded from six
bit-identical CI runs on RECORDING_DATE. A local arm64 reduction lands up to 5% away,
so expect local runs to sit near the band edges.
"""

from __future__ import annotations

import sqlite3
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from soxspipe.commonutils.flux_calibration import _calculate_flux_calibration
from tests.real_data.reporting import report

pytestmark = pytest.mark.real_data

# A LEVEL, RESIDUAL OR RMS MAY DRIFT BY THIS FRACTION ACROSS RUNNER HARDWARE
LEVEL_REL = 0.05
FRACTION_ABS = 0.02


@dataclass(frozen=True)
class Band:
    """An approved value and how far a run may land from it."""

    centre: float
    relTol: float = 0.0
    absTol: float = 0.0


def level(centre: float) -> Band:
    """A level, residual or RMS value, held to `LEVEL_REL`."""
    return Band(centre, relTol=LEVEL_REL)


def fraction(centre: float) -> Band:
    """A fraction, held to `FRACTION_ABS`."""
    return Band(centre, absTol=FRACTION_ABS)


def count(centre: int) -> Band:
    """A count, held exactly."""
    return Band(float(centre))


@dataclass(frozen=True)
class QcBaseline:
    """The approved QC values of one reduced SOF."""

    recipe: str
    sofName: str
    values: Mapping[str, Band]


@dataclass(frozen=True)
class ResponseBaseline:
    """The approved response curve of one standard-star reduction, in flux per count per second."""

    sofName: str
    gridStart: float
    gridStop: float
    median: Band
    # WAVELENGTH IN NM TO THE APPROVED CURVE VALUE THERE
    points: Mapping[float, Band]


QC_BASELINES: tuple[QcBaseline, ...] = (
    QcBaseline(
        recipe="soxs-mbias",
        sofName="20250831T021035_VIS_1X1_1_MBIAS_SOXS.sof",
        values={
            "MBIAS MEDIAN": level(644.4989764939007),
            "MASTER RON": level(1.901401162147522),
        },
    ),
    QcBaseline(
        recipe="soxs-mbias",
        sofName="20250901T100431_VIS_1X1_1_MBIAS_SOXS.sof",
        values={
            "MBIAS MEDIAN": level(643.7105077242613),
            "MASTER RON": level(1.7507745027542114),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20250514T120301_NIR_3_MDARK_5_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-0.94981895210743),
            "HOTPIX FRAC": fraction(0.009302),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20250514T120422_NIR_3_MDARK_10_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-1.61726377464756),
            "HOTPIX FRAC": fraction(0.009088),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20250514T120611_NIR_3_MDARK_15_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-0.6358310977120494),
            "HOTPIX FRAC": fraction(0.008801),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20251128T105414_NIR_3_MDARK_2_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-1.5857252304520588),
            "HOTPIX FRAC": fraction(0.009712),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20251128T105637_NIR_3_MDARK_10_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(0.4038254614750811),
            "HOTPIX FRAC": fraction(0.010009),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20251128T105827_NIR_3_MDARK_15_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-1.7331866021651137),
            "HOTPIX FRAC": fraction(0.009858),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20250514T121150_NIR_3_DSOL_PINHOLE_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.29231),
            "GOODLINES FRAC": fraction(0.835616),
            "DETLINES NUM": count(611),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20250831T030313_VIS_1X1_1_DSOL_PINHOLE_30_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.35245),
            "GOODLINES FRAC": fraction(0.932836),
            "DETLINES NUM": count(130),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20250901T105805_VIS_1X1_1_DSOL_PINHOLE_30_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.43292),
            "GOODLINES FRAC": fraction(0.932836),
            "DETLINES NUM": count(129),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20251128T122146_NIR_3_DSOL_PINHOLE_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.29715),
            "GOODLINES FRAC": fraction(0.844749),
            "DETLINES NUM": count(602),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20250514T091410_NIR_3_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(15),
            "SAMPLES DET FRAC": fraction(0.971),
            "Y RES SD": level(0.043),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20250831T024046_VIS_1X1_1_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(4),
            "SAMPLES DET FRAC": fraction(0.928),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20250901T103524_VIS_1X1_1_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(4),
            "SAMPLES DET FRAC": fraction(0.938),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20251128T121139_NIR_3_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(15),
            "SAMPLES DET FRAC": fraction(0.988),
            "Y RES SD": level(0.044),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20250514T085927_NIR_3_MFLAT_QTH_SLIT1_0_3_75S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.917),
            "COLDPIX FRAC": fraction(0.009302),
            "N LOW SENS": count(0),
            "ORDEXP50": level(6378.843),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20250831T023807_VIS_1X1_1_MFLAT_QTH_SLIT5_0_2_0S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.965),
            "COLDPIX FRAC": fraction(0.005825),
            "ORDEXP50": level(3882.34),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20250901T103245_VIS_1X1_1_MFLAT_QTH_SLIT5_0_2_0S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.964),
            "COLDPIX FRAC": fraction(0.004728),
            "ORDEXP50": level(3904.34),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20251128T114919_NIR_3_MFLAT_QTH_SLIT1_5_2_5S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.928),
            "COLDPIX FRAC": fraction(0.009712),
            "N LOW SENS": count(0),
            "ORDEXP50": level(12719.776),
        },
    ),
    QcBaseline(
        recipe="soxs-spat-solution",
        sofName="20250514T121227_NIR_3_SSOL_MULTPIN_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.22811),
            "GOODLINES FRAC": fraction(0.557594),
            "PINHOLE COUNT MIN": count(9),
        },
    ),
    QcBaseline(
        recipe="soxs-spat-solution",
        sofName="20250831T030421_VIS_1X1_1_SSOL_MULTPIN_30_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.26075),
            "GOODLINES FRAC": fraction(0.854892),
            "PINHOLE COUNT MIN": count(9),
        },
    ),
    QcBaseline(
        recipe="soxs-spat-solution",
        sofName="20251128T122217_NIR_3_SSOL_MULTPIN_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.23285),
            "GOODLINES FRAC": fraction(0.457696),
            "PINHOLE COUNT MIN": count(9),
        },
    ),
    QcBaseline(
        recipe="soxs-nod-std",
        sofName="20250901T030144_VIS_1X1_1_NOD_STD_FLUX_SLIT5_0_300_0S_SOXS.sof",
        values={
            "EFF MEDIAN": level(0.1273),
            "SNR MEDIAN": level(161.41),
            "N ORDERS": count(4),
        },
    ),
    QcBaseline(
        recipe="soxs-nod-std",
        sofName="20251128T030328_NIR_3_NOD_STD_FLUX_SLIT5_0_150_0S_SOXS.sof",
        values={
            "EFF MEDIAN": level(0.1158),
            "SNR MEDIAN": level(50.14),
            "N ORDERS": count(15),
        },
    ),
)

RESPONSE_BASELINES: tuple[ResponseBaseline, ...] = (
    ResponseBaseline(
        sofName="20250901T030144_VIS_1X1_1_NOD_STD_FLUX_SLIT5_0_300_0S_SOXS.sof",
        gridStart=400.0,
        gridStop=800.0,
        median=level(7.075162382086183e-16),
        points={
            450.0: level(2.847147794003331e-16),
            520.0: level(4.990828057279169e-16),
            750.0: level(1.0858371357198532e-15),
        },
    ),
    ResponseBaseline(
        sofName="20251128T030328_NIR_3_NOD_STD_FLUX_SLIT5_0_150_0S_SOXS.sof",
        gridStart=1000.0,
        gridStop=1800.0,
        median=level(2.647135203631308e-16),
        points={
            1050.0: level(2.861237714107415e-16),
            1250.0: level(2.879974378992802e-16),
            1650.0: level(1.8112117393047925e-16),
        },
    ),
)

# THE RECIPES WHOSE SOFS THE TABLE ABOVE MUST COVER IN FULL. THE SCIENCE RECIPES HAVE
# THEIR OWN BASELINE FILES
BASELINED_RECIPES = (
    "soxs-mbias",
    "soxs-mdark",
    "soxs-disp-solution",
    "soxs-order-centres",
    "soxs-mflat",
    "soxs-spat-solution",
    "soxs-nod-std",
)


def _connect_read_only(workspacePath: Path) -> sqlite3.Connection:
    # READ-ONLY, SO A MISSING DATABASE RAISES INSTEAD OF BEING CREATED EMPTY
    return sqlite3.connect(f"{(workspacePath / 'soxspipe.db').as_uri()}?mode=ro", uri=True)


def test_every_calibration_and_standard_star_sof_has_a_baseline(reduced_workspace: Path) -> None:
    connection = _connect_read_only(reduced_workspace)
    try:
        reducedRows = connection.execute("SELECT DISTINCT soxspipe_recipe, sof_name FROM quality_control").fetchall()
    finally:
        connection.close()
    reduced = {(recipe, sofName) for recipe, sofName in reducedRows if recipe in BASELINED_RECIPES}
    baselined = {(baseline.recipe, baseline.sofName) for baseline in QC_BASELINES}
    # A SOF THE GATE STOPS REDUCING, OR A NEW ONE WITHOUT A BASELINE, FAILS HERE
    assert reduced == baselined, (
        f"reduced but not baselined: {reduced - baselined}; baselined but not reduced: {baselined - reduced}"
    )


@pytest.mark.parametrize("baseline", QC_BASELINES, ids=lambda baseline: baseline.sofName)
def test_calibration_qc_matches_approved_baseline(reduced_workspace: Path, baseline: QcBaseline) -> None:
    qcNames = tuple(baseline.values)
    connection = _connect_read_only(reduced_workspace)
    try:
        # ONLY THE WHOLE-PRODUCT VALUE: SOME QCS ALSO HAVE ONE ROW PER ORDER
        sofRows = connection.execute(
            "SELECT qc_name, qc_value FROM quality_control "
            "WHERE soxspipe_recipe = ? AND sof_name = ? AND qc_order = '-1'",
            (baseline.recipe, baseline.sofName),
        ).fetchall()
    finally:
        connection.close()
    qcRows = [(name, value) for name, value in sofRows if name in baseline.values]

    # ONE ROW PER QC NAME: A DUPLICATE OR A MISSING ROW FAILS HERE RATHER THAN BEING SILENTLY COLLAPSED
    assert sorted(name for name, _ in qcRows) == sorted(qcNames)
    measured = {name: float(value) for name, value in qcRows}
    for name in qcNames:
        report(f"{baseline.sofName} {name}", measured[name])
    for name, band in baseline.values.items():
        assert measured[name] == pytest.approx(band.centre, rel=band.relTol, abs=band.absTol), name


def _response_curve(responsePath: Path, wavelengths: np.ndarray) -> np.ndarray:
    """Evaluate the response function on a wavelength grid, the way `flux_calibration` applies it."""
    coefficientTable = Table.read(responsePath, format="fits")
    coefficients = [coefficientTable[f"c{index}"][0] for index in range(int(coefficientTable["polyOrder"][0]) + 1)]
    # ONE COUNT PER SECOND WITH NO EXTINCTION, SO THE RESULT IS THE FLUX ONE COUNT PER SECOND STANDS FOR
    return np.asarray(
        _calculate_flux_calibration(
            wavelengths, counts=1.0, exposureTime=1.0, responseCoefficients=coefficients, extinctionFactors=1.0
        ),
        dtype=float,
    )


@pytest.mark.parametrize("baseline", RESPONSE_BASELINES, ids=lambda baseline: baseline.sofName)
def test_standard_star_response_curve_matches_approved_baseline(
    reduced_workspace: Path, baseline: ResponseBaseline
) -> None:
    stem = baseline.sofName.removesuffix(".sof")
    responsePaths = list(reduced_workspace.rglob(f"{stem}_RESP.fits"))
    assert len(responsePaths) == 1, f"expected one response product for {stem}, found {responsePaths}"

    grid = np.arange(baseline.gridStart, baseline.gridStop + 0.5, 1.0)
    curveMedian = float(np.median(_response_curve(responsePaths[0], grid)))
    pointWavelengths = np.array(sorted(baseline.points))
    pointValues = dict(
        zip(pointWavelengths.tolist(), _response_curve(responsePaths[0], pointWavelengths).tolist(), strict=True)
    )
    report(f"{stem} response median", curveMedian)
    for wavelength, value in pointValues.items():
        report(f"{stem} response at {wavelength:g} nm", value)

    assert curveMedian == pytest.approx(baseline.median.centre, rel=baseline.median.relTol, abs=baseline.median.absTol)
    for wavelength, band in baseline.points.items():
        assert pointValues[wavelength] == pytest.approx(band.centre, rel=band.relTol, abs=band.absTol), wavelength
