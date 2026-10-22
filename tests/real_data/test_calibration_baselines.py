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

The check reads the response curve from the `RESPONSE` table of the response product,
the way flux calibration does. It evaluates the curve on a fixed 1 nm grid inside each
arm with the pipeline's own flux-calibration function, and bands the median and three
fixed wavelengths inside the grid. The grids stop well short of the arm edges, where
the curve is extrapolated. NIR also bands one point outside its grid, at 2020 nm in the
K band. The DY-1647 telluric masks blank that region, so no data constrain the curve
there and the value is the spline's extrapolation. The point guards it anyway, because
the DY-1674 and DY-1679 regressions moved it. Its centre was recorded on 2026-10-22 from
six bit-identical CI runs (39770082366, 39770085713, 39770089302, 39770092950,
39770096520 and 39770099985) of develop 55faf0fd plus a report-only probe (DY-1712).

Bands follow the stare convention: levels, residuals and RMS values to 5% relative,
fractions to ±0.02 absolute, and counts exactly. The centres were re-recorded on
2026-10-14 from six bit-identical CI runs (38566513736, 38566517851, 38566520911,
38566524861, 38566528464 and 38566531516), after the exact dispersion fit (DY-1525),
per-frame stacking masks (DY-1254) and the DY-1518 background mask (DY-1524). Since
DY-1525, local arm64 reductions give the same `XY RES MEDIAN` as CI. The three VIS
dispersion and spatial solution entries were re-recorded on 2026-10-16 from two
identical CI runs (38869883513 and 38869886400), after three blended Th lines were
marked deleted in the VIS line lists (DY-1578). The two response baselines were
re-recorded on 2026-10-18 from six bit-identical CI runs (39165616944, 39165618635,
39165620359, 39165622153, 39165623841 and 39165625784), after the response curve
became a log-space spline (DY-1614). The NIR median fell 8.9%, and every other
response value moved by 4% or less.

The NIR response baseline was re-recorded on 2026-10-19 from six bit-identical CI
runs (39290451388, 39290461991, 39290472137, 39290482158, 39290492461 and
39290503644), after the NIR fit gained wider telluric masks, telluric wing growth
and a smoothing penalty (DY-1647). The NIR median rose 2.8%, 1050 nm rose 2.3%,
1250 nm fell 0.3% and 1650 nm fell 1.3%. The VIS curve moved by about 0.02% and was
left as it is (DY-1648).

The extinction correction became linear instead of a staircase on 2026-10-21 (DY-1684).
In gate run 39616315251 the VIS response median fell 0.44%, 450 nm fell 0.58%, 520 nm
fell 0.25% and 750 nm fell 0.42%. All stay inside their bands and were left as they are.
NIR does not use the extinction correction and is unchanged.

Both standard-star reductions were re-recorded on 2026-10-21 from six bit-identical CI
runs (39680454564, 39680459968, 39680464192, 39680468399, 39680472907 and 39680476407),
after the nod shift-and-add (DY-1674), the variance of the median stack (DY-1675), the
flux-conserving slit shift (DY-1679), the flat-pass sky shifts reused by the
not-flattened pass (DY-1680) and the linear extinction correction (DY-1684). The `SNR
MEDIAN` of both standards rose by about the square root of 2 (VIS 41.2%, NIR 40.1%),
because the stacked nod path now averages A-B and B-A before extraction (DY-1674). Both
standards combined the nod pair (`NOD COMBINED` 1, `NOD COMBINE FALLBACK NUM` 0), so
DY-1675 adds nothing. The VIS response values moved by -0.2% to -0.6% (DY-1684) and the
VIS `EFF MEDIAN` by -0.3%. The NIR response values moved by 0.2% or less: DY-1674 raised
them by up to 0.5% and DY-1679 and DY-1680 took most of that back (DY-1690, DY-1691).
"""

from __future__ import annotations

import logging
import sqlite3
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from soxspipe.commonutils.flux_calibration import _calculate_flux_calibration, load_response_factors
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
            "MBIAS MEDIAN": level(644.498976493901),
            "MASTER RON": level(1.82982501991057),
        },
    ),
    QcBaseline(
        recipe="soxs-mbias",
        sofName="20250901T100431_VIS_1X1_1_MBIAS_SOXS.sof",
        values={
            "MBIAS MEDIAN": level(643.710507724261),
            "MASTER RON": level(1.68993050289757),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20250514T120301_NIR_3_MDARK_5_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-0.94981895210743),
            "HOTPIX FRAC": fraction(0.008368),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20250514T120422_NIR_3_MDARK_10_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-1.61726377464756),
            "HOTPIX FRAC": fraction(0.008138),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20250514T120611_NIR_3_MDARK_15_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-0.635831097712049),
            "HOTPIX FRAC": fraction(0.007809),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20251128T105414_NIR_3_MDARK_2_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-1.58572523045206),
            "HOTPIX FRAC": fraction(0.008151),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20251128T105637_NIR_3_MDARK_10_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(0.403825461475081),
            "HOTPIX FRAC": fraction(0.008577),
        },
    ),
    QcBaseline(
        recipe="soxs-mdark",
        sofName="20251128T105827_NIR_3_MDARK_15_0S_SOXS.sof",
        values={
            "MDARK MEDIAN": level(-1.73318660216511),
            "HOTPIX FRAC": fraction(0.008408),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20250514T121150_NIR_3_DSOL_PINHOLE_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.25505),
            "GOODLINES FRAC": fraction(0.835616),
            "DETLINES NUM": count(611),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20250831T030313_VIS_1X1_1_DSOL_PINHOLE_30_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.23581),
            "GOODLINES FRAC": fraction(0.931298),
            "DETLINES NUM": count(127),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20250901T105805_VIS_1X1_1_DSOL_PINHOLE_30_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.19831),
            "GOODLINES FRAC": fraction(0.931298),
            "DETLINES NUM": count(126),
        },
    ),
    QcBaseline(
        recipe="soxs-disp-solution",
        sofName="20251128T122146_NIR_3_DSOL_PINHOLE_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.29097),
            "GOODLINES FRAC": fraction(0.838661),
            "DETLINES NUM": count(603),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20250514T091410_NIR_3_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(15),
            "SAMPLES DET FRAC": fraction(0.971),
            "Y RES SD": level(0.042),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20250831T024046_VIS_1X1_1_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(4),
            "SAMPLES DET FRAC": fraction(0.935),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20250901T103524_VIS_1X1_1_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(4),
            "SAMPLES DET FRAC": fraction(0.94),
        },
    ),
    QcBaseline(
        recipe="soxs-order-centres",
        sofName="20251128T121139_NIR_3_OLOC_QTH_PINHOLE_10_0S_SOXS.sof",
        values={
            "N ORDERS": count(15),
            "SAMPLES DET FRAC": fraction(0.987),
            "Y RES SD": level(0.044),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20250514T085927_NIR_3_MFLAT_QTH_SLIT1_0_3_75S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.921),
            "COLDPIX FRAC": fraction(0.008368),
            "N LOW SENS": count(0),
            "ORDEXP50": level(6382.758),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20250831T023807_VIS_1X1_1_MFLAT_QTH_SLIT5_0_2_0S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.971),
            "COLDPIX FRAC": fraction(0.005784),
            "ORDEXP50": level(3883.01),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20250901T103245_VIS_1X1_1_MFLAT_QTH_SLIT5_0_2_0S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.97),
            "COLDPIX FRAC": fraction(0.004902),
            "ORDEXP50": level(3902.36),
        },
    ),
    QcBaseline(
        recipe="soxs-mflat",
        sofName="20251128T114919_NIR_3_MFLAT_QTH_SLIT1_5_2_5S_SOXS.sof",
        values={
            "INNER ORDER PIX MEAN": level(0.935),
            "COLDPIX FRAC": fraction(0.008151),
            "N LOW SENS": count(0),
            "ORDEXP50": level(12719.765),
        },
    ),
    QcBaseline(
        recipe="soxs-spat-solution",
        sofName="20250514T121227_NIR_3_SSOL_MULTPIN_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.21581),
            "GOODLINES FRAC": fraction(0.555556),
            "PINHOLE COUNT MIN": count(9),
        },
    ),
    QcBaseline(
        recipe="soxs-spat-solution",
        sofName="20250831T030421_VIS_1X1_1_SSOL_MULTPIN_30_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.18952),
            "GOODLINES FRAC": fraction(0.868533),
            "PINHOLE COUNT MIN": count(9),
        },
    ),
    QcBaseline(
        recipe="soxs-spat-solution",
        sofName="20251128T122217_NIR_3_SSOL_MULTPIN_15_0S_SOXS.sof",
        values={
            "XY RES MEDIAN": level(0.23585),
            "GOODLINES FRAC": fraction(0.458308),
            "PINHOLE COUNT MIN": count(9),
        },
    ),
    QcBaseline(
        recipe="soxs-nod-std",
        sofName="20250901T030144_VIS_1X1_1_NOD_STD_FLUX_SLIT5_0_300_0S_SOXS.sof",
        values={
            "EFF MEDIAN": level(0.1406),
            # 167.53 UNTIL 2026-10-16; 207.34 IN CI RUNS 38915352201 AND 38915356992 AND ON develop
            # ITSELF (38918109330), SO THE MOVE PREDATES DY-1612. 292.84 SINCE 2026-10-21 (DY-1674)
            "SNR MEDIAN": level(292.84),
            "N ORDERS": count(4),
        },
    ),
    QcBaseline(
        recipe="soxs-nod-std",
        sofName="20251128T030328_NIR_3_NOD_STD_FLUX_SLIT5_0_150_0S_SOXS.sof",
        values={
            "EFF MEDIAN": level(0.1181),
            "SNR MEDIAN": level(70.14),
            "N ORDERS": count(15),
        },
    ),
)

RESPONSE_BASELINES: tuple[ResponseBaseline, ...] = (
    ResponseBaseline(
        sofName="20250901T030144_VIS_1X1_1_NOD_STD_FLUX_SLIT5_0_300_0S_SOXS.sof",
        gridStart=400.0,
        gridStop=800.0,
        median=level(6.631187603114011e-16),
        points={
            450.0: level(2.535990712539509e-16),
            520.0: level(4.413322319722507e-16),
            750.0: level(1.007297831677685e-15),
        },
    ),
    ResponseBaseline(
        sofName="20251128T030328_NIR_3_NOD_STD_FLUX_SLIT5_0_150_0S_SOXS.sof",
        gridStart=1000.0,
        gridStop=1800.0,
        median=level(2.3624028979488e-16),
        points={
            1050.0: level(2.966543301877743e-16),
            1250.0: level(2.813704561775645e-16),
            1650.0: level(1.8297532492538825e-16),
            # EXTRAPOLATED, MASKED-REGION GUARD: 2020 NM LIES INSIDE A DY-1647 K-BAND TELLURIC
            # MASK, SO NO DATA CONSTRAIN THE CURVE HERE (DY-1712)
            2020.0: level(1.0080307076018144e-16),
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
    responseFactors = load_response_factors(str(responsePath), wavelengths, logging.getLogger(__name__))
    # ONE COUNT PER SECOND WITH NO EXTINCTION, SO THE RESULT IS THE FLUX ONE COUNT PER SECOND STANDS FOR
    return np.asarray(
        _calculate_flux_calibration(1.0, exposureTime=1.0, responseFactors=responseFactors, extinctionFactors=1.0),
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
