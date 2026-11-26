#!/usr/bin/env python
"""*build the telluric absorption line list that soxspipe ships as `telluric_lines.fits`*

This is a developer tool for DY-1969. It is not part of a reduction and not in the soxspipe command line. `hapi`, the
HITRAN Application Programming Interface, is imported inside `fetch_tables` only. It is a development tool, not a
dependency of the package.

**Usage**

Run it with the `soxspipe` environment interpreter. Put the worktree and a directory that holds a copy of the `hapi`
package (from `hitran-api`) on `PYTHONPATH`. Install nothing into the shared environment. The HITRAN fetches are cached
in `--cache`, so a rerun needs no network. Give the lowest measured per-slit resolving power of each arm, which is the
`SLIT1.5` value (2811 for VIS and 3462 for NIR). Do not give the `SLIT5.0` values (954 and 1229). They are derived from
`SLIT1.0`, not measured, and they give no lines at all: the line fit cannot follow lines that wide.

```bash
PYTHONPATH=/path/to/worktree:/path/to/dir/with/hapi \\
    ~/anaconda/envs/soxspipe/bin/python tools/build_telluric_lines.py \\
    --cache /path/to/hitran/cache --lowest-r-vis 2811 --lowest-r-nir 3462 \\
    --out soxspipe/resources/static_calibrations/soxs/telluric_lines.fits
```

A run prints the number of lines per band and names every band that it dropped. Rerun with the same arguments to
reproduce the shipped file. Its header records the model settings and the selection limits.

HAPI fetches from HITRAN with the host that `hapi.VARIABLES["GLOBAL_HOST"]` names, and its default is plain HTTP
(`http://hitran.org`). The script sets the host to `https://hitran.org` before every fetch. On 2026-11-26 a rebuild
from a fresh cache fetched over HTTPS gave a file byte-identical to the shipped one.

**How the lines are chosen**

1. Source. HITRAN line-by-line parameters from hitran.org, through HAPI. Principal isotopologues only: H2O (molecule 1)
   over the VIS and NIR spans and O2 (7) over the VIS span. HITRAN gives vacuum wavenumbers. The list stores wavelengths
   in air (nm, Ciddor 1996 standard air), the frame of the SOXS wavelength scale and of the sky-line list.
2. Model. Each band is one homogeneous absorbing layer (270 K, 0.75 atm, airmass 1.3) with Voigt profiles. H2O is also
   run at 0.6 and 1.6 times its column, because water vapour varies by a large factor from night to night. The model is
   for line selection only. It is never fitted to data.
3. Simulated measurement. The model is smoothed to a set of resolving powers that runs from the lowest one of the arm
   up to `--highest-r`, and the shipped fit (`soxspipe.commonutils.telluric_shift.fit_line_centroid`) runs on every
   HITRAN line of the band at every resolving power and water column.
4. Selection. A HITRAN line is kept when, in all of the simulated measurements:
   - the fit succeeds and its depth is between `--min-depth` and `--max-depth` (detectable, not saturated);
   - the fitted centre moves by at most `--max-spread` nm (a blend with a neighbour moves the centre as the resolution
     changes);
   - the mean fitted centre is within `--max-offset` nm of the HITRAN wavelength.
   `WAVELENGTH` is the mean fitted centre. `LINE_NM` is the HITRAN air wavelength. Two kept lines whose fitted centres
   are closer than `--merge-nm` are one feature, and the stronger one is kept.
5. Bands. Only lines inside the bands of `BANDS` are kept. A band with fewer than `--min-band-lines` lines is dropped.
"""

from __future__ import annotations

import argparse
import os
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import OptimizeWarning
from scipy.signal import fftconvolve
from scipy.special import wofz

from soxspipe.commonutils import telluric_shift as ts


class Band(NamedTuple):
    """*a telluric band: its name, its wavelength range (nm, air), its arm and the fetched spans that hold its lines*"""

    name: str
    firstNm: float
    lastNm: float
    arm: str
    spanNames: tuple[str, ...]


class Candidate(NamedTuple):
    """*a HITRAN line inside a band: its air wavelength (nm), its strength and its species*"""

    lineNm: float
    strength: float
    species: str


class Feature(NamedTuple):
    """*a kept line: its mean fitted centre (nm), its HITRAN air wavelength (nm), its strength and its species*"""

    fittedNm: float
    lineNm: float
    strength: float
    species: str


class Fetch(NamedTuple):
    """*one HITRAN fetch: the name of its table in the HAPI cache and its wavenumber range (cm-1, vacuum)*"""

    tableName: str
    firstCm: float
    lastCm: float


class Molecule(NamedTuple):
    """*what the model needs of one molecule: name, mass (u), slant column (molecule/cm2), isotopologue, fetches*"""

    name: str
    mass: float
    column: float
    isotopologue: int
    fetches: tuple[Fetch, ...]


BANDS = [
    Band("O2_B", 686.0, 697.0, "VIS", ("H2O_VIS", "O2_VIS")),
    Band("H2O_720", 714.0, 735.0, "VIS", ("H2O_VIS", "O2_VIS")),
    Band("H2O_940", 920.0, 980.0, "NIR", ("H2O_VIS",)),
    Band("H2O_1130", 1110.0, 1160.0, "NIR", ("H2O_NIR",)),
    Band("H2O_1450", 1400.0, 1500.0, "NIR", ("H2O_NIR",)),
    Band("H2O_1800", 1750.0, 1850.0, "NIR", ("H2O_NIR",)),
]

# THE HITRAN MOLECULE IDS. THE TABLE NAME OF A FETCH IS THE KEY OF THE HAPI CACHE. "VIS" IS 10000 TO 14700 CM-1
# (680 TO 1000 NM) AND "NIR" IS 4780 TO 9120 CM-1 (1096 TO 2092 NM). THE SLANT COLUMNS ARE AT THE AIRMASS BELOW, AND THE
# H2O COLUMN IS ABOUT 3 MM OF PRECIPITABLE WATER. A NEW MOLECULE IS ONE MORE ENTRY HERE
WATER_ID = 1
OXYGEN_ID = 7
MOLECULES = {
    WATER_ID: Molecule(
        "H2O", 18.0106, 1.3e22, 1, (Fetch("H2O_VIS", 10000.0, 14700.0), Fetch("H2O_NIR", 4780.0, 9120.0))
    ),
    OXYGEN_ID: Molecule("O2", 31.9898, 4.5e24, 1, (Fetch("O2_VIS", 10000.0, 14700.0),)),
}
HITRAN_FETCH_DATE = "2026-11-26"
# HAPI DEFAULTS TO PLAIN HTTP; EVERY FETCH GOES TO THIS HTTPS HOST INSTEAD
HITRAN_HOST = "https://hitran.org"
HAPI_COLUMNS = ("nu", "sw", "gamma_air", "n_air", "elower", "delta_air")

# SINGLE-LAYER ATMOSPHERE: TEMPERATURE (K), PRESSURE (ATM) AND AIRMASS
TEMPERATURE = 270.0
PRESSURE = 0.75
AIRMASS = 1.3
WATER_COLUMN_FACTORS = (0.6, 1.0, 1.6)

# HITRAN LINES WEAKER THAN THIS (CM-1/(MOLECULE CM-2) AT 296 K) NEVER GIVE A DETECTABLE LINE
MIN_LINE_STRENGTH = 1e-27
# SECOND RADIATION CONSTANT, CM K
SECOND_RADIATION_CONSTANT = 1.4387769
# THE TEMPERATURE OF THE HITRAN LINE STRENGTHS, K
REFERENCE_TEMPERATURE = 296.0
DOPPLER_COEFFICIENT = 3.58e-7
# THE WAVENUMBER GRID STEP, CM-1
GRID_STEP = 0.01
# A LINE REACHES THIS FAR EITHER SIDE OF ITS CENTRE, CM-1
LINE_WING = 3.0
# WAVELENGTH (NM) TIMES WAVENUMBER (CM-1): NM = NM_CM / CM-1
NM_CM = 1e7
# NM IN ONE MICRON, FOR THE CIDDOR AIR FORMULA
NM_PER_MICRON = 1e3
# THE FIXED-POINT STEPS THAT TURN AN AIR WAVELENGTH BACK INTO A VACUUM WAVELENGTH
AIR_INVERSE_ITERATIONS = 3
BAND_PAD_NM = 3.0  # MODEL THIS FAR BEYOND THE BAND, SO THE FIT WINDOWS AT ITS EDGES ARE FILLED
THINNING = 4  # THE FIT RUNS ON EVERY FOURTH GRID POINT
KERNEL_HALF_WIDTHS = 5  # THE CONVOLUTION KERNEL REACHES THIS MANY SIGMA
PREFILTER_FRACTION = 0.5  # A LINE LESS THAN THIS FRACTION OF THE MINIMUM DEPTH IN THE MIDDLE MODEL CANNOT PASS
RESOLVING_POWER_COUNT = 3
DEFAULT_HIGHEST_R = 5000.0
SELECTION_DEFAULTS = {
    "minDepth": 0.04,
    "maxDepth": 0.85,
    "maxSpreadNm": 0.015,
    "maxOffsetNm": 0.05,
    "mergeNm": 0.03,
}
LINE_COLUMNS = ("BAND", "ARM", "WAVELENGTH", "LINE_NM", "STRENGTH", "SPECIES")


@dataclass(frozen=True)
class Selection:
    """*the selection limits and the resolving powers of each arm*"""

    resolvingPowers: Mapping[str, tuple[float, ...]]
    minDepth: float
    maxDepth: float
    maxSpreadNm: float
    maxOffsetNm: float
    mergeNm: float
    minBandLines: int


def resolving_power_grid(lowest: float, highest: float, count: int = RESOLVING_POWER_COUNT) -> tuple[float, ...]:
    """*resolving powers spaced evenly in the logarithm, from the lowest to the highest*

    **Key Arguments:**

    - ``lowest`` -- the lowest resolving power
    - ``highest`` -- the highest resolving power
    - ``count`` -- how many values to return. Default ``RESOLVING_POWER_COUNT``

    **Return:**

    - ``powers`` -- ``count`` values from ``lowest`` to ``highest``, rounded to whole numbers

    **Raises:**

    - ``ValueError`` -- if ``lowest`` is not above 0, or ``highest`` is below ``lowest``, or ``count`` is below 2
    """
    if lowest <= 0 or highest < lowest or count < 2:
        raise ValueError(f"cannot span resolving powers {lowest} to {highest} with {count} values")
    return tuple(float(round(power)) for power in np.geomspace(lowest, highest, count))


def vac_to_air_nm(wavelengthVac: Any) -> Any:
    """*vacuum wavelength to standard-air wavelength (nm), Ciddor 1996*

    **Key Arguments:**

    - ``wavelengthVac`` -- the vacuum wavelength (nm), a number or an array

    **Return:**

    - ``wavelengthAir`` -- the air wavelength (nm)
    """
    sigma2 = (NM_PER_MICRON / np.asarray(wavelengthVac, dtype=float)) ** 2
    refractivity = 0.05792105 / (238.0185 - sigma2) + 0.00167917 / (57.362 - sigma2)
    return wavelengthVac / (1.0 + refractivity)


def vac_to_air_inverse(airNm: float) -> float:
    """*air wavelength (nm) to the matching vacuum wavelength (nm), by fixed-point steps*

    **Key Arguments:**

    - ``airNm`` -- the air wavelength (nm)

    **Return:**

    - ``vacNm`` -- the vacuum wavelength (nm)
    """
    vac = airNm
    for _ in range(AIR_INVERSE_ITERATIONS):
        vac = airNm * (vac / float(vac_to_air_nm(vac)))
    return vac


class HitranFetchError(RuntimeError):
    """*a HITRAN fetch failed, or returned no lines*"""


def drop_cached_table(hapi: Any, cache: str, tableName: str) -> None:
    """*forget a HAPI table and delete its files, so that the next run fetches it again*

    **Key Arguments:**

    - ``hapi`` -- the ``hapi`` module
    - ``cache`` -- the HAPI cache folder
    - ``tableName`` -- the name of the table
    """
    hapi.dropTable(tableName)
    for suffix in (".data", ".header"):
        (Path(cache) / f"{tableName}{suffix}").unlink(missing_ok=True)


def read_table(hapi: Any, cache: str, molecule: Molecule, fetch: Fetch) -> dict[str, Any]:
    """*read the line columns of one fetched table, and refuse a table with no lines*

    **Key Arguments:**

    - ``hapi`` -- the ``hapi`` module
    - ``cache`` -- the HAPI cache folder
    - ``molecule`` -- the molecule of the table
    - ``fetch`` -- the fetch that made the table

    **Return:**

    - ``table`` -- the numpy columns of ``HAPI_COLUMNS``

    **Raises:**

    - ``HitranFetchError`` -- if the table holds no lines. The table is deleted from the cache first
    """
    table = {column: np.array(hapi.getColumn(fetch.tableName, column)) for column in HAPI_COLUMNS}
    if table["nu"].size == 0:
        drop_cached_table(hapi, cache, fetch.tableName)
        raise HitranFetchError(
            f"HITRAN returned no {molecule.name} lines for {fetch.firstCm} to {fetch.lastCm} cm-1 "
            f"(table {fetch.tableName}); the table was deleted from {cache}"
        )
    return table


def fetch_tables(cache: str) -> dict[str, dict[str, Any]]:
    """*fetch every span with HAPI, cached in ``cache``, and return its line columns*

    HAPI queries the host in ``hapi.VARIABLES["GLOBAL_HOST"]``, which defaults to plain HTTP. This function sets it to
    ``HITRAN_HOST`` (HTTPS) before any fetch, so the line data cannot be changed in transit.

    **Key Arguments:**

    - ``cache`` -- the HAPI cache folder. It is created when it does not exist

    **Return:**

    - ``tables`` -- per table name, a dictionary of numpy columns and the HITRAN ``molecule`` id

    **Raises:**

    - ``HitranFetchError`` -- if a fetch fails or returns no lines. The message names the molecule and the range, and
      the partial table is deleted from the cache so that the next run fetches it again
    """
    import hapi

    hapi.VARIABLES["GLOBAL_HOST"] = HITRAN_HOST
    os.makedirs(cache, exist_ok=True)
    hapi.db_begin(cache)
    tables = {}
    for moleculeId, molecule in MOLECULES.items():
        for fetch in molecule.fetches:
            if fetch.tableName not in hapi.tableList():
                try:
                    hapi.fetch(fetch.tableName, moleculeId, molecule.isotopologue, fetch.firstCm, fetch.lastCm)
                except Exception as error:
                    drop_cached_table(hapi, cache, fetch.tableName)
                    raise HitranFetchError(
                        f"the HITRAN fetch of {molecule.name} for {fetch.firstCm} to {fetch.lastCm} cm-1 "
                        f"(table {fetch.tableName}) failed: {error}"
                    ) from error
            tables[fetch.tableName] = {**read_table(hapi, cache, molecule, fetch), "molecule": moleculeId}
    return tables


def voigt_profile(offsets: np.ndarray, doppler: float, lorentz: float) -> np.ndarray:
    """*the area-normalised Voigt profile (cm) at ``offsets`` (cm-1)*

    **Key Arguments:**

    - ``offsets`` -- the distances from the line centre (cm-1)
    - ``doppler`` -- the Doppler half width at half maximum (cm-1)
    - ``lorentz`` -- the Lorentz half width at half maximum (cm-1)

    **Return:**

    - ``profile`` -- the profile at each offset
    """
    sigma = doppler / np.sqrt(2.0 * np.log(2.0))
    z = (offsets + 1j * lorentz) / (sigma * np.sqrt(2.0))
    return wofz(z).real / (sigma * np.sqrt(2.0 * np.pi))


def line_opacity(table: Mapping[str, Any], grid: np.ndarray, columnFactor: float = 1.0) -> np.ndarray:
    """*the optical depth on ``grid`` of every line of one fetched table*

    **Key Arguments:**

    - ``table`` -- one fetched span (``fetch_tables``)
    - ``grid`` -- the wavenumber grid (cm-1), ascending
    - ``columnFactor`` -- the factor on the H2O column. Other molecules ignore it. Default 1.0

    **Return:**

    - ``tau`` -- the optical depth at each grid point
    """
    molecule = table["molecule"]
    nu = table["nu"] + table["delta_air"] * PRESSURE
    strength = table["sw"] * np.exp(
        -SECOND_RADIATION_CONSTANT * table["elower"] * (1.0 / TEMPERATURE - 1.0 / REFERENCE_TEMPERATURE)
    )
    lorentz = table["gamma_air"] * PRESSURE * (REFERENCE_TEMPERATURE / TEMPERATURE) ** table["n_air"]
    doppler = DOPPLER_COEFFICIENT * nu * np.sqrt(TEMPERATURE / MOLECULES[molecule].mass)
    column = MOLECULES[molecule].column * (columnFactor if molecule == WATER_ID else 1.0)
    tau = np.zeros_like(grid)
    for index in np.flatnonzero((nu > grid[0] - LINE_WING) & (nu < grid[-1] + LINE_WING)):
        left = np.searchsorted(grid, nu[index] - LINE_WING)
        right = np.searchsorted(grid, nu[index] + LINE_WING)
        tau[left:right] += (
            column * strength[index] * voigt_profile(grid[left:right] - nu[index], doppler[index], lorentz[index])
        )
    return tau


def convolve_to_resolution(grid: np.ndarray, signal: np.ndarray, resolvingPower: float) -> np.ndarray:
    """*smooth a signal with a Gaussian of constant resolving power, evaluated at the grid centre*

    **Key Arguments:**

    - ``grid`` -- the wavenumber grid (cm-1)
    - ``signal`` -- the signal on the grid
    - ``resolvingPower`` -- the resolving power to smooth to

    **Return:**

    - ``smoothed`` -- the smoothed signal
    """
    centre = 0.5 * (grid[0] + grid[-1])
    sigmaPixels = (centre / resolvingPower) / (2.0 * np.sqrt(2.0 * np.log(2.0))) / GRID_STEP
    half = int(np.ceil(KERNEL_HALF_WIDTHS * sigmaPixels))
    kernel = np.exp(-0.5 * (np.arange(-half, half + 1) / sigmaPixels) ** 2)
    return fftconvolve(signal, kernel / kernel.sum(), mode="same")


def band_models(
    spans: Sequence[Mapping[str, Any]], grid: np.ndarray, resolvingPowers: Sequence[float], factors: Sequence[float]
) -> list[np.ndarray]:
    """*the smoothed transmission on ``grid`` for every water-column factor and resolving power*

    **Key Arguments:**

    - ``spans`` -- the fetched spans that hold the band's lines
    - ``grid`` -- the wavenumber grid (cm-1), ascending
    - ``resolvingPowers`` -- the resolving powers to smooth to
    - ``factors`` -- the factors on the H2O column

    **Return:**

    - ``models`` -- one transmission per factor and resolving power
    """
    models = []
    for factor in factors:
        transmission = np.exp(-sum(line_opacity(span, grid, factor) for span in spans))
        models.extend(convolve_to_resolution(grid, transmission, power) for power in resolvingPowers)
    return models


def simulated_centres(
    waveAir: np.ndarray, models: Sequence[np.ndarray], lineNm: float
) -> tuple[np.ndarray, np.ndarray] | None:
    """*fit one line in every model*

    **Key Arguments:**

    - ``waveAir`` -- the model wavelengths (nm, air), ascending
    - ``models`` -- the transmission models on ``waveAir``
    - ``lineNm`` -- the HITRAN wavelength of the line (nm, air)

    **Return:**

    - ``fits`` -- the fitted centres (nm) and depths of the line in every model, or *None* when any fit fails
    """
    halfWindowNm, sigmaNm = ts.line_window(lineNm)
    centres, depths = [], []
    for transmission in models:
        fit = ts.fit_line_centroid(
            waveAir, transmission, lineNm, halfWindowNm, sigmaGuessNm=sigmaNm, maxShiftNm=ts.MAX_SHIFT_NM
        )
        if fit is None:
            return None
        centres.append(fit.centroidNm)
        depths.append(fit.depth)
    return np.array(centres), np.array(depths)


def merge_close_features(features: Sequence[Feature], mergeNm: float) -> list[Feature]:
    """*keep the strongest of any features whose fitted centres lie within ``mergeNm`` of each other*

    **Key Arguments:**

    - ``features`` -- the kept lines
    - ``mergeNm`` -- the merge distance (nm)

    **Return:**

    - ``kept`` -- the surviving features, in order of wavelength
    """
    kept: list[Feature] = []
    for feature in sorted(features):
        if kept and feature.fittedNm - kept[-1].fittedNm < mergeNm:
            if feature.strength > kept[-1].strength:
                kept[-1] = feature
            continue
        kept.append(feature)
    return kept


def band_candidates(spans: Sequence[Mapping[str, Any]], firstNm: float, lastNm: float) -> list[Candidate]:
    """*the HITRAN lines inside a band, as wavelength in air (nm), strength and species*

    **Key Arguments:**

    - ``spans`` -- the fetched spans that hold the band's lines
    - ``firstNm`` -- the first wavelength of the band (nm, air)
    - ``lastNm`` -- the last wavelength of the band (nm, air)

    **Return:**

    - ``candidates`` -- one ``Candidate`` per line that is strong enough to matter
    """
    candidates = []
    for span in spans:
        wavelength = vac_to_air_nm(NM_CM / span["nu"])
        inBand = (wavelength >= firstNm) & (wavelength <= lastNm) & (span["sw"] > MIN_LINE_STRENGTH)
        candidates.extend(
            Candidate(float(wavelength[index]), float(span["sw"][index]), MOLECULES[span["molecule"]].name)
            for index in np.flatnonzero(inBand)
        )
    return candidates


def passes_selection(centres: np.ndarray, depths: np.ndarray, lineNm: float, selection: Selection) -> bool:
    """*whether the fits of one line in all the models meet the selection limits*

    **Key Arguments:**

    - ``centres`` -- the fitted centres (nm)
    - ``depths`` -- the fitted depths
    - ``lineNm`` -- the HITRAN wavelength of the line (nm, air)
    - ``selection`` -- the limits

    **Return:**

    - ``passes`` -- True if the depths, the spread of the centres and their mean offset are all within the limits
    """
    if depths.min() < selection.minDepth or depths.max() > selection.maxDepth:
        return False
    if centres.max() - centres.min() > selection.maxSpreadNm:
        return False
    return abs(float(centres.mean()) - lineNm) <= selection.maxOffsetNm


def select_band(tables: Mapping[str, Mapping[str, Any]], band: Band, selection: Selection) -> list[Feature]:
    """*the lines kept in one band*

    **Key Arguments:**

    - ``tables`` -- the fetched spans (``fetch_tables``)
    - ``band`` -- an entry of ``BANDS``
    - ``selection`` -- the limits and resolving powers

    **Return:**

    - ``kept`` -- the kept lines
    """
    spans = [tables[name] for name in band.spanNames]
    candidates = band_candidates(spans, band.firstNm, band.lastNm)
    if not candidates:
        return []
    nuLast = NM_CM / vac_to_air_inverse(band.firstNm - BAND_PAD_NM)
    nuFirst = NM_CM / vac_to_air_inverse(band.lastNm + BAND_PAD_NM)
    grid = np.arange(nuFirst, nuLast, GRID_STEP)
    waveAir = vac_to_air_nm(NM_CM / grid)[::-1]
    factors = (
        WATER_COLUMN_FACTORS
        if any(candidate.species == MOLECULES[WATER_ID].name for candidate in candidates)
        else (1.0,)
    )
    models = [model[::-1] for model in band_models(spans, grid, selection.resolvingPowers[band.arm], factors)]
    thin = slice(None, None, THINNING)
    middle = models[len(models) // 2]
    kept = []
    for candidate in candidates:
        if 1.0 - np.interp(candidate.lineNm, waveAir, middle) < PREFILTER_FRACTION * selection.minDepth:
            continue
        result = simulated_centres(waveAir[thin], [model[thin] for model in models], candidate.lineNm)
        if result is not None and passes_selection(*result, candidate.lineNm, selection):
            kept.append(Feature(float(result[0].mean()), candidate.lineNm, candidate.strength, candidate.species))
    return merge_close_features(kept, selection.mergeNm)


def drop_sparse_bands(
    rowsByBand: Mapping[str, list[tuple]], minBandLines: int
) -> tuple[dict[str, list[tuple]], list[str]]:
    """*remove every band that has fewer lines than ``minBandLines``*

    **Key Arguments:**

    - ``rowsByBand`` -- the table rows of each band
    - ``minBandLines`` -- the fewest lines a band keeps

    **Return:**

    - ``kept`` -- the rows of the bands that stay
    - ``dropped`` -- the names of the bands that were removed
    """
    kept = {band: rows for band, rows in rowsByBand.items() if len(rows) >= minBandLines}
    return kept, [band for band in rowsByBand if band not in kept]


def model_header_cards(selection: Selection, dropped: Sequence[str]) -> list[tuple[str, Any, str]]:
    """*the HIERARCH cards that record the model settings, the selection limits and the dropped bands*

    **Key Arguments:**

    - ``selection`` -- the limits and resolving powers
    - ``dropped`` -- the names of the bands that were dropped

    **Return:**

    - ``cards`` -- tuples of keyword, value and comment
    """
    cards = [
        ("HIERARCH MODEL TEMPERATURE", TEMPERATURE, "[K] one-layer atmosphere"),
        ("HIERARCH MODEL PRESSURE", PRESSURE, "[atm]"),
        ("HIERARCH MODEL AIRMASS", AIRMASS, "the columns below are slant columns at this airmass"),
        ("HIERARCH MODEL H2O COLUMN", MOLECULES[WATER_ID].column, "[molecule/cm2] about 3 mm of precipitable water"),
        (
            "HIERARCH MODEL H2O COLUMN FACTORS",
            ",".join(f"{factor:g}" for factor in WATER_COLUMN_FACTORS),
            "times the column",
        ),
        ("HIERARCH MODEL O2 COLUMN", MOLECULES[OXYGEN_ID].column, "[molecule/cm2]"),
    ]
    cards += [
        (f"HIERARCH MODEL RESOLVING POWERS {arm}", ",".join(f"{power:.0f}" for power in powers), "smoothed to these")
        for arm, powers in selection.resolvingPowers.items()
    ]
    cards += [
        ("HIERARCH SELECT MIN DEPTH", selection.minDepth, "fractional line depth"),
        ("HIERARCH SELECT MAX DEPTH", selection.maxDepth, "fractional line depth"),
        ("HIERARCH SELECT MAX SPREAD", selection.maxSpreadNm, "[nm] largest move of the fitted centre"),
        ("HIERARCH SELECT MAX OFFSET", selection.maxOffsetNm, "[nm] mean fitted centre from HITRAN"),
        ("HIERARCH SELECT MERGE DISTANCE", selection.mergeNm, "[nm] closer features are one"),
        ("HIERARCH SELECT MIN BAND LINES", selection.minBandLines, "fewer lines drop the band"),
        ("HIERARCH SELECT DROPPED BANDS", ",".join(dropped) or "none", "bands left with too few lines"),
    ]
    return cards


def build_table(rowsByBand: Mapping[str, list[tuple]]) -> Table:
    """*the line table, with units and descriptions*

    **Key Arguments:**

    - ``rowsByBand`` -- the table rows (``BAND``, ``ARM``, ``WAVELENGTH``, ``LINE_NM``, ``STRENGTH``, ``SPECIES``) of
      each band

    **Return:**

    - ``table`` -- the ``TELLURIC_LINES`` table
    """
    rows = [row for bandRows in rowsByBand.values() for row in bandRows]
    table = Table(rows=rows, names=LINE_COLUMNS, dtype=("U10", "U3", "f8", "f8", "f4", "U4"))
    table["WAVELENGTH"].unit = "nm"
    table["WAVELENGTH"].description = "Reference wavelength in air: mean fitted centre of the simulated measurements"
    table["LINE_NM"].unit = "nm"
    table["LINE_NM"].description = "HITRAN wavelength of the line in air"
    table["STRENGTH"].description = "HITRAN line intensity sw at 296 K, cm-1/(molecule cm-2)"
    table.meta["EXTNAME"] = "TELLURIC_LINES"
    return table


def write_line_list(table: Table, selection: Selection, dropped: Sequence[str], path: str) -> None:
    """*write the line table and the model settings to a FITS file*

    **Key Arguments:**

    - ``table`` -- the ``TELLURIC_LINES`` table
    - ``selection`` -- the limits and resolving powers
    - ``dropped`` -- the names of the bands that were dropped
    - ``path`` -- the FITS file to write
    """
    hdu = fits.table_to_hdu(table)
    hdu.header["HISTORY"] = f"HITRAN via HAPI, fetched {HITRAN_FETCH_DATE}"
    hdu.header["HISTORY"] = "built by tools/build_telluric_lines.py (DY-1969)"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", fits.verify.VerifyWarning)
        for keyword, value, comment in model_header_cards(selection, dropped):
            hdu.header[keyword] = (value, comment)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


def build_line_list(tables: Mapping[str, Mapping[str, Any]], selection: Selection) -> tuple[Table, list[str]]:
    """*select the lines of every band and drop the bands with too few*

    **Key Arguments:**

    - ``tables`` -- the fetched spans (``fetch_tables``)
    - ``selection`` -- the limits and resolving powers

    **Return:**

    - ``table`` -- the ``TELLURIC_LINES`` table
    - ``dropped`` -- the names of the bands that were dropped
    """
    rowsByBand = {}
    for band in BANDS:
        # CURVE_FIT WARNS WHEN IT CANNOT ESTIMATE THE COVARIANCE OF A FLAT WINDOW. THE FIT IS THEN REJECTED, SO IT IS
        # EXPECTED. EVERY OTHER WARNING STAYS VISIBLE
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            features = select_band(tables, band, selection)
        rowsByBand[band.name] = [
            (band.name, band.arm, feature.fittedNm, feature.lineNm, feature.strength, feature.species)
            for feature in features
        ]
    kept, dropped = drop_sparse_bands(rowsByBand, selection.minBandLines)
    return build_table(kept), dropped


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """*read the command line*

    **Key Arguments:**

    - ``argv`` -- the arguments. Default *None*, which reads ``sys.argv``

    **Return:**

    - ``args`` -- the parsed arguments
    """
    parser = argparse.ArgumentParser(description="Build the telluric line list telluric_lines.fits.")
    parser.add_argument("--cache", required=True, help="the HAPI cache folder")
    parser.add_argument("--out", required=True, help="the FITS file to write")
    parser.add_argument(
        "--lowest-r-vis", type=float, required=True, help="lowest measured per-slit resolving power of the VIS arm"
    )
    parser.add_argument(
        "--lowest-r-nir", type=float, required=True, help="lowest measured per-slit resolving power of the NIR arm"
    )
    parser.add_argument("--highest-r", type=float, default=DEFAULT_HIGHEST_R)
    parser.add_argument("--min-depth", type=float, default=SELECTION_DEFAULTS["minDepth"])
    parser.add_argument("--max-depth", type=float, default=SELECTION_DEFAULTS["maxDepth"])
    parser.add_argument("--max-spread", type=float, default=SELECTION_DEFAULTS["maxSpreadNm"], help="nm")
    parser.add_argument("--max-offset", type=float, default=SELECTION_DEFAULTS["maxOffsetNm"], help="nm")
    parser.add_argument("--merge-nm", type=float, default=SELECTION_DEFAULTS["mergeNm"])
    parser.add_argument("--min-band-lines", type=int, default=ts.MIN_BAND_LINES)
    return parser.parse_args(argv)


def selection_from_arguments(args: argparse.Namespace) -> Selection:
    """*the selection that the command line asks for*

    **Key Arguments:**

    - ``args`` -- the parsed arguments (``parse_arguments``)

    **Return:**

    - ``selection`` -- the limits, with each arm's resolving powers running from its lowest value to ``--highest-r``
    """
    powers = {
        "VIS": resolving_power_grid(args.lowest_r_vis, args.highest_r),
        "NIR": resolving_power_grid(args.lowest_r_nir, args.highest_r),
    }
    return Selection(
        powers,
        args.min_depth,
        args.max_depth,
        args.max_spread,
        args.max_offset,
        args.merge_nm,
        args.min_band_lines,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """*build the line list, print the lines per band and the dropped bands, and write the FITS file*

    **Key Arguments:**

    - ``argv`` -- the arguments. Default *None*, which reads ``sys.argv``

    **Return:**

    - ``status`` -- 0
    """
    args = parse_arguments(argv)
    selection = selection_from_arguments(args)
    for arm, powers in selection.resolvingPowers.items():
        print(f"{arm} resolving powers: {', '.join(f'{power:.0f}' for power in powers)}")
    table, dropped = build_line_list(fetch_tables(args.cache), selection)
    for band in BANDS:
        count = int((table["BAND"] == band.name).sum()) if len(table) else 0
        print(f"{band.name:10s} {band.firstNm:7.1f}-{band.lastNm:7.1f} nm  {count:4d} lines")
    print(f"dropped bands: {', '.join(dropped) or 'none'}")
    write_line_list(table, selection, dropped, args.out)
    print(f"wrote {len(table)} lines to {args.out} ({os.path.getsize(args.out)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
