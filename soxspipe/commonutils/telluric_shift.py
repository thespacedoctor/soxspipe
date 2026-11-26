#!/usr/bin/env python
"""
*measure the wavelength shift of a star's telluric absorption lines, as an informational QC value (DY-1969)*

Author
: David Young

Date Created
: November 26, 2026
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

# THE FEWEST FINITE SAMPLES A LINE WINDOW MUST HOLD FOR THE FIT (CENTRE, WIDTH, DEPTH AND 2 CONTINUUM TERMS)
MIN_WINDOW_SAMPLES = 8
# THE RESOLVING POWER THAT SETS THE FIT WINDOW AND THE STARTING LINE WIDTH. THE FIT WIDTH ITSELF IS FREE
REFERENCE_RESOLVING_POWER = 4000.0
HALF_WINDOW_FWHM = 2.5
FWHM_PER_SIGMA = 2.3548200450309493
# A LINE SHALLOWER THAN THIS (FRACTION OF THE CONTINUUM) IS NOT TREATED AS DETECTED
MIN_LINE_DEPTH = 0.03
# THE LARGEST SHIFT (nm) FROM THE CATALOGUE WAVELENGTH THE FIT MAY REACH; 0.3 nm IS 3 A
MAX_SHIFT_NM = 0.3
# A FITTED CENTRE OR WIDTH THIS CLOSE (AS A FRACTION OF ITS RANGE) TO THE EDGE OF THE RANGE IS A FIT THAT FOUND NO LINE
BOUND_MARGIN = 0.02
# THE LINE WIDTH MAY MOVE THIS FACTOR EITHER SIDE OF THE STARTING GUESS
SIGMA_BOUND_FACTOR = 3.0


@dataclass(frozen=True)
class line_fit:
    """*the fitted centre of one absorption line*

    **Key Arguments:**

    - ``centroidNm`` -- the fitted line centre (nm)
    - ``errorNm`` -- the 1 sigma error of the centre (nm)
    - ``depth`` -- the fractional depth of the line below its continuum
    - ``sigmaNm`` -- the fitted Gaussian width of the line (nm)
    """

    centroidNm: float
    errorNm: float
    depth: float
    sigmaNm: float


def _absorption_line(wavelength, centre, sigma, depth, level, slope, reference):
    """*a straight continuum, anchored at ``reference``, times one Gaussian absorption line*

    **Key Arguments:**

    - ``wavelength`` -- the wavelengths (nm)
    - ``centre`` -- the line centre (nm)
    - ``sigma`` -- the Gaussian width (nm)
    - ``depth`` -- the fractional depth of the line
    - ``level`` -- the continuum level at ``reference``
    - ``slope`` -- the continuum slope per nm
    - ``reference`` -- the wavelength at which the continuum equals ``level`` (nm)

    **Return:**

    - ``flux`` -- the model flux at ``wavelength``
    """
    line = 1.0 - depth * np.exp(-0.5 * ((wavelength - centre) / sigma) ** 2)
    return (level + slope * (wavelength - reference)) * line


def _hits_a_bound(params, lower, upper):
    """*whether the fitted centre or width sits on the edge of its allowed range, so the fit found no line there*

    **Key Arguments:**

    - ``params`` -- the fitted parameters, centre first and width second
    - ``lower`` -- the lower bound of each parameter
    - ``upper`` -- the upper bound of each parameter

    **Return:**

    - ``hitsBound`` -- True if the centre or the width is within ``BOUND_MARGIN`` of its range edge
    """
    for index in (0, 1):
        margin = BOUND_MARGIN * (upper[index] - lower[index])
        if params[index] - lower[index] < margin or upper[index] - params[index] < margin:
            return True
    return False


def fit_line_centroid(wavelength, flux, lineNm, halfWindowNm, sigmaGuessNm=0.1, maxShiftNm=None):
    """*fit one absorption line with a Gaussian on a straight continuum and return its centre*

    **Key Arguments:**

    - ``wavelength`` -- the spectrum wavelengths (nm)
    - ``flux`` -- the spectrum flux
    - ``lineNm`` -- the catalogue wavelength of the line (nm), the centre of the fit window
    - ``halfWindowNm`` -- the half-width of the fit window (nm)
    - ``sigmaGuessNm`` -- the starting line width (nm). Default *0.1*
    - ``maxShiftNm`` -- the largest distance (nm) the fitted centre may lie from ``lineNm``. Default *None*
      (the half-window)

    **Return:**

    - ``fit`` -- the ``line_fit``, or *None* when the window has too few samples or the fit fails or has no absorption
    """
    from scipy.optimize import curve_fit

    wavelength = np.asarray(wavelength, dtype=float)
    flux = np.asarray(flux, dtype=float)
    inWindow = np.isfinite(wavelength) & np.isfinite(flux) & (np.abs(wavelength - lineNm) <= halfWindowNm)
    if inWindow.sum() < MIN_WINDOW_SAMPLES:
        return None
    x, y = wavelength[inWindow], flux[inWindow]
    level = float(np.median(y))
    if level <= 0.0:
        return None
    start = [lineNm, sigmaGuessNm, 0.3, level, 0.0]
    reachNm = halfWindowNm if maxShiftNm is None else min(maxShiftNm, halfWindowNm)
    lower = [lineNm - reachNm, sigmaGuessNm / SIGMA_BOUND_FACTOR, 0.0, 0.0, -np.inf]
    upper = [lineNm + reachNm, sigmaGuessNm * SIGMA_BOUND_FACTOR, 1.0, np.inf, np.inf]
    try:
        params, covariance = curve_fit(
            lambda w, c, s, d, lev, sl: _absorption_line(w, c, s, d, lev, sl, lineNm),
            x,
            y,
            p0=start,
            bounds=(lower, upper),
        )
    except (RuntimeError, ValueError):
        return None
    if _hits_a_bound(params, lower, upper):
        return None
    with np.errstate(invalid="ignore"):
        errorNm = float(np.sqrt(covariance[0, 0]))
    if not np.isfinite(errorNm):
        return None
    return line_fit(centroidNm=float(params[0]), errorNm=errorNm, depth=float(params[2]), sigmaNm=float(params[1]))


# A BAND COUNTS AS MEASURED WHEN IT HAS AT LEAST THIS MANY LINES WHOSE CENTROID ERROR IS BELOW THE LIMIT
# (nm; 0.03 nm = 0.3 A). TUNED ON REAL NIR NOD STANDARDS AND VIS FRAMES (DY-1969). THE TICKET SAID 0.2 A, WHICH LEAVES
# THE FRAMES AT SNR 2 TO 13 NOT MEASURED THOUGH THEIR LINES AGREE ON ONE SHIFT
MIN_BAND_LINES = 5
MAX_LINE_ERROR_NM = 0.03
# THE QC VALUE OF A SHIFT THAT COULD NOT BE MEASURED
NOT_MEASURED = "not measured"
# THE FACTOR THAT TURNS A MEDIAN ABSOLUTE DEVIATION INTO A GAUSSIAN SIGMA
MAD_TO_SIGMA = 1.4826


def robust_error(shiftsNm):
    """*the robust scatter of the line shifts divided by the square root of their number*

    **Key Arguments:**

    - ``shiftsNm`` -- the line shifts (nm)

    **Return:**

    - ``errorNm`` -- 1.4826 times the median absolute deviation, over the square root of N
    """
    shiftsNm = np.asarray(shiftsNm, dtype=float)
    scatter = MAD_TO_SIGMA * np.median(np.abs(shiftsNm - np.median(shiftsNm)))
    return float(scatter / np.sqrt(shiftsNm.size))


def good_lines(lineShifts):
    """*the lines whose centroid error is below ``MAX_LINE_ERROR_NM``*

    **Key Arguments:**

    - ``lineShifts`` -- the line table (``band``, ``shiftNm``, ``errorNm``, ...)

    **Return:**

    - ``goodLines`` -- the rows that pass
    """
    return lineShifts.loc[lineShifts["errorNm"] < MAX_LINE_ERROR_NM]


def summarise_bands(lineShifts, bandNames=()):
    """*summarise the line shifts of each band and apply the measurement gate*

    **Key Arguments:**

    - ``lineShifts`` -- the line table (``band``, ``shiftNm``, ``errorNm``, ...)
    - ``bandNames`` -- bands to report even when none of their lines was fitted. Default ()

    **Return:**

    - ``bands`` -- one row per band: ``nLines`` (lines under the error limit), ``shiftNm`` and ``errorNm`` (NaN
      when the band is not measured) and ``isMeasured``
    """
    import pandas as pd

    good = good_lines(lineShifts)
    names = list(dict.fromkeys([*bandNames, *lineShifts["band"]]))
    rows = []
    for name in names:
        shifts = good.loc[good["band"] == name, "shiftNm"].to_numpy(dtype=float)
        isMeasured = shifts.size >= MIN_BAND_LINES
        rows.append(
            {
                "band": name,
                "nLines": int(shifts.size),
                "shiftNm": float(np.median(shifts)) if isMeasured else np.nan,
                "errorNm": robust_error(shifts) if isMeasured else np.nan,
                "isMeasured": bool(isMeasured),
            }
        )
    return pd.DataFrame(rows, columns=["band", "nLines", "shiftNm", "errorNm", "isMeasured"])


@dataclass(frozen=True)
class telluric_shift_result:
    """*the telluric-line shift of one extraction*

    **Key Arguments:**

    - ``shiftNm`` -- the median shift (nm) of the lines in the measured bands, NaN when no band is measured
    - ``errorNm`` -- the robust scatter of those shifts over the square root of their number, NaN when no
      band is measured
    - ``nLines`` -- the number of lines used, 0 when no band is measured
    - ``bands`` -- the per-band summary (``summarise_bands``)
    - ``lines`` -- every fitted line (``band``, ``lineNm``, ``shiftNm``, ``errorNm``, ``depth``, ``sigmaNm``, ``order``)
    """

    shiftNm: float
    errorNm: float
    nLines: int
    bands: pd.DataFrame
    lines: pd.DataFrame

    # THE FRAMES MAKE THE RESULT COMPARABLE BUT NOT HASHABLE
    __hash__ = None


def summarise_telluric_shift(lineShifts, bandNames=()):
    """*combine the fitted lines of an extraction into the overall shift, its error and the per-band table*

    **Key Arguments:**

    - ``lineShifts`` -- the fitted lines (``band``, ``shiftNm``, ``errorNm``, ...)
    - ``bandNames`` -- bands to report even when none of their lines was fitted. Default ()

    **Return:**

    - ``result`` -- the ``telluric_shift_result``
    """
    bands = summarise_bands(lineShifts, bandNames)
    measuredBands = bands.loc[bands["isMeasured"], "band"]
    used = good_lines(lineShifts)
    used = used.loc[used["band"].isin(measuredBands), "shiftNm"].to_numpy(dtype=float)
    if used.size == 0:
        return telluric_shift_result(np.nan, np.nan, 0, bands, lineShifts)
    return telluric_shift_result(float(np.median(used)), robust_error(used), int(used.size), bands, lineShifts)


def _line_window(lineNm):
    """*the fit half-window and the starting line width of a line, from the reference resolving power*

    **Key Arguments:**

    - ``lineNm`` -- the catalogue wavelength of the line (nm)

    **Return:**

    - ``halfWindowNm`` -- the half-width of the fit window (nm)
    - ``sigmaNm`` -- the starting Gaussian width (nm)
    """
    fwhmNm = lineNm / REFERENCE_RESOLVING_POWER
    return HALF_WINDOW_FWHM * fwhmNm, fwhmNm / FWHM_PER_SIGMA


def _fit_lines_in_order(wavelength, flux, lines):
    """*fit every line whose whole window lies inside one order*

    **Key Arguments:**

    - ``wavelength`` -- the wavelengths of the order, ascending (nm)
    - ``flux`` -- the flux of the order
    - ``lines`` -- the line list rows that fall in the order (``BAND``, ``WAVELENGTH``)

    **Return:**

    - ``rows`` -- one dictionary per fitted line
    """
    rows = []
    for band, lineNm in zip(lines["BAND"], lines["WAVELENGTH"], strict=True):
        halfWindowNm, sigmaNm = _line_window(lineNm)
        if lineNm - halfWindowNm < wavelength[0] or lineNm + halfWindowNm > wavelength[-1]:
            continue
        fit = fit_line_centroid(wavelength, flux, lineNm, halfWindowNm, sigmaGuessNm=sigmaNm, maxShiftNm=MAX_SHIFT_NM)
        if fit is None or fit.depth < MIN_LINE_DEPTH:
            continue
        rows.append(
            {
                "band": band,
                "lineNm": float(lineNm),
                "shiftNm": fit.centroidNm - lineNm,
                "errorNm": fit.errorNm,
                "depth": fit.depth,
                "sigmaNm": fit.sigmaNm,
            }
        )
    return rows


def measure_line_shifts(orderSpectra, lines):
    """*fit every catalogue line in the order spectra and return each line's shift from its catalogue wavelength*

    A line that two overlapping orders both cover is reported once, from the order with the smaller centroid error.

    **Key Arguments:**

    - ``orderSpectra`` -- the extracted orders (``order``, ``wavelengthMean``, ``extractedFluxOptimal``)
    - ``lines`` -- the telluric line list (``BAND``, ``WAVELENGTH`` in nm)

    **Return:**

    - ``lineShifts`` -- one row per fitted line (``band``, ``lineNm``, ``shiftNm``, ``errorNm``, ``depth``, ``sigmaNm``,
      ``order``)
    """
    rows = []
    for order, orderDF in orderSpectra.groupby("order", sort=True):
        orderDF = orderDF.sort_values("wavelengthMean", kind="stable")
        wavelength = orderDF["wavelengthMean"].to_numpy(dtype=float)
        flux = orderDF["extractedFluxOptimal"].to_numpy(dtype=float)
        if wavelength.size < MIN_WINDOW_SAMPLES:
            continue
        inOrder = lines.loc[(lines["WAVELENGTH"] > wavelength[0]) & (lines["WAVELENGTH"] < wavelength[-1])]
        rows.extend({**row, "order": order} for row in _fit_lines_in_order(wavelength, flux, inOrder))
    columns = ["band", "lineNm", "shiftNm", "errorNm", "depth", "sigmaNm", "order"]
    if not rows:
        return pd.DataFrame(columns=columns)
    table = pd.DataFrame(rows, columns=columns).sort_values(["lineNm", "errorNm"], kind="stable")
    return table.drop_duplicates("lineNm", keep="first").reset_index(drop=True)


def record_telluric_qcs(qcTable, result, recipeName, dateObs, reductionDateUtc):
    """*append the telluric-shift QCs to a QC table*

    ``TELL SHIFT`` (nm), ``TELL SHIFT ERR`` (nm) and ``TELL NLINES`` go to the FITS header. When no band is
    measured the two shift rows read ``not measured`` and stay out of the header, so a repeated extraction cannot leave
    an older value behind. Each band adds ``TELL SHIFT <band>``,
    ``TELL SHIFT ERR <band>`` and ``TELL NLINES <band>`` to the QC table only. The values are informational and have
    no acceptable range.

    **Key Arguments:**

    - ``qcTable`` -- the QC table to append to (not mutated)
    - ``result`` -- the ``telluric_shift_result``
    - ``recipeName`` -- the recipe name to record against the rows
    - ``dateObs`` -- the observation date (UTC)
    - ``reductionDateUtc`` -- the reduction date (UTC)

    **Return:**

    - ``qcTable`` -- a new QC table with the rows appended
    """
    from soxspipe.commonutils.toolkit import append_qc

    def add(table, name, value, comment, unit, toHeader):
        return append_qc(
            table,
            recipeName=recipeName,
            qcName=name,
            qcValue=value,
            qcComment=comment,
            obsDateUtc=dateObs,
            reductionDateUtc=reductionDateUtc,
            qcUnit=unit,
            toHeader=toHeader,
        )

    isMeasured = bool(result.nLines)
    qcTable = add(
        qcTable,
        "TELL SHIFT",
        round(result.shiftNm, 5) if isMeasured else NOT_MEASURED,
        "Median telluric-line shift of the star",
        "nm",
        isMeasured,
    )
    qcTable = add(
        qcTable,
        "TELL SHIFT ERR",
        round(result.errorNm, 5) if isMeasured else NOT_MEASURED,
        "Robust scatter of the line shifts over sqrt(N)",
        "nm",
        isMeasured,
    )
    qcTable = add(qcTable, "TELL NLINES", result.nLines, "Telluric lines used for the shift", None, True)
    for band in result.bands.itertuples():
        for name, value, comment in (
            (f"TELL SHIFT {band.band}", band.shiftNm, "Telluric-line shift in the band"),
            (f"TELL SHIFT ERR {band.band}", band.errorNm, "Robust scatter over sqrt(N) in the band"),
        ):
            qcTable = add(qcTable, name, round(value, 5) if band.isMeasured else NOT_MEASURED, comment, "nm", False)
        qcTable = add(
            qcTable,
            f"TELL NLINES {band.band}",
            band.nLines,
            "Telluric lines under the error limit in the band",
            None,
            False,
        )
    return qcTable


def load_telluric_lines(log, settings, arm):
    """*load the shipped telluric line list for one arm*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the soxspipe settings dictionary
    - ``arm`` -- the spectrograph arm

    **Return:**

    - ``lines`` -- the lines of the arm (``BAND``, ``WAVELENGTH`` in nm, ``STRENGTH``, ``SPECIES``), sorted by
      wavelength, or *None* when the instrument ships no telluric line list
    """
    from astropy.table import Table

    from soxspipe.commonutils import detector_lookup
    from soxspipe.commonutils.toolkit import get_calibrations_path

    detector = detector_lookup(log=log, settings=settings).get(arm)
    if "telluric-lines" not in detector:
        return None
    path = get_calibrations_path(log=log, settings=settings) + "/" + detector["telluric-lines"]
    lines = Table.read(path, format="fits").to_pandas()
    for column in ("BAND", "ARM", "SPECIES"):
        lines[column] = (
            lines[column].map(lambda value: value.decode() if isinstance(value, bytes) else value).str.strip()
        )
    lines["WAVELENGTH"] = lines["WAVELENGTH"].astype(float)
    inArm = lines["ARM"] == arm.upper()
    return lines.loc[inArm].sort_values("WAVELENGTH").reset_index(drop=True)


def plot_telluric_shift_qc(result, filePath, title):
    """*plot every line's shift against wavelength, coloured by band, with the median marked*

    **Key Arguments:**

    - ``result`` -- the ``telluric_shift_result``
    - ``filePath`` -- the PDF path to write
    - ``title`` -- the figure title

    **Return:**

    - ``filePath`` -- the path the plot was written to

    **Usage:**

    ```python
    from soxspipe.commonutils.telluric_shift import plot_telluric_shift_qc
    plot_telluric_shift_qc(result, "/path/to/plot.pdf", "NIR telluric-line shift")
    ```
    """
    import matplotlib.pyplot as plt

    angstromPerNm = 10.0
    fig, axis = plt.subplots(figsize=(14, 5))
    try:
        measured = set(result.bands.loc[result.bands["isMeasured"], "band"])
        colours = plt.cm.tab10(np.linspace(0, 1, 10))
        for index, (band, rows) in enumerate(result.lines.groupby("band", sort=False)):
            good = rows["errorNm"] < MAX_LINE_ERROR_NM
            label = band if band in measured else f"{band} (not measured)"
            axis.errorbar(
                rows.loc[good, "lineNm"],
                rows.loc[good, "shiftNm"] * angstromPerNm,
                yerr=rows.loc[good, "errorNm"] * angstromPerNm,
                fmt="o",
                ms=4,
                color=colours[index % 10],
                alpha=1.0 if band in measured else 0.4,
                label=label,
            )
            axis.plot(
                rows.loc[~good, "lineNm"],
                rows.loc[~good, "shiftNm"] * angstromPerNm,
                "x",
                color="grey",
                alpha=0.4,
                ms=4,
            )
        axis.axhline(0.0, color="grey", lw=0.5)
        if result.nLines:
            medianA = result.shiftNm * angstromPerNm
            axis.axhline(medianA, color="black", ls="--", lw=1.2, label=f"median {medianA:+.2f} A, N = {result.nLines}")
            axis.axhspan(
                medianA - result.errorNm * angstromPerNm,
                medianA + result.errorNm * angstromPerNm,
                color="black",
                alpha=0.12,
            )
        else:
            axis.text(0.5, 0.9, "no band measured", transform=axis.transAxes, ha="center")
        axis.set_xlabel("catalogue wavelength of the telluric line, in air (nm)")
        axis.set_ylabel("fitted minus catalogue wavelength (A)")
        axis.set_title(title, fontsize=10)
        if axis.get_legend_handles_labels()[0]:
            axis.legend(fontsize=8, loc="best")
        fig.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    return filePath


def measure_telluric_shift(orderSpectra, lines):
    """*measure the telluric-line shift of the extracted orders of a star*

    **Key Arguments:**

    - ``orderSpectra`` -- the extracted orders (``order``, ``wavelengthMean``, ``extractedFluxOptimal``)
    - ``lines`` -- the telluric line list of the arm (``BAND``, ``WAVELENGTH`` in nm)

    **Return:**

    - ``result`` -- the ``telluric_shift_result``. Every band of the line list is in its band table
    """
    lineShifts = measure_line_shifts(orderSpectra, lines)
    return summarise_telluric_shift(lineShifts, list(dict.fromkeys(lines["BAND"])))
