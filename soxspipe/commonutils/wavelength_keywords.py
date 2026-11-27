"""*the wavelength-calibration and resolution cards of a merged spectrum (DY-2358)*

``LAMRMS``, ``LAMNLIN`` and ``SPEC_ERR`` describe how well the dispersion solution fits the arc lines. They are read
from the ``WAVE RMS`` and ``WAVE NLINES`` QCs of the ``DISP_TAB`` that the science frames were reduced with: the recipe
already holds the path of that map, so the merged spectrum always reports the map its extraction used.
``SPEC_ERR = LAMRMS / sqrt(LAMNLIN)``, as the Phase 3 standard says.

``SPEC_RES`` is the resolving power through the science slit. It comes from the nominal ``spec-res`` table of the
detector parameters, by arm and slit name. It is never the pinhole value ``R PIN MEDIAN``, which is for a 0.5 arcsec
aperture and is far above the value for the slit.

``SPEC_SYE`` is the largest of the per-order ``SKY SHIFT O<n> ERR`` QCs, the uncertainty of the sky shift applied to
each order (DY-2307, DY-2359). The extraction writes those QCs, and for a stack each one is already the largest over
the frames or cycles. Only the VIS arm writes them so far, so a NIR spectrum has no ``SPEC_SYE`` yet.

A card that cannot be computed is left out, with one warning in the recipe log. ``SPEC_RES`` is mandatory, so the
warning says the arm fails: the Phase 3 checker fails the missing card, and ``soxs-join`` treats the arm as a failed
arm (DY-2336). ``LAMRMS``, ``LAMNLIN`` and ``SPEC_ERR`` are optional or recommended, so the product ships without them.

X-Shooter has no ``spec-res`` table and its code paths do not change, so these cards are written for SOXS only.

This module imports nothing from a recipe.

Author
: David Young

Date Created
: November 27, 2026
"""

from __future__ import annotations

import copy
import math
import re
from typing import Any

SOXS_INSTRUMENT = "soxs"

LAMRMS_CARD = "LAMRMS"
LAMNLIN_CARD = "LAMNLIN"
SPEC_ERR_CARD = "SPEC_ERR"
SPEC_RES_CARD = "SPEC_RES"
SPEC_SYE_CARD = "SPEC_SYE"
ACCURACY_CARDS = (LAMRMS_CARD, LAMNLIN_CARD, SPEC_ERR_CARD)

# THE PER-ORDER QC THE SKYLINE TUNING WRITES: THE UNCERTAINTY (nm) OF THE SHIFT APPLIED TO AN ORDER
SKY_SHIFT_ERROR_QC = re.compile(r"SKY SHIFT O\d+ ERR")

# THE QC CARDS OF A DISP_TAB THAT `create_dispersion_map.calculate_residuals` WRITES
WAVE_RMS_QC_CARD = "ESO QC WAVE RMS"
WAVE_NLINES_QC_CARD = "ESO QC WAVE NLINES"

# THE SUFFIX THE NIR JH BLOCKING FILTER ADDS TO A SLIT NAME, WHICH THE spec-res TABLE LEAVES OUT
JH_SUFFIX = "JH"

Cards = list[tuple[str, Any, str]]


def _is_soxs(settings: Any) -> bool:
    """*return True unless the settings select another instrument*"""
    return str(settings.get("instrument", SOXS_INSTRUMENT)).lower() == SOXS_INSTRUMENT


def _read_wave_qcs(dispMapPath: Any, log: Any) -> tuple[Any, Any] | None:
    """*read the two wave QC cards of a dispersion map, or None with a warning when it cannot give both*"""
    from astropy.io import fits

    leftOut = f"so {LAMRMS_CARD}, {LAMNLIN_CARD} and {SPEC_ERR_CARD} are left out"
    if not dispMapPath:
        log.warning(f"The spectrum has no dispersion map table to read the wavelength accuracy from, {leftOut}.")
        return None
    try:
        header = fits.getheader(dispMapPath)
    except (OSError, ValueError) as error:
        log.warning(f"The dispersion map {dispMapPath} cannot be read ({error}), {leftOut}.")
        return None
    if WAVE_RMS_QC_CARD not in header or WAVE_NLINES_QC_CARD not in header:
        log.warning(
            f"The dispersion map {dispMapPath} has no {WAVE_RMS_QC_CARD} and {WAVE_NLINES_QC_CARD} QC cards, "
            f"which the map holds only when it was made by a soxspipe version that writes them, {leftOut}."
        )
        return None
    return header[WAVE_RMS_QC_CARD], header[WAVE_NLINES_QC_CARD]


def wavelength_accuracy_cards(dispMapPath: Any, log: Any) -> Cards:
    """*the* ``LAMRMS``, ``LAMNLIN`` *and* ``SPEC_ERR`` *cards of a dispersion map*

    **Key Arguments:**

    - ``dispMapPath`` -- the path to the ``DISP_TAB`` of the science frames
    - ``log`` -- logger

    **Return:**

    - ``cards`` -- ``(keyword, value, comment)`` tuples. Empty, with one warning, when the map is missing, cannot be
      read, has no ``WAVE RMS`` and ``WAVE NLINES`` QCs, or holds a line count below 1 or a value that is not finite

    **Usage:**

    ```python
    from soxspipe.commonutils.wavelength_keywords import wavelength_accuracy_cards
    for keyword, value, comment in wavelength_accuracy_cards(self.dispMap, self.log):
        set_header_card(header, keyword, value, comment)
    ```
    """
    qcs = _read_wave_qcs(dispMapPath, log)
    if qcs is None:
        return []
    try:
        rmsNm, lineCount = float(qcs[0]), int(round(float(qcs[1])))
    except (TypeError, ValueError, OverflowError):
        rmsNm, lineCount = math.nan, 0
    if not math.isfinite(rmsNm) or rmsNm < 0 or lineCount < 1:
        log.warning(
            f"The dispersion map {dispMapPath} gives a {WAVE_RMS_QC_CARD} of {qcs[0]!r} and a {WAVE_NLINES_QC_CARD} "
            f"of {qcs[1]!r}, so {LAMRMS_CARD}, {LAMNLIN_CARD} and {SPEC_ERR_CARD} are left out."
        )
        return []
    return [
        (LAMRMS_CARD, rmsNm, "[nm] RMS of the wavelength solution fit"),
        (LAMNLIN_CARD, lineCount, "Number of lines in the wavelength solution"),
        (SPEC_ERR_CARD, rmsNm / math.sqrt(lineCount), "[nm] Statistical error of wavelengths"),
    ]


def _slit_key(slitName: Any) -> str:
    """*the name of a slit as the* ``spec-res`` *table keys it: stripped, without the NIR JH filter suffix*"""
    return "" if slitName is None else str(slitName).strip().removesuffix(JH_SUFFIX)


def spec_res_cards(arm: str, slitName: Any, settings: Any, log: Any) -> Cards:
    """*the* ``SPEC_RES`` *card of an arm and slit, from the nominal* ``spec-res`` *table*

    **Key Arguments:**

    - ``arm`` -- the arm of the spectrum, ``VIS`` or ``NIR``
    - ``slitName`` -- the slit name of the science frames, for example ``SLIT1.0``. A trailing ``JH`` is ignored
    - ``settings`` -- the settings dictionary, which selects the detector parameters
    - ``log`` -- logger

    **Return:**

    - ``cards`` -- one ``SPEC_RES`` tuple. Empty, with one warning that says the arm fails, when the table has no row
      for the slit: a ``BLANK`` or ``UNDEFINED`` slit, a name that is not in the table, or no slit card at all

    **Usage:**

    ```python
    from soxspipe.commonutils.wavelength_keywords import spec_res_cards
    cards = spec_res_cards("VIS", "SLIT1.0", self.settings, self.log)
    ```
    """
    from soxspipe.commonutils.detector_lookup import detector_lookup

    slit = _slit_key(slitName)
    try:
        value = int(detector_lookup(log=log, settings=settings).get(arm)["spec-res"][slit])
    except (LookupError, AttributeError, TypeError, ValueError) as error:
        log.warning(
            f"The {arm} slit {slitName!r} has no usable row in the spec-res table ({type(error).__name__}: {error}), "
            f"so the spectrum has no {SPEC_RES_CARD} and the arm fails."
        )
        return []
    return [(SPEC_RES_CARD, value, "Nominal resolving power from slit arcs")]


def _largest_sky_shift_error_nm(qcTable: Any) -> float | None:
    """*the largest finite ``SKY SHIFT O<n> ERR`` value in a QC table, or None when it has none*"""
    import pandas as pd

    if not isinstance(qcTable, pd.DataFrame) or not {"qc_name", "qc_value"} <= set(qcTable.columns):
        return None
    rows = qcTable.loc[qcTable["qc_name"].astype(str).str.fullmatch(SKY_SHIFT_ERROR_QC), "qc_value"]
    values = pd.to_numeric(rows, errors="coerce")
    values = values[values.map(math.isfinite)]
    return float(values.max()) if len(values) else None


def spec_sye_cards(qcTable: Any, arm: str, log: Any) -> Cards:
    """*the* ``SPEC_SYE`` *card: the largest* ``SKY SHIFT O<n> ERR`` *QC of the spectrum*

    **Key Arguments:**

    - ``qcTable`` -- the QC table of the recipe, which holds the per-order ``SKY SHIFT O<n> ERR`` rows
    - ``arm`` -- the arm of the spectrum, named in the warning
    - ``log`` -- logger

    **Return:**

    - ``cards`` -- one ``SPEC_SYE`` tuple (nm). Empty, with one warning, when the table holds no finite error row: the
      arm writes none (NIR, until its error formulas exist), or there is no QC table

    **Usage:**

    ```python
    from soxspipe.commonutils.wavelength_keywords import spec_sye_cards
    cards = spec_sye_cards(self.qc, self.arm, self.log)
    ```
    """
    largestNm = _largest_sky_shift_error_nm(qcTable)
    if largestNm is None:
        log.warning(
            f"The {arm} spectrum has no SKY SHIFT O<n> ERR QC to take {SPEC_SYE_CARD} from, so {SPEC_SYE_CARD} is "
            "left out."
        )
        return []
    return [(SPEC_SYE_CARD, largestNm, "[nm] Systematic error of wavelengths")]


def with_wavelength_keywords(
    header: Any, dispMapPath: Any, arm: str, kw: Any, settings: Any, log: Any, qc: Any = None
) -> Any:
    """*give a header the wavelength-calibration and resolution cards, on a copy*

    The cards an earlier stamp left go first, so a card that cannot be computed now does not survive from before. A
    SOXS header gets the cards. Any other instrument gets an unchanged copy and no warning.

    **Key Arguments:**

    - ``header`` -- the header to stamp. It is not changed
    - ``dispMapPath`` -- the path to the ``DISP_TAB`` of the science frames
    - ``arm`` -- the arm of the spectrum
    - ``kw`` -- the keyword lookup (``keyword_lookup(...).get``), which names the slit card of the arm
    - ``settings`` -- the settings dictionary
    - ``log`` -- logger
    - ``qc`` -- the QC table of the recipe, which gives ``SPEC_SYE``. Default *None*, which gives no ``SPEC_SYE``

    **Return:**

    - ``stamped`` -- a copy of ``header`` with ``SPEC_RES``, ``LAMRMS``, ``LAMNLIN``, ``SPEC_ERR`` and ``SPEC_SYE``,
      less any that could not be computed

    **Usage:**

    ```python
    from soxspipe.commonutils.wavelength_keywords import with_wavelength_keywords
    header = with_wavelength_keywords(header, self.dispMap, self.arm, self.kw, self.settings, self.log, self.qc)
    ```
    """
    from soxspipe.commonutils.fits_headers import set_header_card

    stamped = copy.deepcopy(header)
    if not _is_soxs(settings):
        return stamped
    for keyword in (SPEC_RES_CARD, SPEC_SYE_CARD, *ACCURACY_CARDS):
        stamped.pop(keyword, None)
    cards = spec_res_cards(arm, header.get(kw(f"SLIT_{arm}".upper())), settings, log)
    cards += wavelength_accuracy_cards(dispMapPath, log)
    cards += spec_sye_cards(qc, arm, log)
    for keyword, value, comment in cards:
        set_header_card(stamped, keyword, value, comment)
    return stamped


def with_merged_spectrum_cards(
    header: Any, entries: Any, dispMapPath: Any, arm: str, kw: Any, settings: Any, log: Any, qc: Any = None
) -> Any:
    """*give a header every card a merged spectrum takes from its recipe, on a copy*

    The frame-accounting and target position cards of the on-target list, then the wavelength-calibration and
    resolution cards. The three places that write a merged spectrum header use this one composition.

    **Key Arguments:**

    - ``header`` -- the header to stamp. It is not changed
    - ``entries`` -- the on-target list of the recipe, or None
    - ``dispMapPath`` -- the path to the ``DISP_TAB`` of the science frames
    - ``arm`` -- the arm of the spectrum
    - ``kw`` -- the keyword lookup (``keyword_lookup(...).get``)
    - ``settings`` -- the settings dictionary
    - ``log`` -- logger
    - ``qc`` -- the QC table of the recipe, which gives ``SPEC_SYE``. Default *None*, which gives no ``SPEC_SYE``

    **Return:**

    - ``stamped`` -- a stamped copy of ``header``

    **Usage:**

    ```python
    from soxspipe.commonutils.wavelength_keywords import with_merged_spectrum_cards
    header = with_merged_spectrum_cards(
        header, self.onTargetFrames, self.dispMap, self.arm, self.kw, self.settings, self.log, self.qc
    )
    ```
    """
    from soxspipe.commonutils.on_target_frames import with_frame_accounting

    stamped = with_frame_accounting(header, entries, log)
    return with_wavelength_keywords(stamped, dispMapPath, arm, kw, settings, log, qc)
