"""*the wavelength-calibration and resolution cards of a merged spectrum (DY-2358)*

``LAMRMS``, ``LAMNLIN`` and ``SPEC_ERR`` describe how well the dispersion solution fits the arc lines. They are read
from the ``WAVE RMS`` and ``WAVE NLINES`` QCs of the ``DISP_TAB`` that the science frames were reduced with: the recipe
already holds the path of that map, so the merged spectrum always reports the map its extraction used.
``SPEC_ERR = LAMRMS / sqrt(LAMNLIN)``, as the Phase 3 standard says.

``SPEC_RES`` is the resolving power through the science slit. It comes from the nominal ``spec-res`` table of the
detector parameters, by arm and slit name. It is never the pinhole value ``R PIN MEDIAN``, which is for a 0.5 arcsec
aperture and is far above the value for the slit.

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
from typing import Any

SOXS_INSTRUMENT = "soxs"

LAMRMS_CARD = "LAMRMS"
LAMNLIN_CARD = "LAMNLIN"
SPEC_ERR_CARD = "SPEC_ERR"
SPEC_RES_CARD = "SPEC_RES"
ACCURACY_CARDS = (LAMRMS_CARD, LAMNLIN_CARD, SPEC_ERR_CARD)

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

    table = detector_lookup(log=log, settings=settings).get(arm).get("spec-res", {})
    slit = _slit_key(slitName)
    if slit not in table:
        log.warning(
            f"The {arm} slit {slitName!r} has no row in the spec-res table, so the spectrum has no {SPEC_RES_CARD} "
            "and the arm fails."
        )
        return []
    return [(SPEC_RES_CARD, int(table[slit]), "Nominal resolving power from slit arcs")]


def with_wavelength_keywords(header: Any, dispMapPath: Any, arm: str, kw: Any, settings: Any, log: Any) -> Any:
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

    **Return:**

    - ``stamped`` -- a copy of ``header`` with ``SPEC_RES``, ``LAMRMS``, ``LAMNLIN`` and ``SPEC_ERR``, less any that
      could not be computed

    **Usage:**

    ```python
    from soxspipe.commonutils.wavelength_keywords import with_wavelength_keywords
    header = with_wavelength_keywords(header, self.dispMap, self.arm, self.kw, self.settings, self.log)
    ```
    """
    from soxspipe.commonutils.fits_headers import set_header_card

    stamped = copy.deepcopy(header)
    if not _is_soxs(settings):
        return stamped
    for keyword in (SPEC_RES_CARD, *ACCURACY_CARDS):
        stamped.pop(keyword, None)
    cards = spec_res_cards(arm, header.get(kw(f"SLIT_{arm}".upper())), settings, log)
    cards += wavelength_accuracy_cards(dispMapPath, log)
    for keyword, value, comment in cards:
        set_header_card(stamped, keyword, value, comment)
    return stamped
