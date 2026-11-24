#!/usr/bin/env python
"""
*the keyword rules of the Phase 3 standard for a spectrum, and the checks of the headers*

Table 8 of the standard gives each keyword a level. ``KEYWORD_RULES`` holds the `SPECTRUM` column, and the checks here
report a missing keyword, a value of the wrong kind, a value the standard does not allow, the identification and
hygiene cards that ``write_fits_table_to_disk`` sets or drops, and the keywords that must agree between the headers.

Author
: David Young

Date Created
: November 24, 2026
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

from soxspipe.commonutils.phase3 import (
    ORIGIN_VALUE,
    PIPE_ID_PREFIX,
    RENAMED_ESO_CARDS,
    SCRUBBED_ESO_CARDS,
    SPECTRUM_EXTNAME,
    TABLE_DROPPED_CARDS,
    header_card_faults,
)
from soxspipe.commonutils.phase3_check_framework import (
    EXTENSION,
    PRIMARY,
    Finding,
    Phase3Product,
    fail,
    is_phase3_product,
    register_check,
    warn,
)

if TYPE_CHECKING:
    from astropy.io.fits import Header

# THE LEVELS OF TABLE 8 OF THE STANDARD, AND THE WORD EACH HAS IN A MESSAGE
MANDATORY = "M"
RECOMMENDED = "R"
OPTIONAL = "O"
LEVEL_WORDS = {MANDATORY: "mandatory", RECOMMENDED: "recommended", OPTIONAL: "optional"}

# WHERE A RULE COMES FROM: THE STANDARD'S TABLE 8, OR A RULE OF OUR OWN (THE SPEC, DY-2393)
TABLE_8 = "Table 8"
HOUSE_RULE = "house rule, DY-2393"

# HEADER TIMES ARE ROUNDED, SO TEXPTIME MAY EXCEED TELAPSE BY THIS MANY SECONDS WITHOUT BEING A FAULT
TIME_TOLERANCE_S = 1e-3

# THE MULTI-PROGRAMME VALUE OF PROG_ID, WHICH NEEDS A PROGID1 CARD TO LIST THE PROGRAMMES
MULTI_PROGRAMME = "MULTI"


@dataclass(frozen=True)
class KeywordRule:
    """*the standard's rule for one header keyword*

    **Key Arguments:**

    - ``keyword`` -- the keyword, or its stem for an ``indexed`` keyword such as ``TCOMM``
    - ``level`` -- ``MANDATORY``, ``RECOMMENDED`` or ``OPTIONAL``
    - ``hdu`` -- ``PRIMARY`` or ``EXTENSION``, the header the keyword goes in
    - ``kind`` -- the Python type of the value: ``str``, ``float``, ``int`` or ``bool``
    - ``allowEmpty`` -- a string value may be empty
    - ``source`` -- ``TABLE_8`` or ``HOUSE_RULE``, named in the message of a rule that is not from Table 8
    - ``indexed`` -- the rule applies to the keyword with each column index, ``TCOMM1`` to ``TCOMM<TFIELDS>``
    """

    keyword: str
    level: str
    hdu: str
    kind: type
    allowEmpty: bool = False
    source: str = TABLE_8
    indexed: bool = False

    def names(self, header: Header) -> list[str]:
        """*return the keywords this rule covers in a header: its own, or one per column when it is indexed*

        **Key Arguments:**

        - ``header`` -- the header the rule applies to

        **Return:**

        - ``names`` -- the keyword names

        **Usage:**

        ```python
        names = rule.names(product.extension)
        ```
        """
        if not self.indexed:
            return [self.keyword]
        return [f"{self.keyword}{index}" for index in range(1, header.get("TFIELDS", 0) + 1)]

    def label(self, name: str) -> str:
        """*return the keyword for a message, with the source of the rule when it is not Table 8*

        **Key Arguments:**

        - ``name`` -- the keyword name from ``names``

        **Return:**

        - ``label`` -- the name, with ``(house rule, DY-2393)`` after it for a rule of our own

        **Usage:**

        ```python
        text = f"keyword {rule.label(name)} is missing"
        ```
        """
        return name if self.source == TABLE_8 else f"{name} ({self.source})"


def _rules(
    keywords: str, level: str, hdu: str, kind: type, allowEmpty: bool = False, **options: object
) -> list[KeywordRule]:
    return [KeywordRule(keyword, level, hdu, kind, allowEmpty, **options) for keyword in keywords.split()]


def _both(keywords: str, kind: type) -> list[KeywordRule]:
    return _rules(keywords, MANDATORY, PRIMARY, kind) + _rules(keywords, MANDATORY, EXTENSION, kind)


# TABLE 8 OF THE STANDARD, THE `SPECTRUM` COLUMN, WITH THE EXTNAME RULE OF THE SPEC. A KEYWORD THAT TABLE 8 MARKS AS
# APPLICABLE ONLY UNDER A CONDITION (EQUINOX, TIMESYS, ASSONi) HAS ITS OWN CHECK
KEYWORD_RULES = tuple(
    _rules(
        "PRODCATG ORIGIN TELESCOP INSTRUME RADESYS OBSTECH FLUXCAL PROCSOFT PROG_ID PROV1 SPECSYS",
        MANDATORY,
        PRIMARY,
        str,
    )
    + _rules("REFERENC", MANDATORY, PRIMARY, str, allowEmpty=True)
    + _rules(
        "EXPTIME TEXPTIME MJD-OBS MJD-END FLUXERR WAVELMIN WAVELMAX SPEC_BIN SNR SPEC_RES", MANDATORY, PRIMARY, float
    )
    + _rules("OBID1 NCOMBINE", MANDATORY, PRIMARY, int)
    + _rules("CONTNORM TOT_FLUX", MANDATORY, PRIMARY, bool)
    + _rules("SPEC_ERR SPEC_SYE RA_ERR DEC_ERR", RECOMMENDED, PRIMARY, float)
    + _rules("LAMRMS", OPTIONAL, PRIMARY, float)
    + _rules("LAMNLIN", OPTIONAL, PRIMARY, int)
    + _both("OBJECT", str)
    + _both("RA DEC", float)
    + _both("CHECKSUM DATASUM", str)
    + _rules("VOCLASS VOPUB TITLE", MANDATORY, EXTENSION, str)
    + _rules("EXTNAME", MANDATORY, EXTENSION, str, source=HOUSE_RULE)
    + _rules("NELEM", MANDATORY, EXTENSION, int)
    + _rules("APERTURE TELAPSE TMID SPEC_VAL SPEC_BW", MANDATORY, EXTENSION, float)
    + _rules("TCOMM", OPTIONAL, EXTENSION, str, indexed=True)
)


def _has_kind(value: object, kind: type) -> bool:
    """*return True when a header value is of the kind the standard asks for, a bool never passing as a number*"""
    if kind is bool:
        return isinstance(value, bool)
    if isinstance(value, bool):
        return False
    if kind is float:
        return isinstance(value, (int, float)) and math.isfinite(value)
    return isinstance(value, kind)


@register_check
def check_mandatory_and_recommended_keywords(product: Phase3Product) -> list[Finding]:
    """*report each keyword of ``KEYWORD_RULES`` that is missing, or a mandatory string that is empty*

    A missing mandatory keyword is a failure. A missing recommended or optional keyword is a warning.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_mandatory_and_recommended_keywords(product)
    ```
    """
    findings = []
    for rule in KEYWORD_RULES:
        header = product.header(rule.hdu)
        if header is None:
            continue
        for name in rule.names(header):
            where = f"the {rule.hdu} header"
            if name not in header:
                message = f"{LEVEL_WORDS[rule.level]} keyword {rule.label(name)} is missing from {where}"
                findings.append(fail(message) if rule.level == MANDATORY else warn(message))
            elif rule.kind is str and not rule.allowEmpty and rule.level == MANDATORY and header[name] == "":
                findings.append(fail(f"mandatory keyword {rule.label(name)} in {where} is an empty string"))
    return findings


@register_check
def check_keyword_value_types(product: Phase3Product) -> list[Finding]:
    """*report each keyword of ``KEYWORD_RULES`` whose value is not of the kind the standard asks for*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_keyword_value_types(product)
    ```
    """
    findings = []
    for rule in KEYWORD_RULES:
        header = product.header(rule.hdu)
        if header is None:
            continue
        for name in rule.names(header):
            if name not in header or _has_kind(header[name], rule.kind):
                continue
            findings.append(
                fail(
                    f"{name} in the {rule.hdu} header has the type {type(header[name]).__name__}, "
                    f"it must be {rule.kind.__name__}"
                )
            )
    return findings


@dataclass(frozen=True)
class RangeRule:
    """*a range that the standard gives a numeric keyword*

    **Key Arguments:**

    - ``keyword`` -- the keyword
    - ``hdu`` -- ``PRIMARY`` or ``EXTENSION``
    - ``isValid`` -- returns True for a value in the range
    - ``text`` -- the range in words, for the message
    """

    keyword: str
    hdu: str
    isValid: Callable[[float], bool]
    text: str


# A SHORT LIST: THE RANGES SECTION 5 OF THE STANDARD STATES FOR THE SPECTRUM KEYWORDS
RANGE_RULES = (
    RangeRule("EXPTIME", PRIMARY, lambda value: value > 0, "greater than 0"),
    RangeRule("TEXPTIME", PRIMARY, lambda value: value > 0, "greater than 0"),
    RangeRule("TELAPSE", EXTENSION, lambda value: value > 0, "greater than 0"),
    RangeRule("NCOMBINE", PRIMARY, lambda value: value >= 1, "at least 1"),
    RangeRule("SNR", PRIMARY, lambda value: value >= 0, "at least 0"),
    RangeRule("WAVELMIN", PRIMARY, lambda value: value > 0, "greater than 0"),
    RangeRule("SPEC_RES", PRIMARY, lambda value: value > 0, "greater than 0"),
)


@register_check
def check_value_ranges(product: Phase3Product) -> list[Finding]:
    """*report a numeric keyword outside the range the standard gives it, and a WAVELMIN not below WAVELMAX*

    A value that is not a number is left to ``check_keyword_value_types``.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_value_ranges(product)
    ```
    """
    findings = []
    for rule in RANGE_RULES:
        header = product.header(rule.hdu)
        value = header.get(rule.keyword) if header is not None else None
        if _has_kind(value, float) and not rule.isValid(value):
            findings.append(fail(f"{rule.keyword} is {value}, outside its allowed range: it must be {rule.text}"))
    wavelmin, wavelmax = product.primary.get("WAVELMIN"), product.primary.get("WAVELMAX")
    if _has_kind(wavelmin, float) and _has_kind(wavelmax, float) and wavelmin >= wavelmax:
        findings.append(
            fail(f"WAVELMIN ({wavelmin}) is outside its allowed range: it must be below WAVELMAX ({wavelmax})")
        )
    return findings


@register_check
def check_multi_programme_ids(product: Phase3Product) -> list[Finding]:
    """*report a ``PROG_ID`` of ``MULTI`` that has no ``PROGID1`` card to list the programmes*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_multi_programme_ids(product)
    ```
    """
    if product.primary.get("PROG_ID") == MULTI_PROGRAMME and "PROGID1" not in product.primary:
        return [
            fail(f"PROG_ID is {MULTI_PROGRAMME!r}, so PROGID1 and the other PROGIDi cards must list the programmes")
        ]
    return []


# CD MATRIX CARDS ARE NOT ALLOWED IN A 1D SPECTRUM (TABLE 8). THE WRITER DOES NOT DROP THEM, SO THIS RULE EXISTS ONLY
# IN THE CHECKER. (THE BUNIT, SXSPRE AND CCDPROC CARDS ARE DROPPED BY THE WRITER AND LISTED IN `TABLE_DROPPED_CARDS`)
CD_MATRIX_PATTERN = r"CD\d+_\d+"


def _header_hygiene_findings(product: Phase3Product) -> list[Finding]:
    """*return the findings for the cards the Phase 3 writer sets or drops: the identification cards and EXTNAME*"""
    primary, extension = product.primary, product.extension
    findings = []
    name = product.path.name
    if primary.get("ORIGFILE") != name:
        findings.append(fail(f"ORIGFILE is {primary.get('ORIGFILE')!r}, it must be the file's own name {name!r}"))
    if extension is not None and "ORIGFILE" in extension:
        findings.append(fail("ORIGFILE is in the extension header, it belongs in the primary header only"))
    if "EXTNAME" in primary:
        findings.append(fail("EXTNAME is in the primary header, it belongs to the table only"))
    if extension is not None and "EXTNAME" in extension and extension["EXTNAME"] != SPECTRUM_EXTNAME:
        findings.append(
            fail(f"EXTNAME of the table is {extension['EXTNAME']!r}, it must be {SPECTRUM_EXTNAME!r} ({HOUSE_RULE})")
        )
    if "ORIGIN" in primary and primary["ORIGIN"] != ORIGIN_VALUE:
        findings.append(fail(f"ORIGIN is {primary['ORIGIN']!r}, it must be {ORIGIN_VALUE!r}"))
    procsoft = primary.get("PROCSOFT")
    if procsoft is not None and not str(procsoft).startswith(PIPE_ID_PREFIX):
        findings.append(fail(f"PROCSOFT is {procsoft!r}, it must start with {PIPE_ID_PREFIX!r}"))
    return findings


def _scrubbed_card_findings(header: Header, hduName: str) -> list[Finding]:
    """*return the findings for each ESO card that the writer removes or renames, found in a header*"""
    findings = [
        fail(f"{keyword} is in the {hduName} header, the writer removes it from a Phase 3 product")
        for keyword in SCRUBBED_ESO_CARDS
        if keyword in header
    ]
    findings.extend(
        fail(f"{keyword} is in the {hduName} header, the writer renames it to {newName}")
        for keyword, newName in RENAMED_ESO_CARDS.items()
        if keyword in header
    )
    return findings


@register_check
def check_header_hygiene(product: Phase3Product) -> list[Finding]:
    """*report the cards that the Phase 3 writer sets, drops or refuses, and that ESO rejects*

    These are the rules of ``write_fits_table_to_disk``: the identification cards, no ``EXTNAME`` in the primary
    header, and no ``BUNIT``, ``SXSPRE`` or ccdproc card in a 1D table. The checker adds the rule of no ``CD`` matrix.
    For a ``_ph3`` product only, it also fails an ``ESO DPR`` card or the ``ESO TEL TARG EQUINOX`` card that the writer
    removes or renames, a string value over 68 characters and a ``CONTINUE`` card. Every other ``HIERARCH ESO`` card is
    kept in a per-arm product, as in the writer. A joined product drops more ESO cards, which the checker cannot tell
    from the file name, so it does not check that.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_header_hygiene(product)
    ```
    """
    findings = _header_hygiene_findings(product)
    for hduName in (PRIMARY, EXTENSION):
        header = product.header(hduName)
        if header is None:
            continue
        findings.extend(
            fail(f"{keyword} is in the {hduName} header, a 1D spectrum must not carry it")
            for keyword in header
            if keyword in TABLE_DROPPED_CARDS or re.fullmatch(CD_MATRIX_PATTERN, keyword)
        )
        if is_phase3_product(product):
            findings.extend(_scrubbed_card_findings(header, hduName))
            findings.extend(fail(fault) for fault in header_card_faults(header, hduName))
    return findings


# THE VALUES THAT TABLE 8 AND SECTION 5 ALLOW. TIMESYS DEFAULTS TO UTC WHEN ABSENT
SPECTRUM_CATEGORY = "SCIENCE.SPECTRUM"
RADESYS_VALUES = ("ICRS", "FK5")
TIMESYS_VALUES = ("UTC", "TAI")
FLUXCAL_VALUES = ("ABSOLUTE", "UNCALIBRATED")
SPECSYS_VALUES = (
    "TOPOCENT",
    "GEOCENTR",
    "BARYCENT",
    "HELIOCEN",
    "LSRK",
    "LSRD",
    "GALACTOC",
    "LOCALGRP",
    "CMBDIPOL",
    "SOURCE",
)
EQUINOX_J2000 = 2000.0
UNCALIBRATED_FLUXERR = -1


def _not_allowed(header: Header, keyword: str, allowed: tuple[str, ...]) -> list[Finding]:
    """*return a finding when a keyword is in the header with a value outside ``allowed``*"""
    if keyword in header and header[keyword] not in allowed:
        return [fail(f"{keyword} is {header[keyword]!r}, it must be one of {list(allowed)}")]
    return []


@register_check
def check_standard_value_domains(product: Phase3Product) -> list[Finding]:
    """*report a keyword value that the standard does not allow, or a pair of keywords that contradict each other*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_standard_value_domains(product)
    ```
    """
    primary = product.primary
    findings = []
    findings += _not_allowed(primary, "PRODCATG", (SPECTRUM_CATEGORY,))
    findings += _not_allowed(primary, "RADESYS", RADESYS_VALUES)
    findings += _not_allowed(primary, "TIMESYS", TIMESYS_VALUES)
    findings += _not_allowed(primary, "FLUXCAL", FLUXCAL_VALUES)
    findings += _not_allowed(primary, "SPECSYS", SPECSYS_VALUES)
    # EQUINOX IS MANDATORY AND 2000.0 FOR FK5, AND TOLERATED, AT 2000.0, FOR ICRS
    radesys = primary.get("RADESYS")
    if radesys == "FK5" and "EQUINOX" not in primary:
        findings.append(fail("EQUINOX is missing, it is mandatory when RADESYS is 'FK5'"))
    if radesys in RADESYS_VALUES and primary.get("EQUINOX", EQUINOX_J2000) != EQUINOX_J2000:
        findings.append(fail(f"EQUINOX is {primary['EQUINOX']}, it must be {EQUINOX_J2000} with RADESYS {radesys!r}"))
    if (
        primary.get("FLUXCAL") == "UNCALIBRATED"
        and primary.get("FLUXERR", UNCALIBRATED_FLUXERR) != UNCALIBRATED_FLUXERR
    ):
        findings.append(
            fail(f"FLUXERR is {primary['FLUXERR']}, it must be {UNCALIBRATED_FLUXERR} when FLUXCAL is 'UNCALIBRATED'")
        )
    if primary.get("TOT_FLUX") is True and primary.get("FLUXCAL") != "ABSOLUTE":
        findings.append(fail("TOT_FLUX is True, FLUXCAL must be 'ABSOLUTE'"))
    return findings


@register_check
def check_primary_matches_extension(product: Phase3Product) -> list[Finding]:
    """*report an OBJECT, RA or DEC that differs between the primary header and the extension header*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_primary_matches_extension(product)
    ```
    """
    if product.extension is None:
        return []
    return [
        fail(
            f"{keyword} differs: {product.primary[keyword]!r} in the primary header, "
            f"{product.extension[keyword]!r} in the extension"
        )
        for keyword in ("OBJECT", "RA", "DEC")
        if keyword in product.primary
        and keyword in product.extension
        and product.primary[keyword] != product.extension[keyword]
    ]


@register_check
def check_exposure_within_elapsed_time(product: Phase3Product) -> list[Finding]:
    """*report a TEXPTIME that is longer than TELAPSE, the time from MJD-OBS to MJD-END*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine

    **Usage:**

    ```python
    findings = check_exposure_within_elapsed_time(product)
    ```
    """
    texptime = product.primary.get("TEXPTIME")
    telapse = product.extension.get("TELAPSE") if product.extension is not None else None
    if texptime is None or telapse is None:
        return []
    if texptime > telapse + TIME_TOLERANCE_S:
        return [fail(f"TEXPTIME ({texptime} s) is longer than TELAPSE ({telapse} s)")]
    return []
