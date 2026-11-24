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

from dataclasses import dataclass

from soxspipe.commonutils.phase3_check_framework import (
    EXTENSION,
    PRIMARY,
    fail,
    is_phase3_product,
    register_check,
    warn,
)

# THE LEVELS OF TABLE 8 OF THE STANDARD
MANDATORY = "M"
RECOMMENDED = "R"
OPTIONAL = "O"

# HEADER TIMES ARE ROUNDED, SO TEXPTIME MAY EXCEED TELAPSE BY THIS MANY SECONDS WITHOUT BEING A FAULT
TIME_TOLERANCE_S = 1e-3


@dataclass(frozen=True)
class KeywordRule:
    """*the standard's rule for one header keyword*

    **Key Arguments:**

    - ``keyword`` -- the keyword
    - ``level`` -- ``MANDATORY``, ``RECOMMENDED`` or ``OPTIONAL``
    - ``hdu`` -- ``PRIMARY`` or ``EXTENSION``, the header the keyword goes in
    - ``kind`` -- the Python type of the value: ``str``, ``float``, ``int`` or ``bool``
    - ``allowEmpty`` -- a string value may be empty
    """

    keyword: str
    level: str
    hdu: str
    kind: type
    allowEmpty: bool = False


def _rules(keywords, level, hdu, kind, allowEmpty=False):
    return [KeywordRule(keyword, level, hdu, kind, allowEmpty) for keyword in keywords.split()]


def _both(keywords, kind):
    return _rules(keywords, MANDATORY, PRIMARY, kind) + _rules(keywords, MANDATORY, EXTENSION, kind)


# TABLE 8 OF THE STANDARD, THE `SPECTRUM` COLUMN. A KEYWORD THAT TABLE 8 MARKS AS APPLICABLE ONLY UNDER A CONDITION
# (EQUINOX, TIMESYS, ASSONi) HAS ITS OWN CHECK
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
    + _rules("VOCLASS VOPUB TITLE EXTNAME", MANDATORY, EXTENSION, str)
    + _rules("NELEM", MANDATORY, EXTENSION, int)
    + _rules("APERTURE TELAPSE TMID SPEC_VAL SPEC_BW", MANDATORY, EXTENSION, float)
)


def _has_kind(value, kind):
    """*return True when a header value is of the kind the standard asks for, a bool never passing as a number*"""
    import math

    if kind is bool:
        return isinstance(value, bool)
    if isinstance(value, bool):
        return False
    if kind is float:
        return isinstance(value, (int, float)) and math.isfinite(value)
    return isinstance(value, kind)


@register_check
def check_mandatory_and_recommended_keywords(product):
    """*report each keyword of ``KEYWORD_RULES`` that is missing, or a mandatory string that is empty*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
    """
    findings = []
    for rule in KEYWORD_RULES:
        header = product.header(rule.hdu)
        if header is None:
            continue
        where = f"the {rule.hdu} header"
        if rule.keyword not in header:
            if rule.level == MANDATORY:
                findings.append(fail(f"mandatory keyword {rule.keyword} is missing from {where}"))
            else:
                level = "recommended" if rule.level == RECOMMENDED else "optional"
                findings.append(warn(f"{level} keyword {rule.keyword} is missing from {where}"))
        elif rule.kind is str and not rule.allowEmpty and rule.level == MANDATORY and header[rule.keyword] == "":
            findings.append(fail(f"mandatory keyword {rule.keyword} in {where} is an empty string"))
    return findings


@register_check
def check_keyword_value_types(product):
    """*report each keyword of ``KEYWORD_RULES`` whose value is not of the kind the standard asks for*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
    """
    findings = []
    for rule in KEYWORD_RULES:
        header = product.header(rule.hdu)
        if header is None or rule.keyword not in header:
            continue
        value = header[rule.keyword]
        if not _has_kind(value, rule.kind):
            findings.append(
                fail(
                    f"{rule.keyword} in the {rule.hdu} header has the type {type(value).__name__}, "
                    f"it must be {rule.kind.__name__}"
                )
            )
    return findings


# THE CARDS THAT THE PHASE 3 TABLE WRITER DROPS (BUNIT, SXSPRE AND THE CCDPROC CARDS) COME FROM THE WRITER, SO THE
# WRITER AND THE CHECKER CANNOT DISAGREE ABOUT THEM
CD_MATRIX_PATTERN = r"CD\d+_\d+"


def _header_hygiene_findings(product):
    """*return the findings for the cards the Phase 3 writer sets or drops: the identification cards and EXTNAME*"""
    from soxspipe.commonutils.phase3 import ORIGIN_VALUE, PIPE_ID_PREFIX, SPECTRUM_EXTNAME

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
        findings.append(fail(f"EXTNAME of the table is {extension['EXTNAME']!r}, it must be {SPECTRUM_EXTNAME!r}"))
    if "ORIGIN" in primary and primary["ORIGIN"] != ORIGIN_VALUE:
        findings.append(fail(f"ORIGIN is {primary['ORIGIN']!r}, it must be {ORIGIN_VALUE!r}"))
    procsoft = primary.get("PROCSOFT")
    if procsoft is not None and not str(procsoft).startswith(PIPE_ID_PREFIX):
        findings.append(fail(f"PROCSOFT is {procsoft!r}, it must start with {PIPE_ID_PREFIX!r}"))
    return findings


@register_check
def check_header_hygiene(product):
    """*report the cards that the Phase 3 writer sets, drops or refuses, and that ESO rejects*

    These are the rules of ``write_fits_table_to_disk``: the identification cards, no ``EXTNAME`` in the primary
    header, no ``BUNIT``, ``SXSPRE`` or ccdproc card and no ``CD`` matrix in a 1D table, and, for a ``_ph3`` product
    only, no string value over 68 characters and no ``CONTINUE`` card.

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
    """
    import re

    from soxspipe.commonutils.phase3 import TABLE_DROPPED_CARDS, header_card_faults

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


def _not_allowed(header, keyword, allowed):
    """*return a finding when a keyword is in the header with a value outside ``allowed``*"""
    if keyword in header and header[keyword] not in allowed:
        return [fail(f"{keyword} is {header[keyword]!r}, it must be one of {list(allowed)}")]
    return []


@register_check
def check_standard_value_domains(product):
    """*report a keyword value that the standard does not allow, or a pair of keywords that contradict each other*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
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
def check_primary_matches_extension(product):
    """*report an OBJECT, RA or DEC that differs between the primary header and the extension header*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
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
def check_exposure_within_elapsed_time(product):
    """*report a TEXPTIME that is longer than TELAPSE, the time from MJD-OBS to MJD-END*

    **Key Arguments:**

    - ``product`` -- the ``Phase3Product`` to check

    **Return:**

    - ``findings`` -- a list of ``Finding`` records, empty when the product is fine
    """
    texptime = product.primary.get("TEXPTIME")
    telapse = product.extension.get("TELAPSE") if product.extension is not None else None
    if texptime is None or telapse is None:
        return []
    if texptime > telapse + TIME_TOLERANCE_S:
        return [fail(f"TEXPTIME ({texptime} s) is longer than TELAPSE ({telapse} s)")]
    return []
