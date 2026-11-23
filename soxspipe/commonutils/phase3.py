#!/usr/bin/env python
"""
*functions used to make soxs products ESO Phase 3 compliant*

Author
: David Young

Date Created
: June 9, 2026
"""

PER_ARM_MODE = "per-arm"
JOINED_MODE = "joined"
HEADER_MODES = (PER_ARM_MODE, JOINED_MODE)

# CARDS OF THE UV-VIS ARM THAT A JOINED PRODUCT DROPS, AS ASTROPY NAMES THEM (WITHOUT THE HIERARCH PREFIX)
JOINED_DROPPED_PREFIXES = ("ESO DET ", "ESO INS ", "ESO SEQ ", "ESO QC ")

# THE FUNCTION NAMES CCDPROC WRITES AS CARDS, AND THE SHORT NAMES IT WRITES THEIR ARGUMENTS UNDER
CCDPROC_CARDS = (
    "CCD_PROCESS",
    "CCDPROC",
    "COSMICRAY_LACOSMIC",
    "CREATE_DEVIATION",
    "CREATVAR",
    "FLAT_CORRECT",
    "FLATCOR",
    "GAIN_CORRECT",
    "GAINCOR",
    "REBIN",
    "SUBBIAS",
    "SUBDARK",
    "SUBOSCAN",
    "SUBTRACT_BIAS",
    "SUBTRACT_DARK",
    "SUBTRACT_OVERSCAN",
    "TRANSFORM_IMAGE",
    "TRANIM",
    "TRIM_IMAGE",
    "TRIMIM",
    "WCS_PROJECT",
    "WCSPROJ",
)

# CARDS THAT BELONG ON AN IMAGE AND ARE REJECTED ON A 1D TABLE
TABLE_DROPPED_CARDS = ("BUNIT", "SXSPRE", *CCDPROC_CARDS)

# ESO'S LIMIT ON A STRING VALUE, AND THE LENGTH OF ONE FITS CARD. A LONGER CARD IS WRITTEN AS CONTINUE CARDS
MAX_VALUE_LENGTH = 68
CARD_LENGTH = 80


def basic_header_scrubbing(log, settings, header, mode=PER_ARM_MODE, isTable=False):
    """
    *Clean up header keywords/values/comments to make the header ESO Phase 3 compliant*

    The per-arm mode keeps every ``HIERARCH ESO`` card except ``DPR``, and the ``ESO QC`` cards. The joined mode is for
    a header taken from the UV-VIS arm of a joined product: it also drops the ``DET``, ``INS``, ``SEQ`` and ``QC``
    cards and keeps the ``OBS``, ``TPL`` and ``TEL`` cards.

    **Key Arguments:**

    - ``log`` -- the logger object
    - ``settings`` -- the settings dictionary
    - ``header`` -- the header to be scrubbed
    - ``mode`` -- ``PER_ARM_MODE`` or ``JOINED_MODE``. Default *PER_ARM_MODE*
    - ``isTable`` -- the header is for a 1D table product, so it also loses ``BUNIT``, ``SXSPRE`` and the ccdproc
      cards. An image header keeps them. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.phase3 import basic_header_scrubbing
    header = basic_header_scrubbing(log, settings, header)
    ```
    """
    log.debug("starting the ``basic_header_scrubbing`` method")
    from soxspipe.commonutils import keyword_lookup
    from soxspipe.commonutils.fits_headers import set_header_card

    if mode not in HEADER_MODES:
        raise ValueError(f"basic_header_scrubbing: mode must be one of {HEADER_MODES}, not {mode!r}")

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(log=log, settings=settings).get

    removeKw = ["DPR_TECH", "DPR_CATG", "DPR_TYPE"]
    for k in removeKw:
        try:
            header.pop(kw(k))
        except KeyError as e:
            log.debug(f"basic_header_scrubbing: `header.pop(kw(k))` failed, continuing: {e}")

    if "NAXIS" in header and header["NAXIS"] != 0 and "INHERIT" in header:
        del header["INHERIT"]

    # KEYWORDS TO DELETE
    deleteKw = ["ARCFILE"]
    for k in deleteKw:
        try:
            header.pop(k)
        except KeyError as e:
            log.debug(f"basic_header_scrubbing: `header.pop(k)` failed, continuing: {e}")

    # KEYWORDS TO RENAME
    renameKw = {
        "RADECSYS": "RADESYS",
        "HIERARCH ESO TEL TARG EQUINOX": "EQUINOX",
    }
    for oldKw, newKw in renameKw.items():
        if oldKw in header:
            comment = header.comments[oldKw]
            set_header_card(header, newKw, header.pop(oldKw), comment)

    if mode == JOINED_MODE:
        _drop_cards(header, [k for k in header if k.startswith(JOINED_DROPPED_PREFIXES)])
    if isTable:
        _drop_cards(header, TABLE_DROPPED_CARDS)

    log.debug("completed the ``basic_header_scrubbing`` method")
    return header


def _drop_cards(header, keywords):
    """*remove every card of each keyword from a header, whether or not it is there*"""
    for keyword in keywords:
        header.remove(keyword, ignore_missing=True, remove_all=True)


def sort_keywords(log, header):
    """
    *Neatly sort header keywords/values/comments*

    **Key Arguments:**

    - ``header`` -- the header to be sorted

    **Usage:**

    ```python
    from soxspipe.commonutils.phase3 import sort_keywords
    header = sort_keywords(log=log, header=frame.header)
    ```
    """
    from soxspipe.commonutils.fits_headers import set_header_card

    log.debug("starting the ``sort_keywords`` method")

    # NEATLY SORT KEYWORDS
    keywords = [k for k in header if len(k)]
    values = [header[k] for k in header if len(k)]
    comments = [header.comments[k] for k in header if len(k)]
    keywords, values, comments = zip(*sorted(zip(keywords, values, comments)))
    if "COMMENT" not in keywords and "HISTORY" not in keywords:
        header.clear()
        for k, v, c in zip(keywords, values, comments):
            if k == "COMMENT":
                set_header_card(header, k, v)
            else:
                set_header_card(header, k, v, c)

    log.debug("completed the ``sort_keywords`` method")
    return header


def add_qc_to_header(header, qc):
    """*write the QC rows flagged ``to_header`` into a FITS header as ``ESO QC <NAME>`` cards*

    FITS headers cannot hold NaN or infinite floats, so a row with a non-finite value is left out.

    **Key Arguments:**

    - ``header`` -- the astropy FITS header to update in place
    - ``qc`` -- the QC table, with ``qc_name``, ``qc_value``, ``qc_comment`` and ``to_header`` columns

    **Return:**

    - ``header`` -- the updated header
    """
    import numpy as np

    from soxspipe.commonutils.fits_headers import set_header_card

    for n, v, c, h in zip(
        qc["qc_name"].values,
        qc["qc_value"].values,
        qc["qc_comment"].values,
        qc["to_header"].values,
        strict=True,
    ):
        isMissing = isinstance(v, (float, np.floating)) and not np.isfinite(v)
        if h and not isMissing:
            set_header_card(header, f"ESO QC {n}".upper(), v, c)
    return header


def pipe_id():
    """*return the pipeline identification string, as the ``ESO PRO REC1 PIPE ID`` card holds it*"""
    from soxspipe.__version__ import __version__

    return f"soxspipe/v{__version__}"


def validate_header_cards(header, hduName):
    """*raise an error when a header holds a card that ESO rejects*

    ESO rejects a string value over 68 characters, and any ``CONTINUE`` card. A card that astropy would write as
    ``CONTINUE`` cards is the same fault. Commentary cards are not values and are not checked.

    **Key Arguments:**

    - ``header`` -- the header to check
    - ``hduName`` -- the name of the HDU, for the error message

    **Usage:**

    ```python
    from soxspipe.commonutils.phase3 import validate_header_cards
    validate_header_cards(header, "primary")
    ```
    """
    for card in header.cards:
        keyword = card.keyword
        if keyword in ("COMMENT", "HISTORY", ""):
            continue
        if keyword == "CONTINUE":
            raise ValueError(f"The {hduName} header holds a CONTINUE card, which ESO rejects")
        if isinstance(card.value, str) and len(card.value) > MAX_VALUE_LENGTH:
            raise ValueError(
                f"The {hduName} header value of {keyword} is {len(card.value)} characters, over the "
                f"{MAX_VALUE_LENGTH} that ESO accepts"
            )
        if len(card.image) > CARD_LENGTH:
            raise ValueError(
                f"The {hduName} header card {keyword} is too long for one card and would need CONTINUE cards"
            )


def _set_identification_cards(header, filePath):
    """*set the fixed identification cards of a product and remove the primary header's ``EXTNAME``*"""
    import os

    from soxspipe.commonutils.fits_headers import set_header_card

    header.remove("EXTNAME", ignore_missing=True, remove_all=True)
    set_header_card(header, "ORIGFILE", os.path.basename(str(filePath)), "File name of this product")
    set_header_card(header, "ORIGIN", "ESO", "Processing system")
    set_header_card(header, "PROCSOFT", pipe_id(), "Data reduction software")
    set_header_card(header, "REFERENC", "", "Reference publication")


def write_fits_table_to_disk(log, settings, header, tables, filePath, qc=None, mode=PER_ARM_MODE, extname=None):
    """
    *Write a fits table to disk with a header that is ESO Phase 3 compliant*

    The primary header gets ``ORIGFILE`` (the file's own name), ``ORIGIN``, ``PROCSOFT`` and ``REFERENC``, and no
    ``EXTNAME``. Every HDU loses ``BUNIT``, ``SXSPRE`` and the ccdproc cards. The write raises a ``ValueError``, and no
    file is written, when a header holds a value over 68 characters or a ``CONTINUE`` card.

    **Key Arguments:**

    - ``log`` -- the logger object
    - ``settings`` -- the settings dictionary
    - ``header`` -- the header to be written to disk
    - ``tables`` -- the tables to be written to disk, one extension each
    - ``filePath`` -- the path to the output file
    - ``qc`` -- a dataframe containing QC metrics to be added to the header (default: None)
    - ``mode`` -- ``PER_ARM_MODE`` or ``JOINED_MODE``, the header carry-over rules (default: PER_ARM_MODE)
    - ``extname`` -- the ``EXTNAME`` of the first table, e.g. ``SPECTRUM`` (default: None, no ``EXTNAME``)

    **Usage:**

    ```python
    from soxspipe.commonutils.phase3 import write_fits_table_to_disk
    write_fits_table_to_disk(log=log, settings=settings, header=header, tables=[table], filePath=filePath, qc=qc)
    ```
    """
    log.debug("starting the ``write_fits_table_to_disk`` method")
    from astropy.io import fits
    from astropy.table import Table

    from soxspipe.commonutils.atomic_fits import write_hdulist_atomically

    if qc is not None:
        # ADD QC METRICS TO HEADER
        add_qc_to_header(header, qc)

    # CONVERT TO FITS BINARY TABLE HDU
    tableExensions = []
    for table in tables:
        if not isinstance(table, Table):
            table = Table.from_pandas(table)
        hdu = fits.table_to_hdu(table)
        _drop_cards(hdu.header, TABLE_DROPPED_CARDS)
        tableExensions.append(hdu)
    if extname and tableExensions:
        tableExensions[0].name = extname

    header = basic_header_scrubbing(log, settings, header, mode=mode, isTable=True)
    _set_identification_cards(header, filePath)
    header = sort_keywords(log, header)

    validate_header_cards(header, "primary")
    for index, hdu in enumerate(tableExensions, start=1):
        validate_header_cards(hdu.header, f"extension {index}")

    priHDU = fits.PrimaryHDU(header=header)
    hduList = fits.HDUList([priHDU] + tableExensions)
    hduList.verify("fix")
    write_hdulist_atomically(hduList, filePath, checksum=True, overwrite=True)

    log.debug("completed the ``write_fits_table_to_disk`` method")
