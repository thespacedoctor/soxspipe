#!/usr/bin/env python
# encoding: utf-8
"""
*functions used to make soxs products ESO Phase 3 compliant*

Author
: David Young

Date Created
: June 9, 2026
"""


def basic_header_scrubbing(log, settings, header):
    """
    *Clean up header keywords/values/comments to make the header ESO Phase 3 compliant*

    **Key Arguments:**

    - ``log`` -- the logger object
    - ``settings`` -- the settings dictionary
    - ``header`` -- the header to be scrubbed

    **Usage:**

    ```python
    from soxspipe.commonutils.phase3 import basic_header_scrubbing
    header = basic_header_scrubbing(log, settings, header)
    ```
    """
    log.debug("starting the ``basic_header_scrubbing`` method")
    from soxspipe.commonutils import keyword_lookup

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(log=log, settings=settings).get

    removeKw = ["DPR_TECH", "DPR_CATG", "DPR_TYPE"]
    for k in removeKw:
        try:
            header.pop(kw(k))
        except:
            pass

    if "NAXIS" in header and header["NAXIS"] != 0 and "INHERIT" in header:
        del header["INHERIT"]

    # KEYWORDS TO DELETE
    deleteKw = ["ARCFILE"]
    for k in deleteKw:
        try:
            header.pop(k)
        except:
            pass

    # KEYWORDS TO RENAME
    renameKw = {
        "RADECSYS": "RADESYS",
        "HIERARCH ESO TEL TARG EQUINOX": "EQUINOX",
    }
    for oldKw, newKw in renameKw.items():
        if oldKw in header:
            comment = header.comments[oldKw]
            header[newKw] = (header.pop(oldKw), comment)

    log.debug("completed the ``basic_header_scrubbing`` method")
    return header


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
                header[k] = v
            else:
                header[k] = (v, c)

    log.debug("completed the ``sort_keywords`` method")
    return header
