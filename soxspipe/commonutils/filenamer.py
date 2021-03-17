#!/usr/bin/env python
# encoding: utf-8
"""
*Given a FITS object, use the SOXS file-naming scheme to return a filename to be used to save the FITS object to disk*

:Author:
    David Young

:Date Created:
    March  9, 2021
"""
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils import detector_lookup


def filenamer(
        log,
        frame,
        keywordLookup=False,
        detectorLookup=False,
        settings=False):
    """Given a FITS object, use the SOXS file-naming scheme to return a filename to be used to save the FITS object to disk

    **Key Arguments:**
        - ``log`` -- logger
        - ``frame`` -- the CCDData object frame
        - ``keywordLookup`` -- the keyword lookup dictionary (needed if `settings` not provided). Default *False*
        - ``detectorLookup`` -- the detector parameters (needed if `settings` not provided). Default *False*
        - ``settings`` -- the soxspipe settings dictionary (needed if `keywordLookup` and `detectorLookup` not provided). Default *False*

    **Return:**
        - ``filename`` -- stanardised name to for the input frame

    ```python
    frame = CCDData.read(filepath, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
            du_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

    from soxspipe.commonutils import filenamer
    filename = filenamer(
        log=log,
        frame=frame,
        settings=settings
    )
    ```
    """
    log.debug('starting the ``filenamer`` function')

    # GENERATE A FILENAME FOR THE FRAME BASED ON THE FILENAMING
    # CONVENTION
    if keywordLookup:
        kw = keywordLookup
    else:
        kw = keyword_lookup(
            log=log,
            settings=settings
        ).get

    if detectorLookup:
        dp = detectorLookup
    else:
        arm = frame.header[kw("SEQ_ARM")]
        # DETECTOR PARAMETERS LOOKUP OBJECT
        dp = detector_lookup(
            log=log,
            settings=settings
        ).get(arm)

    dateStamp = frame.header[kw("DATE_OBS")].replace(
        "-", "").replace(":", "").split(".")[0]
    obid = frame.header[kw("OBS_ID")]
    arm = frame.header[kw("SEQ_ARM")].lower()
    x = int(dp["binning"][1])
    y = int(dp["binning"][0])
    binning = f"{x}x{y}"

    ttype = None

    obsmode = None

    # DETERMINE THE TYPE
    if frame.header[kw("DPR_TYPE")].upper() == "BIAS":
        if "SXSPRE" in frame.header:
            ttype = "mbias"
        else:
            ttype = "bias"
    elif frame.header[kw("DPR_TYPE")].upper() == "DARK":
        if "SXSPRE" in frame.header:
            ttype = "mdark"
        else:
            ttype = "dark"
    elif "LAMP" in frame.header[kw("DPR_TYPE")].upper() and "FLAT" in frame.header[kw("DPR_TYPE")].upper():
        if "SXSPRE" in frame.header:
            ttype = "mflat"
        else:
            ttype = "flat"
    elif frame.header[kw("DPR_TYPE")].upper() == "LAMP,FMTCHK":
        ttype = "arc"
    elif "LAMP" in frame.header[kw("DPR_TYPE")].upper() and "ORDERDEF" in frame.header[kw("DPR_TYPE")].upper():
        ttype = "flat"

    if ttype is None:
        print(repr(frame.header))
        print()

        print(frame.header[kw("DPR_TYPE")].lower())
        print(frame.header[kw("DPR_TECH")].lower())
        print(frame.header[kw("DPR_CATG")].lower())

        message = "Frame type can't be determined - exiting"
        log.error(message)
        raise TypeError(message)

    filename = f"{dateStamp}_{obid}_{arm}_{binning}_{ttype}"

    maskSlit = None
    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,PINHOLE":
        maskSlit = "onepin"

    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,SLIT" and ttype in ("mflat", "flat"):
        maskSlit = "slit"

    # EXTRA PARAMETERS NEEDED FOR SPECTRUM
    if frame.header[kw("DPR_TECH")].upper() != "IMAGE":

        if maskSlit is None:
            print(repr(frame.header))
            print()

            print(frame.header[kw("DPR_TYPE")].lower())
            print(frame.header[kw("DPR_TECH")].lower())
            print(frame.header[kw("DPR_CATG")].lower())
            message = "Frame mask/slit can't be determined - exiting"
            log.error(message)
            raise TypeError(message)

    if maskSlit:
        filename = f"{filename}_{maskSlit}"

    filename += ".fits"

    log.debug('completed the ``filenamer`` function')
    return filename
