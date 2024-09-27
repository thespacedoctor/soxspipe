#!/usr/bin/env python
# encoding: utf-8
"""
*Given a FITS object, use the SOXS file-naming scheme to return a filename to be used to save the FITS object to disk*

Author
: David Young

Date Created
: March  9, 2021
"""
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


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
        "-", ".").replace(":", ".")
    obid = frame.header[kw("OBS_ID")]
    arm = frame.header[kw("SEQ_ARM")].lower()
    # x = int(dp["binning"][1])
    # y = int(dp["binning"][0])
    if frame.wcs:
        x = int(frame.wcs.to_header(relax=True)["CDELT1"])
        y = int(frame.wcs.to_header(relax=True)["CDELT2"])
        binning = f"_{x}x{y}"
    else:
        binning = ""

    romode = ""

    if kw("DET_READ_SPEED") in frame.header:
        if frame.header[kw("INSTRUME")].strip().upper() == "SOXS":
            romode = "_ro" + str(frame.header[kw("DET_READ_SPEED")])
        else:
            if frame.header[kw("DET_READ_SPEED")] == 1:
                romode = "_rospeed1"
            elif "100k" in frame.header[kw("DET_READ_SPEED")].lower():
                romode = "_slow"
            elif "400k" in frame.header[kw("DET_READ_SPEED")].lower():
                romode = "_fast"
            else:
                log.print(frame.header[kw("DET_READ_SPEED")])
                raise LookupError(f"Cound not parse readout mode")

    filename = f"{dateStamp}_{arm}{binning}{romode}"

    ttype = None
    obsmode = None

    # DETERMINE THE TYPE
    if kw("DPR_TYPE") not in frame.header and kw("PRO_TYPE") in frame.header:
        return None

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
    elif frame.header[kw("DPR_TYPE")].upper() == "LAMP,FMTCHK" or frame.header[kw("DPR_TYPE")].upper() == "LAMP,WAVE":
        ttype = "arc"
    elif "LAMP" in frame.header[kw("DPR_TYPE")].upper() and "ORDERDEF" in frame.header[kw("DPR_TYPE")].upper():
        ttype = "flat"
    elif "OBJECT" in frame.header[kw("DPR_TYPE")].upper() and ("STARE" in frame.header[kw("DPR_TECH")].upper() or "NODDING" in frame.header[kw("DPR_TECH")].upper()):
        object = frame.header[kw("OBJECT")].upper()
        ttype = f"object_stare_{object}".replace(" ", "_").replace("-", "_").replace("__", "_").replace("__", "_")
    elif "STD,FLUX" in frame.header[kw("DPR_TYPE")].upper() and ("STARE" in frame.header[kw("DPR_TECH")].upper() or "NODDING" in frame.header[kw("DPR_TECH")].upper()):
        object = frame.header[kw("OBJECT")].upper()
        ttype = f"std_flux_stare_{object}".replace(" ", "_").replace("-", "_").replace("__", "_").replace("__", "_")

    if ",Q" in frame.header[kw("DPR_TYPE")].upper():
        lamp = "_QLAMP"
    elif ",D" in frame.header[kw("DPR_TYPE")].upper():
        lamp = "_DLAMP"
    else:
        lamp = ""

    if ttype is None:
        print(repr(frame.header))
        print()

        print(frame.header[kw("DPR_TYPE")].lower())
        print(frame.header[kw("DPR_TECH")].lower())
        print(frame.header[kw("DPR_CATG")].lower())

        message = "Frame type can't be determined - exiting"
        log.error(message)
        raise TypeError(message)

    filename = f"{filename}{lamp}_{ttype}"

    maskSlit = None
    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,PINHOLE":
        maskSlit = "onepin"
    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,MULTI-PINHOLE":
        maskSlit = "multipin"

    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,SLIT" and ttype in ("mflat", "flat"):
        maskSlit = "slit"
    if frame.header[kw("DPR_TECH")].upper() in ("ECHELLE,SLIT,STARE", "ECHELLE,SLIT,NODDING") and ("object" in ttype or 'std_flux' in ttype):
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

    filename = filename.upper()
    filename += ".fits"

    log.debug('completed the ``filenamer`` function')
    return filename
