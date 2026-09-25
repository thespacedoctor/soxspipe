#!/usr/bin/env python
"""
*Given a FITS object, use the SOXS file-naming scheme to return a filename to be used to save the FITS object to disk*

Author
: David Young

Date Created
: March  9, 2021
"""

import os

from soxspipe.commonutils import detector_lookup, keyword_lookup

os.environ["TERM"] = "vt100"


def _binning_fragment(frame):
    """build the binning fragment of a filename from the frame WCS

    **Key Arguments:**

    - ``frame`` -- the CCDData object frame

    **Return:**

    - ``binning`` -- ``_<x>x<y>`` from the truncated WCS pixel scales, or an empty string if the frame has no WCS
    """
    if frame.wcs:
        x = int(frame.wcs.to_header(relax=True)["CDELT1"])
        y = int(frame.wcs.to_header(relax=True)["CDELT2"])
        binning = f"_{x}x{y}"
    else:
        binning = ""
    return binning


def _readout_fragment(log, frame, kw):
    """build the read-out mode fragment of a filename

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- the CCDData object frame
    - ``kw`` -- the keyword lookup function

    **Return:**

    - ``romode`` -- the read-out fragment, or an empty string if the frame has no read-speed keyword
    """
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
                raise LookupError("Cound not parse readout mode")
    return romode


def _frame_type(frame, kw):
    """work out the frame type fragment of a filename from the DPR keywords

    **Key Arguments:**

    - ``frame`` -- the CCDData object frame
    - ``kw`` -- the keyword lookup function

    **Return:**

    - ``ttype`` -- the frame type, or *None* if the header matches no known type
    """
    ttype = None

    if frame.header[kw("DPR_TYPE")].upper() == "BIAS":
        ttype = "mbias" if "SXSPRE" in frame.header else "bias"
    elif frame.header[kw("DPR_TYPE")].upper() == "DARK":
        ttype = "mdark" if "SXSPRE" in frame.header else "dark"
    elif (
        "LAMP" in frame.header[kw("DPR_TYPE")].upper()
        and "FLAT" in frame.header[kw("DPR_TYPE")].upper()
    ):
        ttype = "mflat" if "SXSPRE" in frame.header else "flat"
    elif (
        frame.header[kw("DPR_TYPE")].upper() == "LAMP,FMTCHK"
        or frame.header[kw("DPR_TYPE")].upper() == "LAMP,WAVE"
    ):
        ttype = "arc"
    elif (
        "LAMP" in frame.header[kw("DPR_TYPE")].upper()
        and "ORDERDEF" in frame.header[kw("DPR_TYPE")].upper()
    ):
        ttype = "flat"
    elif "OBJECT" in frame.header[kw("DPR_TYPE")].upper() and (
        "STARE" in frame.header[kw("DPR_TECH")].upper()
        or "NODDING" in frame.header[kw("DPR_TECH")].upper()
    ):
        object = frame.header[kw("OBJECT")].upper()
        ttype = (
            f"object_stare_{object}".replace(" ", "_")
            .replace("-", "_")
            .replace("__", "_")
            .replace("__", "_")
        )
    elif "STD,FLUX" in frame.header[kw("DPR_TYPE")].upper() and (
        "STARE" in frame.header[kw("DPR_TECH")].upper()
        or "NODDING" in frame.header[kw("DPR_TECH")].upper()
    ):
        object = frame.header[kw("OBJECT")].upper()
        ttype = (
            f"std_flux_stare_{object}".replace(" ", "_")
            .replace("-", "_")
            .replace("__", "_")
            .replace("__", "_")
        )
    return ttype


def _lamp_fragment(frame, kw):
    """build the lamp fragment of a filename from the DPR type

    **Key Arguments:**

    - ``frame`` -- the CCDData object frame
    - ``kw`` -- the keyword lookup function

    **Return:**

    - ``lamp`` -- ``_QLAMP``, ``_DLAMP`` or an empty string
    """
    if ",Q" in frame.header[kw("DPR_TYPE")].upper():
        lamp = "_QLAMP"
    elif ",D" in frame.header[kw("DPR_TYPE")].upper():
        lamp = "_DLAMP"
    else:
        lamp = ""
    return lamp


def _mask_slit(frame, kw, ttype):
    """work out the mask or slit fragment of a filename from the DPR technique

    **Key Arguments:**

    - ``frame`` -- the CCDData object frame
    - ``kw`` -- the keyword lookup function
    - ``ttype`` -- the frame type from `_frame_type`

    **Return:**

    - ``maskSlit`` -- ``onepin``, ``multipin``, ``slit``, or *None* if the technique matches none
    """
    maskSlit = None
    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,PINHOLE":
        maskSlit = "onepin"
    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,MULTI-PINHOLE":
        maskSlit = "multipin"

    if frame.header[kw("DPR_TECH")].upper() == "ECHELLE,SLIT" and ttype in (
        "mflat",
        "flat",
    ):
        maskSlit = "slit"
    if frame.header[kw("DPR_TECH")].upper() in (
        "ECHELLE,SLIT,STARE",
        "ECHELLE,SLIT,NODDING",
    ) and ("object" in ttype or "std_flux" in ttype):
        maskSlit = "slit"
    return maskSlit


def filenamer(log, frame, keywordLookup=False, detectorLookup=False, settings=False):
    """Given a FITS object, use the SOXS file-naming scheme to return a filename
    to be used to save the FITS object to disk

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- the CCDData object frame
    - ``keywordLookup`` -- the keyword lookup dictionary (needed if `settings` not provided). Default *False*
    - ``detectorLookup`` -- the detector parameters (needed if `settings` not provided). Default *False*
    - ``settings`` -- the soxspipe settings dictionary (needed if `keywordLookup` and `detectorLookup` not provided).
      Default *False*

    **Return:**

    - ``filename`` -- standardised name for the input frame

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
    log.debug("starting the ``filenamer`` function")

    # GENERATE A FILENAME FOR THE FRAME BASED ON THE FILENAMING
    # CONVENTION
    kw = keywordLookup or keyword_lookup(log=log, settings=settings).get

    if detectorLookup:
        dp = detectorLookup
    else:
        arm = frame.header[kw("SEQ_ARM")]
        # DETECTOR PARAMETERS LOOKUP OBJECT. THE RESULT IS UNUSED BUT THE CALL STAYS:
        # IT RAISES LookupError FOR AN UNKNOWN ARM
        dp = detector_lookup(log=log, settings=settings).get(arm)  # noqa: F841

    dateStamp = frame.header[kw("DATE_OBS")].replace("-", ".").replace(":", ".")
    # THE OBSERVATION ID IS NOT IN THE NAME BUT THE READ STAYS: A MISSING KEYWORD RAISES KeyError
    obid = frame.header[kw("OBS_ID")]  # noqa: F841
    arm = frame.header[kw("SEQ_ARM")].lower()
    binning = _binning_fragment(frame)
    romode = _readout_fragment(log, frame, kw)

    filename = f"{dateStamp}_{arm}{binning}{romode}"

    # DETERMINE THE TYPE
    if kw("DPR_TYPE") not in frame.header and kw("PRO_TYPE") in frame.header:
        return None

    ttype = _frame_type(frame, kw)

    lamp = _lamp_fragment(frame, kw)

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

    maskSlit = _mask_slit(frame, kw, ttype)

    # EXTRA PARAMETERS NEEDED FOR SPECTRUM
    if frame.header[kw("DPR_TECH")].upper() != "IMAGE" and maskSlit is None:
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

    log.debug("completed the ``filenamer`` function")
    return filename
