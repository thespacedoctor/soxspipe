#!/usr/bin/env python
"""
*Build the cosmic-ray mask QC product from the pixels L.A.Cosmic flagged in each raw frame*

L.A.Cosmic ORs its hits into the bad-pixel mask of the frame. This module keeps the hits on their own, so a
recipe can write them as one multi-extension FITS file per set of files (SOF).

This module imports nothing from soxspipe, so any module can import it without creating an import cycle.

Author
: David Young

Date Created
: October 22, 2026
"""

import numpy as np
from astropy.io import fits


def new_cosmic_ray_mask(maskBefore, maskAfter):
    """*return the pixels cosmic-ray cleaning newly flagged, with the earlier mask removed*

    **Key Arguments:**

    - ``maskBefore`` -- the frame mask before cleaning (the default bad-pixel mask)
    - ``maskAfter`` -- the frame mask after cleaning

    **Return:**

    - ``crMask`` -- uint8 image, 1 where cleaning newly flagged a pixel and 0 elsewhere

    **Usage:**

    ```python
    crMask = new_cosmic_ray_mask(maskBefore=frame.mask.copy(), maskAfter=cleanedFrame.mask)
    ```
    """
    return (np.asarray(maskAfter, dtype=bool) & ~np.asarray(maskBefore, dtype=bool)).astype(np.uint8)


def cosmic_ray_entry_order(entry):
    """*return the sort key that puts the cosmic-ray mask entries in observation order*

    Frames with a numeric ``mjdObs`` come first, by observation time. Frames without one come after them, by file name.

    **Key Arguments:**

    - ``entry`` -- a cosmic-ray mask entry, as in `cosmic_ray_mask_hdulist`

    **Return:**

    - ``key`` -- the tuple to sort on

    **Usage:**

    ```python
    entries = sorted(self.cosmicRayMasks, key=cosmic_ray_entry_order)
    ```
    """
    if entry["mjdObs"] is None:
        return (1, 0.0, entry["fileName"])
    return (0, entry["mjdObs"], entry["fileName"])


def cosmic_ray_mask_hdulist(entries):
    """*build the multi-extension FITS file of cosmic-ray masks, one image extension per frame*

    **Key Arguments:**

    - ``entries`` -- one dictionary per cleaned frame, in the order to write them, with the keys ``fileName``
      (the raw frame's file name), ``mjdObs`` and ``dateObs`` (the raw frame's observation time, or None) and
      ``mask`` (the uint8 cosmic-ray mask)

    **Return:**

    - ``hduList`` -- an empty primary HDU followed by one ``CRMASK`` image extension per entry. The header of each
      extension names the raw frame (``ORIGFILE``), its ``MJD-OBS`` and ``DATE-OBS`` where known, and the number of
      flagged pixels (``NCRPIX``)

    **Usage:**

    ```python
    hduList = cosmic_ray_mask_hdulist(self.cosmicRayMasks)
    hduList.writeto(filePath, overwrite=True)
    ```
    """
    hduList = fits.HDUList([fits.PrimaryHDU()])
    for entry in entries:
        hdu = fits.ImageHDU(data=entry["mask"], name="CRMASK")
        hdu.header["ORIGFILE"] = (entry["fileName"], "Raw frame the mask was made from")
        if entry["mjdObs"] is not None:
            hdu.header["MJD-OBS"] = (entry["mjdObs"], "Observation start of the raw frame")
        if entry["dateObs"] is not None:
            hdu.header["DATE-OBS"] = (entry["dateObs"], "Observation start of the raw frame")
        hdu.header["NCRPIX"] = (int(entry["mask"].sum()), "Pixels flagged as cosmic rays")
        hduList.append(hdu)
    return hduList
