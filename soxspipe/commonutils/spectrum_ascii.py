#!/usr/bin/env python
"""
*Write the ASCII version of an extracted spectrum from its FITS binary table*

The ASCII file comes from the flux-calibrated table when the spectrum was flux calibrated. Otherwise it comes from the
order-merged (or stacked) table, so a run with no response product still gets an ASCII spectrum, in counts.

Author
: David Young

Date Created
: November 19, 2026
"""

import os

FLUXCAL_DESCRIPTION = "ASCII version of the flux-calibrated extracted spectrum"
COUNTS_DESCRIPTION = "ASCII version of the extracted spectrum (not flux calibrated)"
ASCII_PRODUCT_LABEL = "EXTRACTED_MERGED_ASCII"
NM_TO_ANGSTROM = 10


def ascii_spectrum_path(mergedFitsPath: str, fluxcalFitsPath: str | None = None) -> str:
    """*return the path of the ASCII spectrum*

    **Key Arguments:**

    - ``mergedFitsPath`` -- the path to the order-merged (or stacked) spectrum table
    - ``fluxcalFitsPath`` -- the path to the flux-calibrated table, or None when the spectrum was not flux calibrated

    **Return:**

    - ``asciiPath`` -- the source table's path with a ``.txt`` extension, so ``<stem>_FLUXCAL.txt`` for a
      flux-calibrated spectrum and ``<stem>_EXTRACTED_MERGED.txt`` for the fallback
    """
    return os.path.splitext(fluxcalFitsPath or mergedFitsPath)[0] + ".txt"


def write_ascii_spectrum(sourceFitsPath: str, asciiPath: str) -> None:
    """*write a spectrum table to an ASCII file with the wavelength in Angstrom*

    The FITS table is read, never changed.

    **Key Arguments:**

    - ``sourceFitsPath`` -- the path to the FITS file whose first table extension is written
    - ``asciiPath`` -- the path to the ASCII file. An existing file is replaced
    """
    from astropy.table import Table

    # NAN VALUES (NO SKY MODEL, OUTSIDE THE RESPONSE RANGE) MUST STAY NAN, NOT BECOME MASKED, SO THEY ARE WRITTEN AS NAN
    table = Table.read(sourceFitsPath, format="fits", hdu=1, mask_invalid=False)
    table["WAVE"] = table["WAVE"] * NM_TO_ANGSTROM
    table["WAVE"].format = "{:.2f}"
    table.write(asciiPath, format="ascii", overwrite=True)


def write_ascii_spectrum_product(
    mergedFitsPath: str,
    fluxcalFitsPath: str | None,
    productsTable,
    recipeName: str,
    obsDateUtc: str,
):
    """*write the ASCII spectrum of a pass once and record it as a product*

    **Key Arguments:**

    - ``mergedFitsPath`` -- the path to the order-merged (or stacked) spectrum table
    - ``fluxcalFitsPath`` -- the path to the flux-calibrated table, or None when the pass has none. The ASCII file is
      written from it when it is given
    - ``productsTable`` -- the products table to add the ASCII file to. It is not changed
    - ``recipeName`` -- the recipe name to record against the product
    - ``obsDateUtc`` -- the observation date (UTC) to record against the product

    **Return:**

    - ``productsTable`` -- a new products table with the ``EXTRACTED_MERGED_ASCII`` row added
    """
    from soxspipe.commonutils.toolkit import append_product, utcnow_string

    asciiPath = ascii_spectrum_path(mergedFitsPath, fluxcalFitsPath)
    write_ascii_spectrum(fluxcalFitsPath or mergedFitsPath, asciiPath)
    return append_product(
        productsTable,
        recipeName=recipeName,
        productLabel=ASCII_PRODUCT_LABEL,
        fileName=os.path.basename(asciiPath),
        filePath=asciiPath,
        productDesc=FLUXCAL_DESCRIPTION if fluxcalFitsPath else COUNTS_DESCRIPTION,
        obsDateUtc=obsDateUtc,
        reductionDateUtc=utcnow_string(),
        fileType="TXT",
        label="PROD",
    )
