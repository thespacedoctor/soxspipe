#!/usr/bin/env python
"""
*Write the ASCII version of an extracted spectrum from its FITS binary table*

The ASCII file comes from the ``<SOF stem>_SPECTRUM.fits`` table of the pass. The table holds the flux-calibrated
columns when the spectrum was flux calibrated. Otherwise it holds counts, so a run with no response product still gets
an ASCII spectrum, in counts.

Author
: David Young

Date Created
: November 19, 2026
"""

import os

import pandas as pd

from soxspipe.commonutils.spectrum_names import SPECTRUM_ASCII_LABEL

FLUXCAL_DESCRIPTION = "ASCII version of the flux-calibrated extracted spectrum"
COUNTS_DESCRIPTION = "ASCII version of the extracted spectrum (not flux calibrated)"
NM_TO_ANGSTROM = 10


def ascii_spectrum_path(spectrumFitsPath: str) -> str:
    """*return the path of the ASCII spectrum*

    **Key Arguments:**

    - ``spectrumFitsPath`` -- the path to the spectrum table, ``<SOF stem>_SPECTRUM.fits``

    **Return:**

    - ``asciiPath`` -- the table's path with a ``.txt`` extension, so ``<SOF stem>_SPECTRUM.txt``
    """
    return os.path.splitext(spectrumFitsPath)[0] + ".txt"


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
    # A MASKED CELL IS MISSING DATA; ASTROPY WOULD WRITE IT AS A BLANK, SO IT BECOMES NAN IN A FLOATING-POINT COLUMN
    for name in table.colnames:
        if hasattr(table[name], "mask") and table[name].dtype.kind == "f":
            table[name] = table[name].filled(float("nan"))
    table["WAVE"] = table["WAVE"] * NM_TO_ANGSTROM
    table["WAVE"].format = "{:.2f}"
    table.write(asciiPath, format="ascii", overwrite=True)


def write_ascii_spectrum_product(
    spectrumFitsPath: str,
    isFluxCalibrated: bool,
    productsTable: pd.DataFrame,
    recipeName: str,
    obsDateUtc: str,
) -> pd.DataFrame:
    """*write the ASCII spectrum of a pass once and record it as a product*

    **Key Arguments:**

    - ``spectrumFitsPath`` -- the path to the spectrum table, ``<SOF stem>_SPECTRUM.fits``
    - ``isFluxCalibrated`` -- does the table hold the flux-calibrated columns? It sets the product description
    - ``productsTable`` -- the products table to add the ASCII file to. It is not changed
    - ``recipeName`` -- the recipe name to record against the product
    - ``obsDateUtc`` -- the observation date (UTC) to record against the product

    **Return:**

    - ``productsTable`` -- a new products table with the ``SPECTRUM_ASCII`` row added
    """
    from soxspipe.commonutils.toolkit import append_product, utcnow_string

    asciiPath = ascii_spectrum_path(spectrumFitsPath)
    write_ascii_spectrum(spectrumFitsPath, asciiPath)
    return append_product(
        productsTable,
        recipeName=recipeName,
        productLabel=SPECTRUM_ASCII_LABEL,
        fileName=os.path.basename(asciiPath),
        filePath=asciiPath,
        productDesc=FLUXCAL_DESCRIPTION if isFluxCalibrated else COUNTS_DESCRIPTION,
        obsDateUtc=obsDateUtc,
        reductionDateUtc=utcnow_string(),
        fileType="TXT",
        label="PROD",
    )
