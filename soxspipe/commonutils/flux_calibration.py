#!/usr/bin/env python
"""
*Flux calibrate an extracted science spectrum using an instrument response function*

Author
: David Young

Date Created
: July 28, 2023
"""

import os
from typing import Any

os.environ["TERM"] = "vt100"


def _calculate_flux_calibration(
    counts: Any,
    exposureTime: float,
    responseFactors: Any,
    extinctionFactors: Any,
) -> Any:
    """Return calibrated flux values from the pipeline calibration factors."""
    return counts / exposureTime * extinctionFactors * responseFactors * 10**-17


def _response_factors_from_table(wavelengths: Any, tableWave: Any, tableResponse: Any) -> Any:
    """Interpolate the tabulated response onto the wavelengths in log space; NaN outside the valid table range.

    Table rows with a non-finite wavelength or a response that is not finite and strictly positive are ignored, and
    the remaining rows are sorted by wavelength. A ValueError is raised if no row remains.
    """
    import numpy as np

    wavelengths = np.asarray(wavelengths, dtype=float)
    tableWave = np.asarray(tableWave, dtype=float)
    tableResponse = np.asarray(tableResponse, dtype=float)
    isUsable = np.isfinite(tableWave) & np.isfinite(tableResponse) & (tableResponse > 0)
    if not isUsable.any():
        raise ValueError(
            "The RESPONSE table of the response product has no usable rows: "
            "every row has a non-finite wavelength or a response that is not finite and positive."
        )
    order = np.argsort(tableWave[isUsable], kind="stable")
    logResponse = np.interp(
        wavelengths,
        tableWave[isUsable][order],
        np.log(tableResponse[isUsable][order]),
        left=np.nan,
        right=np.nan,
    )
    return np.exp(logResponse)


def load_response_factors(responseFunction: str, wavelengths: Any, log: Any) -> Any:
    """*read the response product and return the response factor at each wavelength*

    The ``RESPONSE`` curve extension is interpolated in log space; wavelengths outside it get NaN. The factor scales
    flat-fielded electrons per second per unbinned pixel. A product from an older release has no basis card (see
    ``response_basis``) and is refused with a ValueError that asks for the response to be rebuilt.

    **Key Arguments:**

    - ``responseFunction`` -- path to the ``_RESP.fits`` product
    - ``wavelengths`` -- the wavelengths (nm) to evaluate
    - ``log`` -- logger

    **Return:**

    - ``responseFactors`` -- the response factor at each wavelength
    """
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table

    from soxspipe.commonutils.response_basis import require_response_basis

    with fits.open(responseFunction) as hdul:
        require_response_basis(hdul[0].header, responseFunction)

    curve = Table.read(responseFunction, format="fits", hdu="RESPONSE")
    factors = _response_factors_from_table(wavelengths, np.asarray(curve["WAVE"]), np.asarray(curve["RESPONSE"]))
    outsideCount = int(np.count_nonzero(np.isnan(factors) & ~np.isnan(np.asarray(wavelengths, dtype=float))))
    if outsideCount > 0:
        log.info(f"{outsideCount} science samples lie outside the response curve range and were set to NaN.")
    return factors


# OR YOU CAN REMOVE THE CLASS BELOW AND ADD A WORKER FUNCTION ... SNIPPET TRIGGER BELOW
# xt-worker-def


class flux_calibration:
    """
    *The worker class for the flux_calibration module*

    **Key Arguments:**

    - ``log`` -- logger
    - ``responseFunction`` -- the instrument response function.
    - ``extractedSpectrum`` -- the extracted science spectrum
    - ``settings`` -- the settings dictionary

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/main/initialisation.html).

    To initiate a flux_calibration object, use the following:

    :::{todo}
        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``flux_calibration`` to documentation
        - create a blog post about what ``flux_calibration`` does
    :::

    ```python
    usage code
    ```

    """

    # Initialisation
    # 1. @flagged: what are the unique Attributes for each object? Add them
    # to __init__

    def __init__(
        self,
        log,
        responseFunction,
        extractedSpectrum,
        settings=False,
        airmass=1.0,
        exptime=1.0,
        extinctionPath="",
        arm="",
        header=None,
        recipeName="",
        startNightDate="",
        sofName="",
        debug=False,
    ):
        self.log = log
        log.debug("instantiating a new 'flux_calibration' object")
        self.settings = settings
        self.responseFunction = responseFunction
        self.extractedSpectrum = extractedSpectrum
        self.airmass = airmass
        self.exptime = exptime
        self.extinctionPath = extinctionPath
        self.arm = arm
        self.debug = debug
        self.header = header
        self.recipeName = recipeName
        self.startNightDate = startNightDate
        self.sofName = sofName

        import pandas as pd

        from soxspipe.commonutils import keyword_lookup
        from soxspipe.commonutils.toolkit import utility_setup

        self.kw = keyword_lookup(log=self.log, settings=self.settings).get

        self.qcDir, self.productDir = utility_setup(
            log=self.log,
            settings=settings,
            recipeName=recipeName,
            startNightDate=startNightDate,
        )
        self.products = pd.DataFrame()

        return

    def calibrate(self):
        """
        *flux calibrate the science spectrum*

        **Return:**

        - ``flux_calibration``

        **Usage:**

        :::{todo}
            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed
        :::

        ```python
        usage code
        ```
        """
        self.log.debug("starting the ``calibrate`` method")

        import copy
        from contextlib import suppress

        import pandas as pd

        from soxspipe.commonutils.phase3 import write_fits_table_to_disk
        from soxspipe.commonutils.response_basis import counts_per_unbinned_pixel, dispersion_binning
        from soxspipe.commonutils.toolkit import extinction_correction_factor

        flux_calibration = None

        self.log.debug("completed the ``calibrate`` method")
        # STEP TO DO:

        # APPLY EXTINCTION CORRECTION FACTOR
        if self.arm == "UVB" or self.arm == "VIS":
            if self.airmass is None:
                raise ValueError(
                    f"No airmass is known for this {self.arm} spectrum, so it cannot be corrected for extinction; "
                    "check the ESO TEL AIRM START and END cards of its raw frames"
                )
            extinctionCorrectionFactor = extinction_correction_factor(
                self.extractedSpectrum["WAVE"], self.extinctionPath, self.airmass
            )
        else:
            extinctionCorrectionFactor = 1.0

        # APPLY RESPNSE FUNCTION
        responseFactors = load_response_factors(self.responseFunction, self.extractedSpectrum["WAVE"], self.log)

        # THE RESPONSE SCALES FLAT-FIELDED ELECTRONS PER UNBINNED PIXEL, SO A BINNED SPECTRUM IS SCALED DOWN FIRST
        counts = counts_per_unbinned_pixel(self.extractedSpectrum, dispersion_binning(self.header, self.arm, self.kw))
        flux_calibration = _calculate_flux_calibration(
            counts,
            self.exptime,
            responseFactors,
            extinctionCorrectionFactor,
        )

        fluxCalSpectrum = pd.DataFrame(
            {
                "WAVE": self.extractedSpectrum["WAVE"],
                "FLUX_CALIBRATED": flux_calibration,
            }
        )

        if self.debug:
            from matplotlib import pyplot as plt

            plt.plot(self.extractedSpectrum["WAVE"], flux_calibration * 10**-17)
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Flux (erg/cm2/s/Angstrom)")
            plt.show()

        header = copy.deepcopy(self.header)
        with suppress(KeyError):
            header.pop(self.kw("DPR_CATG"))
        with suppress(KeyError):
            header.pop(self.kw("DPR_TYPE"))
        with suppress(KeyError):
            header.pop(self.kw("DET_READ_SPEED"))
        with suppress(KeyError):
            header.pop(self.kw("CONAD"))
        with suppress(KeyError):
            header.pop(self.kw("GAIN"))
        with suppress(KeyError):
            header.pop(self.kw("RON"))

        header[self.kw("SEQ_ARM")] = self.arm
        header["HIERARCH " + self.kw("PRO_TYPE")] = "REDUCED"
        header["HIERARCH " + self.kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{self.arm}".upper()
        if self.airmass is not None:
            header["HIERARCH " + self.kw("PRO_AIRM_MEAN")] = (
                float(self.airmass),
                "mean airmass of the frames in the product",
            )

        filename = f"{self.sofName}_FLUXCAL.fits"
        filePath = f"{self.productDir}/{filename}"

        write_fits_table_to_disk(
            log=self.log, settings=self.settings, header=header, tables=[fluxCalSpectrum], filePath=filePath, qc=None
        )

        from datetime import datetime

        utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

        self.products = pd.concat(
            [
                self.products,
                pd.DataFrame(
                    [
                        {
                            "soxspipe_recipe": self.recipeName,
                            "product_label": "EXTRACTED_FLUXCAL_SPECTRUM",
                            "file_name": filename,
                            "file_type": "FITS",
                            "reduction_date_utc": utcnow,
                            "product_desc": "Flux calibrated extracted spectrum",
                            "file_path": filePath,
                            "obs_date_utc": header["DATE-OBS"],
                            "label": "PROD",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        return filePath, self.products

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
