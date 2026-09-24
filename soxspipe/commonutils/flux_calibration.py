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
    wavelengths: Any,
    counts: Any,
    exposureTime: float,
    responseCoefficients: Any,
    extinctionFactors: Any,
) -> Any:
    """Return calibrated flux values from the pipeline calibration factors."""
    import numpy as np

    responseFactors = np.polyval(responseCoefficients, wavelengths)
    return counts / exposureTime * extinctionFactors * responseFactors * 10**-17


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
        from astropy.table import Table

        from soxspipe.commonutils.phase3 import write_fits_table_to_disk
        from soxspipe.commonutils.toolkit import extinction_correction_factor

        flux_calibration = None

        self.log.debug("completed the ``calibrate`` method")
        # STEP TO DO:

        # APPLY EXTINCTION CORRECTION FACTOR
        if self.arm == "UVB" or self.arm == "VIS":
            extinctionCorrectionFactor = extinction_correction_factor(
                self.extractedSpectrum["WAVE"], self.extinctionPath, self.airmass
            )
        else:
            extinctionCorrectionFactor = 1.0

        # APPLY RESPNSE FUNCTION
        responseFunctionCoeff = Table.read(self.responseFunction, format="fits")
        # NOW UNPACK THE COEFFICIENTS
        polyCoeffs = []
        for idx_coeff in range(0, int(responseFunctionCoeff["polyOrder"]) + 1):
            polyCoeffs.append(responseFunctionCoeff[f"c{idx_coeff}"][0])

        flux_calibration = _calculate_flux_calibration(
            self.extractedSpectrum["WAVE"],
            self.extractedSpectrum["FLUX_COUNTS"],
            self.exptime,
            polyCoeffs,
            extinctionCorrectionFactor,
        )

        fluxCalSpectrum = pd.DataFrame(
            {
                "WAVE": self.extractedSpectrum["WAVE"],
                "FLUX_CALIBRATED": flux_calibration,
            }
        )

        if self.debug:
            import matplotlib

            matplotlib.use("TkAgg")
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
                pd.DataFrame([
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
                ]),
            ],
            ignore_index=True,
        )

        return filePath, self.products

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
