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

# THE COLUMNS OF THE EXTRACTED SPECTRUM THE FLUX CALIBRATION READS, AND THE ONE IT CARRIES WHEN THE SPECTRUM HAS IT
REQUIRED_SPECTRUM_COLUMNS = ("WAVE", "FLUX_COUNTS", "VARIANCE", "SNR")
OPTIONAL_SPECTRUM_COLUMNS = ("SKY_COUNTS",)


def _calculate_flux_calibration(
    counts: Any,
    exposureTime: float,
    responseFactors: Any,
    extinctionFactors: Any,
) -> Any:
    """Return calibrated flux values from the pipeline calibration factors."""
    return counts / exposureTime * extinctionFactors * responseFactors * 10**-17


def _require_spectrum_columns(spectrum: Any) -> None:
    """Raise a ValueError that names the columns of the extracted spectrum the flux calibration needs and lacks."""
    missing = [name for name in REQUIRED_SPECTRUM_COLUMNS if name not in spectrum.columns]
    if missing:
        raise ValueError(
            f"The extracted spectrum is missing the column(s) {', '.join(missing)}, so it cannot be flux calibrated; "
            f"it needs {', '.join(REQUIRED_SPECTRUM_COLUMNS)}."
        )


def _column_values(spectrum: Any, name: str, dtype: Any = None) -> Any:
    """Return a column as a plain array, with masked cells set to NaN.

    ``np.asarray`` on an astropy ``MaskedColumn`` returns the data under the mask, which holds fill values. A masked
    cell is missing data, so it becomes NaN in a floating-point column. Other columns are returned as they are.
    """
    import numpy as np

    column = spectrum[name]
    if np.ma.isMaskedArray(column):
        column = np.ma.asarray(column) if dtype is None else np.ma.asarray(column, dtype=dtype)
        if column.dtype.kind == "f":
            return column.filled(np.nan)
    return np.asarray(column, dtype=dtype)


def _flux_calibrated_columns(spectrum: Any, fluxCalibrated: Any, errorCalibrated: Any) -> dict[str, Any]:
    """Return the columns of the flux-calibrated table, in order.

    The calibrated flux and error come first. The counts, variance, sky and SNR of the extracted spectrum follow
    unchanged, so the sky stays in electrons per pixel. ``SKY_COUNTS`` is carried only when the spectrum has it: the
    nod and offset recipes extract no sky model.
    """
    import numpy as np

    columns = {
        "WAVE": _column_values(spectrum, "WAVE", float),
        "FLUX_CALIBRATED": np.asarray(fluxCalibrated, dtype=float),
        "ERR_CALIBRATED": np.asarray(errorCalibrated, dtype=float),
        "FLUX_COUNTS": _column_values(spectrum, "FLUX_COUNTS"),
        "VARIANCE": _column_values(spectrum, "VARIANCE"),
    }
    for name in OPTIONAL_SPECTRUM_COLUMNS:
        if name in spectrum.columns:
            columns[name] = _column_values(spectrum, name)
    columns["SNR"] = _column_values(spectrum, "SNR")
    return columns


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

        import numpy as np
        import pandas as pd

        from soxspipe.commonutils.fits_headers import set_header_card
        from soxspipe.commonutils.phase3 import PER_ARM_MODE, write_fits_table_to_disk
        from soxspipe.commonutils.response_basis import counts_per_unbinned_pixel, dispersion_binning
        from soxspipe.commonutils.toolkit import extinction_correction_factor

        flux_calibration = None
        _require_spectrum_columns(self.extractedSpectrum)

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
        binning = dispersion_binning(self.header, self.arm, self.kw)
        counts = counts_per_unbinned_pixel(self.extractedSpectrum, binning)
        flux_calibration = _calculate_flux_calibration(
            counts,
            self.exptime,
            responseFactors,
            extinctionCorrectionFactor,
        )
        # THE ERROR IS THE SQUARE ROOT OF THE VARIANCE, PASSED THROUGH THE SAME FACTOR AS THE FLUX
        errorCalibrated = _calculate_flux_calibration(
            np.sqrt(np.asarray(self.extractedSpectrum["VARIANCE"], dtype=float)) / binning,
            self.exptime,
            responseFactors,
            extinctionCorrectionFactor,
        )

        fluxCalSpectrum = pd.DataFrame(
            _flux_calibrated_columns(self.extractedSpectrum, flux_calibration, errorCalibrated)
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

        set_header_card(header, self.kw("SEQ_ARM"), self.arm)
        set_header_card(header, self.kw("PRO_TYPE"), "REDUCED")
        set_header_card(header, self.kw("PRO_CATG"), f"SCI_SLIT_FLUX_{self.arm}".upper())
        if self.airmass is not None:
            set_header_card(
                header,
                self.kw("PRO_AIRM_MEAN"),
                float(self.airmass),
                "mean airmass of the frames in the product",
            )

        filename = f"{self.sofName}_FLUXCAL.fits"
        filePath = f"{self.productDir}/{filename}"

        write_fits_table_to_disk(
            log=self.log,
            settings=self.settings,
            header=header,
            tables=[fluxCalSpectrum],
            filePath=filePath,
            qc=None,
            phase3=PER_ARM_MODE,
            extname="SPECTRUM",
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
