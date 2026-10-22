#!/usr/bin/env python
"""
*Given a standard star extracted spectrum, generate the instrument response function needed to flux calibrate science spectra*

:Author:
    Marco Landoni & David Young

:Date Created:
    July 28, 2023
"""

import os
from typing import Any

from soxspipe.commonutils.response_basis import (
    RESPONSE_BASIS,
    RESPONSE_BASIS_COMMENT,
    RESPONSE_BASIS_KEYWORD,
    counts_per_unbinned_pixel,
    dispersion_binning,
)
from soxspipe.commonutils.response_fit import (
    _ResponseFitConvergenceError,
    fit_response_curve,
    resolve_response_settings,
)
from soxspipe.commonutils.toolkit import extinction_correction_factor

os.environ["TERM"] = "vt100"

__all__ = ["_ResponseFitConvergenceError", "response_function"]


QC_COLUMNS = (
    "soxspipe_recipe",
    "qc_name",
    "qc_value",
    "qc_unit",
    "qc_order",
    "qc_comment",
    "obs_date_utc",
    "reduction_date_utc",
    "to_header",
)
PRODUCT_COLUMNS = (
    "soxspipe_recipe",
    "product_label",
    "file_name",
    "file_type",
    "obs_date_utc",
    "reduction_date_utc",
    "product_desc",
    "file_path",
    "label",
)


def _table_or_empty(table: Any, columns: tuple[str, ...]) -> Any:
    """Return ``table`` if it is a data frame, otherwise an empty data frame with the production ``columns``."""
    import pandas as pd

    if isinstance(table, pd.DataFrame):
        return table
    return pd.DataFrame({column: pd.Series([], dtype="object") for column in columns})


def _sampled_response_table(wavelength: Any, response: Any) -> Any:
    """Return the response at each finite wavelength, sorted by wavelength, as a WAVE/RESPONSE table.

    RESPONSE carries the scale flux calibration applies: flux = counts / exptime * RESPONSE * 1e-17. The counts are
    flat-fielded electrons per unbinned pixel (see ``response_basis``); the unit uses ``ct`` because the FITS unit
    standard has no electron.
    """
    import astropy.units as u
    import numpy as np
    from astropy.table import Table

    wave = np.asarray(wavelength, dtype=float)
    values = np.asarray(response, dtype=float)
    isFinite = np.isfinite(wave)
    order = np.argsort(wave[isFinite], kind="stable")
    responseUnit = u.Unit("1e-17 erg / (cm2 s Angstrom)") / (u.ct / u.s)
    return Table({"WAVE": wave[isFinite][order] * u.nm, "RESPONSE": values[isFinite][order] * responseUnit})


class response_function:
    """
    *Given a standard star extracted spectrum, generate the instrument response function needed to flux calibrate science spectra*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``stdExtractionPath`` -- fits binary table containing the extracted standard spectrum
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
        - ``settings`` --  the pipeline settings
        - ``sofName`` -- name of the originating SOF file
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""
        - ``stdNotFlatExtractionPath`` -- fits binary table containing the extracted standard spectrum without flat correction. Default "".
        - ``orderJoins`` -- a list of tuples indicating the orders to be joined together in the response function fitting. Default None.

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a response_function object, use the following:

    ```python
    from soxspipe.commonutils import response_function
    response = response_function(
        log=log,
        settings=settings,
        recipeName=recipeName,
        sofName=sofName,
        stdExtractionPath=stdExtractionPath
        qcTable=qcTable,
        productsTable=productsTable,
        startNightDate=startNightDate
        stdNotFlatExtractionPath=stdNotFlatExtractionPath,
        orderJoins=orderJoins,
    )
    qcTable, productsTable = response.get()
    ```
    """

    def __init__(
        self,
        log,
        stdExtractionPath,
        recipeName,
        sofName,
        settings=False,
        qcTable=False,
        productsTable=False,
        startNightDate="",
        stdNotFlatExtractionPath="",
        orderJoins=None,
    ):
        self.log = log
        log.debug("instansiating a new 'response_function' object")
        self.settings = settings
        self.stdExtractionPath = stdExtractionPath
        self.recipeSettings = settings["soxs-response"]
        self.instrument = self.settings["instrument"].lower()
        self.qc = _table_or_empty(qcTable, QC_COLUMNS)
        self.products = _table_or_empty(productsTable, PRODUCT_COLUMNS)
        self.recipeName = recipeName
        self.sofName = sofName
        self.orderJoins = orderJoins

        from astropy.io import fits
        from astropy.table import Table

        from soxspipe.commonutils import detector_lookup, keyword_lookup
        from soxspipe.commonutils.toolkit import get_calibrations_path

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(log=self.log, settings=self.settings).get
        kw = self.kw

        hdul = fits.open(self.stdExtractionPath)
        self.header = hdul[0].header
        self.arm = self.header[kw("SEQ_ARM")].strip().upper()  # KW lookup
        self.detectorParams = detector_lookup(log=log, settings=settings).get(self.arm)
        self.dateObs = self.header[kw("DATE_OBS")]
        self.texp = float(self.header[kw("EXPTIME")])
        self.dispersionBinning = dispersion_binning(self.header, self.arm, kw)

        # CONVERTING EXTRACTION (BACK) TO DATAFRAME
        self.stdExtractionDF = Table.read(self.stdExtractionPath, format="fits")
        self.stdExtractionDF = self.stdExtractionDF.to_pandas()

        # SORT BY COLUMN NAME
        self.stdExtractionDF.sort_values(["WAVE"], inplace=True, kind="stable")

        self.calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)

        # GET THE ABSOLUTE STANDARD STAR FLUXES, ASSUMING TO HAVE 1-1 MAPPING BETWEEN OBJECT NAME IN THE FITS HEADER AND DATABASE
        stdAbsFluxDF = Table.read(
            self.calibrationRootPath + "/" + self.detectorParams["flux-standards"],
            format="fits",
        )
        self.stdAbsFluxDF = stdAbsFluxDF.to_pandas()
        # MAKE ALL COLUMNS UPPERCASE
        self.stdAbsFluxDF.columns = [d.upper() for d in self.stdAbsFluxDF.columns]

        # NAME: AKA
        stdAkas = {
            "LTT7987": "CD3017706",
            "EG274": "CD3810980",
            "LTT3218": "CD325613",
            "EG21": "CPD69177",
        }

        self.std_objName = ""
        for ii in kw("OBS_TARG_NAME"), kw("OBJECT"), kw("OBS_NAME"):
            if ii in self.header:
                self.std_objName += self.header[ii].strip().upper().replace(" ", "").replace("-", "").replace("_", "")

        for k, v in stdAkas.items():
            if v in self.std_objName:
                self.std_objName = k
                break

        for ii in self.stdAbsFluxDF.columns:
            if ii in self.std_objName:
                self.std_objName = ii
                break

        if stdNotFlatExtractionPath and len(stdNotFlatExtractionPath) > 1:
            # STD STAR GIVEN, READING THE NON FLAT FIELDED SPECTRUM
            self.stdExtractionNotFlatDF = Table.read(stdNotFlatExtractionPath, format="fits")
            self.stdExtractionNotFlatDF = self.stdExtractionNotFlatDF.to_pandas()

        self.log.print(f"STANDARD-STAR: {self.std_objName}")
        # USING THE AVERAGE AIR MASS
        if self.instrument == "soxs":
            # NEED TO UPDATE THE LOOKUP TABLE EVENTUALLY
            airmass_start = float(self.header["HIERARCH ESO TEL AIRM START"])
            airmass_end = float(self.header["HIERARCH ESO TEL AIRM END"])
        else:
            airmass_start = float(self.header[kw("AIRM_START")])
            airmass_end = float(self.header[kw("AIRM_END")])
        self.airmass = (airmass_start + airmass_end) / 2

        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log,
            settings=settings,
            recipeName=recipeName,
            startNightDate=startNightDate,
        )

        return

    def get(self):
        """
        *get the response_function object*

        **Return:**
            - ``response_function`` -- a set of polynomial coefficients
        """
        self.log.debug("starting the ``get`` method")

        import numpy as np
        import pandas as pd
        from matplotlib import pyplot as plt
        from scipy.interpolate import interp1d
        from scipy.signal import savgol_filter

        response_function = None

        # GET THE EXTRACTED STANDARD STAR'S WAVELENGTH AND FLUX
        stdExtWave = self.stdExtractionDF["WAVE"].values
        # FLAT-FIELDED ELECTRONS PER UNBINNED PIXEL: THE FLAT HAS ALREADY DIVIDED OUT THE PIXEL SCALE
        stdExtFlux = counts_per_unbinned_pixel(self.stdExtractionDF, self.dispersionBinning)
        stdExtWaveNotFlat = self.stdExtractionNotFlatDF["WAVE"].values
        stdExtFluxNotFlat = self.stdExtractionNotFlatDF["FLUX_DENSITY_COUNTS"].values

        if self.std_objName not in self.stdAbsFluxDF.columns:
            message = (
                f"Standard star {self.std_objName} not found in the static calibration database. "
                f"The available STDs are {', '.join(self.stdAbsFluxDF.columns[1:])}"
            )
            self.log.error(message)
            raise LookupError(message)

        # SELECTING ROWS IN THE INTERESTED WAVELENGTH RANGE ADDING A MARGIN TO THE RANGE
        stdAbsFluxDF = self.stdAbsFluxDF
        stdAbsFluxDF = stdAbsFluxDF[
            (stdAbsFluxDF["WAVE"] > np.min(stdExtWaveNotFlat) - 10)
            & (stdAbsFluxDF["WAVE"] < 10 + np.max(stdExtWaveNotFlat))
        ]

        # FLUX IS CONVERTED IN ERG / CM2 / S / ANGs
        self.std_wavelength_to_abs_flux = interp1d(
            np.array(stdAbsFluxDF["WAVE"]),
            np.array(stdAbsFluxDF[self.std_objName]) * 10**17,
            kind="next",
            fill_value="extrapolate",
        )

        # THE NOT-FLATTENED FLUX DENSITY IS ELECTRONS PER NM; DIVIDING BY 10 CONVERTS IT TO ELECTRONS PER ANGSTROM
        stdExtFluxNotFlat = stdExtFluxNotFlat / 10.0

        # NOW DIVIDING FOR THE EXPOSURE TIME
        stdExtFlux = stdExtFlux / self.texp
        stdExtFluxNotFlat = stdExtFluxNotFlat / self.texp

        if False:
            plt.plot(stdExtWave, stdExtFlux)
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Extracted Flux (counts/s)")
            plt.title("Extracted Standard Star Spectrum")
            plt.grid()
            plt.show()

            plt.plot(stdAbsFluxDF["WAVE"], stdAbsFluxDF[self.std_objName] * 10**17)
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Absolute Flux (erg/s/cm²/Å)")
            plt.title("Absolute Standard Star Spectrum")
            plt.grid()
            plt.show()
        # GETTING EFFICIENCY
        stdEfficiencyEstimate = None
        # USING THE STD STAR SPECTRUM THAT IS NOT CORRECTED BY FLAT
        try:
            # CONSTANTS FOR CALCULATIONS
            c = 3 * 10**10  # SPEED OF LIGHT IN CM/S
            h = 6.63 * 10**-27  # PLANCK'S CONSTANT IN ERG·S

            # TELESCOPE AREA BASED ON THE INSTRUMENT
            if self.instrument == "xsh":
                # VLT TELESCOPE AREA (IN CM²)
                area = 400 * 400 * 3.14
            else:
                # NTT TELESCOPE AREA (IN CM²)
                area = 89000

            # EXTRACT WAVELENGTH VALUES AND COMPUTE ABSOLUTE FLUX
            wavelength = stdExtWaveNotFlat
            abs_flux = self.std_wavelength_to_abs_flux(wavelength) * 10**-17
            # CONVERT WAVELENGTH FROM NM TO CM
            wavelength_in_cm = wavelength * 10**-7

            # COMPUTE THE ABSOLUTE PHOTON FLUX OF THE STANDARD STAR
            hc = h * c  # PRECOMPUTE PLANCK'S CONSTANT * SPEED OF LIGHT
            stdAbsPhotonFlux = (abs_flux * wavelength_in_cm / hc) * area

            # COMPUTE THE EFFICIENCY ESTIMATE
            stdEfficiencyEstimate = stdExtFluxNotFlat / stdAbsPhotonFlux
            # Change backend to MacOSX
            if False:
                import matplotlib

                matplotlib.use("MacOSX")
                from matplotlib import pyplot as plt

                print(self.texp)
                plt.plot(
                    wavelength,
                    stdExtFluxNotFlat,
                    label="Extracted Flux (not flat corrected)",
                )
                plt.plot(wavelength, stdAbsPhotonFlux, label="Absolute Photon Flux")
                plt.plot(wavelength, stdEfficiencyEstimate, label="Efficiency Estimate")
                plt.xlabel("Wavelength (nm)")
                plt.ylabel("Efficiency")
                plt.title("Efficiency Estimate")
                plt.grid()
                plt.legend()
                plt.show()

            # APPLY SAVITZKY-GOLAY FILTER TO SMOOTH THE EFFICIENCY ESTIMATE
            stdEfficiencyEstimate = savgol_filter(stdEfficiencyEstimate, 21, 2)
        except Exception as e:
            print(e)
            # LANDING HERE NO STD STAR SPECTRUM WITHOUT FLAT CORRECTION IS PROVIDED
            stdEfficiencyEstimate = None
            pass

        # APPLYING EXTINCTION CORRECTION
        if self.arm == "UVB" or self.arm == "VIS":
            extCorrectionFactor = extinction_correction_factor(
                stdExtWave,
                self.calibrationRootPath + "/" + self.detectorParams["extinction"],
                self.airmass,
            )
            if False:
                from matplotlib import pyplot as plt

                plt.plot(stdExtWave, extCorrectionFactor)
                plt.xlabel("Wavelength (nm)")
                plt.ylabel("Extinction Correction Factor")
                plt.title("Extinction Correction Factor vs Wavelength")
                plt.grid()
                plt.show()
            stdExtFlux = stdExtFlux * extCorrectionFactor

        self.response_function_raw = self.std_wavelength_to_abs_flux(stdExtWave) / stdExtFlux
        fit, fitSettings = self._fit_response(stdExtWave, self.response_function_raw)
        self.response_function_raw = fit.ratio
        qcWarnings = self._add_response_qcs(fit, fitSettings)

        # WRITE RESPONSE FUNCTION TO FITS BINARY TABLE
        self.write_response_function_to_file(wavelength=stdExtWave, response=fit.evaluate(stdExtWave))

        if stdEfficiencyEstimate is not None:
            # CREATE A DATAFRAME FOR EFFICIENCY ESTIMATE
            stdEfficiencyEstimateDF = pd.DataFrame({"WAVE": stdExtWaveNotFlat, "EFFICIENCY": stdEfficiencyEstimate})
            # WRITE THE EFFICIENCY ESTIMATE TO FITS BINARY TABLE
            import copy

            from soxspipe.commonutils.phase3 import write_fits_table_to_disk
            from soxspipe.commonutils.toolkit import add_snr_efficiency_qcs

            filename = f"{self.sofName}_EFFICIENCY.fits"
            filepath = f"{self.productDir}/{filename}"

            header = copy.deepcopy(self.header)
            header[self.kw("SEQ_ARM").upper()] = self.arm
            header[self.kw("PRO_TYPE").upper()] = "REDUCED"
            header[self.kw("PRO_CATG").upper()] = f"EFFICIENCY_TAB_{self.arm}".upper()

            self.qc = add_snr_efficiency_qcs(
                log=self.log,
                spectrumDF=stdEfficiencyEstimateDF,
                qcTable=self.qc,
                orderJoins=self.orderJoins,
                recipeName=self.recipeName,
                dateObs=self.dateObs,
                arm=self.arm,
            )

            write_fits_table_to_disk(
                log=self.log,
                settings=self.settings,
                header=header,
                tables=[stdEfficiencyEstimateDF],
                filePath=filepath,
                qc=self.qc,
            )

            self._register_product("EFFICIENCY", filename, "FITS", "SOXS efficiency estimate", filepath, "QC")

        self.plot_response_curve(
            stdExtWave=stdExtWave,
            stdExtWaveNotFlat=stdExtWaveNotFlat,
            stdExtFlux=stdExtFlux,
            fit=fit,
            fitSettings=fitSettings,
            stdEfficiencyEstimate=stdEfficiencyEstimate,
            qcWarnings=qcWarnings,
        )

        return self.qc, self.products, False

    def _register_product(self, productLabel, fileName, fileType, productDesc, filePath, label):
        """*append one row describing an output product to the products table*

        **Key Arguments:**

        - ``productLabel`` -- the product label, for example ``RESPONSE_FUNC``
        - ``fileName`` -- the product file name
        - ``fileType`` -- the file type, for example ``FITS``
        - ``productDesc`` -- a short description of the product
        - ``filePath`` -- the path to the product file
        - ``label`` -- ``PROD`` for a data product or ``QC`` for a quality-control product
        """
        from datetime import datetime

        import pandas as pd

        self.products = pd.concat(
            [
                self.products,
                pd.DataFrame(
                    [
                        {
                            "soxspipe_recipe": self.recipeName,
                            "product_label": productLabel,
                            "file_name": fileName,
                            "file_type": fileType,
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
                            "product_desc": productDesc,
                            "file_path": filePath,
                            "label": label,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

    def _fit_response(self, stdExtWave, rawRatio):
        """*fit the spline response curve to the raw ratio of catalogue to extracted flux*

        **Key Arguments:**

        - ``stdExtWave`` -- the extracted standard star wavelengths (nm)
        - ``rawRatio`` -- the raw response at each wavelength

        **Return:**

        - ``fit`` -- the ``ResponseFit``
        - ``fitSettings`` -- the complete settings for this arm
        """
        fitSettings = resolve_response_settings(self.recipeSettings, self.arm, self.instrument)
        try:
            fit = fit_response_curve(stdExtWave, rawRatio, fitSettings, log=self.log)
        except _ResponseFitConvergenceError as e:
            self.log.print("fail")
            raise Exception("The fitting of response function did not converge!") from e.originalError
        return fit, fitSettings

    def _add_response_qcs(self, fit, fitSettings):
        """*measure the response fit QCs, add them to the QC table and warn about any out of tolerance*

        **Return:**

        - ``qcWarnings`` -- one message per out-of-tolerance measure. The recipe never fails because of them.
        """
        from soxspipe.commonutils.response_qc import measure_response_qcs, response_qc_rows, response_qc_warnings

        measures = measure_response_qcs(fit, fitSettings["selfcal_bands"])
        self.qc = response_qc_rows(self.qc, measures, self.recipeName, self.dateObs)
        qcWarnings = response_qc_warnings(measures, fitSettings)
        for message in qcWarnings:
            self.log.warning(f"Response QC out of tolerance: {message}")
        return qcWarnings

    def plot_response_curve(
        self,
        stdExtWave,
        stdExtWaveNotFlat,
        stdExtFlux,
        fit,
        fitSettings,
        stdEfficiencyEstimate,
        qcWarnings,
    ):
        """*generate a QC plot for the response curve*

        **Key Arguments:**

        - ``stdExtWave`` -- the extracted standard star wavelength
        - ``stdExtWaveNotFlat`` -- the extracted standard star wavelength (not flattened)
        - ``stdExtFlux`` -- the flat-fielded standard star flux, in e-/s per unbinned pixel and extinction corrected
        - ``fit`` -- the ``ResponseFit``
        - ``fitSettings`` -- the settings for this arm, which hold the masked regions
        - ``stdEfficiencyEstimate`` -- the estimated instrument efficiency
        - ``qcWarnings`` -- out-of-tolerance QC messages, shown in red in the plot title

        **Return:**
            - `plotFilePath` -- the path to the QC plot PDF
        """
        self.log.debug("starting the ``plot_response_curve`` method")

        import matplotlib.pyplot as plt

        from soxspipe.commonutils.response_plot import draw_response_qc_figure

        fig = draw_response_qc_figure(
            fit=fit,
            fitSettings=fitSettings,
            stdExtWave=stdExtWave,
            stdExtFlux=stdExtFlux,
            absFlux=self.std_wavelength_to_abs_flux(stdExtWave),
            stdExtWaveNotFlat=stdExtWaveNotFlat,
            stdEfficiencyEstimate=stdEfficiencyEstimate,
            objName=self.std_objName,
            qcWarnings=qcWarnings,
        )

        plotFilename = self.sofName + "_RESPONSE.pdf"
        plotFilePath = f"{self.qcDir}/{plotFilename}"
        try:
            fig.savefig(plotFilePath, dpi=120, format="pdf", bbox_inches="tight")
        finally:
            plt.close(fig)

        self._register_product("RESPONSE_QC_PLOT", plotFilename, "PDF", "Response curve QC plot.", plotFilePath, "QC")

        if False:
            plt.show()

        self.log.debug("completed the ``plot_response_curve`` method")
        return plotFilePath

    def write_response_function_to_file(self, wavelength: Any, response: Any) -> str:
        """*write out the response curve*

        The one table extension, ``RESPONSE``, holds the fitted spline response at each finite ``wavelength``, sorted
        by wavelength. Flux calibration reads it.

        **Key Arguments:**

        - ``wavelength`` -- the wavelengths (nm) of the extracted standard-star spectrum
        - ``response`` -- the fitted response at each of those wavelengths

        **Return:**

        - ``responseCurvePath`` -- path to the saved file
        """
        self.log.debug("starting the ``write_response_function_to_file`` method")

        import copy

        from astropy.io import fits

        from soxspipe.commonutils.phase3 import add_qc_to_header

        arm = self.arm
        kw = self.kw

        filename = self.sofName + "_RESP.fits"

        header = copy.deepcopy(self.header)

        filePath = f"{self.productDir}/{filename}"

        header[kw("SEQ_ARM").upper()] = arm
        header[kw("PRO_TYPE").upper()] = "REDUCED"
        header[kw("PRO_CATG").upper()] = f"RESP_TAB_{arm}".upper()
        header[RESPONSE_BASIS_KEYWORD] = (RESPONSE_BASIS, RESPONSE_BASIS_COMMENT)

        # WRITE QCs TO HEADERS
        add_qc_to_header(header, self.qc)

        priHDU = fits.PrimaryHDU(header=header)

        curveHDU = fits.table_to_hdu(_sampled_response_table(wavelength, response))
        curveHDU.name = "RESPONSE"

        hduList = fits.HDUList([priHDU, curveHDU])
        hduList.verify("fix")
        hduList.writeto(filePath, checksum=True, overwrite=True)

        self._register_product("RESPONSE_FUNC", filename, "FITS", "Response function curve.", filePath, "PROD")

        self.log.debug("completed the ``write_response_function_to_file`` method")
        return filePath
