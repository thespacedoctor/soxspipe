#!/usr/bin/env python
# encoding: utf-8
"""
*Given a standard star extracted spectrum, generate the instrument response function needed to flux calibrate science spectra*

:Author:
    Marco Landoni & David Young

:Date Created:
    July 28, 2023
"""
import sys
import os
from builtins import object

from soxspipe.commonutils.toolkit import extinction_correction_factor
os.environ['TERM'] = 'vt100'


class response_function(object):
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
            stdNotFlatExtractionPath=""
    ):
        self.log = log
        log.debug("instansiating a new 'response_function' object")
        self.settings = settings
        self.stdExtractionPath = stdExtractionPath
        self.recipeSettings = settings['soxs-response']
        self.instrument = self.settings["instrument"].lower()
        self.qc = qcTable
        self.products = productsTable
        self.recipeName = recipeName
        self.sofName = sofName

        from soxspipe.commonutils.toolkit import get_calibrations_path
        from astropy.table import Table
        from astropy.io import fits

        # CONVERTING EXTRACTION (BACK) TO DATAFRAME
        self.stdExtractionDF = Table.read(
            self.stdExtractionPath, format='fits')
        self.stdExtractionDF = self.stdExtractionDF.to_pandas()

        self.calibrationRootPath = get_calibrations_path(
            log=self.log, settings=self.settings)

        from soxspipe.commonutils import keyword_lookup
        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw

        stdNames = ["LTT7987", "EG274", "LTT3218", "EG21"]
        stdAkas = ["CD-3017706", "CD-3810980", "CD-325613", "CPD-69177"]

        hdul = fits.open(self.stdExtractionPath)
        self.header = hdul[0].header
        self.dateObs = self.header[kw("DATE_OBS")]
        self.texp = float(self.header[kw("EXPTIME")])
        if self.instrument == 'xsh':
            # Name is in the format 'EG 274'
            self.std_objName = self.header[kw("OBS_TARG_NAME")].strip().upper()
        else:
            # Name is in the format 'EG 274'
            self.std_objName = self.header[kw("OBJECT")].strip().upper()

            if "STD," in self.std_objName:
                try:
                    self.std_objName = self.header[kw(
                        "OBS_TARG_NAME")].strip().upper()
                except:
                    pass

        self.std_objName = self.std_objName.split(
            " V")[0].replace(" ", "")  # Hack to reduce xsh data

        #REMOVE SPACES IN NAME 
        self.std_objName = self.std_objName.replace(" ", "")

        
        self.std_objName = self.std_objName.replace("_NOD", "")


        if stdNotFlatExtractionPath and len(stdNotFlatExtractionPath) > 1:
            # STD STAR GIVEN, READING THE NON FLAT FIELDED SPECTRUM
            self.stdExtractionNotFlatDF = Table.read(
                stdNotFlatExtractionPath, format='fits')
            self.stdExtractionNotFlatDF = self.stdExtractionNotFlatDF.to_pandas()
        if self.std_objName in stdAkas:
            for s, a in zip(stdNames, stdAkas):
                if self.std_objName == a:
                    self.std_objName = s
                

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
        self.arm = self.header[kw("SEQ_ARM")].strip().upper()  # KW lookup

        # DETECTOR PARAMETERS LOOKUP OBJECT
        from soxspipe.commonutils import detector_lookup
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        from soxspipe.commonutils.toolkit import utility_setup
        self.qcDir, self.productDir = utility_setup(
            log=self.log, settings=settings, recipeName=recipeName, startNightDate=startNightDate)

        return None

    def get(self):
        """
        *get the response_function object*

        **Return:**
            - ``response_function`` -- a set of polynomial coefficients
        """
        self.log.debug('starting the ``get`` method')

        import pandas as pd
        from scipy.interpolate import interp1d
        import numpy as np
        from scipy.signal import savgol_filter

        from astropy.table import Table

        from matplotlib import pyplot as plt
        # SWITCH BACKEDN TO MACOSX

        response_function = None

        # GET THE EXTRACTED STANDARD STAR'S WAVELENGTH AND FLUX
        stdExtWave = self.stdExtractionDF['WAVE'].values
        stdExtFlux = self.stdExtractionDF['FLUX_COUNTS'].values
        stdExtWaveNotFlat = self.stdExtractionNotFlatDF['WAVE'].values
        stdExtFluxNotFlat = self.stdExtractionNotFlatDF['FLUX_DENSITY_COUNTS'].values

        # GET THE ABSOLUTE STANDARD STAR FLUXES, ASSUMING TO HAVE 1-1 MAPPING BETWEEN OBJECT NAME IN THE FITS HEADER AND DATABASE
        stdAbsFluxDF = Table.read(
            self.calibrationRootPath + "/" + self.detectorParams["flux-standards"], format='fits')
        stdAbsFluxDF = stdAbsFluxDF.to_pandas()
        # MAKE ALL COLUMNS UPPERCASE
        stdAbsFluxDF.columns = [d.upper() for d in stdAbsFluxDF.columns]

        # SELECTING ROWS IN THE INTERESTED WAVELENGTH RANGE ADDING A MARGIN TO THE RANGE
        stdAbsFluxDF = stdAbsFluxDF[(stdAbsFluxDF['WAVE'] > np.min(
            stdExtWaveNotFlat) - 10) & (stdAbsFluxDF['WAVE'] < 10 + np.max(stdExtWaveNotFlat))]

        # FLUX IS CONVERTED IN ERG / CM2 / S / ANG
        try:
            self.std_wavelength_to_abs_flux = interp1d(np.array(stdAbsFluxDF['WAVE']), np.array(
                stdAbsFluxDF[self.std_objName]) * 10**17, kind='next', fill_value="extrapolate")
            # print("YESSTD")
        except Exception as e:
            # print("NOSTD")
            self.log.warning(
                f"Standard star {self.std_objName} not found in the static calibration database. The available STDs are {', '.join(stdAbsFluxDF.columns[1:])}")
            return self.qc, self.products

        # STRONG SKY ABS REGION TO BE EXCLUDED
        if self.arm == 'NIR':
            excludeRegions = [(770,850), (930, 970), (1080, 1200), (1280, 1520), (1785, 1950)]
        else:
            excludeRegions = []
        # INTEGRATING THE EXTRACTED STANDARD IN xx nm BIN-WIDE FILTERS (CONVERTING BACK IN A/PX)
        stdExtFluxNotFlat = stdExtFluxNotFlat / 10.

        # NOW DIVIDING FOR THE EXPOSURE TIME
        stdExtFlux = stdExtFlux / self.texp
        stdExtFluxNotFlat = stdExtFluxNotFlat / self.texp

        if False:
            plt.plot(stdExtWave, stdExtFlux)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Extracted Flux (counts/s)')
            plt.title('Extracted Standard Star Spectrum')
            plt.grid()
            plt.show()


            plt.plot(stdAbsFluxDF['WAVE'], stdAbsFluxDF[self.std_objName]*10**17)
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Absolute Flux (erg/s/cm²/Å)')
            plt.title('Absolute Standard Star Spectrum')
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
            if self.instrument == 'xsh':
                # VLT TELESCOPE AREA (IN CM²)
                area = 400 * 400 * 3.14
            else:
                # NTT TELESCOPE AREA (IN CM²)
                area = 179 * 179 * 3.14

            # EXTRACT WAVELENGTH VALUES AND COMPUTE ABSOLUTE FLUX
            wavelength = stdExtWaveNotFlat
            abs_flux = self.std_wavelength_to_abs_flux(wavelength) * 10**-17
            # CONVERT WAVELENGTH FROM ÅNGSTRÖM TO CM
            wavelength_in_cm = wavelength * 10**-7

            # COMPUTE THE ABSOLUTE PHOTON FLUX OF THE STANDARD STAR
            hc = h * c  # PRECOMPUTE PLANCK'S CONSTANT * SPEED OF LIGHT
            stdAbsPhotonFlux = (abs_flux * wavelength_in_cm / hc) * area

            # COMPUTE THE EFFICIENCY ESTIMATE
            stdEfficiencyEstimate = stdExtFluxNotFlat / stdAbsPhotonFlux
            # Change backend to MacOSX
            if False:
                import matplotlib
                matplotlib.use('MacOSX')
                from matplotlib import pyplot as plt
                print(self.texp)
                plt.plot(wavelength,stdExtFluxNotFlat, label='Extracted Flux (not flat corrected)')
                plt.plot(wavelength,stdAbsPhotonFlux, label='Absolute Photon Flux')
                plt.plot(wavelength,stdEfficiencyEstimate, label='Efficiency Estimate')
                plt.xlabel('Wavelength (nm)')
                plt.ylabel('Efficiency')
                plt.title('Efficiency Estimate')
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
        if self.arm == 'UVB' or self.arm == 'VIS':
            self.log.print('Applying extinction correction')
            extCorrectionFactor = extinction_correction_factor(
                stdExtWave, self.calibrationRootPath + "/" + self.detectorParams["extinction"], self.airmass)
            if False:
                from matplotlib import pyplot as plt
                plt.plot(stdExtWave, extCorrectionFactor)
                plt.xlabel('Wavelength (nm)')
                plt.ylabel('Extinction Correction Factor')
                plt.title('Extinction Correction Factor vs Wavelength')
                plt.grid()
                plt.show()
            stdExtFlux = stdExtFlux * extCorrectionFactor

        self.response_function_raw = self.std_wavelength_to_abs_flux(
            stdExtWave) / stdExtFlux
        wavelength_response = self.stdExtractionDF['WAVE'].values

        polyOrder = int(self.recipeSettings[self.arm.lower()]['poly_order'])
        numIter = 0
        deletedPoints = 1

        # REMOVE EXLUCED REGION FROM WAVELENGTH AND RESPONSE FUNCTION
        for er in excludeRegions:
            elementsToDelete = np.where(
                (wavelength_response >= er[0]) & (wavelength_response <= er[1]))[0]
            wavelength_response = np.delete(
                wavelength_response, elementsToDelete)
            self.response_function_raw = np.delete(
                self.response_function_raw, elementsToDelete)
            
        #SMOOTHING DATA IF NIR
        if self.arm == 'NIR':
            from scipy.ndimage import gaussian_filter1d
            self.response_function_raw = gaussian_filter1d(self.response_function_raw, sigma=5)
        
        # FITTING ITERATIVELY THE DATA WITH A POLYNOMIAL

        while (numIter < int(self.recipeSettings[self.arm.lower()]['max_iteration'])) and (deletedPoints > 0):
            try:
                # FITTING THE DATA
                elementsToDelete = []
                responseFuncCoeffs = np.polyfit(
                    wavelength_response, self.response_function_raw, polyOrder)
                for index, (w, z, zf) in enumerate(zip(wavelength_response, self.response_function_raw, np.polyval(responseFuncCoeffs, wavelength_response))):
                    # if np.abs(np.abs(z)-np.abs(zf)) > 0.05:
                    #ff np.abs(np.abs(z) - np.abs(zf)) / np.abs(z) > 100 or z < 0:

                    if z <0 or np.abs(np.abs(z) - np.abs(zf)) / np.abs(z) > 0.2:
                        elementsToDelete.append(index)

                wavelength_response = np.delete(
                    wavelength_response, elementsToDelete)
                self.response_function_raw = np.delete(
                    self.response_function_raw, elementsToDelete)
                deletedPoints = len(elementsToDelete)
                numIter = numIter + 1
            except Exception as e:
                raise Exception(
                    'The fitting of response function did not converge!')
                sys.exit(1)
        if False:
            plt.plot(wavelength_response, self.response_function_raw)
            plt.show()
        # WRITE RESPONSE FUNCTION TO FITS BINARY TABLE
        self.write_response_function_to_file(
            responseFuncCoeffs=responseFuncCoeffs,
            polyOrder=polyOrder
        )
        print("RESP FUNC WRITTEN")

        if not isinstance(stdEfficiencyEstimate, bool):
            # CREATE A DATAFRAME FOR EFFICIENCY ESTIMATE
            stdEfficiencyEstimateDF = pd.DataFrame({
                'WAVE': stdExtWaveNotFlat,
                'EFFICIENCY': stdEfficiencyEstimate
            })
            # WRITE THE EFFICIENCY ESTIMATE TO FITS BINARY TABLE
            from astropy.table import Table
            import copy
            from astropy.io import fits
            t = Table.from_pandas(stdEfficiencyEstimateDF)
            filename = f"{self.sofName}_EFFICIENCY.fits"
            filepath = f"{self.productDir}/{filename}"

            header = copy.deepcopy(self.header)
            header[self.kw("SEQ_ARM").upper()] = self.arm
            header[self.kw("PRO_TYPE").upper()] = "REDUCED"
            header[self.kw("PRO_CATG").upper()
                   ] = f"EFFICIENCY_TAB_{self.arm}".upper()

            BinTableHDU = fits.table_to_hdu(t)

            priHDU = fits.PrimaryHDU(header=header)

            hduList = fits.HDUList([priHDU, BinTableHDU])
            hduList.writeto(filepath, checksum=True, overwrite=True)

            # ADD TO QC TABLE
            from datetime import datetime
            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": "EFFICIENCY",
                "file_name": filename,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"SOXS efficiency estimate",
                "file_path": filepath,
                "label": "QC"
            }).to_frame().T], ignore_index=True)

        self.plot_response_curve(
            stdExtWave=stdExtWave,
            stdExtWaveNotFlat=stdExtWaveNotFlat,
            stdExtFlux=stdExtFlux,
            binCentreWave=wavelength_response,
            binCentreWaveOriginal=wavelength_response,
            binIntegratedFlux=self.response_function_raw,
            absToExtFluxRatio=None,
            responseFuncCoeffs=responseFuncCoeffs,
            stdEfficiencyEstimate=stdEfficiencyEstimate
        )

        return self.qc, self.products

    def plot_response_curve(
            self,
            stdExtWave,
            stdExtWaveNotFlat,
            stdExtFlux,
            binCentreWave,
            binCentreWaveOriginal,
            binIntegratedFlux,
            absToExtFluxRatio,
            responseFuncCoeffs,
            stdEfficiencyEstimate):
        """*generate a QC plot for the response curve*

        **Key Arguments:**

        - ``stdExtWave`` -- the extracted standard star wavelength
        - ``stdExtFlux`` -- the extracted standard star flux
        - ``stdExtWaveNotFlat`` -- the extracted standard star wavelength (not flattened)
        - ``binCentreWave`` -- binned wavelengths after clipping (during fitting)
        - ``binCentreWaveOriginal`` -- binned wavelengths
        - ``binIntegratedFlux`` -- binned flux
        - ``absToExtFluxRatio`` -- the ratio of the absolute flux vs the extraction flux
        - ``responseFuncCoeffs`` -- the response function coefficients
        - ``stdEfficiencyEstimate`` -- the estimated instrument efficiency

        **Return:**
            - `plotFilePath` -- the path to the QC plot PDF
        """
        self.log.debug('starting the ``plot_response_curve`` method')

        import matplotlib.pyplot as plt
        import numpy as np
        from datetime import datetime
        import pandas as pd

        # WRITE THE QC PLOT TO PDF
        fig = plt.figure(figsize=(6, 22), constrained_layout=True, dpi=120)
        if stdEfficiencyEstimate is not None:
            gs = fig.add_gridspec(30, 4)
        else:
            gs = fig.add_gridspec(25, 4)
        onerow = fig.add_subplot(gs[0:4, :])
        tworow = fig.add_subplot(gs[5:9, :])
        threerow = fig.add_subplot(gs[10:14, :])
        fourrow = fig.add_subplot(gs[15:19, :])
        fiverow = fig.add_subplot(gs[20:24, :])
        if stdEfficiencyEstimate is not None:
            sixrow = fig.add_subplot(gs[25:29, :])

        onerow.plot(stdExtWave, self.std_wavelength_to_abs_flux(
            stdExtWave) * 10**-17, linewidth=0.2)
        onerow.set_title(
            f'{self.std_objName} absolute flux spectrum', fontsize=12)
        onerow.set_xlabel(f"wavelength (nm)", fontsize=9)
        onerow.set_ylabel(
            "flux ($\\mathrm{erg/cm^{2}/s/angstom}$)", fontsize=9)
        onerow.tick_params(axis='both', which='major', labelsize=9)
        # Set y-limits based on the absolute flux spectrum
        abs_flux = self.std_wavelength_to_abs_flux(stdExtWave)
        min_flux, max_flux = np.min(abs_flux), np.max(abs_flux)
        flux_margin = 0.2 * (max_flux - min_flux)  # Add 10% margin
        # onerow.set_ylim(min_flux - flux_margin, max_flux + flux_margin)

        tworow.scatter(binCentreWaveOriginal, binIntegratedFlux,
                       marker='o', s=10, alpha=0.5)
        tworow.set_title('Raw ratio', fontsize=12)
        tworow.set_xlabel(f"wavelength (nm)", fontsize=9)
        tworow.set_ylabel("Ratio $\\frac{F_{\\lambda}}{F_c}$", fontsize=9)
        tworow.tick_params(axis='both', which='major', labelsize=9)

        # threerow.plot(binCentreWave, absToExtFluxRatio, linewidth=0.5)
        threerow.set_title('Fitted Response Function', fontsize=12)
        threerow.plot(stdExtWave, np.polyval(responseFuncCoeffs,
                      stdExtWave), c='red', label="response curve", linewidth=2.5, alpha=0.7)
        threerow.scatter(binCentreWaveOriginal, binIntegratedFlux,
                         marker='o', s=10, alpha=0.2)
        # threerow.set_xlim(min(binCentreWave), max(binCentreWave))
        # threerow.set_ylim(min(absToExtFluxRatio), max(absToExtFluxRatio))
        threerow.set_xlabel(f"wavelength (nm)", fontsize=9)
        threerow.set_ylabel("absolute-extracted flux ratio", fontsize=9)
        threerow.tick_params(axis='both', which='major', labelsize=9)

        # TEST THE X-CALIB

        # FLUX IS ALREADY DIVIDED BY DISPERSION AND CORRECTED FOR THE EXTINCTION !!
        # OTHERWISE, COPY LINES ABOVE
        flux_calib = stdExtFlux * \
            np.polyval(responseFuncCoeffs, stdExtWave)
        flux_calib = flux_calib * 10**-17  # CONVERTING BACK TO PHYS UNITS
        fourrow.plot(stdExtWave, flux_calib, linewidth=0.2)
        fourrow.set_title('Self calibration of std star', fontsize=12)
        fourrow.set_xlabel(f"wavelength (nm)", fontsize=9)
        fourrow.set_ylabel(
            "flux ($\\mathrm{erg/cm^{2}/s/angstom}$)", fontsize=9)
        fourrow.tick_params(axis='both', which='major', labelsize=9)
        # fourrow.set_ylim(min_flux - flux_margin, max_flux + flux_margin)

        fiverow.plot(stdExtWave, (flux_calib * 10**17 - self.std_wavelength_to_abs_flux(
            stdExtWave)) / self.std_wavelength_to_abs_flux(stdExtWave), linewidth=0.2)
        fiverow.set_ylim(-5, 5)
        # plt.plot(np.array(stdAbsFluxDF[0]),np.array(stdAbsFluxDF[4])*10**17,c='red')
        plt.subplots_adjust(hspace=1.0)
        fiverow.set_title('Relative residuals', fontsize=12)
        fiverow.set_xlabel(f"wavelength (nm)", fontsize=9)
        fiverow.set_ylabel("residual", fontsize=9)
        fiverow.tick_params(axis='both', which='major', labelsize=9)

        if stdEfficiencyEstimate is not None:
            sixrow.plot(stdExtWaveNotFlat,
                        stdEfficiencyEstimate, linewidth=0.2)

            sixrow.set_ylim(0, np.min([np.nanmax(stdEfficiencyEstimate), 1.0]))
            # plt.plot(np.array(stdAbsFluxDF[0]),np.array(stdAbsFluxDF[4])*10**17,c='red')
            plt.subplots_adjust(hspace=1.0)
            sixrow.set_title('Efficiency (end-to-end)', fontsize=12)
            sixrow.set_xlabel(f"wavelength (nm)", fontsize=9)
            sixrow.set_ylabel("Efficiency", fontsize=9)
            sixrow.tick_params(axis='both', which='major', labelsize=9)

        plotFilename = self.sofName + "_RESPONSE.pdf"
        plotFilePath = f"{self.qcDir}/{plotFilename}"
        plt.savefig(plotFilePath, dpi=120, bbox_inches='tight')
        plt.close('all')

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "RESPONSE_QC_PLOT",
            "file_name": plotFilename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"Response curve QC plot.",
            "file_path": plotFilePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        if False:
            plt.show()

        self.log.debug('completed the ``plot_response_curve`` method')
        return plotFilePath

    def write_response_function_to_file(
            self,
            responseFuncCoeffs,
            polyOrder):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**

        - ``responseFuncCoeffs`` -- the response curve coefficients
        - ``polyOrder`` -- the order of polynomial used to fit the curve

        **Return:**

        - ``responseCurvePath`` -- path to the saved file
        """
        self.log.debug(
            'starting the ``write_response_function_to_file`` method')

        import pandas as pd
        from astropy.table import Table
        from astropy.io import fits
        from datetime import datetime
        import copy
        arm = self.arm
        kw = self.kw

        coeffDict = {}
        coeffDict["polyOrder"] = polyOrder
        n_coeff = 0
        for i in range(0, polyOrder + 1):
            coeffDict[f'c{i}'] = responseFuncCoeffs[n_coeff]
            n_coeff += 1

        filename = self.sofName + "_RESP.fits"

        header = copy.deepcopy(self.header)

        filePath = f"{self.productDir}/{filename}"
        df = pd.DataFrame([coeffDict])
        t = Table.from_pandas(df)

        BinTableHDU = fits.table_to_hdu(t)

        header[kw("SEQ_ARM").upper()] = arm
        header[kw("PRO_TYPE").upper()] = "REDUCED"
        header[kw("PRO_CATG").upper()] = f"RESP_TAB_{arm}".upper()

        # WRITE QCs TO HEADERS
        for n, v, c, h in zip(self.qc["qc_name"].values, self.qc["qc_value"].values, self.qc["qc_comment"].values, self.qc["to_header"].values):
            if h:
                header[f"ESO QC {n}".upper()] = (v, c)

        priHDU = fits.PrimaryHDU(header=header)

        hduList = fits.HDUList([priHDU, BinTableHDU])
        hduList.writeto(filePath, checksum=True, overwrite=True)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "RESPONSE_FUNC",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"Response function coeffs.",
            "file_path": filePath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

        self.log.debug(
            'completed the ``write_response_function_to_file`` method')
        return filePath
