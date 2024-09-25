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

os.environ['TERM'] = 'vt100'


class response_function(object):
    """
    *Given a standard star extracted spectrum, generate the instrument response function needed to flux calibrate science spectra*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``stdExtractionPath`` -- fits binary table containing the extracted standard spectrum
        - ``settings`` --  the pipeline settings

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a response_function object, use the following:


    ```python
    from soxspipe.commonutils import response_function
    response = response_function(
        log=log,
        settings=settings,
        stdExtractionPath=stdExtractionPath
    )
    response_function = response.get()
    ```
    """

    def __init__(
            self,
            log,
            stdExtractionPath,
            settings=False,
    ):
        self.log = log
        log.debug("instansiating a new 'response_function' object")
        self.settings = settings
        self.stdExtractionPath = stdExtractionPath
        self.recipeSettings = settings['soxs-response']
        self.instrument = self.settings["instrument"].lower()

        from soxspipe.commonutils.toolkit import get_calibrations_path
        from astropy.table import Table
        from astropy.io import fits

        # CONVERTING EXTRACTION TO DATAFRAME
        self.stdExtractionDF = Table.read(self.stdExtractionPath, format='fits')
        self.stdExtractionDF = self.stdExtractionDF.to_pandas()

        self.calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)

        # DAVE: Where should response_function be called? I think in stare and nodding recipes whenever reducing a standard star

        # DAVE: Why would there an issue reading the header?
        # TODO: Use the keyword dictionary to grab header keywords
        try:
            hdul = fits.open(self.stdExtractionPath)
            header = hdul[0].header
            self.texp = float(header['EXPTIME'])
            self.std_objName = str(header['OBJECT']).strip().upper()  # Name is in the format 'EG 274'
            # USING THE AVERAGE AIR MASS
            airmass_start = float(header['HIERARCH ESO TEL AIRM START'])
            airmass_end = float(header['HIERARCH ESO TEL AIRM END'])
            self.airmass = (airmass_start + airmass_end) / 2
            self.arm = str(header['HIERARCH ESO SEQ ARM']).strip().upper()  # KW lookup
        except Exception as e:
            raise Exception('Error reading the FITS header')
            sys.exit(1)

        return None

    def extinction_correction_factor(
            self,
            wave):

        from scipy.interpolate import interp1d
        import numpy as np
        import matplotlib.pyplot as plt
        import pandas as pd

        # DAVE: DO WE NEED PARANAL FOR XSH?
        # TODO: CONVERT extinction_lasilla.dat TO FITS BINARY TABLE

        # READ THE EXTINCTION CURVE FOR THE OBSERVATORY
        # DATA IS ORGANIZED AS FOLLOWS:
        # FIRST COLUMN, WAVELENGTH (IN ANGSTROM), SECOND COLUMN MAG/AIRMASS

        if self.instrument == "soxs":
            extinctionData = pd.read_csv(self.calibrationRootPath + '/extinction_lasilla.dat', sep='\t', header=0)
        else:
            extinctionData = pd.read_csv(self.calibrationRootPath + '/extinction_lasilla.dat', sep='\t', header=0)

        # CONVERT ANG TO NM
        wave_ext = extinctionData['WAVE'] / 10
        # INTERPOLATING ON THE REQUIRED WAVE SCALE
        refitted_ext = interp1d(np.array(wave_ext),
                                np.array(extinctionData['MAG_AIRMASS']), kind='next', fill_value='extrapolate')
        if False:
            fig, ax = plt.subplots(3)
            ax[0].plot(wave, refitted_ext(wave))
            ax[0].plot(wave_ext, extinctionData['MAG_AIRMASS'])
            ax[1].plot(wave_ext, refitted_ext(wave_ext))
            ax[2].plot(wave, 10**(0.4 * refitted_ext(wave) * self.airmass))
            plt.title('Extinction Correction factor')
            plt.show()

        # DAVE: IS THIS CONVERTING TO FLUX SPACE?
        return 10**(0.4 * refitted_ext(wave) * self.airmass)

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
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import figure

        response_function = None

        wave = self.stdExtractionDF['WAVE'].values
        flux = self.stdExtractionDF['FLUX_COUNTS'].values

        # TODO: CONVERT std.dat TO FITS BINARY TABLE

        # READING THE DATA FROM THE DATABASE, ASSUMING TO HAVE 1-1 MAPPING BETWEEN OBJECT NAME IN THE FITS HEADER AND DATABASE
        stdData = pd.read_csv(self.calibrationRootPath + '/std.dat', sep=' ', header=0)

        # SELECTING ROWS IN THE INTERESTED WAVELENGTH RANGE ADDING A MARGIN TO THE RANGE
        selected_rows = stdData[(stdData['WAVE'] > np.min(self.stdExtractionDF['WAVE']) - 10) & (stdData['WAVE'] < 10 + np.max(self.stdExtractionDF['WAVE']))]

        # TODO: ADD A SEPARATE CHECK AND LOGGING FOR std_objName

        # FLUX IS CONVERTED IN ERG / CM2 / S / ANG
        try:
            refitted_flux = interp1d(np.array(selected_rows['WAVE']), np.array(selected_rows[self.std_objName]) * 10 * 10**17, kind='next')
        except Exception as e:
            raise Exception("Standard star %s not found in the static calibration database" % self.std_objName)
            sys.exit(1)

        # STRONG SKY ABS REGION TO BE EXCLUDED

        # TODO: ADD FOR VIS
        # TODO: ADD THIS EXCLUDED REGIONS TO STATIC CALIBRATION ... SAME FOR XSH AND SOXS
        # MIGHT BE BEST TO HAVE A SINGLE LIST ... NO NEED TO DIFFERENTIATE SOXS VS XSH ARMS

        if self.arm == 'UVB':
            exclude_regions = [(200, 400), (590, 600), (405, 416), (426, 440), (460, 475), (563, 574), (478, 495), (528, 538)]
        elif self.arm == 'VIS':
            exclude_regions = [(620, 640), (648, 666), (754, 770), (800, 810), (836, 845)]
        elif self.arm == 'NIR':
            exclude_regions = [(1100, 1190), (1300, 1500), (1800, 1900), (1850, 2700)]
        else:
            exclude_regions = []

        # INTEGRATING THE OBS SPECTRUM IN xx nm BIN-WIDE FILTERS

        # CONVERTING FLUX_COUNTS IN FLUX_COUNTS PER NM DIVIDING BY THE PIXEL SIZE
        dispersion = wave - np.roll(wave, 1)
        dispersion[0] = dispersion[1]  # first element will be out

        flux = flux / dispersion

        # NOW DIVIDING FOR THE EXPOSURE TIME
        flux = flux / self.texp

        # NOT REQUIRED FOR NIR
        if self.arm == 'UVB' or self.arm == 'VIS':
            # TODO: CONVERT PRINT STATEMENTS TO self.log.print
            print('Applying extinction correction')
            factor = self.extinction_correction_factor(wave)
            flux = flux * factor

        bin_width = 3
        bin_starts = np.arange(self.stdExtractionDF['WAVE'].min(), self.stdExtractionDF['WAVE'].max(), bin_width)
        bin_ends = bin_starts + bin_width

        central_wavelengths = []
        integrated_flux = []

        for bin_start, bin_end in zip(bin_starts, bin_ends):
            central_wave = (bin_start + bin_end) / 2
            # EXCLUDE BRIGHT SKY REGIONS
            exclude_bin = any(start <= bin_start <= end or start <= bin_end <= end for start, end in exclude_regions)

            if not exclude_bin:
                mask = (wave >= bin_start) & (wave < bin_end)
                bin_flux = flux[mask]
                bin_integral = np.trapz(bin_flux, wave[mask]) / (bin_end - bin_start)
                # COLLECT THE WAVELENGTHS AND INTERGRATED FLUX
                if not np.isnan(bin_integral) and bin_integral > 0:
                    central_wavelengths.append(central_wave)
                    integrated_flux.append(bin_integral)

        # NOW FINDING THE RESPONSE FUNCTION POINTS AND THEN FIT

        figure(figsize=(8, 12))
        fig, axs = plt.subplots(5)

        axs[0].plot(wave, refitted_flux(wave))
        axs[0].set_title('Tabulated flux')
        axs[1].scatter(central_wavelengths, integrated_flux)
        axs[1].set_title('Passband Photometry')

        # FINDING THE FUNCTION S = F/C
        zp = np.array(refitted_flux(central_wavelengths) / integrated_flux)
        central_wavelengths = np.array(central_wavelengths)

        axs[2].plot(central_wavelengths, zp)
        axs[2].set_title('Response function points vs fit')

        zp_original = zp
        central_wavelengths_original = central_wavelengths
        numIter = 0
        deletedPoints = 1
        order = int(self.recipeSettings['poly_order'])

        # FITTING ITERATIVELY THE DATA WITH A POLYNOMIAL
        while (numIter < int(self.recipeSettings['max_iteration'])) and (deletedPoints > 0):
            try:
                # FITTING THE DATA
                elements_to_delete = []
                coefficients = np.polyfit(central_wavelengths, zp, order)
                for index, (w, z, zf) in enumerate(zip(central_wavelengths, zp, np.polyval(coefficients, central_wavelengths))):
                    # if np.abs(np.abs(z)-np.abs(zf)) > 0.05:
                    if np.abs(np.abs(z) - np.abs(zf)) / np.abs(z) > 0.1:
                        elements_to_delete.append(index)

                central_wavelengths = np.delete(central_wavelengths, elements_to_delete)
                zp = np.delete(zp, elements_to_delete)
                deletedPoints = len(elements_to_delete)
                numIter = numIter + 1
            except Exception as e:
                raise Exception('The fitting of response function did not converge!')
                sys.exit(1)
        axs[2].plot(wave, np.polyval(coefficients, wave), c='red')
        axs[2].set_xlim(min(central_wavelengths), max(central_wavelengths))
        axs[2].set_ylim(min(zp), max(zp))

        response_function = coefficients

        # TEST THE X-CALIB

        # FLUX IS ALREADY DIVIDED BY DISPERSION AND CORRECTED FOR THE EXTINCTION !!
        # OTHERWISE, COPY LINES 192-203
        flux_calib = flux * np.polyval(response_function, self.stdExtractionDF['WAVE'])

        axs[3].plot(self.stdExtractionDF['WAVE'], flux_calib)
        axs[3].set_title('Self calibration of std star')
        # axs[3].set_xlim(0,np.max(flux_calib))

        axs[4].plot(self.stdExtractionDF['WAVE'], (flux_calib - refitted_flux(self.stdExtractionDF['WAVE'])) / refitted_flux(self.stdExtractionDF['WAVE']))
        axs[4].set_ylim(-5, 5)
        # plt.plot(np.array(selected_rows[0]),np.array(selected_rows[4])*10**17,c='red')
        plt.subplots_adjust(hspace=1.0)
        axs[4].set_title('Relative residuals')
        if True:
            plt.show()

        # TODO: ADD THE PLOT AS A QC PRODUCT
        # TODO: ADD response_function AS A FITS BINARY TABLE PRODUCT
        # TODO: ADD response_function AS A PRODUCT PREDICTION IN DATA-ORGANISER

        # DAVE: SHOULD THE STARE & NODDING FOR SCIENCE OBJECT SOF FILES THEN CONTAIN THIS RESPONSE FUNCTION?

        return response_function
