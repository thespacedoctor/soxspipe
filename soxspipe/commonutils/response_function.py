#!/usr/bin/env python
# encoding: utf-8
"""
*Given a standard star extracted spectrum, generate the instrument response function needed to flux calibrate science spectra*

:Author:
    David Young

:Date Created:
    July 28, 2023
"""
import os
from builtins import object
from matplotlib.pyplot import figure
os.environ['TERM'] = 'vt100'
import sys

class response_function(object):
    """
    *The worker class for the response_function module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``stdExtractionPath`` -- fits binary table containing the extracted standard spectrum

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a response_function object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``response_function`` to documentation
        - create a blog post about what ``response_function`` does
    ```

    ```python
    usage code 
    ```

    """
    # Initialisation
    # 1. @flagged: what are the unique attrributes for each object? Add them
    # to __init__

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
        from astropy.table import Table
        from astropy.io import fits


        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # Initial Actions
        # OPEN EXTRACTED SPECTRUM
        # SPEC FORMAT TO PANDAS DATAFRAME

        self.stdExtractionDF = Table.read(self.stdExtractionPath, format='fits')

        self.stdExtractionDF = self.stdExtractionDF.to_pandas()

        #from tabulate import tabulate
        #print(tabulate(stdExtractionDF.head(100), headers='keys', tablefmt='psql'))

        # TODO READING NECESSARY KEYWORDS FROM THE HEADER
        hdul = fits.open(self.stdExtractionPath)
        self.texp = 1.0
        self.std_objName = 'EG 274'
        self.airmass = 1.039
        self.arm = 'UVB'
        #TEXP, OBJECT NAME



    def extinction_correction_factor(self, wave):
        from scipy.interpolate import interp1d
        import numpy as np
        import matplotlib.pyplot as plt
        import pandas as pd

        #READ THE EXTINCTION CURVE FOR THE OBSERVATORY

        extinctionData = pd.read_csv('/Users/mlandoni/Desktop/extinction_lasilla.dat', sep='\t', header=0)

        #DATA IS ORGANIZED AS FOLLOWS:
        #FIRST COLUMN, WAVELENGTH (IN ANGSTROM), SECOND COLUMN MAG/AIRMASS

        wave_ext = extinctionData['WAVE']/10

        #INTERPOLATING ON THE REQUIRED WAVE SCALE

        refitted_ext = interp1d(np.array(wave_ext),
                                 np.array(extinctionData['MAG_AIRMASS']), kind='next',fill_value='extrapolate')
        if False:
            fig, ax = plt.subplots(2)
            ax[0].plot(wave, refitted_ext(wave))
            ax[0].plot(wave_ext, extinctionData['MAG_AIRMASS'])
            ax[1].plot(wave_ext,refitted_ext(wave_ext))
            plt.show()

        return 10**(0.4*refitted_ext(wave)*self.airmass)



        pass
    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def get(self):
        """
        *get the response_function object*

        **Return:**
            - ``response_function``

        **Usage:**

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed
        ```

        ```python
        usage code 
        ```
        """
        self.log.debug('starting the ``get`` method')

        import pandas as pd
        from scipy.interpolate import interp1d
        import numpy as np
        import matplotlib.pyplot as plt

        response_function = None

        wave = np.array(self.stdExtractionDF['WAVE'])
        flux = np.array(self.stdExtractionDF['FLUX_COUNTS'])

        self.extinction_correction_factor(wave)

        # TODO: ADD CORRECT PATH
        # READING THE DATA FROM THE DATABASE, ASSUMING TO HAVE 1-1 MAPPING BETWEEN OBJECT NAME IN THE FITS HEADER AND DATABASE
        stdData = pd.read_csv('/Users/mlandoni/Desktop/std.dat',sep=' ',header=0)



        # SELECTING ROWS IN THE INTERESTED WAVELENGTH RANGE ADDING A MARGIN TO THE RANGE
        selected_rows = stdData[(stdData['WAVE'] > np.min(self.stdExtractionDF['WAVE']) -10 ) & (stdData['WAVE'] < 10 + np.max(self.stdExtractionDF['WAVE']))]

        print(selected_rows['WAVE'])
        #FLUX IS CONVERTED IN ERG / CM2 / S / ANG
        try:
            refitted_flux = interp1d(np.array(selected_rows['WAVE']),np.array(selected_rows[self.std_objName])*10*10**17,kind='next')
        except Exception as e:
            raise Exception("Standard star %s not found in the static calibration database" % self.std_objName)
            sys.exit(1)


        # STRONG SKY ABS REGION TO BE EXCLUDED - 1 TABLE PER ARM
        #TODO Check on FITS Header for the correct arm
        #exclude_regions = [(1100 ,1190), (1300, 1500), (1800, 1900), (1850,2700)]

        exclude_regions = [(200, 400), (590,600), (405, 416), (426,440),(460,475),  (563, 574), (478,495), (528,538)]




        #INTEGRATING THE OBS SPECTRUM IN xx nm BIN-WIDE FILTERS

        # CONVERTING FLUX_COUNTS IN FLUX_COUNTS PER NM DIVIDING BY THE PIXEL SIZE
        dispersion = wave - np.roll(wave, 1)
        dispersion[0] = dispersion[1]  # first element will be out

        flux = flux / dispersion

        #NOW DIVIDING FOR THE EXPOSURE TIME
        flux = flux/self.texp

        if self.arm == 'UVB' or self.arm == 'VIS':
            print('Applying extinction correction')
            factor = self.extinction_correction_factor(wave)
            flux = flux * factor

        bin_width = 3
        bin_starts = np.arange(min(self.stdExtractionDF['WAVE']), max(self.stdExtractionDF['WAVE']), bin_width)
        bin_ends = bin_starts + bin_width

        central_wavelengths = []
        integrated_flux = []




        for bin_start, bin_end in zip(bin_starts, bin_ends):
            central_wave = (bin_start + bin_end) / 2
            exclude_bin = any(start <= bin_start <= end or start <= bin_end <= end for start, end in exclude_regions)

            if not exclude_bin:
                mask = (wave >= bin_start) & (wave < bin_end)
                bin_flux = flux[mask]
                bin_integral = np.trapz(bin_flux, wave[mask])/(bin_end - bin_start)
                if not np.isnan(bin_integral) and bin_integral > 0:
                    central_wavelengths.append(central_wave)
                    integrated_flux.append(bin_integral)

        # NOW FINDING THE RESPONSE FUNCTION POINTS AND THEN FIT


        figure(figsize = (8, 12))
        fig, axs = plt.subplots(5)
        print(wave)

        axs[0].plot(wave,refitted_flux(wave))
        axs[0].set_title('Tabulated flux')
        axs[1].scatter(central_wavelengths,integrated_flux)
        axs[1].set_title('Passband Photometry')


        # FINDING THE FUNCTION S = F/C
        zp = np.array(refitted_flux(central_wavelengths)/integrated_flux)
        central_wavelengths = np.array(central_wavelengths)

        axs[2].plot(central_wavelengths, zp)
        axs[2].set_title('Response function points vs fit')




        zp_original = zp
        central_wavelengths_original = central_wavelengths
        numIter = 0
        deletedPoints = 1
        order = 7


        #FITTING ITERATIVELY THE DATA WITH A POLYNOMIAL
        while (numIter < 5) and (deletedPoints > 0):
            try:
                #FITTING THE DATA
                elements_to_delete = []
                coefficients = np.polyfit(central_wavelengths, zp, order)
                for index, (w,z, zf) in enumerate(zip(central_wavelengths, zp, np.polyval(coefficients, central_wavelengths ))):
                    #if np.abs(np.abs(z)-np.abs(zf)) > 0.05:
                    if np.abs(np.abs(z) - np.abs(zf))/np.abs(z) > 0.1:
                        elements_to_delete.append(index)

                central_wavelengths = np.delete(central_wavelengths, elements_to_delete)
                zp = np.delete(zp, elements_to_delete)
                deletedPoints = len(elements_to_delete)
                numIter = numIter + 1
            except Exception as e:
                raise Exception('The fitting of response function did not converge!')
                sys.exit(1)
        axs[2].plot(wave, np.polyval(coefficients,wave),c='red')
        axs[2].set_xlim(min(central_wavelengths), max(central_wavelengths))
        axs[2].set_ylim(min(zp), max(zp))



        response_function = coefficients
        self.log.debug('completed the ``get`` method')


        #TEST THE X-CALIB

        #FLUX IS ALREADY DIVIDED BY DISPERSION AND CORRECTED FOR THE EXTINCTION !!
        #OTHERWISE, COPY LINES 192-203
        flux_calib = flux*np.polyval(response_function, self.stdExtractionDF['WAVE'])


        axs[3].plot(self.stdExtractionDF['WAVE'], flux_calib)
        axs[3].set_title('Self calibration of std star')
        #axs[3].set_xlim(0,np.max(flux_calib))

        axs[4].plot(self.stdExtractionDF['WAVE'], (flux_calib - refitted_flux(self.stdExtractionDF['WAVE']))/refitted_flux(self.stdExtractionDF['WAVE']))
        axs[4].set_ylim(-5, 5)
        #plt.plot(np.array(selected_rows[0]),np.array(selected_rows[4])*10**17,c='red')
        plt.subplots_adjust(hspace=1.0)
        axs[4].set_title('Relative residuals')
        plt.show()

        return response_function

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
