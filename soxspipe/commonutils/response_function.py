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

os.environ['TERM'] = 'vt100'


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
        print(hdul[1].header)
        self.texp = 1.0
        self.std_objName = ''
        #TEXP, OBJECT NAME




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
        from scipy.interpolate import UnivariateSpline
        import numpy as np
        import matplotlib.pyplot as plt


        response_function = None


        # TODO: read from static library
        stdData = pd.read_csv('/Users/mlandoni/Desktop/std.dat',sep=' ',header=None)



        # SELECTING ROWS IN THE INTERESTED WAVELENGTH RANGE
        selected_rows = stdData[(stdData[0] >= np.min(self.stdExtractionDF['WAVE'])) & (stdData[0] <= np.max(self.stdExtractionDF['WAVE']))]


        # TODO DATA IN THE EXTRACTED FLUX SHALL BE DIVIDED HERE BY THE PIXEL SIZE.
        # TODO THE NUMBER [4] SHALL BE COMPUTED FROM THE OBJECT NAME LIST

        # CONVERT IN AB MAGNITUDE AND FIT TO THE STANDARD STAR FLUX TO A FUNCTION
        spline = UnivariateSpline(selected_rows[0], -2.5*np.log10(selected_rows[4]/ (3.63e-20)) , k=5)
        refitted_flux = spline(selected_rows[0])


        #INTEGRATING THE FLUX IN 10nm BINS

        # STRONG SKY ABS REGION TO BE EXCLUDED - 1 TABLE PER ARM
        #TODO Check on FITS Header for the correct arm
        exclude_regions = [(1100 ,1190), (1300, 1500), (1800, 1900), (1850,2700)]

        #INTEGRATING THE OBS SPECTRUM IN xx nm BIN-WIDE FILTERS
        #TODO DIVIDE BY PIXEL SIZE IN NM
        bin_width = 5
        bin_starts = np.arange(min(self.stdExtractionDF['WAVE']), max(self.stdExtractionDF['WAVE']), bin_width)
        bin_ends = bin_starts + bin_width

        central_wavelengths = []
        integrated_flux = []
        wave = np.array(self.stdExtractionDF['WAVE'])
        flux = np.array(self.stdExtractionDF['FLUX_COUNTS'])
        for bin_start, bin_end in zip(bin_starts, bin_ends):
            central_wave = (bin_start + bin_end) / 2
            exclude_bin = any(start <= bin_start <= end or start <= bin_end <= end for start, end in exclude_regions)

            if not exclude_bin:
                mask = (wave >= bin_start) & (wave < bin_end)
                bin_flux = flux[mask]
                bin_integral = np.trapz(bin_flux, wave[mask])/(bin_end - bin_start)
                if not np.isnan(bin_integral) and bin_integral > 0:
                    central_wavelengths.append(central_wave)
                    integrated_flux.append(bin_integral/bin_width)

        # NOW FINDING THE RESPONSE FUNCTION POINTS AND THEN FIT
        # -2.5log10(counts/texp) + ZP = STDFLUX
        #ZP = STDFLUX + 2.5log10(counts/texp)

        #CONVERTING THE DATA IN AB-MAG

        zp = spline(central_wavelengths) + 2.5*np.log10(integrated_flux)

        #FROM HERE ZP IS THE Y DATA, CENTRAL WAVELENGTH THE X DATA
        zp = np.array(zp)
        central_wavelengths = np.array(central_wavelengths)

        zp_original = zp
        central_wavelengths_original = central_wavelengths
        numIter = 0
        deletedPoints = 1
        order = 4



        while (numIter < 5) and (deletedPoints > 0):
            try:
                #FITTING THE DATA
                elements_to_delete = []
                coefficients = np.polyfit(central_wavelengths, zp, order)
                for index, (w,z, zf) in enumerate(zip(central_wavelengths, zp, np.polyval(coefficients, central_wavelengths ))):
                    #if np.abs(np.abs(z)-np.abs(zf)) > 0.05:
                    if np.abs(np.abs(z) - np.abs(zf))/np.abs(z) > 0.01:
                        elements_to_delete.append(index)

                central_wavelengths = np.delete(central_wavelengths, elements_to_delete)
                zp = np.delete(zp, elements_to_delete)
                deletedPoints = len(elements_to_delete)
                numIter = numIter + 1
            except Exception as e:
                raise Exception('The fitting of response function did not converge!')


        plt.plot(central_wavelengths_original, zp_original, c='blue',alpha=0.4)
        plt.scatter(central_wavelengths_original, np.polyval(coefficients, central_wavelengths_original),c='red')
        plt.show()

        plt.figure()
        wave_range = np.arange(min(self.stdExtractionDF['WAVE']), max(self.stdExtractionDF['WAVE']))
        plt.plot(wave_range, np.polyval(coefficients,wave_range))
        plt.show()

        response_function = coefficients
        self.log.debug('completed the ``get`` method')
        return response_function

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
