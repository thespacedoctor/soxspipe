#!/usr/bin/env python
# encoding: utf-8
"""
*find and fit the continuum in a pinhole flat frame with low-order polynomials. These polynominals are the central loctions of the orders.*

:Author:
    David Young

:Date Created:
    September 10, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils import detector_lookup
from os.path import expanduser
import numpy as np
import unicodecsv as csv
from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial
import matplotlib.pyplot as plt
from astropy.stats import mad_std
from astropy.modeling import models, fitting
from scipy.signal import find_peaks
from random import random
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip, mad_std
from astropy.visualization import hist


class detect_continuum(object):
    """
    *The worker class for the detect_continuum module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``pinholeFlat`` -- calibrationed pinhole flat frame (CCDObject)
        - ``dispersion_map`` -- path to dispersion map csv file containing polynomial fits of the dispersion solution for the frame
        - ``settings`` -- the settings dictionary

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a detect_continuum object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``detect_continuum`` to documentation
        - create a blog post about what ``detect_continuum`` does
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
            pinholeFlat,
            dispersion_map,
            settings=False
    ):
        self.log = log
        log.debug("instansiating a new 'detect_continuum' object")
        self.settings = settings
        self.pinholeFlat = pinholeFlat
        self.dispersion_map = dispersion_map

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        self.arm = pinholeFlat.header[self.kw("SEQ_ARM")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        return None

    def get(self):
        """
        *get the detect_continuum object*

        **Return:**
            - ``detect_continuum``

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

        orderNums, waveLengthMin, waveLengthMax = self.read_spectral_format()
        pixelArrays = self.create_pixel_arrays(
            orderNums,
            waveLengthMin,
            waveLengthMax)

        sliceLength = self.settings[
            "soxs-order-centre"]["slice-length"]
        peakSigmaLimit = self.settings[
            "soxs-order-centre"]["peak-sigma-limit"]

        for order, xy_array in pixelArrays.items():
            found_xy_array = []
            for xy in xy_array:
                found_xy = self.fit_1d_gaussian_to_slice(
                    xy, sliceLength, peakSigmaLimit)
                if found_xy:
                    found_xy_array.append(found_xy)

            percent = 100 * len(found_xy_array) / len(xy_array)
            print(f"{len(found_xy_array)} out of {len(xy_array)} found ({percent:3.0f}%)")

        detect_continuum = None

        self.log.debug('completed the ``get`` method')
        return detect_continuum

    def read_spectral_format(
            self):
        """*read the spectral format table to get some key parameters*

        **Return:**
            - ``orderNums`` -- a list of the order numbers
            - ``waveLengthMin`` -- a list of the maximum wavelengths reached by each order
            - ``waveLengthMax`` -- a list of the minimum wavelengths reached by each order
        """
        self.log.debug('starting the ``read_spectral_format`` method')

        kw = self.kw
        pinholeFlat = self.pinholeFlat
        dp = self.detectorParams

        # READ THE SPECTRAL FORMAT TABLE FILE
        home = expanduser("~")
        calibrationRootPath = self.settings[
            "calibration-data-root"].replace("~", home)
        spectralFormatFile = calibrationRootPath + \
            "/" + dp["spectral format table"]
        spectralFormat = np.genfromtxt(
            spectralFormatFile, delimiter=',', names=True)

        # print(spectralFormat.dtype.names)

        # EXTRACT REQUIRED PARAMETERS
        orderNums = spectralFormat["ORDER"]
        waveLengthMin = spectralFormat["WLMIN"]
        waveLengthMax = spectralFormat["WLMAX"]

        self.log.debug('completed the ``read_spectral_format`` method')
        return orderNums, waveLengthMin, waveLengthMax

    def create_pixel_arrays(
            self,
            orderNums,
            waveLengthMin,
            waveLengthMax):
        """*create a pixel array for the approximate centre of each order*

        **Key Arguments:**
            - ``orderNums`` -- a list of the order numbers
            - ``waveLengthMin`` -- a list of the maximum wavelengths reached by each order
            - ``waveLengthMax`` -- a list of the minimum wavelengths reached by each order

        **Return:**
            - ``pixelArrays`` -- a dictionary of order based pixel arrays
        """
        self.log.debug('starting the ``create_pixel_arrays`` method')

        # READ ORDER SAMPLING RESOLUTION FROM SETTINGS
        sampleCount = self.settings[
            "soxs-order-centre"]["order-sample-count"]

        # CREATE THE WAVELENGTH/ORDER ARRAYS TO BE CONVERTED TO PIXELS
        orderArray = np.asarray([])
        wavelengthArray = np.asarray([])
        for o, wmin, wmax in zip(orderNums, waveLengthMin, waveLengthMax):

            wavelengthArray = np.append(wavelengthArray, np.arange(
                wmin, wmax, (wmax - wmin) / sampleCount))
            orderArray = np.append(orderArray, np.ones(sampleCount) * o)
        order_wave = (orderArray, wavelengthArray)

        # SETUP EMPTY PIXEL ARRAYS
        pixelArrays = {f"o{k:0.0f}": [] for k in orderNums}

        # READ THE FILE
        home = expanduser("~")
        dispersion_map = self.dispersion_map.replace("~", home)

        with open(dispersion_map, 'rb') as csvFile:
            csvReader = csv.DictReader(
                csvFile, dialect='excel', delimiter=',', quotechar='"')
            for row in csvReader:
                axis = row["axis"]
                order_deg = int(row["order-deg"])
                wavelength_deg = int(row["wavelength-deg"])
                coeff = [float(v) for k, v in row.items() if k not in [
                    "axis", "order-deg", "wavelength-deg"]]
                poly = chebyshev_order_wavelength_polynomials(
                    log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg).poly

                if axis == "x":
                    xcoords = poly(order_wave, *coeff)
                if axis == "y":
                    ycoords = poly(order_wave, *coeff)
        csvFile.close()

        # for o, w, x, y in zip(orderArray, wavelengthArray, xcoords, ycoords):
        #     print(o, w, x, y)

        for o, x, y in zip(orderArray, xcoords, ycoords):
            pixelArrays[f"o{o:0.0f}"].append((x, y))

        for order, pixelArray in pixelArrays.items():
            pass

        polyDeg = self.settings[
            "soxs-order-centre"]["poly-deg"]

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        coeff = self.fit_polynomial(
            pixelArray=pixelArray,
            deg=polyDeg
        )

        self.log.debug('completed the ``create_pixel_arrays`` method')
        return pixelArrays

    def fit_1d_gaussian_to_slice(
            self,
            xy,
            sliceLength,
            peakSigmaLimit):
        """*cut a slice from the pinhole flat along the cross-dispersion direction centred on pixel position, fit 1D gaussian and return the peak pixel position*

        **Key Arguments:**
            - ``xy`` -- the x,y pixel coordinate (tuple)
            - ``sliceLength`` -- length of the slice to cut
            - ``peakSigmaLimit`` -- number of simga above median. If peak is not above this level slice will be rejected

        **Return:**
            - ``fit_xy`` -- gaussian fit peak xy position
        """
        self.log.debug('starting the ``fit_1d_gaussian_to_slice`` method')

        # CLIP OUT A SLICE TO INSPECT CENTRED AT POSITION
        halfSlice = sliceLength / 2
        x_fit = xy[0]
        y_fit = xy[1]

        slice = self.pinholeFlat.data[int(y_fit), max(
            0, int(x_fit - halfSlice)):min(2048, int(x_fit + halfSlice))]

        # CHECK THE SLICE POINTS IF NEEDED
        if 1 == 0:
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            plt.show()

        # EVALUATING THE MEAN AND STD-DEV FOR PEAK FINDING - REMOVES SLICE
        # CONTAINING JUST NOISE
        median_r = np.median(slice)
        std_r = mad_std(slice)
        peaks, _ = find_peaks(slice, height=median_r +
                              peakSigmaLimit * std_r, width=1)

        # CHECK PEAK HAS BEEN FOUND
        if peaks is None or len(peaks) <= 0:
            # CHECK THE SLICE POINTS IF NEEDED
            if 1 == 0:
                x = np.arange(0, len(slice))
                plt.figure(figsize=(8, 5))
                plt.plot(x, slice, 'ko')
                plt.xlabel('Position')
                plt.ylabel('Flux')
                plt.show()
            return None

        # FIT THE DATA USING A 1D GAUSSIAN - USING astropy.modeling
        # CENTRE THE GAUSSIAN ON THE PEAK
        g_init = models.Gaussian1D(
            amplitude=1000., mean=peaks[0], stddev=1.)
        # print(f"g_init: {g_init}")
        fit_g = fitting.LevMarLSQFitter()

        # NOW FIT
        g = fit_g(g_init, np.arange(0, len(slice)), slice)
        xpeak = g.mean + max(0, int(x_fit - halfSlice))
        ypeak = y_fit

        # PRINT A FEW PLOTS IF NEEDED - GUASSIAN FIT OVERLAYED
        if 1 == 0 and random() < 0.02:
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            guassx = np.arange(0, max(x), 0.05)
            plt.plot(guassx, g(guassx), label='Gaussian')
            plt.show()

        self.log.debug('completed the ``fit_1d_gaussian_to_slice`` method')
        return (xpeak, ypeak)

    def fit_polynomial(
            self,
            pixelArray,
            deg):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``pixelArray`` -- the array of x,y pixels of measured 1D guassian peak positions
            - ``deg`` -- the degree of the polynomial to fit the order centre trace

        **Return:**
            - ``coeffs`` -- the coefficients of the polynomial fit
        """
        self.log.debug('starting the ``fit_polynomial`` method')

        clippedCount = 1

        poly = chebyshev_xy_polynomial(
            log=self.log, deg=deg).poly

        clippingSigma = self.settings[
            "soxs-order-centre"]["poly-fitting-residual-clipping-sigma"]

        clippingIterationLimit = self.settings[
            "soxs-order-centre"]["clipping-iteration-limit"]

        xcoords = [v[0] for v in pixelArray]
        ycoords = [v[1] for v in pixelArray]

        iteration = 0
        while clippedCount > 0 and iteration < clippingIterationLimit:
            iteration += 1
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLY
            coeff = np.ones((deg + 1))
            coeff, pcov_x = curve_fit(
                poly, xdata=ycoords, ydata=xcoords, p0=coeff)

            residuals, mean_res, std_res, median_res = self.calculate_residuals(
                xcoords=xcoords,
                ycoords=ycoords,
                coeff=coeff,
                deg=deg)

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                residuals, sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc=mad_std)

            # MASK DATA ARRAYS WITH CLIPPED RESIDUAL MASK
            startCount = len(xcoords)
            a = [xcoords, ycoords]
            xcoords, ycoords = [np.ma.compressed(np.ma.masked_array(
                i, masked_residuals.mask)) for i in a]
            clippedCount = startCount - len(xcoords)
            print(f'{clippedCount} pixel positions where clipped in this iteration of fitting an order centre polynomial')

        # PLOT THE RESIDUALS NOW CLIPPING IS COMPLETE
        residuals, mean_res, std_res, median_res = self.calculate_residuals(
            xcoords=xcoords,
            ycoords=ycoords,
            coeff=coeff,
            deg=deg,
            plot=True)

        print(f'\nThe order centre polynomial fitted against the observed 1D gaussian peak positions with a mean residual of {mean_res:2.2f} pixels (stdev = {std_res:2.2f} pixles)')

        self.log.debug('completed the ``fit_polynomials`` method')
        return coeff

    def calculate_residuals(
            self,
            xcoords,
            ycoords,
            coeff,
            deg,
            plot=False):
        """*calculate residuals of the polynomial fits against the observed line postions*

        **Key Arguments:**
            - ``xcoords`` -- the measured x positions of the gaussian peaks
            - ``ycoords`` -- the measurd y positions of the gaussian peaks
            - ``coeff`` -- the coefficients of the fitted polynomial
            - ``deg`` -- degree of the fitted polynomial
            - ``plot`` -- write out a plot to file. Default False.

        **Return:**
            - ``residuals`` -- x residuals
            - ``mean`` -- the mean of the residuals
            - ``std`` -- the stdev of the residuals
            - ``median`` -- the median of the residuals
        """
        self.log.debug('starting the ``calculate_residuals`` method')

        arm = self.arm

        poly = chebyshev_xy_polynomial(
            log=self.log, deg=deg).poly

        # CALCULATE RESIDUALS BETWEEN GAUSSIAN PEAK LINE POSITIONS AND POLY
        # FITTED POSITIONS
        res_x = np.asarray(poly(
            ycoords, *coeff)) - np.asarray(xcoords)

        # CALCULATE COMBINED RESIDUALS AND STATS
        res_mean = np.mean(res_x)
        res_std = np.std(res_x)
        res_median = np.median(res_x)

        # IF PLOT IS REQUIRED:
        if plot:
            fig, ax = plt.subplots(1, 2, figsize=(10, 4))
            plt.subplots_adjust(top=0.85)

            hist(res_x, bins='scott', ax=ax[0], histtype='stepfilled',
                 alpha=0.7, density=True)
            ax[0].set_xlabel('x residuals')
            subtitle = f"mean res: {res_mean:2.3f} pix, res stdev: {res_std:2.3f}"
            fig.suptitle(
                f"residuals of global dispersion solution fitting - single pinhole\n{subtitle}")

            home = expanduser("~")
            outDir = self.settings["intermediate-data-root"].replace("~", home)
            filePath = f"{outDir}/pinhole_flat_{arm}_order_location_residuals.pdf"
            plt.savefig(filePath)

        self.log.debug('completed the ``calculate_residuals`` method')
        return res_x, res_mean, res_std, res_median

    # use the tab-trigger below for new method
    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
