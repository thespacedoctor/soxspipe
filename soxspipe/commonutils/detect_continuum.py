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
import matplotlib.pyplot as plt
from astropy.stats import mad_std
from astropy.modeling import models, fitting
from scipy.signal import find_peaks
from random import random


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
                    "axis", "order_deg", "wavelength_deg"]]
                poly = chebyshev_order_wavelength_polynomials(
                    log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg).poly
                if axis == "x":
                    xcoords = poly(order_wave, *coeff)
                if axis == "y":
                    ycoords = poly(order_wave, *coeff)
        csvFile.close()

        for o, x, y in zip(orderArray, xcoords, ycoords):
            pixelArrays[f"o{o:0.0f}"].append((x, y))

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

    # use the tab-trigger below for new method
    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
