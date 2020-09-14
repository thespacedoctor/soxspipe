#!/usr/bin/env python
# encoding: utf-8
"""
*find and fit the continuum in a pinhole flat frame with low-order polynomials. These polynominals are the central loctions of the orders*

:Author:
    David Young & Marco Landoni

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
import collections
from fundamentals.renderer import list_of_dictionaries


class detect_continuum(object):
    """
    *find and fit the continuum in a pinhole flat frame with low-order polynomials. These polynominals are the central loctions of the orders*

    **Key Arguments:**
        - ``log`` -- logger
        - ``pinholeFlat`` -- calibrationed pinhole flat frame (CCDObject)
        - ``dispersion_map`` -- path to dispersion map csv file containing polynomial fits of the dispersion solution for the frame
        - ``settings`` -- the settings dictionary

    **Usage:**

    To use the ``detect_continuum`` object, use the following:

    ```python
    from soxspipe.commonutils import detect_continuum
    detector = detect_continuum(
        log=log,
        pinholeFlat=pinholeFlat,
        dispersion_map=dispersion_map,
        settings=settings
    )
    order_table_path = detector.get()
    ```

    """

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

        # DEG OF THE POLYNOMIALS TO FIT THE ORDER CENTRE LOCATIONS
        self.polyDeg = self.settings[
            "soxs-order-centre"]["poly-deg"]

        return None

    def get(self):
        """
        *return the order centre table filepath*

        **Return:**
            - ``order_table_path`` -- file path to the order centre table giving polynomial coeffs to each order fit
        """
        self.log.debug('starting the ``get`` method')

        arm = self.arm

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = self.read_spectral_format()

        # CONVERT WAVELENGTH TO PIXEL POSTIONS AND RETURN ARRAY OF POSITIONS TO
        # SAMPLE THE TRACES
        pixelArrays = self.create_pixel_arrays(
            orderNums,
            waveLengthMin,
            waveLengthMax)

        # SLICE LENGTH TO SAMPLE TRACES IN THE CROSS-DISPERSION DIRECTION
        sliceLength = self.settings[
            "soxs-order-centre"]["slice-length"]
        peakSigmaLimit = self.settings[
            "soxs-order-centre"]["peak-sigma-limit"]

        # FOR EACH ORDER, FOR EACH PIXEL POSITION SAMPLE, FIT A 1D GAUSSIAN IN
        # CROSS-DISPERSION DIRECTTION. RETURN PEAK POSTIONS
        guassianPixelArrays = {}
        for order, xy_array in pixelArrays.items():
            found_xy_array = []
            for xy in xy_array:
                found_xy = self.fit_1d_gaussian_to_slice(
                    xy, sliceLength, peakSigmaLimit)
                if found_xy:
                    found_xy_array.append(found_xy)
            guassianPixelArrays[order] = found_xy_array

            percent = 100 * len(found_xy_array) / len(xy_array)
            print(f"{len(found_xy_array)} out of {len(xy_array)} found ({percent:3.0f}%)")

        orderLoctions = {}
        allResiduals = []
        allXfit = []
        allXcoords = []
        allYcoords = []
        for order, pixelArray in guassianPixelArrays.items():
            # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
            coeff, residuals, xfit, xcoords, ycoords = self.fit_polynomial(
                pixelArray=pixelArray
            )
            allResiduals.extend(residuals)
            allXfit.extend(xfit)
            allXcoords.extend(xcoords)
            allYcoords.extend(ycoords)
            orderLoctions[order] = coeff

        self.plot_results(
            allResiduals=allResiduals,
            allXfit=allXfit,
            allXcoords=allXcoords,
            allYcoords=allYcoords,
            orderLoctions=orderLoctions
        )

        mean_res = np.mean(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        print(f'\nThe order centre polynomial fitted against the observed 1D gaussian peak positions with a mean residual of {mean_res:2.2f} pixels (stdev = {std_res:2.2f} pixles)')

        # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE
        order_table_path = self.write_order_table_to_file(orderLoctions)

        self.log.debug('completed the ``get`` method')
        return order_table_path

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
        waveLengthMin = spectralFormat["WLMINFUL"]
        waveLengthMax = spectralFormat["WLMAXFUL"]

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

        xcoords, ycoords, orderArray = zip(
            *[(x, y, o) for x, y, o in zip(xcoords, ycoords, orderArray) if (x > 0 and y > 0)])

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

        try:
            slice = self.pinholeFlat.data[int(y_fit), max(
                0, int(x_fit - halfSlice)):min(2048, int(x_fit + halfSlice))]
        except:
            return None

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
            pixelArray):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``pixelArray`` -- the array of x,y pixels of measured 1D guassian peak positions

        **Return:**
            - ``coeffs`` -- the coefficients of the polynomial fit
            - ``residuals`` -- the residuals of the fit compared to original guassian peak positions
            - ``xfit`` -- the polynomial fits x-positions
            - ``xcoords`` -- the clean guassian peak x-coordinate list (post-clipping)
            - ``ycoords`` -- the clean guassian peak x-coordinate list (post-clipping)
        """
        self.log.debug('starting the ``fit_polynomial`` method')

        arm = self.arm

        clippedCount = 1

        poly = chebyshev_xy_polynomial(
            log=self.log, deg=self.polyDeg).poly

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
            coeff = np.ones((self.polyDeg + 1))
            coeff, pcov_x = curve_fit(
                poly, xdata=ycoords, ydata=xcoords, p0=coeff)

            residuals, mean_res, std_res, median_res, xfit = self.calculate_residuals(
                xcoords=xcoords,
                ycoords=ycoords,
                coeff=coeff)

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                residuals, sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc=mad_std)

            # MASK DATA ARRAYS WITH CLIPPED RESIDUAL MASK
            startCount = len(xcoords)
            a = [xcoords, ycoords, residuals]
            xcoords, ycoords, residuals = [np.ma.compressed(np.ma.masked_array(
                i, masked_residuals.mask)) for i in a]
            clippedCount = startCount - len(xcoords)
            print(f'{clippedCount} pixel positions where clipped in this iteration of fitting an order centre polynomial')

        self.log.debug('completed the ``fit_polynomials`` method')
        return coeff, residuals, xfit, xcoords, ycoords

    def calculate_residuals(
            self,
            xcoords,
            ycoords,
            coeff):
        """*calculate residuals of the polynomial fits against the observed line postions*

        **Key Arguments:**
            - ``xcoords`` -- the measured x positions of the gaussian peaks
            - ``ycoords`` -- the measurd y positions of the gaussian peaks
            - ``coeff`` -- the coefficients of the fitted polynomial

        **Return:**
            - ``residuals`` -- x residuals
            - ``mean`` -- the mean of the residuals
            - ``std`` -- the stdev of the residuals
            - ``median`` -- the median of the residuals
        """
        self.log.debug('starting the ``calculate_residuals`` method')

        arm = self.arm

        poly = chebyshev_xy_polynomial(
            log=self.log, deg=self.polyDeg).poly

        # CALCULATE RESIDUALS BETWEEN GAUSSIAN PEAK LINE POSITIONS AND POLY
        # FITTED POSITIONS
        xfit = poly(
            ycoords, *coeff)
        res_x = np.asarray(xfit) - np.asarray(xcoords)

        # CALCULATE COMBINED RESIDUALS AND STATS
        res_mean = np.mean(res_x)
        res_std = np.std(res_x)
        res_median = np.median(res_x)

        self.log.debug('completed the ``calculate_residuals`` method')
        return res_x, res_mean, res_std, res_median, xfit

    def plot_results(
            self,
            allResiduals,
            allXfit,
            allXcoords,
            allYcoords,
            orderLoctions):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**
            - ``allResiduals`` -- list of all residuals
            - ``allXfit`` -- list of all fitted x-positions 
            - ``allXcoords`` -- cleaned list of all guassian x-pixel positions
            - ``allYcoords`` -- cleaned list of all guassian y-pixel positions
            - ``orderLoctions`` -- dictionary of order-location polynomial coeff

        **Return:**
            - ``filePath`` -- path to the plot pdf
        """
        self.log.debug('starting the ``plot_results`` method')

        arm = self.arm

        # a = plt.figure(figsize=(40, 15))
        if arm == "UVB":
            fig = plt.figure(figsize=(6, 13.5), constrained_layout=True)
        else:
            fig = plt.figure(figsize=(6, 11), constrained_layout=True)
        gs = fig.add_gridspec(6, 4)

        # CREATE THE GID OF AXES
        toprow = fig.add_subplot(gs[0:2, :])
        midrow = fig.add_subplot(gs[2:4, :])
        bottomleft = fig.add_subplot(gs[4:, 0:2])
        bottomright = fig.add_subplot(gs[4:, 2:])

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.rot90(self.pinholeFlat.data, 1)
        toprow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        toprow.set_title(
            "1D guassian peak positions (post-clipping)", fontsize=10)
        x = np.ones(len(allXcoords)) * \
            self.pinholeFlat.data.shape[1] - allXcoords
        toprow.scatter(allYcoords, x, marker='x', c='red', s=4)
        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])
        toprow.set_ylabel("x-axis", fontsize=8)
        toprow.set_xlabel("y-axis", fontsize=8)
        toprow.tick_params(axis='both', which='major', labelsize=9)

        midrow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        midrow.set_title(
            "order-location fit solutions", fontsize=10)
        ylinelist = np.arange(0, self.pinholeFlat.data.shape[0], 3)
        poly = chebyshev_xy_polynomial(
            log=self.log, deg=self.polyDeg).poly
        for o, coeff in orderLoctions.items():
            xfit = poly(ylinelist, *coeff)
            xfit = np.ones(len(xfit)) * \
                self.pinholeFlat.data.shape[1] - xfit
            xfit, ylinelist = zip(
                *[(x, y) for x, y in zip(xfit, ylinelist) if x > 0 and x < (self.pinholeFlat.data.shape[1]) - 10])
            midrow.plot(ylinelist, xfit)

        # xfit = np.ones(len(xfit)) * \
        #     self.pinholeFrame.data.shape[1] - xfit
        # midrow.scatter(yfit, xfit, marker='x', c='blue', s=4)
        # midrow.set_yticklabels([])
        # midrow.set_xticklabels([])
        midrow.set_ylabel("x-axis", fontsize=8)
        midrow.set_xlabel("y-axis", fontsize=8)
        midrow.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomleft.scatter(allXcoords, allResiduals, alpha=0.6, s=2)
        bottomleft.set_xlabel('x pixel position')
        bottomleft.set_ylabel('x residual')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomright.scatter(allYcoords, allResiduals, alpha=0.6, s=2)
        bottomright.set_xlabel('y pixel position')
        bottomright.tick_params(axis='both', which='major', labelsize=9)
        # bottomright.set_ylabel('x residual')
        bottomright.set_yticklabels([])

        mean_res = np.mean(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        fig.suptitle(f"traces of order-centre locations - pinhole flat-frame\n{subtitle}", fontsize=12)

        # plt.show()
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        filePath = f"{outDir}/pinhole_flat_{arm}_order_location_residuals.pdf"
        plt.savefig(filePath)

        self.log.debug('completed the ``plot_results`` method')
        return filePath

    def write_order_table_to_file(
            self,
            orderLoctions):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**
            - ``orderLoctions`` -- dictionary of the order coefficients

        **Return:**
            - ``order_table_path`` -- path to the order table file
        """
        self.log.debug('starting the ``write_order_table_to_file`` method')

        arm = self.arm

        # SORT COEFFICIENT OUTPUT TO WRITE TO FILE
        listOfDictionaries = []
        for k, v in orderLoctions.items():
            orderDict = collections.OrderedDict(sorted({}.items()))
            orderDict["order"] = int(k.replace("o", ""))
            orderDict["degy"] = self.polyDeg
            n_coeff = 0
            for i in range(0, self.polyDeg + 1):
                orderDict[f'CENT_c{i}'] = v[n_coeff]
                n_coeff += 1
            listOfDictionaries.append(orderDict)

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        order_table_path = f"{outDir}/order_centre_{arm}_locations.csv"
        dataSet = list_of_dictionaries(
            log=self.log,
            listOfDictionaries=listOfDictionaries
        )
        csvData = dataSet.csv(filepath=order_table_path)

        self.log.debug('completed the ``write_order_table_to_file`` method')
        return order_table_path

    # use the tab-trigger below for new method
    # xt-class-method
