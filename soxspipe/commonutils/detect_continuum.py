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
from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial
from soxspipe.commonutils.filenamer import filenamer
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
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
import pandas as pd
from soxspipe.commonutils.toolkit import cut_image_slice


class _base_detect(object):

    def fit_polynomial(
            self,
            pixelList,
            order,
            xCol,
            yCol):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``pixelList`` -- data-frame group containing x,y pixel array
            - ``order`` -- the order to fit
            - ``xCol`` -- name of x-pixel column
            - ``yCol`` -- name of y-pixel column

        **Return:**
            - ``coeffs`` -- the coefficients of the polynomial fit
            - ``pixelList`` -- the pixel list but now with fits and residuals included
        """
        self.log.debug('starting the ``fit_polynomial`` method')

        arm = self.arm

        clippedCount = 1

        poly = chebyshev_xy_polynomial(
            log=self.log, deg=self.polyDeg).poly
        clippingSigma = self.recipeSettings[
            "poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.recipeSettings[
            "clipping-iteration-limit"]

        iteration = 0
        mask = (pixelList['order'] == order)
        pixelListFiltered = pixelList.loc[mask]
        while clippedCount > 0 and iteration < clippingIterationLimit:
            pixelListFiltered = pixelList.loc[mask]

            startCount = len(pixelListFiltered.index)
            iteration += 1
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLY
            coeff = np.ones((self.polyDeg + 1))
            # NOTE X AND Y COLUMN ARE CORRECLY IN xdata AND ydata - WANT TO
            # FIND X (UNKNOWN) WRT Y (KNOWNN)
            try:
                coeff, pcov_x = curve_fit(
                    poly, xdata=pixelListFiltered[yCol].values, ydata=pixelListFiltered[xCol].values, p0=coeff)
            except TypeError as e:
                # REMOVE THIS ORDER FROM PIXEL LIST
                pixelList.drop(index=pixelList[mask].index, inplace=True)
                coeff = None
                return coeff, pixelList
            except Exception as e:
                raise e

            res, res_mean, res_std, res_median, xfit = self.calculate_residuals(
                orderPixelTable=pixelListFiltered,
                coeff=coeff,
                xCol=xCol,
                yCol=yCol)

            pixelList.loc[mask, "x_fit_res"] = res
            pixelList.loc[mask, "x_fit"] = xfit

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                res, sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc=mad_std)
            pixelList.loc[mask, "mask"] = masked_residuals.mask

            # REMOVE FILTERED ROWS FROM DATA FRAME
            removeMask = (pixelList["mask"] == True)
            pixelList.drop(index=pixelList[removeMask].index, inplace=True)
            pixelListFiltered = pixelList.loc[mask]
            clippedCount = startCount - len(pixelListFiltered.index)

            # MASK DATA ARRAYS WITH CLIPPED RESIDUAL MASK
            print(f'{clippedCount} pixel positions where clipped in this iteration of fitting an order centre polynomial')

        self.log.debug('completed the ``fit_polynomials`` method')
        return coeff, pixelList

    def calculate_residuals(
            self,
            orderPixelTable,
            coeff,
            xCol,
            yCol):
        """*calculate residuals of the polynomial fits against the observed line postions*

        **Key Arguments:**
            - ``orderPixelTable`` -- data-frame containing pixel list for given order
            - ``coeff`` -- the coefficients of the fitted polynomial
            - ``xCol`` -- name of x-pixel column
            - ``yCol`` -- name of y-pixel column

        **Return:**
            - ``res`` -- x residuals
            - ``mean`` -- the mean of the residuals
            - ``std`` -- the stdev of the residuals
            - ``median`` -- the median of the residuals
            - ``xfit`` -- fitted x values
        """
        self.log.debug('starting the ``calculate_residuals`` method')

        arm = self.arm

        poly = chebyshev_xy_polynomial(
            log=self.log, deg=self.polyDeg).poly

        # CALCULATE RESIDUALS BETWEEN GAUSSIAN PEAK LINE POSITIONS AND POLY
        # FITTED POSITIONS
        xfit = poly(
            orderPixelTable[yCol].values, *coeff)
        res = xfit - orderPixelTable[xCol].values

        # CALCULATE COMBINED RESIDUALS AND STATS
        res_mean = np.mean(res)
        res_std = np.std(res)
        res_median = np.median(res)

        self.log.debug('completed the ``calculate_residuals`` method')
        return res, res_mean, res_std, res_median, xfit

    def write_order_table_to_file(
            self,
            frame,
            orderPolyTable):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**
            - ``frame`` -- the calibration frame used to generate order location data
            - ``orderPolyTable`` -- data-frames containing centre location coefficients (and possibly also order edge coeffs)

        **Return:**
            - ``order_table_path`` -- path to the order table file
        """
        self.log.debug('starting the ``write_order_table_to_file`` method')

        arm = self.arm

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)

        filename = filenamer(
            log=self.log,
            frame=frame,
            settings=self.settings
        )
        filename = filename.replace("MFLAT", "FLAT")
        filename = filename.split("FLAT")[0] + "ORDER_LOCATIONS.csv"

        order_table_path = f"{outDir}/{filename}"
        orderPolyTable.to_csv(order_table_path, index=False)

        self.log.debug('completed the ``write_order_table_to_file`` method')
        return order_table_path


class detect_continuum(_base_detect):
    """
    *find and fit the continuum in a pinhole flat frame with low-order polynomials. These polynominals are the central loctions of the orders*

    **Key Arguments:**
        - ``log`` -- logger
        - ``pinholeFlat`` -- calibrationed pinhole flat frame (CCDObject)
        - ``dispersion_map`` -- path to dispersion map csv file containing polynomial fits of the dispersion solution for the frame
        - ``settings`` -- the recipe settings dictionary
        - ``recipeName`` -- the recipe name as given in the settings dictionary

    **Usage:**

    To use the ``detect_continuum`` object, use the following:

    ```python
    from soxspipe.commonutils import detect_continuum
    detector = detect_continuum(
        log=log,
        pinholeFlat=pinholeFlat,
        dispersion_map=dispersion_map,
        settings=settings,
        recipeName="soxs-order-centre"
    )
    order_table_path = detector.get()
    ```
    """

    def __init__(
            self,
            log,
            pinholeFlat,
            dispersion_map,
            settings=False,
            recipeName=False
    ):
        self.log = log
        log.debug("instansiating a new 'detect_continuum' object")
        self.settings = settings
        if recipeName:
            self.recipeSettings = settings[recipeName]
        else:
            self.recipeSettings = False
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
        self.polyDeg = self.recipeSettings["poly-deg"]

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
        orderPixelTable = self.create_pixel_arrays(
            orderNums,
            waveLengthMin,
            waveLengthMax)
        # SLICE LENGTH TO SAMPLE TRACES IN THE CROSS-DISPERSION DIRECTION
        self.sliceLength = self.recipeSettings["slice-length"]
        self.peakSigmaLimit = self.recipeSettings["peak-sigma-limit"]

        # PREP LISTS WITH NAN VALUE IN CONT_X AND CONT_Y BEFORE FITTING
        orderPixelTable['cont_x'] = np.nan
        orderPixelTable['cont_y'] = np.nan

        # FOR EACH ORDER, FOR EACH PIXEL POSITION SAMPLE, FIT A 1D GAUSSIAN IN
        # CROSS-DISPERSION DIRECTTION. RETURN PEAK POSTIONS
        orderPixelTable = orderPixelTable.apply(
            self.fit_1d_gaussian_to_slice, axis=1)
        allLines = len(orderPixelTable.index)
        # DROP ROWS WITH NAN VALUES
        orderPixelTable.dropna(axis='index', how='any',
                               subset=['cont_x'], inplace=True)
        foundLines = len(orderPixelTable.index)
        percent = 100 * foundLines / allLines
        print(f"{foundLines} out of {allLines} found ({percent:3.0f}%)")

        # GET UNIQUE VALUES IN COLUMN
        uniqueOrders = orderPixelTable['order'].unique()

        orderLocations = {}
        orderPixelTable['x_fit'] = np.nan
        orderPixelTable['x_fit_res'] = np.nan

        for o in uniqueOrders:
            # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
            coeff, orderPixelTable = self.fit_polynomial(
                pixelList=orderPixelTable,
                order=o,
                xCol="cont_x",
                yCol="cont_y"
            )
            orderLocations[o] = coeff

        # SORT CENTRE TRACE COEFFICIENT OUTPUT TO PANDAS DATAFRAME
        columnsNames = ["order", "degy_cent", "ymin", "ymax"]
        coeffColumns = [f'cent_c{i}' for i in range(0, self.polyDeg + 1)]
        columnsNames.extend(coeffColumns)
        myDict = {k: [] for k in columnsNames}
        for k, v in orderLocations.items():
            myDict["order"].append(k)
            myDict["degy_cent"].append(self.polyDeg)
            n_coeff = 0
            for i in range(0, self.polyDeg + 1):
                myDict[f'cent_c{i}'].append(v[n_coeff])
                n_coeff += 1

            myDict["ymin"].append(
                np.min(orderPixelTable.loc[(orderPixelTable['order'] == k)]["fit_y"].values))
            myDict["ymax"].append(
                np.max(orderPixelTable.loc[(orderPixelTable['order'] == k)]["fit_y"].values))
        orderPolyTable = pd.DataFrame(myDict)

        # HERE IS THE LINE LIST IF NEEDED FOR QC
        orderPixelTable.drop(columns=['mask'], inplace=True)

        plotPath = self.plot_results(
            orderPixelTable=orderPixelTable,
            orderPolyTable=orderPolyTable
        )

        mean_res = np.mean(np.abs(orderPixelTable['x_fit_res'].values))
        std_res = np.std(np.abs(orderPixelTable['x_fit_res'].values))

        print(f'\nThe order centre polynomial fitted against the observed 1D gaussian peak positions with a mean residual of {mean_res:2.2f} pixels (stdev = {std_res:2.2f} pixels)')

        # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE
        order_table_path = self.write_order_table_to_file(
            frame=self.pinholeFlat, orderPolyTable=orderPolyTable)

        print(f'\nFind results of the order centre fitting here: {plotPath}')

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

        # READ CSV FILE TO PANDAS DATAFRAME
        specFormatTable = pd.read_csv(
            spectralFormatFile, index_col=False, na_values=['NA', 'MISSING'])
        # # DATAFRAME INFO

        # EXTRACT REQUIRED PARAMETERS
        orderNums = specFormatTable["ORDER"].values
        waveLengthMin = specFormatTable["WLMINFUL"].values
        waveLengthMax = specFormatTable["WLMAXFUL"].values

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
            - ``orderPixelTable`` -- a data-frame containing lines and associated pixel locations
        """
        self.log.debug('starting the ``create_pixel_arrays`` method')

        # READ ORDER SAMPLING RESOLUTION FROM SETTINGS
        sampleCount = self.settings[
            "soxs-order-centre"]["order-sample-count"]

        # CREATE THE WAVELENGTH/ORDER ARRAYS TO BE CONVERTED TO PIXELS
        myDict = {
            "order": np.asarray([]),
            "wavelength": np.asarray([]),
            "slit_position": np.asarray([])
        }
        for o, wmin, wmax in zip(orderNums, waveLengthMin, waveLengthMax):

            wlArray = np.arange(
                wmin, wmax, (wmax - wmin) / sampleCount)
            myDict["wavelength"] = np.append(myDict["wavelength"], wlArray)
            myDict["order"] = np.append(
                myDict["order"], np.ones(len(wlArray)) * o)
            myDict["slit_position"] = np.append(
                myDict["slit_position"], np.zeros(len(wlArray)))

        orderPixelTable = pd.DataFrame(myDict)
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispersion_map,
            orderPixelTable=orderPixelTable
        )

        self.log.debug('completed the ``create_pixel_arrays`` method')
        return orderPixelTable

    def fit_1d_gaussian_to_slice(
            self,
            pixelPostion):
        """*cut a slice from the pinhole flat along the cross-dispersion direction centred on pixel position, fit 1D gaussian and return the peak pixel position*

        **Key Arguments:**
            - ``pixelPostion`` -- the x,y pixel coordinate from orderPixelTable data-frame (series)

        **Return:**
            - ``pixelPostion`` -- now including gaussian fit peak xy position
        """
        self.log.debug('starting the ``fit_1d_gaussian_to_slice`` method')

        # CLIP OUT A SLICE TO INSPECT CENTRED AT POSITION
        halfSlice = self.sliceLength / 2

        slice = cut_image_slice(log=self.log, frame=self.pinholeFlat,
                                width=1, length=self.sliceLength, x=pixelPostion["fit_x"], y=pixelPostion["fit_y"], median=True, plot=False)

        if slice is None:
            pixelPostion["cont_x"] = np.nan
            pixelPostion["cont_y"] = np.nan
            return pixelPostion

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

        if not median_r:
            pixelPostion["cont_x"] = np.nan
            pixelPostion["cont_y"] = np.nan
            return pixelPostion

        peaks, _ = find_peaks(slice, height=median_r +
                              self.peakSigmaLimit * std_r, width=1)

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
            pixelPostion["cont_x"] = np.nan
            pixelPostion["cont_y"] = np.nan
            return pixelPostion

        # FIT THE DATA USING A 1D GAUSSIAN - USING astropy.modeling
        # CENTRE THE GAUSSIAN ON THE PEAK
        g_init = models.Gaussian1D(
            amplitude=1000., mean=peaks[0], stddev=1.)
        # print(f"g_init: {g_init}")
        fit_g = fitting.LevMarLSQFitter()

        # NOW FIT
        g = fit_g(g_init, np.arange(0, len(slice)), slice)
        pixelPostion["cont_x"] = g.mean + \
            max(0, int(pixelPostion["fit_x"] - halfSlice))
        pixelPostion["cont_y"] = pixelPostion["fit_y"]

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
        return pixelPostion

    def plot_results(
            self,
            orderPixelTable,
            orderPolyTable):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**
            - ``orderPixelTable`` -- the pixel table with residuals of fits
            - ``orderPolyTable`` -- data-frame of order-location polynomial coeff

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
        x = np.ones(len(orderPixelTable.index)) * \
            self.pinholeFlat.data.shape[1] - orderPixelTable['cont_x'].values
        toprow.scatter(orderPixelTable[
                       'cont_y'].values, x, marker='x', c='red', s=4)
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

        # ONLY DO THIS FOR SMALL DATAFRAMES - THIS IS AN ANTIPATTERN
        for index, row in orderPolyTable.iterrows():
            o = row["order"]
            coeff = [float(v) for k, v in row.items() if "cent_" in k]
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
        bottomleft.scatter(orderPixelTable['cont_x'].values, orderPixelTable[
                           'x_fit_res'].values, alpha=0.2, s=1)
        bottomleft.set_xlabel('x pixel position')
        bottomleft.set_ylabel('x residual')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomright.scatter(orderPixelTable['cont_y'].values, orderPixelTable[
                            'x_fit_res'].values, alpha=0.2, s=1)
        bottomright.set_xlabel('y pixel position')
        bottomright.tick_params(axis='both', which='major', labelsize=9)
        # bottomright.set_ylabel('x residual')
        bottomright.set_yticklabels([])

        mean_res = np.mean(np.abs(orderPixelTable['x_fit_res'].values))
        std_res = np.std(np.abs(orderPixelTable['x_fit_res'].values))

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        fig.suptitle(f"traces of order-centre locations - pinhole flat-frame\n{subtitle}", fontsize=12)

        # plt.show()
        filename = filenamer(
            log=self.log,
            frame=self.pinholeFlat,
            settings=self.settings
        )
        filename = filename.split("FLAT")[0] + "ORDER_CENTRES_residuals.pdf"

        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        filePath = f"{outDir}/{filename}"
        plt.savefig(filePath)

        self.log.debug('completed the ``plot_results`` method')
        return filePath

    # use the tab-trigger below for new method
    # xt-class-method
