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
from astropy.table import Table
from astropy.io import fits
from soxspipe.commonutils.toolkit import read_spectral_format
from soxspipe.commonutils.toolkit import cut_image_slice
import pandas as pd
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
from fundamentals.renderer import list_of_dictionaries
import collections
from astropy.visualization import hist
from astropy.stats import sigma_clip, mad_std
from scipy.optimize import curve_fit
from random import random
from scipy.signal import find_peaks
from astropy.modeling import models, fitting
from astropy.stats import mad_std
import matplotlib.pyplot as plt
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial, chebyshev_order_xy_polynomials
import numpy as np
from os.path import expanduser
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from fundamentals import tools
from builtins import object
import sys
import os
from io import StringIO
import copy
from contextlib import suppress
os.environ['TERM'] = 'vt100'


class _base_detect(object):

    def fit_order_polynomial(
            self,
            pixelList,
            order,
            y_deg,
            xCol,
            yCol,
            exponents_included=False):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``pixelList`` -- data-frame group containing x,y pixel array
            - ``order`` -- the order to fit
            - ``y_deg`` -- degree for polynomial to fit
            - ``xCol`` -- name of x column
            - ``yCol`` -- name of y column
            - ``exponents_included`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

        **Return:**
            - ``coeffs`` -- the coefficients of the polynomial fit
            - ``pixelList`` -- the pixel list but now with fits and residuals included
        """
        self.log.debug('starting the ``fit_order_polynomial`` method')

        arm = self.arm
        self.y_deg = y_deg

        clippedCount = 1

        poly = chebyshev_xy_polynomial(
            log=self.log, yCol=yCol, y_deg=y_deg, exponents_included=exponents_included).poly

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
            coeff = np.ones((y_deg + 1))
            # NOTE X AND Y COLUMN ARE CORRECLY IN xdata AND ydata - WANT TO
            # FIND X (UNKNOWN) WRT Y (KNOWNN)
            try:
                coeff, pcov_x = curve_fit(
                    poly, xdata=pixelListFiltered, ydata=pixelListFiltered[xCol].values, p0=coeff)
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
                y_deg=y_deg,
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

            sys.stdout.write("\x1b[1A\x1b[2K")
            print(f'\t\tORDER {order:0.0f}: {clippedCount} pixel positions where clipped in iteration {iteration} of fitting the polynomial')

        self.log.debug('completed the ``fit_order_polynomial`` method')
        return coeff, pixelList

    def fit_global_polynomial(
            self,
            pixelList,
            y_deg,
            order_deg,
            xCol="cont_x",
            yCol="cont_y",
            orderCol="order",
            exponents_included=False):
        """*iteratively fit the global polynomial to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``pixelList`` -- data-frame group containing x,y pixel array
            - ``y_deg`` -- degree for polynomial to fit y-values
            - ``order_deg`` -- degree for polynomial to fit order-values
            - ``exponents_included`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

        **Return:**
            - ``coeffs`` -- the coefficients of the polynomial fit
            - ``pixelList`` -- the pixel list but now with fits and residuals included
        """
        self.log.debug('starting the ``fit_global_polynomial`` method')

        arm = self.arm

        clippedCount = 1

        poly = chebyshev_order_xy_polynomials(log=self.log, yCol=yCol, orderCol=orderCol, order_deg=order_deg, y_deg=y_deg, exponents_included=exponents_included).poly

        clippingSigma = self.recipeSettings[
            "poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.recipeSettings[
            "clipping-iteration-limit"]

        iteration = 0

        while clippedCount > 0 and iteration < clippingIterationLimit:
            startCount = len(pixelList.index)
            iteration += 1
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLY
            coeff = np.ones((self.yDeg + 1) * (self.orderDeg + 1))
            try:
                coeff, pcov_x = curve_fit(
                    poly, xdata=pixelList, ydata=pixelList[xCol].values, p0=coeff)
            except TypeError as e:
                # REMOVE THIS ORDER FROM PIXEL LIST
                coeff = None
                return coeff, pixelList
            except Exception as e:
                raise e

            res, res_mean, res_std, res_median, xfit = self.calculate_residuals(
                orderPixelTable=pixelList,
                coeff=coeff,
                y_deg=y_deg,
                order_deg=order_deg,
                orderCol=orderCol,
                xCol=xCol,
                yCol=yCol)

            pixelList["x_fit_res"] = res
            pixelList["x_fit"] = xfit

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                res, sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc=mad_std)
            pixelList["mask"] = masked_residuals.mask

            # REMOVE FILTERED ROWS FROM DATA FRAME
            removeMask = (pixelList["mask"] == True)
            pixelList.drop(index=pixelList[removeMask].index, inplace=True)
            clippedCount = startCount - len(pixelList.index)

            sys.stdout.write("\x1b[1A\x1b[2K")
            print(f'\t\tGLOBAL FIT: {clippedCount} pixel positions where clipped in iteration {iteration} of fitting the polynomial')

        self.log.debug('completed the ``fit_global_polynomial`` method')
        return coeff, pixelList

    def calculate_residuals(
            self,
            orderPixelTable,
            coeff,
            y_deg,
            xCol,
            yCol,
            orderCol=False,
            order_deg=False):
        """*calculate residuals of the polynomial fits against the observed line postions*

        **Key Arguments:**
            - ``orderPixelTable`` -- data-frame containing pixel list for given order
            - ``coeff`` -- the coefficients of the fitted polynomial
            - ``y_deg`` -- degree for polynomial to fit y-values
            - ``xCol`` -- name of x-pixel column
            - ``yCol`` -- name of y-pixel column
            - ``orderCol`` -- name of the order column (global fits only)
            - ``order_deg`` -- degree for polynomial to fit order-values (global fits only)

        **Return:**
            - ``res`` -- x residuals
            - ``mean`` -- the mean of the residuals
            - ``std`` -- the stdev of the residuals
            - ``median`` -- the median of the residuals
            - ``xfit`` -- fitted x values
        """
        self.log.debug('starting the ``calculate_residuals`` method')

        arm = self.arm

        poly = chebyshev_order_xy_polynomials(
            log=self.log, yCol=yCol, orderCol=orderCol, order_deg=self.orderDeg, y_deg=self.yDeg).poly

        # CALCULATE RESIDUALS BETWEEN GAUSSIAN PEAK LINE POSITIONS AND POLY
        # FITTED POSITIONS
        xfit = poly(
            orderPixelTable, *coeff)
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
            orderPolyTable,
            orderMetaTable):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**
            - ``frame`` -- the calibration frame used to generate order location data
            - ``orderPolyTable`` -- data-frames containing centre location coefficients (and possibly also order edge coeffs)
            - ``orderMetaTable`` -- extra order meta data to be added in an extra FITS extension

        **Return:**
            - ``order_table_path`` -- path to the order table file
        """
        from astropy.table import Table
        self.log.debug('starting the ``write_order_table_to_file`` method')

        arm = self.arm
        kw = self.kw

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)

        filename = filenamer(
            log=self.log,
            frame=frame,
            settings=self.settings
        )
        filename = filename.replace("MFLAT", "FLAT")
        filename = filename.split("FLAT")[0] + "ORDER_LOCATIONS.fits"

        order_table_path = f"{outDir}/{filename}"

        header = copy.deepcopy(frame.header)
        header.pop(kw("DPR_TECH"))
        header.pop(kw("DPR_CATG"))
        header.pop(kw("DPR_TYPE"))

        with suppress(KeyError):
            header.pop(kw("DET_READ_SPEED"))
        with suppress(KeyError):
            header.pop(kw("CONAD"))
        with suppress(KeyError):
            header.pop(kw("GAIN"))
        with suppress(KeyError):
            header.pop(kw("RON"))

        header[kw("PRO_TECH")] = "ECHELLE,SLIT"

        orderPolyTable = Table.from_pandas(orderPolyTable)
        BinTableHDU = fits.table_to_hdu(orderPolyTable)
        orderMetaTable = Table.from_pandas(orderMetaTable)
        BinTableHDU2 = fits.table_to_hdu(orderMetaTable)

        header[kw("SEQ_ARM")] = arm
        header[kw("PRO_TYPE")] = "REDUCED"
        header[kw("PRO_CATG")] = f"ORDER_TAB_{arm}".upper()
        priHDU = fits.PrimaryHDU(header=header)

        hduList = fits.HDUList([priHDU, BinTableHDU, BinTableHDU2])
        hduList.writeto(order_table_path, checksum=True, overwrite=True)

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
        - ``qcTable`` -- the data frame to collect measured QC metrics 
        - ``productsTable`` -- the data frame to collect output products

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
            recipeName=False,
            qcTable=False,
            productsTable=False
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
        self.qc = qcTable
        self.products = productsTable

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        self.arm = pinholeFlat.header[self.kw("SEQ_ARM")]
        self.dateObs = pinholeFlat.header[self.kw("DATE_OBS")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        # DEG OF THE POLYNOMIALS TO FIT THE ORDER CENTRE LOCATIONS
        self.yDeg = self.recipeSettings["y-deg"]
        self.orderDeg = self.recipeSettings["order-deg"]

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
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=arm)

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

        # SETUP EXPONENTS AHEAD OF TIME - SAVES TIME ON POLY FITTING
        for i in range(0, self.yDeg + 1):
            orderPixelTable[f"y_pow_{i}"] = orderPixelTable["cont_y"].pow(i)
        for i in range(0, self.orderDeg + 1):
            orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)
        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        coeff, orderPixelTable = self.fit_global_polynomial(
            pixelList=orderPixelTable,
            y_deg=self.yDeg,
            order_deg=self.orderDeg,
            exponents_included=True
        )

        # orderLocations[o] = coeff
        coeff_dict = {"degorder_cent": self.orderDeg,
                      "degy_cent": self.yDeg}
        n_coeff = 0
        for i in range(0, self.orderDeg + 1):
            for j in range(0, self.yDeg + 1):
                coeff_dict[f'cent_{i}{j}'] = coeff[n_coeff]
                n_coeff += 1
        coeffColumns = coeff_dict.keys()
        dataSet = list_of_dictionaries(
            log=self.log,
            listOfDictionaries=[coeff_dict]
        )

        # WRITE CSV DATA TO PANDAS DATAFRAME TO ASTROPY TABLE TO FITS
        fakeFile = StringIO(dataSet.csv())
        orderPolyTable = pd.read_csv(fakeFile, index_col=False, na_values=['NA', 'MISSING'])
        fakeFile.close()

        # HERE IS THE LINE LIST IF NEEDED FOR QC
        orderPixelTable.drop(columns=['mask'], inplace=True)

        plotPath, orderMetaTable = self.plot_results(
            orderPixelTable=orderPixelTable,
            orderPolyTable=orderPolyTable
        )

        mean_res = np.mean(np.abs(orderPixelTable['x_fit_res'].values))
        std_res = np.std(np.abs(orderPixelTable['x_fit_res'].values))

        print(f'\nThe order centre polynomial fitted against the observed 1D gaussian peak positions with a mean residual of {mean_res:2.2f} pixels (stdev = {std_res:2.2f} pixels)')

        # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE
        order_table_path = self.write_order_table_to_file(
            frame=self.pinholeFlat, orderPolyTable=orderPolyTable, orderMetaTable=orderMetaTable)

        print(f'\nFind results of the order centre fitting here: {plotPath}')

        self.log.debug('completed the ``get`` method')
        return order_table_path

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
            - ``orderMetaTable`` -- dataframe of useful order fit metadata
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

        poly = chebyshev_order_xy_polynomials(
            log=self.log, yCol="y", orderCol="order", order_deg=self.orderDeg, y_deg=self.yDeg).poly
        for index, row in orderPolyTable.iterrows():
            coeff = [float(v) for k, v in row.items() if "cent_" in k]
        uniqueOrders = orderPixelTable['order'].unique()
        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {"y": ylinelist}
        df = pd.DataFrame(myDict)
        ymin = []
        ymax = []
        xmin = []
        xmax = []
        for o in uniqueOrders:
            df["order"] = o
            xfit = poly(df, *coeff)
            xfit = np.ones(len(xfit)) * \
                self.pinholeFlat.data.shape[1] - xfit
            xfit, yfit = zip(
                *[(x, y) for x, y in zip(xfit, ylinelist) if x > 0 and x < (self.pinholeFlat.data.shape[1]) - 10])
            l = midrow.plot(yfit, xfit)
            ymin.append(min(yfit))
            ymax.append(max(yfit))
            xmin.append(self.pinholeFlat.data.shape[1] - max(xfit))
            xmax.append(self.pinholeFlat.data.shape[1] - min(xfit))
            midrow.text(yfit[10], xfit[10] - 20, int(o), fontsize=6, c="white", verticalalignment='bottom')

        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        orderMetaTable = {
            "order": uniqueOrders,
            "ymin": ymin,
            "ymax": ymax,
            "xmin": xmin,
            "xmax": xmax,
        }
        orderMetaTable = pd.DataFrame(orderMetaTable)

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
        return filePath, orderMetaTable

    # use the tab-trigger below for new method
    # xt-class-method
