#!/usr/bin/env python
# encoding: utf-8
"""
*detect arc-lines on a pinhole frame to generate a dispersion solution*

:Author:
    Marco Landoni & David Young

:Date Created:
    September  1, 2020
"""
################# GLOBAL IMPORTS ####################
from fundamentals import fmultiprocess
from tabulate import tabulate
import pandas as pd
from ccdproc import Combiner
from soxspipe.commonutils.toolkit import unpack_order_table, read_spectral_format
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
from astropy.visualization import hist
import warnings
from photutils.utils import NoDetectionsWarning
from astropy.nddata import CCDData
import math
from soxspipe.commonutils.toolkit import get_calibrations_path
import numpy as np
from astropy.stats import sigma_clip, mad_std
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib import cm
from os.path import expanduser
from fundamentals.renderer import list_of_dictionaries
from scipy.optimize import curve_fit
from photutils import DAOStarFinder
from photutils import datasets
from astropy.stats import sigma_clipped_stats
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from fundamentals import tools
from builtins import object
import sys
from astropy.table import Table
import os
from io import StringIO
from contextlib import suppress
from astropy.io import fits
import copy

os.environ['TERM'] = 'vt100'


class create_dispersion_map(object):
    """
    *detect arc-lines on a pinhole frame to generate a dispersion solution*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``pinholeFrame`` -- the calibrated pinhole frame (single or multi)
        - ``firstGuessMap`` -- the first guess dispersion map from the `soxs_disp_solution` recipe (needed in `soxs_spat_solution` recipe). Default *False*.
        - ``orderTable`` -- the order geometry table
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products

    **Usage:**

    ```python
    from soxspipe.commonutils import create_dispersion_map
    mapPath, mapImagePath = create_dispersion_map(
        log=log,
        settings=settings,
        pinholeFrame=frame,
        firstGuessMap=False
    ).get()
    ```
    """

    def __init__(
            self,
            log,
            settings,
            pinholeFrame,
            firstGuessMap=False,
            orderTable=False,
            qcTable=False,
            productsTable=False
    ):
        self.log = log
        log.debug("instantiating a new 'create_dispersion_map' object")
        self.settings = settings
        self.pinholeFrame = pinholeFrame
        self.firstGuessMap = firstGuessMap
        self.orderTable = orderTable
        self.qc = qcTable
        self.products = productsTable

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        self.kw = kw
        self.arm = pinholeFrame.header[kw("SEQ_ARM")]
        self.dateObs = pinholeFrame.header[kw("DATE_OBS")]

        # WHICH RECIPE ARE WE WORKING WITH?
        if self.firstGuessMap:
            self.recipeSettings = self.settings["soxs-spatial-solution"]
        else:
            self.recipeSettings = self.settings["soxs-disp-solution"]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        warnings.simplefilter('ignore', NoDetectionsWarning)

        return None

    def get(self):
        """
        *generate the dispersion map*

        **Return:**
            - ``mapPath`` -- path to the file containing the coefficients of the x,y polynomials of the global dispersion map fit
        """
        self.log.debug('starting the ``get`` method')

        # WHICH RECIPE ARE WE WORKING WITH?
        if self.firstGuessMap:
            slit_deg = self.recipeSettings["slit-deg"]
        else:
            slit_deg = 0

        # READ PREDICTED LINE POSITIONS FROM FILE - RETURNED AS DATAFRAME
        orderPixelTable = self.get_predicted_line_list()

        # GET THE WINDOW SIZE FOR ATTEMPTING TO DETECT LINES ON FRAME
        windowSize = self.recipeSettings["pixel-window-size"]
        self.windowHalf = int(windowSize / 2)

        # DETECT THE LINES ON THE PINHILE FRAME AND
        # ADD OBSERVED LINES TO DATAFRAME
        orderPixelTable = orderPixelTable.apply(
            self.detect_pinhole_arc_line, axis=1)

        # DROP MISSING VALUES
        orderPixelTable.dropna(axis='index', how='any', subset=[
            'observed_x'], inplace=True)

        order_deg = self.recipeSettings["order-deg"]
        wavelength_deg = self.recipeSettings["wavelength-deg"]

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        popt_x, popt_y = self.fit_polynomials(
            orderPixelTable=orderPixelTable,
            wavelength_deg=wavelength_deg,
            order_deg=order_deg,
            slit_deg=slit_deg,
        )

        # WRITE THE MAP TO FILE
        mapPath = self.write_map_to_file(
            popt_x, popt_y, order_deg, wavelength_deg, slit_deg)

        if self.firstGuessMap and self.orderTable:
            mapImagePath = self.map_to_image(dispersionMapPath=mapPath)
            return mapPath, mapImagePath

        self.log.debug('completed the ``get`` method')
        return mapPath, None

    def get_predicted_line_list(
            self):
        """*lift the predicted line list from the static calibrations*

        **Return:**
            - ``orderPixelTable`` -- a panda's data-frame containing wavelength,order,slit_index,slit_position,detector_x,detector_y
        """
        self.log.debug('starting the ``get_predicted_line_list`` method')

        kw = self.kw
        pinholeFrame = self.pinholeFrame
        dp = self.detectorParams

        # WHICH TYPE OF PINHOLE FRAME DO WE HAVE - SINGLE OR MULTI
        if self.pinholeFrame.header[kw("DPR_TECH")] == "ECHELLE,PINHOLE":
            frameTech = "single"
        elif self.pinholeFrame.header[kw("DPR_TECH")] == "ECHELLE,MULTI-PINHOLE":
            frameTech = "multi"
        else:
            raise TypeError(
                "The input frame needs to be a calibrated single- or multi-pinhole arc lamp frame")

        # FIND THE APPROPRIATE PREDICTED LINE-LIST
        arm = self.arm
        if arm != "NIR" and kw('WIN_BINX') in pinholeFrame.header:
            binx = int(self.pinholeFrame.header[kw('WIN_BINX')])
            biny = int(self.pinholeFrame.header[kw('WIN_BINY')])
        else:
            binx = 1
            biny = 1

        # READ THE FILE
        home = expanduser("~")

        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)
        predictedLinesFile = calibrationRootPath + "/" + dp["predicted pinhole lines"][frameTech][f"{binx}x{biny}"]

        # LINE LIST TO PANDAS DATAFRAME
        print(predictedLinesFile)
        dat = Table.read(predictedLinesFile, format='fits')
        orderPixelTable = dat.to_pandas()

        # RENAME ALL COLUMNS FOR CONSISTENCY
        listName = []
        listName[:] = [l if l else l for l in listName]
        orderPixelTable.columns = [d.lower() if d.lower() in [
            "order", "wavelength"] else d for d in orderPixelTable.columns]

        # WANT TO DETERMINE SYSTEMATIC SHIFT IF FIRST GUESS SOLUTION PRESENT
        if self.firstGuessMap:
            # ADD SOME EXTRA COLUMNS TO DATAFRAME

            # FILTER THE PREDICTED LINES TO ONLY SLIT POSITION INCLUDED IN
            # SINGLE PINHOLE FRAMES
            slitIndex = int(dp["mid_slit_index"])

            # GET THE OBSERVED PIXELS VALUES
            orderPixelTable = dispersion_map_to_pixel_arrays(
                log=self.log,
                dispersionMapPath=self.firstGuessMap,
                orderPixelTable=orderPixelTable
            )

            # CREATE A COPY OF THE DATA-FRAME TO DETERMINE SHIFTS
            tmpList = orderPixelTable.copy()

            mask = (tmpList['slit_index'] == slitIndex)
            tmpList.loc[mask, 'shift_x'] = tmpList.loc[
                mask, 'detector_x'].values - tmpList.loc[mask, 'fit_x'].values
            tmpList.loc[mask, 'shift_y'] = tmpList.loc[
                mask, 'detector_y'].values - tmpList.loc[mask, 'fit_y'].values

            # MERGING SHIFTS INTO MAIN DATAFRAME
            tmpList = tmpList.loc[tmpList['shift_x'].notnull(
            ), ['wavelength', 'order', 'shift_x', 'shift_y']]
            orderPixelTable = orderPixelTable.merge(tmpList, on=[
                'wavelength', 'order'], how='outer')

            # DROP ROWS WITH MISSING SHIFTS
            orderPixelTable.dropna(axis='index', how='any', subset=[
                'shift_x'], inplace=True)

            # SHIFT DETECTOR LINE PIXEL POSITIONS BY SHIFTS
            # UPDATE FILTERED VALUES
            orderPixelTable.loc[
                :, 'detector_x'] -= orderPixelTable.loc[:, 'shift_x']
            orderPixelTable.loc[
                :, 'detector_y'] -= orderPixelTable.loc[:, 'shift_y']

            # DROP HELPER COLUMNS
            orderPixelTable.drop(columns=['fit_x', 'fit_y',
                                          'shift_x', 'shift_y'], inplace=True)

        self.log.debug('completed the ``get_predicted_line_list`` method')
        return orderPixelTable

    def detect_pinhole_arc_line(
            self,
            predictedLine):
        """*detect the observed position of an arc-line given the predicted pixel positions*

        **Key Arguments:**
            - ``predictedLine`` -- single predicted line coordinates from predicted line-list

        **Return:**
            - ``predictedLine`` -- the line with the observed pixel coordinates appended (if detected, otherwise nan)
        """
        self.log.debug('starting the ``detect_pinhole_arc_line`` method')

        pinholeFrame = self.pinholeFrame
        windowHalf = self.windowHalf
        x = predictedLine['detector_x']
        y = predictedLine['detector_y']

        # CLIP A STAMP FROM IMAGE AROUNDS PREDICTED POSITION
        xlow = int(np.max([x - windowHalf, 0]))
        xup = int(np.min([x + windowHalf, pinholeFrame.shape[1]]))
        ylow = int(np.max([y - windowHalf, 0]))
        yup = int(np.min([y + windowHalf, pinholeFrame.shape[0]]))
        stamp = pinholeFrame[ylow:yup, xlow:xup]
        # CONVERT TO MASKED ARRAY
        stamp = np.asanyarray(stamp)

        # USE DAOStarFinder TO FIND LINES WITH 2D GUASSIAN FITTING
        mean, median, std = sigma_clipped_stats(stamp, sigma=3.0)

        daofind = DAOStarFinder(
            fwhm=2.0, threshold=5. * std, roundlo=-3.0, roundhi=3.0, sharplo=-3.0, sharphi=3.0)
        sources = daofind(stamp - median)

        # plt.clf()
        # plt.imshow(stamp)
        old_resid = windowHalf * 4
        if sources:
            # FIND SOURCE CLOSEST TO CENTRE
            if len(sources) > 1:
                for source in sources:
                    tmp_x = source['xcentroid']
                    tmp_y = source['ycentroid']
                    new_resid = ((windowHalf - tmp_x)**2 +
                                 (windowHalf - tmp_y)**2)**0.5
                    if new_resid < old_resid:
                        observed_x = tmp_x + xlow
                        observed_y = tmp_y + ylow
                        old_resid = new_resid
            else:
                observed_x = sources[0]['xcentroid'] + xlow
                observed_y = sources[0]['ycentroid'] + ylow
            # plt.scatter(observed_x - xlow, observed_y -
            #             ylow, marker='x', s=30)
            # plt.show()
        else:
            observed_x = np.nan
            observed_y = np.nan
        # plt.show()

        predictedLine['observed_x'] = observed_x
        predictedLine['observed_y'] = observed_y

        self.log.debug('completed the ``detect_pinhole_arc_line`` method')
        return predictedLine

    def write_map_to_file(
            self,
            xcoeff,
            ycoeff,
            order_deg,
            wavelength_deg,
            slit_deg):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**
            - ``xcoeff`` -- the x-coefficients
            - ``ycoeff`` -- the y-coefficients
            - ``order_deg`` -- degree of the order fitting
            - ``wavelength_deg`` -- degree of wavelength fitting
            - ``slit_deg`` -- degree of the slit fitting (False for single pinhole)

        **Return:**
            - ``disp_map_path`` -- path to the saved file
        """
        self.log.debug('starting the ``write_map_to_file`` method')

        arm = self.arm
        kw = self.kw

        # SORT X COEFFICIENT OUTPUT TO WRITE TO FILE
        coeff_dict_x = {}
        coeff_dict_x["axis"] = "x"
        coeff_dict_x["order-deg"] = order_deg
        coeff_dict_x["wavelength-deg"] = wavelength_deg
        coeff_dict_x["slit-deg"] = slit_deg
        n_coeff = 0
        for i in range(0, order_deg + 1):
            for j in range(0, wavelength_deg + 1):
                for k in range(0, slit_deg + 1):
                    coeff_dict_x[f'c{i}{j}{k}'] = xcoeff[n_coeff]
                    n_coeff += 1

        # SORT Y COEFFICIENT OUTPUT TO WRITE TO FILE
        coeff_dict_y = {}
        coeff_dict_y["axis"] = "y"
        coeff_dict_y["order-deg"] = order_deg
        coeff_dict_y["wavelength-deg"] = wavelength_deg
        coeff_dict_y["slit-deg"] = slit_deg
        n_coeff = 0
        for i in range(0, order_deg + 1):
            for j in range(0, wavelength_deg + 1):
                for k in range(0, slit_deg + 1):
                    coeff_dict_y[f'c{i}{j}{k}'] = ycoeff[n_coeff]
                    n_coeff += 1

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)

        filename = filenamer(
            log=self.log,
            frame=self.pinholeFrame,
            settings=self.settings
        )

        header = copy.deepcopy(self.pinholeFrame.header)
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

        if slit_deg == 0:
            filename = filename.split("ARC")[0] + "DISP_MAP.fits"
            header[kw("PRO_TECH")] = "ECHELLE,PINHOLE"
        else:
            filename = filename.split("ARC")[0] + "2D_MAP.fits"
            header[kw("PRO_TECH")] = "ECHELLE,MULTI-PINHOLE"
        filePath = f"{outDir}/{filename}"
        dataSet = list_of_dictionaries(
            log=self.log,
            listOfDictionaries=[coeff_dict_x, coeff_dict_y]
        )

        # WRITE CSV DATA TO PANDAS DATAFRAME TO ASTROPY TABLE TO FITS
        fakeFile = StringIO(dataSet.csv())
        df = pd.read_csv(fakeFile, index_col=False, na_values=['NA', 'MISSING'])
        fakeFile.close()
        t = Table.from_pandas(df)
        BinTableHDU = fits.table_to_hdu(t)

        header[kw("SEQ_ARM")] = arm
        header[kw("PRO_TYPE")] = "REDUCED"
        header[kw("PRO_CATG")] = f"DISP_TAB_{arm}".upper()
        priHDU = fits.PrimaryHDU(header=header)

        hduList = fits.HDUList([priHDU, BinTableHDU])
        hduList.writeto(filePath, checksum=True, overwrite=True)

        self.log.debug('completed the ``write_map_to_file`` method')
        return filePath

    def calculate_residuals(
            self,
            orderPixelTable,
            xcoeff,
            ycoeff,
            order_deg,
            wavelength_deg,
            slit_deg):
        """*calculate residuals of the polynomial fits against the observed line positions*

        **Key Arguments:**

            - ``orderPixelTable`` -- the predicted line list as a data frame
            - ``xcoeff`` -- the x-coefficients
            - ``ycoeff`` -- the y-coefficients
            - ``order_deg`` -- degree of the order fitting
            - ``wavelength_deg`` -- degree of wavelength fitting
            - ``slit_deg`` -- degree of the slit fitting (False for single pinhole)

        **Return:**
            - ``residuals`` -- combined x-y residuals
            - ``mean`` -- the mean of the combine residuals
            - ``std`` -- the stdev of the combine residuals
            - ``median`` -- the median of the combine residuals
        """
        self.log.debug('starting the ``calculate_residuals`` method')

        arm = self.arm

        # POLY FUNCTION NEEDS A DATAFRAME AS INPUT
        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg, slit_deg=slit_deg, exponents_included=True).poly

        # CALCULATE X & Y RESIDUALS BETWEEN OBSERVED LINE POSITIONS AND POLY
        # FITTED POSITIONS
        orderPixelTable["fit_x"] = poly(orderPixelTable, *xcoeff)
        orderPixelTable["fit_y"] = poly(orderPixelTable, *ycoeff)
        orderPixelTable["residuals_x"] = orderPixelTable[
            "fit_x"] - orderPixelTable["observed_x"]
        orderPixelTable["residuals_y"] = orderPixelTable[
            "fit_y"] - orderPixelTable["observed_y"]

        # CALCULATE COMBINED RESIDUALS AND STATS
        orderPixelTable["residuals_xy"] = np.sqrt(np.square(
            orderPixelTable["residuals_x"]) + np.square(orderPixelTable["residuals_y"]))
        combined_res_mean = np.mean(orderPixelTable["residuals_xy"])
        combined_res_std = np.std(orderPixelTable["residuals_xy"])
        combined_res_median = np.median(orderPixelTable["residuals_xy"])

        self.log.debug('completed the ``calculate_residuals`` method')
        return combined_res_mean, combined_res_std, combined_res_median, orderPixelTable

    def fit_polynomials(
            self,
            orderPixelTable,
            wavelength_deg,
            order_deg,
            slit_deg):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``orderPixelTable`` -- data frame containing order, wavelengths, slit positions and observed pixel positions
            - ``wavelength_deg`` -- degree of wavelength fitting
            - ``order_deg`` -- degree of the order fitting
            - ``slit_deg`` -- degree of the slit fitting (0 for single pinhole)

        **Return:**
            - ``xcoeff`` -- the x-coefficients post clipping
            - ``ycoeff`` -- the y-coefficients post clipping
        """
        self.log.debug('starting the ``fit_polynomials`` method')

        arm = self.arm

        if self.firstGuessMap:
            recipe = "soxs-spatial-solution"
        else:
            recipe = "soxs-disp-solution"

        clippedCount = 1

        # ADD EXPONENTS TO ORDERTABLE UP-FRONT
        for i in range(0, order_deg + 1):
            orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)
        for j in range(0, wavelength_deg + 1):
            orderPixelTable[f"wavelength_pow_{j}"] = orderPixelTable["wavelength"].pow(j)
        for k in range(0, slit_deg + 1):
            orderPixelTable[f"slit_position_pow_{k}"] = orderPixelTable["slit_position"].pow(k)

        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg, slit_deg=slit_deg, exponents_included=True).poly

        clippingSigma = self.settings[
            recipe]["poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.settings[
            recipe]["clipping-iteration-limit"]

        iteration = 0
        while clippedCount > 0 and iteration < clippingIterationLimit:
            iteration += 1
            observed_x = orderPixelTable["observed_x"].to_numpy()
            observed_y = orderPixelTable["observed_y"].to_numpy()
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLYS
            # FIRST X
            coeff = np.ones((order_deg + 1) *
                            (wavelength_deg + 1) * (slit_deg + 1))
            self.log.info("""curvefit x""" % locals())

            xcoeff, pcov_x = curve_fit(
                poly, xdata=orderPixelTable, ydata=observed_x, p0=coeff)

            # NOW Y
            self.log.info("""curvefit y""" % locals())
            ycoeff, pcov_y = curve_fit(
                poly, xdata=orderPixelTable, ydata=observed_y, p0=coeff)

            self.log.info("""calculate_residuals""" % locals())
            mean_res, std_res, median_res, orderPixelTable = self.calculate_residuals(
                orderPixelTable=orderPixelTable,
                xcoeff=xcoeff,
                ycoeff=ycoeff,
                order_deg=order_deg,
                wavelength_deg=wavelength_deg,
                slit_deg=slit_deg)

            # SIGMA-CLIP THE DATA
            self.log.info("""sigma_clip""" % locals())
            masked_residuals = sigma_clip(
                orderPixelTable["residuals_xy"], sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc=mad_std)
            orderPixelTable["residuals_masked"] = masked_residuals.mask
            # RETURN BREAKDOWN OF COLUMN VALUE COUNT
            valCounts = orderPixelTable[
                'residuals_masked'].value_counts(normalize=False)
            if True in valCounts:
                clippedCount = valCounts[True]
            else:
                clippedCount = 0
            print(f'{clippedCount} arc lines where clipped in this iteration of fitting a global dispersion map')

            # REMOVE FILTERED ROWS FROM DATA FRAME
            mask = (orderPixelTable['residuals_masked'] == True)
            orderPixelTable.drop(index=orderPixelTable[
                                 mask].index, inplace=True)

        # a = plt.figure(figsize=(40, 15))
        if arm == "UVB":
            fig = plt.figure(figsize=(6, 13.5), constrained_layout=True)
        else:
            fig = plt.figure(figsize=(6, 11), constrained_layout=True)
        gs = fig.add_gridspec(6, 4)

        # CREATE THE GRID OF AXES
        toprow = fig.add_subplot(gs[0:2, :])
        midrow = fig.add_subplot(gs[2:4, :])
        bottomleft = fig.add_subplot(gs[4:, 0:2])
        bottomright = fig.add_subplot(gs[4:, 2:])

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.flipud(np.rot90(self.pinholeFrame.data, 1))
        toprow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        toprow.set_title(
            "observed arc-line positions (post-clipping)", fontsize=10)

        x = orderPixelTable["observed_x"]
        # x = np.ones(orderPixelTable.shape[0]) * \
        #     self.pinholeFrame.data.shape[1] - orderPixelTable["observed_x"]
        toprow.scatter(orderPixelTable["observed_y"],
                       x, marker='x', c='red', s=4)

        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])
        toprow.set_ylabel("x-axis", fontsize=8)
        toprow.set_xlabel("y-axis", fontsize=8)
        toprow.tick_params(axis='both', which='major', labelsize=9)
        toprow.invert_yaxis()

        midrow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        midrow.set_title(
            "global dispersion solution", fontsize=10)

        xfit = orderPixelTable["fit_x"]
        # xfit = np.ones(orderPixelTable.shape[0]) * \
        #     self.pinholeFrame.data.shape[1] - orderPixelTable["fit_x"]
        midrow.scatter(orderPixelTable["fit_y"],
                       xfit, marker='x', c='blue', s=4)

        # midrow.set_yticklabels([])
        # midrow.set_xticklabels([])
        midrow.set_ylabel("x-axis", fontsize=8)
        midrow.set_xlabel("y-axis", fontsize=8)
        midrow.tick_params(axis='both', which='major', labelsize=9)
        midrow.invert_yaxis()

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomleft.scatter(orderPixelTable["residuals_x"], orderPixelTable[
                           "residuals_y"], alpha=0.4)
        bottomleft.set_xlabel('x residual')
        bottomleft.set_ylabel('y residual')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        hist(orderPixelTable["residuals_xy"], bins='scott', ax=bottomright, histtype='stepfilled',
             alpha=0.7, density=True)
        bottomright.set_xlabel('xy residual')
        bottomright.tick_params(axis='both', which='major', labelsize=9)
        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        fig.suptitle(f"residuals of global dispersion solution fitting - single pinhole\n{subtitle}", fontsize=12)

        # GET FILENAME FOR THE RESIDUAL PLOT
        filename = filenamer(
            log=self.log,
            frame=self.pinholeFrame,
            settings=self.settings
        )
        if self.firstGuessMap:
            filename = filename.split("ARC")[0] + "2D_MAP_RESIDUALS.pdf"
        else:
            filename = filename.split("ARC")[0] + "DISP_MAP_RESIDUALS.pdf"

        # plt.show()
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        filePath = f"{outDir}/{filename}"
        plt.savefig(filePath)

        print(f'\nThe dispersion maps fitted against the observed arc-line positions with a mean residual of {mean_res:2.2f} pixels (stdev = {std_res:2.2f} pixels)')

        self.log.debug('completed the ``fit_polynomials`` method')
        return xcoeff, ycoeff

    def create_placeholder_images(
            self,
            order=False,
            plot=False,
            reverse=False):
        """*create CCDData objects as placeholders to host the 2D images of the wavelength and spatial solutions from dispersion solution map*

        **Key Arguments:**
            - ``order`` -- specific order to generate the placeholder pixels for. Inner-order pixels set to NaN, else set to 0. Default *False* (generate all orders)
            - ``plot`` -- generate plots of placeholder images (for debugging). Default *False*.
            - ``reverse`` -- Inner-order pixels set to 0, else set to NaN (reverse of default output).

        **Return:**
            - ``slitMap`` -- placeholder image to add pixel slit positions to
            - ``wlMap`` -- placeholder image to add pixel wavelength values to

        **Usage:**

        ```python
        slitMap, wlMap = self._create_placeholder_images(order=order)
        ```
        """
        self.log.debug('starting the ``create_placeholder_images`` method')

        kw = self.kw
        dp = self.detectorParams

        # UNPACK THE ORDER TABLE
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTable, extend=0.)

        # CREATE THE IMAGE SAME SIZE AS DETECTOR - NAN INSIDE ORDERS, 0 OUTSIDE
        science_pixels = dp["science-pixels"]
        xlen = science_pixels["columns"]["end"] - \
            science_pixels["columns"]["start"]
        ylen = science_pixels["rows"]["end"] - science_pixels["rows"]["start"]
        xlen, ylen
        if reverse:
            seedArray = np.empty((ylen, xlen))
            seedArray[:] = np.nan
        else:
            seedArray = np.zeros((ylen, xlen))
        wlMap = CCDData(seedArray, unit="adu")
        orderMap = wlMap.copy()
        uniqueOrders = orderPixelTable['order'].unique()
        expandEdges = 0
        for o in uniqueOrders:
            if order and o != order:
                continue
            ycoord = orderPixelTable.loc[
                (orderPixelTable["order"] == o)]["ycoord"]
            xcoord_edgeup = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                "xcoord_edgeup"] + expandEdges
            xcoord_edgelow = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                "xcoord_edgelow"] - expandEdges
            xcoord_edgelow, xcoord_edgeup, ycoord = zip(*[(x1, x2, y) for x1, x2, y in zip(xcoord_edgelow, xcoord_edgeup, ycoord) if x1 > 0 and x1 < wlMap.data.shape[
                                                        1] and x2 > 0 and x2 < wlMap.data.shape[1] and y > 0 and y < wlMap.data.shape[0]])
            if reverse:
                for y, u, l in zip(ycoord, np.ceil(xcoord_edgeup).astype(int), np.floor(xcoord_edgelow).astype(int)):
                    wlMap.data[y, l:u] = 0
                    orderMap.data[y, l:u] = o
            else:
                for y, u, l in zip(ycoord, np.ceil(xcoord_edgeup).astype(int), np.floor(xcoord_edgelow).astype(int)):
                    wlMap.data[y, l:u] = np.NaN
                    orderMap.data[y, l:u] = np.NaN

        # SLIT MAP PLACEHOLDER SAME AS WAVELENGTH MAP PLACEHOLDER
        slitMap = wlMap.copy()

        # PLOT CCDDATA OBJECT
        if plot:
            import matplotlib.pyplot as plt
            rotatedImg = np.rot90(slitMap.data, 1)
            std = np.nanstd(slitMap.data)
            mean = np.nanmean(slitMap.data)
            vmax = mean + 3 * std
            vmin = mean - 3 * std
            plt.figure(figsize=(12, 5))
            plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
                       cmap='gray', alpha=1)
            plt.colorbar()
            plt.xlabel(
                "y-axis", fontsize=10)
            plt.ylabel(
                "x-axis", fontsize=10)
            plt.show()

        self.log.debug('completed the ``create_placeholder_images`` method')
        return slitMap, wlMap, orderMap

    def map_to_image(
            self,
            dispersionMapPath):
        """*convert the dispersion map to images in the detector format showing pixel wavelength values and slit positions*

        **Key Arguments:**
            - ``dispersionMapPath`` -- path to the full dispersion map to convert to images

        **Return:**
            - ``dispersion_image_filePath`` -- path to the FITS image with an extension for wavelength values and another for slit positions

        **Usage:**

        ```python
        mapImagePath = self.map_to_image(dispersionMapPath=mapPath)
        ```
        """
        self.log.debug('starting the ``map_to_image`` method')

        self.dispersionMapPath = dispersionMapPath
        kw = self.kw
        dp = self.detectorParams

        self.map_to_image_displacement_threshold = self.recipeSettings[
            "map_to_image_displacement_threshold"]

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm)

        # GENERATE SLIT ARRAY VALUES - THIS GRID WILL BE THE SAME FOR ALL
        # ORDERS
        slitLength = dp["slit_length"]
        grid_res_slit = self.recipeSettings["grid_res_slit"]
        halfGrid = (slitLength / 2) * 1.1
        self.slitArray = np.arange(-halfGrid, halfGrid +
                                   grid_res_slit, grid_res_slit)

        # CREATE GRIDS FOR ZOOM-IN STAMPS TO FIND CENTRE OF PIXELS
        self.gridSize = self.recipeSettings["zoom_grid_size"]
        grid = np.arange(self.gridSize)
        self.gridSlit = np.tile(grid, (1, self.gridSize))[0]
        self.gridWl = np.repeat(grid, self.gridSize)

        combinedSlitImage = False
        combinedWlImage = False

        # DEFINE AN INPUT ARRAY
        inputArray = [(order, minWl, maxWl) for order, minWl,
                      maxWl in zip(orderNums, waveLengthMin, waveLengthMax)]
        results = fmultiprocess(log=self.log, function=self.order_to_image,
                                inputArray=inputArray, poolSize=False, timeout=3600, turnOffMP=False)

        slitImages = [r[0] for r in results]
        wlImages = [r[1] for r in results]

        # NOTE TO SELF: if having issue with multiprocessing stalling, try and
        # import required modules into the method/function running this
        # fmultiprocess function instead of at the module level

        slitMap, wlMap, orderMap = self.create_placeholder_images(reverse=True)

        combinedSlitImage = Combiner(slitImages)
        combinedSlitImage = combinedSlitImage.sum_combine()
        combinedWlImage = Combiner(wlImages)
        combinedWlImage = combinedWlImage.sum_combine()

        combinedWlImage.data += wlMap.data
        combinedSlitImage.data += wlMap.data

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)

        # GET THE EXTENSION (WITH DOT PREFIX)
        extension = os.path.splitext(dispersionMapPath)[1]
        filename = os.path.basename(
            dispersionMapPath).replace(f"_MAP{extension}", "_MAP_IMAGE.fits")

        dispersion_image_filePath = f"{outDir}/{filename}"
        # WRITE CCDDATA OBJECT TO FILE
        from astropy.io import fits
        primary_hdu = fits.PrimaryHDU(combinedWlImage.data)
        primary_hdu.header['EXTNAME'] = 'WAVELENGTH'
        image_hdu = fits.ImageHDU(combinedSlitImage.data)
        image_hdu.header['EXTNAME'] = 'SLIT'
        image_hdu2 = fits.ImageHDU(orderMap.data)
        image_hdu2.header['EXTNAME'] = 'ORDER'
        hdul = fits.HDUList([primary_hdu, image_hdu, image_hdu2])
        hdul.writeto(dispersion_image_filePath, output_verify='exception',
                     overwrite=True, checksum=False)

        self.log.debug('completed the ``map_to_image`` method')
        return dispersion_image_filePath

    def order_to_image(
            self,
            orderInfo):
        """*convert a single order in the dispersion map to wavelength and slit position images*

        **Key Arguments:**
            - ``orderInfo`` -- tuple containing the order number to generate the images for, the minimum wavelength to consider (from format table) and maximum wavelength to consider (from format table).

        **Return:**
            - ``slitMap`` -- the slit map with order values filled
            - ``wlMap`` -- the wavelengths map with order values filled

        **Usage:**

        ```python
        slitMap, wlMap = self.order_to_image(order=order,minWl=minWl, maxWl=maxWl)
        ```
        """
        self.log.debug('starting the ``order_to_image`` method')

        slitArray = self.slitArray
        gridSlit = self.gridSlit
        gridWl = self.gridWl

        (order, minWl, maxWl) = orderInfo

        slitMap, wlMap, orderMap = self.create_placeholder_images(order=order)

        # GENERATE INITIAL FULL-ORDER WAVELENGTH ARRAY FOR PARTICULAR ORDER
        # (ADD A LITTLE WRIGGLE ROOM AT EACH SIDE OF RANGE)
        wlArray = np.arange(minWl - 20, maxWl + 20, self.recipeSettings[
                            "grid_res_wavelength"])
        # ONE SINGLE-VALUE SLIT ARRAY FOR EVERY WAVELENGTH ARRAY
        bigSlitArray = np.concatenate(
            [np.ones(wlArray.shape[0]) * slitArray[i] for i in range(0, slitArray.shape[0])])
        # NOW THE BIG WAVELEGTH ARRAY
        bigWlArray = np.tile(wlArray, np.shape(slitArray)[0])

        remainingPixels = 1
        iteration = 0
        iterationLimit = 20
        remainingCount = 1
        while remainingPixels and remainingCount and iteration < iterationLimit:
            iteration += 1

            orderPixelTable, remainingCount = self.convert_and_fit(
                order=order, bigWlArray=bigWlArray, bigSlitArray=bigSlitArray, slitMap=slitMap, wlMap=wlMap, iteration=iteration, plots=False)

            if not remainingCount:
                continue

            g = orderPixelTable.groupby(['pixel_x', 'pixel_y', 'order'])
            g = g.agg(["first", np.sum, np.mean, np.std])

            estimatedValues = g.reset_index()

            # SET LOWER LIMIT TO SLIT/WAVELENGTH STD
            limit = self.map_to_image_displacement_threshold / 100
            estimatedValues["wavelength_std"] = np.where(estimatedValues["wavelength"][
                "std"] <= limit, limit, estimatedValues["wavelength"]["std"])
            estimatedValues["slit_position_std"] = np.where(estimatedValues["slit_position"][
                "std"] <= limit, limit, estimatedValues["slit_position"]["std"])

            # SEED GRID ARRAYS ADDED FOR EACH PIXEL
            estimatedValues["gridMeshSlit"] = list(
                np.tile(gridSlit, (len(estimatedValues.index), 1)))
            estimatedValues["gridMeshWl"] = list(
                np.tile(gridWl, (len(estimatedValues.index), 1)))

            # CALCULATE THE DIMENSIONS NEEDED FOR EACH REMAINING PIXEL GRIDS
            estimatedValues["mean_offset_x"] = estimatedValues[
                "fit_x"]["mean"] - (estimatedValues["pixel_x"] + 0.5)
            estimatedValues["mean_offset_y"] = estimatedValues[
                "fit_y"]["mean"] - (estimatedValues["pixel_y"] + 0.5)
            estimatedValues["mean_offset_xy"] = np.sqrt(np.square(
                estimatedValues["mean_offset_x"]) + np.square(estimatedValues["mean_offset_y"]))
            estimatedValues['guess_wavelength'] = np.where(estimatedValues["mean_offset_xy"] <= estimatedValues["residual_xy"]["first"], estimatedValues["wavelength"][
                "mean"], estimatedValues["wavelength"][
                "first"])
            estimatedValues['guess_slit_position'] = np.where(estimatedValues["mean_offset_xy"] <= estimatedValues[
                                                              "residual_xy"]["first"], estimatedValues["slit_position"]["mean"], estimatedValues["slit_position"]["first"])
            estimatedValues['best_offset_x'] = np.where(estimatedValues["mean_offset_xy"] <= estimatedValues[
                "residual_xy"]["first"], abs(estimatedValues["mean_offset_x"]), abs(estimatedValues["residual_x"]["first"]))
            estimatedValues['best_offset_y'] = np.where(estimatedValues["mean_offset_xy"] <= estimatedValues[
                "residual_xy"]["first"], abs(estimatedValues["mean_offset_y"]), abs(estimatedValues["residual_y"]["first"]))
            estimatedValues['offset_std_ratio_x'] = estimatedValues[
                'best_offset_x'] * 2 / estimatedValues["fit_x"]["std"]
            estimatedValues['offset_std_ratio_y'] = estimatedValues[
                'best_offset_y'] * 2 / estimatedValues["fit_y"]["std"]

            estimatedValues["wlArrayMin"] = estimatedValues[
                'guess_wavelength'] - estimatedValues["wavelength_std"] * estimatedValues['offset_std_ratio_y']
            estimatedValues["wlArrayMax"] = estimatedValues[
                'guess_wavelength'] + estimatedValues["wavelength_std"] * estimatedValues['offset_std_ratio_y']
            estimatedValues["slArrayMin"] = estimatedValues[
                'guess_slit_position'] - estimatedValues["slit_position_std"] * estimatedValues['offset_std_ratio_x']
            estimatedValues["slArrayMax"] = estimatedValues[
                'guess_slit_position'] + estimatedValues["slit_position_std"] * estimatedValues['offset_std_ratio_x']
            estimatedValues["wlArray"] = estimatedValues["wlArrayMin"] + estimatedValues[
                "gridMeshWl"] * (estimatedValues["wavelength_std"] * (estimatedValues['offset_std_ratio_y'] * 2) / (self.gridSize - 1))
            estimatedValues["slitArray"] = estimatedValues["slArrayMin"] + estimatedValues[
                "gridMeshSlit"] * (estimatedValues["slit_position_std"] * (estimatedValues['offset_std_ratio_x'] * 2) / (self.gridSize - 1))

            # import sqlite3 as sql
            # # CONNECT TO THE DATABASE
            # conn = sql.connect("/tmp/pandas_export.db")
            # # SEND TO DATABASE
            # estimatedValues.to_sql('my_export_table', con=conn,
            #                        index=False, if_exists='replace')

            # COMBINE ALL PIXEL ARRAYS INTO 2 BIG ARRAYS
            bigWlArray = np.concatenate(estimatedValues["wlArray"].values)
            bigSlitArray = np.concatenate(estimatedValues["slitArray"].values)

            remainingPixels = np.count_nonzero(np.isnan(wlMap.data))

        self.log.debug('completed the ``order_to_image`` method')
        return slitMap, wlMap

    def convert_and_fit(
            self,
            order,
            bigWlArray,
            bigSlitArray,
            slitMap,
            wlMap,
            iteration,
            plots=False):
        """*convert wavelength and slit position grids to pixels*

        **Key Arguments:**
            - ``order`` -- the order being considered
            - ``bigWlArray`` -- 1D array of all wavelengths to be converted
            - ``bigSlitArray`` -- 1D array of all split-positions to be converted (same length as `bigWlArray`)
            - ``slitMap`` -- place-holder image hosting fitted pixel slit-position values
            - ``wlMap`` -- place-holder image hosting fitted pixel wavelength values
            - ``iteration`` -- the iteration index (used for CL reporting)

        **Return:**
            - ``orderPixelTable`` -- dataframe containing unfitted pixel info
            - ``remainingCount`` -- number of remaining pixels in orderTable

        **Usage:**

        ```python
        orderPixelTable = self.convert_and_fit(
                order=order, bigWlArray=bigWlArray, bigSlitArray=bigSlitArray, slitMap=slitMap, wlMap=wlMap)
        ```
        """
        self.log.debug('starting the ``convert_and_fit`` method')

        # CREATE PANDAS DATAFRAME WITH LARGE ARRAYS - ONE ROW PER
        # WAVELENGTH-SLIT GRID CELL
        myDict = {
            "order": np.ones(bigWlArray.shape[0]) * order,
            "wavelength": bigWlArray,
            "slit_position": bigSlitArray
        }
        orderPixelTable = pd.DataFrame(myDict)

        # GET DETECTOR PIXEL POSITIONS FOR ALL WAVELENGTH-SLIT GRID CELLS
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispersionMapPath,
            orderPixelTable=orderPixelTable
        )

        # INTEGER PIXEL VALUES & FIT DISPLACEMENTS FROM PIXEL CENTRES
        orderPixelTable["pixel_x"] = np.floor(orderPixelTable["fit_x"].values)
        orderPixelTable["pixel_y"] = np.floor(orderPixelTable["fit_y"].values)
        orderPixelTable["residual_x"] = orderPixelTable[
            "fit_x"] - (orderPixelTable["pixel_x"] + 0.5)
        orderPixelTable["residual_y"] = orderPixelTable[
            "fit_y"] - (orderPixelTable["pixel_y"] + 0.5)
        orderPixelTable["residual_xy"] = np.sqrt(np.square(
            orderPixelTable["residual_x"]) + np.square(orderPixelTable["residual_y"]))

        # ADD A COUNT COLUMN FOR THE NUMBER OF SMALL SLIT/WL PIXELS FALLING IN
        # LARGE DETECTOR PIXELS
        count = orderPixelTable.groupby(
            ['pixel_x', 'pixel_y']).size().reset_index(name='count')
        orderPixelTable = pd.merge(orderPixelTable, count, how='left', left_on=[
                                   'pixel_x', 'pixel_y'], right_on=['pixel_x', 'pixel_y'])
        orderPixelTable = orderPixelTable.sort_values(
            ['order', 'pixel_x', 'pixel_y', 'residual_xy'])

        # REMVOE LOW COUNT PIXELS AT EDGES OF ORDER
        mask = (orderPixelTable['count'] > 2)
        orderPixelTable = orderPixelTable.loc[mask]

        if plots:
            # PLOT CENTRAL PIXEL
            medX = orderPixelTable.iloc[
                int(len(orderPixelTable.index) / 2)]['pixel_x']
            medY = orderPixelTable.iloc[
                int(len(orderPixelTable.index) / 2)]['pixel_y']
            mask = (orderPixelTable['pixel_x'] == medX) & (
                orderPixelTable['pixel_y'] == medY)
            filteredDf = orderPixelTable.loc[mask]
            fit_x = filteredDf['fit_x'].values - medX - 0.5
            fit_y = filteredDf['fit_y'].values - medY - 0.5

            count = int(len(filteredDf.index))
            windowSize = self.gridSize**2 / (2 * count)
            mean_x = np.mean(filteredDf['fit_x'].values)
            mean_y = np.mean(filteredDf['fit_y'].values)
            mean_x_offset = abs(mean_x - (medX + 0.5))
            mean_y_offset = abs(mean_y - (medY + 0.5))
            mean_xy_offset = np.sqrt(np.square(
                mean_x_offset) + np.square(mean_y_offset))

            if mean_xy_offset < filteredDf['residual_xy'].values[0]:
                best_offset_x = abs(mean_x_offset)
                best_offset_y = abs(mean_y_offset)
            else:
                best_offset_x = abs(filteredDf['residual_x'].values[0])
                best_offset_y = abs(filteredDf['residual_y'].values[0])

            offsetStdRatioX = best_offset_x * 2. / np.std(fit_x)
            offsetStdRatioY = best_offset_y * 2. / np.std(fit_y)

            if filteredDf['residual_x'].min() > 0 or filteredDf['residual_x'].max() < 0 or filteredDf['residual_y'].min() > 0 or filteredDf['residual_y'].max() < 0:
                centred = False
            else:
                centred = True

            from tabulate import tabulate
            print(tabulate(filteredDf, headers='keys', tablefmt='psql'))
            print(f"PIXEL LOCATION STDEV: {np.std(fit_x)}, {np.std(fit_y)}")
            print(f"OFFSET RATIOS: {offsetStdRatioX}, {offsetStdRatioY}")
            print(f"WAVELENGTH/SLIT STDEV: {filteredDf['wavelength'].std()}, {filteredDf['wavelength'].std()}")
            print(f"OFFSET RATIOS: {offsetStdRatioX}, {offsetStdRatioY}")

            # print(filteredDf)
            print(f"Pixel ({medX}, {medY})")
            plt.grid()
            plt.xlim([-0.5, 0.5])
            plt.ylim([-0.5, 0.5])
            plt.plot(fit_x, fit_y, 'o', color='black', label=f"sample (count = {count})", ms=2)
            plt.plot(np.mean(fit_x), np.mean(fit_y), 'x',
                     color='cyan', label=f"mean (centred = {centred})", ms=3)
            plt.plot(fit_x[0], fit_y[0], 'o', color='green', label=f"closest {filteredDf['residual_xy'].values[0]:0.3f}", ms=2)
            plt.plot(0, 0, 'x', color='red', label="pixel centre", ms=2)
            plt.legend(loc=2, prop={'size': 8})
            if mean_xy_offset > filteredDf['residual_xy'].values[0]:
                plt.gca().add_patch(Rectangle((fit_x[0] - np.std(fit_x) * offsetStdRatioX, fit_y[0] - np.std(fit_y) * offsetStdRatioY), np.std(fit_x) * offsetStdRatioX * 2, np.std(fit_y) * offsetStdRatioY * 2,
                                              edgecolor='red',
                                              facecolor='none',
                                              lw=1))
            else:
                plt.gca().add_patch(Rectangle((np.mean(fit_x) - np.std(fit_x) * offsetStdRatioX, np.mean(fit_y) - np.std(fit_y) * offsetStdRatioY), np.std(fit_x) * offsetStdRatioX * 2, np.std(fit_y) * offsetStdRatioY * 2,
                                              edgecolor='cyan',
                                              facecolor='none',
                                              lw=1))

        # FILTER TO WL/SLIT POSITION CLOSE ENOUGH TO CENTRE OF PIXEL
        mask = (orderPixelTable['residual_xy'] <
                self.map_to_image_displacement_threshold)
        # KEEP ONLY VALUES CLOSEST TO CENTRE OF PIXEL
        newPixelValue = orderPixelTable.loc[mask].drop_duplicates(
            subset=['pixel_x', 'pixel_y'], keep="first")
        # REMOVE PIXELS FOUND IN newPixelValue FROM orderPixelTable
        orderPixelTable = newPixelValue[['pixel_x', 'pixel_y']].merge(orderPixelTable, on=[
            'pixel_x', 'pixel_y'], how='right', indicator=True).query('_merge == "right_only"').drop('_merge', 1)
        remainingCount = orderPixelTable.drop_duplicates(
            subset=['pixel_x', 'pixel_y'], keep="first")

        # ADD FITTED PIXELS TO PLACE HOLDER IMAGES
        for xx, yy, wavelength, slit_position in zip(newPixelValue["pixel_x"].values.astype(int), newPixelValue["pixel_y"].values.astype(int), newPixelValue["wavelength"].values, newPixelValue["slit_position"].values):
            try:
                wlMap.data[yy, xx] = np.where(
                    np.isnan(wlMap.data[yy, xx]), wavelength, wlMap.data[yy, xx])
                slitMap.data[yy, xx] = np.where(
                    np.isnan(slitMap.data[yy, xx]), slit_position, slitMap.data[yy, xx])
            except (IndexError):
                # PIXELS OUTSIDE OF DETECTOR EDGES - IGNORE
                pass

        print(f"ORDER {order:02d}, iteration {iteration:02d}. Fit found for {len(newPixelValue.index)} new pixels, {len(remainingCount.index)} image pixel remain to be constrained ({np.count_nonzero(np.isnan(wlMap.data))} nans in place-holder image)")

        if plots:
            # PLOT CCDDATA OBJECT
            rotatedImg = np.rot90(slitMap.data, 3)
            std = np.nanstd(slitMap.data)
            mean = np.nanmean(slitMap.data)
            cmap = cm.gray
            cmap.set_bad(color='#ADD8E6')
            vmax = np.nanmax(slitMap.data)
            vmin = np.nanmin(slitMap.data)
            plt.figure(figsize=(24, 10))
            plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
                       cmap=cmap, alpha=1)
            plt.gca().invert_yaxis()
            plt.colorbar()
            plt.xlabel(
                "y-axis", fontsize=10)
            plt.ylabel(
                "x-axis", fontsize=10)
            plt.show()

        remainingCount = len(remainingCount.index)

        self.log.debug('completed the ``convert_and_fit`` method')
        return orderPixelTable, remainingCount
