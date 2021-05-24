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
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils import detector_lookup
from astropy.stats import sigma_clipped_stats
from photutils import datasets
from photutils import DAOStarFinder
from scipy.optimize import curve_fit
from fundamentals.renderer import list_of_dictionaries
from os.path import expanduser
import matplotlib.pyplot as plt
from os.path import expanduser
from astropy.stats import sigma_clip, mad_std
import numpy as np
import math
from photutils.utils import NoDetectionsWarning
import warnings
from astropy.visualization import hist
from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
import pandas as pd
from tabulate import tabulate


class create_dispersion_map(object):
    """
    *detect arc-lines on a pinhole frame to generate a dispersion solution*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``pinholeFrame`` -- the calibrated pinhole frame (single or multi)
        - ``firstGuessMap`` -- the first guess dispersion map from the `soxs_disp_solution` recipe (needed in `soxs_spat_solution` recipe). Default *False*.

    **Usage:**

    ```python
    from soxspipe.commonutils import create_dispersion_map
    mapPath = create_dispersion_map(
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
            firstGuessMap=False
    ):
        self.log = log
        log.debug("instantiating a new 'create_dispersion_map' object")
        self.settings = settings
        self.pinholeFrame = pinholeFrame
        self.firstGuessMap = firstGuessMap

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        self.kw = kw
        self.arm = pinholeFrame.header[kw("SEQ_ARM")]

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
            recipe = "soxs-spatial-solution"
            slit_deg = self.settings[recipe]["slit-deg"]
        else:
            recipe = "soxs-disp-solution"
            slit_deg = 0

        # READ PREDICTED LINE POSITIONS FROM FILE - RETURNED AS DATAFRAME
        orderPixelTable = self.get_predicted_line_list()

        # GET THE WINDOW SIZE FOR ATTEMPTING TO DETECT LINES ON FRAME
        windowSize = self.settings[recipe]["pixel-window-size"]
        self.windowHalf = int(windowSize / 2)

        # DETECT THE LINES ON THE PINHILE FRAME AND
        # ADD OBSERVED LINES TO DATAFRAME
        orderPixelTable = orderPixelTable.apply(
            self.detect_pinhole_arc_line, axis=1)

        # DROP MISSING VALUES
        orderPixelTable.dropna(axis='index', how='any', subset=[
            'observed_x'], inplace=True)

        order_deg = self.settings[recipe]["order-deg"]
        wavelength_deg = self.settings[
            recipe]["wavelength-deg"]

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

        self.log.debug('completed the ``get`` method')
        return mapPath

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
        if kw('WIN_BINX') in pinholeFrame.header:
            binx = int(self.pinholeFrame.header[kw('WIN_BINX')])
            biny = int(self.pinholeFrame.header[kw('WIN_BINY')])
        else:
            binx = 1
            biny = 1

        # READ THE FILE
        home = expanduser("~")
        calibrationRootPath = self.settings[
            "calibration-data-root"].replace("~", home)
        predictedLinesFile = calibrationRootPath + "/" + dp["predicted pinhole lines"][frameTech][f"{binx}x{biny}"]

        # READ CSV FILE TO PANDAS DATAFRAME
        orderPixelTable = pd.read_csv(predictedLinesFile)

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
        if slit_deg == 0:
            filename = filename.split("ARC")[0] + "DISP_MAP.csv"
        else:
            filename = filename.split("ARC")[0] + "2D_MAP.csv"
        filePath = f"{outDir}/{filename}"
        dataSet = list_of_dictionaries(
            log=self.log,
            listOfDictionaries=[coeff_dict_x, coeff_dict_y]
        )
        csvData = dataSet.csv(filepath=filePath)

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
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg, slit_deg=slit_deg).poly

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

        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg, slit_deg=slit_deg).poly

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

        # CREATE THE GID OF AXES
        toprow = fig.add_subplot(gs[0:2, :])
        midrow = fig.add_subplot(gs[2:4, :])
        bottomleft = fig.add_subplot(gs[4:, 0:2])
        bottomright = fig.add_subplot(gs[4:, 2:])

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.flipud(np.rot90(self.pinholeFrame.data, 1))
        toprow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        toprow.set_title(
            "observed arc-line positions (post-clipping)", fontsize=10)

        x = np.ones(orderPixelTable.shape[0]) * \
            self.pinholeFrame.data.shape[1] - orderPixelTable["observed_x"]
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

        xfit = np.ones(orderPixelTable.shape[0]) * \
            self.pinholeFrame.data.shape[1] - orderPixelTable["fit_x"]
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

        print(f'\nThe dispersion maps fitted against the observed arc-line positions with a mean residual of {mean_res:2.2f} pixels (stdev = {std_res:2.2f} pixles)')

        self.log.debug('completed the ``fit_polynomials`` method')
        return xcoeff, ycoeff

    # use the tab-trigger below for new method
    # xt-class-method
