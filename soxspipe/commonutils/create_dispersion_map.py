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


class create_dispersion_map(object):
    """
    *detect arc-lines on a pinhole frame to generate a dispersion solution*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``pinholeFrame`` -- the calibrated pinhole frame (single or multi)
        - ``firstGuessMap`` -- the first guess dispersion map from the `soxs_disp_solution` receipe (needed in `soxs_spat_solution` recipe). Deefault *False*.

    **Usage:**

    ```eval_rst
    .. todo::

        - add a tutorial about ``create_dispersion_map`` to documentation
    ```

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
        log.debug("instansiating a new 'create_dispersion_map' object")
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

        **Usage:**


        ```eval_rst
        .. todo::

            - update the package tutorial if needed
        ```
        """
        self.log.debug('starting the ``get`` method')

        # READ PREDICTED LINE POSITIONS FROM FILE
        predictedLines = self.get_predicted_line_list()

        # ROUND OFF THE PIXEL POSITIONS
        observed_x = []
        observed_y = []
        windowSize = self.settings["soxs-disp-solution"]["pixel-window-size"]
        windowHalf = int(windowSize / 2)

        # DETECT THE LINES ON THE PINHILE FRAME
        for x, y in zip(predictedLines["detector_x"].astype(int), predictedLines["detector_y"].astype(int)):
            found_x, found_y = self.detect_pinhole_arc_line(x, y, windowHalf)
            observed_x.append(found_x)
            observed_y.append(found_y)

        # COLLECT INFO FOR DETECTED LINES
        wavelengths = [w for w, x in zip(
            predictedLines["Wavelength"], observed_x) if x != None]
        orders = [w for w, x in zip(
            predictedLines["Order"], observed_x) if x != None]
        observed_x = [x for x in observed_x if x != None]
        observed_y = [y for y in observed_y if y != None]

        order_deg = self.settings["soxs-disp-solution"]["order-deg"]
        wavelength_deg = self.settings["soxs-disp-solution"]["wavelength-deg"]

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        popt_x, popt_y = self.fit_polynomials(
            observed_x=observed_x,
            observed_y=observed_y,
            wavelengths=wavelengths,
            orders=orders,
            order_deg=order_deg,
            wavelength_deg=wavelength_deg
        )

        # WRITE THE MAP TO FILE
        mapPath = self.write_map_to_file(
            popt_x, popt_y, order_deg, wavelength_deg)

        self.log.debug('completed the ``get`` method')
        return mapPath

    def get_predicted_line_list(
            self):
        """*lift the predicted line list from the static calcibrations*

        **Return:**
            - ``predictedLines`` -- a dictionary of lists detailing Wavelength,Order,slit_index,slit_position,detector_x,detector_y
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
        predictedLines = np.genfromtxt(
            predictedLinesFile, delimiter=',', names=True)

        # WANT TO DETERMINE SYSTEMATIC SHIFT IF FIRST GUESS SOLUTION PRESENT
        if self.firstGuessMap:
            slitIndex = int(dp["mid_slit_index"])

            # FILTER THE PREDICTED LINES TO MAKE A DICTIONARY OF PREDICTED LINE
            # POSITIONS PER ORDER
            fPredictedLines = predictedLines[np.where(
                predictedLines["slit_index"] == slitIndex)]
            orders = np.unique(fPredictedLines["Order"])
            slitindexes = np.unique(predictedLines["slit_index"])
            orderWavelengthDict = {
                o: fPredictedLines[np.where(
                    fPredictedLines["Order"] == o)]["Wavelength"] for o in orders}
            orderPredictedXDict = {
                o: fPredictedLines[np.where(
                    fPredictedLines["Order"] == o)]["detector_x"] for o in orders}
            orderPredictedYDict = {
                o: fPredictedLines[np.where(
                    fPredictedLines["Order"] == o)]["detector_y"] for o in orders}

            pixelArrays = dispersion_map_to_pixel_arrays(
                log=self.log,
                dispersionMapPath=self.firstGuessMap,
                orderWavelengthDict=orderWavelengthDict
            )

            # PACK THE PIXEL SHIFTS INTO ORDER ARRAYS
            for o in orders:
                xcoords, ycoords = zip(*[(p[0], p[1]) for p in pixelArrays[o]])
                xshift = orderPredictedXDict[o] - np.array(xcoords)
                yshift = orderPredictedYDict[o] - np.array(ycoords)
                # NOW SHIFT ALL LINES FROM PREDICTED LINE LIST
                for s in slitindexes:
                    # this = predictedLines[np.where(predictedLines["slit_index"] == s)][
                    # np.where(predictedLines["Order"] == o)]["detector_x"] -
                    # xshift
                    predictedLines[
                        np.where((predictedLines["slit_index"] == s) & (predictedLines["Order"] == o))]["detector_x"] = predictedLines[
                        np.where((predictedLines["slit_index"] == s) & (predictedLines["Order"] == o))] - xshift

        self.log.debug('completed the ``get_predicted_line_list`` method')
        return predictedLines

    def detect_pinhole_arc_line(
            self,
            x,
            y,
            windowHalf):
        """*detect the observed position of an arc-line given the predicted pixel positions*

        **Key Arguments:**
            - ``x`` -- the predicted x-positions
            - ``y`` -- the predicted y-positions
            - ``windowHalf`` -- half the width of the image stmap to cut and search (in pixels)

        **Return:**
            - ``observed_x`` -- the observed x-position of the line (None if not found)
            - ``observed_y`` -- the observed y-position of the line (None if not found)
        """
        self.log.debug('starting the ``detect_pinhole_arc_line`` method')

        pinholeFrame = self.pinholeFrame

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
            observed_x = None
            observed_y = None
        # plt.show()

        self.log.debug('completed the ``detect_pinhole_arc_line`` method')
        return observed_x, observed_y

    def write_map_to_file(
            self,
            xcoeff,
            ycoeff,
            order_deg,
            wavelength_deg,
            slit_deg=False):
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
        n_coeff = 0
        for i in range(0, order_deg + 1):
            for j in range(0, wavelength_deg + 1):
                coeff_dict_x[f'c{i}{j}'] = xcoeff[n_coeff]
                n_coeff += 1

        # SORT Y COEFFICIENT OUTPUT TO WRITE TO FILE
        coeff_dict_y = {}
        coeff_dict_y["axis"] = "y"
        coeff_dict_y["order-deg"] = order_deg
        coeff_dict_y["wavelength-deg"] = wavelength_deg
        n_coeff = 0
        for i in range(0, order_deg + 1):
            for j in range(0, wavelength_deg + 1):
                coeff_dict_y[f'c{i}{j}'] = ycoeff[n_coeff]
                n_coeff += 1

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)

        filename = filenamer(
            log=self.log,
            frame=self.pinholeFrame,
            settings=self.settings
        )
        filename = filename.split("ARC")[0] + "DISP_MAP.csv"
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
            observed_x,
            observed_y,
            wavelengths,
            orders,
            xcoeff,
            ycoeff,
            order_deg,
            wavelength_deg,
            slit_deg=False):
        """*calculate residuals of the polynomial fits against the observed line postions*

        **Key Arguments:**
            - ``observed_x`` -- the measured x positions of the lines
            - ``observed_x`` -- the measurd y positions of the lines
            - ``wavelengths`` -- the wavelengths of the lines
            - ``orders`` -- the orders of the lines
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

        # POLY FUNCTION NEEDS A TUPLE AS INPUT
        order_wave = (orders, wavelengths)

        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg).poly

        # CALCULATE X & Y RESIDUALS BETWEEN OBSERVED LINE POSITIONS AND POLY
        # FITTED POSITIONS
        xfit = poly(order_wave, *xcoeff)
        yfit = poly(order_wave, *ycoeff)
        residuals_x = np.asarray(xfit) - np.asarray(observed_x)
        residuals_y = np.asarray(yfit) - np.asarray(observed_y)

        # CALCULATE COMBINED RESIDUALS AND STATS
        combined_res = np.sqrt(np.square(residuals_x) + np.square(residuals_y))
        combined_res_mean = np.mean(combined_res)
        combined_res_std = np.std(combined_res)
        combined_res_median = np.median(combined_res)

        self.log.debug('completed the ``calculate_residuals`` method')
        return combined_res, combined_res_mean, combined_res_std, combined_res_median, residuals_x, residuals_y, xfit, yfit

    def fit_polynomials(
            self,
            observed_x,
            observed_y,
            wavelengths,
            orders,
            order_deg,
            wavelength_deg,
            slit_deg=False):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``observed_x`` -- the measured x positions of the lines
            - ``observed_x`` -- the measurd y positions of the lines
            - ``wavelengths`` -- the wavelengths of the lines
            - ``orders`` -- the orders of the lines
            - ``order_deg`` -- degree of the order fitting
            - ``wavelength_deg`` -- degree of wavelength fitting
            - ``slit_deg`` -- degree of the slit fitting (False for single pinhole)

        **Return:**
            - ``xcoeff`` -- the x-coefficients post clipping
            - ``ycoeff`` -- the y-coefficients post clipping
        """
        self.log.debug('starting the ``fit_polynomials`` method')

        arm = self.arm

        clippedCount = 1

        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg).poly

        clippingSigma = self.settings[
            "soxs-disp-solution"]["poly-fitting-residual-clipping-sigma"]

        clippingIterationLimit = self.settings[
            "soxs-disp-solution"]["clipping-iteration-limit"]

        iteration = 0
        while clippedCount > 0 and iteration < clippingIterationLimit:
            iteration += 1
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLYS
            # FIRST X
            coeff = np.ones((order_deg + 1) * (wavelength_deg + 1))
            xcoeff, pcov_x = curve_fit(poly, xdata=(
                orders, wavelengths), ydata=observed_x, p0=coeff)

            # NOW Y
            ycoeff, pcov_y = curve_fit(poly, xdata=(
                orders, wavelengths), ydata=observed_y, p0=coeff)

            residuals, mean_res, std_res, median_res, residuals_x, residuals_y, xfit, yfit = self.calculate_residuals(
                observed_x=observed_x,
                observed_y=observed_y,
                wavelengths=wavelengths,
                orders=orders,
                xcoeff=xcoeff,
                ycoeff=ycoeff,
                order_deg=order_deg,
                wavelength_deg=wavelength_deg,
                slit_deg=False)

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                residuals, sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc=mad_std)

            # MASK DATA ARRAYS WITH CLIPPED RESIDUAL MASK
            startCount = len(observed_x)
            a = [observed_x, observed_y, wavelengths, orders]
            observed_x, observed_y, wavelengths, orders = [np.ma.compressed(np.ma.masked_array(
                i, masked_residuals.mask)) for i in a]
            clippedCount = startCount - len(observed_x)
            print(f'{clippedCount} arc lines where clipped in this iteration of fitting a global dispersion map')

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
        rotatedImg = np.rot90(self.pinholeFrame.data, 1)
        toprow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        toprow.set_title(
            "observed arc-line positions (post-clipping)", fontsize=10)
        x = np.ones(observed_x.shape) * \
            self.pinholeFrame.data.shape[1] - observed_x
        toprow.scatter(observed_y, x, marker='x', c='red', s=4)
        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])
        toprow.set_ylabel("x-axis", fontsize=8)
        toprow.set_xlabel("y-axis", fontsize=8)
        toprow.tick_params(axis='both', which='major', labelsize=9)

        midrow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        midrow.set_title(
            "global dispersion solution", fontsize=10)
        xfit = np.ones(len(xfit)) * \
            self.pinholeFrame.data.shape[1] - xfit
        midrow.scatter(yfit, xfit, marker='x', c='blue', s=4)
        # midrow.set_yticklabels([])
        # midrow.set_xticklabels([])
        midrow.set_ylabel("x-axis", fontsize=8)
        midrow.set_xlabel("y-axis", fontsize=8)
        midrow.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomleft.scatter(residuals_x, residuals_y, alpha=0.4)
        bottomleft.set_xlabel('x residual')
        bottomleft.set_ylabel('y residual')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        hist(residuals, bins='scott', ax=bottomright, histtype='stepfilled',
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
