#!/usr/bin/env python
# encoding: utf-8
"""
*using a fully-illuminated slit flat frame detect and record the order-edges*

:Author:
    David Young & Marco Landoni

:Date Created:
    September 18, 2020
"""
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
import numpy as np
import matplotlib.pyplot as plt
from soxspipe.commonutils.toolkit import cut_image_slice
from soxspipe.commonutils import _base_detect
from scipy.signal import medfilt
import matplotlib.pyplot as plt
import random
from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial
from os.path import expanduser
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils.toolkit import unpack_order_table
import collections
import unicodecsv as csv
from fundamentals.renderer import list_of_dictionaries


class detect_order_edges(_base_detect):
    """
    *using a fully-illuminated slit flat frame detect and record the order-edges*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``flatFrame`` -- the flat frame to detect the order edges on
        - ``orderCentreTable`` -- the order centre table

    **Usage:**

    ```eval_rst
    .. todo::

        - add a tutorial about ``detect_order_edges`` to documentation
    ```

    ```python
    from soxspipe.commonutils import detect_order_edges
    edges = detect_order_edges(
        log=log,
        flatFrame=flatFrame,
        orderCentreTable=orderCentreTable,
        settings=settings,
    )
    edges.get()
    ```
    """

    def __init__(
            self,
            log,
            flatFrame,
            orderCentreTable,
            settings=False,

    ):
        self.log = log
        log.debug("instansiating a new 'detect_order_edges' object")
        self.settings = settings
        self.orderCentreTable = orderCentreTable
        self.flatFrame = flatFrame
        # xt-self-arg-tmpx

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = flatFrame.header[kw("SEQ_ARM")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        # DEG OF THE POLYNOMIALS TO FIT THE ORDER CENTRE LOCATIONS
        self.polyDeg = self.settings[
            "soxs-mflat"]["poly-deg"]

        return None

    def get(self):
        """
        *get the detect_order_edges object*

        **Return:**
            - ``orderTablePath`` -- path to the new order table
        """
        self.log.debug('starting the ``get`` method')

        orderTablePath = None

        # GET PARAMETERS FROM SETTINGS
        sliceLength = int(
            self.settings["soxs-mflat"]["slice-length-for-edge-detection"])
        sliceWidth = int(
            self.settings["soxs-mflat"]["slice-width-for-edge-detection"])
        minThresholdPercenage = int(
            self.settings["soxs-mflat"]["min-percentage-threshold-for-edge-detection"]) / 100
        maxThresholdPercenage = int(
            self.settings["soxs-mflat"]["max-percentage-threshold-for-edge-detection"]) / 100

        # UNPACK THE ORDER TABLE (CENTRE LOCATION ONLY AT THIS STAGE)
        orders, orderCentres, orderLimits = unpack_order_table(
            log=self.log, orderTablePath=self.orderCentreTable)

        # index = -1
        # orders, orderCentres, orderLimits = [orders[index]], [
        #     orderCentres[index]], [orderLimits[index]]

        # GENERATE ARRAYS OF PIXEL LOCATIONS OF UPPER AND LOWER ORDER EDGES
        orderEdgeMax = []
        orderEdgeMin = []
        halfSlice = sliceLength / 2
        for xcoords, ycoords in orderCentres:

            minThreshold, maxThreshold = self.determine_order_flux_threshold(
                xcoords, ycoords, minThresholdPercenage, maxThresholdPercenage, sliceWidth, sliceLength)

            # FOR EACH PIXEL COORDINATE
            xmaxArray = []
            xminArray = []
            for x, y in zip(xcoords, ycoords):
                xmax, xmin = self.determine_lower_upper_edge_pixel_positions(
                    x, y, minThreshold, maxThreshold, sliceWidth, sliceLength)
                if xmax == None:
                    continue
                xmaxArray.append((xmax + int(x - halfSlice), y))
                xminArray.append((xmin + int(x - halfSlice), y))
            orderEdgeMax.append(xmaxArray)
            orderEdgeMin.append(xminArray)

        # PARAMETERS & VARIABLES FOR FITTING EDGES
        orderMaxLocations = {}
        orderMinLocations = {}
        allResiduals = []
        allXfit = []
        allXcoords = []
        allYcoords = []
        ylims = []
        clippingSigma = self.settings[
            "soxs-mflat"]["poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.settings[
            "soxs-mflat"]["clipping-iteration-limit"]

        # FIRST FIT UPPER EDGE PIXEL-POSTIONS WITH POLYNOMIAL
        for order, pixelArray in zip(orders, orderEdgeMax):
            # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
            coeff, residuals, xfit, xcoords, ycoords = self.fit_polynomial(
                pixelArray=pixelArray,
                clippingSigma=clippingSigma,
                clippingIterationLimit=clippingIterationLimit
            )
            allResiduals.extend(residuals)
            allXfit.extend(xfit)
            allXcoords.extend(xcoords)
            allYcoords.extend(ycoords)
            orderMaxLocations[order] = coeff
            print()

        # NOW FIT LOWER EDGE PIXEL-POSTIONS WITH POLYNOMIAL
        for order, pixelArray in zip(orders, orderEdgeMin):
            # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
            coeff, residuals, xfit, xcoords, ycoords = self.fit_polynomial(
                pixelArray=pixelArray,
                clippingSigma=clippingSigma,
                clippingIterationLimit=clippingIterationLimit
            )
            allResiduals.extend(residuals)
            allXfit.extend(xfit)
            allXcoords.extend(xcoords)
            allYcoords.extend(ycoords)
            orderMinLocations[order] = coeff
            print()

        # GENERATE AN OUTPUT PLOT OF RESULTS AND FITTING RESIDUALS
        self.plot_results(
            allResiduals=allResiduals,
            allXfit=allXfit,
            allXcoords=allXcoords,
            allYcoords=allYcoords,
            orderMaxLocations=orderMaxLocations,
            orderMinLocations=orderMinLocations,
            ylims=orderLimits
        )

        # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE
        order_table_path = self.write_order_table_to_file(
            upperEdges=orderMaxLocations, lowerEdges=orderMinLocations,  ylims=orderLimits)

        self.log.debug('completed the ``get`` method')
        return orderTablePath

    def plot_results(
            self,
            allResiduals,
            allXfit,
            allXcoords,
            allYcoords,
            orderMaxLocations,
            orderMinLocations,
            ylims):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**
            - ``allResiduals`` -- list of all residuals
            - ``allXfit`` -- list of all fitted x-positions
            - ``allXcoords`` -- cleaned list of all guassian x-pixel positions
            - ``allYcoords`` -- cleaned list of all guassian y-pixel positions
            - ``orderMaxLocations`` -- dictionary of upper order-location polynomial coeff
            - ``orderMinLocations`` -- dictionary of lower order-location polynomial coeff
            - ``ylims`` -- the limits of the fit

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
        rotatedImg = np.flipud(np.rot90(self.flatFrame.data, 1))
        # rotatedImg = self.flatFrame.data
        std = np.std(self.flatFrame.data)
        mean = np.mean(self.flatFrame.data)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        toprow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=1)
        toprow.invert_yaxis()
        toprow.set_title(
            "upper and lower order edge detections", fontsize=10)
        toprow.scatter(allYcoords, allXcoords, marker='x', c='red', s=0.2)
        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])
        toprow.set_ylabel("x-axis", fontsize=8)
        toprow.set_xlabel("y-axis", fontsize=8)
        toprow.tick_params(axis='both', which='major', labelsize=9)

        midrow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=0.9)
        midrow.invert_yaxis()
        midrow.set_title(
            "order-location fit solutions", fontsize=10)
        ylinelist = np.arange(0, self.flatFrame.data.shape[0], 3)
        poly = chebyshev_xy_polynomial(
            log=self.log, deg=self.polyDeg).poly

        for ylim, mincoeff, maxcoeff in zip(ylims, orderMinLocations.values(), orderMaxLocations.values()):
            ylinelist = np.arange(ylim[0], ylim[1], 3)
            xfitmin = poly(ylinelist, *mincoeff)
            xfitmax = poly(ylinelist, *maxcoeff)
            xfitmin, xfitmax, ylinelist = zip(
                *[(x1, x2, y) for x1, x2, y in zip(xfitmin, xfitmax, ylinelist) if x1 > 0 and x1 < (self.flatFrame.data.shape[1]) - 10 and x2 > 0 and x2 < (self.flatFrame.data.shape[1]) - 10])
            # midrow.plot(ylinelist, xfit)
            midrow.fill_between(ylinelist, xfitmin, xfitmax, alpha=0.4)

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
        bottomleft.scatter(allXcoords, allResiduals, alpha=0.2, s=0.05)
        bottomleft.set_xlabel('x pixel position')
        bottomleft.set_ylabel('x residual')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomright.scatter(allYcoords, allResiduals, alpha=0.2, s=0.2)
        bottomright.set_xlabel('y pixel position')
        bottomright.tick_params(axis='both', which='major', labelsize=9)
        # bottomright.set_ylabel('x residual')
        bottomright.set_yticklabels([])

        mean_res = np.mean(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        fig.suptitle(f"detection of order-edge locations - flat-frame\n{subtitle}", fontsize=12)

        # plt.show()
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        filePath = f"{outDir}/flat_{arm}_order_edge_residuals.pdf"
        plt.savefig(filePath, dpi=720)

        self.log.debug('completed the ``plot_results`` method')
        return filePath

    def determine_order_flux_threshold(
            self,
            xcoords,
            ycoords,
            minThresholdPercenage,
            maxThresholdPercenage,
            sliceWidth,
            sliceLength):
        """*determine the flux threshold at the central column of each order*

        **Key Arguments:**
            - ``xcoords`` -- the x-coordinate array for the central trace of the order
            - ``ycoords`` -- the y-coordinate array for the central trace of the order
            - ``minThresholdPercenage`` -- the lowest percentage threshold to try and detect order edges at
            - ``maxThresholdPercenage`` -- the highest percentage threshold to try and detect order edges at
            - ``sliceWidth`` -- the width of the slice to take from the image (pixels)
            - ``sliceLength`` -- the length of the slice to take from the image, in the cross-dispersion direction (pixels)

        **Return:**
            - ``minThreshold`` -- minimum threshold value in flux
            - ``maxThreshold`` -- maximum threshold value in flux
        """
        self.log.debug(
            'starting the ``determine_order_flux_threshold`` method')

        # DETERMINE THE FLUX THRESHOLD FROM THE CENTRAL COLUMN
        # CUT A MEDIAN COLLAPSED SLICE
        index = int(len(xcoords) / 2)
        x = xcoords[index]
        y = ycoords[index]
        slice = cut_image_slice(log=self.log, frame=self.flatFrame,
                                width=sliceWidth, length=sliceLength, x=x, y=y, median=True, plot=False)
        # SMOOTH WITH A MEDIAN FILTER
        medSlide = medfilt(slice, 9)
        # DETERMINE THRESHOLD FLUX VALUE
        maxvalue = np.max(
            medSlide[int(len(medSlide) / 2 - 8):int(len(medSlide) / 2 + 8)])
        minvalue = np.min(medSlide)
        minThreshold = minvalue + \
            (maxvalue - minvalue) * minThresholdPercenage
        maxThreshold = minvalue + \
            (maxvalue - minvalue) * maxThresholdPercenage

        # SANITY CHECK PLOT OF CROSS-SECTION
        if 1 == 0:
            # CHECK THE SLICE POINTS IF NEEDED
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko', alpha=0.5)
            plt.plot(x, medSlide, 'rx', alpha=0.8)
            # plt.hlines(maxvalue, 0, len(slice), label='max')
            # plt.hlines(minvalue, 0, len(slice), label='min')
            plt.hlines(minThreshold, 0, len(slice),
                       label=f'threshold {minThreshold:0.2f},  {maxThreshold:0.2f}',  colors='red')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            plt.legend()
            plt.show()

        self.log.debug(
            'completed the ``determine_order_flux_threshold`` method')
        return minThreshold, maxThreshold

    def determine_lower_upper_edge_pixel_positions(
            self,
            x,
            y,
            minThreshold,
            maxThreshold,
            sliceWidth,
            sliceLength):
        """*from a pixel postion somewhere on the trace of the order centre, return the lower and upper edges of the order*

        **Key Arguments:**
            - ``x`` -- the pixel postion x value
            - ``y`` -- the pixel position y value
            - ``minThreshold`` -- minimum threshold value in flux
            - ``maxThreshold`` -- maximum threshold value in flux
            - ``sliceWidth`` -- the width of the slice to take from the image (pixels)
            - ``sliceLength`` -- the length of the slice to take from the image, in the cross-dispersion direction (pixels)

        **Return:**
            - ``xmin`` -- the lower edge pixel position x value
            - ``xmax`` -- the upper edge pixel position x value
        """
        self.log.debug(
            'starting the ``determine_lower_upper_edge_limits`` method')

        threshold = minThreshold

        # CUT A MEDIAN COLLAPSED SLICE
        slice = cut_image_slice(log=self.log, frame=self.flatFrame,
                                width=sliceWidth, length=sliceLength, x=x, y=y, median=True, plot=False)
        if slice is None:
            return None, None

        # SMOOTH WITH A MEDIAN FILTER
        medSlide = medfilt(slice, 9)

        # FIND FIRST TIME VALUE DROPS BELOW THRESHOLD WORKING OUT FROM
        # MIDDLE OF SLICE
        middle = int(len(medSlide) / 2)
        firstHalf = medSlide[0:middle]
        secondHalf = medSlide[middle:]

        # ITERATE UP TO MAX THRESHOLD
        hit = False
        while hit == False and threshold < maxThreshold:
            try:
                xmaxguess = np.where(secondHalf < threshold)[
                    0][0] + middle
                xminguess = np.where(firstHalf < threshold)[0][-1]
                hit = True
            except:
                threshold = threshold * 1.1

        # IF WE STILL DIDN'T GET A HIT THEN REJECT
        if hit == False:
            return None, None

        # REPORT THE EXACT PIXEL POSTION AT THE FLUX THRESHOLD
        xmax = xmaxguess - \
            (threshold - medSlide[xmaxguess]) / \
            (medSlide[xmaxguess - 1] - medSlide[xmaxguess])

        xmin = xminguess + \
            (threshold - medSlide[xminguess]) / \
            (medSlide[xminguess + 1] - medSlide[xminguess])

        # IF THE WIDTH BETWEEN MIN AND MAX IS TOO SMALL THEN REJECT
        if xmax - xmin < 10:
            return None, None

        # SANITY CHECK PLOT OF CROSS-SECTION
        if 1 == 0 and random.randint(1, 101) < 2:
            # CHECK THE SLICE POINTS IF NEEDED
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko', alpha=0.5)
            plt.plot(x, medSlide, 'rx', alpha=0.8)
            plt.plot(xmin, threshold, 'ro',
                     alpha=0.8, label="order edge")
            plt.plot(xmax, threshold, 'ro', alpha=0.8)
            # plt.hlines(maxvalue, 0, len(slice), label='max')
            # plt.hlines(minvalue, 0, len(slice), label='min')
            plt.hlines(threshold, 0, len(slice),
                       label='threshold',  colors='red')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            plt.legend()
            plt.show()

        self.log.debug(
            'completed the ``determine_lower_upper_edge_limits`` method')
        return xmin, xmax

    def write_order_table_to_file(
            self,
            upperEdges,
            lowerEdges,
            ylims):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**
            - ``upperEdges`` -- dictionary of the order upper-edge location coefficients
            - ``lowerEdges`` -- dictionary of the order lower-edge location coefficients
            - ``ylims`` -- the y-min, y-max limits of the order

        **Return:**
            - ``order_table_path`` -- path to the order table file
        """
        self.log.debug('starting the ``write_order_table_to_file`` method')

        arm = self.arm

        # SORT UPPER COEFFICIENT OUTPUT TO WRITE TO FILE
        listOfDictionaries = []
        for k, v in upperEdges.items():
            orderDict = collections.OrderedDict(sorted({}.items()))
            orderDict["order"] = k
            orderDict["degy_edge"] = self.polyDeg
            n_coeff = 0
            for i in range(0, self.polyDeg + 1):
                orderDict[f'EDGEUP_c{i}'] = v[n_coeff]
                n_coeff += 1
            listOfDictionaries.append(orderDict)

        # SORT LOWER COEFFICIENT OUTPUT TO WRITE TO FILE
        for k, v in lowerEdges.items():
            orderDict = collections.OrderedDict(sorted({}.items()))
            order = k
            for l in listOfDictionaries:
                if l["order"] == order:
                    n_coeff = 0
                    for i in range(0, self.polyDeg + 1):
                        l[f'EDGELOW_c{i}'] = v[n_coeff]
                        n_coeff += 1

        # RE-ADD CENTRAL ORDERS
        with open(self.orderCentreTable, 'rb') as csvFile:
            csvReader = csv.DictReader(
                csvFile, dialect='excel', delimiter=',', quotechar='"')
            for row in csvReader:
                order = int(row["order"])
                for l in listOfDictionaries:
                    if l["order"] == order:
                        l["degy_cent"] = row["degy_cent"]
                        l["ymin"] = row["ymin"]
                        l["ymax"] = row["ymax"]
                        for k, v in row.items():
                            if "CENT_" in k:
                                l[k] = row[k]

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        order_table_path = f"{outDir}/order_locations_{arm}.csv"
        dataSet = list_of_dictionaries(
            log=self.log,
            listOfDictionaries=listOfDictionaries
        )
        csvData = dataSet.csv(filepath=order_table_path)

        self.log.debug('completed the ``write_order_table_to_file`` method')
        return order_table_path

    # use the tab-trigger below for new method
    # xt-class-method
