#!/usr/bin/env python
# encoding: utf-8
"""
*using a fully-illuminated slit flat frame detect and record the order-edges*

:Author:
    David Young & Marco Landoni

:Date Created:
    September 18, 2020
"""
from datetime import datetime, date, time
import pandas as pd
from soxspipe.commonutils.filenamer import filenamer
from fundamentals.renderer import list_of_dictionaries
import unicodecsv as csv
import collections
from soxspipe.commonutils.toolkit import unpack_order_table
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from os.path import expanduser
from soxspipe.commonutils.polynomials import chebyshev_order_xy_polynomials
import random
from scipy.signal import medfilt
from soxspipe.commonutils import _base_detect
from soxspipe.commonutils.toolkit import cut_image_slice
import matplotlib.pyplot as plt
import numpy as np
from fundamentals import tools
from builtins import object
import sys
import os
from io import StringIO
os.environ['TERM'] = 'vt100'


class detect_order_edges(_base_detect):
    """
    *using a fully-illuminated slit flat frame detect and record the order-edges*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``flatFrame`` -- the flat frame to detect the order edges on
        - ``orderCentreTable`` -- the order centre table
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
        - ``verbose`` -- verbose. True or False. Default *False*
        - ``qcTable`` -- the data frame to collect measured QC metrics 
        - ``productsTable`` -- the data frame to collect output products

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
        recipeName="soxs-mflat",
        verbose=False,
        qcTable=False,
        productsTable=False
    )
    productsTable, qcTable = edges.get()
    ```
    """

    def __init__(
            self,
            log,
            flatFrame,
            orderCentreTable,
            settings=False,
            recipeName="soxs-mflat",
            verbose=False,
            qcTable=False,
            productsTable=False
    ):
        self.log = log
        log.debug("instansiating a new 'detect_order_edges' object")
        self.settings = settings
        self.recipeName = recipeName
        if recipeName:
            self.recipeSettings = settings[recipeName]
        else:
            self.recipeSettings = False

        self.orderCentreTable = orderCentreTable
        self.flatFrame = flatFrame
        self.verbose = verbose
        self.qc = qcTable
        self.products = productsTable

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = flatFrame.header[kw("SEQ_ARM")]
        self.dateObs = flatFrame.header[kw("DATE_OBS")]

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
        *get the detect_order_edges object*

        **Return:**
            - ``orderTablePath`` -- path to the new order table
        """
        self.log.debug('starting the ``get`` method')

        print("\n# DETECTING THE ORDER EDGES FROM MASTER-FLAT FRAME")

        orderTablePath = None

        # GET PARAMETERS FROM SETTINGS
        self.sliceLength = int(
            self.recipeSettings["slice-length-for-edge-detection"])
        self.sliceWidth = int(
            self.recipeSettings["slice-width-for-edge-detection"])
        self.minThresholdPercenage = int(
            self.recipeSettings["min-percentage-threshold-for-edge-detection"]) / 100
        self.maxThresholdPercenage = int(
            self.recipeSettings["max-percentage-threshold-for-edge-detection"]) / 100

        # UNPACK THE ORDER TABLE (CENTRE LOCATION ONLY AT THIS STAGE)
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderCentreTable)

        # ADD MIN AND MAX FLUX THRESHOLDS TO ORDER TABLE
        print("\tDETERMINING ORDER FLUX THRESHOLDS")
        orderMetaTable["maxThreshold"] = np.nan
        orderMetaTable["minThreshold"] = np.nan
        orderMetaTable = orderMetaTable.apply(
            self.determine_order_flux_threshold, axis=1, orderPixelTable=orderPixelTable)

        # ADD THRESHOLDS TO orderPixelTable
        orderPixelTable["maxThreshold"] = np.nan
        orderPixelTable["minThreshold"] = np.nan
        uniqueOrders = orderMetaTable['order'].unique()
        for o in uniqueOrders:
            orderPixelTable.loc[(orderPixelTable["order"] == o), ["minThreshold", "maxThreshold"]] = orderMetaTable.loc[
                (orderMetaTable["order"] == o), ["minThreshold", "maxThreshold"]].values

        print("\tMEASURING PIXEL-POSITIONS AT ORDER-EDGES WHERE FLUX THRESHOLDS ARE MET")
        orderPixelTable["xcoord_upper"] = np.nan
        orderPixelTable["xcoord_lower"] = np.nan
        orderPixelTable = orderPixelTable.apply(
            self.determine_lower_upper_edge_pixel_positions, axis=1)

        # DROP ROWS WITH NAN VALUES
        orderPixelTable.dropna(axis='index', how='any',
                               subset=['xcoord_upper'], inplace=True)
        orderPixelTable.dropna(axis='index', how='any',
                               subset=['xcoord_lower'], inplace=True)

        for o in uniqueOrders:
            orderMetaTable.loc[(orderMetaTable["order"] == o), "ymin"] = np.nanmin(orderPixelTable.loc[(orderPixelTable["order"] == o), ["ycoord"]].values)
            orderMetaTable.loc[(orderMetaTable["order"] == o), "ymax"] = np.nanmax(orderPixelTable.loc[(orderPixelTable["order"] == o), ["ycoord"]].values)

        # REDEFINE UNIQUE ORDERS IN CASE ONE OR MORE IS COMPLETELY MISSING
        uniqueOrders = orderPixelTable['order'].unique()

        # SETUP EXPONENTS AHEAD OF TIME - SAVES TIME ON POLY FITTING
        for i in range(0, self.yDeg + 1):
            orderPixelTable[f"y_pow_{i}"] = orderPixelTable["ycoord"].pow(i)
        for i in range(0, self.orderDeg + 1):
            orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        print("\tFITTING POLYNOMIALS TO MEASURED PIXEL-POSITIONS AT UPPER ORDER-EDGES\n")
        upperCoeff, orderPixelTable = self.fit_global_polynomial(
            pixelList=orderPixelTable,
            yCol="ycoord",
            xCol="xcoord_upper",
            orderCol="order",
            y_deg=self.yDeg,
            order_deg=self.orderDeg,
            exponents_included=True
        )

        # RENAME SOME INDIVIDUALLY
        orderPixelTable.rename(columns={
            "x_fit": "xcoord_upper_fit", "x_fit_res": "xcoord_upper_fit_res"}, inplace=True)

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        print("\tFITTING POLYNOMIALS TO MEASURED PIXEL-POSITIONS AT LOWER ORDER-EDGES\n")
        lowerCoeff, orderPixelTable = self.fit_global_polynomial(
            pixelList=orderPixelTable,
            yCol="ycoord",
            xCol="xcoord_lower",
            orderCol="order",
            y_deg=self.yDeg,
            order_deg=self.orderDeg,
            exponents_included=True
        )

        # RENAME SOME INDIVIDUALLY
        orderPixelTable.rename(columns={
            "x_fit": "xcoord_lower_fit", "x_fit_res": "xcoord_lower_fit_res"}, inplace=True)

        # orderLocations[o] = coeff
        coeff_dict = {
            "degorder_edgeup": self.orderDeg,
            "degy_edgeup": self.yDeg,
            "degorder_edgelow": self.orderDeg,
            "degy_edgelow": self.yDeg
        }
        n_coeff = 0
        for i in range(0, self.orderDeg + 1):
            for j in range(0, self.yDeg + 1):
                coeff_dict[f'edgeup_c{i}{j}'] = upperCoeff[n_coeff]
                coeff_dict[f'edgelow_c{i}{j}'] = lowerCoeff[n_coeff]
                n_coeff += 1
        coeffColumns = coeff_dict.keys()
        dataSet = list_of_dictionaries(
            log=self.log,
            listOfDictionaries=[coeff_dict]
        )

        # WRITE CSV DATA TO PANDAS DATAFRAME TO ASTROPY TABLE TO FITS
        fakeFile = StringIO(dataSet.csv())
        orderEdgePolyTable = pd.read_csv(fakeFile, index_col=False, na_values=['NA', 'MISSING'])
        fakeFile.close()

        # MERGE DATAFRAMES
        cols_to_use = orderEdgePolyTable.columns.difference(orderPolyTable.columns)
        orderPolyTable = orderPolyTable.join(
            orderEdgePolyTable[cols_to_use])

        # GENERATE AN OUTPUT PLOT OF RESULTS AND FITTING RESIDUALS
        print("\tMEASURING AND PLOTTING RESIDUALS OF FITS")
        allResiduals = np.concatenate((orderPixelTable[
            'xcoord_lower_fit_res'], orderPixelTable['xcoord_upper_fit_res']))
        plotPath = self.plot_results(
            orderPixelTable=orderPixelTable,
            orderPolyTable=orderPolyTable,
            orderMetaTable=orderMetaTable
        )

        # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE
        orderTablePath = self.write_order_table_to_file(
            frame=self.flatFrame, orderPolyTable=orderPolyTable, orderMetaTable=orderMetaTable)
        mean_res = np.mean(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        orderTablePath = os.path.abspath(orderTablePath)
        plotPath = os.path.abspath(plotPath)
        orderTableName = os.path.basename(orderTablePath)
        plotName = os.path.basename(plotPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # RECORD PRODUCTS AND QCs
        if not isinstance(self.products, bool):
            self.products = self.products.append({
                "soxspipe_recipe": self.recipeName,
                "product_label": "ORDER_LOC",
                "product_desc": "table of coefficients from polynomial fits to order locations",
                "file_name": orderTableName,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "file_path": orderTablePath
            }, ignore_index=True)
            self.products = self.products.append({
                "soxspipe_recipe": self.recipeName,
                "product_label": "",
                "product_desc": "visualisation of goodness of order edge fitting",
                "file_name": plotName,
                "file_type": "PDF",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "file_path": plotPath
            }, ignore_index=True)
        if not isinstance(self.qc, bool):
            self.qc = self.qc.append({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "ORDER EDGE RES MEAN",
                "qc_value": mean_res,
                "qc_comment": "mean residual of order-edge fit",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }, ignore_index=True)
            self.qc = self.qc.append({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "ORDER EDGE RES STDEV",
                "qc_value": std_res,
                "qc_comment": "stdev of residuals to order edge fit",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }, ignore_index=True)

        self.log.debug('completed the ``get`` method')

        return self.products, self.qc

    def plot_results(
            self,
            orderPixelTable,
            orderPolyTable,
            orderMetaTable):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**
            - ``orderPixelTable`` -- the pixel table with residuals of fits
            - ``orderPolyTable`` -- data-frame of order-location polynomial coeff
            - ``orderMetaTable`` -- data-frame containing the limits of the fit

        **Return:**
            - ``filePath`` -- path to the plot pdf
        """
        self.log.debug('starting the ``plot_results`` method')

        allResiduals = np.concatenate((orderPixelTable[
            'xcoord_lower_fit_res'], orderPixelTable['xcoord_upper_fit_res']))
        allXcoords = np.concatenate((
            orderPixelTable['xcoord_lower'], orderPixelTable['xcoord_upper'])),
        allYcoords = np.concatenate((
            orderPixelTable['ycoord'], orderPixelTable['ycoord'])),

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
        std = np.nanstd(self.flatFrame.data)
        mean = np.nanmean(self.flatFrame.data)
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

        poly = chebyshev_order_xy_polynomials(
            log=self.log, yCol="y", orderCol="order", order_deg=self.orderDeg, y_deg=self.yDeg).poly

        # UPPER
        for index, row in orderPolyTable.iterrows():
            coeffupper = [float(v) for k, v in row.items() if "edgeup_" in k]
            coefflower = [float(v) for k, v in row.items() if "edgelow_" in k]

        uniqueOrders = orderPixelTable['order'].unique()
        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {"y": ylinelist}
        df = pd.DataFrame(myDict)

        from tabulate import tabulate

        for o in uniqueOrders:
            o = int(o)
            ymin = orderMetaTable.loc[(orderMetaTable["order"] == o), "ymin"].values[0]
            ymax = orderMetaTable.loc[(orderMetaTable["order"] == o), "ymax"].values[0]

            df["order"] = o
            xfitup = poly(df, *coeffupper)
            xfitlow = poly(df, *coefflower)
            # xfit = np.ones(len(xfit)) * \
            #     self.flatFrame.data.shape[1] - xfit
            xfitup, yfitup = zip(
                *[(x, y) for x, y in zip(xfitup, ylinelist) if x > 0 and x < (self.flatFrame.data.shape[1]) - 10 and y >= ymin and y <= ymax])
            xfitlow, yfitlow = zip(
                *[(x, y) for x, y in zip(xfitlow, ylinelist) if x > 0 and x < (self.flatFrame.data.shape[1]) - 10 and y >= ymin and y <= ymax])
            l = midrow.plot(yfitlow, xfitlow)
            u = midrow.plot(yfitup, xfitup, c=l[0].get_color())
            midrow.fill_between(yfitlow, xfitlow, xfitup, alpha=0.4, fc=l[0].get_color())
            midrow.text(yfitlow[10], xfitlow[10] + 5, int(o), fontsize=6, c='white', verticalalignment='bottom')

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

        filename = filenamer(
            log=self.log,
            frame=self.flatFrame,
            settings=self.settings
        )
        filename = filename.split("SLIT")[0] + "ORDER_EDGES_residuals.pdf"
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        filePath = f"{outDir}/{filename}"
        plt.savefig(filePath, dpi=720)

        self.log.debug('completed the ``plot_results`` method')
        return filePath

    def determine_order_flux_threshold(
            self,
            orderData,
            orderPixelTable):
        """*determine the flux threshold at the central column of each order*

        **Key Arguments:**
            - ``orderData`` -- one row in the orderTable
            - ``orderPixelTable`` the order table containing pixel arrays

        **Return:**
            - ``orderData`` -- orderData with min and max flux thresholds added
        """
        self.log.debug(
            'starting the ``determine_order_flux_threshold`` method')

        minThresholdPercenage = self.minThresholdPercenage
        maxThresholdPercenage = self.maxThresholdPercenage
        sliceWidth = self.sliceWidth
        sliceLength = self.sliceLength
        order = orderData["order"]

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable['order'] == order)
        xcoords = orderPixelTable.loc[mask, "xcoord_centre"].values
        ycoords = orderPixelTable.loc[mask, "ycoord"].values

        # xpd-update-filter-dataframe-column-values

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
        orderData["minThreshold"] = minvalue + \
            (maxvalue - minvalue) * minThresholdPercenage
        orderData["maxThreshold"] = minvalue + \
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
            plt.hlines(orderData["minThreshold"], 0, len(slice),
                       label=f'threshold {orderData["minThreshold"]:0.2f},  {orderData["maxThreshold"]:0.2f}', colors='red')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            plt.legend()
            plt.show()

        self.log.debug(
            'completed the ``determine_order_flux_threshold`` method')
        return orderData

    def determine_lower_upper_edge_pixel_positions(
            self,
            orderData):
        """*from a pixel postion somewhere on the trace of the order centre, return the lower and upper edges of the order*

        **Key Arguments:**
            - ``orderData`` -- one row in the orderTable

        **Return:**
            - ``orderData`` -- orderData with upper and lower edge xcoord arrays added
        """
        self.log.debug(
            'starting the ``determine_lower_upper_edge_limits`` method')

        sliceWidth = self.sliceWidth
        sliceLength = self.sliceLength
        halfSlice = sliceLength / 2
        minThreshold = orderData["minThreshold"]
        maxThreshold = orderData["maxThreshold"]
        x = orderData["xcoord_centre"]
        y = orderData["ycoord"]

        threshold = minThreshold

        # CUT A MEDIAN COLLAPSED SLICE
        slice = cut_image_slice(log=self.log, frame=self.flatFrame,
                                width=sliceWidth, length=sliceLength, x=x, y=y, median=True, plot=False)
        if slice is None:
            return orderData

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
            return orderData

        # REPORT THE EXACT PIXEL POSTION AT THE FLUX THRESHOLD
        xmax = xmaxguess - \
            (threshold - medSlide[xmaxguess]) / \
            (medSlide[xmaxguess - 1] - medSlide[xmaxguess]) - 2

        xmin = xminguess + \
            (threshold - medSlide[xminguess]) / \
            (medSlide[xminguess + 1] - medSlide[xminguess]) + 2

        # IF THE WIDTH BETWEEN MIN AND MAX IS TOO SMALL THEN REJECT
        if xmax - xmin < 10:
            return orderData
        else:
            orderData["xcoord_upper"] = xmax + int(x - halfSlice)
            orderData["xcoord_lower"] = xmin + int(x - halfSlice)

        # SANITY CHECK PLOT OF CROSS-SECTION
        if 1 == 0 and random.randint(1, 5001) < 2:
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
                       label='threshold', colors='red')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            plt.legend()
            plt.show()

        self.log.debug(
            'completed the ``determine_lower_upper_edge_limits`` method')
        return orderData
