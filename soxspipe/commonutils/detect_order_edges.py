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

from soxspipe.commonutils.filenamer import filenamer
import unicodecsv as csv
import collections
from soxspipe.commonutils.toolkit import unpack_order_table
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from os.path import expanduser
from soxspipe.commonutils.polynomials import chebyshev_order_xy_polynomials
import random

from soxspipe.commonutils import _base_detect
from soxspipe.commonutils.toolkit import cut_image_slice


from fundamentals import tools
from builtins import object
import random
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
        - ``tag`` -- e.g. '_DLAMP' to differentiate between UV-VIS lamps
        - ``sofName`` -- name of the originating SOF file
        - ``binx`` -- binning in x-axis
        - ``biny`` -- binning in y-axis

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
    productsTable, qcTable, orderDetectionCounts = edges.get()
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
            productsTable=False,
            tag="",
            sofName=False,
            binx=1,
            biny=1
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
        self.tag = tag
        self.sofName = sofName

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = flatFrame.header[kw("SEQ_ARM")]
        self.dateObs = flatFrame.header[kw("DATE_OBS")]

        self.binx = binx
        self.biny = biny

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        # DEG OF THE POLYNOMIALS TO FIT THE ORDER CENTRE LOCATIONS
        self.axisBDeg = self.recipeSettings["disp-axis-deg"]
        self.orderDeg = self.recipeSettings["order-deg"]

        self.inst = flatFrame.header[kw("INSTRUME")]
        # SET IMAGE ORIENTATION
        if self.inst == "SOXS":
            self.axisA = "y"
            self.axisB = "x"
            self.axisAbin = self.biny
            self.axisBbin = self.binx
        elif self.inst == "XSHOOTER":
            self.axisA = "x"
            self.axisB = "y"
            self.axisAbin = self.binx
            self.axisBbin = self.biny

        return None

    def get(self):
        """
        *get the detect_order_edges object*

        **Return:**
            - ``orderTablePath`` -- path to the new order table
        """
        self.log.debug('starting the ``get`` method')

        import numpy as np
        import pandas as pd

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
            log=self.log, orderTablePath=self.orderCentreTable, binx=self.binx, biny=self.biny, pixelDelta=9)

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
            orderPixelTable.loc[(orderPixelTable["order"] == o), ["minThreshold", "maxThreshold"]] = orderMetaTable.loc[(orderMetaTable["order"] == o), ["minThreshold", "maxThreshold"]].values

        print("\tMEASURING PIXEL-POSITIONS AT ORDER-EDGES WHERE FLUX THRESHOLDS ARE MET")
        orderPixelTable[f"{self.axisA}coord_upper"] = np.nan
        orderPixelTable[f"{self.axisA}coord_lower"] = np.nan
        orderPixelTable = orderPixelTable.apply(
            self.determine_lower_upper_edge_pixel_positions, axis=1)

        # DROP ROWS WITH NAN VALUES
        orderPixelTable.dropna(axis='index', how='any',
                               subset=[f"{self.axisA}coord_upper"], inplace=True)
        orderPixelTable.dropna(axis='index', how='any',
                               subset=[f"{self.axisA}coord_lower"], inplace=True)
        if self.axisAbin > 1:
            orderPixelTable[f"{self.axisA}coord_upper"] *= self.axisAbin
        if self.axisBbin > 1:
            orderPixelTable[f"{self.axisB}coord"] *= self.axisBbin

        # REDEFINE UNIQUE ORDERS IN CASE ONE OR MORE IS COMPLETELY MISSING
        uniqueOrders = orderPixelTable['order'].unique()

        for o in uniqueOrders:
            orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}min"] = np.nanmin(orderPixelTable.loc[(orderPixelTable["order"] == o), [f"{self.axisB}coord"]].values)
            orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}max"] = np.nanmax(orderPixelTable.loc[(orderPixelTable["order"] == o), [f"{self.axisB}coord"]].values)

        # REDEFINE UNIQUE ORDERS IN CASE ONE OR MORE IS COMPLETELY MISSING
        uniqueOrders = orderPixelTable['order'].unique()

        # SETUP EXPONENTS AHEAD OF TIME - SAVES TIME ON POLY FITTING
        for i in range(0, self.axisBDeg + 1):
            orderPixelTable[f"{self.axisB}_pow_{i}"] = orderPixelTable[f"{self.axisB}coord"].pow(i)
        for i in range(0, self.orderDeg + 1):
            orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        upperCoeff, orderPixelTable, clippedUpper = self.fit_global_polynomial(
            pixelList=orderPixelTable,
            axisBCol=f"{self.axisB}coord",
            axisACol=f"{self.axisA}coord_upper",
            orderCol="order",
            exponentsIncluded=True
        )

        # RENAME SOME INDIVIDUALLY
        orderPixelTable.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_upper_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_upper_fit_res"}, inplace=True)
        clippedUpper.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_upper_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_upper_fit_res"}, inplace=True)

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        print("\tFITTING POLYNOMIALS TO MEASURED PIXEL-POSITIONS AT LOWER ORDER-EDGES\n")
        lowerCoeff, orderPixelTable, clippedLower = self.fit_global_polynomial(
            pixelList=orderPixelTable,
            axisBCol=f"{self.axisB}coord",
            axisACol=f"{self.axisA}coord_lower",
            orderCol="order",
            exponentsIncluded=True
        )

        # RENAME SOME INDIVIDUALLY
        orderPixelTable.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_lower_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_lower_fit_res"}, inplace=True)
        clippedLower.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_lower_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_lower_fit_res"}, inplace=True)

        # orderLocations[o] = coeff
        coeff_dict = {
            "degorder_edgeup": self.orderDeg,
            f"deg{self.axisB}_edgeup": self.axisBDeg,
            "degorder_edgelow": self.orderDeg,
            f"deg{self.axisB}_edgelow": self.axisBDeg
        }
        n_coeff = 0
        for i in range(0, self.orderDeg + 1):
            for j in range(0, self.axisBDeg + 1):
                coeff_dict[f'edgeup_c{i}{j}'] = upperCoeff[n_coeff]
                coeff_dict[f'edgelow_c{i}{j}'] = lowerCoeff[n_coeff]
                n_coeff += 1
        coeffColumns = coeff_dict.keys()

        # RETURN BREAKDOWN OF ORDER EDGE DETECTION POSITION COUNTS (NEEDED TO COMPARE D and Q LAMPS)
        orderDetectionCounts = orderPixelTable['order'].value_counts(normalize=False).sort_index().to_frame()

        orderEdgePolyTable = pd.DataFrame([coeff_dict])

        # MERGE DATAFRAMES
        cols_to_use = orderEdgePolyTable.columns.difference(orderPolyTable.columns)
        orderPolyTable = orderPolyTable.join(
            orderEdgePolyTable[cols_to_use])

        # GENERATE AN OUTPUT PLOT OF RESULTS AND FITTING RESIDUALS
        print("\tMEASURING AND PLOTTING RESIDUALS OF FITS")
        allResiduals = np.concatenate((orderPixelTable[
            f"{self.axisA}coord_lower_fit_res"], orderPixelTable[f"{self.axisA}coord_upper_fit_res"]))
        plotPath = self.plot_results(
            orderPixelTable=orderPixelTable,
            orderPolyTable=orderPolyTable,
            orderMetaTable=orderMetaTable,
            clippedDataUpper=clippedUpper,
            clippedDataLower=clippedLower
        )

        # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE

        orderTablePath = self.write_order_table_to_file(
            frame=self.flatFrame, orderPolyTable=orderPolyTable, orderMetaTable=orderMetaTable)
        mean_res = np.mean(np.abs(allResiduals))
        min_res = np.min(np.abs(allResiduals))
        max_res = np.max(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        orderTablePath = os.path.abspath(orderTablePath)
        plotPath = os.path.abspath(plotPath)
        orderTableName = os.path.basename(orderTablePath)
        plotName = os.path.basename(plotPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # RECORD PRODUCTS AND QCs
        if not isinstance(self.products, bool):
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": f"ORDER_LOC{self.tag}",
                "product_desc": "table of coefficients from polynomial fits to order locations",
                "file_name": orderTableName,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "file_path": orderTablePath
            }).to_frame().T], ignore_index=True)
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": f"ORDER_LOC_RES{self.tag}",
                "product_desc": "visualisation of goodness of order edge fitting",
                "file_name": plotName,
                "file_type": "PDF",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "file_path": plotPath
            }).to_frame().T], ignore_index=True)
        if not isinstance(self.qc, bool):
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESMIN",
                "qc_value": min_res,
                "qc_comment": "[px] Minimum residual in order edge fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESMAX",
                "qc_value": max_res,
                "qc_comment": "[px] Maximum residual in order edge fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESRMS",
                "qc_value": std_res,
                "qc_comment": "[px] Std-dev of residual order edge fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``get`` method')

        return self.products, self.qc, orderDetectionCounts

    def plot_results(
            self,
            orderPixelTable,
            orderPolyTable,
            orderMetaTable,
            clippedDataUpper,
            clippedDataLower):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**
            - ``orderPixelTable`` -- the pixel table with residuals of fits
            - ``orderPolyTable`` -- data-frame of order-location polynomial coeff
            - ``orderMetaTable`` -- data-frame containing the limits of the fit
            - ``clippedDataUpper`` -- the sigma-clipped data from upper edge
            - ``clippedDataLower`` -- the sigma-clipped data from lower edge

        **Return:**
            - ``filePath`` -- path to the plot pdf
        """
        self.log.debug('starting the ``plot_results`` method')

        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt

        allResiduals = np.concatenate((orderPixelTable[
            f"{self.axisA}coord_lower_fit_res"], orderPixelTable[f"{self.axisA}coord_upper_fit_res"]))
        allAxisACoords = np.concatenate((
            orderPixelTable[f"{self.axisA}coord_lower"], orderPixelTable[f"{self.axisA}coord_upper"]))
        allAxisBCoords = np.concatenate((
            orderPixelTable[f"{self.axisB}coord"], orderPixelTable[f"{self.axisB}coord"]))
        allAxisACoordsClipped = np.concatenate((
            clippedDataLower[f"{self.axisA}coord_lower"], clippedDataUpper[f"{self.axisA}coord_upper"]))
        allAxisBCoordsClipped = np.concatenate((
            clippedDataLower[f"{self.axisB}coord"], clippedDataUpper[f"{self.axisB}coord"]))

        if self.axisBbin > 1:
            allAxisBCoords /= self.axisBbin
        if self.axisAbin > 1:
            allAxisACoords /= self.axisAbin

        arm = self.arm

        # a = plt.figure(figsize=(40, 15))
        if arm == "UVB" or self.inst == "SOXS":
            fig = plt.figure(figsize=(5, 13.5), constrained_layout=True)
            # CREATE THE GID OF AXES
            gs = fig.add_gridspec(6, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:, 0:2])
            bottomright = fig.add_subplot(gs[4:, 2:])
        else:
            fig = plt.figure(figsize=(6, 11), constrained_layout=True)
            # CREATE THE GID OF AXES
            gs = fig.add_gridspec(6, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:, 0:2])
            bottomright = fig.add_subplot(gs[4:, 2:])

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = self.flatFrame.data
        if self.axisA == "x":
            rotatedImg = np.rot90(rotatedImg, 1)
            rotatedImg = np.flipud(rotatedImg)
        # rotatedImg = self.flatFrame.data
        std = np.nanstd(self.flatFrame.data)
        mean = np.nanmean(self.flatFrame.data)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        toprow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=1)
        if self.axisA == "x":
            toprow.invert_yaxis()
        toprow.set_title(
            "upper and lower order edge detections", fontsize=10)
        toprow.scatter(allAxisBCoordsClipped, allAxisACoordsClipped, marker='x', c='red', s=4, alpha=0.6, linewidths=0.5)
        toprow.scatter(allAxisBCoords, allAxisACoords, marker='o', c='yellow', s=0.3, alpha=0.6)
        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])
        toprow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        toprow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        toprow.tick_params(axis='both', which='major', labelsize=9)

        midrow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=0.9)
        if self.axisA == "x":
            midrow.invert_yaxis()
        midrow.set_title(
            "order-location fit solutions", fontsize=10)
        if self.axisA == "x":
            axisALength = self.flatFrame.data.shape[1]
            axisBLength = self.flatFrame.data.shape[0]
        else:
            axisALength = self.flatFrame.data.shape[0]
            axisBLength = self.flatFrame.data.shape[1]

        if self.axisBbin > 1:
            axisBLength *= self.axisBbin
        if self.axisAbin > 1:
            axisALength *= self.axisAbin
        axisBlinelist = np.arange(0, axisBLength, 3)

        poly = chebyshev_order_xy_polynomials(
            log=self.log, axisBCol=self.axisB, orderCol="order", orderDeg=self.orderDeg, axisBDeg=self.axisBDeg).poly

        # UPPER
        for index, row in orderPolyTable.iterrows():
            coeffupper = [float(v) for k, v in row.items() if "edgeup_" in k]
            coefflower = [float(v) for k, v in row.items() if "edgelow_" in k]

        uniqueOrders = orderPixelTable['order'].unique()
        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {f"{self.axisB}": axisBlinelist}
        df = pd.DataFrame(myDict)

        colors = []
        for o in uniqueOrders:
            o = int(o)
            axisBmin = orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}min"].values[0]
            axisBmax = orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}max"].values[0]

            df["order"] = o
            axisAfitupStart = poly(df, *coeffupper)
            axisAfitlowStart = poly(df, *coefflower)

            # xfit = np.ones(len(xfit)) * \
            #     self.flatFrame.data.shape[1] - xfit
            axisAfitup, axisBfitup = zip(
                *[(a, b) for a, b in zip(axisAfitupStart, axisBlinelist) if a > 0 and a < (axisALength) - 10])
            axisAfitlow, axisBfitlow = zip(
                *[(a, b) for a, b in zip(axisAfitlowStart, axisBlinelist) if a > 0 and a < (axisALength) - 10])
            if len(axisBfitlow) < len(axisBfitup):
                half = int(len(axisAfitlowStart) / 2)
                try:
                    axisAfitlowExtra, axisBfitlowExtra = zip(
                        *[(0, b) for a, b in zip(axisAfitlowStart[:half], axisBlinelist[:half]) if (a < 0 or a > (axisALength) - 10) and b in axisBfitup])
                    axisAfitlow = axisAfitlowExtra + axisAfitlow
                    axisBfitlow = axisBfitlowExtra + axisBfitlow
                except:
                    pass
                try:
                    axisAfitlowExtra, axisBfitlowExtra = zip(
                        *[(0, b) for a, b in zip(axisAfitlowStart[half:], axisBlinelist[half:]) if (a < 0 or a > (axisALength) - 10) and b in axisBfitup])
                    axisAfitlow += axisAfitlowExtra
                    axisBfitlow += axisBfitlowExtra
                except:
                    pass

            if self.axisAbin > 1:
                axisAfitup = np.array(axisAfitup) / self.axisAbin
                axisAfitlow = np.array(axisAfitlow) / self.axisAbin
            if self.axisBbin > 1:
                axisBfitup = np.array(axisBfitup) / self.axisBbin
                axisBfitlow = np.array(axisBfitlow) / self.axisBbin

            l = midrow.plot(axisBfitlow, axisAfitlow)
            colors.append(l[0].get_color())
            u = midrow.plot(axisBfitup, axisAfitup, c=l[0].get_color())
            try:
                midrow.fill_between(axisBfitlow, axisAfitlow, axisAfitup, alpha=0.4, fc=l[0].get_color())
            except:
                pass
            midrow.text(axisBfitlow[10], axisAfitlow[10] + 5, int(o), fontsize=6, c="white", verticalalignment='bottom')

        # xfit = np.ones(len(xfit)) * \
        #     self.pinholeFrame.data.shape[1] - xfit
        # midrow.scatter(yfit, xfit, marker='x', c='blue', s=4)
        # midrow.set_yticklabels([])
        # midrow.set_xticklabels([])
        midrow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        midrow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        midrow.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        for o, c in zip(uniqueOrders, colors):
            mask = (orderPixelTable['order'] == o)
            orderAxisACoords = np.concatenate((orderPixelTable.loc[mask][f"{self.axisA}coord_lower"], orderPixelTable.loc[mask][f"{self.axisA}coord_upper"]))
            orderAxisBCoords = np.concatenate((orderPixelTable.loc[mask][f"{self.axisB}coord"], orderPixelTable.loc[mask][f"{self.axisB}coord"]))
            orderResiduals = np.concatenate((orderPixelTable.loc[mask][
                f"{self.axisA}coord_lower_fit_res"], orderPixelTable.loc[mask][f"{self.axisA}coord_upper_fit_res"]))
            bottomleft.scatter(orderAxisACoords, orderResiduals, alpha=0.6, s=0.2, c=c)
            bottomleft.text(orderAxisACoords[10], orderResiduals[10], int(o), fontsize=6, c=c, verticalalignment='bottom')
            bottomright.scatter(orderAxisBCoords, orderResiduals, alpha=0.6, s=0.2, c=c)
            bottomright.text(orderAxisBCoords[10], orderResiduals[10], int(o), fontsize=6, c=c, verticalalignment='bottom')
        bottomleft.set_xlabel(f'{self.axisA} pixel position')
        bottomleft.set_ylabel(f'{self.axisA} residual')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomright.set_xlabel(f'{self.axisB} pixel position')
        bottomright.tick_params(axis='both', which='major', labelsize=9)
        # bottomright.set_ylabel(f'{self.axisA} residual')
        bottomright.set_yticklabels([])

        mean_res = np.mean(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        fig.suptitle(f"detection of order-edge locations - flat-frame\n{subtitle}", fontsize=12)

        filename = filenamer(
            log=self.log,
            frame=self.flatFrame,
            settings=self.settings
        )
        filename = filename.split("SLIT")[0] + "ORDER_EDGES_residuals.pdf"
        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home) + "/qc/pdf"
        filePath = f"{outDir}/{filename}"
        plt.tight_layout()
        # plt.show()
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

        import numpy as np
        from scipy.signal import medfilt

        minThresholdPercenage = self.minThresholdPercenage
        maxThresholdPercenage = self.maxThresholdPercenage
        sliceWidth = self.sliceWidth
        sliceLength = self.sliceLength
        order = orderData["order"]

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable['order'] == order)
        axisAcoords = orderPixelTable.loc[mask, f"{self.axisA}coord_centre"].values
        axisBcoords = orderPixelTable.loc[mask, f"{self.axisB}coord"].values

        # xpd-update-filter-dataframe-column-values

        # DETERMINE THE FLUX THRESHOLD FROM THE CENTRAL COLUMN
        # CUT A MEDIAN COLLAPSED SLICE
        index = int(len(axisAcoords) / 2)
        if self.axisA == "x":
            x = axisAcoords[index]
            y = axisBcoords[index]
        else:
            y = axisAcoords[index]
            x = axisBcoords[index]

        slice, slice_length_offset, slice_width_centre = cut_image_slice(log=self.log, frame=self.flatFrame,
                                                                         width=sliceWidth, length=sliceLength, x=x, y=y, sliceAxis=self.axisA, median=True, plot=False)
        # SMOOTH WITH A MEDIAN FILTER
        medSlide = medfilt(slice, 9)
        # DETERMINE THRESHOLD FLUX VALUE
        maxvalue = np.max(
            medSlide[int(len(medSlide) / 2 - 8):int(len(medSlide) / 2 + 8)])
        minvalue = np.min(medSlide)

        orderData["minThreshold"] = minvalue + (maxvalue - minvalue) * minThresholdPercenage
        orderData["maxThreshold"] = minvalue + (maxvalue - minvalue) * maxThresholdPercenage

        # SANITY CHECK PLOT OF CROSS-SECTION
        if 1 == 0:
            import matplotlib.pyplot as plt
            # CHECK THE SLICE POINTS IF NEEDED
            print(order)
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

        import numpy as np
        from scipy.signal import medfilt

        sliceWidth = self.sliceWidth
        sliceLength = self.sliceLength
        halfSlice = sliceLength / 2
        minThreshold = orderData["minThreshold"]
        maxThreshold = orderData["maxThreshold"]

        if self.axisA == "x":
            x = orderData["xcoord_centre"]
            y = orderData["ycoord"]
            axisACoord = x
        else:
            x = orderData["xcoord"]
            y = orderData["ycoord_centre"]
            axisACoord = y

        threshold = minThreshold
        thresholdRange = maxThreshold - minThreshold

        # CUT A MEDIAN COLLAPSED SLICE
        slice, slice_length_offset, slice_width_centre = cut_image_slice(log=self.log, frame=self.flatFrame,
                                                                         width=sliceWidth, length=sliceLength, x=x, y=y, median=True, sliceAxis=self.axisA, plot=False)
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
                axisAmaxguess = np.where(secondHalf < threshold)[
                    0][0] + middle
                axisAminguess = np.where(firstHalf < threshold)[0][-1]
                hit = True
            except:
                threshold = threshold + thresholdRange * 0.1

        # IF WE STILL DIDN'T GET A HIT THEN REJECT
        if hit == False:
            return orderData

        # REPORT THE EXACT PIXEL POSTION AT THE FLUX THRESHOLD
        axisAmax = axisAmaxguess - \
            (threshold - medSlide[axisAmaxguess]) / \
            (medSlide[axisAmaxguess - 1] - medSlide[axisAmaxguess]) - 2

        axisAmin = axisAminguess + \
            (threshold - medSlide[axisAminguess]) / \
            (medSlide[axisAminguess + 1] - medSlide[axisAminguess]) + 2

        # IF THE WIDTH BETWEEN MIN AND MAX IS TOO SMALL THEN REJECT
        if axisAmax - axisAmin < 10:
            return orderData
        else:
            orderData[f"{self.axisA}coord_upper"] = axisAmax + int(axisACoord - halfSlice)
            orderData[f"{self.axisA}coord_lower"] = axisAmin + int(axisACoord - halfSlice)

        # SANITY CHECK PLOT OF CROSS-SECTION
        if 1 == 0 and random.randint(1, 5001) < 2:
            import matplotlib.pyplot as plt
            # CHECK THE SLICE POINTS IF NEEDED
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko', alpha=0.5)
            plt.plot(x, medSlide, 'rx', alpha=0.8)
            plt.plot(axisAmin, threshold, 'ro',
                     alpha=0.8, label="order edge")
            plt.plot(axisAmax, threshold, 'ro', alpha=0.8)
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
