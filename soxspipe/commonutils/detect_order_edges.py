#!/usr/bin/env python
# encoding: utf-8
"""
*using a fully-illuminated slit flat frame detect and record the order-edges*

Author
: David Young & Marco Landoni

Date Created
: September 18, 2020
"""
from datetime import datetime, date, time
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.toolkit import unpack_order_table, get_calibration_lamp
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from os.path import expanduser
from soxspipe.commonutils.polynomials import chebyshev_order_xy_polynomials
from soxspipe.commonutils import _base_detect
from soxspipe.commonutils.toolkit import cut_image_slice
from fundamentals import tools
from builtins import object
import sys
import os


os.environ['TERM'] = 'vt100'


class detect_order_edges(_base_detect):
    """
    *using a fully-illuminated slit flat frame detect and record the order-edges*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
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
        - ``extendToEdges`` -- if true, extend the order edge tracing to the edges of the frame (Default *True*)
        - ``lampTag`` -- add this tag to the end of the product filename (Default *False*)

    **Usage:**

    ```python
    from soxspipe.commonutils import detect_order_edges
    edges = detect_order_edges(
        log=log,
        flatFrame=flatFrame,
        orderCentreTable=orderCentreTable,
        settings=settings,
        recipeSettings=recipeSettings,
        recipeName="soxs-mflat",
        verbose=False,
        qcTable=False,
        productsTable=False,
        extendToEdges=True,
        lampTag=False
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
            recipeSettings=False,
            recipeName="soxs-mflat",
            verbose=False,
            qcTable=False,
            productsTable=False,
            tag="",
            sofName=False,
            binx=1,
            biny=1,
            extendToEdges=True,
            lampTag=False
    ):
        self.log = log
        log.debug("instantiating a new 'detect_order_edges' object")
        self.settings = settings
        self.recipeName = recipeName
        self.recipeSettings = recipeSettings

        self.orderCentreTable = orderCentreTable
        self.flatFrame = flatFrame
        self.verbose = verbose
        self.qc = qcTable
        self.products = productsTable
        self.tag = tag
        self.sofName = sofName
        self.extendToEdges = extendToEdges
        self.lampTag = lampTag

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = flatFrame.header[kw("SEQ_ARM")]
        self.dateObs = flatFrame.header[kw("DATE_OBS")]
        self.exptime = flatFrame.header[kw("EXPTIME")]
        try:
            self.slit = flatFrame.header[kw(f"SLIT_{self.arm}".upper())]
        except:
            self.log.warning(kw(f"SLIT_{self.arm}".upper()) + " keyword not found")
            self.slit = ""

        # if self.exptime < 59 or self.exptime > 61:
        #     raise Exception("too short")

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

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
            self.axisAbin = self.binx
            self.axisBbin = self.biny
        else:
            self.axisA = "y"
            self.axisB = "x"
            self.axisAbin = self.biny
            self.axisBbin = self.binx

        # self.lamp = get_calibration_lamp(log=log, frame=flatFrame, kw=self.kw)

        home = expanduser("~")
        self.qcDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}/"
        self.qcDir = self.qcDir.replace("//", "/")
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.qcDir):
            os.makedirs(self.qcDir)

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=self.flatFrame, show=False, ext='data', stdWindow=3, title=False, surfacePlot=True)

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
        from astropy.stats import mad_std

        self.log.print("\n# DETECTING THE ORDER EDGES FROM MASTER-FLAT FRAME")

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
            log=self.log, orderTablePath=self.orderCentreTable, binx=self.binx, biny=self.biny, pixelDelta=25)

        # REMOVE TOP 2 ORDERS IF BLOCKING FILTER USED
        if "JH" in self.slit:
            orderPixelTable = orderPixelTable.loc[(orderPixelTable["order"] > 12)]
            orderMetaTable = orderMetaTable.loc[(orderMetaTable["order"] > 12)]

        # ADD MIN AND MAX FLUX THRESHOLDS TO ORDER TABLE
        self.log.print("\tDETERMINING ORDER FLUX THRESHOLDS")
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

        self.log.print("\tMEASURING PIXEL-POSITIONS AT ORDER-EDGES WHERE FLUX THRESHOLDS ARE MET")
        orderPixelTable[f"{self.axisA}coord_upper"] = np.nan
        orderPixelTable[f"{self.axisA}coord_lower"] = np.nan
        orderPixelTable = orderPixelTable.apply(
            self.determine_lower_upper_edge_pixel_positions, axis=1)

        # DROP ROWS WITH NAN VALUES
        orderPixelTable.dropna(axis='index', how='any',
                               subset=[f"{self.axisA}coord_upper"], inplace=True)
        orderPixelTable.dropna(axis='index', how='any',
                               subset=[f"{self.axisA}coord_lower"], inplace=True)

        # MEASURE ORDER HEIGHT, FIND MEDIAN AND RE-ADJUST LOWER FROM UPPER EDGE
        orderPixelTable["order_height_upper"] = orderPixelTable[f"{self.axisA}coord_upper"] - orderPixelTable[f"{self.axisA}coord_centre"]
        orderPixelTable["order_height_lower"] = orderPixelTable[f"{self.axisA}coord_centre"] - orderPixelTable[f"{self.axisA}coord_lower"]

        # REDEFINE UNIQUE ORDERS IN CASE ONE OR MORE IS COMPLETELY MISSING
        uniqueOrders = orderPixelTable['order'].unique()

        # COUPLE THE UPPER AND LOWER EDGES FOR HOMOGENEOUS HEIGHT
        for o in uniqueOrders:
            mask = (orderPixelTable['order'] == o)
            medianHeight = orderPixelTable.loc[mask]['order_height_upper'].median() + mad_std(orderPixelTable.loc[mask]['order_height_upper'] * 3.)
            orderPixelTable.loc[mask, f"{self.axisA}coord_upper"] = orderPixelTable.loc[mask, f"{self.axisA}coord_centre"] + medianHeight
            medianHeight = orderPixelTable.loc[mask]['order_height_lower'].median() + mad_std(orderPixelTable.loc[mask]['order_height_lower'] * 3.)
            orderPixelTable.loc[mask, f"{self.axisA}coord_lower"] = orderPixelTable.loc[mask, f"{self.axisA}coord_centre"] - medianHeight

        if self.axisAbin > 1:
            orderPixelTable[f"{self.axisA}coord_upper"] *= self.axisAbin
        if self.axisAbin > 1:
            orderPixelTable[f"{self.axisA}coord_lower"] *= self.axisAbin
        if self.axisBbin > 1:
            orderPixelTable[f"{self.axisB}coord"] *= self.axisBbin

        if self.axisA == "x":
            minAxis = 1
            maxAxis = 0
        else:
            minAxis = 0
            maxAxis = 1

        for o in uniqueOrders:
            if self.extendToEdges:
                orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}min"] = 0
                orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}max"] = self.flatFrame.data.shape[maxAxis]
            else:
                orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}min"] = np.nanmin(orderPixelTable.loc[(orderPixelTable["order"] == o), [f"{self.axisB}coord"]].values)
                orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}max"] = np.nanmax(orderPixelTable.loc[(orderPixelTable["order"] == o), [f"{self.axisB}coord"]].values)

        # CONVERT COLUMN TYPE
        orderPixelTable[f"{self.axisB}coord"] = orderPixelTable[f"{self.axisB}coord"].astype(float)
        orderPixelTable[f"{self.axisA}coord_upper"] = orderPixelTable[f"{self.axisA}coord_upper"].astype(float)
        orderPixelTable[f"{self.axisA}coord_lower"] = orderPixelTable[f"{self.axisA}coord_lower"].astype(float)

        # SETUP EXPONENTS AHEAD OF TIME - SAVES TIME ON POLY FITTING
        for i in range(0, self.axisBDeg + 1):
            orderPixelTable[f"{self.axisB}_pow_{i}"] = orderPixelTable[f"{self.axisB}coord"].pow(i)
        for i in range(0, self.orderDeg + 1):
            orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        self.log.print("\tFITTING POLYNOMIALS TO MEASURED PIXEL-POSITIONS AT UPPER ORDER-EDGES\n")
        orderPixelTableUpper = orderPixelTable.dropna(axis='index', how='any',
                                                      subset=[f"{self.axisA}coord_upper"])

        upperCoeff, orderPixelTableUpper, clippedUpper = self.fit_global_polynomial(
            pixelList=orderPixelTableUpper,
            axisBCol=f"{self.axisB}coord",
            axisACol=f"{self.axisA}coord_upper",
            orderCol="order",
            exponentsIncluded=True
        )

        # RENAME COLUMNS
        orderPixelTableUpper.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_upper_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_upper_fit_res"}, inplace=True)
        clippedUpper.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_upper_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_upper_fit_res"}, inplace=True)

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        self.log.print("\tFITTING POLYNOMIALS TO MEASURED PIXEL-POSITIONS AT LOWER ORDER-EDGES\n")
        orderPixelTableLower = orderPixelTable.dropna(axis='index', how='any',
                                                      subset=[f"{self.axisA}coord_lower"])

        lowerCoeff, orderPixelTableLower, clippedLower = self.fit_global_polynomial(
            pixelList=orderPixelTable,
            axisBCol=f"{self.axisB}coord",
            axisACol=f"{self.axisA}coord_lower",
            orderCol="order",
            exponentsIncluded=True
        )

        # RENAME SOME INDIVIDUALLY
        orderPixelTableLower.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_lower_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_lower_fit_res"}, inplace=True)
        clippedLower.rename(columns={
            f"{self.axisA}_fit": f"{self.axisA}coord_lower_fit", f"{self.axisA}_fit_res": f"{self.axisA}coord_lower_fit_res"}, inplace=True)

        if isinstance(lowerCoeff, type(None)) or isinstance(upperCoeff, type(None)):
            raise Exception("Pipeline failed to determine the edges of the orders in this lamp-flat")

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

        # WRITE OUT THE FITS TO THE ORDER LOCATION TABLE
        orderTablePath = self.write_order_table_to_file(
            frame=self.flatFrame, orderPolyTable=orderPolyTable, orderMetaTable=orderMetaTable)

        allResiduals = np.concatenate((orderPixelTableLower[
            f"{self.axisA}coord_lower_fit_res"], orderPixelTableUpper[f"{self.axisA}coord_upper_fit_res"]))
        mean_res = np.mean(np.abs(allResiduals))
        min_res = np.min(np.abs(allResiduals))
        max_res = np.max(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # RECORD QCs
        if not isinstance(self.qc, bool):
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESMIN",
                "qc_value": f"{min_res:0.2E}",
                "qc_comment": "[px] Minimum residual in order edge fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESMAX",
                "qc_value": f"{max_res:0.2E}",
                "qc_comment": "[px] Maximum residual in order edge fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESRMS",
                "qc_value": f"{std_res:0.2E}",
                "qc_comment": "[px] Std-dev of residual order edge fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

        # GENERATE AN OUTPUT PLOT OF RESULTS AND FITTING RESIDUALS
        self.log.print("\tMEASURING AND PLOTTING RESIDUALS OF FITS")

        plotPath = self.plot_results(
            orderPixelTableUpper=orderPixelTableUpper,
            orderPixelTableLower=orderPixelTableLower,
            orderPolyTable=orderPolyTable,
            orderMetaTable=orderMetaTable,
            clippedDataUpper=clippedUpper,
            clippedDataLower=clippedLower
        )

        orderTablePath = os.path.abspath(orderTablePath)
        plotPath = os.path.abspath(plotPath)
        orderTableName = os.path.basename(orderTablePath)
        plotName = os.path.basename(plotPath)

        # RECORD PRODUCTS
        if not isinstance(self.products, bool):
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": f"ORDER_LOC{self.tag}",
                "product_desc": "table of coefficients from polynomial fits to order locations",
                "file_name": orderTableName,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "file_path": orderTablePath,
                "label": "PROD"
            }).to_frame().T], ignore_index=True)
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": f"ORDER_LOC_RES{self.tag}",
                "product_desc": "visualisation of goodness of order edge fitting",
                "file_name": plotName,
                "file_type": "PDF",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "file_path": plotPath,
                "label": "QC"
            }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``get`` method')

        return self.products, self.qc, orderDetectionCounts

    def plot_results(
            self,
            orderPixelTableUpper,
            orderPixelTableLower,
            orderPolyTable,
            orderMetaTable,
            clippedDataUpper,
            clippedDataLower):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**

        - ``orderPixelTableUpper`` -- the pixel table with residuals of fits for the upper edges
        - ``orderPixelTableLower`` -- the pixel table with residuals of fits for the lower edges
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

        allResiduals = np.concatenate((orderPixelTableLower[
            f"{self.axisA}coord_lower_fit_res"], orderPixelTableUpper[f"{self.axisA}coord_upper_fit_res"]))
        allAxisACoords = np.concatenate((
            orderPixelTableLower[f"{self.axisA}coord_lower"], orderPixelTableUpper[f"{self.axisA}coord_upper"]))
        allAxisBCoords = np.concatenate((
            orderPixelTableLower[f"{self.axisB}coord"], orderPixelTableUpper[f"{self.axisB}coord"]))
        allAxisACoordsClipped = np.concatenate((
            clippedDataLower[f"{self.axisA}coord_lower"], clippedDataUpper[f"{self.axisA}coord_upper"]))
        allAxisBCoordsClipped = np.concatenate((
            clippedDataLower[f"{self.axisB}coord"], clippedDataUpper[f"{self.axisB}coord"]))

        if self.axisBbin > 1:
            allAxisBCoords /= self.axisBbin
            allAxisBCoordsClipped /= self.axisBbin
        if self.axisAbin > 1:
            allAxisACoords /= self.axisAbin
            allAxisACoordsClipped /= self.axisAbin

        arm = self.arm

        rotateImage = self.detectorParams["rotate-qc-plot"]
        flipImage = self.detectorParams["flip-qc-plot"]

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = self.flatFrame.data
        if rotateImage:
            rotatedImg = np.rot90(rotatedImg, rotateImage / 90)
        if flipImage:
            rotatedImg = np.flipud(rotatedImg)
            if not rotateImage:
                aLen = rotatedImg.shape[0]
                aLen = rotatedImg.shape[0]
                allAxisACoords = aLen - allAxisACoords
                allAxisACoordsClipped = aLen - allAxisACoordsClipped

        if rotatedImg.shape[0] / rotatedImg.shape[1] > 0.8:
            fig = plt.figure(figsize=(6, 12))
            # CREATE THE GID OF AXES
            gs = fig.add_gridspec(6, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            if False:
                bottomleft = fig.add_subplot(gs[4:, 0:2])
                bottomright = fig.add_subplot(gs[4:, 2:])
            settingsAx = fig.add_subplot(gs[4:, 2:])
            qcAx = fig.add_subplot(gs[4:, 0:2])
        else:
            fig = plt.figure(figsize=(6, 10))
            # CREATE THE GID OF AXES
            gs = fig.add_gridspec(6, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            if False:
                bottomleft = fig.add_subplot(gs[4:, 0:2])
                bottomright = fig.add_subplot(gs[4:, 2:])
            settingsAx = fig.add_subplot(gs[4:, 2:])
            qcAx = fig.add_subplot(gs[4:, 0:2])

        # rotatedImg = self.flatFrame.data
        std = np.nanstd(self.flatFrame.data)
        mean = np.nanmean(self.flatFrame.data)
        vmax = mean + 2 * std
        vmin = mean - 1 * std

        toprow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=1)
        midrow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=0.9)

        toprow.set_title(
            "upper and lower order edge detections", fontsize=10)
        toprow.scatter(allAxisBCoords, allAxisACoords, marker='o', c='green', s=0.3, alpha=0.6, label="detected order edge location")
        toprow.scatter(allAxisBCoordsClipped, allAxisACoordsClipped, marker='x', c='red', s=4, alpha=0.6, linewidths=0.5, label="locations clipped during edge fitting")
        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])
        toprow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        toprow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        if arm.upper() == "UVB":
            toprow.xaxis.set_label_coords(0.2, -0.13)
        else:
            toprow.xaxis.set_label_coords(0.4, -0.13)
        toprow.tick_params(axis='both', which='major', labelsize=9)
        toprow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.13),
                      fontsize=4)

        toprow.set_xlim([0, rotatedImg.shape[1]])
        if self.axisA == "x":
            toprow.invert_yaxis()
            midrow.invert_yaxis()

        toprow.set_ylim([0, rotatedImg.shape[0]])
        midrow.set_ylim([0, rotatedImg.shape[0]])

        midrow.set_title(
            "order-location fit solutions", fontsize=10)
        if self.axisB == "y":
            axisALength = self.flatFrame.data.shape[1]
            axisBLength = self.flatFrame.data.shape[0]
        elif self.axisB == "x":
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

        uniqueOrders = orderPixelTableLower['order'].unique()
        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {f"{self.axisB}": axisBlinelist}
        df = pd.DataFrame(myDict)

        colors = []
        labelAdded = None
        for o in uniqueOrders:
            o = int(o)
            axisBmin = orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}min"].values[0]
            axisBmax = orderMetaTable.loc[(orderMetaTable["order"] == o), f"{self.axisB}max"].values[0]

            df["order"] = o
            axisAfitupStart = poly(df, *coeffupper)
            axisAfitlowStart = poly(df, *coefflower)

            if flipImage and not rotateImage:
                axisAfitupStart = aLen - axisAfitupStart
                axisAfitlowStart = aLen - axisAfitlowStart

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

            l = midrow.plot(axisBfitlow, axisAfitlow, linewidth=0.1)
            colors.append(l[0].get_color())

            if labelAdded == None:
                label1 = "order edge"
                label2 = "order region"
                labelAdded = True
            else:
                label1 = None
                label2 = None

            u = midrow.plot(axisBfitup, axisAfitup, c=l[0].get_color(), label=label1, linewidth=0.1)
            try:
                midrow.fill_between(axisBfitlow, axisAfitlow, axisAfitup, alpha=0.4, fc=l[0].get_color(), label=label2)
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
        midrow.xaxis.set_label_coords(0.5, -0.12)
        midrow.tick_params(axis='both', which='major', labelsize=9)

        midrow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.12),
                      fontsize=4)

        # PLOT THE FINAL RESULTS:
        if False:
            plt.subplots_adjust(top=0.92)
            for o, c in zip(uniqueOrders, colors):
                maskLower = (orderPixelTableLower['order'] == o)
                maskUpper = (orderPixelTableUpper['order'] == o)
                orderAxisACoords = np.concatenate((orderPixelTableLower.loc[maskLower][f"{self.axisA}coord_lower"], orderPixelTableUpper.loc[maskUpper][f"{self.axisA}coord_upper"]))
                orderAxisBCoords = np.concatenate((orderPixelTableLower.loc[maskLower][f"{self.axisB}coord"], orderPixelTableUpper.loc[maskUpper][f"{self.axisB}coord"]))
                orderResiduals = np.concatenate((orderPixelTableLower.loc[maskLower][
                    f"{self.axisA}coord_lower_fit_res"], orderPixelTableUpper.loc[maskUpper][f"{self.axisA}coord_upper_fit_res"]))
                bottomleft.scatter(orderAxisACoords, orderResiduals, alpha=0.6, s=0.2, c=c)
                try:
                    bottomleft.text(orderAxisACoords[10], orderResiduals[10], int(o), fontsize=6, c=c, verticalalignment='bottom')
                except:
                    pass
                bottomright.scatter(orderAxisBCoords, orderResiduals, alpha=0.6, s=0.2, c=c)
                try:
                    bottomright.text(orderAxisBCoords[10], orderResiduals[10], int(o), fontsize=6, c=c, verticalalignment='bottom')
                except:
                    pass

            bottomleft.set_xlabel(f'{self.axisA} pixel position')
            bottomleft.set_ylabel(f'{self.axisA} residual')
            bottomleft.tick_params(axis='both', which='major', labelsize=9)

            # PLOT THE FINAL RESULTS:
            plt.subplots_adjust(top=0.92)
            bottomright.set_xlabel(f'{self.axisB} pixel position')
            bottomright.tick_params(axis='both', which='major', labelsize=9)
            # bottomright.set_ylabel(f'{self.axisA} residual')
            bottomright.set_yticklabels([])

        from soxspipe.commonutils.toolkit import qc_settings_plot_tables
        qc_settings_plot_tables(log=self.log, qc=self.qc, qcAx=qcAx, settings={**self.recipeSettings, **{"exptime": self.exptime}}, settingsAx=settingsAx)

        mean_res = np.mean(np.abs(allResiduals))
        std_res = np.std(np.abs(allResiduals))

        lamp = ""
        if self.tag:
            lamp = " " + self.tag.replace("_", "")
        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        slitWidth = ""
        if self.slit:
            slitWidth = f" {self.slit.replace('x11','')}\""
        fig.suptitle(f"detection of order-edge locations - {arm}{lamp}{slitWidth} flat-frame\n{subtitle}", fontsize=12)

        if self.sofName:
            filename = self.sofName + f"_ORDER_LOCATIONS{self.tag}.pdf"
        else:
            filename = filenamer(
                log=self.log,
                frame=self.flatFrame,
                settings=self.settings
            )
            filename = filename.split("SLIT")[0] + "ORDER_EDGES_residuals.pdf"

        filePath = f"{self.qcDir}/{filename}"
        plt.tight_layout()
        # plt.show()
        plt.savefig(filePath, dpi=720)
        plt.close()

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

        if slice is None:
            return orderData

        # SMOOTH WITH A MEDIAN FILTER
        medSlide = medfilt(slice, 9)
        # DETERMINE THRESHOLD FLUX VALUE
        maxvalue = np.max(
            medSlide[int(len(medSlide) / 2 - 8):int(len(medSlide) / 2 + 8)])
        minvalue = np.min(medSlide)

        orderData["minThreshold"] = minvalue + (maxvalue - minvalue) * minThresholdPercenage
        orderData["maxThreshold"] = minvalue + (maxvalue - minvalue) * maxThresholdPercenage

        # SANITY CHECK PLOT OF CROSS-SECTION
        if False:
            import matplotlib.pyplot as plt
            # CHECK THE SLICE POINTS IF NEEDED
            self.log.print(order)
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko', alpha=0.5)
            plt.plot(x, medSlide, 'rx', alpha=0.8)
            plt.hlines(maxvalue, 0, len(slice), label='max')
            plt.hlines(minvalue, 0, len(slice), label='min')
            order = orderData["order"]
            plt.hlines(orderData["minThreshold"], 0, len(slice),
                       label=f'threshold {orderData["minThreshold"]:0.2f},  {orderData["maxThreshold"]:0.2f}', colors='red')
            plt.title(f'Order {order}, centre = {axisBcoords[index]}')
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

        minThresholdPercenage = self.minThresholdPercenage
        maxThresholdPercenage = self.maxThresholdPercenage

        import numpy as np
        from scipy.signal import medfilt
        import random
        sliceWidth = self.sliceWidth
        sliceLength = self.sliceLength
        halfSlice = sliceLength / 2

        orderData[f"{self.axisA}coord_upper"] = np.nan
        orderData[f"{self.axisA}coord_lower"] = np.nan

        if self.axisA == "x":
            x = orderData["xcoord_centre"]
            y = orderData["ycoord"]
            axisACoord = x
            axisBCoord = y
        else:
            x = orderData["xcoord"]
            y = orderData["ycoord_centre"]
            axisACoord = y
            axisBCoord = x

        # CUT A MEDIAN COLLAPSED SLICE
        slice, slice_length_offset, slice_width_centre = cut_image_slice(log=self.log, frame=self.flatFrame,
                                                                         width=sliceWidth, length=sliceLength, x=x, y=y, median=True, sliceAxis=self.axisA, plot=False)
        if slice is None:
            return orderData

        # SMOOTH WITH A MEDIAN FILTER
        medSlide = medfilt(slice, 9)

        # DETERMINE THRESHOLD FLUX VALUE
        maxvalue = np.max(
            medSlide[int(len(medSlide) / 2 - 8):int(len(medSlide) / 2 + 8)])
        minvalue = np.min(medSlide)

        orderData["minThreshold"] = minvalue + (maxvalue - minvalue) * minThresholdPercenage
        orderData["maxThreshold"] = minvalue + (maxvalue - minvalue) * maxThresholdPercenage
        minThreshold = orderData["minThreshold"]
        maxThreshold = orderData["maxThreshold"]
        threshold = minThreshold
        thresholdRange = maxThreshold - minThreshold

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

        # IF THE WIDTH BETWEEN MIN AND MAX IS TOO SMALL THEN REJECT
        if (medSlide[axisAminguess + 1] - medSlide[axisAminguess] == 0) or (medSlide[axisAmaxguess - 1] - medSlide[axisAmaxguess] == 0):
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
            # SANITY CHECK PLOT OF CROSS-SECTION
            if False:
                import matplotlib.pyplot as plt
                # CHECK THE SLICE POINTS IF NEEDED
                print(minThreshold, middle)
                print(axisAminguess, axisAmaxguess, threshold)
                print(axisAmin, axisAmax)
                print(len(slice))
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
                order = orderData["order"]
                plt.title(f'Order {order}, {self.axisB} = {axisBCoord}')
                plt.legend()
                plt.show()
            return orderData
        else:
            orderData[f"{self.axisA}coord_upper"] = axisAmax + int(axisACoord - halfSlice) + 1
            orderData[f"{self.axisA}coord_lower"] = axisAmin + int(axisACoord - halfSlice) - 1

        # SANITY CHECK PLOT OF CROSS-SECTION
        if False and random.randint(1, 5001) < 2000:
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
            order = orderData["order"]
            plt.title(f'Order {order}, {self.axisB} = {axisBCoord}')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            plt.legend()
            plt.show()

        self.log.debug(
            'completed the ``determine_lower_upper_edge_limits`` method')
        return orderData
