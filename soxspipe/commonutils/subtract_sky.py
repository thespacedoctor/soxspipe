#!/usr/bin/env python
# encoding: utf-8
"""
*Subtract the sky background using the Kelson Method*

Author
: David Young

Date Created
: April 14, 2022
"""

from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
from fundamentals import tools
from builtins import object
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils.toolkit import read_spectral_format
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
import sys
import os
from datetime import datetime
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.toolkit import quicklook_image
from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
from os.path import expanduser

os.environ["TERM"] = "vt100"


class subtract_sky(object):
    """
    *Subtract the sky background from a science image using the Kelson Method*

    A model of the sky-background is created using a method similar to that described in Kelson, D. (2003), *Optimal Techniques in Two-dimensional Spectroscopy: Background Subtraction for the 21st Century* (http://dx.doi.org/10.1086/375502). This model-background is then subtracted from the original science image to leave only non-sky flux.

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the soxspipe settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
    - ``objectFrame`` -- the image frame in need of sky subtraction
    - ``twoDMap`` -- 2D dispersion map image path
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products
    - ``dispMap`` -- path to dispersion map. Default *False*
    - ``sofName`` -- name of the originating SOF file. Default *False*
    - ``recipeName`` -- name of the recipe as it appears in the settings dictionary. Default *soxs-stare*
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""

    **Usage:**

    To setup your logger and settings, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/master/initialisation.html).

    To initiate a `subtract_sky` object, use the following:

    ```python
    from soxspipe.commonutils import subtract_sky
    skymodel = subtract_sky(
        log=log,
        settings=settings,
        recipeSettings=recipeSettings,
        objectFrame=objectFrame,
        twoDMap=twoDMap,
        qcTable=qc,
        productsTable=products,
        dispMap=dispMap
    )
    skymodelCCDData, skySubtractedCCDData, qcTable, productsTable = skymodel.subtract()
    ```

    """

    def __init__(
        self,
        log,
        settings,
        recipeSettings,
        objectFrame,
        twoDMap,
        qcTable,
        productsTable,
        dispMap=False,
        sofName=False,
        recipeName="soxs-stare",
        startNightDate="",
        debug=False,
    ):
        self.log = log
        log.debug("instantiating a new 'subtract_sky' object")
        self.settings = settings
        self.objectFrame = objectFrame
        self.twoDMap = twoDMap
        self.dispMap = dispMap
        self.qc = qcTable
        self.products = productsTable
        self.sofName = sofName
        self.recipeName = recipeName
        self.recipeSettings = recipeSettings
        self.startNightDate = startNightDate
        self.debug = debug

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(log=self.log, settings=self.settings).get
        kw = self.kw
        self.arm = objectFrame.header[kw("SEQ_ARM")]
        self.inst = objectFrame.header[kw("INSTRUME")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(log=log, settings=settings).get(self.arm)
        dp = self.detectorParams

        # UNPACK THE 2D DISP IMAGE MAP AND THE OBJECT IMAGE TO GIVE A
        # DATA FRAME CONTAINING ONE ROW FOR EACH PIXEL WITH COLUMNS X, Y, FLUX, WAVELENGTH, SLIT-POSITION, ORDER
        self.mapDF, self.interOrderMask = twoD_disp_map_image_to_dataframe(
            log=self.log,
            slit_length=dp["slit_length"],
            twoDMapPath=twoDMap,
            associatedFrame=self.objectFrame,
            kw=kw,
            dispAxis=self.detectorParams["dispersion-axis"],
        )

        # DETERMINE SLIT
        self.slit = objectFrame.header[kw(f"SLIT_{self.arm}".upper())]
        # ACCOUNT FOR BLOCKING FILTER
        if "JH" in self.slit:
            self.mapDF = self.mapDF.loc[(self.mapDF["order"] > 12)]

        quicklook_image(
            log=self.log,
            CCDObject=self.objectFrame,
            show=self.debug,
            ext=False,
            stdWindow=0.1,
            title="science frame awaiting sky-subtraction",
            surfacePlot=False,
            dispMap=dispMap,
            dispMapImage=twoDMap,
            settings=self.settings,
            skylines=True,
        )

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

        self.dateObs = objectFrame.header[kw("DATE_OBS")]

        # GET A TEMPLATE FILENAME USED TO NAME PRODUCTS
        if self.sofName:
            self.filenameTemplate = self.sofName + ".fits"
        else:
            self.filenameTemplate = filenamer(log=self.log, frame=self.objectFrame, settings=self.settings)

        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log, settings=settings, recipeName=recipeName, startNightDate=startNightDate
        )

        if self.arm != "NIR" and kw("WIN_BINX") in self.objectFrame.header:
            self.binx = int(self.objectFrame.header[kw("WIN_BINX")])
            self.biny = int(self.objectFrame.header[kw("WIN_BINY")])
        else:
            self.binx = 1
            self.biny = 1
        return

    def subtract(self):
        """
        *generate and subtract a sky-model from the input frame*

        **Return:**

        - ``skymodelCCDData`` -- CCDData object containing the model sky frame
        - ``skySubtractedCCDData`` -- CCDData object containing the sky-subtracted frame
        - ``qcTable`` -- the data frame containing measured QC metrics
        - ``productsTable`` -- the data frame containing collected output products
        """
        self.log.debug("starting the ``get`` method")

        import numpy as np
        import pandas as pd

        pd.options.mode.chained_assignment = None

        self.log.print(f"\n# MODELLING SKY BACKGROUND AND REMOVING FROM SCIENCE FRAME")

        # THESE PLACEHOLDERS ARE INITIALLY BLANK AND AWAITING PIXEL VALUES TO BE ADDED
        skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData = self.create_placeholder_images()

        uniqueOrders = self.mapDF["order"].unique()
        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # THE BSPLINE ORDER TO FIT SKY WITH
        self.bspline_order = self.recipeSettings["sky-subtraction"]["bspline_order"]

        # SELECT A SINGLE ORDER TO GENERATE QC PLOTS FOR
        self.qcPlotOrder = int(np.median(uniqueOrders)) - 1
        # uniqueOrders = [qcPlotOrder]

        allimageMapOrder = []
        allimageMapOrderWithObject = []

        # SPLIT ORDERS INTO THEIR OWN DATAFRAMES
        imageMapOrders = []
        for o in uniqueOrders:
            # SELECT DATAFRAME CONTAINING ONLY A SINGLE ORDER
            imageMapOrders.append(self.mapDF[self.mapDF["order"] == o])

        # GET OVER SAMPLED SKY & SKY+OBJECT AS LISTS OF DATAFRAMES
        self.log.print(f"\n  ## CLIPPING DEVIANT PIXELS AND PIXELS WITH OBJECT FLUX\n")

        # NOTE MULTIPROCESSING THIS BLOCK RESULTS IN SLOWER PERFORMANCE
        for o in uniqueOrders:
            # SELECT ONLY A DATAFRAME CONTAINING ONLY A SINGLE ORDER
            imageMapOrder = self.mapDF[self.mapDF["order"] == o]

            # MASK OUTLYING PIXELS (imageMapOrderWithObject) AND ALSO THEN THE OBJECT PIXELS (imageMapOrderSkyOnly)
            imageMapOrder = self.get_over_sampled_sky_from_order(
                imageMapOrder,
                clipBPs=True,
                clipSlitEdge=self.recipeSettings["sky-subtraction"]["clip-slit-edge-fraction"],
            )
            allimageMapOrder.append(imageMapOrder)

        # MASK OUT OBJECT PIXELS
        allimageMapOrder = self.clip_object_slit_positions(
            allimageMapOrder, aggressive=self.recipeSettings["sky-subtraction"]["aggressive_object_masking"]
        )

        self.log.print(
            f"\n  ## FITTING SKY-FLUX WITH A BSPLINE (WAVELENGTH) AND LOW-ORDER POLY (SLIT-ILLUMINATION PROFILE)\n"
        )

        # NOTE MULTIPROCESSING THIS BLOCK RESULTS IN SLOWER PERFORMANCE
        newAllimageMapOrder = []
        allFluxErrorRatios = []
        allResidualFloor = []
        alltck = []
        allKnots = []
        totalKnots = 0
        for o, imageMapOrder in zip(uniqueOrders, allimageMapOrder):
            imageMapOrder, tck, knots, flux_error_ratio, residualFloor = self.fit_bspline_curve_to_sky(imageMapOrder)
            totalKnots += len(knots)
            newAllimageMapOrder.append(imageMapOrder)
            allFluxErrorRatios.append(flux_error_ratio)
            allResidualFloor.append(residualFloor)
            alltck.append(tck)
            allKnots.append(knots)
        allimageMapOrder = newAllimageMapOrder

        allFluxErrorRatios = np.concatenate(allFluxErrorRatios)
        self.log.print(
            f"\n\tFULL FRAME SKY-MODEL FLUX TO ERROR METRICS: MEAN {allFluxErrorRatios.mean():0.3f}, STD {allFluxErrorRatios.std():0.3f}, MEDIAN {np.median(allFluxErrorRatios):0.3f}, MEAN RES FLOOR: {np.mean(allResidualFloor):0.3f}, KNOT COUNT: {totalKnots}"
        )
        # self.log.print(f'\t{allFluxErrorRatios.mean():0.3f},  {allFluxErrorRatios.std():0.3f},  {np.median(allFluxErrorRatios):0.3f},  {allFluxErrorRatios.max():0.3f},  {allFluxErrorRatios.min():0.3f}, {allFluxErrorRatios.max()-allFluxErrorRatios.min():0.3f},{np.mean(allResidualFloor):0.3f},{totalKnots}')

        for o, imageMapOrder, tck, knots in zip(uniqueOrders, allimageMapOrder, alltck, allKnots):
            if isinstance(imageMapOrder, pd.core.frame.DataFrame):
                # INJECT THE PIXEL VALUES BACK INTO THE PLACEHOLDER IMAGES
                skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData = (
                    self.add_data_to_placeholder_images(
                        imageMapOrder, skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData
                    )
                )
                if self.debug or (
                    o == self.qcPlotOrder and self.recipeSettings["sky-subtraction"]["sky_model_qc_plot"]
                ):
                    qc_plot_path = self.plot_sky_sampling(
                        order=o, imageMapOrderDF=imageMapOrder, tck=tck, knotLocations=knots
                    )
                    basename = os.path.basename(qc_plot_path)
                    self.products = pd.concat(
                        [
                            self.products,
                            pd.Series(
                                {
                                    "soxspipe_recipe": "soxs-stare",
                                    "product_label": "SKY_MODEL_QC_PLOTS",
                                    "file_name": basename,
                                    "file_type": "PDF",
                                    "obs_date_utc": self.dateObs,
                                    "reduction_date_utc": utcnow,
                                    "product_desc": f"QC plots for the sky-background modelling",
                                    "file_path": qc_plot_path,
                                    "label": "QC",
                                }
                            )
                            .to_frame()
                            .T,
                        ],
                        ignore_index=True,
                    )

        # SET NANs TO 0
        skymodelCCDData.data[np.isnan(skymodelCCDData.data)] = 0
        skymodelCCDData.uncertainty.array[np.isnan(skymodelCCDData.uncertainty.array)] = 0
        skySubtractedCCDData.data[np.isnan(skySubtractedCCDData.data)] = 0
        skySubtractedCCDData.uncertainty.array[np.isnan(skySubtractedCCDData.uncertainty.array)] = 0

        if self.recipeSettings["sky-subtraction"]["sky_model_qc_plot"]:
            comparisonPdf = self.plot_image_comparison(self.objectFrame, skymodelCCDData, skySubtractedCCDData)

            filename = os.path.basename(comparisonPdf)
            self.products = pd.concat(
                [
                    self.products,
                    pd.Series(
                        {
                            "soxspipe_recipe": "soxs-stare",
                            "product_label": "SKY SUBTRACTION QUICKLOOK",
                            "file_name": filename,
                            "file_type": "PDF",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "product_desc": f"Sky-subtraction quicklook",
                            "file_path": comparisonPdf,
                            "label": "QC",
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )

        self.log.debug("completed the ``get`` method")
        return skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData, self.qc, self.products

    def get_over_sampled_sky_from_order(self, imageMapOrder, clipBPs=True, clipSlitEdge=False):
        """*unpack the over sampled sky from an order*

        **Key Arguments:**

        - ``imageMapOrder`` -- single order dataframe from object image and 2D map
        - ``clipBPs`` -- clip bad-pixels? Default *True*
        - ``clipSlitEdge`` -- clip the slit edges. Percentage of slit width to clip. Default *False*

        **Return:**

        - `imageMapOrder` -- input order dataframe with outlying pixels masked AND object pixels masked

        **Usage:**

        ```python
        imageMapOrderWithObject, imageMapOrderSkyOnly = skymodel.get_over_sampled_sky_from_order(imageMapOrder, o, ignoreBP=False, clipSlitEdge=0.00)
        ```
        """
        self.log.debug("starting the ``get_over_sampled_sky_from_order`` method")

        from astropy.stats import sigma_clip, mad_std

        # COLLECT SETTINGS
        percentile_clipping_sigma = self.recipeSettings["sky-subtraction"]["percentile_clipping_sigma"]
        percentile_clipping_iterations = self.recipeSettings["sky-subtraction"]["percentile_clipping_iterations"]
        percentile_rolling_window_size = self.recipeSettings["sky-subtraction"]["percentile_rolling_window_size"]
        self.rollingWindowSize = int(percentile_rolling_window_size)

        imageMapOrder["flagged_all_clipped"] = False
        imageMapOrder["flagged_object_clipped"] = False
        imageMapOrder["flagged_bspline_clipped"] = False
        imageMapOrder["flagged_edge_clipped"] = False
        imageMapOrder["flagged_bad_pixel_clipped"] = False

        # ASSIGN ORDER-EDGE PIXELS A 'clipped' FLAG
        if clipSlitEdge:
            slitRange = imageMapOrder["slit_position"].max() - imageMapOrder["slit_position"].min()
            clipSlitEdge *= slitRange
            mask = (imageMapOrder["slit_position"] > imageMapOrder["slit_position"].max() - clipSlitEdge) | (
                imageMapOrder["slit_position"] < imageMapOrder["slit_position"].min() + clipSlitEdge
            )
            imageMapOrder.loc[mask, "flagged_edge_clipped"] = True
            imageMapOrder.loc[mask, "flagged_all_clipped"] = True

        # ASSIGN BAD-PIXELS A 'clipped' FLAG?
        if clipBPs:
            mask = imageMapOrder["mask"] == True
            imageMapOrder.loc[mask, "flagged_bad_pixel_clipped"] = True
            imageMapOrder.loc[mask, "flagged_all_clipped"] = True

        # CLIP THE MOST DEVIANT PIXELS WITHIN A WAVELENGTH ROLLING WINDOW - BAD-PIXELS AND CRHs AND OBJECTS
        imageMapOrder = self.rolling_window_clipping(
            imageMapOrderDF=imageMapOrder,
            windowSize=self.rollingWindowSize,
            sigma_clip_limit=percentile_clipping_sigma,
            max_iterations=percentile_clipping_iterations,
        )

        self.log.debug("completed the ``get_over_sampled_sky_from_order`` method")
        return imageMapOrder

    def plot_sky_sampling(self, order, imageMapOrderDF, tck=False, knotLocations=False):
        """*generate a plot of sky sampling*

        **Key Arguments:**

        - ``order`` -- the order number.
        - ``imageMapOrderDF`` -- dataframe with various processed data for order
        - ``tck`` -- spline parameters to replot
        - ``knotLocations`` -- wavelength locations of all knots used in the fit

        **Return:**

        - ``filePath`` -- path to the generated QC plots PDF

        **Usage:**

        ```python
        self.plot_sky_sampling(
            order=myOrder,
            imageMapOrderDF=imageMapOrder
        )
        ```
        """
        self.log.debug("starting the ``plot_sky_sampling`` method")

        import numpy as np
        import scipy.interpolate as ip
        import numpy.ma as ma
        from copy import copy

        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib import colors

        # SET COLOURS FOR VARIOUS STAGES
        red = "#dc322f"
        blue = "#268bd2"
        black = "#002b36"
        grey = "#93a1a1"
        green = "green"
        orange = "#cb4b16"
        violet = "#6c71c4"
        purple = "purple"

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotateImage = self.detectorParams["rotate-qc-plot"]
        flipImage = self.detectorParams["flip-qc-plot"]

        # MAKE A COPY OF THE FRAME TO NOT ALTER ORIGINAL DATA
        frame = self.objectFrame.copy()

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=frame,
            show=False,
            ext="data",
            stdWindow=3,
            title=False,
            surfacePlot=True,
            saveToPath=False,
        )

        # SETUP THE PLOT SUB-PANELS
        # print("FIX ME")
        if True:
            fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=180)
        else:
            # REMOVE ME
            fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=180)

        gs = fig.add_gridspec(11, 4)
        # CREATE THE GID OF AXES
        onerow = fig.add_subplot(gs[1:2, :])
        tworow = fig.add_subplot(gs[2:4, :])
        threerow = fig.add_subplot(gs[4:5:, :])
        fourrow = fig.add_subplot(gs[5:6:, :])
        fiverow = fig.add_subplot(gs[6:7:, :])
        sixrow = fig.add_subplot(gs[7:8:, :])
        sevenrow = fig.add_subplot(gs[8:9:, :])
        eightrow = fig.add_subplot(gs[9:10:, :])
        ninerow = fig.add_subplot(gs[10:11:, :])

        # FIND ORDER PIXELS - MASK THE REST
        nonOrderMask = np.ones_like(frame.data)
        for x, y in zip(imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB]):
            if self.detectorParams["dispersion-axis"] == "x":
                nonOrderMask[y][x] = 0
            else:
                nonOrderMask[x][y] = 0

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        nonOrderMask = ma.make_mask(nonOrderMask)
        combinedMask = (nonOrderMask == 1) | (frame.mask == 1)
        frame.mask = nonOrderMask == 1

        # RAW IMAGE PANEL
        # ROTATE THE IMAGE FOR BETTER LAYOUT
        if rotateImage:
            rotatedImg = np.flipud(np.rot90(frame, 1))
        else:
            rotatedImg = frame
        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        maskedDataArray = np.ma.array(frame.data, mask=combinedMask)
        maskedDataValues = np.array(maskedDataArray.filled(np.nan), dtype=float, copy=True)
        std = np.nanstd(maskedDataValues)
        mean = np.nanmean(maskedDataValues)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        im = onerow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap="gray", alpha=1)
        medianValue = np.median(rotatedImg.data.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="unprocessed frame")]

        onerow.set_title("Object & Sky Frame", fontsize=10)
        onerow.set_xlabel("y-axis", fontsize=10)
        onerow.set_ylabel("x-axis", fontsize=10)
        ylimMinImage = imageMapOrderDF[self.axisB].min() - 10
        ylimMaxImage = imageMapOrderDF[self.axisB].max() + 10
        onerow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        onerow.set_xlim(ylimMinImage, ylimMaxImage)

        # ORIGINAL DATA AND PERCENTILE SMOOTHED WAVELENGTH VS FLUX
        tworow.plot(
            imageMapOrderDF["wavelength"].values,
            imageMapOrderDF["flux"].values,
            label="unprocessed",
            alpha=0.2,
            c=grey,
            zorder=0,
        )
        tworow.set_title(
            "STEP 1. Identify and clip outlying pixels (CRHs etc) and pixels containing object flux.", fontsize=10
        )

        # RAW MARKERS
        tworow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "flux"].values,
            label="unclipped",
            s=0.5,
            c=black,
            alpha=1.0,
            zorder=1,
        )

        columnName = [
            "flagged_object_clipped",
            "flagged_bad_pixel_clipped",
            "flagged_edge_clipped",
            "flagged_bspline_clipped",
        ]
        colours = [blue, orange, green, purple]
        alphas = [0.2, 0.3, 0.8, 1.0]
        labels = ["flagged_object_clipped", "bad pixels", "order edges", "bspline clipped"]

        for cn, cl, lb, al in zip(columnName, colours, labels, alphas):
            tworow.scatter(
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "wavelength"].values,
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "flux"].values,
                label=lb,
                s=8,
                marker="x",
                c=cl,
                zorder=3,
                alpha=al,
            )

        # PERCENTILE LINE
        tworow.plot(
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False)
                & (imageMapOrderDF["flagged_object_clipped"] == False),
                "wavelength",
            ].values,
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False)
                & (imageMapOrderDF["flagged_object_clipped"] == False),
                "flux_percentile_smoothed",
            ].values,
            label="percentile-smoothed",
            c=blue,
            zorder=3,
        )

        # SIGMA RESIDUAL
        weights = tworow.plot(
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "residual_windowed_std"].values
            - imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "residual_windowed_std"].max() * 1.2,
            label="$\\sigma$ residual scatter (shifted)",
            c=black,
        )
        ylimmin = (
            -imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "residual_windowed_std"].max() * 1.3
        )

        if ylimmin < -3000:
            ylimmin = -300

        from astropy.stats import sigma_clip, mad_std

        # SIGMA-CLIP THE DATA
        masked = sigma_clip(
            imageMapOrderDF["flux"], sigma_lower=30, sigma_upper=30, maxiters=1, cenfunc="median", stdfunc=mad_std
        )

        tworow.set_ylim(ylimmin, masked.max())

        tworow.set_ylabel("flux ($e^{-}$)", fontsize=10)
        tworow.set_xlabel("wavelength", fontsize=10)
        tworow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.0)
        # tworow.set_xticks([], [])

        # SLIT-POSITION RESIDUAL PANEL (SHOWING OBJECT)
        std = imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == True, "residual_global_sigma"].std()
        median = imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == True, "residual_global_sigma"].median()

        threerow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "slit_position"].values,
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "residual_global_sigma"].values,
            label="deviations",
            s=0.5,
            alpha=0.5,
            c=grey,
            zorder=1,
        )

        columnName = ["flagged_object_clipped", "flagged_bad_pixel_clipped", "flagged_bspline_clipped"]
        colours = [blue, orange, purple]
        alphas = [0.2, 1, 1]
        labels = ["flagged_object_clipped", "bad pixels", "bspline clipped"]

        for cn, cl, lb, al in zip(columnName, colours, labels, alphas):
            threerow.scatter(
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "slit_position"].values,
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "residual_global_sigma"].values,
                label=lb,
                s=3,
                marker="x",
                c=cl,
                zorder=3,
                alpha=al,
            )

        threerow.set_ylim(median - 3 * std, median + 7 * std)
        threerow.set_xlabel("slit-position relative to slit centre (arcsec)", fontsize=10)
        threerow.set_ylabel("flux minus smoothed flux residual ($\\sigma$)", fontsize=10)

        threerow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.0)

        # IMAGE SHOWING CLIPPED PIXEL MASK
        im = fourrow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap="gray", alpha=1)

        columnName = [
            "flagged_object_clipped",
            "flagged_bad_pixel_clipped",
            "flagged_edge_clipped",
            "flagged_bspline_clipped",
        ]
        colours = [blue, orange, green, purple]
        alphas = [1, 1, 1, 1]
        labels = ["flagged_object_clipped", "bad pixels", "order edges", "bspline clipped"]
        patches = []

        for cn, cl, lb, al in zip(columnName, colours, labels, alphas):
            clippedMask = nonOrderMask
            clippedMask = np.zeros_like(frame.data)
            for x, y in zip(
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, self.axisA].values,
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, self.axisB].values,
            ):
                clippedMask[y][x] = 1
            clippedMask = ma.make_mask(clippedMask)
            imageMask = np.ma.array(np.ones_like(frame.data), mask=~clippedMask)
            # MAKE A COLOR MAP OF FIXED COLORS
            cmap = colors.ListedColormap([cl, cl])
            bounds = [0, 5, 10]
            norm = colors.BoundaryNorm(bounds, cmap.N)
            cmap.set_bad(cl, 0.0)
            fourrow.imshow(np.flipud(np.rot90(imageMask, 1)), cmap=cmap, norm=norm, alpha=al, interpolation="nearest")
            patches.append(mpatches.Patch(color=cl, label=lb))

        fourrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

        nonOrderMask = nonOrderMask == 0
        imageMask = np.ma.array(np.ones_like(frame.data), mask=nonOrderMask)
        cmap = copy(cm.gray)
        cmap.set_bad("green", 0.0)
        fourrow.imshow(np.flipud(np.rot90(imageMask, 1)), vmin=-10, vmax=-9, cmap=cmap, alpha=1.0)
        fourrow.set_xlabel("y-axis", fontsize=10)
        fourrow.set_ylabel("x-axis", fontsize=10)
        fourrow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        fourrow.set_xlim(ylimMinImage, ylimMaxImage)
        # fourrow.invert_xaxis()

        # PLOT WAVELENGTH VS FLUX SKY MODEL
        fiverow.set_title("STEP 2. Fit a univariate bspline to sky-flux as a function of wavelength", fontsize=10)
        fiverow.scatter(
            imageMapOrderDF.loc[(imageMapOrderDF["flagged_all_clipped"] == False), "wavelength"].values,
            imageMapOrderDF.loc[(imageMapOrderDF["flagged_all_clipped"] == False), "flux"].values,
            s=1,
            c=black,
            alpha=0.1,
            zorder=1,
            label="unclipped pixels",
        )

        if False:
            fiverow.scatter(
                imageMapOrderDF.loc[imageMapOrderDF["flagged_bspline_clipped"] == True, "wavelength"].values,
                imageMapOrderDF.loc[imageMapOrderDF["flagged_bspline_clipped"] == True, "flux"].values,
                s=8,
                c=violet,
                marker="x",
                alpha=1.0,
                zorder=1,
                label="flagged_all_clipped",
            )

        if tck:
            wl = np.linspace(imageMapOrderDF["wavelength"].min(), imageMapOrderDF["wavelength"].max(), 1000000)
            sky = ip.splev(wl, tck)
            # sky = ip.splev(wl, tck, der=1)
            # sky = ip.splev(wl, tck, der=2)
            if not isinstance(knotLocations, bool) and len(knotLocations):
                knotSky = ip.splev(knotLocations, tck)
                fiverow.scatter(knotLocations, knotSky, marker=7, s=5, alpha=0.2, c=red, zorder=3, label="knots")
            skymodel = fiverow.plot(wl, sky, label="sky model", c=blue, zorder=3)
        else:
            skymodel = fiverow.plot(
                imageMapOrderDF["wavelength"].values,
                imageMapOrderDF["sky_model"].values,
                label="sky model",
                c=blue,
                zorder=3,
            )
        if ylimmin < -3000:
            ylimmin = -300
        fiverow.set_ylim(
            imageMapOrderDF["flux_percentile_smoothed"].min() - 30,
            imageMapOrderDF["flux_percentile_smoothed"].max() * 1.2,
        )
        fiverow.set_ylabel("counts", fontsize=10)
        fiverow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.0)

        # BUILD IMAGE OF SKY MODEL
        skyModelImage = np.zeros_like(frame.data)
        for x, y, skypixel in zip(
            imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_model"]
        ):
            skyModelImage[y][x] = skypixel
        nonOrderMask = nonOrderMask == 0
        skyModelImage = np.ma.array(skyModelImage, mask=nonOrderMask)
        cmap = copy(cm.gray)
        std = np.nanstd(skyModelImage)
        mean = np.nanmean(skyModelImage)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        im = sixrow.imshow(np.flipud(np.rot90(skyModelImage, 1)), vmin=vmin, vmax=vmax, cmap=cmap, alpha=1.0)
        sixrow.set_ylabel("x-axis", fontsize=10)
        sixrow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        sixrow.set_xlim(ylimMinImage, ylimMaxImage)
        # sixrow.invert_xaxis()
        medianValue = np.median(skyModelImage.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="sky model")]
        sixrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
        sixrow.set_xticks([], [])

        # BUILD SKY-SUBTRACTED IMAGE
        skySubImage = np.zeros_like(frame.data)
        for x, y, skypixel in zip(
            imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_subtracted_flux"]
        ):
            skySubImage[y][x] = skypixel
        skySubMask = nonOrderMask == 1
        skySubImage = np.ma.array(skySubImage, mask=skySubMask)
        cmap = copy(cm.gray)
        std = np.nanstd(skySubImage)
        mean = np.nanmedian(skySubImage)
        vmax = mean + 0.2 * std
        vmin = mean - 0.2 * std
        im = sevenrow.imshow(np.flipud(np.rot90(skySubImage, 1)), vmin=0, vmax=50, cmap=cmap, alpha=1.0)
        sevenrow.set_title("STEP 3. Subtract the sky-model from the original data.", fontsize=10)
        sevenrow.set_xlabel("y-axis", fontsize=10)
        sevenrow.set_ylabel("x-axis", fontsize=10)
        sevenrow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        sevenrow.set_xlim(ylimMinImage, ylimMaxImage)
        # sevenrow.invert_xaxis()
        medianValue = np.median(skySubImage.data.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="sky-subtracted frame")]
        sevenrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

        # SUBTRACTED SKY RESIDUAL PANEL
        eightrow.scatter(
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "wavelength"
            ].values,
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] > 0),
                "sky_subtracted_flux",
            ].values,
            s=3,
            alpha=0.2,
            c="orange",
            zorder=1,
            label="slit position > 0",
        )
        eightrow.scatter(
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "wavelength"
            ].values,
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] < 0),
                "sky_subtracted_flux",
            ].values,
            s=3,
            alpha=0.2,
            c=blue,
            zorder=1,
            label="slit position < 0",
        )
        eightrow.scatter(
            knotLocations, np.zeros_like(knotLocations), marker=7, s=15, alpha=0.7, c=red, zorder=3, label="knots"
        )
        eightrow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.0)

        mean = np.absolute(
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "sky_subtracted_flux"]
        ).mean()
        std = np.absolute(
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "sky_subtracted_flux"]
        ).std()

        eightrow.set_ylim(-1000, 1000)
        eightrow.set_xlabel("wavelength (nm)", fontsize=10)
        eightrow.set_ylabel("residual", fontsize=10)

        # SUBTRACTED SKY RESIDUAL/ERROR PANEL
        ninerow.scatter(
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "wavelength"
            ].values,
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] > 0),
                "sky_subtracted_flux_weighted",
            ].values,
            s=3,
            alpha=0.2,
            c="orange",
            zorder=1,
        )
        ninerow.scatter(
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "wavelength"
            ].values,
            imageMapOrderDF.loc[
                (imageMapOrderDF["flagged_all_clipped"] == False) & (imageMapOrderDF["slit_position"] < 0),
                "sky_subtracted_flux_weighted",
            ].values,
            s=3,
            alpha=0.2,
            c=blue,
            zorder=1,
        )
        ninerow.scatter(knotLocations, np.zeros_like(knotLocations), marker=7, s=15, alpha=0.7, c=red, zorder=3)

        mean = np.absolute(
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "sky_subtracted_flux_weighted"]
        )[100:-100].mean()
        std = np.absolute(
            imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "sky_subtracted_flux_weighted"]
        )[100:-100].std()

        try:
            ninerow.set_ylim(mean - 10 * std, mean + 10 * std)
        except:
            pass

        ninerow.set_xlabel("wavelength (nm)", fontsize=10)
        ninerow.set_ylabel("residual (weighted)", fontsize=10)

        fig.suptitle(f"{self.arm} sky model: order {order}", fontsize=12, y=0.97)

        filename = self.filenameTemplate.replace(".fits", f"_SKYMODEL_QC_PLOTS_ORDER_{int(order)}.pdf")

        filePath = f"{self.qcDir}/{filename}"

        plt.show()
        plt.savefig(filePath, dpi=120, format="pdf")
        plt.close("all")

        self.log.debug("completed the ``plot_sky_sampling`` method")
        return filePath

    def rolling_window_clipping(self, imageMapOrderDF, windowSize, sigma_clip_limit=5, max_iterations=10):
        """*clip pixels in a rolling wavelength window*

        **Key Arguments:**

        - ``imageMapOrderDF`` --  dataframe with various processed data for a given order
        - ``windowSize`` -- the window-size used to perform rolling window clipping (number of data-points)
        - ``sigma_clip_limit`` -- clip data values straying beyond this sigma limit. Default *5*
        - ``max_iterations`` -- maximum number of iterations when clipping

        **Return:**

        - ``imageMapOrderDF`` -- image order dataframe with 'clipped' == True for those pixels that have been clipped via rolling window clipping

        **Usage:**

        ```python
        imageMapOrder = self.rolling_window_clipping(
            imageMapOrderDF=imageMapOrder,
            windowSize=23,
            sigma_clip_limit=4,
            max_iterations=10
        )
        ```
        """
        self.log.debug("starting the ``rolling_window_clipping`` method")

        import numpy as np
        from astropy.stats import sigma_clip, mad_std

        allPixels = len(imageMapOrderDF.index)
        order = imageMapOrderDF["order"].values[0]

        iteration = 0
        quantile = 0.3
        lastClipped = -1
        while iteration < max_iterations:
            iteration += 1
            # CALCULATE PERCENTILE SMOOTH DATA & RESIDUALS
            mask_clipped = imageMapOrderDF["flagged_all_clipped"] == True
            mask_object = imageMapOrderDF["flagged_object_clipped"] == True

            # RESETS
            imageMapOrderDF["flux_std"] = np.nan
            imageMapOrderDF["flux_upper_limit"] = np.nan
            imageMapOrderDF["flux_minus_smoothed_residual"] = np.nan
            imageMapOrderDF["flux_minus_smoothed_residual_std"] = np.nan
            imageMapOrderDF["flux_minus_smoothed_residual_upper_limit"] = np.nan

            if iteration <= 2:
                imageMapOrderDF["flux_percentile_smoothed"] = np.nan
                imageMapOrderDF.loc[~mask_clipped, "flux_percentile_smoothed"] = (
                    imageMapOrderDF.loc[~mask_clipped, "flux"]
                    .rolling(window=windowSize, center=True, closed="both")
                    .quantile(quantile)
                )

            imageMapOrderDF.loc[~mask_clipped, "flux_minus_smoothed_residual"] = (
                imageMapOrderDF.loc[~mask_clipped, "flux"]
                - imageMapOrderDF.loc[~mask_clipped, "flux_percentile_smoothed"]
            )

            # imageMapOrderDF.loc[~mask_clipped & ~mask_object, "flux_minus_smoothed_residual_std"] = (
            #     imageMapOrderDF.loc[~mask_clipped & ~mask_object, "flux_minus_smoothed_residual"]
            #     .rolling(window=windowSize * 10, center=True, min_periods=windowSize, closed="both")
            #     .std()
            # )

            imageMapOrderDF.loc[~mask_clipped, "flux_minus_smoothed_residual_std"] = (
                imageMapOrderDF.loc[~mask_clipped, "flux_minus_smoothed_residual"]
                .rolling(window=windowSize, center=True, closed="both")
                .std()
            )

            imageMapOrderDF.loc[~mask_clipped, "flux_minus_smoothed_residual_upper_limit"] = (
                imageMapOrderDF.loc[~mask_clipped, "flux_minus_smoothed_residual_std"] * sigma_clip_limit
            )

            imageMapOrderDF.loc[~mask_clipped, "flux_upper_limit"] = (
                imageMapOrderDF.loc[~mask_clipped, "flux_percentile_smoothed"]
                + imageMapOrderDF.loc[~mask_clipped, "flux_minus_smoothed_residual_upper_limit"]
            )

            # NOW CLIP THE OBJECT PIXELS
            masked = imageMapOrderDF.loc[~mask_clipped, "flux"] > imageMapOrderDF.loc[~mask_clipped, "flux_upper_limit"]
            imageMapOrderDF.loc[(~mask_clipped & ~mask_object & masked), "flagged_object_clipped"] = True

            # RESET MASKS
            mask_object = imageMapOrderDF["flagged_object_clipped"] == True
            imageMapOrderDF.loc[(mask_object), "flagged_all_clipped"] = True
            mask_clipped = imageMapOrderDF["flagged_all_clipped"] == True

            # masked = (
            #     imageMapOrderDF.loc[~mask_clipped & ~mask_object, "flux"]
            #     < imageMapOrderDF.loc[~mask_clipped & ~mask_object, "flux_percentile_smoothed"]
            #     - 1.0 * imageMapOrderDF.loc[~mask_clipped & ~mask_object, "flux_std"]
            # )
            # imageMapOrderDF.loc[(~mask_clipped & ~mask_object & masked), "flagged_object_clipped"] = True

            totalClipped = len(imageMapOrderDF.loc[(imageMapOrderDF["flagged_object_clipped"] == True)].index)
            percent = (float(totalClipped) / float(allPixels)) * 100.0
            # print(f"ORDER {order}: iteration {iteration} - {totalClipped} pixels clipped in total = {percent:1.1f}%")

            if totalClipped == lastClipped:
                # print("NO CHANGE IN CLIPPED PIXELS - STOPPING ITERATIONS")
                break
            lastClipped = totalClipped

            if self.debug and False:
                # PLOT THE CLIPPED PIXELS IN EACH ITERATION
                self.plot_order_skymodel_fitting_quicklook(
                    imageMapOrderDF,
                    None,
                    title=f"clipping pixels containing object flux\niteration {iteration} - {percent:1.1f}% clipped",
                )

            if self.arm.upper() in ("VIS"):
                if iteration == 5:
                    totalClipped = len(imageMapOrderDF.loc[(imageMapOrderDF["flagged_object_clipped"] == True)].index)
                    percent = (float(totalClipped) / float(allPixels)) * 100.0
                    if percent < 2:
                        print("FIXING PERCENT < 2")
                        imageMapOrderDF["flagged_object_clipped"] = False
                        sigma_clip_limit -= 0.1
                        if quantile > 0.1:
                            quantile -= 0.05
                        iteration = 0
                # if iteration == max_iterations:
                #     totalClipped = len(imageMapOrderDF.loc[(imageMapOrderDF["flagged_object_clipped"] == True)].index)
                #     percent = (float(totalClipped) / float(allPixels)) * 100.0
                #     if percent > 70:
                #         print("FIXING PERCENT > 70")
                #         imageMapOrderDF["flagged_object_clipped"] = False
                #         windowSize = windowSize + 1
                #         iteration = 0

        totalClipped = len(imageMapOrderDF.loc[(imageMapOrderDF["flagged_object_clipped"] == True)].index)
        percent = (float(totalClipped) / float(allPixels)) * 100.0

        # sys.stdout.flush()
        # sys.stdout.write("\x1b[1A\x1b[2K")
        percent = (float(totalClipped) / float(allPixels)) * 100.0
        self.log.print(f"\tORDER {order}: {totalClipped} pixels clipped in total = {percent:1.1f}%)")

        if percent > 85.0:
            imageMapOrderDF["flagged_object_clipped"] = False
            self.log.warning(
                f"ORDER {order}: More than 85% of pixels flagged to be clipped ({percent:1.1f}%). Clipping 0% instead."
            )

        std = imageMapOrderDF.loc[imageMapOrderDF["flagged_all_clipped"] == False, "flux_minus_smoothed_residual"].std()
        imageMapOrderDF["residual_global_sigma"] = imageMapOrderDF["flux_minus_smoothed_residual"] / std

        self.log.debug("completed the ``rolling_window_clipping`` method")
        return imageMapOrderDF

    def fit_bspline_curve_to_sky(self, imageMapOrder):
        """*fit a single-order univariate bspline to the unclipped sky pixels (wavelength vs flux)*

        **Key Arguments:**

        - ``imageMapOrder`` -- single order dataframe, containing sky flux with object(s) and CRHs removed
        - ``order`` -- the order number

        **Return:**

        - ``imageMapOrder`` -- same `imageMapOrder` as input but now with `sky_model` (bspline fit of the sky) and `sky_subtracted_flux` columns
        - ``tck`` -- the fitted bspline components. t for knots, c of coefficients, k for order

        **Usage:**

        ```python
        imageMapOrder, tck = self.fit_bspline_curve_to_sky(
            imageMapOrder
        )
        ```

        """
        self.log.debug("starting the ``fit_bspline_curve_to_sky`` method")

        import numpy as np
        import scipy.interpolate as ip
        import pandas as pd
        from astropy.stats import sigma_clip

        # CAN NOT ADD ANOTHER KNOT TO A GROUP OF DATA POINTS SMALLER THAN min_points_per_knot
        min_points_per_knot = self.recipeSettings["sky-subtraction"]["min_points_per_knot"]
        # WHEN REFINING THE BSPLINE FIT, USE THIS SIGMA FOR CLIPPING
        bsplineSigma = self.recipeSettings["sky-subtraction"]["bspline_fitting_residual_clipping_sigma"]
        # NUMBER OF ITERATIONS USED TO FIT THE BSPLINE
        bsplineIterations = self.recipeSettings["sky-subtraction"]["bspline_iteration_limit"]
        # USE THIS PERCENTILE TO DETERMINE THE RESIDUAL FLOOR
        residual_floor_percentile = self.recipeSettings["sky-subtraction"]["residual_floor_percentile"]

        # SORT BY COLUMN NAME
        imageMapOrder.sort_values(by=["wavelength"], inplace=True)
        order = imageMapOrder["order"].values[0]

        # CREATE ARRAYS NEEDED FOR BSPLINE FITTING
        mask_all_clipped = imageMapOrder["flagged_all_clipped"] == True
        goodWl = imageMapOrder.loc[~mask_all_clipped, "wavelength"]

        # REPLACE NANS IN residual_windowed_std
        data = imageMapOrder["residual_windowed_std"].values
        mask = np.isnan(data)
        data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
        imageMapOrder["residual_windowed_std"] = data

        imageMapOrder.loc[~mask_all_clipped, "weights"] = (
            1 / imageMapOrder.loc[~mask_all_clipped, "residual_windowed_std"].values
        )

        imageMapOrder.loc[~mask_all_clipped, "weights2"] = imageMapOrder.loc[~mask_all_clipped, "flux"].values + 1000
        imageMapOrder.loc[imageMapOrder["weights2"] < 0.01, "weights2"] = 0.01

        described_weights = imageMapOrder.loc[~mask_all_clipped, "weights2"].describe()
        if self.debug:
            print("weights2 description", described_weights)

        imageMapOrder.loc[~mask_all_clipped, "weights2"] = imageMapOrder.loc[~mask_all_clipped, "weights2"].values / (
            imageMapOrder.loc[~mask_all_clipped, "residual_windowed_std"].values * 2.0
        )

        described_weights = imageMapOrder.loc[~mask_all_clipped, "weights2"].describe()
        if self.debug:
            print("weights2 description", described_weights)

        imageMapOrder.loc[~mask_all_clipped, "weights"] = np.pow(
            imageMapOrder.loc[~mask_all_clipped, "weights2"].values, 1.2
        )
        imageMapOrder.loc[~mask_all_clipped, "weights2"] = np.pow(
            imageMapOrder.loc[~mask_all_clipped, "weights2"].values, 1.2
        )

        described_weights = imageMapOrder.loc[~mask_all_clipped, "weights2"].describe()
        if self.debug:
            print("weights2 description", described_weights)
        # imageMapOrder["weights"] = 1 / imageMapOrder["error"].values

        # WE WILL UPDATE THIS VALUE LATER IN WORKFLOW WITH SLIT-ILLUMINATION CORRECTION
        imageMapOrder["slit_normalisation_ratio"] = 1

        # FOR WEIGHTED BSPLINES WE ONLY NEED *INTERIOR* KNOTS (DON'T GO BEYOND RANGE OF DATA)
        # CAN'T HAVE MORE KNOTS THAN DATA POINTS
        # NUMBER OF 'DEFAULT' KNOTS
        defaultPointsPerKnot = self.recipeSettings["sky-subtraction"]["starting_points_per_knot"]

        # STARTER KNOTS USED TO MEASURE THE RESIDUAL FLOOR BEFORE REVERTING TO defaultPointsPerKnot
        if self.arm.upper() == "NIR":
            starterPointsPerKnot = 750
        else:
            starterPointsPerKnot = 1000

        if self.binx > 1:
            defaultPointsPerKnot /= self.binx
            starterPointsPerKnot /= self.binx
        if self.biny > 1:
            defaultPointsPerKnot /= self.biny
            starterPointsPerKnot /= self.biny

        if self.debug:
            print(f"defaultPointsPerKnot: {defaultPointsPerKnot}, order: {order}")
            print(f"min_points_per_knot: {min_points_per_knot}, order: {order}")
            print(f"starterPointsPerKnot: {starterPointsPerKnot}, order: {order}")

        defaultPointsPerKnot = int(defaultPointsPerKnot)
        starterPointsPerKnot = int(starterPointsPerKnot)

        if defaultPointsPerKnot:
            n_interior_knots = int(goodWl.values.shape[0] / defaultPointsPerKnot)
            # QUANTILE SPACES - i.e. PERCENTAGE VALUES TO PLACE THE KNOTS, FROM 0-1, ALONG WAVELENGTH RANGE
            try:
                qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
            except:
                qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
            defaultKnots = np.quantile(goodWl, qs)
        else:
            defaultKnots = np.array([])

        if self.debug:
            print("default knot count", len(defaultKnots))

        n_interior_knots = int(goodWl.values.shape[0] / starterPointsPerKnot)
        # QUANTILE SPACES - i.e. PERCENTAGE VALUES TO PLACE THE KNOTS, FROM 0-1, ALONG WAVELENGTH RANGE
        qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
        starterKnots = np.quantile(goodWl, qs)

        if self.debug:
            print("starter knot count", len(starterKnots))

        extraKnots = np.array([])
        iterationCount = -3
        residualFloor = False

        slitIlluminationCorrectionIteration = 4
        tiltAdjustmentIteration = 4

        slitCorrectIterationLimit = 2
        slitCorrectIterations = 0

        # CLIP NAN FLUX
        imageMapOrder.loc[imageMapOrder["flux"].isnull(), "flagged_all_clipped"] = True
        mask_all_clipped = imageMapOrder["flagged_all_clipped"] == True

        lastExtraKnotCount = -1
        while iterationCount < bsplineIterations:
            iterationCount += 1

            if iterationCount == 2:
                mask_noisy = imageMapOrder["flagged_noisy_region"] == True
                imageMapOrder.loc[mask_noisy, "weights2"] = (
                    imageMapOrder.loc[mask_noisy, "weights2"] / imageMapOrder.loc[mask_noisy, "flux"].abs() * 0.1
                )
                imageMapOrder.loc[mask_noisy, "weights"] = (
                    imageMapOrder.loc[mask_noisy, "weights"] / imageMapOrder.loc[mask_noisy, "flux"].abs() * 0.1
                )

            # CREATE ARRAYS NEEDED FOR BSPLINE FITTING
            goodWl = imageMapOrder.loc[~mask_all_clipped, "wavelength"]
            goodFlux = (
                imageMapOrder.loc[~mask_all_clipped, "flux"]
                / imageMapOrder.loc[~mask_all_clipped, "slit_normalisation_ratio"]
            )
            if iterationCount < 5:
                goodWeights = imageMapOrder.loc[~mask_all_clipped, "weights"]
            else:
                goodWeights = imageMapOrder.loc[~mask_all_clipped, "weights2"]

            baseFlux = np.median(goodFlux.values)
            goodFlux = goodFlux.values
            goodWeights = goodWeights.values
            goodFlux[0] = baseFlux
            goodFlux[-1] = baseFlux
            goodWeights[0] = 10e4
            goodWeights[-1] = 10e4

            if iterationCount < 5:
                baseKnots = starterKnots

            if iterationCount == 5:
                baseKnots = defaultKnots

            if self.debug:
                print("base knot count", len(baseKnots), iterationCount)
                print("extra knot count", len(extraKnots), iterationCount)

            allKnots = np.sort(np.concatenate((extraKnots, baseKnots)))

            if self.debug:
                print("total knot count", len(allKnots), iterationCount)

            if slitCorrectIterations < slitCorrectIterationLimit:

                # SLIT CORRECT ACTUALLY HELPS
                if iterationCount == slitIlluminationCorrectionIteration and False:
                    # FIT SLIT-ILLUMINATION PROFILE
                    imageMapOrder = self.cross_dispersion_flux_normaliser(imageMapOrder)
                    extraKnots = np.array([])
                    iterationCount = -2
                    slitIlluminationCorrectionIteration = -99
                    tiltAdjustmentIteration = 4

                # FLUX SHUFFLING MAKE EXTRACTION WORSE
                if iterationCount == tiltAdjustmentIteration and self.arm.upper() in ["NIR"] and False:
                    # FIT SLIT-ILLUMINATION PROFILE
                    imageMapOrder = self.adjust_tilt(imageMapOrder, tck)
                    extraKnots = np.array([])
                    iterationCount = -2
                    # tiltAdjustmentIteration = -99
                    slitCorrectIterations += 1

            if iterationCount > 1:
                # POTENTIAL NEW KNOTS PLACED HALF WAY BETWEEN ADJACENT CURRENT KNOTS
                meanResiduals = []

                # GROUP ALL DATA POINTS BETWEEN KNOTS
                nosiyRegionMask = imageMapOrder["flagged_noisy_region"] == True
                df = imageMapOrder.loc[~mask_all_clipped & ~nosiyRegionMask]
                ind = np.digitize(df["wavelength"], allKnots)

                ## GROUP ALL UNCLIPPED DATA POINTS BETWEEN KNOTS WHERE THERE ARE SKYLINES
                group = imageMapOrder.loc[~mask_all_clipped & ~nosiyRegionMask].groupby(ind)

                ## ANY GROUPS WITH A MEAN > 0 WILL CONTAIN HIGH RESIDUALS SOMEWHERE
                meanResiduals = group["sky_local_vs_global"].mean()

                counts = group.size()
                potentialNewKnots = group["wavelength"].mean()
                potentialNewKnots2 = group["wavelength"].mean() - group["wavelength"].std()
                potentialNewKnots3 = group["wavelength"].mean() + group["wavelength"].std()
                mask = counts < min_points_per_knot
                meanResiduals[mask] = -100000

                meanResiduals = np.array(meanResiduals)

                mask = np.ma.masked_where(meanResiduals > 0, meanResiduals).mask
                # ELSE ADD NEW KNOT IF ABOVE FLOOT
                for nk in [potentialNewKnots2, potentialNewKnots3]:
                    nk = np.ma.compressed(np.ma.masked_array(np.array(nk), ~mask))
                    allKnots = np.sort(np.concatenate((nk, allKnots)))
                    extraKnots = np.sort(np.concatenate((nk, extraKnots)))

                # NOW ADD KNOTS AT SKYLINE NODES
                skylineNodes = ~mask_all_clipped & (imageMapOrder["flagged_sky_line"] == "node")
                df = imageMapOrder.loc[skylineNodes]
                ind = np.digitize(df["wavelength"], allKnots)
                group = imageMapOrder.loc[skylineNodes].groupby(ind)
                counts = group.size()
                potentialNewKnots = group["wavelength"].mean()
                mask = counts < min_points_per_knot
                newKnots = np.ma.compressed(np.ma.masked_array(potentialNewKnots, mask))
                allKnots = np.sort(np.concatenate((newKnots, allKnots)))
                extraKnots = np.sort(np.concatenate((newKnots, extraKnots)))

                # if order == self.qcPlotOrder:
                #     print(f"EXTRA KNOTS: {len(newKnots)} .... {len(allKnots)} ... {iterationCount}")

            try:

                tck, fp, ier, msg = ip.splrep(
                    goodWl, goodFlux, t=allKnots, k=self.bspline_order, w=goodWeights, full_output=True
                )
            except:
                raise ValueError(
                    f"BSpline fit failed for order {order} on iteration {iterationCount}. Possibly too many knots ({len(allKnots)}) for the number of data points ({goodWl.values.shape[0]})."
                )
            t, c, k = tck

            if ier in (10, 30):
                self.log.info(
                    f"\t\tpoor fit on iteration {iterationCount} for order {imageMapOrder['order'].values[0]}. Reverting to last iteration.\n"
                )
                tck = tck_previous
                break
            else:
                tck_previous = tck

            if iterationCount >= -1:
                # FIRST PASS SIGMA CLIPPING OF BSPLINE
                for _ in range(3):
                    mask_all_clipped = imageMapOrder["flagged_all_clipped"] == True

                    # ## ROLLING MEDIAN CLIPPING
                    # imageMapOrder.loc[~mask_all_clipped, "sky_flux_rolling_median"] = (
                    #     imageMapOrder.loc[~mask_all_clipped, "flux"]
                    #     .rolling(35, center=True, min_periods=3, closed="both")
                    #     .median()
                    # )
                    # imageMapOrder.loc[~mask_all_clipped, "sky_flux_rolling_std"] = (
                    #     imageMapOrder.loc[~mask_all_clipped, "flux"]
                    #     .rolling(35, center=True, min_periods=3, closed="both")
                    #     .std()
                    # )
                    # mask_rolling_median_clipping = ~mask_all_clipped & (
                    #     (
                    #         imageMapOrder["flux"]
                    #         > imageMapOrder["sky_flux_rolling_median"] + imageMapOrder["sky_flux_rolling_std"] * 3
                    #     )
                    #     | (
                    #         imageMapOrder["flux"]
                    #         < imageMapOrder["sky_flux_rolling_median"] - imageMapOrder["sky_flux_rolling_std"] * 2
                    #     )
                    # )
                    # imageMapOrder.loc[mask_rolling_median_clipping, "flagged_all_clipped"] = True
                    # imageMapOrder.loc[mask_rolling_median_clipping, "flagged_bspline_clipped"] = True

                    # residuals = imageMapOrder.loc[~mask_all_clipped, "sky_subtracted_flux"]
                    # imageMapOrder.loc[~mask_all_clipped, "bspline_sky_residual_windowed_std"] = (
                    #     imageMapOrder.loc[~mask_all_clipped, "sky_subtracted_flux"]
                    #     .rolling(15, center=True, min_periods=3, closed="both")
                    #     .std()
                    # )

                    # mask_residual_clipping = ~mask_all_clipped & (
                    #     (
                    #         imageMapOrder["sky_subtracted_flux"]
                    #         > imageMapOrder["bspline_sky_residual_windowed_std"] * bsplineSigma
                    #     )
                    #     | (
                    #         imageMapOrder["sky_subtracted_flux"]
                    #         < -imageMapOrder["bspline_sky_residual_windowed_std"] * bsplineSigma
                    #     )
                    # )
                    # imageMapOrder.loc[mask_residual_clipping, "flagged_bspline_clipped"] = True
                    # imageMapOrder.loc[mask_residual_clipping, "flagged_all_clipped"] = True

            if iterationCount > 0:
                imageMapOrder, residualFloor = self.determine_residual_floor(imageMapOrder, tck, iterationCount)

            if iterationCount > 5:
                window = 15
            else:
                window = 25

            if self.arm.upper() in ["NIR"]:
                window = 50

            imageMapOrder.loc[~mask_all_clipped, "sky_residuals"] = imageMapOrder.loc[
                ~mask_all_clipped, "flux"
            ].values - ip.splev(imageMapOrder.loc[~mask_all_clipped, "wavelength"].values, tck)

            def sliding_median(arr, window):
                if window % 2 == 0:
                    window += 1

                medianArry = np.median(np.lib.stride_tricks.sliding_window_view(arr, (window,)), axis=1)
                start = np.ones(window // 2) * medianArry[0]
                end = np.ones(window // 2) * medianArry[-1]
                medianArry = np.concatenate((start, medianArry, end))
                return medianArry

            ## USE ROLLING MEAN TO ESTIMATE THE LOCAL RESIDUALS, WHICH CAN BE USED TO ADD NEW KNOTS IN HIGH-RESIDUAL AREAS
            imageMapOrder.loc[~mask_all_clipped, "sky_residual_floor_local"] = sliding_median(
                imageMapOrder.loc[~mask_all_clipped, "sky_residuals"].values, window=window
            )

            if self.debug and iterationCount > 4:
                self.plot_order_skymodel_fitting_quicklook(
                    imageMapOrder, tck, title=f"Fitting the sky model\niteration {iterationCount}", knots=allKnots
                )

            # GENERATE SKY-MODEL FROM BSPLINE
            imageMapOrder["sky_model_wl"] = ip.splev(imageMapOrder["wavelength"].values, tck)
            imageMapOrder["sky_model_wl_derivative"] = ip.splev(imageMapOrder["wavelength"].values, tck, der=1)
            imageMapOrder["sky_subtracted_flux"] = (
                imageMapOrder["flux"] / imageMapOrder["slit_normalisation_ratio"]
            ) - imageMapOrder["sky_model_wl"]
            imageMapOrder["sky_subtracted_flux_weighted"] = (
                imageMapOrder["sky_subtracted_flux"]
                * imageMapOrder["sky_model_wl_derivative"].abs()
                / imageMapOrder["residual_windowed_std"]
            )
            imageMapOrder["sky_subtracted_flux_weighted_abs"] = imageMapOrder["sky_subtracted_flux_weighted"].abs()

            flux_error_ratio = imageMapOrder.loc[
                imageMapOrder["flagged_all_clipped"] == False, "sky_subtracted_flux_weighted"
            ].values

            if flux_error_ratio[1000:-1000].shape[0]:
                flux_error_ratio = flux_error_ratio[1000:-1000]

            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            try:
                self.log.print(
                    f"\tOrder: {order}, Iteration {iterationCount}, RES {flux_error_ratio.mean():0.3f}, STD {flux_error_ratio.std():0.3f}, MEDIAN {np.median(flux_error_ratio):0.3f}, MAX {flux_error_ratio.max():0.3f}, MIN {flux_error_ratio.min():0.3f}"
                )
            except:
                pass

            if iterationCount >= 5:
                if lastExtraKnotCount == len(extraKnots):
                    self.log.info(f"\t\tNo new knots added on iteration {iterationCount}. Stopping iterations.\n")
                    break
            lastExtraKnotCount = len(extraKnots)

        if not lastExtraKnotCount:
            imageMapOrder["sky_model_wl"] = baseFlux
            imageMapOrder["sky_model_wl_derivative"] = 1
            imageMapOrder["sky_model"] = baseFlux

            imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model"]
            imageMapOrder["sky_subtracted_flux_weighted"] = 1
            imageMapOrder["sky_subtracted_flux_weighted_abs"] = imageMapOrder["sky_subtracted_flux_weighted"].abs()
        else:
            imageMapOrder["sky_model_wl"] = ip.splev(imageMapOrder["wavelength"].values, tck)
            imageMapOrder["sky_model_wl_derivative"] = ip.splev(imageMapOrder["wavelength"].values, tck, der=1)
            imageMapOrder["sky_model"] = imageMapOrder["sky_model_wl"] * imageMapOrder["slit_normalisation_ratio"]
            # REPLACE VALUES LESS THAN ZERO IN COLUMN WITH ZERO
            imageMapOrder["sky_model"] = imageMapOrder["sky_model"].apply(lambda x: max(0, x))

            imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model"]
            imageMapOrder["sky_subtracted_flux_weighted"] = (
                imageMapOrder["sky_subtracted_flux"]
                * imageMapOrder["sky_model_wl_derivative"].abs()
                / (imageMapOrder["residual_windowed_std"] * 10)
            )

        imageMapOrder["sky_subtracted_flux_weighted_abs"] = imageMapOrder["sky_subtracted_flux_weighted"].abs()
        flux_error_ratio = imageMapOrder.loc[
            imageMapOrder["flagged_all_clipped"] == False, "sky_subtracted_flux_weighted"
        ].values
        if flux_error_ratio[1000:-1000].shape[0]:
            flux_error_ratio = flux_error_ratio[1000:-1000]

        self.log.debug("completed the ``fit_bspline_curve_to_sky`` method")
        return imageMapOrder, tck, allKnots, flux_error_ratio, residualFloor

    def create_placeholder_images(self):
        """*create placeholder images for the sky model and sky-subtracted frame*

        **Return:**

        - ``skymodelCCDData`` -- placeholder for sky model image
        - ``skySubtractedCCDData`` -- placeholder for sky-subtracted image

        **Usage:**

        ```python
        skymodelCCDData, skySubtractedCCDData = self.create_placeholder_images()
        ```
        """
        self.log.debug("starting the ``create_placeholder_images`` method")

        import numpy as np

        # CREATE AN IMAGE ARRAY TO HOST WAVELENGTH AND SLIT-POSITIONS
        skymodelCCDData = self.objectFrame.copy()
        skymodelCCDData.data[:] = np.nan
        skySubtractedCCDData = skymodelCCDData.copy()
        skySubtractedResidualsCCDData = skymodelCCDData.copy()

        self.log.debug("completed the ``create_placeholder_images`` method")
        return skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData

    def add_data_to_placeholder_images(
        self, imageMapOrderDF, skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData
    ):
        """*add sky-model and sky-subtracted data to placeholder images*

        **Key Arguments:**

        - ``imageMapOrderDF`` -- single order dataframe from object image and 2D map
        - ``skymodelCCDData`` -- the sky model
        - ``skySubtractedCCDData`` -- the sky-subtracted data

        **Usage:**

        ```python
        self.add_data_to_placeholder_images(imageMapOrderSkyOnly, skymodelCCDData, skySubtractedCCDData)
        ```
        """
        self.log.debug("starting the ``add_data_to_placeholder_images`` method")

        for x, y, skypixel in zip(
            imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_model"]
        ):
            if self.detectorParams["dispersion-axis"] == "x":
                skymodelCCDData.data[y][x] = skypixel
            else:
                skymodelCCDData.data[x][y] = skypixel

        for x, y, skypixel in zip(
            imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_subtracted_flux"]
        ):
            if self.detectorParams["dispersion-axis"] == "x":
                skySubtractedCCDData.data[y][x] = skypixel
            else:
                skySubtractedCCDData.data[x][y] = skypixel
        for x, y, skypixel in zip(
            imageMapOrderDF[self.axisA],
            imageMapOrderDF[self.axisB],
            imageMapOrderDF["sky_subtracted_flux"] / imageMapOrderDF["error"],
        ):
            if self.detectorParams["dispersion-axis"] == "x":
                skySubtractedResidualsCCDData.data[y][x] = skypixel
            else:
                skySubtractedResidualsCCDData.data[x][y] = skypixel

        self.log.debug("completed the ``add_data_to_placeholder_images`` method")
        return skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData

    def plot_image_comparison(self, objectFrame, skyModelFrame, skySubFrame):
        """*generate a plot of original image, sky-model and sky-subtraction image*

        **Key Arguments:**

        - ``objectFrame`` -- object frame
        - ``skyModelFrame`` -- sky model frame
        - ``skySubFrame`` -- sky subtracted frame

        **Return:**

        - ``filePath`` -- path to the plot pdf
        """
        self.log.debug("starting the ``plot_results`` method")

        import matplotlib.pyplot as plt
        import numpy.ma as ma
        import numpy as np

        arm = self.arm

        # a = plt.figure(figsize=(40, 15))
        if arm == "UVB":
            fig = plt.figure(figsize=(6, 13.5), constrained_layout=True)
        else:
            fig = plt.figure(figsize=(6, 11), constrained_layout=True)
        gs = fig.add_gridspec(6, 4)

        # CREATE THE GRID OF AXES
        toprow = fig.add_subplot(gs[0:2, :])
        midrow = fig.add_subplot(gs[2:4, :])
        bottomrow = fig.add_subplot(gs[4:6, :])

        # FIND ORDER PIXELS - MASK THE REST
        nonOrderMask = np.ones_like(objectFrame.data)
        for x, y in zip(self.mapDF[self.axisA], self.mapDF[self.axisB]):
            nonOrderMask[y][x] = 0

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        nonOrderMask = ma.make_mask(nonOrderMask)
        combinedMask = (nonOrderMask == 1) | (objectFrame.mask == 1)

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.rot90(objectFrame.data, 1)
        maskedDataArray = np.ma.array(objectFrame.data, mask=combinedMask)
        maskedDataValues = np.array(maskedDataArray.filled(np.nan), dtype=float, copy=True)

        try:
            std = np.nanstd(maskedDataValues)
            mean = np.nanmean(maskedDataValues)
        except:
            std = np.std(maskedDataValues)
            mean = np.mean(maskedDataValues)
        vmax = mean + 1 * std
        vmin = mean - 0.1 * std
        toprow.imshow(rotatedImg, vmin=0, vmax=100, cmap="gray", alpha=1.0)
        toprow.set_title(f"Original {arm} Frame", fontsize=10)
        toprow.set_ylabel("x-axis", fontsize=8)
        toprow.set_xlabel("y-axis", fontsize=8)
        toprow.tick_params(axis="both", which="major", labelsize=9)

        rotatedImg = np.rot90(skyModelFrame.data, 1)
        maskedDataArray = np.ma.array(skyModelFrame.data, mask=combinedMask)
        maskedDataValues = np.array(maskedDataArray.filled(np.nan), dtype=float, copy=True)
        std = np.nanstd(maskedDataValues)
        mean = np.nanmean(maskedDataValues)
        vmax = mean + 1 * std
        vmin = mean - 1 * std
        midrow.imshow(rotatedImg, vmin=0, vmax=100, cmap="gray", alpha=1.0)
        midrow.set_title(f"Sky-model for {arm} Frame", fontsize=10)
        midrow.set_ylabel("x-axis", fontsize=8)
        midrow.set_xlabel("y-axis", fontsize=8)
        midrow.tick_params(axis="both", which="major", labelsize=9)

        rotatedImg = np.rot90(skySubFrame.data, 1)
        maskedDataArray = np.ma.array(skySubFrame.data, mask=combinedMask)
        maskedDataValues = np.array(maskedDataArray.filled(np.nan), dtype=float, copy=True)

        try:
            std = np.nanstd(maskedDataValues)
            mean = np.nanmean(maskedDataValues)
        except:
            std = np.std(maskedDataValues)
            mean = np.mean(maskedDataValues)

        vmax = 0 + std
        vmin = 0
        bottomrow.imshow(rotatedImg, vmin=vmin, vmax=30, cmap="gray", alpha=1.0)
        bottomrow.set_title(f"Sky-subtracted {arm} Frame", fontsize=10)
        bottomrow.set_ylabel("x-axis", fontsize=8)
        bottomrow.set_xlabel("y-axis", fontsize=8)
        bottomrow.tick_params(axis="both", which="major", labelsize=9)

        # plt.show()
        filename = self.filenameTemplate.replace(".fits", "_skysub_quicklook.pdf")

        filePath = f"{self.qcDir}/{filename}"
        plt.savefig(filePath, dpi=120, format="pdf")
        plt.close("all")

        self.log.debug("completed the ``plot_results`` method")
        return filePath

    def _rectify_order(self, order, imageMapOrder, remove_clipped=False, conserve_flux=False):
        """*rectify order on a fine slit-position, wavelength grid*

        **Key Arguments:**

        - ``order`` -- order to be rectified
        - ``imageMapOrder`` -- the image map for this order (wavelength, slit-position and flux for each physical pixel
        - ``conserve_flux`` -- conserve the flux budget across the entire image

        **Return:**

        - None

        **Usage:**

        ```python
        usage code
        ```
        """
        self.log.debug("starting the ``rectify_order`` method")

        import numpy as np
        import pandas as pd

        dispMap = self.dispMap
        kw = self.kw
        dp = self.detectorParams
        arm = self.arm

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm
        )

        for o, minWl, maxWl in zip(orderNums, waveLengthMin, waveLengthMax):
            if o == order:
                orderInfo = (order, minWl, maxWl)
        order, minWl, maxWl = orderInfo

        minWl = minWl - 5
        maxWl = minWl + 5

        # DYANIMICALLY DETERMINE SIZE OF SUB-PIXELS
        slit_pixel_range = imageMapOrder[self.axisA].max() - imageMapOrder[self.axisA].min()
        wl_pixel_range = imageMapOrder[self.axisB].max() - imageMapOrder[self.axisB].min()

        wl_range = maxWl - minWl
        slitLength = dp["slit_length"]
        slitLength = 4
        sl_range = dp["slit_length"]
        sl_range = 4

        straighten_grid_res_wavelength = 2 * (wl_range / wl_pixel_range)  # in nm
        straighten_grid_res_slit = 2 * (sl_range / slit_pixel_range)  # in arcsec

        halfGrid = slitLength / 2
        slitArray = np.arange(-halfGrid, halfGrid + straighten_grid_res_slit, straighten_grid_res_slit)

        wlArray = np.arange(minWl, maxWl, straighten_grid_res_wavelength)

        # ONE SINGLE-VALUE SLIT ARRAY FOR EVERY WAVELENGTH ARRAY
        bigSlitArray = np.concatenate([np.ones(wlArray.shape[0]) * slitArray[i] for i in range(0, slitArray.shape[0])])
        # NOW THE BIG WAVELEGTH ARRAY
        bigWlArray = np.tile(wlArray, np.shape(slitArray)[0])

        # CREATE PANDAS DATAFRAME WITH LARGE ARRAYS - ONE ROW PER
        # WAVELENGTH-SLIT GRID CELL
        myDict = {
            "order": np.ones(bigWlArray.shape[0]) * order,
            "wavelength": bigWlArray,
            "slit_position": bigSlitArray,
        }
        orderPixelTable = pd.DataFrame(myDict)

        # GET DETECTOR PIXEL POSITIONS FOR ALL WAVELENGTH-SLIT GRID CELLS
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispMap,
            orderPixelTable=orderPixelTable,
            removeOffDetectorLocation=False,
        )
        # INTEGER PIXEL VALUES & FIT DISPLACEMENTS FROM PIXEL CENTRES
        orderPixelTable["pixel_x"] = np.floor(orderPixelTable["fit_x"].values)
        orderPixelTable["pixel_y"] = np.floor(orderPixelTable["fit_y"].values)

        # xpd-update-filter-dataframe-column-values

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        # mask = (orderPixelTable["pixel_x"] < self.objectFrame.shape[1]) & (orderPixelTable["pixel_y"] < self.objectFrame.shape[0])
        # orderPixelTable = orderPixelTable.loc[mask]

        # xpd-update-filter-dataframe-column-values

        pixel_x = orderPixelTable["pixel_x"].values.astype(int)
        pixel_y = orderPixelTable["pixel_y"].values.astype(int)

        # fluxValues = self.objectFrame.data[pixel_y, pixel_x].byteswap().newbyteorder()
        # try:
        #     orderPixelTable["flux"] = fluxValues.byteswap().newbyteorder()
        #     orderPixelTable.sort_values(['slit_position', 'wavelength'])
        # except:
        #     orderPixelTable["flux"] = fluxValues
        #     orderPixelTable.sort_values(['slit_position', 'wavelength'])

        orderPixelTable = pd.merge(
            orderPixelTable,
            imageMapOrder[["x", "y", "flux", "flagged_all_clipped"]],
            how="left",
            left_on=["pixel_x", "pixel_y"],
            right_on=["x", "y"],
        )

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = orderPixelTable["flux"].isnull()
        self.log.print(orderPixelTable.loc[~mask, "wavelength"].min())

        # DROP MISSING VALUES
        # orderPixelTable.dropna(axis='index', how='any', subset=['x'], inplace=True)

        # orderPixelTable = orderPixelTable[['order', 'wavelength', 'slit_position', 'fit_x', 'fit_y', 'flux', 'clipped']]
        # orderPixelTable['weight'] = 100

        if conserve_flux:
            # ADD A COUNT COLUMN FOR THE NUMBER OF SMALL SLIT/WL PIXELS FALLING IN LARGE DETECTOR PIXELS
            count = orderPixelTable.groupby(["pixel_x", "pixel_y"]).size().reset_index(name="count")
            orderPixelTable = pd.merge(
                orderPixelTable, count, how="left", left_on=["pixel_x", "pixel_y"], right_on=["pixel_x", "pixel_y"]
            )

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        if remove_clipped:
            mask = orderPixelTable["flagged_all_clipped"] == True
            orderPixelTable.loc[mask, "flux"] = np.nan

        # RESTRUCTURE FLUXES INTO A STRAIGHTENED IMAGE
        imageArray = np.array([])
        for index, slit in enumerate(slitArray):
            rowFlux = orderPixelTable[(orderPixelTable["slit_position"] == slit)]["flux"].values
            if index == 0:
                imageArray = rowFlux
            else:
                imageArray = np.vstack((imageArray, rowFlux))

        imageArray[imageArray > 80000] = np.nan
        imageArray[imageArray < -7000] = np.nan

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=imageArray,
            show=False,
            ext="data",
            stdWindow=3,
            title=False,
            surfacePlot=True,
            inst="dummy",
        )

        self.log.debug("completed the ``rectify_order`` method")
        return imageArray

    def calculate_residuals(self, skyPixelsDF, fluxcoeff, orderDeg, wavelengthDeg, slitDeg, writeQCs=False):
        """*calculate residuals of the polynomial fits against the observed line positions*

        **Key Arguments:**

        - ``skyPixelsDF`` -- the predicted line list as a data frame
        - ``fluxcoeff`` -- the flux-coefficients
        - ``orderDeg`` -- degree of the order fitting
        - ``wavelengthDeg`` -- degree of wavelength fitting
        - ``slitDeg`` -- degree of the slit fitting (False for single pinhole)
        - ``writeQCs`` -- write the QCs to dataframe? Default *False*

        **Return:**

        - ``residuals`` -- combined x-y residuals
        - ``mean`` -- the mean of the combine residuals
        - ``std`` -- the stdev of the combine residuals
        - ``median`` -- the median of the combine residuals
        """
        self.log.debug("starting the ``calculate_residuals`` method")

        import numpy as np

        arm = self.arm

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # POLY FUNCTION NEEDS A DATAFRAME AS INPUT
        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDeg, wavelengthDeg=wavelengthDeg, slitDeg=slitDeg, exponentsIncluded=True
        ).poly

        # CALCULATE RESIDUALS BETWEEN MEASURED FLUX AND POLY
        # FITTED FLUX
        skyPixelsDF["fit_sky_subtracted_flux"] = poly(skyPixelsDF, *fluxcoeff)
        skyPixelsDF["residuals_sky_subtracted_flux"] = (
            skyPixelsDF["fit_sky_subtracted_flux"] - skyPixelsDF["sky_subtracted_flux"]
        )

        # CALCULATE COMBINED RESIDUALS AND STATS
        res_mean = np.mean(skyPixelsDF["residuals_sky_subtracted_flux"])
        res_std = np.std(skyPixelsDF["residuals_sky_subtracted_flux"])
        res_median = np.median(skyPixelsDF["residuals_sky_subtracted_flux"])

        self.log.debug("completed the ``calculate_residuals`` method")
        return res_mean, res_std, res_median, skyPixelsDF

    def clip_object_slit_positions(self, order_dataframes, aggressive=False):
        """*clip out pixels flagged as an object*

        **Key Arguments:**

        - ``order_dataframes`` -- a list of order data-frames with pixels potentially containing the object flagged.

        **Return:**

        - ``order_dataframes`` -- the order dataframes with the object(s) slit-ranges clipped
        - ``sky_only_dataframes`` -- dataframes with object removed

        **Usage:**

        ```python
        allimageMapOrder = self.clip_object_slit_positions(
            allimageMapOrder,
            aggressive=True
        )
        ```
        """
        self.log.debug("starting the ``clip_object_slit_positions`` method")

        import numpy as np
        import pandas as pd

        # COMBINE ALL ORDERS AND KEEP ONLY PIXELS FLAGGED AS POTENTIAL OBJECT
        allimageMapOrder = pd.concat(order_dataframes)
        mask = allimageMapOrder["flagged_object_clipped"] == True
        allimageMapOrder = allimageMapOrder.loc[mask]

        percentile_rolling_window_size = self.recipeSettings["sky-subtraction"]["percentile_rolling_window_size"]
        noise_rolling_window_size = self.recipeSettings["sky-subtraction"]["noise_rolling_window_size"]

        if aggressive:

            # BIN FLAGGED PIXEL COUNTS INTO DISCRETE SLIT-POSITION RANGES
            nbins = 100
            minsp = allimageMapOrder["slit_position"].min()
            maxsp = allimageMapOrder["slit_position"].max()
            bins = np.linspace(minsp, maxsp, nbins)
            result = allimageMapOrder["slit_position"].value_counts(bins=bins, sort=False, normalize=True) * nbins

            # REMOVE MEDIAN x 3 -- ONLY OBJECTS SHOULD REMAIN POSITIVE IN COUNTS
            result -= result.median()
            result -= result.abs().median()
            # result -= result.abs().median()

            # NEED 3 POSITIVE BINS IN A ROW TO BE SELECTED AS AN OBJECT
            object_ranges = []
            postiveCount = 0
            # AVOID EDGES WHEN SELECTING OBJECT SLIT-POSITIONS
            edges = int(nbins / 20)
            lower = False
            for sp, count in zip(bins[edges:-edges], result[edges:-edges]):
                if count > 0:
                    postiveCount += 1
                    upper = sp
                    if count > 0.05:
                        record_range = True
                else:
                    if postiveCount > 4 and record_range:
                        object_ranges.append([lower, upper])
                    postiveCount = 0
                    lower = sp
                    upper = False
                    record_range = False
            if postiveCount > 4:
                object_ranges.append([lower, upper])

            if 1 == 0:
                import matplotlib.pyplot as plt

                self.log.print(object_ranges)
                width = (maxsp - minsp) / nbins
                fig, ax = plt.subplots()
                bins = bins[:-1]
                rects1 = ax.bar(bins - width / 2, result, width, label="count")
                fig.tight_layout()
                plt.show()

        # NOW FOR EACH OBJECT SLIT-RANGE, FLAG AS CLIPPED IN ORIGINAL ORDER DATAFRAMES
        for df in order_dataframes:
            if aggressive:
                for objectt in object_ranges:

                    df.loc[(df["slit_position"].between(objectt[0], objectt[1])), "flagged_all_clipped"] = True
                    df.loc[(df["slit_position"].between(objectt[0], objectt[1])), "flagged_object_clipped"] = True
                    df.loc[((df["flagged_object_clipped"] == True)), "flagged_all_clipped"] = True
            else:
                # df.loc[((df['slit_position'].between(object[0], object[1])) & (df['object'] == True)), "flagged_all_clipped"] = True
                df.loc[((df["flagged_object_clipped"] == True)), "flagged_all_clipped"] = True
            # df.loc[
            #     ((df["flagged_all_clipped"] == False) & (df["flagged_object_clipped"] == True)),
            #     "flagged_object_clipped",
            # ] = False

            notClippedOrObjectMask = (df["flagged_all_clipped"] == False) & (df["flagged_object_clipped"] == False)

            df.loc[notClippedOrObjectMask, "residual_windowed_std"] = (
                df.loc[notClippedOrObjectMask, "flux_minus_smoothed_residual"]
                .rolling(percentile_rolling_window_size, center=True, min_periods=3, closed="both")
                .std()
            )

            # GENERATE A LARGER ROLLING WINDOWED STD TO BE USED IN DETERMINING NOISE
            df.loc[notClippedOrObjectMask, "residual_windowed_long_median"] = (
                df.loc[notClippedOrObjectMask, "flux_minus_smoothed_residual"]
                .rolling(noise_rolling_window_size, center=True, min_periods=30, closed="both")
                .median()
            )

            # GENERATE A LARGER ROLLING WINDOWED STD TO BE USED IN DETERMINING NOISE
            df.loc[notClippedOrObjectMask, "flux_windowed_long_median"] = (
                df.loc[notClippedOrObjectMask, "flux_percentile_smoothed"]
                .rolling(noise_rolling_window_size, center=True, min_periods=30, closed="both")
                .median()
            )

        self.log.debug("completed the ``clip_object_slit_positions`` method")
        return order_dataframes

    def cross_dispersion_flux_normaliser(self, orderDF):
        """*measure and normalise the flux in the cross-dispersion direction*

        **Key Arguments:**

        - ``orderDF`` -- a single order dataframe containing sky-subtraction flux residuals used to determine and remove a slit-illumination correction

        **Return:**

        - `correctedOrderDF` -- dataframe with slit-illumination correction factor added (flux-normaliser)

        **Usage:**

        ```python
        correctedOrderDF = self.cross_dispersion_flux_normaliser(orderDF)
        ```
        """
        self.log.debug("starting the ``cross_dispersion_flux_normaliser`` method")

        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt
        import pandas as pd
        import scipy.interpolate as ip

        slit_illumination_order = self.recipeSettings["sky-subtraction"]["slit_illumination_order"]
        order = orderDF["order"].values[0]
        orderDF["slit_normalisation_ratio"] = 1.0

        # COLLECT THE SKYLINE PIXELS
        mask = orderDF["flagged_all_clipped"] == False
        # mask = ((orderDF["flagged_all_clipped"] == False) & (orderDF["flagged_sky_line"] == False))
        thisOrder = orderDF.loc[mask]

        # SORT BY COLUMN NAME
        thisOrder.sort_values(["slit_position"], inplace=True)

        # GROUP INTO DISCRETE WAVELENGTH BINS
        # DEFINE THE BINS FOR COLUMN 'wavelength'
        lower = float("%0.*f" % (1, thisOrder["slit_position"].min()))
        upper = float("%0.*f" % (1, thisOrder["slit_position"].max()))

        # DO AN INITIAL CLIPPING OF THE DATA
        # masked_residuals = sigma_clip(
        #     thisOrder["sky_subtracted_flux"], sigma_lower=3, sigma_upper=3, maxiters=5, cenfunc='mean', stdfunc='std')
        # thisOrder["flagged_all_clipped"] = masked_residuals.mask

        # FIT THE DATA WITH A POLYMONIAL
        mask = thisOrder["flagged_all_clipped"] == False
        sp = thisOrder.loc[mask]["slit_position"].values
        fx = thisOrder.loc[mask]["sky_subtracted_flux"].values
        wt = 1 / thisOrder.loc[mask]["residual_windowed_std"].values

        # STARTER KNOTS USED TO MEASURE THE RESIDUAL FLOOR BEFORE REVERTING TO defaultPointsPerKnot
        starterPointsPerKnot = 500

        if self.binx > 1:

            starterPointsPerKnot /= self.binx
        if self.biny > 1:
            starterPointsPerKnot /= self.biny

        starterPointsPerKnot = int(starterPointsPerKnot)

        n_interior_knots = int(sp.shape[0] / starterPointsPerKnot)
        # QUANTILE SPACES - i.e. PERCENTAGE VALUES TO PLACE THE KNOTS, FROM 0-1,  WAVELENGTH RANGE
        qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
        starterKnots = np.quantile(sp, qs)

        tck, fp, ier, msg = ip.splrep(sp, fx, k=1, w=wt, t=starterKnots, full_output=True)

        iteration = 0
        spCopy = np.copy(sp)
        fxCopy = np.copy(fx)
        while iteration < 3 and True:
            iteration += 1
            coeff = np.polyfit(spCopy, fxCopy, deg=slit_illumination_order)
            residuals = fxCopy - np.polyval(coeff, spCopy)
            masked_residuals = sigma_clip(
                residuals, sigma_lower=3, sigma_upper=3, maxiters=1, cenfunc="mean", stdfunc="std"
            )
            # REDUCE ARRAYS TO NON-MASKED VALUES
            a = [spCopy, fxCopy]
            spCopy, fxCopy = [np.ma.compressed(np.ma.masked_array(i, masked_residuals.mask)) for i in a]

        # FIX ME
        if self.debug:
            fig = plt.figure(figsize=(15, 4))
            plt.title("Slit Illumination Profile")

            plt.plot(sp, fx, ".", ms=1, color="red", alpha=0.3)

            xp = np.linspace(lower, upper, 100)
            fitFx = np.polyval(coeff, xp)
            plt.plot(xp, fitFx, color="blue")

            slitPos = np.linspace(lower, upper, 10000)
            sky = ip.splev(slitPos, tck)
            plt.plot(slitPos, sky, label="sky model", c="#268bd2", zorder=3)

            plt.xlabel("slit-position", fontsize=12)
            plt.ylabel("normalised flux", fontsize=12)
            plt.show()
            plt.close("all")

        # USE THE SLIT-ILLUMINATION FUNCTION TO CREATE A FLUX-NORMALISATION

        if False:
            orderDF["flux"] -= np.polyval(coeff, orderDF["slit_position"].values)
        # orderDF['flux'] -= ip.splev(orderDF['slit_position'].values, tck)
        orderDF["slit_normalisation_ratio"] = 1

        self.log.debug("completed the ``cross_dispersion_flux_normaliser`` method")
        return orderDF

    def adjust_tilt(self, orderDF, tck):
        """*correct the tilt of the slit by attempting to minimise the residuals of the bspline fit while shifting the tilt angle*

        **Key Arguments:**

        - ``orderDF`` -- a single order dataframe containing sky-subtraction flux residuals used to determine and remove a slit-illumination correction
        - ``tck`` -- spline parameters to replot

        **Return:**

        - `orderDF` -- dataframe with wavelengths adjusted for a corrected tilt

        **Usage:**

        ```python
        orderDF = self.adjust_tilt(orderDF, tck)
        ```
        """
        self.log.debug("starting the ``adjust_tilt`` method")

        import scipy.interpolate
        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt
        import pandas as pd
        import scipy.interpolate as ip

        slit_illumination_order = self.recipeSettings["sky-subtraction"]["slit_illumination_order"]

        # COLLECT THE SKYLINE PIXELS
        # mask = ((orderDF["flagged_all_clipped"] == False) & (orderDF["flagged_sky_line"] != False))
        mask = orderDF["flagged_all_clipped"] == False
        thisOrder = orderDF.loc[mask]
        order = orderDF["order"].values[0]

        # GROUP INTO DISCRETE WAVELENGTH BINS
        # DEFINE THE BINS FOR COLUMN 'wavelength'
        lower = float("%0.*f" % (1, thisOrder["wavelength"].min()))
        upper = float("%0.*f" % (1, thisOrder["wavelength"].max()))
        wlBins = np.linspace(lower, upper, 1000)
        # BIN RESULTS
        thisOrder["wl_bins"] = pd.cut(thisOrder["wavelength"], bins=wlBins)
        wlGroups = thisOrder[["pixelScale", "wl_bins"]].groupby(["wl_bins"])
        thisOrder["pixelScale"] = wlGroups["pixelScale"].transform(lambda x: x.mean())

        minimum = 1000000
        bestShift = 0.0

        if order == self.qcPlotOrder or True:
            orderEdge = int(len(thisOrder.index) / 10)
            flux = thisOrder["flux"][orderEdge:-orderEdge]
            wl = thisOrder["wavelength"][orderEdge:-orderEdge]
            sp = thisOrder["slit_position"][orderEdge:-orderEdge]
            ps = thisOrder["pixelScale"][orderEdge:-orderEdge]
            scatter = thisOrder["residual_windowed_std"][orderEdge:-orderEdge]
            for shift in np.arange(-0.01, 0.01, 0.00001):
                # shiftedWl = wl - sp * ps * shift
                shiftedWl = wl - sp * shift
                sky = ip.splev(shiftedWl, tck)
                der = ip.splev(shiftedWl, tck, der=1)

                skySubtractedFlux = np.abs((flux - sky) * der / scatter)
                mean = np.mean(skySubtractedFlux)
                # median = np.median(skySubtractedFlux)
                # p80 = np.percentile(skySubtractedFlux, 80)
                # p90 = np.percentile(skySubtractedFlux, 90)
                # p95 = np.percentile(skySubtractedFlux, 95)
                # p99 = np.percentile(skySubtractedFlux, 99)

                if mean < minimum and mean > 0:
                    minimum = mean
                    bestShift = shift

            # if order == self.qcPlotOrder:
            #     print(f"{shift:0.5f}, {mean:0.3f}, {median:0.3f}, {p80:0.3f}, {p90:0.3f}, {p95:0.3f}, {p99:0.3f}")

        # print(f"BEST SHIFT: {bestShift:0.5f}. MIN: {minimum:0.2f}")

        orderDF["wavelength"] = orderDF["wavelength"] - orderDF["slit_position"] * bestShift

        # orderDF['wavelength'] = orderDF["wavelength"] - orderDF["slit_position"] * orderDF["pixelScale"] * bestShift

        orderDF.sort_values(by=["wavelength"], inplace=True)

        self.log.debug("completed the ``adjust_tilt`` method")
        return orderDF

    def determine_residual_floor(self, imageMapOrder, tck, iteration):
        """*determine residual floor and flag sky-lines*

        **Key Arguments:**
            - ``imageMapOrderDF`` --  dataframe with various processed data for a given order
            - ``tck`` -- the fitted bspline components. t for knots, c of coefficients, k for order
            - ``iteration`` -- the iteration number of the sky-subtraction loop, used to determine how aggressive to be in flagging sky-lines and determining the residual floor

        **Return:**

        - `imageMapOrder` -- same dataframe but now with sky-line locations flagged
        - `residualFloor` -- the residual floor determined within regions containing no skylines.

        **Usage:**

        ```python
        imageMapOrder, residualFloor = self.determine_residual_floor(imageMapOrder, tck)
        ```
        """
        self.log.debug("starting the ``determine_residual_floor`` method")

        import matplotlib.pyplot as plt
        import numpy as np
        import scipy.interpolate as ip
        from astropy.stats import sigma_clip
        from astropy.stats import sigma_clipped_stats

        order = imageMapOrder["order"].values[0]

        # USE THIS PERCENTILE TO DETERMINE THE RESIDUAL FLOOR
        residual_floor_percentile = self.recipeSettings["sky-subtraction"]["residual_floor_percentile"]
        # SIGNIFICANCE OF THE FIRST DERIVATIVE OF THE SKY (DO WE HAVE A LINE)
        line_significance = self.recipeSettings["sky-subtraction"]["skyline_significance"] / 10.0
        # SIGNIFICANCE OF THE SECOND DERIVATIVE OF THE SKY, i.e. A RISE OR FALL IN THE SKY SPECTRUM ALONG THE DISPERSION AXIS (LINE WINGS)
        node_significance = self.recipeSettings["sky-subtraction"]["skyline_significance"] / 10.0
        # SIGNIFICANCE OF THE LOCAL ROLLING WINDOW RESIDUAL VALUE COMPARED TO THE GLOBAL RESIDUALS - TO DETERMINE IF WE ARE JUST LOOKING AT NOISE
        noise_significance = self.recipeSettings["sky-subtraction"]["noise_sigma"]

        # DO SOME CLIPPING ON THE INITIAL SKY SUBTRACTION RESIDUALS
        imageMapOrder["sky_residuals"] = np.nan
        imageMapOrder["sky_residuals_clipped"] = False
        if "flagged_sky_line" not in imageMapOrder.columns:
            imageMapOrder["flagged_sky_line"] = False

        mask_unclipped = imageMapOrder["flagged_all_clipped"] == False
        imageMapOrder.loc[mask_unclipped, "sky_residuals"] = imageMapOrder.loc[
            mask_unclipped, "flux"
        ].values - ip.splev(imageMapOrder.loc[mask_unclipped, "wavelength"].values, tck)

        # FUDGE FOR NON-DARK SUBTRACTED DATA
        imageMapOrder.replace([np.inf, -np.inf], 1, inplace=True)

        ## DETERMINE THE RESIDUAL FLOOR WITHIN UNCLIPPED REGIONS - LONG WINDOW ROLLING QUANTILE
        window = int(20000 / iteration)
        if window < 1000:
            window = 1000

        # window = 5000

        # def sliding_median(arr, window):
        #         if window % 2 == 0:
        #             window += 1

        #         medianArry = np.median(np.lib.stride_tricks.sliding_window_view(arr, (window,)), axis=1)
        #         start = np.ones(window // 2) * medianArry[0]
        #         end = np.ones(window // 2) * medianArry[-1]
        #         medianArry = np.concatenate((start, medianArry, end))
        #         return medianArry

        #     ## USE ROLLING MEAN TO ESTIMATE THE LOCAL RESIDUALS, WHICH CAN BE USED TO ADD NEW KNOTS IN HIGH-RESIDUAL AREAS
        #     imageMapOrder.loc[~mask_all_clipped, "sky_residual_floor"] = sliding_median(
        #         imageMapOrder.loc[~mask_all_clipped, "sky_residuals"].values, window=window
        #     )

        imageMapOrder.loc[mask_unclipped, "sky_residual_floor"] = (
            imageMapOrder.loc[mask_unclipped, "sky_residuals"]
            .rolling(window=window, center=True, closed="both", min_periods=25)
            .quantile(residual_floor_percentile / 100.0)
        )

        # RESET NOISY FLAG
        imageMapOrder["flagged_noisy_region"] = False

        mean, median, std = sigma_clipped_stats(
            imageMapOrder.loc[(mask_unclipped), "residual_windowed_long_median"],
            sigma=2.0,
            stdfunc="mad_std",
            cenfunc="median",
            maxiters=10,
        )

        if std < 1.0:
            std = 1.0

        mask_noise_limit = imageMapOrder["residual_windowed_long_median"] > median + noise_significance * std

        ## FLAG SECTIONS OF NEGATIVE FLUX AS NOISE
        mean, median, std = sigma_clipped_stats(
            imageMapOrder.loc[(mask_unclipped), "flux_windowed_long_median"],
            sigma=2.0,
            stdfunc="std",
            cenfunc="mean",
            maxiters=5,
        )

        # if std < 1.0:
        #     std = 1.0

        mask_noise_limit = mask_noise_limit | (
            imageMapOrder["flux_windowed_long_median"] < median - noise_significance * std
        )
        imageMapOrder.loc[(mask_unclipped & mask_noise_limit), "flagged_noisy_region"] = True

        # MARK RISING SKYLINES
        mask_unclipped = imageMapOrder["flagged_all_clipped"] == False
        mask_unnoisy = imageMapOrder["flagged_noisy_region"] == False

        imageMapOrder.loc[(mask_unclipped), "sky_local_vs_global"] = (
            imageMapOrder.loc[(mask_unclipped), "sky_residual_floor_local"]
            - imageMapOrder.loc[(mask_unclipped), "sky_residual_floor"]
        )

        mask_skyline = imageMapOrder.loc[(mask_unclipped), "sky_local_vs_global"] > 0
        mask_negative = imageMapOrder.loc[(mask_unclipped), "sky_local_vs_global"] < 0
        imageMapOrder.loc[
            (mask_unclipped & mask_negative),
            "sky_local_vs_global",
        ] = 0

        imageMapOrder.loc[
            (mask_unclipped & mask_skyline),
            "flagged_sky_line",
        ] = "line"

        # # print(defaultPointsPerKnot)
        # imageMapOrder.loc[(mask_unclipped & mask_unnoisy), "sky_d1"] = (
        #     imageMapOrder.loc[(mask_unclipped & mask_unnoisy), "sky_d1"]
        #     .rolling(25, center=True, min_periods=3)
        #     .median()
        # )
        # # print(defaultPointsPerKnot)
        # imageMapOrder.loc[(mask_unclipped & mask_unnoisy), "sky_d2"] = (
        #     imageMapOrder.loc[(mask_unclipped & mask_unnoisy), "sky_d2"].rolling(9, center=True, min_periods=3).median()
        # )

        # ## FLAG SECTIONS OF NEGATIVE FLUX AS NOISE
        # mean, median, std = sigma_clipped_stats(
        #     imageMapOrder.loc[(mask_unclipped & mask_unnoisy), "sky_d1"],
        #     sigma=15,
        #     stdfunc="mad_std",
        #     cenfunc="median",
        #     maxiters=3,
        # )

        # imageMapOrder.loc[
        #     (mask_unclipped & mask_unnoisy)
        #     & (imageMapOrder["sky_d1"] > mean + line_significance * std)
        #     & (imageMapOrder["sky_d0"] > mean + line_significance * std),
        #     "flagged_sky_line",
        # ] = "rise"
        # imageMapOrder.loc[
        #     (mask_unclipped & mask_unnoisy) & (imageMapOrder["sky_d1"] < mean - line_significance * std),
        #     "flagged_sky_line",
        # ] = "fall"

        # MARK SKY INFECTION POINTS
        ## FLAG SECTIONS OF NEGATIVE FLUX AS NOISE
        # mask_skyline = imageMapOrder["flagged_sky_line"].isin(["rise", "fall", "line"])
        # if iteration > 4:
        #     mean, median, std = sigma_clipped_stats(
        #         imageMapOrder.loc[(mask_unclipped & mask_unnoisy), "sky_d2"],
        #         sigma=100,
        #         stdfunc="mad_std",
        #         cenfunc="median",
        #         maxiters=1,
        #     )

        #     imageMapOrder.loc[
        #         (mask_unclipped & mask_unnoisy) & (imageMapOrder["sky_d2"] > mean + node_significance * std),
        #         "flagged_sky_line",
        #     ] = "node"
        #     imageMapOrder.loc[
        #         (mask_unclipped & mask_unnoisy) & (imageMapOrder["sky_d2"] < mean - node_significance * std),
        #         "flagged_sky_line",
        #     ] = "node"

        # CALCULATE THE RESIDUAL FLOOR
        mask_nonSkyline = imageMapOrder["flagged_sky_line"] == False

        # DETERMINE SKYLINES PEAKS
        skyFlux = imageMapOrder.loc[(mask_unclipped & mask_unnoisy), "flux"].values
        skyFluxMean, skyFluxStd = np.median(skyFlux), skyFlux.std()
        highFluxMask = imageMapOrder["flux"] > skyFluxMean + 1.5 * skyFluxStd
        imageMapOrder.loc[(mask_unclipped & mask_unnoisy & mask_nonSkyline & highFluxMask), "flagged_sky_line"] = "peak"

        self.log.debug("completed the ``determine_residual_floor`` method")
        # print(residualFloor)
        return imageMapOrder, 5

    def plot_order_skymodel_fitting_quicklook(self, imageMapOrder, tck, title=None, knots=False):
        """Quick-look diagnostic plot of the sky-model fit for a single order."""
        from soxspipe.commonutils.toolkit import get_calibrations_path
        from astropy.table import Table
        import matplotlib.pyplot as plt
        import numpy as np
        import scipy.interpolate as ip
        from astropy.stats import sigma_clipped_stats

        # GET THE SKYLINES TO PLOT FOR REFERENCE
        dp = detector_lookup(log=self.log, settings=self.settings).get(self.arm)
        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)
        skylines = calibrationRootPath + "/" + dp["skylines"]

        dat = Table.read(skylines, format="fits")
        skylinesDF = dat.to_pandas()
        # FILTER TO STRONG SKY LINES ONLY FOR PLOTTING
        skylinesDF["WAVELENGTH"] = skylinesDF["WAVELENGTH"].astype(float)
        skylinesDF["FLUX"] = skylinesDF["FLUX"].astype(float)
        if self.arm == "VIS":
            mask = skylinesDF["FLUX"] > 5
        else:
            mask = skylinesDF["FLUX"] > 100
        skylinesDF = skylinesDF.loc[mask]

        # GENERATE THE FIGURE FOR THE PLOT
        from matplotlib.gridspec import GridSpec

        fig = plt.figure(num=None, figsize=(24, 10), dpi=150, facecolor=None, edgecolor=None, frameon=True)
        gs = GridSpec(2, 1, figure=fig, height_ratios=[3, 1], hspace=0.4)
        ax = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax.set_xlabel("wavelength (nm)")
        ax.set_ylabel("Counts")
        ax2.set_title("sky subtraction residuals", fontsize=12)
        ax2.set_xlabel("wavelength (nm)")
        ax2.set_ylabel("Flux")
        if title:
            ax.set_title(title, fontsize=14)

        def refresh_and_plot(pause=None):
            for axis in [ax, ax2]:
                handles, labels = axis.get_legend_handles_labels()
                unique_handles = []
                unique_labels = []

                for handle, label in zip(handles, labels):
                    if not label or label.startswith("_") or label in unique_labels:
                        continue

                    unique_handles.append(handle)
                    unique_labels.append(label)

                existing_legend = getattr(axis, "legend_", None)
                if existing_legend is not None:
                    existing_legend.remove()

                if unique_handles:
                    lgnd = axis.legend(
                        unique_handles,
                        unique_labels,
                        loc="upper right",
                        bbox_to_anchor=(1.2, 1.0),
                        prop={"size": 8},
                    )
                    for legobj in lgnd.legend_handles:
                        try:
                            legobj.set_sizes([30])
                            legobj.set_alpha(0.7)
                        except:
                            pass
            fig.canvas.draw()
            fig.canvas.flush_events()
            if pause:
                plt.pause(pause)
            else:
                plot.show()

        # BUILD THE MASKS FOR PLOTTING VARIOUS DATA POINTS
        mask_everything = imageMapOrder["wavelength"] > 0.0
        mask_object = imageMapOrder["flagged_object_clipped"] == True
        mask_all_clipped = (imageMapOrder["flagged_all_clipped"] == True) | (
            imageMapOrder["flagged_object_clipped"] == True
        )
        mask_bad_pixel = imageMapOrder["flagged_bad_pixel_clipped"] == True
        mask_order_edge = imageMapOrder["flagged_edge_clipped"] == True

        # LINES
        wl_not_clipped = imageMapOrder.loc[~mask_all_clipped, "wavelength"].values
        flux_not_clipped = imageMapOrder.loc[~mask_all_clipped, "flux"].values
        flux_percentile_smoothed = imageMapOrder.loc[~mask_all_clipped, "flux_percentile_smoothed"].values
        flux_minus_smoothed_residuals = imageMapOrder.loc[~mask_all_clipped, "flux_minus_smoothed_residual"].values
        flux_minus_smoothed_residual_upper_limit = imageMapOrder.loc[
            ~mask_all_clipped, "flux_minus_smoothed_residual_upper_limit"
        ].values
        flux_upper_limit = imageMapOrder.loc[~mask_all_clipped, "flux_upper_limit"].values

        if tck:
            mask_noisy = imageMapOrder["flagged_noisy_region"] == True
            mask_nonSkyline = imageMapOrder["flagged_sky_line"] == False
            mask_skyline = imageMapOrder["flagged_sky_line"].isin(["rise", "fall", "line"])
            mask_skynode = imageMapOrder["flagged_sky_line"].isin(["node"])
            mask_bspline_clipped = imageMapOrder["flagged_bspline_clipped"] == True
            flux_bspline_sky_residuals = imageMapOrder.loc[~mask_all_clipped, "sky_residuals"].values
            # sky_flux_rolling_median = imageMapOrder.loc[~mask_all_clipped, "sky_flux_rolling_median"].values
            sky_residual_floor = imageMapOrder.loc[~mask_all_clipped, "sky_residual_floor"].values
            sky_residual_floor_local = imageMapOrder.loc[~mask_all_clipped, "sky_residual_floor_local"].values

        # SET PLOT LIMITS
        mean, median, std = sigma_clipped_stats(
            imageMapOrder.loc[~mask_all_clipped & imageMapOrder["flux"] > -50, "flux"].values, sigma=3.0, maxiters=3
        )
        range_sigma = 3
        ax.set_ylim(mean - std, mean + range_sigma * std)

        # PLOT SKY LINES AS VERTICAL LINES ON SKY PANEL
        label = "catalogue skyline"
        for _, row in skylinesDF.iterrows():
            ax.axvline(row["WAVELENGTH"], color="gray", linestyle="-", linewidth=0.5, alpha=0.5, zorder=0, label=label)
            ax2.axvline(row["WAVELENGTH"], color="gray", linestyle="-", linewidth=0.5, alpha=0.5, zorder=0, label=label)
            if label:
                label = None

        mmin = imageMapOrder.loc[(~mask_object), "wavelength"].values.min()
        mmax = imageMapOrder.loc[(~mask_object), "wavelength"].values.max()
        ax.set_xlim(mmin, mmax)
        ax2.set_xlim(mmin, mmax)

        # Put a legend on plot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])
        box2 = ax2.get_position()
        ax2.set_position([box2.x0, box2.y0, box2.width * 0.95, box2.height])

        # FLUX - SMOOTHED RESIDUALS
        if True and not tck:
            mean, median, std = sigma_clipped_stats(
                imageMapOrder.loc[
                    ~mask_all_clipped & imageMapOrder["flux_minus_smoothed_residual"] > -50,
                    "flux_minus_smoothed_residual",
                ].values,
                sigma=3.0,
                maxiters=3,
            )
            range_sigma = 3
            ax2.set_ylim(mean - std, mean + range_sigma * std)
            ax2.scatter(
                wl_not_clipped,
                flux_minus_smoothed_residuals,
                s=1,
                c="grey",
                alpha=0.2,
                zorder=1,
                label="flux residuals (windowed)",
            )

            if False:
                ax2.scatter(
                    wl_not_clipped,
                    flux_minus_smoothed_residual_upper_limit,
                    s=1,
                    c="red",
                    alpha=0.01,
                    zorder=1,
                    label="flux residuals upper limit",
                )
        elif True:
            mean, median, std = sigma_clipped_stats(
                imageMapOrder.loc[
                    ~mask_all_clipped & imageMapOrder["sky_residuals"] > -50,
                    "sky_residuals",
                ].values,
                sigma=3.0,
                maxiters=3,
            )
            range_sigma = 7
            ax2.set_ylim(mean - std, mean + range_sigma * std)
            ax2.scatter(
                wl_not_clipped,
                flux_bspline_sky_residuals,
                s=1,
                c="grey",
                alpha=0.2,
                zorder=1,
                label="flux residuals (windowed)",
            )
            ax2.plot(
                wl_not_clipped,
                sky_residual_floor,
                c="blue",
                alpha=0.2,
                zorder=1,
                label="sky residual floor",
            )
            ax2.plot(
                wl_not_clipped,
                sky_residual_floor_local,
                c="green",
                alpha=0.2,
                zorder=1,
                label="sky residual floor local",
            )
            if not isinstance(knots, bool):
                ax2.scatter(
                    knots,
                    np.zeros_like(knots),
                    marker="v",
                    s=1,
                    c="red",
                    alpha=0.5,
                    zorder=10,
                    label="bspline knots",
                )

        ## UNCLIPPED
        if True:
            ax.scatter(
                wl_not_clipped,
                flux_not_clipped,
                s=1,
                c="green",
                alpha=0.1,
                zorder=1,
                label="unclipped",
            )
        if False:
            ax.scatter(
                wl_not_clipped,
                flux_upper_limit,
                s=1,
                c="red",
                alpha=0.01,
                zorder=1,
                label="flux upper limit (windowed)",
            )
        if True and not tck:
            # SHOW ROLLING SCATTER
            ax.scatter(
                wl_not_clipped,
                flux_percentile_smoothed,
                s=1,
                c="blue",
                alpha=0.01,
                zorder=5,
                label="flux percentile smoothed (windowed)",
            )

        ## EDGES
        if False and not tck:
            ax.scatter(
                imageMapOrder.loc[mask_order_edge, "wavelength"].values,
                imageMapOrder.loc[mask_order_edge, "flux"].values,
                s=1,
                marker="x",
                c="green",
                alpha=0.1,
                zorder=2,
                label="edge clipped",
            )

        ## BAD PIXELS
        if True and not tck:
            ax.scatter(
                imageMapOrder.loc[mask_bad_pixel, "wavelength"].values,
                imageMapOrder.loc[mask_bad_pixel, "flux"].values,
                s=1,
                marker="v",
                c="red",
                alpha=0.1,
                zorder=2,
                label="bad pixels clipped",
            )
        # refresh_and_plot(0.2)

        ## OBJECT
        if True and not tck:
            ax.scatter(
                imageMapOrder.loc[mask_object, "wavelength"].values,
                imageMapOrder.loc[mask_object, "flux"].values,
                s=1,
                c="cyan",
                alpha=0.1,
                zorder=2,
                label="object",
            )
            # refresh_and_plot(0.2)

        if not tck:
            plt.pause(1)
            plt.close("all")  # close the figure
            return

        ## SHOT-NOISE
        if False:
            ax.scatter(
                imageMapOrder.loc[(~mask_all_clipped & mask_nonSkyline), "wavelength"].values,
                imageMapOrder.loc[(~mask_all_clipped & mask_nonSkyline), "flux"].values,
                s=1,
                c="green",
                alpha=0.1,
                zorder=1,
                label="shot-noise",
            )
            refresh_and_plot(0.2)

        # PREDICTED SKYLINES - ORANGE
        if True:
            ax.scatter(
                imageMapOrder.loc[(~mask_all_clipped & mask_skyline), "wavelength"].values,
                imageMapOrder.loc[(~mask_all_clipped & mask_skyline), "flux"].values,
                s=1,
                c="orange",
                alpha=0.3,
                zorder=2,
                label="predicted skyline",
            )
            # refresh_and_plot(0.2)

        # PREDICTED SKYLINE WINGS - RED
        if True:
            ax.scatter(
                imageMapOrder.loc[(~mask_all_clipped & mask_skynode), "wavelength"].values,
                imageMapOrder.loc[(~mask_all_clipped & mask_skynode), "flux"].values,
                s=1,
                c="red",
                alpha=0.3,
                zorder=3,
                label="predicted skyline wings/peak",
            )
            # refresh_and_plot(0.2)

        # BSPLINE CLIPPED
        if True:
            # SHOW BSPLINE CLIPPED FROM PREVIOUS ITERATION
            ax.scatter(
                imageMapOrder.loc[(mask_bspline_clipped), "wavelength"].values,
                imageMapOrder.loc[(mask_bspline_clipped), "flux"].values,
                s=1,
                c="blue",
                alpha=1,
                zorder=10,
                label="bspline clipped",
            )
            # refresh_and_plot(0.2)

        # NOISY DATA REGION
        if True:

            ax.scatter(
                imageMapOrder.loc[(mask_noisy), "wavelength"].values,
                imageMapOrder.loc[(mask_noisy), "flux"].values,
                s=2,
                marker="x",
                c="grey",
                alpha=0.3,
                zorder=4,
                label="noisy region",
            )
            # refresh_and_plot(0.2)

        # PLOT SKY RESIDUALS
        if False:
            ax.scatter(
                imageMapOrder.loc[mask_all_clipped, "wavelength"].values,
                imageMapOrder.loc[mask_all_clipped, "sky_residuals"].values,
                s=1,
                c="green",
                alpha=1,
                zorder=13,
                label="sky residuals",
            )
            refresh_and_plot(0.2)

        # CLIPPED RESIDUALS
        if False:
            ax.scatter(
                imageMapOrder.loc[(imageMapOrder["sky_residuals_clipped"] == True), "wavelength"].values,
                imageMapOrder.loc[(imageMapOrder["sky_residuals_clipped"] == True), "sky_residuals"].values,
                s=5,
                c="red",
                marker="x",
                alpha=1,
                zorder=13,
                label="clipped residuals",
            )
            refresh_and_plot(0.2)

        wl = np.linspace(imageMapOrder["wavelength"].min(), imageMapOrder["wavelength"].max(), 1000000)
        sky = ip.splev(wl, tck)
        # sky_d1 = ip.splev(wl, tck, der=1)
        # sky_d2 = ip.splev(wl, tck, der=2)
        # sky_d3 = ip.splev(wl, tck, der=3)
        # ax.plot(wl, sky, label="sky model", c="#268bd2", zorder=10, alpha=0.7)

        ## BSPLINE SKY
        if True:
            ax.plot(
                wl,
                sky,
                label="sky model bspline",
                c="purple",
                zorder=3,
                alpha=0.5,
            )
            # refresh_and_plot(0.2)
            # ax.plot(wl, sky_d2 / 100, label="sky model 2nd derivative / 100", c="violet", zorder=3, alpha=0.2)

        # ROLLING MEDIAN OF THE SKY SUBTRACTED FLUX
        if False:
            ax.plot(
                wl_not_clipped,
                sky_flux_rolling_median,
                label="sky flux_rolling_median",
                c="purple",
                zorder=3,
                alpha=0.5,
            )
            # ax.plot(wl, sky_d2 / 100, label="sky model 2nd derivative / 100", c="violet", zorder=3, alpha=0.2)

        if False:
            ax.plot(
                wl,
                sky_d1,
                label="sky model derivative",
                c="pink",
                zorder=5,
                alpha=0.5,
            )
            refresh_and_plot(0.2)
            # ax.plot(wl, sky_d2 / 100, label="sky model 2nd derivative / 100", c="violet", zorder=3, alpha=0.2)

        if False:
            ax.plot(
                wl,
                sky_d3,
                label="sky model derivative 3",
                c="orange",
                zorder=5,
                alpha=0.5,
            )
            refresh_and_plot(0.2)

        if False:
            # SHOW ROLLING SCATTER
            ax.plot(
                imageMapOrder.loc[(mask_all_clipped), "wavelength"].values,
                imageMapOrder.loc[(mask_all_clipped), "residual_windowed_long_median"].values - 2 * std,
                c="black",
                alpha=0.3,
                zorder=0,
                label="residual windowed std",
            )
            refresh_and_plot(0.2)

        plt.pause(0.1)
        plt.show()
        plt.close("all")  # close the figure

    # use the tab-trigger below for new method
    # xt-class-method
