#!/usr/bin/env python
# encoding: utf-8
"""
*Subtract the sky background using the Kelson Method*

:Author:
    David Young

:Date Created:
    April 14, 2022
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
os.environ['TERM'] = 'vt100'


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

    **Usage:**

    To setup your logger and settings, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a `subtract_sky` object, use the following:

    ```eval_rst
    .. todo::

        - add a tutorial about ``subtract_sky`` to documentation
    ```

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
            recipeName="soxs-stare"
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

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = objectFrame.header[kw("SEQ_ARM")]
        self.inst = objectFrame.header[kw("INSTRUME")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)
        dp = self.detectorParams

        # UNPACK THE 2D DISP IMAGE MAP AND THE OBJECT IMAGE TO GIVE A
        # DATA FRAME CONTAINING ONE ROW FOR EACH PIXEL WITH COLUMNS X, Y, FLUX, WAVELENGTH, SLIT-POSITION, ORDER
        self.mapDF, self.interOrderMask = twoD_disp_map_image_to_dataframe(log=self.log, slit_length=dp["slit_length"], twoDMapPath=twoDMap, associatedFrame=self.objectFrame, kw=kw)

        # DETERMINE SLIT
        self.slit = objectFrame.header[kw(f"SLIT_{self.arm}".upper())]
        # ACCOUNT FOR BLOCKING FILTER
        if "JH" in self.slit:
            self.mapDF = self.mapDF.loc[(self.mapDF["order"] > 12)]

        quicklook_image(
            log=self.log, CCDObject=self.objectFrame, show=False, ext=False, stdWindow=0.1, title="science frame awaiting sky-subtraction", surfacePlot=False, dispMap=dispMap, dispMapImage=twoDMap, settings=self.settings, skylines=True)

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
            self.filenameTemplate = filenamer(
                log=self.log,
                frame=self.objectFrame,
                settings=self.settings
            )

        home = expanduser("~")
        self.qcDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}/"
        self.qcDir = self.qcDir.replace("//", "/")
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.qcDir):
            os.makedirs(self.qcDir)

        if self.arm != "NIR" and kw('WIN_BINX') in self.objectFrame.header:
            self.binx = int(self.objectFrame.header[kw('WIN_BINX')])
            self.biny = int(self.objectFrame.header[kw('WIN_BINY')])
        else:
            self.binx = 1
            self.biny = 1
        return

    def subtract(self):
        """
        *generate and subtract a sky-model from the input frame*

        **Return:**
            - ``skymodelCCDData`` -- CCDData object containing the model sky frame
            - ``skySubtractedCCDData`` -- CCDData object containing the sky-subtacted frame
            - ``qcTable`` -- the data frame containing measured QC metrics
            - ``productsTable`` -- the data frame containing collected output products
        """
        self.log.debug('starting the ``get`` method')

        import numpy as np
        import pandas as pd
        pd.options.mode.chained_assignment = None

        self.log.print(f'\n# MODELLING SKY BACKGROUND AND REMOVING FROM SCIENCE FRAME')

        # THESE PLACEHOLDERS ARE INITAILLY BLANK AND AWAITING PIXEL VALUES TO BE ADDED
        skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData = self.create_placeholder_images()

        uniqueOrders = self.mapDF['order'].unique()

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # THE BSPLINE ORDER TO FIT SKY WITH
        bspline_order = self.recipeSettings["sky-subtraction"]["bspline_order"]

        # SELECT A SINGLE ORDER TO GENERATE QC PLOTS FOR
        qcPlotOrder = int(np.median(uniqueOrders)) - 1
        print("FIX ME")
        qcPlotOrder = 12

        allimageMapOrder = []
        allimageMapOrderWithObject = []

        # GET OVER SAMPLED SKY & SKY+OBJECT AS LISTS OF DATAFRAMES
        self.log.print(f"\n  ## CLIPPING DEVIANT PIXELS AND PIXELS WITH OBJECT FLUX\n")

        imageMapOrders = []
        for o in uniqueOrders:
            # SELECT DATAFRAME CONTAINING ONLY A SINGLE ORDER
            imageMapOrders.append(self.mapDF[self.mapDF["order"] == o])

        # NOTE MULTIPROCESSING THIS BLOCK RESULTS IN SLOWER PERFORMANCE
        for o in uniqueOrders:
            # SELECT ONLY A DATAFRAME CONTAINING ONLY A SINGLE ORDER
            imageMapOrder = self.mapDF[self.mapDF["order"] == o]

            # MASK OUTLYING PIXELS (imageMapOrderWithObject) AND ALSO THEN THE OBJECT PIXELS (imageMapOrderSkyOnly)
            imageMapOrder = self.get_over_sampled_sky_from_order(imageMapOrder, clipBPs=True, clipSlitEdge=self.recipeSettings["sky-subtraction"]["clip-slit-edge-fraction"])
            allimageMapOrder.append(imageMapOrder)

        # MASK OUT OBJECT PIXELS
        allimageMapOrder = self.clip_object_slit_positions(allimageMapOrder, aggressive=self.recipeSettings["sky-subtraction"]["aggressive_object_masking"])

        self.log.print(f"\n  ## FITTING SKY-FLUX WITH A BSPLINE (WAVELENGTH) AND LOW-ORDER POLY (SLIT-ILLUMINATION PROFILE)\n")

        # NOTE MULTIPROCESSING THIS BLOCK RESULTS IN SLOWER PERFORMANCE
        newAllimageMapOrder = []
        allFluxErrorRatios = []
        allResidualFloor = []
        alltck = []
        allKnots = []
        totalKnots = 0
        for o, imageMapOrder in zip(uniqueOrders, allimageMapOrder):

            imageMapOrder, tck, knots, flux_error_ratio, residualFloor = self.fit_bspline_curve_to_sky(imageMapOrder, bspline_order)
            totalKnots += len(knots)
            newAllimageMapOrder.append(imageMapOrder)
            allFluxErrorRatios.append(flux_error_ratio)
            allResidualFloor.append(residualFloor)
            alltck.append(tck)
            allKnots.append(knots)
        allimageMapOrder = newAllimageMapOrder

        allFluxErrorRatios = np.concatenate(allFluxErrorRatios)
        self.log.print(f'\n\tFULL FRAME SKY-MODEL FLUX TO ERROR METRICS: MEAN {allFluxErrorRatios.mean():0.3f}, STD {allFluxErrorRatios.std():0.3f}, MEDIAN {np.median(allFluxErrorRatios):0.3f}, MAX {allFluxErrorRatios.max():0.3f}, MIN {allFluxErrorRatios.min():0.3f}, RANGE {allFluxErrorRatios.max()-allFluxErrorRatios.min():0.3f}, MEAN RES FLOOR: {np.mean(allResidualFloor):0.3f}, KNOT COUNT: {totalKnots}')
        self.log.print(f'\t{allFluxErrorRatios.mean():0.3f},  {allFluxErrorRatios.std():0.3f},  {np.median(allFluxErrorRatios):0.3f},  {allFluxErrorRatios.max():0.3f},  {allFluxErrorRatios.min():0.3f}, {allFluxErrorRatios.max()-allFluxErrorRatios.min():0.3f},{np.mean(allResidualFloor):0.3f},{totalKnots}')

        for o, imageMapOrder, tck, knots in zip(uniqueOrders, allimageMapOrder, alltck, allKnots):
            if isinstance(imageMapOrder, pd.core.frame.DataFrame):
                # INJECT THE PIXEL VALUES BACK INTO THE PLACEHOLDER IMAGES
                skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData = self.add_data_to_placeholder_images(imageMapOrder, skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData)
                if o == qcPlotOrder and True:
                    qc_plot_path = self.plot_sky_sampling(order=o, imageMapOrderDF=imageMapOrder, tck=tck, knotLocations=knots)
                    basename = os.path.basename(qc_plot_path)
                    self.products = pd.concat([self.products, pd.Series({
                        "soxspipe_recipe": "soxs-stare",
                        "product_label": "SKY_MODEL_QC_PLOTS",
                        "file_name": basename,
                        "file_type": "PDF",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "product_desc": f"QC plots for the sky-background modelling",
                        "file_path": qc_plot_path,
                        "label": "QC"
                    }).to_frame().T], ignore_index=True)

        # SET NANs TO 0
        skymodelCCDData.data[np.isnan(skymodelCCDData.data)] = 0
        skymodelCCDData.uncertainty.array[np.isnan(skymodelCCDData.uncertainty.array)] = 0
        skySubtractedCCDData.data[np.isnan(skySubtractedCCDData.data)] = 0
        skySubtractedCCDData.uncertainty.array[np.isnan(skySubtractedCCDData.uncertainty.array)] = 0

        comparisonPdf = self.plot_image_comparison(self.objectFrame, skymodelCCDData, skySubtractedCCDData)

        filename = os.path.basename(comparisonPdf)
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "SKY SUBTRACTION QUICKLOOK",
            "file_name": filename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"Sky-subtraction quicklook",
            "file_path": comparisonPdf,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``get`` method')
        return skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData, self.qc, self.products

    def get_over_sampled_sky_from_order(
            self,
            imageMapOrder,
            clipBPs=True,
            clipSlitEdge=False):
        """*unpack the over sampled sky from an order*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe from object image and 2D map
            - ``clipBPs`` -- clip bad-pixels? Deafult *True*
            - ``clipSlitEdge`` -- clip the slit edges. Percentage of slit width to clip. Default *False*

        **Return:**
            - `imageMapOrder` -- input order dataframe with outlying pixels masked AND object pixels masked

        **Usage:**

        ```python
        imageMapOrderWithObject, imageMapOrderSkyOnly = skymodel.get_over_sampled_sky_from_order(imageMapOrder, o, ignoreBP=False, clipSlitEdge=0.00)
        ```
        """
        self.log.debug('starting the ``get_over_sampled_sky_from_order`` method')

        from astropy.stats import sigma_clip, mad_std

        # COLLECT SETTINGS
        percentile_clipping_sigma = self.recipeSettings["sky-subtraction"]["percentile_clipping_sigma"]
        percentile_clipping_iterations = self.recipeSettings["sky-subtraction"]["percentile_clipping_iterations"]
        percentile_rolling_window_size = self.recipeSettings["sky-subtraction"]["percentile_rolling_window_size"]
        self.rollingWindowSize = int(percentile_rolling_window_size)

        imageMapOrder["clipped"] = False
        imageMapOrder["object"] = False
        imageMapOrder["bspline_clipped"] = False
        imageMapOrder["edge_clipped"] = False
        imageMapOrder["bad_pixel_clipped"] = False

        # ASSIGN ORDER-EDGE PIXELS A 'clipped' FLAG
        if clipSlitEdge:
            slitRange = imageMapOrder["slit_position"].max() - imageMapOrder["slit_position"].min()
            clipSlitEdge *= slitRange
            mask = ((imageMapOrder['slit_position'] > imageMapOrder["slit_position"].max() - clipSlitEdge) | (imageMapOrder['slit_position'] < imageMapOrder["slit_position"].min() + clipSlitEdge))
            imageMapOrder.loc[mask, "edge_clipped"] = True
            imageMapOrder.loc[mask, "clipped"] = True

        # ASSIGN BAD-PIXELS A 'clipped' FLAG?
        if clipBPs:
            mask = (imageMapOrder['mask'] == True)
            imageMapOrder.loc[mask, "bad_pixel_clipped"] = True
            imageMapOrder.loc[mask, "clipped"] = True

        # CLIP THE MOST DEVIANT PIXELS WITHIN A WAVELENGTH ROLLING WINDOW - BAD-PIXELS AND CRHs AND OBJECTS
        imageMapOrder = self.rolling_window_clipping(imageMapOrderDF=imageMapOrder, windowSize=self.rollingWindowSize, sigma_clip_limit=percentile_clipping_sigma, max_iterations=percentile_clipping_iterations)

        self.log.debug('completed the ``get_over_sampled_sky_from_order`` method')
        return imageMapOrder

    def plot_sky_sampling(
            self,
            order,
            imageMapOrderDF,
            tck=False,
            knotLocations=False):
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
        self.log.debug('starting the ``plot_sky_sampling`` method')

        import numpy as np
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt
        import scipy.interpolate as ip
        import numpy.ma as ma
        from matplotlib import cm
        from matplotlib import colors
        from copy import copy

        # SET COLOURS FOR VARIOUS STAGES
        red = "#dc322f"
        blue = "#268bd2"
        black = "#002b36"
        grey = "#93a1a1"
        green = "green"
        orange = "#cb4b16"
        violet = "#6c71c4"

        # MAKE A COPY OF THE FRAME TO NOT ALTER ORIGINAL DATA
        frame = self.objectFrame.copy()

        # SETUP THE PLOT SUB-PANELS
        print("FIX ME")
        if False:
            fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=320)
        else:
            # REMOVE ME
            fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=100)

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
            nonOrderMask[y][x] = 0

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        nonOrderMask = ma.make_mask(nonOrderMask)
        combinedMask = (nonOrderMask == 1) | (frame.mask == 1)
        frame.mask = (nonOrderMask == 1)

        # RAW IMAGE PANEL
        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.flipud(np.rot90(frame, 1))
        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        maskedDataArray = np.ma.array(frame.data, mask=combinedMask)
        std = np.nanstd(maskedDataArray)
        mean = np.nanmean(maskedDataArray)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        im = onerow.imshow(rotatedImg, vmin=0, vmax=100, cmap='gray', alpha=1)
        medianValue = np.median(rotatedImg.data.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="unprocessed frame")]

        onerow.set_title("Object & Sky Frame", fontsize=10)
        onerow.set_xlabel(
            "y-axis", fontsize=10)
        onerow.set_ylabel(
            "x-axis", fontsize=10)
        ylimMinImage = imageMapOrderDF[self.axisB].min() - 10
        ylimMaxImage = imageMapOrderDF[self.axisB].max() + 10
        onerow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        onerow.set_xlim(ylimMinImage, ylimMaxImage)

        # ORIGINAL DATA AND PERCENTILE SMOOTHED WAVELENGTH VS FLUX
        tworow.plot(
            imageMapOrderDF["wavelength"].values,
            imageMapOrderDF["flux"].values, label='unprocessed', alpha=0.2, c=grey, zorder=0)
        tworow.set_title("STEP 1. Identify and clip outlying pixels (CRHs etc) and pixels containing object flux.", fontsize=10)

        # RAW MARKERS
        tworow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].values, label='unclipped', s=0.5, c=black, alpha=1., zorder=1)

        columnName = ["object", "bad_pixel_clipped", "edge_clipped", "bspline_clipped"]
        colours = [blue, orange, green, violet]
        alphas = [.2, .3, .8, .8]
        labels = ["object", "bad pixels", "order edges", "bspline clipped"]

        for cn, cl, lb, al in zip(columnName, colours, labels, alphas):
            tworow.scatter(
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "wavelength"].values,
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "flux"].values, label=lb, s=8, marker="x", c=cl, zorder=3, alpha=al)

        # PERCENTILE LINE
        tworow.plot(
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False), "wavelength"].values,
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False), "flux_smoothed"].values, label='percentile-smoothed', c=blue, zorder=3)

        # SIGMA RESIDUAL
        weights = tworow.plot(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"].values - imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"].max() * 1.2, label='$\sigma$ residual scatter (shifted)', c=black)
        ylimmin = -imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"].max() * 1.3

        if ylimmin < -3000:
            ylimmin = -300

        from astropy.stats import sigma_clip, mad_std
        # SIGMA-CLIP THE DATA
        masked = sigma_clip(
            imageMapOrderDF["flux"], sigma_lower=30, sigma_upper=30, maxiters=1, cenfunc='median', stdfunc=mad_std)

        tworow.set_ylim(ylimmin, masked.max())

        tworow.set_ylabel(
            "flux ($e^{-}$)", fontsize=10)
        tworow.set_xlabel(
            "wavelength", fontsize=10)
        tworow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
        # tworow.set_xticks([], [])

        # SLIT-POSITION RESIDUAL PANEL (SHOWING OBJECT)
        std = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "residual_global_sigma"].std()
        median = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "residual_global_sigma"].median()

        threerow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "slit_position"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_global_sigma"].values, label='deviations', s=0.5, alpha=0.5, c=grey, zorder=1)

        columnName = ["object", "bad_pixel_clipped", "bspline_clipped"]
        colours = [blue, orange, violet]
        alphas = [.2, 1, .4]
        labels = ["object", "bad pixels", "bspline clipped"]

        for cn, cl, lb, al in zip(columnName, colours, labels, alphas):
            threerow.scatter(
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "slit_position"].values,
                imageMapOrderDF.loc[imageMapOrderDF[cn] == True, "residual_global_sigma"].values, label=lb, s=3, marker="x", c=cl, zorder=3, alpha=al)

        threerow.set_ylim(median - 3 * std, median + 7 * std)
        threerow.set_xlabel(
            "slit-position relative to slit centre (arcsec)", fontsize=10)
        threerow.set_ylabel("flux minus smoothed flux residual ($\sigma$)", fontsize=10)

        threerow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        # IMAGE SHOWING CLIPPED PIXEL MASK
        im = fourrow.imshow(rotatedImg, vmin=0, vmax=100, cmap='gray', alpha=1)

        columnName = ["object", "bad_pixel_clipped", "edge_clipped", "bspline_clipped"]
        colours = [blue, orange, green, violet]
        alphas = [1, 1, 1, 1]
        labels = ["object", "bad pixels", "order edges", "bspline clipped"]
        patches = []

        for cn, cl, lb, al in zip(columnName, colours, labels, alphas):
            clippedMask = nonOrderMask
            clippedMask = np.zeros_like(frame.data)
            for x, y in zip(imageMapOrderDF.loc[imageMapOrderDF[cn] == True, self.axisA].values, imageMapOrderDF.loc[imageMapOrderDF[cn] == True, self.axisB].values):
                clippedMask[y][x] = 1
            clippedMask = ma.make_mask(clippedMask)
            imageMask = np.ma.array(np.ones_like(frame.data), mask=~clippedMask)
            # MAKE A COLOR MAP OF FIXED COLORS
            cmap = colors.ListedColormap([cl, cl])
            bounds = [0, 5, 10]
            norm = colors.BoundaryNorm(bounds, cmap.N)
            cmap.set_bad(cl, 0.)
            fourrow.imshow(np.flipud(np.rot90(imageMask, 1)), cmap=cmap, norm=norm, alpha=al, interpolation='nearest')
            patches.append(mpatches.Patch(color=cl, label=lb))

        fourrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        nonOrderMask = (nonOrderMask == 0)
        imageMask = np.ma.array(np.ones_like(frame.data), mask=nonOrderMask)
        cmap = copy(cm.gray)
        cmap.set_bad("green", 0.0)
        fourrow.imshow(np.flipud(np.rot90(imageMask, 1)), vmin=-10, vmax=-9, cmap=cmap, alpha=1.)
        fourrow.set_xlabel(
            "y-axis", fontsize=10)
        fourrow.set_ylabel(
            "x-axis", fontsize=10)
        fourrow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        fourrow.set_xlim(ylimMinImage, ylimMaxImage)
        # fourrow.invert_xaxis()

        # PLOT WAVELENGTH VS FLUX SKY MODEL
        fiverow.set_title("STEP 2. Fit a univariate bspline to sky-flux as a function of wavelength", fontsize=10)
        fiverow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].values, s=3, c=black, alpha=0.2, zorder=1, label="unclipped")
        fiverow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["bspline_clipped"] == True, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["bspline_clipped"] == True, "flux"].values, s=8, c=violet, marker="x", alpha=1., zorder=1, label="clipped")

        print(imageMapOrderDF.loc[imageMapOrderDF["bspline_clipped"] == True, "wavelength"].values.shape)

        if tck:
            wl = np.linspace(imageMapOrderDF["wavelength"].min(), imageMapOrderDF["wavelength"].max(), 1000000)
            sky = ip.splev(wl, tck)
            knotSky = ip.splev(knotLocations, tck)
            skymodel = fiverow.plot(
                wl, sky, label='sky model', c=blue, zorder=3)
            fiverow.scatter(knotLocations, knotSky, marker=7, s=15, alpha=0.7, c=red, zorder=3, label='knots')
        else:
            skymodel = fiverow.plot(
                imageMapOrderDF["wavelength"].values,
                imageMapOrderDF["sky_model"].values, label='sky model', c=blue, zorder=3)
        if ylimmin < -3000:
            ylimmin = -300
        fiverow.set_ylim(imageMapOrderDF["flux_smoothed"].min() - 30, imageMapOrderDF["flux_smoothed"].max() * 1.2)
        fiverow.set_ylabel(
            "counts", fontsize=10)
        fiverow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        # BUILD IMAGE OF SKY MODEL
        skyModelImage = np.zeros_like(frame.data)
        for x, y, skypixel in zip(imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_model"]):
            skyModelImage[y][x] = skypixel
        nonOrderMask = (nonOrderMask == 0)
        skyModelImage = np.ma.array(skyModelImage, mask=nonOrderMask)
        cmap = copy(cm.gray)
        std = np.nanstd(skyModelImage)
        mean = np.nanmean(skyModelImage)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        im = sixrow.imshow(np.flipud(np.rot90(skyModelImage, 1)), vmin=0, vmax=100, cmap=cmap, alpha=1.)
        sixrow.set_ylabel(
            "x-axis", fontsize=10)
        sixrow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        sixrow.set_xlim(ylimMinImage, ylimMaxImage)
        # sixrow.invert_xaxis()
        medianValue = np.median(skyModelImage.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="sky model")]
        sixrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        sixrow.set_xticks([], [])

        # BUILD SKY-SUBTRACTED IMAGE
        skySubImage = np.zeros_like(frame.data)
        for x, y, skypixel in zip(imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_subtracted_flux"]):
            skySubImage[y][x] = skypixel
        skySubMask = (nonOrderMask == 1)
        skySubImage = np.ma.array(skySubImage, mask=skySubMask)
        cmap = copy(cm.gray)
        std = np.nanstd(skySubImage)
        mean = np.nanmedian(skySubImage)
        vmax = mean + 0.2 * std
        vmin = mean - 0.2 * std
        im = sevenrow.imshow(np.flipud(np.rot90(skySubImage, 1)), vmin=0, vmax=50, cmap=cmap, alpha=1.)
        sevenrow.set_title("STEP 3. Subtract the sky-model from the original data.", fontsize=10)
        sevenrow.set_xlabel(
            "y-axis", fontsize=10)
        sevenrow.set_ylabel(
            "x-axis", fontsize=10)
        sevenrow.set_ylim(imageMapOrderDF[self.axisA].min() - 10, imageMapOrderDF[self.axisA].max() + 10)
        sevenrow.set_xlim(ylimMinImage, ylimMaxImage)
        # sevenrow.invert_xaxis()
        medianValue = np.median(skySubImage.data.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="sky-subtracted frame")]
        sevenrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        # SUBTRACTED SKY RESIDUAL PANEL
        eightrow.scatter(imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "wavelength"].values, imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "sky_subtracted_flux"].values, s=3, alpha=0.2, c="orange", zorder=1, label="slit position > 0")
        eightrow.scatter(imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "wavelength"].values, imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "sky_subtracted_flux"].values, s=3, alpha=0.2, c=blue, zorder=1, label="slit position < 0")
        eightrow.scatter(knotLocations, np.zeros_like(knotLocations), marker=7, s=15, alpha=0.7, c=red, zorder=3, label="knots")
        eightrow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        mean = np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux"]).mean()
        std = np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux"]).std()

        eightrow.set_ylim(-1000, 1000)
        eightrow.set_xlabel(
            "wavelength (nm)", fontsize=10)
        eightrow.set_ylabel("residual", fontsize=10)

        # SUBTRACTED SKY RESIDUAL/ERROR PANEL
        ninerow.scatter(imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "wavelength"].values, imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "sky_subtracted_flux_weighted"].values, s=3, alpha=0.2, c="orange", zorder=1)
        ninerow.scatter(imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "wavelength"].values, imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "sky_subtracted_flux_weighted"].values, s=3, alpha=0.2, c=blue, zorder=1)
        ninerow.scatter(knotLocations, np.zeros_like(knotLocations), marker=7, s=15, alpha=0.7, c=red, zorder=3)

        mean = np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux_weighted"])[100:-100].mean()
        std = np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux_weighted"])[100:-100].std()

        ninerow.set_ylim(mean - 10 * std, mean + 10 * std)
        ninerow.set_xlabel(
            "wavelength (nm)", fontsize=10)
        ninerow.set_ylabel("residual (weighted)", fontsize=10)

        fig.suptitle(f"{self.arm} sky model: order {order}", fontsize=12, y=0.97)

        filename = self.filenameTemplate.replace(".fits", f"_SKYMODEL_QC_PLOTS_ORDER_{int(order)}.pdf")

        filePath = f"{self.qcDir}/{filename}"

        plt.show()
        plt.savefig(filePath, dpi='figure')
        plt.close()

        self.log.debug('completed the ``plot_sky_sampling`` method')
        return filePath

    def rolling_window_clipping(
            self,
            imageMapOrderDF,
            windowSize,
            sigma_clip_limit=5,
            max_iterations=10):
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
        self.log.debug('starting the ``rolling_window_clipping`` method')

        import numpy as np
        from astropy.stats import sigma_clip, mad_std

        allPixels = len(imageMapOrderDF.index)
        order = imageMapOrderDF["order"].values[0]

        # CALCULATE PERCENTILE SMOOTH DATA & RESIDUALS
        imageMapOrderDF["flux_smoothed"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].rolling(window=windowSize, center=True).quantile(.30)

        imageMapOrderDF["flux_minus_smoothed_residual"] = imageMapOrderDF["flux"] - imageMapOrderDF["flux_smoothed"]

        notClippedOrObjectMask = ((imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False))
        residuals = imageMapOrderDF.loc[notClippedOrObjectMask, "flux"] - imageMapOrderDF.loc[notClippedOrObjectMask, "flux_smoothed"]

        # SIGMA-CLIP THE DATA
        masked_residuals = sigma_clip(
            residuals, sigma_lower=3000, sigma_upper=sigma_clip_limit, maxiters=max_iterations, cenfunc='median', stdfunc=mad_std)
        imageMapOrderDF.loc[notClippedOrObjectMask, "object"] = masked_residuals.mask

        totalClipped = len(imageMapOrderDF.loc[(imageMapOrderDF["object"] == True)].index)
        percent = (float(totalClipped) / float(allPixels)) * 100.

        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        percent = (float(totalClipped) / float(allPixels)) * 100.
        self.log.print(f'\tORDER {order}: {totalClipped} pixels clipped in total = {percent:1.1f}%)')

        std = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"].std()
        imageMapOrderDF["residual_global_sigma"] = imageMapOrderDF["flux_minus_smoothed_residual"] / std

        self.log.debug('completed the ``rolling_window_clipping`` method')
        return imageMapOrderDF

    def fit_bspline_curve_to_sky(
            self,
            imageMapOrder,
            bspline_order):
        """*fit a single-order univariate bspline to the unclipped sky pixels (wavelength vs flux)*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe, containing sky flux with object(s) and CRHs removed
            - ``order`` -- the order number
            - ``bspline_order`` -- the order of the bspline to fit

        **Return:**
            - ``imageMapOrder`` -- same `imageMapOrder` as input but now with `sky_model` (bspline fit of the sky) and `sky_subtracted_flux` columns
            - ``tck`` -- the fitted bspline components. t for knots, c of coefficients, k for order

        **Usage:**

        ```python
        imageMapOrder, tck = self.fit_bspline_curve_to_sky(
            imageMapOrder,
            bspline_order
        )
        ```

        """
        self.log.debug('starting the ``fit_bspline_curve_to_sky`` method')

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

        residual_floor_percentile = self.recipeSettings["sky-subtraction"]["residual_floor_percentile"]

        # SORT BY COLUMN NAME - RUN ON A COPY SO NOT TO CONFUSE THE REST OF CODE
        imageMapOrder.sort_values(by=['wavelength'], inplace=True)
        order = imageMapOrder['order'].values[0]

        # CREATE ARRAYS NEEDED FOR BSPLINE FITTING
        goodWl = imageMapOrder.loc[imageMapOrder["clipped"] == False, "wavelength"]
        print("FIX WEIGHTS")
        imageMapOrder["weights"] = 1 / imageMapOrder["error"]

        # THESE WEIGHTS ARE NOT WORKING AS WELL AS 1/ERROR
        # imageMapOrder["weights"] = 1 / imageMapOrder["residual_windowed_sigma"].abs()
        # imageMapOrder["weights"] = imageMapOrder["weights"].replace(np.nan, 0.00000000)
        # imageMapOrder["weights"] = imageMapOrder["weights"].apply(lambda x: min(1, x))
        # imageMapOrder["weights"] = imageMapOrder["weights"].apply(lambda x: max(0.1, x))

        # WE WILL UPDATE THIS VALUE LATER IN WORKFLOW WITH SLIT-ILLUMINATION CORRECTION
        imageMapOrder["sky_model_slit"] = 0

        # FOR WEIGHTED BSPLINES WE ONLY NEED *INTERIOR* KNOTS (DON'T GO BEYOND RANGE OF DATA)
        # CAN'T HAVE MORE KNOTS THAN DATA POINTS
        # NUMBER OF 'DEFAULT' KNOTS
        defaultPointsPerKnot = self.recipeSettings["sky-subtraction"]["starting_points_per_knot"]

        if self.binx > 1:
            defaultPointsPerKnot /= self.binx
        if self.biny > 1:
            defaultPointsPerKnot /= self.biny
        defaultPointsPerKnot = int(defaultPointsPerKnot)

        n_interior_knots = int(goodWl.values.shape[0] / defaultPointsPerKnot)
        # QUANTILE SPACES - i.e. PERCENTAGE VALUES TO PLACE THE KNOTS, FROM 0-1, ALONGS WAVELENGTH RANGE
        qs = np.linspace(0, 1, n_interior_knots + 2)[1: -1]
        allKnots = np.quantile(goodWl, qs)
        extraKnots = np.array([])
        iterationCount = 0
        residualFloor = False
        residualFloorIterationLimit = 5
        slitIlluminationCorrectionIteration = 3

        while iterationCount < bsplineIterations:
            iterationCount += 1

            # CREATE ARRAYS NEEDED FOR BSPLINE FITTING
            goodWl = imageMapOrder.loc[imageMapOrder["clipped"] == False, "wavelength"]
            goodFlux = imageMapOrder.loc[imageMapOrder["clipped"] == False, "flux"] - imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_model_slit"]
            goodWeights = imageMapOrder.loc[imageMapOrder["clipped"] == False, "weights"]

            if iterationCount == slitIlluminationCorrectionIteration:
                # FIT SLIT-ILUMINATION PROFILE
                imageMapOrder = self.cross_dispersion_flux_normaliser(imageMapOrder)

            if iterationCount > 1:
                # POTENTIAL NEW KNOTS PLACED HALF WAY BETWEEN ADJACENT CURRENT KNOTS
                meanResiduals = []
                df = imageMapOrder.loc[imageMapOrder["clipped"] == False]

                ind = np.digitize(df["wavelength"], allKnots)
                group = imageMapOrder.loc[imageMapOrder["clipped"] == False].groupby(ind)

                # sys.exit(0)
                meanResiduals = group["sky_subtracted_flux_weighted_abs"].median()

                counts = group.size()
                potentialNewKnots = group["wavelength"].median()
                mask = counts < min_points_per_knot
                meanResiduals[mask] = residualFloor - 1

                meanResiduals = np.array(meanResiduals)
                potentialNewKnots = np.array(potentialNewKnots)
                # IF MEAN RESIDUAL BELOW PRE-AGREED FLOOR THEN SKIP NEW KNOT
                mask = np.ma.masked_where(meanResiduals < residualFloor, meanResiduals).mask
                # ELSE ADD NEW KNOT IF ABOVE FLOOT
                extraKnots = np.ma.compressed(np.ma.masked_array(potentialNewKnots, mask))
                meanResiduals = np.ma.compressed(np.ma.masked_array(meanResiduals, mask))

                allKnots = np.sort(np.concatenate((extraKnots, allKnots)))

                if order == 20:
                    print(residualFloor)
                    print(f"EXTRA KNOTS: {len(extraKnots)} .... {len(allKnots)} ... {iterationCount}")

            tck, fp, ier, msg = ip.splrep(goodWl, goodFlux, t=allKnots, k=bspline_order, w=goodWeights, full_output=True)
            t, c, k = tck

            if ier == 10:
                self.log.info(f"\t\tpoor fit on iteration {iterationCount} for order {imageMapOrder['order'].values[0]}. Reverting to last iteration.\n")
                tck = tck_previous
                break
            else:
                tck_previous = tck

            # GENERATE SKY-MODEL FROM BSPLINE
            imageMapOrder["sky_model_wl"] = ip.splev(imageMapOrder["wavelength"].values, tck)
            imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model_wl"] - imageMapOrder["sky_model_slit"]
            imageMapOrder["sky_subtracted_flux_weighted"] = imageMapOrder["sky_subtracted_flux"] / imageMapOrder["residual_windowed_std"]
            imageMapOrder["sky_subtracted_flux_weighted_abs"] = imageMapOrder["sky_subtracted_flux_weighted"].abs()

            if iterationCount <= residualFloorIterationLimit:
                # RECALUATE THE RESIDUAL FLOOR WE ARE CONVERGING TO
                allResiduals = imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_subtracted_flux_weighted_abs"]
                meanResidual = np.mean(allResiduals[1000:-1000])
                medianResidual = np.median(allResiduals[1000:-1000])
                p90Residual = np.nanpercentile(allResiduals[1000:-1000], residual_floor_percentile)
                residualFloor = p90Residual
                print(np.nanpercentile(allResiduals[1000:-1000], 90), np.nanpercentile(allResiduals[1000:-1000], 80), np.nanpercentile(allResiduals[1000:-1000], 70))
            else:
                imageMapOrder["sky_subtracted_flux_weighted"] = imageMapOrder["sky_subtracted_flux"] / (imageMapOrder["residual_windowed_std"] * 4.)
                imageMapOrder["sky_subtracted_flux_weighted_abs"] = imageMapOrder["sky_subtracted_flux_weighted"].abs()
                residuals = imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_subtracted_flux_weighted"]

                if iterationCount > bsplineIterations - 3:
                    # SIGMA-CLIP THE DATA
                    flux_error_ratio = imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_subtracted_flux_weighted"].values
                    masked_residuals = sigma_clip(
                        residuals, sigma_lower=bsplineSigma, sigma_upper=bsplineSigma, maxiters=7, cenfunc='mean', stdfunc='std')
                    imageMapOrder.loc[imageMapOrder["clipped"] == False, "bspline_clipped"] = masked_residuals.mask
                    imageMapOrder.loc[imageMapOrder["clipped"] == False, "clipped"] = masked_residuals.mask

            flux_error_ratio = imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_subtracted_flux_weighted"].values

            if flux_error_ratio[1000:-1000].shape[0]:
                flux_error_ratio = flux_error_ratio[1000:-1000]

            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print(f'\tOrder: {order}, Iteration {iterationCount}, RES {flux_error_ratio.mean():0.3f}, STD {flux_error_ratio.std():0.3f}, MEDIAN {np.median(flux_error_ratio):0.3f}, MAX {flux_error_ratio.max():0.3f}, MIN {flux_error_ratio.min():0.3f}')
            # self.log.print(fp, ier, msg)

        imageMapOrder["sky_model_wl"] = ip.splev(imageMapOrder["wavelength"].values, tck)
        imageMapOrder["sky_model"] = imageMapOrder["sky_model_wl"] + imageMapOrder["sky_model_slit"]
        # REPLACE VALUES LESS THAN ZERO IN COLUMN WITH ZERO
        imageMapOrder["sky_model"] = imageMapOrder["sky_model"].apply(lambda x: max(0, x))

        imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model"]
        imageMapOrder["sky_subtracted_flux_weighted"] = imageMapOrder["sky_subtracted_flux"] / imageMapOrder["residual_windowed_std"]

        print(defaultPointsPerKnot)
        imageMapOrder["sky_subtracted_flux_rolling_median"] = imageMapOrder["sky_subtracted_flux"].abs().rolling(defaultPointsPerKnot).median()
        flux_error_ratio = imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_subtracted_flux_weighted"].values
        if flux_error_ratio[1000:-1000].shape[0]:
            flux_error_ratio = flux_error_ratio[1000:-1000]

        self.log.debug('completed the ``fit_bspline_curve_to_sky`` method')
        return imageMapOrder, tck, allKnots, flux_error_ratio, residualFloor

    def create_placeholder_images(
            self):
        """*create placeholder images for the sky model and sky-subtracted frame*

        **Return:**
            - ``skymodelCCDData`` -- placeholder for sky model image
            - ``skySubtractedCCDData`` -- placeholder for sky-subtracted image

        **Usage:**

        ```python
        skymodelCCDData, skySubtractedCCDData = self.create_placeholder_images()
        ```
        """
        self.log.debug('starting the ``create_placeholder_images`` method')

        import numpy as np

        # CREATE AN IMAGE ARRAY TO HOST WAVELENGTH AND SLIT-POSITIONS
        skymodelCCDData = self.objectFrame.copy()
        skymodelCCDData.data[:] = np.nan
        skySubtractedCCDData = skymodelCCDData.copy()
        skySubtractedResidualsCCDData = skymodelCCDData.copy()

        self.log.debug('completed the ``create_placeholder_images`` method')
        return skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData

    def add_data_to_placeholder_images(
            self,
            imageMapOrderDF,
            skymodelCCDData,
            skySubtractedCCDData,
            skySubtractedResidualsCCDData):
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
        self.log.debug('starting the ``add_data_to_placeholder_images`` method')

        for x, y, skypixel in zip(imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_model"]):
            skymodelCCDData.data[y][x] = skypixel

        for x, y, skypixel in zip(imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_subtracted_flux"]):
            skySubtractedCCDData.data[y][x] = skypixel
        for x, y, skypixel in zip(imageMapOrderDF[self.axisA], imageMapOrderDF[self.axisB], imageMapOrderDF["sky_subtracted_flux"] / imageMapOrderDF["error"]):
            skySubtractedResidualsCCDData.data[y][x] = skypixel

        self.log.debug('completed the ``add_data_to_placeholder_images`` method')
        return skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData

    def plot_image_comparison(
            self,
            objectFrame,
            skyModelFrame,
            skySubFrame):
        """*generate a plot of original image, sky-model and sky-subtraction image*

        **Key Arguments:**
            - ``objectFrame`` -- object frame
            - ``skyModelFrame`` -- sky model frame
            - ``skySubFrame`` -- sky subtracted frame

        **Return:**
            - ``filePath`` -- path to the plot pdf
        """
        self.log.debug('starting the ``plot_results`` method')

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

        try:
            std = np.nanstd(maskedDataArray)
            mean = np.nanmean(maskedDataArray)
        except:
            std = np.std(maskedDataArray)
            mean = np.mean(maskedDataArray)
        vmax = mean + 1 * std
        vmin = mean - 0.1 * std
        toprow.imshow(rotatedImg, vmin=0, vmax=100, cmap='gray', alpha=1.)
        toprow.set_title(
            f"Original {arm} Frame", fontsize=10)
        toprow.set_ylabel("x-axis", fontsize=8)
        toprow.set_xlabel("y-axis", fontsize=8)
        toprow.tick_params(axis='both', which='major', labelsize=9)

        rotatedImg = np.rot90(skyModelFrame.data, 1)
        maskedDataArray = np.ma.array(skyModelFrame.data, mask=combinedMask)
        std = np.nanstd(maskedDataArray)
        mean = np.nanmean(maskedDataArray)
        vmax = mean + 1 * std
        vmin = mean - 1 * std
        midrow.imshow(rotatedImg, vmin=0, vmax=100, cmap='gray', alpha=1.)
        midrow.set_title(
            f"Sky-model for {arm} Frame", fontsize=10)
        midrow.set_ylabel("x-axis", fontsize=8)
        midrow.set_xlabel("y-axis", fontsize=8)
        midrow.tick_params(axis='both', which='major', labelsize=9)

        rotatedImg = np.rot90(skySubFrame.data, 1)
        maskedDataArray = np.ma.array(skySubFrame.data, mask=combinedMask)

        try:
            std = np.nanstd(maskedDataArray)
            mean = np.nanmean(maskedDataArray)
        except:
            std = np.std(maskedDataArray)
            mean = np.mean(maskedDataArray)

        vmax = 0 + std
        vmin = 0
        bottomrow.imshow(rotatedImg, vmin=vmin, vmax=30, cmap='gray', alpha=1.)
        bottomrow.set_title(
            f"Sky-subtracted {arm} Frame", fontsize=10)
        bottomrow.set_ylabel("x-axis", fontsize=8)
        bottomrow.set_xlabel("y-axis", fontsize=8)
        bottomrow.tick_params(axis='both', which='major', labelsize=9)

        # plt.show()
        filename = self.filenameTemplate.replace(".fits", "_skysub_quicklook.pdf")

        filePath = f"{self.qcDir}/{filename}"
        plt.savefig(filePath, dpi=720)
        plt.close()

        self.log.debug('completed the ``plot_results`` method')
        return filePath

    def rectify_order(
            self,
            order,
            imageMapOrder,
            remove_clipped=False,
            conserve_flux=False):
        """*rectify order on a fine slit-postion, wavelength grid*

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

        ---

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - write a command-line tool for this method
            - update package tutorial with command-line tool info if needed
        ```
        """
        self.log.debug('starting the ``rectify_order`` method')

        import numpy as np
        import pandas as pd

        dispMap = self.dispMap
        kw = self.kw
        dp = self.detectorParams
        arm = self.arm

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm)

        for o, minWl, maxWl in zip(orderNums, waveLengthMin, waveLengthMax):
            if o == order:
                orderInfo = (order, minWl, maxWl)
        (order, minWl, maxWl) = orderInfo

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

        halfGrid = (slitLength / 2)
        slitArray = np.arange(-halfGrid, halfGrid +
                              straighten_grid_res_slit, straighten_grid_res_slit)

        wlArray = np.arange(minWl, maxWl, straighten_grid_res_wavelength)

        # ONE SINGLE-VALUE SLIT ARRAY FOR EVERY WAVELENGTH ARRAY
        bigSlitArray = np.concatenate(
            [np.ones(wlArray.shape[0]) * slitArray[i] for i in range(0, slitArray.shape[0])])
        # NOW THE BIG WAVELEGTH ARRAY
        bigWlArray = np.tile(wlArray, np.shape(slitArray)[0])

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
            dispersionMapPath=self.dispMap,
            orderPixelTable=orderPixelTable,
            removeOffDetectorLocation=False
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

        orderPixelTable = pd.merge(orderPixelTable, imageMapOrder[['x', 'y', 'flux', 'clipped']], how='left', left_on=[
            'pixel_x', 'pixel_y'], right_on=['x', 'y'])

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable['flux'].isnull())
        self.log.print(orderPixelTable.loc[~mask, "wavelength"].min())

        # DROP MISSING VALUES
        # orderPixelTable.dropna(axis='index', how='any', subset=['x'], inplace=True)

        # orderPixelTable = orderPixelTable[['order', 'wavelength', 'slit_position', 'fit_x', 'fit_y', 'flux', 'clipped']]
        # orderPixelTable['weight'] = 100

        if conserve_flux:
            # ADD A COUNT COLUMN FOR THE NUMBER OF SMALL SLIT/WL PIXELS FALLING IN LARGE DETECTOR PIXELS
            count = orderPixelTable.groupby(['pixel_x', 'pixel_y']).size().reset_index(name='count')
            orderPixelTable = pd.merge(orderPixelTable, count, how='left', left_on=['pixel_x', 'pixel_y'], right_on=['pixel_x', 'pixel_y'])

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        if remove_clipped:
            mask = (orderPixelTable['clipped'] == True)
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
            log=self.log, CCDObject=imageArray, show=False, ext='data', stdWindow=3, title=False, surfacePlot=True, inst="dummy")

        self.log.debug('completed the ``rectify_order`` method')
        return imageArray

    def calculate_residuals(
            self,
            skyPixelsDF,
            fluxcoeff,
            orderDeg,
            wavelengthDeg,
            slitDeg,
            writeQCs=False):
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
        self.log.debug('starting the ``calculate_residuals`` method')

        import numpy as np

        arm = self.arm

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # POLY FUNCTION NEEDS A DATAFRAME AS INPUT
        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDeg, wavelengthDeg=wavelengthDeg, slitDeg=slitDeg, exponentsIncluded=True).poly

        # CALCULATE RESIDUALS BETWEEN MEASURED FLUX AND POLY
        # FITTED FLUX
        skyPixelsDF["fit_sky_subtracted_flux"] = poly(skyPixelsDF, *fluxcoeff)
        skyPixelsDF["residuals_sky_subtracted_flux"] = skyPixelsDF[
            "fit_sky_subtracted_flux"] - skyPixelsDF["sky_subtracted_flux"]

        # CALCULATE COMBINED RESIDUALS AND STATS
        res_mean = np.mean(skyPixelsDF["residuals_sky_subtracted_flux"])
        res_std = np.std(skyPixelsDF["residuals_sky_subtracted_flux"])
        res_median = np.median(skyPixelsDF["residuals_sky_subtracted_flux"])

        self.log.debug('completed the ``calculate_residuals`` method')
        return res_mean, res_std, res_median, skyPixelsDF

    def clip_object_slit_positions(
            self,
            order_dataframes,
            aggressive=False):
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
        self.log.debug('starting the ``clip_object_slit_positions`` method')

        import numpy as np
        import pandas as pd

        # COMBINE ALL ORDERS AND KEEP ONLY PIXELS FLAGGED AS POTENTIAL OBJECT
        allimageMapOrder = pd.concat(order_dataframes)
        mask = (allimageMapOrder['object'] == True)
        allimageMapOrder = allimageMapOrder.loc[mask]

        if aggressive:

            # BIN FLAGGED PIXEL COUNTS INTO DISCRETE SLIT-POSTION RANGES
            nbins = 100
            minsp = allimageMapOrder['slit_position'].min()
            maxsp = allimageMapOrder['slit_position'].max()
            bins = np.linspace(minsp, maxsp, nbins)
            result = allimageMapOrder['slit_position'].value_counts(bins=bins, sort=False, normalize=True) * nbins

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
                rects1 = ax.bar(bins - width / 2, result, width, label='count')
                fig.tight_layout()
                plt.show()

        # NOW FOR EACH OBJECT SLIT-RANGE, FLAG AS CLIPPED IN ORIGINAL ORDER DATAFRAMES
        for df in order_dataframes:
            if aggressive:
                for object in object_ranges:
                    df.loc[(df['slit_position'].between(object[0], object[1])), "clipped"] = True
                    df.loc[(df['slit_position'].between(object[0], object[1])), "object"] = True
            else:
                # df.loc[((df['slit_position'].between(object[0], object[1])) & (df['object'] == True)), "clipped"] = True
                df.loc[((df['object'] == True)), "clipped"] = True
            df.loc[((df['clipped'] == False) & (df['object'] == True)), "object"] = False

            notClippedOrObjectMask = ((df["clipped"] == False) & (df["object"] == False))
            df.loc[notClippedOrObjectMask, "residual_windowed_std"] = df.loc[notClippedOrObjectMask, "flux_minus_smoothed_residual"].rolling(self.rollingWindowSize, center=True).std()
            df.loc[notClippedOrObjectMask, "residual_windowed_sigma"] = df.loc[notClippedOrObjectMask, "flux_minus_smoothed_residual"] / df.loc[notClippedOrObjectMask, "residual_windowed_std"]

        self.log.debug('completed the ``clip_object_slit_positions`` method')
        return order_dataframes

    def cross_dispersion_flux_normaliser(
            self,
            orderDF):
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
        self.log.debug('starting the ``slit-profile-flux-normaliser`` method')

        import scipy.interpolate
        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt

        slit_illumination_order = self.recipeSettings["sky-subtraction"]["slit_illumination_order"]

        # 3 IMAGE DIMENSIONS - WAVELENGTH, SLIT-POSTIION AND FLUX
        mask = ((orderDF["clipped"] == False) & (orderDF["object"] == False))
        thisOrder = orderDF.loc[mask]
        wl = thisOrder['wavelength'].values
        sp = thisOrder['slit_position'].values
        fx = thisOrder['sky_subtracted_flux'].values

        # SOME LIGHT CLIPPING
        masked_residuals = sigma_clip(fx, sigma_lower=15, sigma_upper=15, maxiters=3, cenfunc='mean', stdfunc='std')
        a = [wl, sp, fx]
        wl, sp, fx = [np.ma.compressed(np.ma.masked_array(
            i, masked_residuals.mask)) for i in a]

        # ITERATIVELY FIT SLIT-ILLUMINATION
        iteration = 1
        while iteration < 10:
            iteration += 1
            coeff = np.polyfit(sp, fx, deg=slit_illumination_order)
            residuals = fx - np.polyval(coeff, sp)
            masked_residuals = sigma_clip(residuals, sigma_lower=3, sigma_upper=3, maxiters=1, cenfunc='mean', stdfunc='std')
            # REDUCE ARRAYS TO NON-MASKED VALUES
            a = [sp, fx]
            sp, fx = [np.ma.compressed(np.ma.masked_array(
                i, masked_residuals.mask)) for i in a]

        if False:
            fig = plt.figure(figsize=(15, 4))
            plt.title("Slit Illumination Profile")
            plt.plot(sp, fx, '.', ms=0.3, color="blue")
            plt.xlabel('slit-position', fontsize=8)
            plt.ylabel('sky_subtracted_flux', fontsize=8)
            xp = np.linspace(sp.min(), sp.max(), 100)
            fitFx = np.polyval(coeff, xp)
            plt.plot(xp, fitFx, color="red")
            plt.show()

        # USE THE SLIT-ILLUMINATION FUNCTION TO CREATE A FLUX-NORMALISATION
        thisOrder['sky_model_slit'] = np.polyval(coeff, thisOrder['slit_position'].values)
        orderDF['sky_model_slit'] = np.polyval(coeff, orderDF['slit_position'].values)

        if False:
            # CHECKING RESULTS - SHOULD NOW HAVE A STRAIGHT LINE PROFILE
            thisOrder['sky_subtracted_flux'] = thisOrder['flux'] - thisOrder['sky_model_wl'] - thisOrder['sky_model_slit']
            sp = thisOrder['slit_position'].values
            fx = thisOrder['sky_subtracted_flux'].values
            xp = np.linspace(sp.min(), sp.max(), 100)

            fig = plt.figure(figsize=(15, 4))
            plt.title("Original Sky Model Corrected for Slit Illumination")
            plt.plot(sp, fx, '.', ms=0.3, color="blue")
            plt.xlabel('slit-position', fontsize=8)
            plt.ylabel('sky_subtracted_flux (corrected)', fontsize=8)
            fitFx = np.polyval(coeff, xp)
            plt.plot(xp, fitFx, color="red")
            plt.show()

        self.log.debug('completed the ``slit-profile-flux-normaliser`` method')
        return orderDF

    # use the tab-trigger below for new method
    # xt-class-method
