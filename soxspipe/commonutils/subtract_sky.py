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
        self.mapDF, self.interOrderMask = twoD_disp_map_image_to_dataframe(log=self.log, slit_length=dp["slit_length"], twoDMapPath=twoDMap, associatedFrame=self.objectFrame, kw=kw, dispAxis=self.detectorParams["dispersion-axis"])

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
        - ``skySubtractedCCDData`` -- CCDData object containing the sky-subtracted frame
        - ``qcTable`` -- the data frame containing measured QC metrics
        - ``productsTable`` -- the data frame containing collected output products
        """
        self.log.debug('starting the ``get`` method')

        import numpy as np
        import pandas as pd
        pd.options.mode.chained_assignment = None

        self.log.print(f'\n# MODELLING SKY BACKGROUND AND REMOVING FROM SCIENCE FRAME')

        # THESE PLACEHOLDERS ARE INITIALLY BLANK AND AWAITING PIXEL VALUES TO BE ADDED
        skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData = self.create_placeholder_images()

        uniqueOrders = self.mapDF['order'].unique()
        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # THE BSPLINE ORDER TO FIT SKY WITH
        self.bspline_order = self.recipeSettings["sky-subtraction"]["bspline_order"]

        # SELECT A SINGLE ORDER TO GENERATE QC PLOTS FOR
        #print("FIX ME")
        self.qcPlotOrder = int(np.median(uniqueOrders)) - 1
        #self.qcPlotOrder = 16
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
            imageMapOrder, tck, knots, flux_error_ratio, residualFloor = self.fit_bspline_curve_to_sky(imageMapOrder)
            totalKnots += len(knots)
            newAllimageMapOrder.append(imageMapOrder)
            allFluxErrorRatios.append(flux_error_ratio)
            allResidualFloor.append(residualFloor)
            alltck.append(tck)
            allKnots.append(knots)
        allimageMapOrder = newAllimageMapOrder

        allFluxErrorRatios = np.concatenate(allFluxErrorRatios)
        self.log.print(f'\n\tFULL FRAME SKY-MODEL FLUX TO ERROR METRICS: MEAN {allFluxErrorRatios.mean():0.3f}, STD {allFluxErrorRatios.std():0.3f}, MEDIAN {np.median(allFluxErrorRatios):0.3f}, MEAN RES FLOOR: {np.mean(allResidualFloor):0.3f}, KNOT COUNT: {totalKnots}')
        # self.log.print(f'\t{allFluxErrorRatios.mean():0.3f},  {allFluxErrorRatios.std():0.3f},  {np.median(allFluxErrorRatios):0.3f},  {allFluxErrorRatios.max():0.3f},  {allFluxErrorRatios.min():0.3f}, {allFluxErrorRatios.max()-allFluxErrorRatios.min():0.3f},{np.mean(allResidualFloor):0.3f},{totalKnots}')

        for o, imageMapOrder, tck, knots in zip(uniqueOrders, allimageMapOrder, alltck, allKnots):
            if isinstance(imageMapOrder, pd.core.frame.DataFrame):
                # INJECT THE PIXEL VALUES BACK INTO THE PLACEHOLDER IMAGES
                skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData = self.add_data_to_placeholder_images(imageMapOrder, skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData)
                if o == self.qcPlotOrder and True:
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
        - ``clipBPs`` -- clip bad-pixels? Default *True*
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
        purple = "purple"

        # MAKE A COPY OF THE FRAME TO NOT ALTER ORIGINAL DATA
        frame = self.objectFrame.copy()

        # SETUP THE PLOT SUB-PANELS
        # print("FIX ME")
        if True:
            fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=320)
        else:
            # REMOVE ME
            fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=150)

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
        im = onerow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=1)
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
        colours = [blue, orange, green, purple]
        alphas = [.2, .3, .8, 1.]
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
        colours = [blue, orange, purple]
        alphas = [.2, 1, 1]
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
        im = fourrow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=1)

        columnName = ["object", "bad_pixel_clipped", "edge_clipped", "bspline_clipped"]
        colours = [blue, orange, green, purple]
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
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "wavelength"].values,
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] > 0), "flux"].values, s=3, c=orange, alpha=1, zorder=1, label="unclipped slit position > 0")
        fiverow.scatter(
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "wavelength"].values,
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["slit_position"] < 0), "flux"].values, s=3, c=blue, alpha=1, zorder=1, label="unclipped slit position < 0")

        # fiverow.scatter(
        #     imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False), "wavelength"].values,
        #     imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False), "flux"].values, s=1, c=black, alpha=0.5, zorder=1, label="unclipped pixels")

        if False:
            fiverow.scatter(
                imageMapOrderDF.loc[imageMapOrderDF["bspline_clipped"] == True, "wavelength"].values,
                imageMapOrderDF.loc[imageMapOrderDF["bspline_clipped"] == True, "flux"].values, s=8, c=violet, marker="x", alpha=1., zorder=1, label="clipped")

        if tck:
            wl = np.linspace(imageMapOrderDF["wavelength"].min(), imageMapOrderDF["wavelength"].max(), 1000000)
            sky = ip.splev(wl, tck)
            # sky = ip.splev(wl, tck, der=1)
            # sky = ip.splev(wl, tck, der=2)
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
        im = sixrow.imshow(np.flipud(np.rot90(skyModelImage, 1)), vmin=vmin, vmax=vmax, cmap=cmap, alpha=1.)
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

        # plt.show()
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

        iteration = 0

        while iteration < max_iterations:
            iteration += 1
            # CALCULATE PERCENTILE SMOOTH DATA & RESIDUALS
            notClippedOrObjectMask = ((imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False))
            imageMapOrderDF.loc[notClippedOrObjectMask, "flux_smoothed"] = imageMapOrderDF.loc[notClippedOrObjectMask, "flux"].rolling(window=windowSize, center=True).quantile(.30)
            imageMapOrderDF.loc[notClippedOrObjectMask, "flux_std"] = imageMapOrderDF.loc[notClippedOrObjectMask, "flux"].rolling(window=windowSize, center=True).std()
            imageMapOrderDF["flux_minus_smoothed_residual"] = imageMapOrderDF["flux"] - imageMapOrderDF["flux_smoothed"]
            notClippedOrObjectMask = ((imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False))
            masked = (imageMapOrderDF.loc[notClippedOrObjectMask, "flux"] > imageMapOrderDF.loc[notClippedOrObjectMask, "flux_smoothed"] + sigma_clip_limit * imageMapOrderDF.loc[notClippedOrObjectMask, "flux_std"])
            imageMapOrderDF.loc[(notClippedOrObjectMask & masked), "object"] = True
            notClippedOrObjectMask = ((imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False))
            masked = (imageMapOrderDF.loc[notClippedOrObjectMask, "flux"] < imageMapOrderDF.loc[notClippedOrObjectMask, "flux_smoothed"] - 1. * imageMapOrderDF.loc[notClippedOrObjectMask, "flux_std"])
            imageMapOrderDF.loc[(notClippedOrObjectMask & masked), "object"] = True

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
            imageMapOrder):
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

        # SORT BY COLUMN NAME
        imageMapOrder.sort_values(by=['wavelength'], inplace=True)
        order = imageMapOrder['order'].values[0]

        # CREATE ARRAYS NEEDED FOR BSPLINE FITTING
        goodWl = imageMapOrder.loc[imageMapOrder["clipped"] == False, "wavelength"]

        # REPLACE NANS IN residual_windowed_std
        data = imageMapOrder["residual_windowed_std"].values
        mask = np.isnan(data)
        data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
        imageMapOrder["residual_windowed_std"] = data

        # imageMapOrder["weights"] = 1 / imageMapOrder["residual_windowed_std"].values
        imageMapOrder["weights"] = 1 / imageMapOrder["error"].values

        # WE WILL UPDATE THIS VALUE LATER IN WORKFLOW WITH SLIT-ILLUMINATION CORRECTION
        imageMapOrder["slit_normalisation_ratio"] = 1

        # FOR WEIGHTED BSPLINES WE ONLY NEED *INTERIOR* KNOTS (DON'T GO BEYOND RANGE OF DATA)
        # CAN'T HAVE MORE KNOTS THAN DATA POINTS
        # NUMBER OF 'DEFAULT' KNOTS
        defaultPointsPerKnot = self.recipeSettings["sky-subtraction"]["starting_points_per_knot"]

        # STARTER KNOTS USED TO MEASURE THE RESIDUAL FLOOR BEFORE REVERTING TO defaultPointsPerKnot
        if self.arm.upper == "NIR":
            starterPointsPerKnot = 25
        else:
            starterPointsPerKnot = 100

        if self.binx > 1:
            defaultPointsPerKnot /= self.binx
            starterPointsPerKnot /= self.binx
        if self.biny > 1:
            defaultPointsPerKnot /= self.biny
            starterPointsPerKnot /= self.biny

        defaultPointsPerKnot = int(defaultPointsPerKnot)
        starterPointsPerKnot = int(starterPointsPerKnot)

        n_interior_knots = int(goodWl.values.shape[0] / defaultPointsPerKnot)
        # QUANTILE SPACES - i.e. PERCENTAGE VALUES TO PLACE THE KNOTS, FROM 0-1, ALONGS WAVELENGTH RANGE
        qs = np.linspace(0, 1, n_interior_knots + 2)[1: -1]
        defaultKnots = np.quantile(goodWl, qs)

        n_interior_knots = int(goodWl.values.shape[0] / starterPointsPerKnot)
        # QUANTILE SPACES - i.e. PERCENTAGE VALUES TO PLACE THE KNOTS, FROM 0-1, ALONGS WAVELENGTH RANGE
        qs = np.linspace(0, 1, n_interior_knots + 2)[1: -1]
        starterKnots = np.quantile(goodWl, qs)

        extraKnots = np.array([])
        iterationCount = -3
        residualFloor = False

        slitIlluminationCorrectionIteration = 4
        tiltAdjustmentIteration = 4

        slitCorrectIterationLimit = 2
        slitCorrectIterations = 0

        while iterationCount < bsplineIterations:
            iterationCount += 1

            # CREATE ARRAYS NEEDED FOR BSPLINE FITTING
            goodWl = imageMapOrder.loc[imageMapOrder["clipped"] == False, "wavelength"]
            goodFlux = imageMapOrder.loc[imageMapOrder["clipped"] == False, "flux"] / imageMapOrder.loc[imageMapOrder["clipped"] == False, "slit_normalisation_ratio"]
            goodWeights = imageMapOrder.loc[imageMapOrder["clipped"] == False, "weights"]

            if iterationCount < 5:
                baseKnots = starterKnots
            if iterationCount == 5:
                baseKnots = defaultKnots

            allKnots = np.sort(np.concatenate((extraKnots, baseKnots)))

            if slitCorrectIterations < slitCorrectIterationLimit:

                # SLIT CORRECT ACTUALLY HELPS
                if iterationCount == slitIlluminationCorrectionIteration and True:
                    # FIT SLIT-ILUMINATION PROFILE
                    imageMapOrder = self.cross_dispersion_flux_normaliser(imageMapOrder)
                    extraKnots = np.array([])
                    iterationCount = -2
                    slitIlluminationCorrectionIteration = -99
                    tiltAdjustmentIteration = 4

                # FLUX SHUFFLING MAKE EXTRACTION WORSE
                if iterationCount == tiltAdjustmentIteration and self.arm.upper() in ["NIR"] and False:
                    # FIT SLIT-ILUMINATION PROFILE
                    imageMapOrder = self.adjust_tilt(imageMapOrder, tck)
                    extraKnots = np.array([])
                    iterationCount = -2
                    # tiltAdjustmentIteration = -99
                    slitCorrectIterations += 1

            if iterationCount > 1:
                # POTENTIAL NEW KNOTS PLACED HALF WAY BETWEEN ADJACENT CURRENT KNOTS
                meanResiduals = []

                unclippedAndSkyline = (imageMapOrder["clipped"] == False) & (imageMapOrder["sky_line"].isin(["rise", "fall"]))
                df = imageMapOrder.loc[unclippedAndSkyline]
                ind = np.digitize(df["wavelength"], allKnots)

                imageMapOrder["sky_subtracted_flux_abs"] = imageMapOrder["sky_subtracted_flux"].abs()
                group = imageMapOrder.loc[unclippedAndSkyline].groupby(ind)
                meanResiduals = group["sky_subtracted_flux_abs"].mean()

                counts = group.size()
                potentialNewKnots = group["wavelength"].mean()
                mask = counts < min_points_per_knot
                meanResiduals[mask] = residualFloor - 1

                meanResiduals = np.array(meanResiduals)
                potentialNewKnots = np.array(potentialNewKnots)
                mask = np.ma.masked_where(meanResiduals < residualFloor, meanResiduals).mask
                # ELSE ADD NEW KNOT IF ABOVE FLOOT
                newKnots = np.ma.compressed(np.ma.masked_array(potentialNewKnots, mask))
                meanResiduals = np.ma.compressed(np.ma.masked_array(meanResiduals, mask))

                allKnots = np.sort(np.concatenate((newKnots, allKnots)))
                extraKnots = np.sort(np.concatenate((newKnots, extraKnots)))

                # NOW ADD KNOTS AT SKYLINE NODES
                skylineNodes = (imageMapOrder["clipped"] == False) & (imageMapOrder["sky_line"] == "node")
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

            tck, fp, ier, msg = ip.splrep(goodWl, goodFlux, t=allKnots, k=self.bspline_order, w=goodWeights, full_output=True)
            t, c, k = tck

            if ier == 10:
                self.log.info(f"\t\tpoor fit on iteration {iterationCount} for order {imageMapOrder['order'].values[0]}. Reverting to last iteration.\n")
                tck = tck_previous
                break
            else:
                tck_previous = tck

            if iterationCount == -1:
                # FIRST PASS SIGMA CLIPPING OF BSPLINE
                residuals = imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_subtracted_flux"]
                masked_residuals = sigma_clip(
                    residuals, sigma_lower=bsplineSigma, sigma_upper=bsplineSigma, maxiters=3, cenfunc='mean', stdfunc='std')
                imageMapOrder.loc[imageMapOrder["clipped"] == False, "bspline_clipped"] = masked_residuals.mask
                imageMapOrder.loc[imageMapOrder["clipped"] == False, "clipped"] = masked_residuals.mask

            if iterationCount == 0 and residualFloor == False:
                imageMapOrder, residualFloor = self.determine_residual_floor(imageMapOrder, tck)

            # GENERATE SKY-MODEL FROM BSPLINE
            imageMapOrder["sky_model_wl"] = ip.splev(imageMapOrder["wavelength"].values, tck)
            imageMapOrder["sky_model_wl_derivative"] = ip.splev(imageMapOrder["wavelength"].values, tck, der=1)
            imageMapOrder["sky_subtracted_flux"] = (imageMapOrder["flux"] / imageMapOrder["slit_normalisation_ratio"]) - imageMapOrder["sky_model_wl"]
            imageMapOrder["sky_subtracted_flux_weighted"] = imageMapOrder["sky_subtracted_flux"] * imageMapOrder["sky_model_wl_derivative"].abs() / imageMapOrder["residual_windowed_std"]
            imageMapOrder["sky_subtracted_flux_weighted_abs"] = imageMapOrder["sky_subtracted_flux_weighted"].abs()

            flux_error_ratio = imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_subtracted_flux_weighted"].values

            if flux_error_ratio[1000:-1000].shape[0]:
                flux_error_ratio = flux_error_ratio[1000:-1000]

            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print(f'\tOrder: {order}, Iteration {iterationCount}, RES {flux_error_ratio.mean():0.3f}, STD {flux_error_ratio.std():0.3f}, MEDIAN {np.median(flux_error_ratio):0.3f}, MAX {flux_error_ratio.max():0.3f}, MIN {flux_error_ratio.min():0.3f}')
            # self.log.print(fp, ier, msg)

        imageMapOrder["sky_model_wl"] = ip.splev(imageMapOrder["wavelength"].values, tck)
        imageMapOrder["sky_model_wl_derivative"] = ip.splev(imageMapOrder["wavelength"].values, tck, der=1)
        imageMapOrder["sky_model"] = imageMapOrder["sky_model_wl"] * imageMapOrder["slit_normalisation_ratio"]
        # REPLACE VALUES LESS THAN ZERO IN COLUMN WITH ZERO
        imageMapOrder["sky_model"] = imageMapOrder["sky_model"].apply(lambda x: max(0, x))

        imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model"]
        imageMapOrder["sky_subtracted_flux_weighted"] = imageMapOrder["sky_subtracted_flux"] * imageMapOrder["sky_model_wl_derivative"].abs() / (imageMapOrder["residual_windowed_std"] * 10)
        imageMapOrder["sky_subtracted_flux_weighted_abs"] = imageMapOrder["sky_subtracted_flux_weighted"].abs()

        # print(defaultPointsPerKnot)
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

    def _rectify_order(
            self,
            order,
            imageMapOrder,
            remove_clipped=False,
            conserve_flux=False):
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

        percentile_rolling_window_size = self.recipeSettings["sky-subtraction"]["percentile_rolling_window_size"]

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
                    # df.loc[(df['slit_position'].between(object[0], object[1])), "object"] = True
                    df.loc[((df['object'] == True)), "clipped"] = True
            else:
                # df.loc[((df['slit_position'].between(object[0], object[1])) & (df['object'] == True)), "clipped"] = True
                df.loc[((df['object'] == True)), "clipped"] = True
            df.loc[((df['clipped'] == False) & (df['object'] == True)), "object"] = False

            notClippedOrObjectMask = ((df["clipped"] == False) & (df["object"] == False))

            df.loc[notClippedOrObjectMask, "residual_windowed_std"] = df.loc[notClippedOrObjectMask, "flux_minus_smoothed_residual"].rolling(percentile_rolling_window_size, center=True).std()

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
        self.log.debug('starting the ``cross_dispersion_flux_normaliser`` method')

        import scipy.interpolate
        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt
        import pandas as pd
        import scipy.interpolate as ip

        slit_illumination_order = self.recipeSettings["sky-subtraction"]["slit_illumination_order"]
        order = orderDF['order'].values[0]
        orderDF['slit_normalisation_ratio'] = 1

        # COLLECT THE SKYLINE PIXELS
        mask = (orderDF["clipped"] == False)
        # mask = ((orderDF["clipped"] == False) & (orderDF["sky_line"] == False))
        thisOrder = orderDF.loc[mask]

        # SORT BY COLUMN NAME
        thisOrder.sort_values(['slit_position'], inplace=True)

        # GROUP INTO DISCRETE WAVELENGTH BINS
        # DEFINE THE BINS FOR COLUMN 'wavelength'
        lower = float("%0.*f" % (1, thisOrder['slit_position'].min()))
        upper = float("%0.*f" % (1, thisOrder['slit_position'].max()))

        # DO AN INITIAL CLIPPING OF THE DATA
        # masked_residuals = sigma_clip(
        #     thisOrder["sky_subtracted_flux"], sigma_lower=3, sigma_upper=3, maxiters=5, cenfunc='mean', stdfunc='std')
        # thisOrder["clipped"] = masked_residuals.mask

        # FIT THE DATA WITH A POLYMONIAL
        mask = ((thisOrder["clipped"] == False))
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
        # QUANTILE SPACES - i.e. PERCENTAGE VALUES TO PLACE THE KNOTS, FROM 0-1, ALONGS WAVELENGTH RANGE
        qs = np.linspace(0, 1, n_interior_knots + 2)[1: -1]
        starterKnots = np.quantile(sp, qs)

        tck, fp, ier, msg = ip.splrep(sp, fx, k=1, w=wt, t=starterKnots, full_output=True)

        iteration = 0
        spCopy = np.copy(sp)
        fxCopy = np.copy(fx)
        while iteration < 3:
            iteration += 1
            coeff = np.polyfit(spCopy, fxCopy, deg=slit_illumination_order)
            residuals = fxCopy - np.polyval(coeff, spCopy)
            masked_residuals = sigma_clip(residuals, sigma_lower=3, sigma_upper=3, maxiters=1, cenfunc='mean', stdfunc='std')
            # REDUCE ARRAYS TO NON-MASKED VALUES
            a = [spCopy, fxCopy]
            spCopy, fxCopy = [np.ma.compressed(np.ma.masked_array(
                i, masked_residuals.mask)) for i in a]

        # FIX ME
        if order == self.qcPlotOrder and False:
            fig = plt.figure(figsize=(15, 4))
            plt.title("Slit Illumination Profile")

            plt.plot(sp, fx, '.', ms=1, color="red", alpha=0.3)

            xp = np.linspace(lower, upper, 100)
            fitFx = np.polyval(coeff, xp)
            plt.plot(xp, fitFx, color="blue")

            slitPos = np.linspace(lower, upper, 10000)
            sky = ip.splev(slitPos, tck)
            plt.plot(
                slitPos, sky, label='sky model', c="#268bd2", zorder=3)

            plt.xlabel('slit-position', fontsize=12)
            plt.ylabel('normalised flux', fontsize=12)
            plt.show()

        # USE THE SLIT-ILLUMINATION FUNCTION TO CREATE A FLUX-NORMALISATION

        orderDF['flux'] -= np.polyval(coeff, orderDF['slit_position'].values)
        # orderDF['flux'] -= ip.splev(orderDF['slit_position'].values, tck)
        orderDF['slit_normalisation_ratio'] = 1

        self.log.debug('completed the ``cross_dispersion_flux_normaliser`` method')
        return orderDF

    def adjust_tilt(
            self,
            orderDF,
            tck):
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
        self.log.debug('starting the ``adjust_tilt`` method')

        import scipy.interpolate
        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt
        import pandas as pd
        import scipy.interpolate as ip

        slit_illumination_order = self.recipeSettings["sky-subtraction"]["slit_illumination_order"]

        # COLLECT THE SKYLINE PIXELS
        # mask = ((orderDF["clipped"] == False) & (orderDF["sky_line"] != False))
        mask = ((orderDF["clipped"] == False))
        thisOrder = orderDF.loc[mask]
        order = orderDF['order'].values[0]

        # GROUP INTO DISCRETE WAVELENGTH BINS
        # DEFINE THE BINS FOR COLUMN 'wavelength'
        lower = float("%0.*f" % (1, thisOrder['wavelength'].min()))
        upper = float("%0.*f" % (1, thisOrder['wavelength'].max()))
        wlBins = np.linspace(lower, upper, 1000)
        # BIN RESULTS
        thisOrder['wl_bins'] = pd.cut(thisOrder['wavelength'], bins=wlBins)
        wlGroups = thisOrder[["pixelScale", "wl_bins"]].groupby(['wl_bins'])
        thisOrder['pixelScale'] = wlGroups['pixelScale'].transform(lambda x: x.mean())

        minimum = 1000000
        bestShift = 0.

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

        orderDF['wavelength'] = orderDF["wavelength"] - orderDF["slit_position"] * bestShift

        # orderDF['wavelength'] = orderDF["wavelength"] - orderDF["slit_position"] * orderDF["pixelScale"] * bestShift

        orderDF.sort_values(by=['wavelength'], inplace=True)

        self.log.debug('completed the ``adjust_tilt`` method')
        return orderDF

    def determine_residual_floor(
            self,
            imageMapOrder,
            tck):
        """*determine residual floor and flag sky-lines*

        **Key Arguments:**
            - ``imageMapOrderDF`` --  dataframe with various processed data for a given order
            - ``tck`` -- the fitted bspline components. t for knots, c of coefficients, k for order

        **Return:**

        - `imageMapOrder` -- same dataframe but now with sky-line locations flagged
        - `residualFloor` -- the residual floor determined within regions containing no skylines.

        **Usage:**

        ```python
        imageMapOrder, residualFloor = self.determine_residual_floor(imageMapOrder, tck)
        ```
        """
        self.log.debug('starting the ``determine_residual_floor`` method')

        import matplotlib.pyplot as plt
        import numpy as np
        import scipy.interpolate as ip
        from astropy.stats import sigma_clip
        from astropy.stats import sigma_clipped_stats

        order = imageMapOrder["order"].values[0]

        # USE THIS PERCENTILE TO DETERMINE THE RESIDUAL FLOOR
        residual_floor_percentile = self.recipeSettings["sky-subtraction"]["residual_floor_percentile"]

        # DO SOME CLIPPING ON THE INITIAL SKY SUBTRACTION RESIDUALS

        imageMapOrder["sky_residuals"] = np.nan
        imageMapOrder["sky_residuals_clipped"] = False
        imageMapOrder["sky_line"] = False
        imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_residuals"] = (imageMapOrder.loc[imageMapOrder["clipped"] == False, "flux"].values - ip.splev(imageMapOrder.loc[imageMapOrder["clipped"] == False, "wavelength"].values, tck)) / ip.splev(imageMapOrder.loc[imageMapOrder["clipped"] == False, "wavelength"].values, tck, der=1)

        # FUDGE FOR NON-DARK SUBTRACTED DATA
        imageMapOrder.replace([np.inf, -np.inf], 1, inplace=True)

        # SIGMA-CLIP THE DATA
        masked_residuals = sigma_clip(
            imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_residuals"], sigma_lower=25, sigma_upper=25, maxiters=1, cenfunc='median', stdfunc="mad_std")
        imageMapOrder.loc[imageMapOrder["clipped"] == False, "sky_residuals_clipped"] = masked_residuals.mask
        imageMapOrder.loc[imageMapOrder["clipped"] == False, "clipped"] = masked_residuals.mask

        # MARK REGIONS AS NOISY
        unclipped = (imageMapOrder["clipped"] == False)
        imageMapOrder["noisy"] = False

        mean, median, std = sigma_clipped_stats(imageMapOrder.loc[(unclipped), "residual_windowed_std"], sigma=50.0, stdfunc="std", cenfunc="mean", maxiters=2)

        # I STOPPED MASKING NOISY AREAS OF THE DATA AS THIS WAS CLIPPING OUT SKYLINES IN NON-NOISY AREAS
        # if std > 1.5:
        #     noiseLimit = (imageMapOrder["residual_windowed_std"] > mean + 5 * std)
        #     imageMapOrder.loc[(unclipped & noiseLimit), "noisy"] = True

        # MARK RISING SKYLINES
        unclipped = (imageMapOrder["clipped"] == False)
        unnoisy = (imageMapOrder["noisy"] == False)

        imageMapOrder.loc[(unclipped), "sky_d1"] = ip.splev(imageMapOrder.loc[unclipped, "wavelength"].values, tck, der=1)
        imageMapOrder.loc[(unclipped), "sky_d2"] = ip.splev(imageMapOrder.loc[unclipped, "wavelength"].values, tck, der=2)

        std = imageMapOrder.loc[(unclipped & unnoisy), "sky_d1"].std()
        mean = imageMapOrder.loc[(unclipped & unnoisy), "sky_d1"].mean()

        if self.arm.upper == "NIR":
            lineMultiplier = 0.2
            nodeMultiplier = 0.3
        else:
            lineMultiplier = 1.5
            nodeMultiplier = 2

        imageMapOrder.loc[(unclipped & unnoisy) & (imageMapOrder["sky_d1"] > mean + lineMultiplier * std), "sky_line"] = "rise"
        imageMapOrder.loc[(unclipped & unnoisy) & (imageMapOrder["sky_d1"] < mean - lineMultiplier * std), "sky_line"] = "fall"

        # MARK SKY INFECTION POINTS
        std = imageMapOrder.loc[(unclipped & unnoisy), "sky_d2"].std()
        mean = imageMapOrder.loc[(unclipped & unnoisy), "sky_d2"].mean()
        skyline = (imageMapOrder["sky_line"].isin(["rise", "fall"]))
        imageMapOrder.loc[(unclipped & unnoisy & skyline) & (imageMapOrder["sky_d2"] > mean + nodeMultiplier * std), "sky_line"] = "node"
        imageMapOrder.loc[(unclipped & unnoisy & skyline) & (imageMapOrder["sky_d2"] < mean - nodeMultiplier * std), "sky_line"] = "node"

        # CALCULATE THE RESIDUAL FLOOR
        nonSkyline = (imageMapOrder["sky_line"] == False)
        nonSkyResiduals = imageMapOrder.loc[(unclipped & nonSkyline & unnoisy), "flux"].values - ip.splev(imageMapOrder.loc[(unclipped & nonSkyline & unnoisy), "wavelength"].values, tck)
        residualFloor = np.percentile(np.abs(nonSkyResiduals), residual_floor_percentile)

        # DETERMINE SKYLINES PEAKS
        skyFlux = imageMapOrder.loc[(unclipped & unnoisy), "flux"].values
        skyFluxMean, skyFluxStd = np.median(skyFlux), skyFlux.std()
        highFluxMask = (imageMapOrder["flux"] > skyFluxMean + 1.5 * skyFluxStd)
        imageMapOrder.loc[(unclipped & unnoisy & nonSkyline & highFluxMask), "sky_line"] = "peak"

        if self.qcPlotOrder == order and False:

            # generate the figure for the plot
            fig = plt.figure(
                num=None,
                figsize=(24, 8),
                dpi=150,
                facecolor=None,
                edgecolor=None,
                frameon=True)

            # add the axes
            ax = fig.add_axes(
                [0.05, 0.1, 0.99, 0.8],
                polar=False,
                frameon=True)
            ax.set_xlabel('wavelength (nm)')
            ax.set_ylabel('Counts')
            # ax.set_xlim(0, 10)
            # ax.set_ylim(0, 12)

            unclipped = (imageMapOrder["clipped"] == False)
            noisy = (imageMapOrder["noisy"] == True)
            nonSkyline = (imageMapOrder["sky_line"] == False)
            skyline = (imageMapOrder["sky_line"].isin(["rise", "fall"]))
            skynode = (imageMapOrder["sky_line"].isin(["node"]))
            bsplineClipped = (imageMapOrder["bspline_clipped"] == True)

            # # ORIGINAL DATA POINTS
            # ax.scatter(
            #     imageMapOrder.loc[(unclipped), "wavelength"].values,
            #     imageMapOrder.loc[(unclipped), "flux"].values * 100 / (imageMapOrder.loc[(unclipped), "residual_windowed_std"].values), s=1, c="green", alpha=0.3, zorder=10, label="no sky")

            ax.scatter(
                imageMapOrder.loc[(unclipped & nonSkyline), "wavelength"].values,
                imageMapOrder.loc[(unclipped & nonSkyline), "flux"].values, s=1, c="green", alpha=0.3, zorder=10, label="no sky")
            # SKYLINES
            ax.scatter(
                imageMapOrder.loc[(unclipped & skyline), "wavelength"].values,
                imageMapOrder.loc[(unclipped & skyline), "flux"].values, s=1, c="orange", alpha=0.3, zorder=10, label="skyline")
            # # SKYLINES NODES
            ax.scatter(
                imageMapOrder.loc[(unclipped & skynode), "wavelength"].values,
                imageMapOrder.loc[(unclipped & skynode), "flux"].values, s=1, c="red", alpha=0.3, zorder=10, label="skyline node")

            if False:
                # SHOW BSPLINE CLIPPED FROM PREVIOUS ITERATION
                ax.scatter(
                    imageMapOrder.loc[(bsplineClipped), "wavelength"].values,
                    imageMapOrder.loc[(bsplineClipped), "flux"].values, s=1, c="blue", alpha=1, zorder=10, label="bspline clipped")

            if True:
                # SHOW NOISY DATA
                ax.scatter(
                    imageMapOrder.loc[(noisy), "wavelength"].values,
                    imageMapOrder.loc[(noisy), "flux"].values, s=1, c="blue", alpha=0.1, zorder=10, label="noisy data")

            if False:
                # PLOT SKY RESIDUALS AND CLIPPED RESIDUALS
                ax.scatter(
                    imageMapOrder.loc[unclipped, "wavelength"].values,
                    imageMapOrder.loc[unclipped, "sky_residuals"].values, s=1, c="green", alpha=1, zorder=13, label="unclipped")
                ax.scatter(
                    imageMapOrder.loc[(imageMapOrder["sky_residuals_clipped"] == True), "wavelength"].values,
                    imageMapOrder.loc[(imageMapOrder["sky_residuals_clipped"] == True), "sky_residuals"].values, s=5, c="red", marker="x", alpha=1, zorder=13, label="unclipped")

            if False:
                # SHOW ROLLING SCATTER

                ax.scatter(
                    imageMapOrder.loc[(unclipped & noiseLimit), "wavelength"].values,
                    imageMapOrder.loc[(unclipped & noiseLimit), "residual_windowed_std"].values, s=1, c="blue", alpha=1, zorder=10, label="bspline clipped")

            wl = np.linspace(imageMapOrder["wavelength"].min(), imageMapOrder["wavelength"].max(), 1000000)
            sky = ip.splev(wl, tck)
            sky_d1 = ip.splev(wl, tck, der=1)
            sky_d2 = ip.splev(wl, tck, der=2)
            ax.plot(
                wl, sky, label='sky model', c="#268bd2", zorder=3)

            if True:
                ax.plot(
                    wl, sky_d1, label='sky model derivative', c="pink", zorder=3)
                ax.plot(
                    wl, sky_d2 / 100, label='sky model 2nd derivative / 100', c="violet", zorder=3)

            # Put a legend on plot
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='upper right', bbox_to_anchor=(1.1, 0.5), prop={'size': 8})

            title = "label"
            plt.show()
            plt.clf()  # clear figure

            # sys.exit(0)

        self.log.debug('completed the ``determine_residual_floor`` method')
        return imageMapOrder, residualFloor

    # use the tab-trigger below for new method
    # xt-class-method
