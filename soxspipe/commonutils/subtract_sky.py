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
from copy import copy
from datetime import datetime

from soxspipe.commonutils import keyword_lookup

from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.toolkit import quicklook_image
from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
from os.path import expanduser
os.environ['TERM'] = 'vt100'


class subtract_sky(object):
    """
    *Subtract the sky background using the Kelson Method*

    A model of the sky-background is created using a method similar to that described in Kelson, D. (2003), *Optimal Techniques in Two-dimensional Spectroscopy: Background Subtraction for the 21st Century (http://dx.doi.org/10.1086/375502). This model-background is then subtracted from the object spectrum.

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``objectFrame`` -- the image frame in need of sky subtraction
        - ``twoDMap`` -- 2D dispersion map image path
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``dispMap`` -- path to dispersion map. Default * False*
        - ``sofName`` -- name of the originating SOF file
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a subtract_sky object, use the following:

    ```eval_rst
    .. todo::

        - create cl-util for this class
        - add a tutorial about ``subtract_sky`` to documentation
    ```

    ```python
    from soxspipe.commonutils import subtract_sky
    skymodel = subtract_sky(
        log=log,
        settings=settings,
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
            twoDMap,
            objectFrame,
            qcTable,
            productsTable,
            dispMap=False,
            settings=False,
            sofName=False,
            recipeName="soxs-stare"
    ):
        self.log = log
        log.debug("instansiating a new 'subtract_sky' object")
        self.settings = settings
        self.objectFrame = objectFrame
        self.twoDMap = twoDMap
        self.dispMap = dispMap
        self.qc = qcTable
        self.products = productsTable
        self.sofName = sofName
        self.recipeName = recipeName

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

        # DATA FRAME CONTAINING ALL PIXELS - X, Y, FLUX, WAVELENGTH, SLIT-POSITION, ORDERE
        self.mapDF = twoD_disp_map_image_to_dataframe(log=self.log, slit_length=dp["slit_length"], twoDMapPath=twoDMap, assosiatedFrame=self.objectFrame)

        if self.inst == "SOXS":
            self.axisA = "y"
            self.axisB = "x"
        elif self.inst == "XSHOOTER":
            self.axisA = "x"
            self.axisB = "y"

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

        return None

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

        print(f'\n# MODELLING SKY BACKGROUND AND REMOVING FROM SCIENCE FRAME\n')

        skymodelCCDData, skySubtractedCCDData = self.create_placeholder_images()

        uniqueOrders = self.mapDF['order'].unique()
        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # BSPLINE ORDER TO FIT SKY WITH
        bspline_order = self.settings["sky-subtraction"]["bspline_order"]

        qcPlotOrder = int(np.median(uniqueOrders)) - 1
        allimageMapOrder = []
        allimageMapOrderWithObject = []
        for o in uniqueOrders:
            # SELECT ONLY ONE ORDER OF DATA
            imageMapOrder = self.mapDF[self.mapDF["order"] == o]
            imageMapOrderWithObject, imageMapOrder = self.get_over_sampled_sky_from_order(imageMapOrder, o, ignoreBP=False, clipSlitEdge=0.05)
            allimageMapOrder.append(imageMapOrder)
            allimageMapOrderWithObject.append(imageMapOrderWithObject)

        allimageMapOrder = self.clip_object_slit_positions(allimageMapOrder, aggressive=self.settings["sky-subtraction"]["aggressive_object_masking"])

        for o, imageMapOrder, imageMapOrderWithObject in zip(uniqueOrders, allimageMapOrder, allimageMapOrderWithObject):
            # if o != qcPlotOrder:
            #     continue

            imageMapOrder, tck = self.fit_csaps_curve(imageMapOrder, o, bspline_order)

            # CLIP EXTREMES OF ORDERS
            wlrange = imageMapOrder["wavelength"].max() - imageMapOrder["wavelength"].max()
            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            mask = (imageMapOrder["wavelength"] < imageMapOrder["wavelength"].max() - wlrange * 0.05) & (imageMapOrder["wavelength"] < imageMapOrder["wavelength"].min() + wlrange * 0.05)
            imageMapOrder.loc[mask, "sky_model"] = np.nan
            imageMapOrder.loc[mask, "sky_subtracted_flux"] = np.nan

            # tck = False
            # imageMapOrder, tck = self.fit_bspline_curve_to_sky(imageMapOrder, o, bspline_order)

            # residuals = imageMapOrder.loc[imageMapOrder["clipped"] == False]["sky_subtracted_flux"]
            # slit = imageMapOrder.loc[imageMapOrder["clipped"] == False]["slit_position"]

            # import matplotlib.pyplot as plt
            # plt.scatter(residuals, slit, s=0.5, alpha=0.1)
            # plt.xlim([-1000, 1000])
            # plt.show()

            # print(o)
            # imageMapOrder = self.fit_surface_to_sky_residuals(imageMapOrder)

            if isinstance(imageMapOrder, pd.core.frame.DataFrame):
                # REMOVE FILTERED ROWS FROM DATA FRAME
                # residuals = imageMapOrder.loc[imageMapOrder["clipped"] == False]["fit_sky_subtracted_flux"]
                # slit = imageMapOrder.loc[imageMapOrder["clipped"] == False]["slit_position"]

                # import matplotlib.pyplot as plt
                # plt.scatter(residuals, slit, s=0.5, alpha=0.1)
                # plt.xlim([-1000, 1000])
                # plt.show()

                # skyPixelsDF.append(imageMapOrder)

                # imageMapOrder = self.fit_surface_to_sky(imageMapOrder, o, bspline_order)
                skymodelCCDData, skySubtractedCCDData = self.add_data_to_placeholder_images(imageMapOrder, skymodelCCDData, skySubtractedCCDData)
                if o == qcPlotOrder:
                    qc_plot_path = self.plot_sky_sampling(order=o, imageMapOrderWithObjectDF=imageMapOrderWithObject, imageMapOrderDF=imageMapOrder, tck=tck)
                    basename = os.path.basename(qc_plot_path)
                    self.products = self.products.append({
                        "soxspipe_recipe": "soxs-stare",
                        "product_label": "SKY_MODEL_QC_PLOTS",
                        "file_name": basename,
                        "file_type": "PDF",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "product_desc": f"QC plots for the sky-background modelling",
                        "file_path": qc_plot_path
                    }, ignore_index=True)

        # plt.show()

        filename = self.filenameTemplate.replace(".fits", "_SKYMODEL.fits")
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home) + f"/product/{self.recipeName}"
        outDir = outDir.replace("//", "/")
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        filePath = f"{outDir}/{filename}"
        self.products = self.products.append({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "SKY_MODEL",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"The sky background model",
            "file_path": filePath
        }, ignore_index=True)

        # WRITE CCDDATA OBJECT TO FILE
        HDUList = skymodelCCDData.to_hdu(
            hdu_mask='QUAL', hdu_uncertainty='ERRS', hdu_flags=None)
        HDUList[0].name = "FLUX"
        HDUList.writeto(filePath, output_verify='exception',
                        overwrite=True, checksum=True)

        filename = self.filenameTemplate.replace(".fits", "_SKYSUB.fits")
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home) + f"/product/{self.recipeName}"
        filePath = f"{outDir}/{filename}"
        self.products = self.products.append({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "SKY_SUBTRACTED_OBJECT",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"The sky-subtracted object",
            "file_path": filePath
        }, ignore_index=True)

        # WRITE CCDDATA OBJECT TO FILE
        HDUList = skySubtractedCCDData.to_hdu(
            hdu_mask='QUAL', hdu_uncertainty='ERRS', hdu_flags=None)
        HDUList[0].name = "FLUX"
        HDUList.writeto(filePath, output_verify='exception',
                        overwrite=True, checksum=True)

        comparisonPdf = self.plot_image_comparison(self.objectFrame, skymodelCCDData, skySubtractedCCDData)

        filename = os.path.basename(comparisonPdf)
        self.products = self.products.append({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "SKY SUBTRACTION QUICKLOOK",
            "file_name": filename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"Sky-subtraction quicklook",
            "file_path": comparisonPdf
        }, ignore_index=True)

        self.log.debug('completed the ``get`` method')
        return skymodelCCDData, skySubtractedCCDData, self.qc, self.products

    def get_over_sampled_sky_from_order(
            self,
            imageMapOrder,
            order,
            ignoreBP=True,
            clipSlitEdge=False):
        """*unpack the over sampled sky from an order*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe from object image and 2D map
            - ``order`` -- the order number
            - ``ignoreBP`` -- ignore bad-pixels? Deafult *True*
            - ``clipSlitEdge`` -- clip the slit edges. Percentage of slit width to clip

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
        self.log.debug('starting the ``get_over_sampled_sky_from_order`` method')

        # MEDIAN CLIPPING FIRST USED TO CLIP MOST DEVIANT PIXELS (BAD AND CRHs)
        median_clipping_sigma = self.settings["sky-subtraction"]["median_clipping_sigma"]
        median_clipping_iterations = self.settings["sky-subtraction"]["median_clipping_iterations"]
        # ROLLING WINDOW LENGTH IN DATA POINTS
        median_rolling_window_size = self.settings["sky-subtraction"]["median_rolling_window_size"]
        # PECENTILE CLIPPING USED TO CLIP THE OBJECT(S) BEFORE FITTING A SKY MODEL
        percential_clipping_sigma = self.settings["sky-subtraction"]["percential_clipping_sigma"]
        percential_clipping_iterations = self.settings["sky-subtraction"]["percential_clipping_iterations"]
        percential_rolling_window_size = self.settings["sky-subtraction"]["percential_rolling_window_size"]

        # FINDING A DYNAMIC SIZE FOR PERCENTILE FILTERING WINDOW
        windowSize = int(len(imageMapOrder.loc[imageMapOrder["y"] == imageMapOrder["y"].median()].index))

        imageMapOrder["clipped"] = False
        imageMapOrder["object"] = False

        if clipSlitEdge:
            slitRange = imageMapOrder["slit_position"].max() - imageMapOrder["slit_position"].min()
            clipSlitEdge *= slitRange
            mask = ((imageMapOrder['slit_position'] > imageMapOrder["slit_position"].max() - clipSlitEdge) | (imageMapOrder['slit_position'] < imageMapOrder["slit_position"].min() + clipSlitEdge))
            imageMapOrder.loc[mask, "clipped"] = True

        if not ignoreBP:
            # REMOVE FILTERED ROWS FROM DATA FRAME
            mask = (imageMapOrder['mask'] == True)
            imageMapOrder.drop(index=imageMapOrder[mask].index, inplace=True)

         # CLIP -VE FLUX
        # mask = (imageMapOrder['flux'] < 0)
        # imageMapOrder.loc[mask, "clipped"] = True

        # CLIP THE MOST DEVIANT PIXELS WITHIN A WAVELENGTH ROLLING WINDOW - BAD-PIXELS AND CRHs
        # print(f'# Clipping extremely deviant pixels via a rolling wavelength window (typically bad-pixels and CRHs)')
        imageMapOrder = self.rolling_window_clipping(imageMapOrderDF=imageMapOrder, windowSize=int(median_rolling_window_size), sigma_clip_limit=median_clipping_sigma, max_iterations=median_clipping_iterations, median_centre_func=True)
        imageMapOrderWithObject = imageMapOrder.copy()

        # NOW SOME MORE ROBUST CLIPPING WITHIN A WAVELENGTH ROLLING WINDOW TO ALSO REMOVE OBJECT(S)
        # print(f'# Robustly clipping deviant pixels via a rolling wavelength window (now including object(s))')
        imageMapOrder["residual_global_sigma_old"] = imageMapOrder["residual_global_sigma"]
        imageMapOrder = self.rolling_window_clipping(imageMapOrderDF=imageMapOrder, windowSize=int(percential_rolling_window_size), sigma_clip_limit=percential_clipping_sigma, max_iterations=percential_clipping_iterations)

        self.log.debug('completed the ``get_over_sampled_sky_from_order`` method')
        return imageMapOrderWithObject, imageMapOrder

    def plot_sky_sampling(
            self,
            order,
            imageMapOrderWithObjectDF,
            imageMapOrderDF,
            tck=False):
        """*generate a plot of sky sampling*

        **Key Arguments:**
            - ``order`` -- the order number.
            - ``imageMapOrderWithObjectDF`` -- dataframe with various processed data without object clipped
            - ``imageMapOrderDF`` -- dataframe with various processed data for order
            - ``tck`` -- spline parameters to replot

        **Return:**
            - ``filePath`` -- path to the generated QC plots PDF

        **Usage:**

        ```python
        self.plot_sky_sampling(
            order=myOrder,
            imageMapOrderWithObjectDF=imageMapOrderWithObject,
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

        # SET COLOURS FOR VARIOUS STAGES
        medianColor = "#C21820"
        percentileColor = "#C21820"
        skyColor = "purple"
        rawColor = "#93a1a1"
        rawColor = "#002b36"
        # SET PLOT LAYER ORDERS
        medianZ = 3
        skyZ = 3
        percentileZ = 2
        unclippedZ = 1
        # SET MARKER SIZES
        rawMS = 0.5
        medianMS = 3
        percentileMS = 3

        # MAKE A COPY OF THE FRAME TO NOT ALTER ORIGINAL DATA
        frame = self.objectFrame.copy()

        # SETUP THE PLOT SUB-PANELS
        fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=320)
        # REMOVE ME
        # fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=100)
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
        for x, y in zip(imageMapOrderWithObjectDF["x"], imageMapOrderWithObjectDF["y"]):
            nonOrderMask[y][x] = 0

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        nonOrderMask = ma.make_mask(nonOrderMask)
        combinedMask = (nonOrderMask == 1) | (frame.mask == 1)
        frame.mask = (nonOrderMask == 1)

        # RAW IMAGE PANEL
        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.flipud(np.rot90(frame, 1))
        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        dataArray = np.asarray(frame)
        maskedDataArray = np.ma.array(frame.data, mask=combinedMask)
        std = np.nanstd(maskedDataArray)
        mean = np.nanmean(maskedDataArray)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        im = onerow.imshow(rotatedImg, vmin=0, vmax=100, cmap='gray', alpha=1)
        medianValue = np.median(rotatedImg.data.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="unprocessed frame")]
        onerow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        onerow.set_xlabel(
            "y-axis", fontsize=10)
        onerow.set_ylabel(
            "x-axis", fontsize=10)
        ylimMinImage = imageMapOrderWithObjectDF["y"].min() - 10
        ylimMaxImage = imageMapOrderWithObjectDF["y"].max() + 10
        onerow.set_ylim(imageMapOrderWithObjectDF["x"].min() - 10, imageMapOrderWithObjectDF["x"].max() + 10)
        onerow.set_xlim(ylimMinImage, ylimMaxImage)
        onerow.invert_xaxis()

        # ORIGINAL DATA AND PERCENTILE SMOOTHED WAVELENGTH VS FLUX
        raw = tworow.plot(
            imageMapOrderWithObjectDF["wavelength"].values,
            imageMapOrderWithObjectDF["flux"].values, label='unprocessed (unp)', c=rawColor, zorder=0)
        # RAW MARKERS
        tworow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].values, label='unclipped', s=rawMS, c=rawColor, alpha=1., zorder=unclippedZ)
        # ROBUSTLY CLIPPED
        tworow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "flux"].values, label='clipped', s=percentileMS, marker="x", c=percentileColor, zorder=percentileZ, alpha=.5)
        # MEDIAN CLIPPED
        # label='median clipped'
        label = None
        tworow.scatter(
            imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "wavelength"].values,
            imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "flux"].values, label=label, s=medianMS, marker="x", c=medianColor, zorder=medianZ, alpha=.5)
        # MEDIAN LINE
        # label='median-smoothed (ms)'
        label = None
        # tworow.plot(
        #     imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "wavelength"].values,
        #     imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "flux_smoothed"].values, label=label, c=medianColor, zorder=medianZ)
        # PERCENTILE LINE
        tworow.plot(
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False), "wavelength"].values,
            imageMapOrderDF.loc[(imageMapOrderDF["clipped"] == False) & (imageMapOrderDF["object"] == False), "flux_smoothed"].values, label='percentile-smoothed', c="#C21820", zorder=percentileZ)

        # SIGMA RESIDUAL
        weights = tworow.plot(
            imageMapOrderWithObjectDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderWithObjectDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"].values * 5 - imageMapOrderWithObjectDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"].max() * 1.2, label='5$\sigma$ residual scatter (shifted)', c="#002b36")
        ylimmin = -imageMapOrderWithObjectDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"].max() * 1.3
        if ylimmin < -3000:
            ylimmin = -300
        tworow.set_ylim(ylimmin, imageMapOrderWithObjectDF["flux_smoothed"].max() * 1.2)

        tworow.set_ylabel(
            "counts", fontsize=10)
        tworow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.)
        # REMOVE ME
        # tworow.set_xlim(1504, 1507)
        tworow.set_xticks([], [])

        # WAVELENGTH RESIDUAL PANEL
        threerow.scatter(imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "wavelength"].values, imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "residual_global_sigma"].values, s=rawMS, alpha=1., c=rawColor, zorder=unclippedZ)
        threerow.scatter(imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "wavelength"].values, imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "residual_global_sigma"].values, s=medianMS, marker="x", c=medianColor, zorder=medianZ, alpha=0.5)
        threerow.scatter(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "wavelength"].values, imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "residual_global_sigma_old"].values, s=percentileMS, marker="x", c=percentileColor, zorder=percentileZ, alpha=0.5)
        std = imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "residual_global_sigma"].std()
        median = imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "residual_global_sigma"].median()
        threerow.set_ylim(median - 3 * std, median + 7 * std)
        threerow.set_xlabel(
            "wavelength (nm)", fontsize=10)
        # REMOVE ME
        # threerow.set_xlim(1504, 1507)
        threerow.set_ylabel("residual ($\sigma$)", fontsize=10)

        # SLIT-POSITION RESIDUAL PANEL (SHOWING OBJECT)
        fourrow.scatter(
            imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "slit_position"].values,
            imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "residual_global_sigma"].values, label='deviations', s=medianMS, marker="x", c=medianColor, zorder=medianZ, alpha=0.2)
        fourrow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "slit_position"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "residual_global_sigma_old"].values, label='deviations', s=percentileMS, marker="x", c=percentileColor, zorder=percentileZ, alpha=0.2)
        fourrow.set_ylim(median - 3 * std, median + 7 * std)
        fourrow.set_xlabel(
            "slit-position relative to slit centre (arcsec)", fontsize=10)
        fourrow.set_ylabel("residual ($\sigma$)", fontsize=10)
        fourrow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "slit_position"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_global_sigma_old"].values, label='deviations', s=rawMS, alpha=0.5, c=rawColor, zorder=unclippedZ)

        # IMAGE SHOWING CLIPPED PIXEL MASK
        im = fiverow.imshow(rotatedImg, vmin=0, vmax=100, cmap='gray', alpha=1)
        # medianValue = np.median(rotatedImg.data.ravel())
        # color = im.cmap(im.norm(medianValue))
        # patches = [mpatches.Patch(color=color, label="clipped frame")]
        # fiverow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        percentileClipMask = nonOrderMask
        percentileClipMask = np.zeros_like(frame.data)
        for x, y in zip(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "x"].values, imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "y"].values):
            percentileClipMask[y][x] = 1
        percentileClipMask = ma.make_mask(percentileClipMask)
        imageMask = np.ma.array(np.ones_like(frame.data), mask=~percentileClipMask)
        # MAKE A COLOR MAP OF FIXED COLORS
        cmap = colors.ListedColormap([percentileColor, percentileColor])
        bounds = [0, 5, 10]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        cmap.set_bad(medianColor, 0.)
        fiverow.imshow(np.flipud(np.rot90(imageMask, 1)), cmap=cmap, norm=norm, alpha=1., interpolation='nearest')
        medianClipMask = np.zeros_like(frame.data)
        for x, y in zip(imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "x"].values, imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "y"].values):
            medianClipMask[y][x] = 1
        medianClipMask = ma.make_mask(medianClipMask)
        imageMask = np.ma.array(np.ones_like(frame.data), mask=~medianClipMask)
        # MAKE A COLOR MAP OF FIXED COLORS
        cmap = colors.ListedColormap([medianColor, medianColor])
        bounds = [0, 5, 10]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        cmap.set_bad(medianColor, 0.)
        fiverow.imshow(np.flipud(np.rot90(imageMask, 1)), cmap=cmap, norm=norm, alpha=1., interpolation='nearest')

        nonOrderMask = (nonOrderMask == 0)
        imageMask = np.ma.array(np.ones_like(frame.data), mask=nonOrderMask)
        cmap = copy(cm.gray)
        cmap.set_bad("green", 0.0)
        fiverow.imshow(np.flipud(np.rot90(imageMask, 1)), vmin=-10, vmax=-9, cmap=cmap, alpha=1.)
        fiverow.set_xlabel(
            "y-axis", fontsize=10)
        fiverow.set_ylabel(
            "x-axis", fontsize=10)
        fiverow.set_ylim(imageMapOrderWithObjectDF["x"].min() - 10, imageMapOrderWithObjectDF["x"].max() + 10)
        fiverow.set_xlim(ylimMinImage, ylimMaxImage)
        # REMOVE ME
        # fiverow.set_xlim(1510, 1600)
        # fiverow.set_ylim(570, 630)
        fiverow.invert_xaxis()

        # PLOT WAVELENGTH VS FLUX SKY MODEL
        sixrow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].values, label='unclipped', s=rawMS, c=rawColor, alpha=0.5, zorder=unclippedZ)
        if tck:
            wl = np.linspace(imageMapOrderDF["wavelength"].min(), imageMapOrderDF["wavelength"].max(), 1000000)
            sky = ip.splev(wl, tck)
            skymodel = sixrow.plot(
                wl, sky, label='sky model', c=skyColor, zorder=skyZ)
        else:
            skymodel = sixrow.plot(
                imageMapOrderDF["wavelength"].values,
                imageMapOrderDF["sky_model"].values, label='sky model', c=skyColor, zorder=skyZ)
        if ylimmin < -3000:
            ylimmin = -300
        sixrow.set_ylim(ylimmin, imageMapOrderWithObjectDF["flux_smoothed"].max() * 1.2)
        # REMOVE ME
        # sixrow.set_xlim(1504, 1507)
        sixrow.set_ylabel(
            "counts", fontsize=10)
        sixrow.legend(loc=2, fontsize=8, bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        # BUILD IMAGE OF SKY MODEL
        skyModelImage = np.zeros_like(frame.data)
        for x, y, skypixel in zip(imageMapOrderDF["x"], imageMapOrderDF["y"], imageMapOrderDF["sky_model"]):
            skyModelImage[y][x] = skypixel
        nonOrderMask = (nonOrderMask == 0)
        skyModelImage = np.ma.array(skyModelImage, mask=nonOrderMask)
        cmap = copy(cm.gray)
        std = np.nanstd(skyModelImage)
        mean = np.nanmean(skyModelImage)
        vmax = mean + 2 * std
        vmin = mean - 1 * std
        im = sevenrow.imshow(np.flipud(np.rot90(skyModelImage, 1)), vmin=0, vmax=100, cmap=cmap, alpha=1.)
        # sevenrow.set_xlabel(
        #     "y-axis", fontsize=10)
        sevenrow.set_ylabel(
            "x-axis", fontsize=10)
        sevenrow.set_ylim(imageMapOrderWithObjectDF["x"].min() - 10, imageMapOrderWithObjectDF["x"].max() + 10)
        sevenrow.set_xlim(ylimMinImage, ylimMaxImage)
        # REMOVE ME
        # sevenrow.set_xlim(1510, 1600)
        # sevenrow.set_ylim(570, 630)
        sevenrow.invert_xaxis()
        medianValue = np.median(skyModelImage.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="sky model")]
        sevenrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        sevenrow.set_xticks([], [])

        # BUILD SKY-SUBTRACTED IMAGE
        skySubImage = np.zeros_like(frame.data)
        for x, y, skypixel in zip(imageMapOrderDF["x"], imageMapOrderDF["y"], imageMapOrderDF["sky_subtracted_flux"]):
            skySubImage[y][x] = skypixel
        skySubMask = (nonOrderMask == 1) | (medianClipMask == 1)
        skySubImage = np.ma.array(skySubImage, mask=skySubMask)
        cmap = copy(cm.gray)
        std = np.nanstd(skySubImage)
        mean = np.nanmedian(skySubImage)
        vmax = mean + 0.2 * std
        vmin = mean - 0.2 * std
        im = eightrow.imshow(np.flipud(np.rot90(skySubImage, 1)), vmin=0, vmax=50, cmap=cmap, alpha=1.)
        eightrow.set_xlabel(
            "y-axis", fontsize=10)
        eightrow.set_ylabel(
            "x-axis", fontsize=10)
        eightrow.set_ylim(imageMapOrderWithObjectDF["x"].min() - 10, imageMapOrderWithObjectDF["x"].max() + 10)
        eightrow.set_xlim(ylimMinImage, ylimMaxImage)
        # REMOVE ME
        # eightrow.set_xlim(1510, 1600)
        # eightrow.set_ylim(570, 630)
        eightrow.invert_xaxis()
        medianValue = np.median(skySubImage.data.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="sky-subtracted frame")]
        eightrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        # SUBTRACTED SKY RESIDUAL PANEL
        ninerow.scatter(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values, np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux"].values), s=rawMS, alpha=1., c=rawColor, zorder=unclippedZ)
        ninerow.scatter(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "wavelength"].values, np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "sky_subtracted_flux"].values), s=medianMS, marker="x", c=medianColor, zorder=medianZ, alpha=.5)
        # ninerow.scatter(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "wavelength"].values, imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "residual_global_sigma_old"].values, s=percentileMS, marker="x", c=percentileColor, zorder=percentileZ, alpha=1.)
        # std = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux"].std()

        mean = np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux"]).mean()
        std = np.absolute(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "sky_subtracted_flux"]).std()

        ninerow.set_ylim(-100, 1000)
        ninerow.set_xlabel(
            "wavelength (nm)", fontsize=10)
        # REMOVE ME
        # ninerow.set_xlim(1504, 1507)
        ninerow.set_ylabel("residual", fontsize=10)

        fig.suptitle(f"{self.arm} sky model: order {order}", fontsize=12, y=0.97)

        filename = self.filenameTemplate.replace(".fits", f"_SKYMODEL_QC_PLOTS_ORDER_{int(order)}.pdf")
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home) + "/qc/pdf"
        filePath = f"{outDir}/{filename}"
        # REMOVE ME
        # plt.show()
        plt.savefig(filePath, dpi='figure')

        self.log.debug('completed the ``plot_sky_sampling`` method')
        return filePath

    def rolling_window_clipping(
            self,
            imageMapOrderDF,
            windowSize,
            sigma_clip_limit=5,
            max_iterations=10,
            median_centre_func=False):
        """*clip pixels in a rolling wavelength window*

        **Key Arguments:**
            - ``imageMapOrderDF`` --  dataframe with various processed data for a given order
            - ``windowSize`` -- the window-size used to perform rolling window clipping (number of data-points)
            - ``sigma_clip_limit`` -- clip data values straying beyond this sigma limit. Default *5*
            - ``max_iterations`` -- maximum number of iterations when clipping
            - ``median_centre_func`` -- use a median centre function for rolling window instead of quantile (use to clip most deviate pixels only). Default *False*

        **Return:**
            - ``imageMapOrderDF`` -- image order dataframe with 'clipped' == True for those pixels that have been clipped via rolling window clipping

        **Usage:**

        ```python
        imageMapOrder = self.rolling_window_clipping(
            imageMapOrderDF=imageMapOrder,
            windowSize=23,
            sigma_clip_limit=4,
            max_iterations=10,
            median_centre_func=True
        )
        ```
        """
        self.log.debug('starting the ``rolling_window_clipping`` method')

        i = 1
        newlyClipped = -1
        allPixels = len(imageMapOrderDF.index)

        while i <= max_iterations:

            # CALCULATE PERCENTILE SMOOTH DATA & RESIDUALS
            if median_centre_func:
                imageMapOrderDF["flux_smoothed"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].rolling(window=windowSize, center=True).median()
            elif i == 1:
                imageMapOrderDF["flux_smoothed"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].rolling(window=windowSize, center=True).quantile(.40)

            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"] - imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_smoothed"]
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"].rolling(windowSize).std()
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_sigma"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"] / imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"]

            if median_centre_func:
                # CLIP ABOVE AND BELOW
                imageMapOrderDF.loc[imageMapOrderDF["residual_windowed_sigma"].abs() > sigma_clip_limit, "clipped"] = True
            else:
                # CLIP ONLY HIGH VALUES
                # imageMapOrderDF.loc[((imageMapOrderDF["residual_windowed_sigma"] > sigma_clip_limit) & (imageMapOrderDF["residual_global_sigma_old"] > -10.0)), "clipped"] = True
                imageMapOrderDF.loc[((imageMapOrderDF["residual_windowed_sigma"] > sigma_clip_limit) & (imageMapOrderDF["residual_global_sigma_old"] > -10.0)), "object"] = True

            if newlyClipped == -1:
                totalClipped = len(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True].index)
                newlyClipped = totalClipped
            else:
                newlyClipped = len(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True].index) - totalClipped
                totalClipped = len(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True].index)

            # Cursor up one line and clear line
            sys.stdout.write("\x1b[1A\x1b[2K")
            percent = (float(totalClipped) / float(allPixels)) * 100.
            print(f'\tITERATION {i}: {newlyClipped} deviant pixels have been newly clipped within a {windowSize} data-point rolling window ({totalClipped} pixels clipped in total = {percent:1.1f}%)')
            if newlyClipped == 0:
                break
            i += 1

        imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_scatter_windowed_std"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].rolling(windowSize).std()
        std = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"].std()
        imageMapOrderDF["residual_global_sigma"] = imageMapOrderDF["flux_minus_smoothed_residual"] / std

        self.log.debug('completed the ``rolling_window_clipping`` method')
        return imageMapOrderDF

    def fit_bspline_curve_to_sky(
            self,
            imageMapOrder,
            order,
            bspline_order):
        """*fit a bspline to the unclipped sky pixels (wavelength vs flux)*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe, containing sky flux with object(s) and CRHs removed
            - ``order`` -- the order number

        **Return:**
            - ``imageMapOrder`` -- same `imageMapOrder` as input but now with `sky_model` (bspline fit of the sky) and `sky_subtracted_flux` columns

        **Usage:**

        ```python
        imageMapOrder = self.fit_bspline_curve_to_sky(
            imageMapOrder,
            myOrder
        )
        ```

        """
        self.log.debug('starting the ``fit_bspline_curve_to_sky`` method')

        import numpy as np
        import scipy.interpolate as ip

        # VARIABLES TO PLAY WITH
        # 1. WEIGHTS ... IF USING WEIGHTS YOU NEED s=0. If the errors in the y values have standard-deviation given by the vector d, then w should be 1/d. Default is ones(len(x)). Note weights are relative to one another (don't have to by between 0-1)
        # 2. KNOT PLACEMENT. The knots needed for task=-1. If given then task is automatically set to -1.
        # 3. DEGREE OF FIT. It is recommended to use cubic splines. Even values of k should be avoided especially with small s values. 1 <= k <= 5
        # 4. TASK: If task==0 find t and c for a given smoothing factor, s. If task=-1 find the weighted least square spline for a given set of knots, t. These should be interior knots as knots on the ends will be added automatically.
        # 5. Smoothing.  If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range (m-sqrt(2*m),m+sqrt(2*m)) where m is the number of datapoints in x, y, and w. default : s=m-sqrt(2*m) if weights are supplied. s = 0.0 (interpolating) if no weights are supplied.

        # RETURNS
        # tck: A tuple (t,c,k) containing the vector of knots, the B-spline coefficients, and the degree of the spline.
        # fp: The weighted sum of squared residuals of the spline approximation.

        # REMOVE ME
        bspline_order = 3

        # SORT BY COLUMN NAME
        df = imageMapOrder.copy()
        df.sort_values(by=['wavelength'], inplace=True)
        # df.drop_duplicates(subset=['wavelength'], inplace=True)

        goodWl = df.loc[df["clipped"] == False]["wavelength"]
        goodFlux = df.loc[df["clipped"] == False]["flux"]
        df["weights"] = 1 / df["error"]
        df["weights"] = df["weights"].replace(np.nan, 0.00000001)

        goodWeights = df.loc[df["clipped"] == False, "weights"]

        # N = int(goodWl.shape[0] / 25)
        # seedKnots = np.linspace(xmin, xmax, N)

        # t for knots
        # c of coefficients
        # k for order

        # FOR WEIGHTED BSPLINES WE ONLY NEED *INTERIOR* KNOT (DON'T GO BEYOND RANGE OF DATA)
        # CAN'T HAVE MORE KNOTS THAN DATA POINTS
        n_interior_knots = int(goodWl.values.shape[0] / 2)
        qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
        # qs = np.append(np.array([0, 0, 0]), qs)
        # qs = np.append(qs, np.array([1, 1, 1]))

        print(qs.min(), qs.max(), qs.mean(), qs.shape[0])

        knots = np.quantile(goodWl, qs)
        print(knots)

        print(qs.shape[0])
        print(knots.min(), knots.max(), knots.mean(), knots.shape[0])
        print(goodWl.min(), goodWl.max(), goodWl.mean(), goodWl.shape[0])

        # tck = ip.splrep(goodWl, goodFlux, t=knots, k=3)

        tck, fp, ier, msg = ip.splrep(goodWl, goodFlux, t=knots, k=bspline_order, w=goodWeights, full_output=True)
        print("GOODNESS:")
        print(fp / 4163287709537)

        # sys.exit(0)
        t, c, k = tck
        # sky_model = ip.BSpline(t, c, k)(imageMapOrder["wavelength"].values)
        sky_model = ip.splev(imageMapOrder["wavelength"].values, tck)

        imageMapOrder["sky_model"] = sky_model

        # t, c, k = splrep(goodWl, goodFlux, t=seedKnots[1:-1], w=goodWeights, s=0.0, k=rowFitOrder, task=-1)
        # spline = BSpline(t, c, k, extrapolate=True)

        # t for knots
        # c of coefficients
        # k for order
        # t, c, k = splrep(goodWl, goodFlux, w=goodWeights, s=0.0, k=rowFitOrder)
        # spline = BSpline(t, c, k, extrapolate=True)

        # spl = splrep(goodWl, goodFlux)
        # imageMapOrder["sky_model"] = splev(imageMapOrder["wavelength"].values, spl)

        # t, c, k = ip.splrep(goodWl, goodFlux, s=0.0, k=bspline_order)
        # print(t)
        # print(len(t))
        # print(len(goodWl))
        # spline = ip.BSpline(t, c, k, extrapolate=False)

        # imageMapOrder["sky_model"] = spline(imageMapOrder["wavelength"].values)
        imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model"]

        self.log.debug('completed the ``fit_bspline_curve_to_sky`` method')
        return imageMapOrder, tck

    def fit_csaps_curve(
            self,
            imageMapOrder,
            order,
            bspline_order):
        """*fit a bspline to the unclipped sky pixels (wavelength vs flux)*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe, containing sky flux with object(s) and CRHs removed
            - ``order`` -- the order number

        **Return:**
            - ``imageMapOrder`` -- same `imageMapOrder` as input but now with `sky_model` (bspline fit of the sky) and `sky_subtracted_flux` columns

        **Usage:**

        ```python
        imageMapOrder = self.fit_bspline_curve_to_sky(
            imageMapOrder,
            myOrder
        )
        ```

        """
        self.log.debug('starting the ``fit_bspline_curve_to_sky`` method')

        # VARIABLES TO PLAY WITH
        # 1. WEIGHTS ... IF USING WEIGHTS YOU NEED s=0. If the errors in the y values have standard-deviation given by the vector d, then w should be 1/d. Default is ones(len(x)). Note weights are relative to one another (don't have to by between 0-1)
        # 2. KNOT PLACEMENT. The knots needed for task=-1. If given then task is automatically set to -1.
        # 3. DEGREE OF FIT. It is recommended to use cubic splines. Even values of k should be avoided especially with small s values. 1 <= k <= 5
        # 4. TASK: If task==0 find t and c for a given smoothing factor, s. If task=-1 find the weighted least square spline for a given set of knots, t. These should be interior knots as knots on the ends will be added automatically.
        # 5. Smoothing.  If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range (m-sqrt(2*m),m+sqrt(2*m)) where m is the number of datapoints in x, y, and w. default : s=m-sqrt(2*m) if weights are supplied. s = 0.0 (interpolating) if no weights are supplied.

        # RETURNS
        # tck: A tuple (t,c,k) containing the vector of knots, the B-spline coefficients, and the degree of the spline.
        # fp: The weighted sum of squared residuals of the spline approximation.

        import numpy as np
        import matplotlib.pyplot as plt
        from csaps import csaps
        from astropy.stats import sigma_clip

        # SORT BY COLUMN NAME
        df = imageMapOrder.copy()
        df.sort_values(by=['wavelength'], inplace=True)
        df.drop_duplicates(subset=['wavelength'], inplace=True)

        goodWl = df.loc[df["clipped"] == False]["wavelength"]
        goodFlux = df.loc[df["clipped"] == False]["flux"]
        goodSlit = df.loc[df["clipped"] == False]["slit_position"]
        df["weights"] = 1 / df["error"]
        df["weights"] = df["weights"].replace(np.nan, 0.00000001)
        goodWeights = df.loc[df["clipped"] == False, "weights"]
        median_rolling_window_size = int(self.settings["sky-subtraction"]["median_rolling_window_size"])

        fittingDF = df.copy()
        mask = (df["clipped"] == False)
        fittingDF = df.loc[mask]
        smooth = 0.99999999999

        iteration = 0
        while iteration < 3:
            iteration += 1
            fittingDF["residuals"] = csaps(fittingDF["wavelength"].values, fittingDF["flux"].values, fittingDF["wavelength"].values, weights=fittingDF["weights"].values, smooth=smooth)
            fittingDF["residuals"] = fittingDF["residuals"].abs()

            # print(len(fittingDF.index))
            # print(fittingDF["residuals"].max())
            # print(fittingDF["residuals"].mean())
            # print(fittingDF["residuals"].min())
            # print()
            # SIGMA-CLIP THE DATA
            fittingDF["residuals_smoothed"] = fittingDF["residuals"].rolling(window=median_rolling_window_size, center=True).median()
            fittingDF["residuals_minus_smoothed_residual"] = fittingDF["residuals"] - fittingDF["residuals_smoothed"]
            std = fittingDF["residuals_minus_smoothed_residual"].std()
            fittingDF["residual_sigma"] = fittingDF["residuals_minus_smoothed_residual"] / std
            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            mask = (fittingDF["residual_sigma"] < 5.0)
            fittingDF = fittingDF.loc[mask]

        fittingDF["residuals"] = csaps(fittingDF["wavelength"].values, fittingDF["flux"].values, fittingDF["wavelength"].values, weights=fittingDF["weights"].values, smooth=smooth)

        sky_model = csaps(fittingDF["wavelength"].values, fittingDF["flux"].values, imageMapOrder["wavelength"].values, weights=fittingDF["weights"].values, smooth=smooth)

        imageMapOrder["sky_model"] = sky_model
        imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model"]

        self.log.debug('completed the ``fit_bspline_curve_to_sky`` method')
        return imageMapOrder, False

    def fit_bspline_surface_to_sky(
            self,
            imageMapOrder,
            order,
            bspline_order):
        """*fit a bspline to the unclipped sky pixels (wavelength vs flux)*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe, containing sky flux with object(s) and CRHs removed
            - ``order`` -- the order number

        **Return:**
            - ``imageMapOrder`` -- same `imageMapOrder` as input but now with `sky_model` (bspline fit of the sky) and `sky_subtracted_flux` columns

        **Usage:**

        ```python
        imageMapOrder = self.fit_bspline_curve_to_sky(
            imageMapOrder,
            myOrder
        )
        ```

        """
        self.log.debug('starting the ``fit_bspline_curve_to_sky`` method')

        # VARIABLES TO PLAY WITH
        # 1. WEIGHTS ... IF USING WEIGHTS YOU NEED s=0. If the errors in the y values have standard-deviation given by the vector d, then w should be 1/d. Default is ones(len(x)). Note weights are relative to one another (don't have to by between 0-1)
        # 2. KNOT PLACEMENT. The knots needed for task=-1. If given then task is automatically set to -1.
        # 3. DEGREE OF FIT. It is recommended to use cubic splines. Even values of k should be avoided especially with small s values. 1 <= k <= 5
        # 4. TASK: If task==0 find t and c for a given smoothing factor, s. If task=-1 find the weighted least square spline for a given set of knots, t. These should be interior knots as knots on the ends will be added automatically.
        # 5. Smoothing.  If the weights represent the inverse of the standard-deviation of y, then a good s value should be found in the range (m-sqrt(2*m),m+sqrt(2*m)) where m is the number of datapoints in x, y, and w. default : s=m-sqrt(2*m) if weights are supplied. s = 0.0 (interpolating) if no weights are supplied.

        # RETURNS
        # tck: A tuple (t,c,k) containing the vector of knots, the B-spline coefficients, and the degree of the spline.
        # fp: The weighted sum of squared residuals of the spline approximation.

        import numpy as np
        import matplotlib.pyplot as plt
        import scipy

        # SORT BY COLUMN NAME
        df = imageMapOrder.copy()
        df.sort_values(by=['wavelength'], inplace=True)
        df.drop_duplicates(subset=['wavelength'], inplace=True)

        goodWl = df.loc[df["clipped"] == False]["wavelength"]
        goodFlux = df.loc[df["clipped"] == False]["flux"]
        goodSlit = df.loc[df["clipped"] == False]["slit_position"]
        # goodWl = df["wavelength"]
        # goodFlux = df["flux"]
        # goodSlit = df["slit_position"]
        df["weights"] = 1 / df["error"]
        df["weights"] = df["weights"].replace(np.nan, 0.00000001)
        goodWeights = df.loc[df["clipped"] == False, "weights"]

        print(goodWl[:100])
        print(goodSlit[:100])
        print(goodFlux[:100])

        from astropy.stats import sigma_clip

        iteration = 0
        while iteration < 3:
            iteration += 1
            interp_surface = scipy.interpolate.SmoothBivariateSpline(goodWl, goodSlit, goodFlux, ky=2, kx=5, s=100000000000)
            tmp = interp_surface(goodWl, goodSlit, grid=False)
            diff = goodFlux - tmp
            print(f"MEAN diff: {diff.mean()}")
            # SIGMA-CLIP THE DATA
            # masked_diff = sigma_clip(diff, sigma_lower=3, sigma_upper=3, maxiters=1, cenfunc='median', stdfunc='mad_std')
            # REDUCE ARRAYS TO NON-MASKED VALUES
            # aa = [goodWl, goodSlit, goodFlux]
            # goodWl, goodSlit, goodFlux = [np.ma.compressed(np.ma.masked_array(i, masked_diff.mask)) for i in aa]

        test_x = np.arange(df['wavelength'].min() + 5, df['wavelength'].max() - 5, 0.1)
        test_y = np.arange(df['slit_position'].min() + 1, df['slit_position'].max() - 1, 0.1)

        surface = interp_surface(test_x, test_y)
        surface = np.rot90(surface)

        fig, axes = plt.subplots(1, 1, figsize=(16, 2))
        extent = [test_x[0], test_x[-1], test_y[0], test_y[-1]]

        std = surface.std()
        mean = surface.mean()

        im = axes.imshow(surface, aspect='auto', cmap='nipy_spectral', extent=extent, vmin=mean - 2 * std, vmax=mean + 2 * std)
        fig.colorbar(im, ax=axes)
        axes.plot(goodWl, goodSlit, 'w.', ms=0.8, alpha=0.4)

        plt.show()

        # N = int(goodWl.shape[0] / 25)
        # seedKnots = np.linspace(xmin, xmax, N)

        # t for knots
        # c of coefficients
        # k for order

        # FOR WEIGHTED BSPLINES WE ONLY NEED *INTERIOR* KNOT (DON'T GO BEYOND RANGE OF DATA)
        # CAN'T HAVE MORE KNOTS THAN DATA POINTS
        n_interior_knots = int(goodWl.values.shape[0] / 2)
        qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
        # qs = np.append(np.array([0, 0, 0]), qs)
        # qs = np.append(qs, np.array([1, 1, 1]))

        print(qs.min(), qs.max(), qs.mean(), qs.shape[0])

        knots = np.quantile(goodWl, qs)
        print(knots)

        print(qs.shape[0])
        print(knots.min(), knots.max(), knots.mean(), knots.shape[0])
        print(goodWl.min(), goodWl.max(), goodWl.mean(), goodWl.shape[0])

        # tck = ip.splrep(goodWl, goodFlux, t=knots, k=3)
        print("start")
        from csaps import csaps
        sky_model = csaps(goodWl, goodFlux, imageMapOrder["wavelength"].values, smooth=0.9999999)

        # tck, fp, ier, msg = ip.splrep(goodWl, goodFlux, t=knots, k=bspline_order, w=goodWeights, full_output=True)
        # print("GOODNESS:")
        # print(fp / 4163287709537)

        imageMapOrder["sky_model"] = sky_model

        # t, c, k = splrep(goodWl, goodFlux, t=seedKnots[1:-1], w=goodWeights, s=0.0, k=rowFitOrder, task=-1)
        # spline = BSpline(t, c, k, extrapolate=True)

        # t for knots
        # c of coefficients
        # k for order
        # t, c, k = splrep(goodWl, goodFlux, w=goodWeights, s=0.0, k=rowFitOrder)
        # spline = BSpline(t, c, k, extrapolate=True)

        # spl = splrep(goodWl, goodFlux)
        # imageMapOrder["sky_model"] = splev(imageMapOrder["wavelength"].values, spl)

        # t, c, k = ip.splrep(goodWl, goodFlux, s=0.0, k=bspline_order)
        # print(t)
        # print(len(t))
        # print(len(goodWl))
        # spline = ip.BSpline(t, c, k, extrapolate=False)

        # imageMapOrder["sky_model"] = spline(imageMapOrder["wavelength"].values)
        imageMapOrder["sky_subtracted_flux"] = imageMapOrder["flux"] - imageMapOrder["sky_model"]

        self.log.debug('completed the ``fit_bspline_curve_to_sky`` method')
        return imageMapOrder, False

    def create_placeholder_images(
            self):
        """*create placeholder images for the sky model and sky-subtracted frame*

        **Key Arguments:**
            # -

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

        self.log.debug('completed the ``create_placeholder_images`` method')
        return skymodelCCDData, skySubtractedCCDData

    def add_data_to_placeholder_images(
            self,
            imageMapOrderDF,
            skymodelCCDData,
            skySubtractedCCDData):
        """*add sky-model and sky-subtracted data to placeholder images*

        **Key Arguments:**
            - ``imageMapOrderDF`` -- single order dataframe from object image and 2D map

        **Usage:**

        ```python
        self.add_data_to_placeholder_images(imageMapOrder)
        ```
        """
        self.log.debug('starting the ``add_data_to_placeholder_images`` method')

        for x, y, skypixel in zip(imageMapOrderDF["x"], imageMapOrderDF["y"], imageMapOrderDF["sky_model"]):
            skymodelCCDData.data[y][x] = skypixel
        for x, y, skypixel in zip(imageMapOrderDF["x"], imageMapOrderDF["y"], imageMapOrderDF["sky_subtracted_flux"]):
            skySubtractedCCDData.data[y][x] = skypixel

        self.log.debug('completed the ``add_data_to_placeholder_images`` method')
        return skymodelCCDData, skySubtractedCCDData

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

        # CREATE THE GID OF AXES
        toprow = fig.add_subplot(gs[0:2, :])
        midrow = fig.add_subplot(gs[2:4, :])
        bottomrow = fig.add_subplot(gs[4:6, :])

        # FIND ORDER PIXELS - MASK THE REST
        nonOrderMask = np.ones_like(objectFrame.data)
        for x, y in zip(self.mapDF["x"], self.mapDF["y"]):
            nonOrderMask[y][x] = 0

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        nonOrderMask = ma.make_mask(nonOrderMask)
        combinedMask = (nonOrderMask == 1) | (objectFrame.mask == 1)

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.rot90(objectFrame.data, 1)
        maskedDataArray = np.ma.array(objectFrame.data, mask=combinedMask)
        std = np.nanstd(maskedDataArray)
        mean = np.nanmean(maskedDataArray)
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
        std = np.nanstd(maskedDataArray)
        mean = np.nanmean(maskedDataArray)
        vmax = 0 + std
        vmin = 0
        bottomrow.imshow(rotatedImg, vmin=vmin, vmax=30, cmap='gray', alpha=1.)
        bottomrow.set_title(
            f"Sky-subtracted {arm} Frame", fontsize=10)
        bottomrow.set_ylabel("x-axis", fontsize=8)
        bottomrow.set_xlabel("y-axis", fontsize=8)
        bottomrow.tick_params(axis='both', which='major', labelsize=9)
        # subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        # fig.suptitle(f"traces of order-centre locations - pinhole flat-frame\n{subtitle}", fontsize=12)

        # plt.show()
        filename = self.filenameTemplate.replace(".fits", "_skysub_quicklook.pdf")

        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home) + "/qc/pdf"
        filePath = f"{outDir}/{filename}"
        plt.savefig(filePath, dpi=720)

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
        print(orderPixelTable.loc[~mask, "wavelength"].min())

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
            log=self.log, CCDObject=imageArray, show=True, ext='data', stdWindow=3, title=False, surfacePlot=True, inst="dummy")

        self.log.debug('completed the ``rectify_order`` method')
        return imageArray

    # use the tab-trigger below for new method
    def fit_surface_to_sky_residuals(
            self,
            skyPixelsDF):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``skyPixelsDF`` -- data frame containing non-clipped pixels containing only sky flux

        **Return:**
            - ``xcoeff`` -- the x-coefficients post clipping
            - ``ycoeff`` -- the y-coefficients post clipping
            - ``res_plots`` -- plot of fit residuals
        """
        self.log.debug('starting the ``fit_surface_to_sky`` method')

        import numpy as np
        from astropy.stats import sigma_clip

        # FIX ME - ADD TO SETTINGS FILE
        orderDeg = 0
        wavelengthDeg = 2
        slitDeg = 3

        # REMOVE FILTERED ROWS FROM DATA FRAME
        skyPixelsDFCopy = skyPixelsDF.copy()
        mask = ((skyPixelsDFCopy['clipped'] == True) | (skyPixelsDFCopy['mask'] == True))
        skyPixelsDFCopy.drop(index=skyPixelsDFCopy[
            mask].index, inplace=True)
        # skyPixelsDFCopy = skyPixelsDFCopy[skyPixelsDFCopy.index % 10 == 0]

        # ADD EXPONENTS TO ORDERTABLE UP-FRONT
        for i in range(0, orderDeg + 1):
            skyPixelsDFCopy[f"order_pow_{i}"] = skyPixelsDFCopy["order"].pow(i)
        for j in range(0, wavelengthDeg + 1):
            skyPixelsDFCopy[f"wavelength_pow_{j}"] = skyPixelsDFCopy["wavelength"].pow(j)
        for k in range(0, slitDeg + 1):
            skyPixelsDFCopy[f"slit_position_pow_{k}"] = skyPixelsDFCopy["slit_position"].pow(k)

        poly = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDeg, wavelengthDeg=wavelengthDeg, slitDeg=slitDeg, exponentsIncluded=True).poly

        clippingSigma = 4
        clippingIterationLimit = 5

        print("\n# FINDING DISPERSION SOLUTION\n")

        iteration = 0
        clippedCount = 1
        from scipy.optimize import curve_fit
        while clippedCount > 0 and iteration < clippingIterationLimit:
            iteration += 1
            sky_subtracted_flux = skyPixelsDFCopy["sky_subtracted_flux"].to_numpy()
            if sky_subtracted_flux.shape[0] < 20:
                return None
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLYS
            # FIRST X
            coeff = np.ones((orderDeg + 1) *
                            (wavelengthDeg + 1) * (slitDeg + 1))
            self.log.info("""curvefit x""" % locals())

            fluxcoeff, pcov_x = curve_fit(
                poly, xdata=skyPixelsDFCopy, ydata=sky_subtracted_flux, p0=coeff)

            self.log.info("""calculate_residuals""" % locals())
            mean_res, std_res, median_res, skyPixelsDFCopy = self.calculate_residuals(
                skyPixelsDF=skyPixelsDFCopy,
                fluxcoeff=fluxcoeff,
                orderDeg=orderDeg,
                wavelengthDeg=wavelengthDeg,
                slitDeg=slitDeg,
                writeQCs=False)

            # SIGMA-CLIP THE DATA
            self.log.info("""sigma_clip""" % locals())
            masked_residuals = sigma_clip(
                skyPixelsDFCopy["residuals_sky_subtracted_flux"], sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
            skyPixelsDFCopy["residuals_masked"] = masked_residuals.mask
            # RETURN BREAKDOWN OF COLUMN VALUE COUNT
            valCounts = skyPixelsDFCopy[
                'residuals_masked'].value_counts(normalize=False)
            if True in valCounts:
                clippedCount = valCounts[True]
            else:
                clippedCount = 0

            if iteration > 1:
                # Cursor up one line and clear line
                sys.stdout.write("\x1b[1A\x1b[2K")

            print(f'ITERATION {iteration:02d}: {clippedCount} pixels clipped in this iteration of fitting the sky subtracted flux residuals')

            # REMOVE FILTERED ROWS FROM DATA FRAME
            mask = (skyPixelsDFCopy['residuals_masked'] == True)
            skyPixelsDFCopy.drop(index=skyPixelsDFCopy[
                mask].index, inplace=True)

            import sqlite3 as sql
            # CONNECT TO THE DATABASE
            conn = sql.connect("pandas_export.db")
            # SEND TO DATABASE
            skyPixelsDFCopy.to_sql('my_export_table', con=conn,
                                   index=False, if_exists='replace')

        self.log.debug('completed the ``fit_surface_to_sky`` method')
        return skyPixelsDFCopy

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
        """*find object slit ranges*

        **Key Arguments:**
            - ``order_dataframes`` -- a list of order data-frames with pixels potentially containing the object flagged.

        **Return:**
            - ``order_dataframes`` -- the order dataframes with the object(s) slit-ranges clipped
            - ``aggressive`` -- mask entire slit range where an object is expected to lie. Default *False*

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

        # BIN FLAGGED PIXEL COUNTS INTO DISCRETE SLIT-POSTION RANGES
        nbins = 100
        minsp = allimageMapOrder['slit_position'].min()
        maxsp = allimageMapOrder['slit_position'].max()
        bins = np.linspace(minsp, maxsp, nbins)
        result = allimageMapOrder['slit_position'].value_counts(bins=bins, sort=False, normalize=True) * nbins

        # REMOVE MEDIAN x 3 -- ONLY OBJECTS SHOULD REMAIN POSITIVE IN COUNTS
        result -= result.median()
        result -= result.abs().median()
        result -= result.abs().median()

        # NEED 3 POSTIVE BINS IN A ROW TO BE SELECTED AS AN OBJECT
        object_ranges = []
        postiveCount = 0
        for wl, count in zip(bins, result):
            if count > 0:
                postiveCount += 1
                upper = wl
                if count > 0.05:
                    record_range = True
            else:
                if postiveCount > 4 and record_range:
                    object_ranges.append([lower, upper])
                postiveCount = 0
                lower = wl
                upper = False
                record_range = False
        if postiveCount > 4:
            object_ranges.append([lower, upper])

        if 1 == 0:
            import matplotlib.pyplot as plt
            print(object_ranges)
            width = (maxsp - minsp) / nbins
            fig, ax = plt.subplots()
            bins = bins[:-1]
            rects1 = ax.bar(bins - width / 2, result, width, label='count')
            fig.tight_layout()
            plt.show()

        # NOW FOR EACH OBJECT SLIT-RANGE, FLAG AS CLIPPED IN ORIGINAL ORDER DATAFRAMES
        for df in order_dataframes:
            for object in object_ranges:
                if aggressive:
                    df.loc[(df['slit_position'].between(object[0], object[1])), "clipped"] = True
                else:
                    df.loc[(df['slit_position'].between(object[0], object[1])) & (df['object'] == True), "clipped"] = True

        self.log.debug('completed the ``clip_object_slit_positions`` method')
        return order_dataframes

    # use the tab-trigger below for new method
    # xt-class-method
