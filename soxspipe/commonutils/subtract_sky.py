#!/usr/bin/env python
# encoding: utf-8
"""
*Subtract the sky background using the Kelson Method*

:Author:
    David Young

:Date Created:
    April 14, 2022
"""
import pandas as pd
import numpy.ma as ma
import numpy as np
from fundamentals import tools
from builtins import object
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils.toolkit import read_spectral_format
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays

import sys
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm
import scipy
from copy import copy
from datetime import datetime
from matplotlib import colors
from soxspipe.commonutils import keyword_lookup
import scipy.interpolate as ip
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.toolkit import quicklook_image
from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
from os.path import expanduser
os.environ['TERM'] = 'vt100'
pd.options.mode.chained_assignment = None


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
    ):
        self.log = log
        log.debug("instansiating a new 'subtract_sky' object")
        self.settings = settings
        self.objectFrame = objectFrame
        self.twoDMap = twoDMap
        self.dispMap = dispMap
        self.qc = qcTable
        self.products = productsTable

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

        # INITIAL ACTIONS
        self.mapDF = twoD_disp_map_image_to_dataframe(log=self.log, slit_length=dp["slit_length"], twoDMapPath=twoDMap, assosiatedFrame=self.objectFrame)
        quicklook_image(
            log=self.log, CCDObject=objectFrame, show=False, ext=False, stdWindow=5, title="Object Frame with dispersion solution grid", surfacePlot=True, dispMap=dispMap, dispMapDF=self.mapDF)

        if self.inst == "SOXS":
            self.axisA = "y"
            self.axisB = "x"
        elif self.inst == "XSHOOTER":
            self.axisA = "x"
            self.axisB = "y"

        self.dateObs = objectFrame.header[kw("DATE_OBS")]

        # GET A TEMPLATE FILENAME USED TO NAME PRODUCTS
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

        print(f'\n# MODELLING SKY BACKGROUND AND REMOVING FROM SCIENCE FRAME')

        skymodelCCDData, skySubtractedCCDData = self.create_placeholder_images()

        uniqueOrders = self.mapDF['order'].unique()
        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # BSPLINE ORDER TO FIT SKY WITH
        bspline_order = self.settings["sky-subtraction"]["bspline_order"]

        qcPlotOrder = int(np.median(uniqueOrders)) - 1
        # REMOVE ME

        for o in uniqueOrders:
            imageMapOrder = self.mapDF[self.mapDF["order"] == o]
            imageMapOrderWithObject, imageMapOrder = self.get_over_sampled_sky_from_order(imageMapOrder, o, ignoreBP=True)
            # REMOVE ME
            imageMapOrder = self.fit_bspline_to_sky_dead(imageMapOrder, o, bspline_order)
            # self.rectify_order(order=o, imageMapOrder=imageMapOrder)

            imageMapOrder = self.fit_bspline_to_sky(imageMapOrder, o, bspline_order)
            # imageMapOrder = self.fit_surface_to_sky(imageMapOrder, o, bspline_order)
            skymodelCCDData, skySubtractedCCDData = self.add_data_to_placeholder_images(imageMapOrder, skymodelCCDData, skySubtractedCCDData)
            if o == qcPlotOrder:
                qc_plot_path = self.plot_sky_sampling(order=o, imageMapOrderWithObjectDF=imageMapOrderWithObject, imageMapOrderDF=imageMapOrder)
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

        filename = self.filenameTemplate.replace("_SLIT", "_SKYMODEL")
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
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

        filename = self.filenameTemplate.replace("_SLIT", "_SKYSUB")
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
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
            ignoreBP=True):
        """*unpack the over sampled sky from an order*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe from object image and 2D map
            - ``order`` -- the order number
            - ``ignoreBP`` -- ignore bad-pixels? Deafult *True*

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

        imageMapOrder.sort_values("wavelength", inplace=True)
        imageMapOrder["clipped"] = False

        if not ignoreBP:
            # REMOVE FILTERED ROWS FROM DATA FRAME
            mask = (imageMapOrder['mask'] == True)
            imageMapOrder.drop(index=imageMapOrder[mask].index, inplace=True)

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
            imageMapOrderDF):
        """*generate a plot of sky sampling*

        **Key Arguments:**
            - ``order`` -- the order number.
            - ``imageMapOrderWithObjectDF`` -- dataframe with various processed data without object clipped
            - ``imageMapOrderDF`` -- dataframe with various processed data for order

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

        # SET COLOURS FOR VARIOUS STAGES
        medianColor = "blue"
        percentileColor = "blue"
        skyColor = "purple"
        rawColor = "#93a1a1"
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
        # fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=320)
        # REMOVE ME
        fig = plt.figure(figsize=(8, 9), constrained_layout=True, dpi=100)
        gs = fig.add_gridspec(10, 4)
        # CREATE THE GID OF AXES
        onerow = fig.add_subplot(gs[1:2, :])
        tworow = fig.add_subplot(gs[2:4, :])
        threerow = fig.add_subplot(gs[4:5:, :])
        fourrow = fig.add_subplot(gs[5:6:, :])
        fiverow = fig.add_subplot(gs[6:7:, :])
        sixrow = fig.add_subplot(gs[7:8:, :])
        sevenrow = fig.add_subplot(gs[8:9:, :])
        eightrow = fig.add_subplot(gs[9:10:, :])

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
        # onerow.text(5, 5, 'your legend', bbox={'facecolor': 'white', 'pad': 10})
        onerow.invert_xaxis()

        # ORIGINAL DATA AND PERCENTILE SMOOTHED WAVELENGTH VS FLUX
        raw = tworow.plot(
            imageMapOrderWithObjectDF["wavelength"].values,
            imageMapOrderWithObjectDF["flux"].values, label='unprocessed (unp)', c=rawColor, zorder=0)
        # RAW MARKERS
        tworow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].values, label='unclipped', s=rawMS, c=rawColor, alpha=0.5, zorder=unclippedZ)
        # ROBUSTLY CLIPPED
        tworow.scatter(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "flux"].values, label='clipped', s=percentileMS, marker="x", c=percentileColor, zorder=percentileZ, alpha=0.2)
        # MEDIAN CLIPPED
        # label='median clipped'
        label = None
        tworow.scatter(
            imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "wavelength"].values,
            imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "flux"].values, label=label, s=medianMS, marker="x", c=medianColor, zorder=medianZ, alpha=0.2)
        # MEDIAN LINE
        # label='median-smoothed (ms)'
        label = None
        # tworow.plot(
        #     imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "wavelength"].values,
        #     imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "flux_smoothed"].values, label=label, c=medianColor, zorder=medianZ)
        # PERCENTILE LINE
        tworow.plot(
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "wavelength"].values,
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_smoothed"].values, label='percentile-smoothed', c=percentileColor, zorder=percentileZ)

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
        tworow.set_xticks([], [])

        # WAVELENGTH RESIDUAL PANEL
        threerow.scatter(imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "wavelength"].values, imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "residual_global_sigma"].values, s=rawMS, alpha=0.5, c=rawColor, zorder=unclippedZ)
        threerow.scatter(imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "wavelength"].values, imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == True, "residual_global_sigma"].values, s=medianMS, marker="x", c=medianColor, zorder=medianZ, alpha=0.2)
        threerow.scatter(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "wavelength"].values, imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True, "residual_global_sigma_old"].values, s=percentileMS, marker="x", c=percentileColor, zorder=percentileZ, alpha=0.2)
        std = imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "residual_global_sigma"].std()
        median = imageMapOrderWithObjectDF.loc[imageMapOrderWithObjectDF["clipped"] == False, "residual_global_sigma"].median()
        threerow.set_ylim(median - 3 * std, median + 7 * std)
        threerow.set_xlabel(
            "wavelength (nm)", fontsize=10)
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
        fiverow.invert_xaxis()

        # PLOT WAVELENGTH VS FLUX SKY MODEL
        sixrow.scatter(
            imageMapOrderDF["wavelength"].values,
            imageMapOrderDF["flux"].values, label='unclipped', s=rawMS, c=rawColor, alpha=0.5, zorder=unclippedZ)
        skymodel = sixrow.plot(
            imageMapOrderDF["wavelength"].values,
            imageMapOrderDF["sky_model"].values, label='sky model', c=skyColor, zorder=skyZ)
        if ylimmin < -3000:
            ylimmin = -300
        sixrow.set_ylim(ylimmin, imageMapOrderWithObjectDF["flux_smoothed"].max() * 1.2)
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
        eightrow.invert_xaxis()
        medianValue = np.median(skySubImage.data.ravel())
        color = im.cmap(im.norm(medianValue))
        patches = [mpatches.Patch(color=color, label="sky-subtracted frame")]
        eightrow.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        fig.suptitle(f"{self.arm} sky model: order {order}", fontsize=12, y=0.97)

        filename = self.filenameTemplate.split("SLIT")[0] + f"SKYMODEL_QC_PLOTS_ORDER_{int(order)}.pdf"
        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        filePath = f"{outDir}/{filename}"
        # REMOVE ME
        plt.show()
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
                imageMapOrderDF["flux_smoothed"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"].rolling(window=windowSize, center=True).quantile(.35)

            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux"] - imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_smoothed"]
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"].rolling(windowSize).std()
            imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_sigma"] = imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "flux_minus_smoothed_residual"] / imageMapOrderDF.loc[imageMapOrderDF["clipped"] == False, "residual_windowed_std"]

            if median_centre_func:
                # CLIP ABOVE AND BELOW
                imageMapOrderDF.loc[imageMapOrderDF["residual_windowed_sigma"].abs() > sigma_clip_limit, "clipped"] = True
            else:
                # CLIP ONLY HIGH VALUES
                imageMapOrderDF.loc[((imageMapOrderDF["residual_windowed_sigma"] > sigma_clip_limit) & (imageMapOrderDF["residual_global_sigma_old"] > -10.0)), "clipped"] = True

            if newlyClipped == -1:
                totalClipped = len(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True].index)
                newlyClipped = totalClipped
            else:
                newlyClipped = len(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True].index) - totalClipped
                totalClipped = len(imageMapOrderDF.loc[imageMapOrderDF["clipped"] == True].index)
            if i > 1:
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

    def fit_bspline_to_sky_dead(
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
        imageMapOrder = self.fit_bspline_to_sky(
            imageMapOrder,
            myOrder
        )
        ```

        """
        self.log.debug('starting the ``fit_bspline_to_sky`` method')

        # SORT BY COLUMN NAME
        df = imageMapOrder.copy()
        df.sort_values(by=['wavelength'], inplace=True)
        df.drop_duplicates(subset=['wavelength'], inplace=True)

        goodWl = df.loc[df["clipped"] == False]["wavelength"]
        goodFlux = df.loc[df["clipped"] == False]["flux"]
        goodSlit = df.loc[df["clipped"] == False]["slit_position"]
        df["weights"] = 1 / df["flux_scatter_windowed_std"].abs()
        df["weights"] = df["weights"].replace(np.nan, 0.0000000001)
        goodWeights = df.loc[df["clipped"] == False, "weights"]

        goodWl = goodWl.values
        goodFlux = goodFlux.values
        goodSlit = goodSlit.values
        goodWeights = goodWeights.values

        import scipy
        interp_func = scipy.interpolate.SmoothBivariateSpline(goodWl, goodSlit, goodFlux, goodWeights)
        seedx = np.arange(int(goodWl.min()), int(goodWl.max()), 0.1)
        seedy = np.arange(int(goodSlit.min()), int(goodSlit.max()), 0.05)
        grid_x, grid_y = np.meshgrid(seedx, seedy)
        grid_x = grid_x.flatten()
        grid_y = grid_y.flatten()

        interp_result = interp_func(seedx, seedy)
        interp_result = np.rot90(np.flipud(interp_result), k=-1)

        extent = [seedx.min(), seedx.max(), seedy.min(), seedy.max()]
        print(extent)

        fig, ax = plt.subplots()
        std = np.nanstd(goodFlux)
        mean = np.nanmean(goodFlux)
        vmax = mean + 3 * std
        vmin = mean - 3 * std
        shw = ax.imshow(interp_result, aspect='auto', cmap='nipy_spectral',
                        extent=extent, vmin=vmin, vmax=vmax)
        # plt.plot(grid_x, grid_y, 'w.')
        fig.colorbar(shw)
        plt.show()

        extent = [test_x[0], test_x[-1], test_y[0], test_y[-1]]
        im = axes[0].imshow(perfect_result, aspect='auto', cmap='nipy_spectral', extent=extent, vmin=-1.5, vmax=2.5)
        fig.colorbar(im, ax=axes[0])
        axes[0].plot(train_x, train_y, 'w.')
        axes[0].set_title('Perfect result, sampled function')
        im = axes[1].imshow(interp_result, aspect='auto', cmap='nipy_spectral', extent=extent, vmin=-1.5, vmax=2.5)
        fig.colorbar(im, ax=axes[1])
        axes[1].plot(train_x, train_y, 'w.')
        axes[1].set_title('SmoothBivariateSpline')
        plt.show()

        # some preparation, loading data and stuff
        # all my data is stored in 'data'

        # Create the knots (10 knots in each direction, making 100 total
        wlKnot = numpy.linspace(5, data.shape[0] - 5, 10)
        slKnot = numpy.linspace(5, data.shape[1] - 5, 10)

        # Create all weights, and set them to 0 when the data is NaN
        weights = numpy.ones(data.shape)
        weights[numpy.isnan(data)] = 1e-15  # weights must be >0

        # LSQBivariateSpline needs x and y coordinates as 1-D arrays
        x, y = numpy.indices(data.shape)
        spline_fit = scipy.interpolate.LSQBivariateSpline(x.ravel(), y.ravel(), data.ravel(),
                                                          xcoord, ycoord,
                                                          w=weights.ravel(),
                                                          bbox=[None, None, None, None],
                                                          kx=2, ky=2)
        # N = int(goodWl.shape[0] / 25)
        # seedKnots = np.linspace(xmin, xmax, N)

        # t for knots
        # c of coefficients
        # k for order

        # Fit
        n_interior_knots = int(imageMapOrder["wavelength"].values.shape[0] / 9)
        qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
        knots = np.quantile(goodWl, qs)
        # tck = ip.splrep(goodWl, goodFlux, t=knots, k=3)
        tck = ip.splrep(goodWl, goodFlux, t=knots, k=bspline_order)
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

        self.log.debug('completed the ``fit_bspline_to_sky`` method')
        return imageMapOrder

    def fit_bspline_to_sky(
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
        imageMapOrder = self.fit_bspline_to_sky(
            imageMapOrder,
            myOrder
        )
        ```

        """
        self.log.debug('starting the ``fit_bspline_to_sky`` method')

        # SORT BY COLUMN NAME
        df = imageMapOrder.copy()
        df.sort_values(by=['wavelength'], inplace=True)
        df.drop_duplicates(subset=['wavelength'], inplace=True)

        goodWl = df.loc[df["clipped"] == False]["wavelength"]
        goodFlux = df.loc[df["clipped"] == False]["flux"]
        df["weights"] = 1 / df["flux_scatter_windowed_std"].abs()
        df["weights"] = df["weights"].replace(np.nan, 0)
        goodWeights = df.loc[df["clipped"] == False, "weights"]

        # N = int(goodWl.shape[0] / 25)
        # seedKnots = np.linspace(xmin, xmax, N)

        # t for knots
        # c of coefficients
        # k for order

        # Fit
        n_interior_knots = int(imageMapOrder["wavelength"].values.shape[0] / 9)
        qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
        knots = np.quantile(goodWl, qs)
        # tck = ip.splrep(goodWl, goodFlux, t=knots, k=3)
        tck = ip.splrep(goodWl, goodFlux, t=knots, k=bspline_order)
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

        self.log.debug('completed the ``fit_bspline_to_sky`` method')
        return imageMapOrder

    def fit_surface_to_sky(
            self,
            imageMapOrder,
            order,
            bspline_order):
        """*fit a nurbs surface to the unclipped sky pixels (wavelength vs flux)*

        **Key Arguments:**
            - ``imageMapOrder`` -- single order dataframe, containing sky flux with object(s) and CRHs removed
            - ``order`` -- the order number

        **Return:**
            - ``imageMapOrder`` -- same `imageMapOrder` as input but now with `sky_model` (bspline fit of the sky) and `sky_subtracted_flux` columns

        **Usage:**

        ```python
        imageMapOrder = self.fit_surface_to_sky(
            imageMapOrder,
            myOrder
        )
        ```

        """
        self.log.debug('starting the ``fit_surface_to_sky`` method')

        dispMap = self.dispMap
        kw = self.kw
        dp = self.detectorParams
        arm = self.arm

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm)

        grid_res_slit = 0.5
        grid_res_wl = 1
        slitLength = dp["slit_length"]
        halfGrid = (slitLength / 2) * 1.1
        slitArray = np.arange(-halfGrid, halfGrid +
                              grid_res_slit, grid_res_slit)

        for o, minWl, maxWl in zip(orderNums, waveLengthMin, waveLengthMax):
            if o == order:
                orderInfo = (order, minWl, maxWl)
        (order, minWl, maxWl) = orderInfo
        wlArray = np.arange(minWl - 20, maxWl + 20, grid_res_wl)

        # SIZES NEEDED LATER FOR NURBS
        wl_size = wlArray.shape[0]
        sp_size = slitArray.shape[0]

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

        from tabulate import tabulate

        # print(tabulate(orderPixelTable.head(100), headers='keys', tablefmt='psql'))
        # print(tabulate(imageMapOrder.head(100), headers='keys', tablefmt='psql'))

        orderPixelTable = pd.merge(orderPixelTable, imageMapOrder[['x', 'y', 'flux', 'clipped']], how='left', left_on=[
            'pixel_x', 'pixel_y'], right_on=['x', 'y'])
        # DROP MISSING VALUES
        # orderPixelTable.dropna(axis='index', how='any', subset=['x'], inplace=True)

        orderPixelTable = orderPixelTable[['order', 'wavelength', 'slit_position', 'fit_x', 'fit_y', 'flux', 'clipped']]
        orderPixelTable['weight'] = 100
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable['flux'].isnull())
        # UPDATE FILTERED VALUES
        orderPixelTable.loc[mask, 'weight'] = 1
        orderPixelTable.loc[mask, 'flux'] = 1
        mask = (orderPixelTable['clipped'] == True)
        orderPixelTable.loc[mask, 'weight'] = 1
        mask = (orderPixelTable['weight'] > 1)

        orderPixelTable['flux'] /= 1000

        orderPixelTable['ctrlpts'] = orderPixelTable[['fit_x', 'fit_y', 'flux']].to_numpy().tolist()
        # orderPixelTable['ctrlpts'] = orderPixelTable.apply(lambda row: [row['fit_x'], row['fit_y'], row['flux'], row['weight']], axis=1)
        ctrlpts = orderPixelTable['ctrlpts'].values
        weights = orderPixelTable['weight'].values

        from geomdl import NURBS
        from geomdl import utilities as utils
        from geomdl import compatibility as compat
        from geomdl.visualization import VisMPL

        degree_sp = 2
        degree_wl = 3

        # Combine weights vector with the control points list
        t_ctrlptsw = compat.combine_ctrlpts_weights(ctrlpts, weights)
        # Since NURBS-Python uses v-row order, we need to convert the exported ones
        n_ctrlptsw = compat.flip_ctrlpts_u(t_ctrlptsw, sp_size, wl_size)
        # Since we have no information on knot vectors, let's auto-generate them
        n_knotvector_sp = utils.generate_knot_vector(degree_sp, sp_size)
        n_knotvector_wl = utils.generate_knot_vector(degree_wl, wl_size)

        # Create a NURBS surface instance
        surf = NURBS.Surface()

        # Fill the surface object
        surf.degree_u = degree_sp
        surf.degree_v = degree_wl
        surf.set_ctrlpts(n_ctrlptsw, sp_size, wl_size)
        surf.knotvector_u = n_knotvector_sp
        surf.knotvector_v = n_knotvector_wl

        # Set evaluation delta
        surf.delta = 0.05

        # Set visualization component
        vis_config = VisMPL.VisConfig(alpha=0.75, ctrlpts=False)
        vis_comp = VisMPL.VisSurface(vis_config)
        surf.vis = vis_comp

        # Render the surface
        surf.render()

        sys.exit(0)

        # SORT BY COLUMN NAME
        df = imageMapOrder.copy()
        df.sort_values(by=['wavelength'], inplace=True)
        df.drop_duplicates(subset=['wavelength'], inplace=True)

        # data = df.loc[df["clipped"] == False][["x", "y", "wavelength", "slit_position", "flux"]]
        xlen = data["x"].max() - data["x"].min()
        ylen = data["y"].max() - data["y"].min()
        xmin = data["x"].min()
        ymin = data["y"].min()

        ctrlpts = []
        weights = []
        badWeight = 0.00000001
        for x in range(xlen + 1):
            row = []
            for y in range(ylen + 1):
                row.append([x * badWeight, y * badWeight, 0, badWeight])
            ctrlpts.append(row)

        # sys.exit(0)

        for x, y, wl, sp, fx, clipped in zip(df["x"], df["y"], df["wavelength"], df["slit_position"], df["flux"], df["clipped"]):
            if fx > 0 and fx < 1000:
                ctrlpts[x - xmin][y - ymin][2] = fx
                ctrlpts[x - xmin][y - ymin][3] = 1

        # Create a BSpline surface
        surf = NURBS.Surface()

        # Set degrees
        surf.degree_u = 3
        surf.degree_v = 3

        # Set control points
        surf.ctrlpts2d = ctrlpts
        # surf_set_ctrlpts(n_ctrlptsw, p_size_u, p_size_v)

        # Set knot vectors
        from geomdl import utilities as utils
        n_knotvector_u = utils.generate_knot_vector(3, xlen + 1)
        n_knotvector_v = utils.generate_knot_vector(3, ylen + 1)
        surf.knotvector_u = n_knotvector_u
        surf.knotvector_v = n_knotvector_v

        # Set evaluation delta
        surf.delta = 0.025

        # Evaluate surface points
        surf.evaluate()

        # Import and use Matplotlib's colormaps
        from matplotlib import cm

        # Plot the control points grid and the evaluated surface
        vis_config = VisMPL.VisConfig(alpha=0.75, ctrlpts=True)
        vis_obj = VisMPL.VisSurface(vis_config)
        surf.vis = vis_obj
        surf.render(colormap=cm.cool)

        sys.exit(0)

        goodWl = df.loc[df["clipped"] == False]["wavelength"]
        goodFlux = df.loc[df["clipped"] == False]["flux"]
        df["weights"] = 1 / df["flux_scatter_windowed_std"].abs()
        df["weights"] = df["weights"].replace(np.nan, 0)
        goodWeights = df.loc[df["clipped"] == False, "weights"]

        # N = int(goodWl.shape[0] / 25)
        # seedKnots = np.linspace(xmin, xmax, N)

        # t for knots
        # c of coefficients
        # k for order

        # Fit
        n_interior_knots = int(imageMapOrder["wavelength"].values.shape[0] / 9)
        qs = np.linspace(0, 1, n_interior_knots + 2)[1:-1]
        knots = np.quantile(goodWl, qs)
        # tck = ip.splrep(goodWl, goodFlux, t=knots, k=3)
        tck = ip.splrep(goodWl, goodFlux, t=knots, k=bspline_order)
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

        self.log.debug('completed the ``fit_surface_to_sky`` method')
        return imageMapOrder

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
        filename = self.filenameTemplate.split("SLIT")[0] + "skysub_quicklook.pdf"

        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        filePath = f"{outDir}/{filename}"
        plt.savefig(filePath, dpi=720)

        self.log.debug('completed the ``plot_results`` method')
        return filePath

    # use the tab-trigger below for new method
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

        minWl = minWl + 20
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

        # SET THE WEIGHTS
        orderPixelTable["weight"] = 100
        mask = (orderPixelTable['clipped'] == True)
        orderPixelTable.loc[mask, "weight"] = 1
        mask = (orderPixelTable['flux'].isnull())
        orderPixelTable.loc[mask, "weight"] = 1
        mask = (orderPixelTable['clipped'] == True)
        orderPixelTable.loc[mask, "weight"] = 1
        mask = (orderPixelTable['weight'] == 1)
        orderPixelTable.loc[mask, "flux"] = 0

        orderPixelTable["flux"] /= 1000

        wl_array = orderPixelTable["wavelength"].values
        sl_array = orderPixelTable["slit_position"].values
        fx_array = orderPixelTable["flux"].values
        weight_array = orderPixelTable["weight"].values

        from geomdl import NURBS
        from geomdl import utilities as utils
        from geomdl import compatibility as compat
        from geomdl.visualization import VisMPL

        degree_sp = 2
        degree_wl = 3
        sp_size = len(slitArray)
        wl_size = len(wlArray)

        orderPixelTable['ctrlpts'] = orderPixelTable[['fit_x', 'fit_y', 'flux']].to_numpy().tolist()
        ctrlpts = orderPixelTable['ctrlpts'].values
        # Combine weights vector with the control points list
        t_ctrlptsw = compat.combine_ctrlpts_weights(ctrlpts, weight_array)
        # Since NURBS-Python uses v-row order, we need to convert the exported ones
        n_ctrlptsw = compat.flip_ctrlpts_u(t_ctrlptsw, sp_size, wl_size)
        # Since we have no information on knot vectors, let's auto-generate them
        n_knotvector_sp = utils.generate_knot_vector(degree_sp, sp_size)
        n_knotvector_wl = utils.generate_knot_vector(degree_wl, wl_size)

        # Create a NURBS surface instance
        surf = NURBS.Surface()

        # Fill the surface object
        surf.degree_u = degree_sp
        surf.degree_v = degree_wl
        surf.set_ctrlpts(n_ctrlptsw, sp_size, wl_size)
        surf.knotvector_u = n_knotvector_sp
        surf.knotvector_v = n_knotvector_wl

        # Set evaluation delta
        surf.delta = 0.05

        # Set visualization component
        vis_config = VisMPL.VisConfig(alpha=0.75, ctrlpts=False)
        vis_comp = VisMPL.VisSurface(vis_config)
        surf.vis = vis_comp

        # Render the surface
        surf.render()

        # print("FITTING")
        # print(len(wl_array))
        # interp_func = scipy.interpolate.SmoothBivariateSpline(wl_array, sl_array, fx_array, weight_array, s=0.0)
        # seedx = np.arange(int(wl_array.min()), int(wl_array.max()), 0.1)
        # seedy = np.arange(int(sl_array.min()), int(sl_array.max()), 0.05)
        # grid_x, grid_y = np.meshgrid(seedx, seedy)
        # grid_x = grid_x.flatten()
        # grid_y = grid_y.flatten()

        # interp_result = interp_func(seedx, seedy)
        # interp_result = np.rot90(np.flipud(interp_result), k=-1)

        # extent = [seedx.min(), seedx.max(), seedy.min(), seedy.max()]
        # print(extent)

        xcoord = np.linspace(wl_array.min(), wl_array.max(), 10)
        ycoord = np.linspace(sl_array.min(), sl_array.max(), 10)
        interp_func = scipy.interpolate.LSQBivariateSpline(wl_array, sl_array, fx_array,
                                                           wl_array, sl_array,
                                                           w=weight_array,
                                                           kx=2, ky=2)

        seedx = np.arange(int(wl_array.min()), int(wl_array.max()), 0.1)
        seedy = np.arange(int(sl_array.min()), int(sl_array.max()), 0.05)
        grid_x, grid_y = np.meshgrid(seedx, seedy)
        grid_x = grid_x.flatten()
        grid_y = grid_y.flatten()

        interp_result = interp_func(seedx, seedy)
        interp_result = np.rot90(np.flipud(interp_result), k=-1)

        quicklook_image(
            log=self.log, CCDObject=interp_result, show=True, ext='data', stdWindow=3, title=False, surfacePlot=True, inst="dummy")

        self.log.debug('completed the ``rectify_order`` method')
        return imageArray

    # use the tab-trigger below for new method
    # xt-class-method
