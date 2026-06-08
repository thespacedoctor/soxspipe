#!/usr/bin/env python
# encoding: utf-8
"""
*perform optimal source extraction using the Horne method (Horne 1986)*

Author
: Marco Landoni & David Young

Date Created
: May 17, 2023
"""

from fundamentals import tools
from builtins import object
import sys
import os
from line_profiler import profile

os.environ["TERM"] = "vt100"


# TODO: include the BPM in create_cross_dispersion_slice (at least)
# TODO: find a more robust solution for when horneDenominatorSum == 0 (all pixels in a slice do not pass variance cuts). See where fudged == true
# TODO: revisit how the wavelength for each slice is calculated ... take from the continuum, or the central 3-5 pixels?


class horne_extraction(object):
    """
    *perform optimal source extraction using the Horne method (Horne 1986)*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
    - ``skySubtractedFrame`` -- path to sky subtracted frame
    - ``unflattenedFrame`` -- path to unflattened frame
    - ``subtractedFrame`` -- path to the frame that has been subtracted in the AB or BA sequence. This is use to extract a sky spectrum. Default *False* (i.e. not used)
    - ``twoDMapPath`` -- path to 2D dispersion map image path
    - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products (if False no products are saved to file)
    - ``dispersionMap`` -- the FITS binary table containing dispersion map polynomial
    - ``sofName`` -- the set-of-files filename
    - ``locationSetIndex`` -- the index of the AB cycle locations (nodding mode only). Default *False*
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""
    - ``notFlattened`` -- flag to indicate if the frame is flattened or not. Default *False*
    - ``debug`` -- flag to indicate if debug mode is on (shows plots). Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/main/initialisation.html).

    To initiate a horne_extraction object, use the following:

    ```python
    from soxspipe.commonutils import horne_extraction
    optimalExtractor = horne_extraction(
        log=log,
        skySubtractedFrame=skySubtractedFrame,
        unflattenedFrame=unflattenedFrame,
        subtractedFrame=subtractedFrame,
        twoDMapPath=twoDMap,
        settings=settings,
        recipeName="soxs-stare",
        qcTable=qc,
        productsTable=products,
        dispersionMap=dispMap,
        sofName=sofName,
        locationSetIndex=locationSetIndex,
        startNightDate=startNightDate,
        debug=debug,
        turnOffMP=turnOffMP

    )
    qc, products = optimalExtractor.extract()
    ```

    """

    def __init__(
        self,
        log,
        settings,
        recipeSettings,
        skySubtractedFrame,
        unflattenedFrame,
        twoDMapPath,
        subtractedFrame=False,
        recipeName=False,
        qcTable=False,
        productsTable=False,
        dispersionMap=False,
        sofName=False,
        locationSetIndex=False,
        startNightDate="",
        notFlattened=False,
        debug=False,
        turnOffMP=False,
    ):
        import numpy as np
        import pandas as pd
        from astropy.io import fits
        from astropy.nddata import CCDData
        from astropy import units as u
        from os.path import expanduser
        from soxspipe.commonutils import keyword_lookup
        from soxspipe.commonutils import detect_continuum
        from soxspipe.commonutils.toolkit import unpack_order_table
        from soxspipe.commonutils import detector_lookup
        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
        import matplotlib.pyplot as plt
        from soxspipe.commonutils.toolkit import get_skylines_dataframe

        self.log = log
        log.debug("instantiating a new 'horne_extraction' object")
        self.dispersionMap = dispersionMap
        self.twoDMapPath = twoDMapPath
        self.settings = settings
        self.products = productsTable
        self.qc = qcTable
        self.recipeName = recipeName
        self.sofName = sofName
        self.noddingSequence = ""
        self.recipeSettings = recipeSettings
        self.startNightDate = startNightDate
        self.debug = debug
        self.unflattenedFrame = unflattenedFrame
        self.turnOffMP = turnOffMP
        self.subtractedFrame = subtractedFrame

        if notFlattened:
            self.notFlattened = "_NOTFLAT"
        else:
            self.notFlattened = ""

        if notFlattened:
            self.skySubtractedFrame = unflattenedFrame

        # DETECTING SEQUENCE AUTOMATICALLY
        try:
            self.noddingSequence = "_A" if int(skySubtractedFrame.header["HIERARCH ESO SEQ CUMOFF Y"] > 0) else "_B"
            if locationSetIndex:
                self.noddingSequence += str(locationSetIndex)
        except:
            self.noddingSequence = ""

        # COLLECT SETTINGS FROM SETTINGS FILE
        self.slitHalfLength = int(self.recipeSettings["horne-extraction-slit-length"] / 2)
        self.clippingSigma = self.recipeSettings["horne-extraction-profile-clipping-sigma"]
        self.clippingIterationLimit = self.recipeSettings["horne-extraction-profile-clipping-iteration-count"]
        self.globalClippingSigma = self.recipeSettings["horne-extraction-profile-global-clipping-sigma"]

        # TODO: replace this value with true value from FITS header object
        self.ron = 3.0

        # OPEN THE SKY-SUBTRACTED FRAME
        if isinstance(skySubtractedFrame, CCDData):
            self.skySubtractedFrame = skySubtractedFrame
        else:
            self.skySubtractedFrame = CCDData.read(
                skySubtractedFrame,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(log=self.log, settings=self.settings).get
        kw = self.kw
        self.arm = self.skySubtractedFrame.header[kw("SEQ_ARM")]
        self.dateObs = self.skySubtractedFrame.header[kw("DATE_OBS")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

        # GET SKYLINES DATAFRAME
        self.skylinesDF = get_skylines_dataframe(
            self.log, self.settings, self.arm, minBrightnessVIS=1, minBrightnessNIR=0.5
        )

        if self.twoDMapPath:
            self.mapDF, self.interOrderMask = twoD_disp_map_image_to_dataframe(
                log=self.log,
                slit_length=self.detectorParams["slit_length"],
                twoDMapPath=twoDMapPath,
                associatedFrame=self.skySubtractedFrame,
                kw=kw,
                dispAxis=self.detectorParams["dispersion-axis"],
            )

            self.skySubtractedFrame.data[self.interOrderMask == 1] = np.nan

        # USE SUBTRACTED FRAME (NODDING ONLY) TO EXTRACT A SKY SPECTRUM
        if subtractedFrame == False:
            self.subtractedFrame = None
        else:
            # OPEN THE SKY-MODEL FRAME
            if isinstance(subtractedFrame, CCDData):
                self.subtractedFrame = subtractedFrame
            else:
                self.subtractedFrame = CCDData.read(
                    subtractedFrame,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

        # GET A TEMPLATE FILENAME USED TO NAME PRODUCTS
        if self.sofName:
            self.filenameTemplate = self.sofName + ".fits"
        else:
            self.filenameTemplate = filenamer(log=self.log, frame=self.skySubtractedFrame, settings=self.settings)

        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log,
            settings=settings,
            recipeName=recipeName,
            startNightDate=startNightDate,
        )

        # OPEN AND UNPACK THE 2D IMAGE MAP
        self.twoDMap = fits.open(twoDMapPath)

        try:
            dpBinx = self.twoDMap[0].header[kw("WIN_BINX")]
            dpBiny = self.twoDMap[0].header[kw("WIN_BINY")]
        except:
            dpBinx = 1
            dpBiny = 1

        # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        binx = 1
        biny = 1
        try:
            binx = int(self.skySubtractedFrame.header[kw("WIN_BINX")])
            biny = int(self.skySubtractedFrame.header[kw("WIN_BINY")])
        except:
            pass

        # ADJUST SLIT HEIGHT IF BINNING IN USE
        if binx > 1 or biny > 1:
            if self.detectorParams["dispersion-axis"] == "x":
                self.slitHalfLength /= binx
            else:
                self.slitHalfLength /= biny
            self.slitHalfLength = round(self.slitHalfLength)

        binxRatio = binx / dpBinx
        binyRatio = biny / dpBiny

        xdim = int(self.twoDMap[0].data.shape[1] / binxRatio)
        ydim = int(self.twoDMap[0].data.shape[0] / binyRatio)
        xarray = np.tile(np.arange(0, xdim), ydim)
        yarray = np.repeat(np.arange(0, ydim), xdim)

        # self.skySubtractedFrame.data[self.skySubtractedFrame.data == 0] = np.nan

        if binxRatio > 1 or binyRatio > 1:
            from astropy.nddata import block_reduce

            self.twoDMap["WAVELENGTH"].data = block_reduce(
                self.twoDMap["WAVELENGTH"].data, (binyRatio, binxRatio), func=np.mean
            )
            self.twoDMap["SLIT"].data = block_reduce(self.twoDMap["SLIT"].data, (binyRatio, binxRatio), func=np.mean)
            self.twoDMap["ORDER"].data = block_reduce(self.twoDMap["ORDER"].data, (binyRatio, binxRatio), func=np.mean)

        self.imageMap = pd.DataFrame.from_dict(
            {
                "x": xarray,
                "y": yarray,
                "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().astype(np.float32),
                "slit_position": self.twoDMap["SLIT"].data.flatten().astype(np.float32),
                "order": self.twoDMap["ORDER"].data.flatten().astype(np.float32),
                "flux": self.skySubtractedFrame.data.flatten().astype(np.float32),
            }
        )
        self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)

        # REMOVE IF THE ABOVE .astype(float) CONVERSION IS WORKING
        # try:
        #     self.imageMap = pd.DataFrame.from_dict({
        #         "x": xarray,
        #         "y": yarray,
        #         "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().astype(float),
        #         "slit_position": self.twoDMap["SLIT"].data.flatten().astype(float),
        #         "order": self.twoDMap["ORDER"].data.flatten().astype(int),
        #         "flux": self.skySubtractedFrame.data.flatten().astype(float)
        #     })
        #     self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)
        # except:
        #     try:
        #         self.imageMap = pd.DataFrame.from_dict({
        #             "x": xarray,
        #             "y": yarray,
        #             "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().byteswap().newbyteorder(),
        #             "slit_position": self.twoDMap["SLIT"].data.flatten().byteswap().newbyteorder(),
        #             "order": self.twoDMap["ORDER"].data.flatten().byteswap().newbyteorder(),
        #             "flux": self.skySubtractedFrame.data.flatten().byteswap().newbyteorder()
        #         })
        #         self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)
        #     except:
        #         self.imageMap = pd.DataFrame.from_dict({
        #             "x": xarray,
        #             "y": yarray,
        #             "wavelength": self.twoDMap["WAVELENGTH"].data.flatten(),
        #             "slit_position": self.twoDMap["SLIT"].data.flatten(),
        #             "order": self.twoDMap["ORDER"].data.flatten(),
        #             "flux": self.skySubtractedFrame.data.flatten()
        #         })
        #         self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)

        # REMOVE ZEROS
        mask = (self.imageMap["wavelength"] == 0) & (self.imageMap["slit_position"] == 0)
        self.imageMap = self.imageMap.loc[~mask]

        # FIND THE OBJECT TRACE IN EACH ORDER
        detector = detect_continuum(
            log=self.log,
            traceFrame=self.unflattenedFrame,
            dispersion_map=self.dispersionMap,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            sofName=self.sofName,
            recipeName=self.recipeName,
            qcTable=self.qc,
            productsTable=self.products,
            locationSetIndex=locationSetIndex,
            startNightDate=self.startNightDate,
            debug=self.debug,
        )
        (
            productPath,
            self.qc,
            self.products,
            orderPolyTable,
            self.orderPixelTable,
            orderMetaTable,
        ) = detector.get()

        # FAILED TO FIND THE TRACE
        if productPath is None:
            return None

        # UNPACK THE ORDER TABLE
        orderPolyTable, self.orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=productPath
        )

        # ORDER CHECK IN CASE OF POOR CONTINUUM FITTING
        # GET UNIQUE VALUES IN COLUMN
        self.inst = self.skySubtractedFrame.header[self.kw("INSTRUME")]
        if self.arm.upper() == "VIS" and self.inst.upper() == "SOXS":
            keepOrders = []
            uniqueOrders = self.orderPixelTable["order"].unique()
            for o in uniqueOrders:
                mask = self.orderPixelTable["order"] == o
                if self.orderPixelTable.loc[mask][f"{self.axisA}coord_centre"].std() < 100:
                    keepOrders.append(o)
                else:
                    self.log.warning(f"Bad continuum fit to order {o}; this order will not be extracted")
            mask = self.orderPixelTable["order"].isin(keepOrders)
            self.orderPixelTable = self.orderPixelTable.loc[mask]

        # xpd-update-filter-dataframe-column-values

    def extract(self):
        """*extract the full spectrum order-by-order and return FITS Binary table containing order-merged spectrum*

        **Return:**

        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``mergedSpectumDF`` -- path to the FITS binary table containing the merged spectrum
        """
        self.log.debug("starting the ``extract`` method")

        if self.orderPixelTable is None:
            self.log.error("No trace found in the data, optimal extraction cannot be performed.")
            return self.qc, self.products, None, None, None

        import matplotlib.pyplot as plt
        import pandas as pd
        from astropy.table import Table
        import copy
        from contextlib import suppress
        from astropy.io import fits

        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser
        from datetime import datetime
        from soxspipe.commonutils.toolkit import (
            read_spectral_format,
            add_snr_efficiency_qcs,
        )
        from soxspipe.commonutils import dispersion_map_to_pixel_arrays
        import numpy as np
        import scipy.ndimage
        from astropy.stats import sigma_clip
        import skimage.transform as skt

        kw = self.kw
        arm = self.arm

        uniqueOrders = self.orderPixelTable["order"].unique()
        extractions = []

        self.log.print("\n# PERFORMING OPTIMAL SOURCE EXTRACTION (Horne Method)")

        # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        binx = 1
        biny = 1
        try:
            binx = int(self.skySubtractedFrame.header[kw("WIN_BINX")])
            biny = int(self.skySubtractedFrame.header[kw("WIN_BINY")])
        except:
            pass

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax, amins, amaxs = read_spectral_format(
            log=self.log,
            settings=self.settings,
            arm=self.arm,
            dispersionMap=self.dispersionMap,
            extended=False,
            binx=binx,
            biny=biny,
        )

        self.log.print("\tBuilding wavelength, slit-position, flux, error and bad-pixel arrays")

        # ZOOM IN ON THE DATA TO GET MORE PIXELS ACROSS THE SLIT (NEAREST
        # NEIGHBOUR INTERPOLATION)
        zoomFactor = 11
        if self.detectorParams["dispersion-axis"] == "x":
            zoomTuple = (1, zoomFactor)
        else:
            zoomTuple = (zoomFactor, 1)
        if self.detectorParams["dispersion-axis"] == "x":
            output_shape = (
                self.twoDMap["WAVELENGTH"].data.shape[0],
                self.twoDMap["WAVELENGTH"].data.shape[1] * zoomFactor,
            )
        else:
            output_shape = (
                self.twoDMap["WAVELENGTH"].data.shape[0] * zoomFactor,
                self.twoDMap["WAVELENGTH"].data.shape[1],
            )

        # NEAREST NEIGHBOUR INTERPOLATION USING NUMPY REPEAT (FASTER THAN skt.resize)
        if self.detectorParams["dispersion-axis"] == "x":

            def _zoom(arr):
                return np.repeat(arr, zoomFactor, axis=1)

        else:

            def _zoom(arr):
                return np.repeat(arr, zoomFactor, axis=0)

        wlZoom = _zoom(self.twoDMap["WAVELENGTH"].data)
        slitZoom = _zoom(self.twoDMap["SLIT"].data)
        ## THIS IS (D-S) IN HORNE 1986, I.E. THE SKY-SUBTRACTED FLUX
        rawFluxZoom = _zoom(self.skySubtractedFrame.data)
        if self.subtractedFrame:
            skyZoom = _zoom(self.subtractedFrame.data)
        errorZoom = _zoom(self.skySubtractedFrame.uncertainty.array)
        bpmZoom = _zoom(self.skySubtractedFrame.mask)

        def rebin(arr, binx, biny):
            """REBIN 2D ARRAY ARR TO SHAPE NEW_SHAPE BY AVERAGING."""
            new_shape = (arr.shape[0] // binx, binx, arr.shape[1] // biny, biny)
            return arr.reshape(new_shape).mean(axis=(1, 3))

        def initial_sigma_clipping(rawFluxArray, bpmArray, order=None):
            """RUN SOME INITIAL SIGMA CLIPPING TO CATCH MORE BAD PIXELS"""

            import numpy as np
            from astropy.stats import sigma_clip

            newBpm = []

            bpmArray[bpmArray > 0] = 1

            rawFluxMasked = np.ma.array(rawFluxArray, mask=bpmArray)
            rawFluxMasked = sigma_clip(
                rawFluxMasked,
                sigma_lower=2,
                sigma_upper=7,
                maxiters=1,
                cenfunc="mean",
                stdfunc="std",
                axis=1,
            )
            newBpm = rawFluxMasked.mask

            self.log.print(
                f"\t\t{sum(sum(newBpm)) - sum(sum(bpmArray))} additional bad pixels found from initial sigma clipping (order:{order})",
            )

            return newBpm

        # ADD SOME DATA TO THE SLICES
        orderSlices = []
        wlMinMax = []
        # uniqueOrders = [16]
        for order, amin, amax, wlmin, wlmax in zip(orderNums, amins, amaxs, waveLengthMin, waveLengthMax):
            if order in uniqueOrders:
                orderTable = self.orderPixelTable.loc[
                    (self.orderPixelTable["order"] == order)
                    & (self.orderPixelTable[f"{self.axisB}coord"] > amin)
                    & (self.orderPixelTable[f"{self.axisB}coord"] < amax)
                ]

                self.axisAstart = (
                    np.round((orderTable[f"{self.axisA}coord_centre"] * zoomFactor)).astype(int)
                    - self.slitHalfLength * zoomFactor
                )
                self.axisAstop = (
                    np.round((orderTable[f"{self.axisA}coord_centre"] * zoomFactor)).astype(int)
                    + self.slitHalfLength * zoomFactor
                )
                self.axisBcoord = orderTable[f"{self.axisB}coord"].round().astype(int)
                self.axisAcoords = list(
                    map(
                        lambda x: list(range(x[0], x[1])),
                        zip(self.axisAstart, self.axisAstop),
                    )
                )
                self.axisBcoords = list(
                    map(
                        lambda x: [x] * self.slitHalfLength * 2 * zoomFactor,
                        self.axisBcoord,
                    )
                )

                if self.detectorParams["dispersion-axis"] == "x":
                    orderTable["wavelength"] = list(
                        rebin(
                            wlZoom[self.axisBcoords, self.axisAcoords],
                            zoomTuple[0],
                            zoomTuple[1],
                        )
                    )
                    rawFluxRebinned = rebin(
                        rawFluxZoom[self.axisBcoords, self.axisAcoords],
                        zoomTuple[0],
                        zoomTuple[1],
                    )
                    orderTable["sliceRawFlux"] = list(rawFluxRebinned)
                    orderTable["slit"] = list(
                        rebin(
                            slitZoom[self.axisBcoords, self.axisAcoords],
                            zoomTuple[0],
                            zoomTuple[1],
                        )
                    )

                    if self.subtractedFrame:
                        orderTable["sliceSky"] = list(
                            rebin(
                                skyZoom[self.axisBcoords, self.axisAcoords],
                                zoomTuple[0],
                                zoomTuple[1],
                            )
                        )
                    else:
                        orderTable["sliceSky"] = list([0] * len(self.axisBcoords))
                    orderTable["sliceError"] = list(
                        rebin(
                            errorZoom[self.axisBcoords, self.axisAcoords],
                            zoomTuple[0],
                            zoomTuple[1],
                        )
                    )

                    bpmRebinned = rebin(
                        bpmZoom[self.axisBcoords, self.axisAcoords],
                        zoomTuple[0],
                        zoomTuple[1],
                    )
                    newBpm = initial_sigma_clipping(rawFluxRebinned, bpmRebinned, order=order)

                    orderTable["bpMask"] = list(newBpm)

                else:
                    orderTable["wavelength"] = list(
                        rebin(
                            wlZoom[self.axisAcoords, self.axisBcoords],
                            zoomTuple[1],
                            zoomTuple[0],
                        )
                    )
                    rawFluxRebinned = rebin(
                        rawFluxZoom[self.axisAcoords, self.axisBcoords],
                        zoomTuple[1],
                        zoomTuple[0],
                    )
                    orderTable["sliceRawFlux"] = list(rawFluxRebinned)
                    orderTable["slit"] = list(
                        rebin(
                            slitZoom[self.axisAcoords, self.axisBcoords],
                            zoomTuple[1],
                            zoomTuple[0],
                        )
                    )
                    if self.subtractedFrame:
                        orderTable["sliceSky"] = list(
                            rebin(
                                skyZoom[self.axisAcoords, self.axisBcoords],
                                zoomTuple[1],
                                zoomTuple[0],
                            )
                        )
                    else:
                        orderTable["sliceSky"] = list([0] * len(self.axisAcoords))
                    orderTable["sliceError"] = list(
                        rebin(
                            errorZoom[self.axisAcoords, self.axisBcoords],
                            zoomTuple[1],
                            zoomTuple[0],
                        )
                    )
                    bpmRebinned = rebin(
                        bpmZoom[self.axisAcoords, self.axisBcoords],
                        zoomTuple[1],
                        zoomTuple[0],
                    )

                    newBpm = initial_sigma_clipping(rawFluxRebinned, bpmRebinned, order=order)
                    orderTable["bpMask"] = list(newBpm)

                orderSlices.append(orderTable)
                wlMinMax.append((wlmin, wlmax))

        # REDUCE MEMORY USAGE
        del slitZoom, wlZoom, rawFluxZoom, errorZoom, bpmZoom
        try:
            del skyZoom
        except:
            pass

        self.log.print(f"\tExtracting {len(orderSlices)} orders\n\n")

        from fundamentals import fmultiprocess

        if self.debug or self.turnOffMP:
            turnOffMP = True
        else:
            turnOffMP = False

        extractions = fmultiprocess(
            log=self.log,
            function=extract_single_order,
            inputArray=orderSlices,
            poolSize=False,
            timeout=300,
            funclog=self.log,
            ron=self.ron,
            slitHalfLength=self.slitHalfLength,
            clippingSigma=self.clippingSigma,
            clippingIterationLimit=self.clippingIterationLimit,
            globalClippingSigma=self.globalClippingSigma,
            axisA=self.axisA,
            axisB=self.axisB,
            hornePolyOrder=self.recipeSettings["horne-extraction-profile-poly-order"],
            debug=self.debug,
            turnOffMP=turnOffMP,
            progressBar=True,
        )

        updatedExtractions = []
        for e, wlTuple in zip(extractions, wlMinMax):
            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            if e is not None:
                mask = (e["wavelengthMean"] > wlTuple[0]) & (e["wavelengthMean"] < wlTuple[1])
                e = e.loc[mask]
                updatedExtractions.append(e)

        extractions = updatedExtractions

        self.plot_extracted_spectrum_qc(extractions=extractions)

        # MERGE THE ORDER SPECTRA
        extractedOrdersDF = pd.concat(extractions, ignore_index=True)
        if False:
            extractedOrdersDF = self.tune_wavelength_calibration_to_skylines(extractedOrdersDF, arm=arm)

        mergedSpectumDF, orderJoins = self.merge_extracted_orders(extractedOrdersDF)

        self.qc = add_snr_efficiency_qcs(
            log=self.log,
            spectrumDF=mergedSpectumDF,
            qcTable=self.qc,
            orderJoins=orderJoins,
            recipeName=self.recipeName,
            dateObs=self.dateObs,
        )

        if not isinstance(self.products, bool):

            # CONVERT TO FITS BINARY TABLE
            header = copy.deepcopy(self.skySubtractedFrame.header)
            with suppress(KeyError):
                header.pop(kw("DPR_CATG"))
            with suppress(KeyError):
                header.pop(kw("DPR_TYPE"))
            with suppress(KeyError):
                header.pop(kw("DET_READ_SPEED"))
            with suppress(KeyError):
                header.pop(kw("CONAD"))
            with suppress(KeyError):
                header.pop(kw("GAIN"))
            with suppress(KeyError):
                header.pop(kw("RON"))

            # header["HIERARCH " + kw("PRO_TECH")] = header.pop(kw("DPR_TECH"))
            extractedOrdersTable = Table.from_pandas(extractedOrdersDF)
            BinTableHDU = fits.table_to_hdu(extractedOrdersTable)
            header[kw("SEQ_ARM")] = arm
            header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
            header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{arm}".upper()
            priHDU = fits.PrimaryHDU(header=header)
            hduList = fits.HDUList([priHDU, BinTableHDU])

            # DISCRETE ORDERS
            filename = self.filenameTemplate.replace(
                ".fits",
                f"_EXTRACTED_ORDERS{self.noddingSequence}{self.notFlattened}.fits",
            )
            filePath = f"{self.productDir}/{filename}"
            hduList.verify("fix")
            hduList.writeto(filePath, checksum=True, overwrite=True)

            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

            if not self.notFlattened:
                self.products = pd.concat(
                    [
                        self.products,
                        pd.Series(
                            {
                                "soxspipe_recipe": "soxs-stare",
                                "product_label": f"EXTRACTED_ORDERS_TABLE{self.noddingSequence}{self.notFlattened}",
                                "file_name": filename,
                                "file_type": "FITS",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "product_desc": f"Table of the extracted source in each order",
                                "file_path": filePath,
                                "label": "PROD",
                            }
                        )
                        .to_frame()
                        .T,
                    ],
                    ignore_index=True,
                )

            # NOW MERGED SPECTRUM
            filename = self.filenameTemplate.replace(
                ".fits",
                f"_EXTRACTED_MERGED{self.noddingSequence}{self.notFlattened}.fits",
            )
            filePath = f"{self.productDir}/{filename}"

            def _round_scalar_or_quantity(value, decimals):
                if hasattr(value, "value"):
                    value = value.value
                return round(float(value), decimals)

            for col, decimals in [
                ("WAVE", 2),
                ("FLUX_COUNTS", 3),
                ("SNR", 2),
                ("FLUX_DENSITY_COUNTS", 3),
            ]:
                mergedSpectumDF[col] = mergedSpectumDF[col].apply(lambda x: _round_scalar_or_quantity(x, decimals))

            mergedTable = Table.from_pandas(mergedSpectumDF)
            BinTableHDU = fits.table_to_hdu(mergedTable)
            hduList = fits.HDUList([priHDU, BinTableHDU])
            hduList.verify("fix")
            hduList.writeto(filePath, checksum=True, overwrite=True)

            # EXPORTING SPECTRUM IN ASCII FORMAT

            # CHECKING IF WE ARE IN A NODDING SEQUENCE
            if self.noddingSequence or len(self.noddingSequence) > 0 and not self.notFlattened:
                pass
            elif not self.notFlattened:
                # SAVE THE MERGED ASTROPY TABLE TO TXT FILE
                # SAVE THE TABLE stackedSpectrum TO DISK IN ASCII FORMAT
                asciiFilepath = filePath.replace(".fits", f".txt")
                asciiFilename = filename.replace(".fits", f".txt")
                mergedTable["WAVE"] = mergedTable["WAVE"] * 10  # CONVERTING TO ANGSTROMS
                mergedTable["WAVE"].format = "{:.2f}"
                mergedTable.write(asciiFilepath, format="ascii", overwrite=True)

                self.products = pd.concat(
                    [
                        self.products,
                        pd.Series(
                            {
                                "soxspipe_recipe": self.recipeName,
                                "product_label": f"EXTRACTED_MERGED_ASCII{self.notFlattened}",
                                "file_name": asciiFilename,
                                "file_type": "TXT",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "product_desc": f"Ascii version of extracted source spectrum",
                                "file_path": asciiFilepath,
                                "label": "PROD",
                            }
                        )
                        .to_frame()
                        .T,
                    ],
                    ignore_index=True,
                )

                self.products = pd.concat(
                    [
                        self.products,
                        pd.Series(
                            {
                                "soxspipe_recipe": "soxs-stare",
                                "product_label": f"EXTRACTED_MERGED_TABLE{self.noddingSequence}{self.notFlattened}",
                                "file_name": filename,
                                "file_type": "FITS",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "product_desc": f"Table of the extracted, order-merged",
                                "file_path": filePath,
                                "label": "PROD",
                            }
                        )
                        .to_frame()
                        .T,
                    ],
                    ignore_index=True,
                )
            extractionFilepath = filePath
        else:
            extractionFilepath = False

        self.log.debug("completed the ``extract`` method")
        return self.qc, self.products, mergedSpectumDF, orderJoins, extractionFilepath

    def weighted_average(self, group):
        import numpy as np

        group["FLUX_COUNTS"] = (group["flux_resampled"] * np.abs(group["snr"])).sum() / np.abs(group["snr"]).sum()
        # group['wf'] = group['flux_resampled'].mean()

        return group["FLUX_COUNTS"]

    def residual_merge(self, group):
        import numpy as np

        # residual = np.abs(group[0]['flux_resampled']) -  np.abs(group[1]['flux_resampled'])
        # if residual <= 0:
        #     group['wf'] = group[0]['flux_resampled']
        # else:
        #     group['wf'] = group[1]['flux_resampled']
        if len(group) > 1:
            if np.abs(group.iloc[0]["flux_resampled"]) >= np.abs(group.iloc[1]["flux_resampled"]):
                group["FLUX_COUNTS"] = group.iloc[1]["flux_resampled"]
            else:
                group["FLUX_COUNTS"] = group.iloc[0]["flux_resampled"]
        else:
            group["FLUX_COUNTS"] = group.iloc[0]["flux_resampled"]

        return group["FLUX_COUNTS"]

    def tune_wavelength_calibration_to_skylines(self, extractedOrdersDF, arm):
        """Skyline-based wavelength calibration tuning helper.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders
        - ``arm`` -- the arm of the instrument

        **Return:**

        - ``extractedOrdersDF`` -- unchanged input dataframe
        """

        self.log.debug("starting the ``tune_wavelength_calibration_to_skylines`` method")

        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from datetime import datetime

        skylineWave = pd.to_numeric(self.skylinesDF["WAVELENGTH"], errors="coerce").dropna().to_numpy()

        uniqueOrders = np.sort(extractedOrdersDF["order"].unique())
        orderShifts = []
        for o in uniqueOrders:
            mask = extractedOrdersDF["order"] == o
            orderDF = extractedOrdersDF.loc[mask]
            orderDF["wavelengthMean_shifted"] = orderDF["wavelengthMean"]
            if orderDF.empty:
                continue

            for iteration in range(2):

                wave = pd.to_numeric(orderDF["wavelengthMean_shifted"], errors="coerce")
                sky = pd.to_numeric(orderDF["skyFlux"], errors="coerce")
                valid = wave.notna() & sky.notna()
                if not valid.any():
                    continue

                # FIND AND MARK THE PEAKS IN THE SKY SPECTRUM
                from scipy.signal import find_peaks, savgol_filter

                skySmoothed = savgol_filter(sky[valid], window_length=21, polyorder=3)

                # Find peaks using vectorized operations
                peaks, _ = find_peaks(
                    skySmoothed,
                    height=np.nanmedian(skySmoothed) + 1 * np.nanstd(skySmoothed),
                    distance=50,
                )

                waveVals = wave.loc[valid].to_numpy()
                skyVals = skySmoothed
                wmin = np.nanmin(waveVals)
                wmax = np.nanmax(waveVals)
                localSkylines = skylineWave[(skylineWave >= wmin) & (skylineWave <= wmax)]

                # MATCH THE PEAKS TO THE SKYLINES
                peakWave = waveVals[peaks]
                matched_skyline_waves = []
                matched_shifts = []
                for ww in localSkylines:
                    if len(peakWave) == 0:
                        continue
                    idx = np.where(np.abs(peakWave - ww) < 3)[0]
                    if idx.size:
                        for i in idx:
                            matched_skyline_waves.append(ww)
                            matched_shifts.append(ww - peakWave[i])

                from astropy.stats import sigma_clip, mad_std

                # SIGMA-CLIP THE DATA
                masked_matchedShifts = sigma_clip(
                    matched_shifts, sigma_lower=1.5, sigma_upper=1.5, maxiters=3, cenfunc="mean", stdfunc="std"
                )

                clippedWave = np.asarray(matched_skyline_waves)[np.asarray(masked_matchedShifts.mask)]
                clippedShifts = np.asarray(masked_matchedShifts.data)[np.asarray(masked_matchedShifts.mask)]
                goodWave = np.asarray(matched_skyline_waves)[~np.asarray(masked_matchedShifts.mask)]
                goodShifts = np.asarray(masked_matchedShifts.data)[~np.asarray(masked_matchedShifts.mask)]
                medianShift = np.median(goodShifts) if len(goodShifts) > 0 else 0
                if np.isnan(medianShift):
                    medianShift = 0
                orderDF["wavelengthMean_shifted"] = orderDF["wavelengthMean_shifted"] + medianShift

                if iteration == 0:
                    print(medianShift, "nm shift applied to order", o)
                    orderShifts.append(medianShift)
                    utcnow = datetime.utcnow()
                    utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

                    self.qc = pd.concat(
                        [
                            self.qc,
                            pd.Series(
                                {
                                    "soxspipe_recipe": self.recipeName,
                                    "qc_name": f"WL SKY SHIFT O{int(o)}",
                                    "qc_value": round(medianShift, 3),
                                    "qc_comment": "Shift applied to wavelength solution based on skyline matches",
                                    "qc_order": int(o),
                                    "qc_unit": "nm",
                                    "obs_date_utc": self.dateObs,
                                    "reduction_date_utc": utcnow,
                                    "to_header": True,
                                }
                            )
                            .to_frame()
                            .T,
                        ],
                        ignore_index=True,
                    )

                # GENERATE A PLOT TO VISUALISE THE POTENTIAL SHIFTS (SHIFT VS WAVELENGTH OF THE MATCHED SKYLINE)
                if self.debug:

                    fig, (ax1, ax2) = plt.subplots(
                        2,
                        1,
                        figsize=(12, 6),
                        dpi=180,
                        gridspec_kw={"height_ratios": [2, 1]},
                    )

                    # Top: sky spectrum with peaks and skyline markers
                    if iteration == 0:
                        color = "tab:orange"
                    else:
                        color = "tab:green"
                    ax1.plot(waveVals, skyVals, color=color, linewidth=0.7, label="skyFlux")
                    for ww in localSkylines:
                        ax1.axvline(ww, color="grey", alpha=0.25, linewidth=0.6)
                    if len(localSkylines):
                        ax1.plot([], [], color="grey", alpha=0.7, linewidth=0.8, label="skylines")
                    ax1.plot(waveVals[peaks], skyVals[peaks], "x", color="red", label="skyFlux peaks")
                    ax1.set_ylabel("sky flux ($e^{-}$)")
                    ax1.legend(loc="best", fontsize=8)

                    # Bottom: shifts (skyline_wavelength vs shift)
                    ax2.axhline(0, color="k", linestyle="--", linewidth=0.7)
                    if matched_shifts:
                        ax2.scatter(matched_skyline_waves, matched_shifts, c="tab:green", s=30)
                        ax2.scatter(clippedWave, clippedShifts, c="tab:red", s=30, marker="x")
                        # show median shift
                        ax2.axhline(
                            medianShift,
                            color="tab:orange",
                            linestyle=":",
                            linewidth=0.8,
                            label=f"median={medianShift:.3f} nm",
                        )
                        ax2.legend(loc="best", fontsize=8)
                    ax2.set_xlabel("wavelength (nm)")
                    ax2.set_ylabel("shift (nm)")

                    fig.suptitle(f"\nOrder {int(o)} skyFlux and skyline shifts")

                    plt.show()
                    plt.clf()

                    # plt.close(fig)

        for o, s in zip(uniqueOrders, orderShifts):
            if o == 1:
                rShift = s

        for o, s in zip(uniqueOrders, orderShifts):
            mask = extractedOrdersDF["order"] == o
            extractedOrdersDF.loc[mask, "wavelengthMean"] = extractedOrdersDF.loc[mask, "wavelengthMean"] + s
            if s == 0 and arm.upper() == "VIS" and o in [2, 3]:
                extractedOrdersDF.loc[mask, "wavelengthMean"] = extractedOrdersDF.loc[mask, "wavelengthMean"] + rShift
                # orderDF["wavelengthMean"] -= rShift * 2.5
            if orderDF.empty:
                continue

        return extractedOrdersDF

    def merge_extracted_orders(self, extractedOrdersDF):
        """*merge the extracted order spectra in one continuous spectrum*

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders

        **Return:**

        - None
        """

        self.log.debug("starting the ``merge_extracted_orders`` method")

        # THINGS TO TRY
        # - experiment with FluxConserving, Linear and Spline resampling
        # - need to dynamically define the ends of each order (possibly from order centre traces) and don't extract beyong these points
        # - set a S/N threshold, below which the data point is ignored
        # - run some kind of median smoothing to remove obvious spikes (with higher resolution than spectrograph)

        self.log.print(f"\n# MERGING ORDERS INTO SINGLE SPECTRUM")

        import numpy as np
        from astropy.table import Table
        from specutils.manipulation import (
            FluxConservingResampler,
            LinearInterpolatedResampler,
        )
        from specutils import Spectrum1D
        import astropy.units as u
        import pandas as pd
        from astropy.io import fits
        import matplotlib.pyplot as plt

        import matplotlib
        from datetime import datetime
        from astropy.nddata import VarianceUncertainty
        from specutils.manipulation import median_smooth
        from astropy.stats import sigma_clipped_stats
        from scipy.interpolate import interp1d
        from soxspipe.commonutils.toolkit import calculate_rolling_snr
        from soxspipe.commonutils.toolkit import plot_merged_spectrum_qc

        # ASTROPY HAS RESET LOGGING LEVEL -- FIX
        import logging

        logging.getLogger().setLevel(logging.INFO + 5)

        kw = self.kw
        arm = self.arm

        # PARAMETERS FROM INPUT FILE
        # THIS IS THE STEP SIZE IN NM (0.06 nm IS SIMILAR TO XSHOOTER EXTRACTION)
        if self.arm.upper() == "NIR":
            stepWavelengthOrderMerge = 0.06
        elif self.arm.upper() == "UVB":
            stepWavelengthOrderMerge = 0.02
        elif self.arm.upper() == "VIS":
            stepWavelengthOrderMerge = 0.02

        ratio = 1 / stepWavelengthOrderMerge
        order_list = []

        # A FIX FOR SORTING OF SOXS VIS ORDERS
        if self.arm.upper() in ["VIS"]:
            mask = extractedOrdersDF["order"] == 4
            extractedOrdersDF.loc[mask, "order"] = 11
            mask = extractedOrdersDF["order"] == 3
            extractedOrdersDF.loc[mask, "order"] = 13
            mask = extractedOrdersDF["order"] == 2
            extractedOrdersDF.loc[mask, "order"] = 14
            mask = extractedOrdersDF["order"] == 1
            extractedOrdersDF.loc[mask, "order"] = 12
            extractedOrdersDF["order"] = extractedOrdersDF["order"] - 10

        uniqueOrders = np.sort(extractedOrdersDF["order"].unique())

        lastOrderMin = False
        orderJoins = {}
        orderGaps = {}
        for o in uniqueOrders:
            o = int(o)
            mask = extractedOrdersDF["order"] == o
            if lastOrderMin:
                order_join_wl = (extractedOrdersDF.loc[mask]["wavelengthMean"].max() + lastOrderMin) / 2.0
                orderJoins[f"{o-1}{o}"] = order_join_wl
                gap = extractedOrdersDF.loc[mask]["wavelengthMean"].max() - lastOrderMin
                orderGaps[f"{o-1}{o}"] = gap
                # print(f"ORDER: {o}, JOIN: {order_join_wl}, GAP: {gap}")
            lastOrderMin = extractedOrdersDF.loc[mask]["wavelengthMean"].min()

        # A FIX FOR THE UV XSHOOTER DATA (FLATS TAKEN WITH DIFFERENT LAMPS)
        # GET UNIQUE VALUES IN COLUMN
        if self.arm.upper() == "UVB":
            orderJoins["2021"] = 363.5

        stepRatio = 20

        for o in uniqueOrders:
            o = int(o)
            thisKey = f"{o-1}{o}"
            if thisKey in orderJoins.keys():
                mask = extractedOrdersDF["order"] == o - 1
                gap = orderGaps[f"{o-1}{o}"]
                if gap > stepWavelengthOrderMerge * stepRatio * 2.1:
                    maxwl = orderJoins[thisKey] - stepWavelengthOrderMerge * stepRatio
                    mask = (extractedOrdersDF["order"] == o - 1) & (
                        extractedOrdersDF["wavelengthMean"]
                        <= orderJoins[thisKey] - stepWavelengthOrderMerge * stepRatio
                    )
                    extractedOrdersDF = extractedOrdersDF.loc[~mask]

                if f"{o}{o+1}" in orderGaps.keys():
                    gap = orderGaps[f"{o-1}{o}"]
                if gap > stepWavelengthOrderMerge * stepRatio * 2.1:
                    minwl = orderJoins[thisKey] + stepWavelengthOrderMerge * stepRatio
                    mask = (extractedOrdersDF["order"] == o) & (
                        extractedOrdersDF["wavelengthMean"] > orderJoins[thisKey] + stepWavelengthOrderMerge * stepRatio
                    )
                    extractedOrdersDF = extractedOrdersDF.loc[~mask]

        if self.arm.upper() == "UVB":
            # CLIP DICHROICH REGION
            mask = extractedOrdersDF["wavelengthMean"] > 556
            extractedOrdersDF = extractedOrdersDF.loc[~mask]

        # SORT BY COLUMN NAME
        extractedOrdersDF.sort_values(["wavelengthMean"], ascending=[True], inplace=True)

        # DEFINE THE WAVELENGTH ARRAY
        # ENSURE THE COLUMN IS NOT EMPTY AND CONTAINS VALID NUMERIC VALUES
        if extractedOrdersDF["wavelengthMean"].notna().any():
            min_wavelength = np.min(extractedOrdersDF["wavelengthMean"])
            max_wavelength = np.max(extractedOrdersDF["wavelengthMean"])

            # ENSURE THE RANGE IS VALID
            if min_wavelength < max_wavelength and stepWavelengthOrderMerge > 0:
                start = float(format(min_wavelength * ratio, ".0f")) / ratio
                stop = float(format(max_wavelength * ratio, ".0f")) / ratio

                # DEFINE THE WAVELENGTH ARRAY
                wave_resample_grid = np.arange(start, stop, step=stepWavelengthOrderMerge)
            else:
                raise ValueError("Invalid range or step size for wavelength resampling.")
        else:
            raise ValueError("The 'wavelengthMean' column is empty or contains invalid values.")
        # wave_resample_grid = np.arange(float(format(np.min(extractedOrdersDF['wavelengthMean']) * ratio, '.0f')) / ratio, float(
        #     format(np.max(extractedOrdersDF['wavelengthMean']) * ratio, '.0f')) / ratio, step=stepWavelengthOrderMerge)

        # ADD UNITS TO THE VARIOUS COLUMNS
        extractedOrdersDF["extractedFluxOptimal"] = extractedOrdersDF["extractedFluxOptimal"].values * u.electron
        extractedOrdersDF["extractedFluxBoxcar"] = extractedOrdersDF["extractedFluxBoxcar"].values * u.electron
        extractedOrdersDF["extractedFluxBoxcarRobust"] = (
            extractedOrdersDF["extractedFluxBoxcarRobust"].values * u.electron
        )
        extractedOrdersDF["pixelScaleNm"] = extractedOrdersDF["pixelScaleNm"].values * u.nm
        extractedOrdersDF["wavelengthMean"] = extractedOrdersDF["wavelengthMean"].values * u.nm

        # SAVE THE FLUX DENSITY COUNTS BY DIVING FOR THE PROPER PIXEL SIZE IN WAVELENGTH SPACE - COUNTS ARE THEN E-/NM
        extractedOrdersDF["extractedFluxDensityOptimal"] = (
            extractedOrdersDF["extractedFluxOptimal"].values / extractedOrdersDF["pixelScaleNm"].values
        )
        extractedOrdersDF["extractedFluxDensityBoxcar"] = (
            extractedOrdersDF["extractedFluxBoxcar"].values / extractedOrdersDF["pixelScaleNm"].values
        )
        extractedOrdersDF["extractedFluxDensityBoxcarRobust"] = (
            extractedOrdersDF["extractedFluxBoxcarRobust"].values / extractedOrdersDF["pixelScaleNm"].values
        )

        flux_orig = extractedOrdersDF["extractedFluxOptimal"].values
        spectrum_orig = Spectrum1D(
            flux=flux_orig,
            spectral_axis=extractedOrdersDF["wavelengthMean"].values,
            uncertainty=VarianceUncertainty(extractedOrdersDF["varianceSpectrum"].values),
            bin_specification="center",
        )

        fluxDensity_orig = extractedOrdersDF["extractedFluxDensityOptimal"].values
        spectrumFD_orig = Spectrum1D(
            flux=fluxDensity_orig,
            spectral_axis=extractedOrdersDF["wavelengthMean"].values,
            bin_specification="center",
        )

        sky_orig = extractedOrdersDF["skyFlux"].values * u.electron
        spectrumSky_orig = Spectrum1D(
            flux=sky_orig,
            spectral_axis=extractedOrdersDF["wavelengthMean"].values,
            bin_specification="center",
        )

        # Fast linear resampling — np.interp is O(N log M) vs FluxConservingResampler's
        # O(N·M), and is accurate for output grids with similar resolution to the input.
        wave_in = spectrum_orig.spectral_axis.to_value(u.nm)
        flux_in = spectrum_orig.flux.to_value(u.electron)
        var_in = spectrum_orig.uncertainty.array
        sky_in = spectrumSky_orig.flux.to_value(u.electron)
        fluxDensity_in = spectrumFD_orig.flux.to_value(u.electron / u.nm)

        # Deduplicate wavelengths that may overlap at order join boundaries
        _, unique_idx = np.unique(wave_in, return_index=True)
        flux_resampled = np.interp(
            wave_resample_grid,
            wave_in[unique_idx],
            flux_in[unique_idx],
            left=0.0,
            right=0.0,
        ).astype(np.float32)
        fluxDensity_resampled = np.interp(
            wave_resample_grid,
            wave_in[unique_idx],
            fluxDensity_in[unique_idx],
            left=0.0,
            right=0.0,
        ).astype(np.float32)
        variance_resampled = np.interp(
            wave_resample_grid,
            wave_in[unique_idx],
            var_in[unique_idx],
            left=0.0,
            right=0.0,
        ).astype(np.float32)
        sky_resampled = np.interp(
            wave_resample_grid,
            wave_in[unique_idx],
            sky_in[unique_idx],
            left=0.0,
            right=0.0,
        ).astype(np.float32)

        # flux_resampled = median_smooth(flux_resampled, width=3)
        merged_orders = pd.DataFrame()
        merged_orders["WAVE"] = wave_resample_grid * u.nm
        merged_orders["FLUX_COUNTS"] = flux_resampled * u.electron
        merged_orders["VARIANCE"] = variance_resampled * u.electron
        merged_orders["SKY_COUNTS"] = sky_resampled * u.electron
        merged_orders["SNR"] = merged_orders["FLUX_COUNTS"].values / np.sqrt(merged_orders["VARIANCE"].values)

        # merged_orders = calculate_rolling_snr(dataframe=merged_orders, flux_column="FLUX_COUNTS", window_size=300)

        merged_orders["FLUX_DENSITY_COUNTS"] = fluxDensity_resampled * u.electron / u.nm

        if False:
            self.products, filePath = plot_merged_spectrum_qc(
                merged_orders=merged_orders,
                products=self.products,
                log=self.log,
                qcDir=self.qcDir,
                filenameTemplate=self.filenameTemplate,
                noddingSequence=self.noddingSequence,
                dateObs=self.dateObs,
                arm=self.arm,
                recipeName=self.recipeName,
                orderJoins=orderJoins,
                debug=self.debug,
            )

        return merged_orders, orderJoins

    def plot_extracted_spectrum_qc(self, extractions):
        """*plot extracted spectrum QC plot*

        **Key Arguments:**

        - ``extractions`` -- dataframes hosting order extractions

        **Usage:**

        ```python
        optimalExtractor.plot_extracted_spectrum_qc(extractions)
        ```

        """
        self.log.debug("starting the ``plot_extracted_spectrum_qc`` method")

        from astropy.stats import sigma_clip

        # DO NOT PLOT IF PRODUCT TABLE HAS NOT BEEN PASSED
        # if isinstance(self.products, bool) and not self.products:
        #     return

        import matplotlib.pyplot as plt

        from datetime import datetime
        import pandas as pd
        from astropy.stats import sigma_clipped_stats

        fig = plt.figure(figsize=(14, 12), constrained_layout=True, dpi=180)
        gs = fig.add_gridspec(4, 1)
        toprow = fig.add_subplot(gs[0:1, :])
        secondrow = fig.add_subplot(gs[1:2, :], sharex=toprow)
        thirdrow = fig.add_subplot(gs[2:3, :], sharex=toprow)
        fourthrow = fig.add_subplot(gs[3:4, :], sharex=toprow)
        addedLegend = True

        allExtractions = pd.concat(extractions, ignore_index=True)

        # SIGMA-CLIP THE DATA TO FIND PLOT AXIS LIMITS
        arrayMask = sigma_clip(
            allExtractions["extractedFluxBoxcarRobust"],
            sigma_lower=30000,
            sigma_upper=5.0,
            maxiters=3,
            cenfunc="mean",
            stdfunc="std",
        )
        mean, median, std = sigma_clipped_stats(
            allExtractions["extractedFluxBoxcarRobust"],
            sigma=5.0,
            stdfunc="std",
            cenfunc="mean",
            maxiters=3,
        )
        maxFlux = arrayMask.max() + 0.5 * std

        for df in extractions:

            if not len(df["order"].values):
                continue
            o = df["order"].values[0]

            extracted_wave_spectrum = df["wavelengthMean"]
            extracted_spectrum = df["extractedFluxOptimal"]
            extracted_spectrum_nonopt = df["extractedFluxBoxcarRobust"]
            extracted_variance_spectrum = df["varianceSpectrum"]
            extracted_sky_spectrum = df["skyFlux"]
            extracted_snr = df["snr"]

            try:
                if "PAE" in self.settings and self.settings["PAE"] and True:
                    line = toprow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum[10:-10],
                        zorder=2,
                        linewidth=0.1,
                    )
                    secondrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum_nonopt[10:-10],
                        color=line[0].get_color(),
                        zorder=2,
                        linewidth=0.1,
                    )
                    thirdrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_snr[10:-10],
                        color=line[0].get_color(),
                        zorder=2,
                        linewidth=0.1,
                    )
                    fourthrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_sky_spectrum[10:-10],
                        color=line[0].get_color(),
                        zorder=2,
                        linewidth=0.1,
                    )
                    toprow.text(
                        extracted_wave_spectrum[10:-10].mean(),
                        extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std(),
                        int(o),
                        fontsize=10,
                        c=line[0].get_color(),
                        verticalalignment="bottom",
                    )

                else:
                    if addedLegend:
                        label_opt = "Optimal Extraction"
                        label_nonopt = "Robust Boxcar Extraction"
                        label_snr = "SNR"
                        label_sky = "Sky Flux"
                    else:
                        label_opt = None
                        label_nonopt = None
                        label_snr = None
                        label_sky = None
                    line = toprow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum[10:-10],
                        alpha=0.8,
                        zorder=2,
                        label=label_opt,
                        linewidth=0.5,
                    )
                    secondrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum_nonopt[10:-10],
                        color=line[0].get_color(),
                        alpha=0.8,
                        zorder=1,
                        label=label_nonopt,
                        linewidth=0.5,
                    )
                    thirdrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_snr[10:-10],
                        color=line[0].get_color(),
                        alpha=0.8,
                        zorder=2,
                        label=label_snr,
                        linewidth=0.5,
                    )
                    fourthrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_sky_spectrum[10:-10],
                        color=line[0].get_color(),
                        alpha=0.8,
                        zorder=2,
                        label=label_sky,
                        linewidth=0.5,
                    )
                    if extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std() < maxFlux:
                        toprow.text(
                            extracted_wave_spectrum[10:-10].mean(),
                            extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std(),
                            int(o),
                            fontsize=10,
                            c=line[0].get_color(),
                            verticalalignment="bottom",
                        )
                    addedLegend = False
            except:
                self.log.warning(f"Order skipped: {o}")

        if (
            self.skylinesDF is not None
            and not self.skylinesDF.empty
            and {"WAVELENGTH", "FLUX"}.issubset(self.skylinesDF.columns)
        ):
            xlim = fourthrow.get_xlim()
            for wavelength in self.skylinesDF["WAVELENGTH"].to_numpy():
                fourthrow.axvline(
                    wavelength,
                    ymin=0,
                    ymax=1,
                    color="grey",
                    alpha=0.2,
                    linewidth=0.5,
                    zorder=0,
                )
            fourthrow.set_xlim(xlim)

        toprow.set_title(f"Extracted Object Spectrum QC ({self.arm.upper()})", fontsize=11)
        toprow.set_ylabel("flux ($e^{-}$)", fontsize=10)
        secondrow.set_ylabel("boxcar flux ($e^{-}$)", fontsize=10)
        thirdrow.set_ylabel("SNR", fontsize=10)
        fourthrow.set_ylabel("sky flux ($e^{-}$)", fontsize=10)

        # toprow.set_yscale("log")
        # secondrow.set_yscale("log")
        fourthrow.set_xlabel(f"wavelength (nm)", fontsize=10)

        toprow.legend(fontsize=8, loc="best")
        secondrow.legend(fontsize=8, loc="best")
        thirdrow.legend(fontsize=8, loc="best")
        fourthrow.legend(fontsize=8, loc="best")

        toprow.set_ylim(-200, maxFlux)
        secondrow.set_ylim(-200, maxFlux)
        filename = self.filenameTemplate.replace(
            ".fits",
            f"_EXTRACTED_ORDERS_QC_PLOT{self.noddingSequence}{self.notFlattened}.pdf",
        )
        filePath = f"{self.qcDir}/{filename}"
        # plt.tight_layout()
        # import matplotlib

        # matplotlib.use("MacOSX")
        # plt.show()
        plt.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
        plt.close(fig)
        plt.close("all")

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        if not self.notFlattened and not isinstance(self.products, bool):
            self.products = pd.concat(
                [
                    self.products,
                    pd.Series(
                        {
                            "soxspipe_recipe": "soxs-stare",
                            "product_label": f"EXTRACTED_ORDERS_QC_PLOT{self.noddingSequence}{self.notFlattened}",
                            "file_name": filename,
                            "file_type": "PDF",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "product_desc": f"QC plot of extracted source",
                            "file_path": filePath,
                            "label": "QC",
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )

        self.log.debug("completed the ``plot_extracted_spectrum_qc`` method")
        return None


def extract_single_order(
    crossDispersionSlices,
    funclog,
    ron,
    slitHalfLength,
    clippingSigma,
    clippingIterationLimit,
    globalClippingSigma,
    axisA,
    axisB,
    hornePolyOrder=3,
    debug=False,
):
    """
    *extract the object spectrum for a single order*

    **Return:**

    - ``crossDispersionSlices`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    """
    # log.debug('starting the ``extract_single_order`` method')

    import yaml
    import pandas as pd
    import numpy as np
    from astropy.stats import sigma_clip
    import matplotlib.pyplot as plt

    # WE ARE BUILDING A SET OF CROSS-SLIT OBJECT PROFILES
    # ALONG THE DISPERSION AXIS
    crossSlitProfiles = []

    # 1) SELECTING THE ORDER FROM THE ORDER PIXEL TABLE - THIS IS THE CONTINUUM OF THE OBJECT
    order = crossDispersionSlices["order"].values[0]

    # SLICE SINGLE IMAGE INTO CROSS-DISPERSION SLICES
    crossDispersionSlices = create_cross_dispersion_slices(crossDispersionSlices=crossDispersionSlices)

    # RETURN IF NO SLICES WERE CREATED
    if not len(crossDispersionSlices.index):
        return None

    # VERTICALLY STACK THE SLICES INTO PSEUDO-RECTIFIED IMAGE
    keys_to_stack = ["sliceRawFlux", "sliceFluxNormalised", "sliceMask"]

    # PREALLOCATE A DICTIONARY FOR THE STACKED ARRAYS
    stacked_images = {key: np.vstack(crossDispersionSlices[key]) for key in keys_to_stack}

    # ASSIGN THE STACKED ARRAYS TO VARIABLES
    fluxRawImage = stacked_images["sliceRawFlux"]
    # THIS IS A FIRST GUESS AT P IN HORNE'S ALGORITHM (EQUATION 5) - WE WILL ITERATIVELY FIT THIS PROFILE AND CLIP OUTLIERS TO IMPROVE IT
    fluxNormalisedImage = stacked_images["sliceFluxNormalised"]
    maskImage = stacked_images["sliceMask"]
    # errorImage = stacked_images["sliceError"]
    # bpMaskImage = stacked_images["bpMask"]
    # wavelengthImage = stacked_images["wavelength"]
    # slitImage = stacked_images["slit"]

    ## SKY-SUBTRACTED FLUX RECTIFIED IMAGE
    fluxRawImageMasked = np.ma.masked_array(fluxRawImage, mask=maskImage)
    fluxRawImageMasked = fluxRawImageMasked.filled(np.nan)

    # PLOT THE RECTIFIED IMAGES
    if debug:

        from astropy.io import fits

        hdu = fits.PrimaryHDU(data=fluxRawImageMasked.T)

        # Add Wavelength and Slit Position to the header
        hdu.header["WAVELENGTH"] = (300, "Starting Wavelength (Angstroms)")
        hdu.header["SLITPOS"] = (-5, "Starting Slit Position")

        # Save the FITS file
        hdu.writeto(f"/tmp/spectrum_image_{order}.fits", overwrite=True)

        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True,
        )
        fig.suptitle(f"Raw Flux for Order {order}", fontsize=16)
        plt.imshow(fluxRawImage.T, interpolation="none", aspect="auto")
        plt.show()

        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True,
        )
        fig.suptitle(f"Normalise Flux for Order {order}", fontsize=16)
        plt.imshow(fluxNormalisedImage.T, interpolation="none", aspect="auto")
        plt.show()

        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True,
        )
        fig.suptitle(f"Mask for Order {order}", fontsize=16)
        plt.imshow(maskImage.T, interpolation="none", aspect="auto")
        plt.show()

        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True,
        )
        fig.suptitle(f"Raw Flux & Mask for Order {order}", fontsize=16)
        plt.imshow(fluxRawImageMasked.T, interpolation="none", aspect="auto")
        plt.show()

    # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
    # ITERATE FIRST PIXEL IN EACH SLICE AND THEN MOVE TO NEXT
    crossDispersionSlices = fit_profiles_along_dispersion_axis_for_each_slt_position(
        crossDispersionSlices=crossDispersionSlices,
        fluxNormalisedImage=fluxNormalisedImage,
        maskImage=maskImage,
        slitHalfLength=slitHalfLength,
        clippingSigma=clippingSigma,
        clippingIterationLimit=clippingIterationLimit,
        hornePolyOrder=hornePolyOrder,
        axisB=axisB,
        order=order,
        debug=debug,
        plt=plt,
    )

    # VARIANCE AS PROPAGATED THROUGH THE REDUCTION CASCADE

    # VARIANCE REJECTION NUMBER FROM HORNE 86 PAPER -- XXX THIS IS NOT CLIPPING IN THE CORRECT PLACE - NEEDS TO CLIP AFTER OPTIMAL FLUX CALCULATED
    if False:
        crossDispersionSlices["sliceRawFlux"] * crossDispersionSlices["sliceFittedProfileNormalised"]
        fitted_flux = (
            crossDispersionSlices["sliceRawFluxMaskedSum"] * crossDispersionSlices["sliceFittedProfileNormalised"]
        )
        residuals = crossDispersionSlices["sliceRawFlux"] - fitted_flux
        sliceRejection = np.square(residuals) / crossDispersionSlices["sliceVariance"]
        # CREATE A MASK FOR PIXELS WHERE VARIANCE IS TOO HIGH
        # XXX DO I NEED THIS MASK?
        mask = np.zeros_like(np.stack(sliceRejection, axis=0))
        mask[np.stack(sliceRejection, axis=0) > globalClippingSigma] = 1
        print(sum(mask.flatten()), "pixels masked in global clipping")
        flipMask = 1 - mask
        # COUNTS HOW MANY UNMASKED PIXELS REMAIN IN EACH ROW AFTER CLIPPING
        goodRowCounts = flipMask.sum(axis=1).astype(int)

        # 1D ARRAY OF GOOD VALUES RETURNED, & THEN RESHAPE INTO SAME LENGTH AS DATAFRAME (SPLIT BACK INTO CROSS-DISPERSION SLICES)
        # APPLY THE SAME BAD-PIXEL MASK TO OTHER STACKED ARRAYS TO CREATE "GOOD" VERSIONS OF THESE ARRAYS
        oneDGood = np.ma.masked_array(
            np.stack(crossDispersionSlices["sliceFittedProfileNormalised"], axis=0),
            mask,
        ).compressed()
        crossDispersionSlices["sliceFittedProfileNormalisedGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])
        oneDGood = np.ma.masked_array(np.stack(crossDispersionSlices["sliceRawFlux"], axis=0), mask).compressed()
        crossDispersionSlices["sliceRawFluxGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])
        oneDGood = np.ma.masked_array(np.stack(crossDispersionSlices["sliceVariance"], axis=0), mask).compressed()
        crossDispersionSlices["sliceVarianceGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])
        oneDGood = np.ma.masked_array(np.stack(crossDispersionSlices["wavelength"], axis=0), mask).compressed()
        crossDispersionSlices["wavelengthGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])

        ## RENORMALISE THE GOOD FLUXES AGAIN
        sliceFittedProfileSums = [x.sum() for x in crossDispersionSlices["sliceFittedProfileNormalised"]]
        crossDispersionSlices["sliceFittedProfileNormalised"] = (
            crossDispersionSlices["sliceFittedProfileNormalised"] / sliceFittedProfileSums
        )

    # CALCULATE HORNE 86 NUMERATOR (EQU 8)
    crossDispersionSlices["horneNumerator"] = (
        crossDispersionSlices["sliceRawFlux"]
        * crossDispersionSlices["sliceFittedProfileNormalised"]
        / crossDispersionSlices["sliceVariance"]
    )
    crossDispersionSlices["horneNumeratorSum"] = [x.sum() for x in crossDispersionSlices["horneNumerator"]]

    # CALCULATE HORNE 86 DENOMINATOR (EQU 8)
    crossDispersionSlices["horneDenominator"] = (
        np.power(crossDispersionSlices["sliceFittedProfileNormalised"], 2) / crossDispersionSlices["sliceVariance"]
    )
    crossDispersionSlices["horneDenominatorSum"] = [x.sum() for x in crossDispersionSlices["horneDenominator"]]
    wavelengthMasks = np.stack(crossDispersionSlices["wavelengthMask"].values)
    crossDispersionSlices["wavelengthMean"] = np.ma.mean(wavelengthMasks, axis=1)

    # CALCULATE THE FINAL EXTRACTED SPECTRA
    crossDispersionSlices["varianceSpectrum"] = 1 / crossDispersionSlices["horneDenominatorSum"]
    crossDispersionSlices["extractedFluxOptimal"] = (
        crossDispersionSlices["horneNumeratorSum"] / crossDispersionSlices["horneDenominatorSum"]
    )
    fluxStack = np.vstack(crossDispersionSlices["sliceRawFlux"])
    crossDispersionSlices["extractedFluxBoxcar"] = fluxStack.sum(axis=1)
    skyFluxStack = np.vstack(crossDispersionSlices["sliceSky"])
    crossDispersionSlices["skyFlux"] = np.median(skyFluxStack, axis=1)
    crossDispersionSlices["extractedFluxBoxcarRobust"] = crossDispersionSlices["sliceRawFluxMaskedSum"]
    crossDispersionSlices["snr"] = crossDispersionSlices["extractedFluxOptimal"] / np.power(
        crossDispersionSlices["varianceSpectrum"], 0.5
    )

    # SORT BY COLUMN NAME
    crossDispersionSlices.sort_values(["wavelengthMean"], ascending=[True], inplace=True)

    # REMOVE 0 WAVELENGTH
    crossDispersionSlices = crossDispersionSlices.loc[crossDispersionSlices["wavelengthMean"] > 0]

    crossDispersionSlices["extractedFluxBoxcarRobust"] = crossDispersionSlices["extractedFluxBoxcarRobust"].astype(
        float
    )
    crossDispersionSlices.dropna(
        how="any",
        subset=[
            "pixelScaleNm",
            "wavelengthMean",
            "extractedFluxOptimal",
            "snr",
            "varianceSpectrum",
            "extractedFluxBoxcarRobust",
        ],
        inplace=True,
    )

    # FILTER DATA FRAME
    # FIRST CREATE THE MASK
    mask = crossDispersionSlices["fullColumnMask"] == False
    crossDispersionSlices = crossDispersionSlices.loc[mask]

    # log.debug('completed the ``extract_single_order`` method')

    return crossDispersionSlices[
        [
            "order",
            f"{axisA}coord_centre",
            f"{axisB}coord",
            "wavelengthMean",
            "pixelScaleNm",
            "varianceSpectrum",
            "snr",
            "extractedFluxOptimal",
            "extractedFluxBoxcar",
            "extractedFluxBoxcarRobust",
            "skyFlux",
        ]
    ]


def create_cross_dispersion_slices(crossDispersionSlices):
    """This function is used to create a single, 1-pixel wide cross-dispersion slices of object data. When applied to the dataframe, a single slice is created for each discrete pixel position in the dispersion direction

    **Key Arguments:**

    - ``crossDispersionSlices`` -- the seed dataframe

    **Return:**

    - ``crossDispersionSlices`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    """

    import numpy as np
    from astropy.stats import sigma_clip

    # VERTICALLY STACK THE SLICES INTO PSEUDO-RECTIFIED IMAGE
    bpMask = np.vstack(crossDispersionSlices["bpMask"])
    sliceRawFlux = np.vstack(crossDispersionSlices["sliceRawFlux"])
    crossDispersionSlices["sliceVariance"] = crossDispersionSlices["sliceError"] ** 2
    sliceVariance = np.vstack(crossDispersionSlices["sliceVariance"])

    # CALCULATE THE PIXEL SCALE BELOW ANY CLIPPED OCCURS
    wlTmp = np.vstack(crossDispersionSlices["wavelength"])
    crossDispersionSlices["pixelScaleNm"] = np.ma.mean(wlTmp, axis=1)
    this = (crossDispersionSlices["pixelScaleNm"].values[2:] - crossDispersionSlices["pixelScaleNm"].values[:-2]) / 2
    this = np.insert(this, 0, np.nan)
    this = np.append(this, np.nan)
    crossDispersionSlices["pixelScaleNm"] = np.abs(this)

    ## REMOVE BAD PIXELS AND COSMIC RAYS FROM THE FLUX ARRAY
    maskedArrays, fullColumnMask, sliceRawFluxMasked, sliceVarianceMasked = _mask_slice_raw_flux(
        sliceRawFlux, sliceVariance, bpMask
    )

    crossDispersionSlices["sliceMask"] = [x.mask for x in maskedArrays]
    crossDispersionSlices["sliceRawFluxMasked"] = maskedArrays
    crossDispersionSlices["fullColumnMask"] = fullColumnMask
    sliceRawFluxMaskedSum = sliceRawFluxMasked.sum(axis=1)
    sliceVarianceMaskedSum = sliceVarianceMasked.sum(axis=1)
    ## THIS IS EQUATION (2) IN HORNE 1986 PAPER - THIS IS THE BOXCAR EXTRACTION FLUX
    crossDispersionSlices["sliceRawFluxMaskedSum"] = [row for row in sliceRawFluxMaskedSum]
    ## THIS IS EQUATION (3) IN HORNE 1986 PAPER - THIS IS THE BOXCAR EXTRACTION VARIANCE
    crossDispersionSlices["sliceVarianceMaskedSum"] = [row for row in sliceVarianceMaskedSum]

    ## THIS IS THE NORMALISED FLUX USED FOR FITTING THE DISPERSION PROFILES - THIS IS THE FRACTIONAL FLUX IN HORNE 1986 PAPER
    sliceFluxNormalised = sliceRawFluxMasked / sliceRawFluxMaskedSum[:, np.newaxis]
    crossDispersionSlices["sliceFluxNormalised"] = [row for row in sliceFluxNormalised]

    # FAIL SAFE FOR BAD WAVELENGTH VALUES - SOME ODD WAVELENGTHS FROM DISPERSION SOLUTION CAN CAUSE PROBLEMS WITH PROFILE FITTING
    wlMask = np.vstack(crossDispersionSlices["wavelength"])
    wlMask[wlMask > 0] = 3
    wlMask[wlMask < 0.1] = 1
    wlMask[wlMask > 2] = 0
    wlArray = np.vstack(crossDispersionSlices["wavelength"])
    maskedImage = np.ma.array(wlArray, mask=wlMask)
    maskedImage = sigma_clip(
        maskedImage,
        sigma_lower=1,
        sigma_upper=1,
        maxiters=3,
        cenfunc="mean",
        stdfunc="std",
        axis=1,
    )
    maskedArrays = [np.ma.masked_array(row, mask=maskedImage.mask[i]) for i, row in enumerate(maskedImage.data)]
    crossDispersionSlices["wavelengthMask"] = maskedArrays

    # REMOVE SLICES CONTAINING ALL NANS - THIS IS THE SKY-SUBTRACTED FLUX. (D-S) IN HORNE 1986 PAPER
    crossDispersionSlices["sliceRawFlux"] = [
        np.nan if np.isnan(x).all() else x for x in crossDispersionSlices["sliceRawFlux"]
    ]
    crossDispersionSlices.dropna(axis="index", how="any", subset=["sliceRawFlux"], inplace=True)

    # REMOVE SLICES WITH FULLY MASKED WAVELENGTH
    crossDispersionSlices["wavelength"] = [
        np.nan if x.mask.sum() > x.mask.shape[0] / 1.1 else x for x in crossDispersionSlices["wavelengthMask"]
    ]
    crossDispersionSlices.dropna(axis="index", how="any", subset=["wavelength"], inplace=True)

    return crossDispersionSlices


def fit_profiles_along_dispersion_axis_for_each_slt_position(
    crossDispersionSlices,
    fluxNormalisedImage,
    maskImage,
    slitHalfLength,
    clippingSigma,
    clippingIterationLimit,
    hornePolyOrder,
    axisB,
    order,
    debug,
    plt,
):
    import numpy as np
    from astropy.stats import sigma_clip

    crossSlitProfiles = []
    wave_coords = crossDispersionSlices[f"{axisB}coord"]

    # DETERMINE LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS
    for slitPixelIndex in range(0, slitHalfLength * 2):

        iteration = 1
        clipped_count = 1

        fractions = fluxNormalisedImage[:, slitPixelIndex]
        wave_px = wave_coords
        mask = maskImage[:, slitPixelIndex]

        # fractions WILL STILL CONTAIN BAD-PIXEL/CRHs SO EXCLUDE PIXELS MASKED IN STEP 1 ABOVE
        a = [fractions, wave_px]
        fractions, wave_px = [np.ma.compressed(np.ma.masked_array(i, mask)) for i in a]

        startCount = len(fractions)
        coeff = []

        while (iteration < clippingIterationLimit) and (clipped_count > 0):
            # FIT A POLYNOMIAL TO THE FRACTIONAL FLUXES
            if len(wave_px):
                coeff = np.polyfit(wave_px, fractions, deg=hornePolyOrder)
            else:
                coeff = []
                break
            residuals = fractions - np.polyval(coeff, wave_px)

            # REMOVE REMAINING OUTLIERS
            masked_residuals = sigma_clip(
                residuals,
                sigma_lower=clippingSigma,
                sigma_upper=clippingSigma,
                maxiters=1,
                cenfunc="mean",
                stdfunc="std",
            )
            # REDUCE ARRAYS TO NON-MASKED VALUES
            a = [fractions, wave_px]
            fractions, wave_px = [np.ma.compressed(np.ma.masked_array(i, masked_residuals.mask)) for i in a]
            clipped_count = startCount - len(fractions)
            percent = (float(clipped_count) / float(startCount)) * 100.0
            # self.log.print(f"\tProfile fitting iteration {iteration}, slice index {slitPixelIndex+1}/{slitHalfLength * 2}. {clipped_count} clipped ({percent:0.2f}%) - ORDER {order}")
            iteration = iteration + 1

        # GENERATE THE FINAL FITTING PROFILE FOR THIS SLIT POSITION
        if len(coeff):
            profile = np.polyval(coeff, wave_coords)
            profile[profile < 0] = 0
        else:
            profile = np.zeros_like(wave_coords)
        crossSlitProfiles.append(profile)

        if debug:
            plt.scatter(wave_px, fractions, alpha=0.2)
            plt.plot(wave_px, np.polyval(coeff, wave_px), color="red")
            plt.title(f"Fitted Profile for Order {order}")
            plt.ylim([-1, 1])
            plt.show()

    crossSlitProfiles = np.array(crossSlitProfiles)
    transposedProfiles = crossSlitProfiles.T.tolist()
    crossDispersionSlices["sliceFittedProfile"] = [np.array(t) for t in transposedProfiles]

    sliceFittedProfileSums = [x.sum() for x in crossDispersionSlices["sliceFittedProfile"]]
    crossDispersionSlices["sliceFittedProfileNormalised"] = (
        crossDispersionSlices["sliceFittedProfile"] / sliceFittedProfileSums
    )

    if debug:
        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True,
        )
        fig.suptitle(f"CROSS SLIT PROFILE {order}", fontsize=16)
        plt.imshow(
            np.vstack(crossDispersionSlices["sliceFittedProfileNormalised"]).T,
            interpolation="none",
            aspect="auto",
        )
        plt.show()

    return crossDispersionSlices


def _mask_slice_raw_flux(sliceRawFlux, sliceVariance, bpMask):
    import numpy as np
    from astropy.stats import sigma_clip

    # BUILDS A 2D DISTANCE-FROM-CENTER WEIGHTING ARRAY .. ONLY USED FOR SIGMA-CLIPPING THE FLUX ARRAY, NOT FOR WEIGHTING THE PROFILE FITTING
    sliceIndexArray = np.tile(
        np.abs(np.arange(sliceRawFlux.shape[1]) - sliceRawFlux.shape[1] // 2) * 50 + 1,
        (sliceRawFlux.shape[0], 1),
    )
    bpMask[bpMask > 1] = 1
    sliceRawFluxMasked = np.ma.array(sliceRawFlux, mask=np.isnan(sliceRawFlux) | bpMask)
    sliceIndexArrayMasked = sigma_clip(
        sliceRawFluxMasked * sliceIndexArray,
        sigma_lower=3000,
        sigma_upper=7,
        maxiters=3,
        cenfunc="mean",
        stdfunc="std",
        axis=1,
    )
    sliceRawFluxMasked.data[sliceIndexArrayMasked.mask] = 0
    sliceVarianceMasked = np.ma.array(sliceVariance, mask=sliceIndexArrayMasked.mask)
    maskedArrays = [
        np.ma.masked_array(row, mask=sliceIndexArrayMasked.mask[i]) for i, row in enumerate(sliceRawFluxMasked.data)
    ]

    # FULL SLICE MASK IF MORE THAN 1(?) PIXEL CLIPPED
    fullColumnMask = []
    for ma in maskedArrays:
        if np.ma.count_masked(ma) > 1:
            ma.mask = True
            fullColumnMask.append(True)
        else:
            fullColumnMask.append(False)

    return maskedArrays, fullColumnMask, sliceRawFluxMasked, sliceVarianceMasked
