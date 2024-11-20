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


os.environ['TERM'] = 'vt100'


# TODO: replace RON value with true value from FITS header
# TODO: include the BPM in create_cross_dispersion_slice (at least)
# TODO: find a more robust solution for when horneDenominatorSum == 0 (all pixels in a slice do not pass variance cuts). See where fudged == true
# TODO: delete weight collection
# TODO: revisit how the wavelength for each slice is calculated ... take from the continuum, or the central 3-5 pixels?
# TODO: replace 'gain' with true gain from fits headers

class horne_extraction(object):
    """
    *perform optimal source extraction using the Horne method (Horne 1986)*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
    - ``skyModelFrame`` -- path to sky model frame
    - ``skySubtractedFrame`` -- path to sky subtracted frame
    - ``twoDMapPath`` -- path to 2D dispersion map image path
    - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products (if False no products are saved to file)
    - ``dispersionMap`` -- the FITS binary table containing dispersion map polynomial
        - ``sofName`` -- the set-of-files filename
        - ``locationSetIndex`` -- the index of the AB cycle locations (nodding mode only). Default *False*

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/master/initialisation.html).

    To initiate a horne_extraction object, use the following:

    ```python
    from soxspipe.commonutils import horne_extraction
    optimalExtractor = horne_extraction(
        log=log,
        skyModelFrame=skyModelFrame,
        skySubtractedFrame=skySubtractedFrame,
        twoDMapPath=twoDMap,
        settings=settings,
        recipeName="soxs-stare",
        qcTable=qc,
        productsTable=products,
        dispersionMap=dispMap,
        sofName=sofName
    )
    qc, products = optimalExtractor.extract()
    ```

    """

    def __init__(
            self,
            log,
            settings,
            recipeSettings,
            skyModelFrame,
            skySubtractedFrame,
            twoDMapPath,
            recipeName=False,
            qcTable=False,
            productsTable=False,
            dispersionMap=False,
            sofName=False,
            locationSetIndex=False

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
        from ccdproc import cosmicray_lacosmic

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

        # DETECTING SEQUENCE AUTOMATICALLY
        try:
            self.noddingSequence = "_A" if int(skySubtractedFrame.header['HIERARCH ESO SEQ CUMOFF Y'] > 0) else "_B"
            if locationSetIndex:
                self.noddingSequence += str(locationSetIndex)
        except:
            self.noddingSequence = ""

        home = expanduser("~")
        self.outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"

        # COLLECT SETTINGS FROM SETTINGS FILE
        self.slitHalfLength = int(self.recipeSettings["horne-extraction-slit-length"] / 2)
        self.clippingSigma = self.recipeSettings["horne-extraction-profile-clipping-sigma"]
        self.clippingIterationLimit = self.recipeSettings["horne-extraction-profile-clipping-iteration-count"]
        self.globalClippingSigma = self.recipeSettings["horne-extraction-profile-global-clipping-sigma"]

        # TODO: replace this value with true value from FITS header
        self.ron = 3.0

        # OPEN THE SKY-SUBTRACTED FRAME
        if isinstance(skySubtractedFrame, CCDData):
            self.skySubtractedFrame = skySubtractedFrame
        else:
            self.skySubtractedFrame = CCDData.read(skySubtractedFrame, hdu=0, unit=u.electron,
                                                   hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS',
                                                   key_uncertainty_type='UTYPE')

        if True and self.recipeSettings["use_lacosmic"]:
            oldCount = self.skySubtractedFrame.mask.sum()
            oldMask = self.skySubtractedFrame.mask.copy()
            self.skySubtractedFrame = cosmicray_lacosmic(self.skySubtractedFrame, sigclip=4.0, gain_apply=False, niter=3, cleantype="meanmask")
            newCount = self.skySubtractedFrame.mask.sum()
            self.skySubtractedFrame.mask = oldMask

        # CHECK SKY MODEL FRAME IS USED (ONLY IN STARE MODE)
        if skyModelFrame == False:
            self.skyModelFrame = None
        else:
            # OPEN THE SKY-MODEL FRAME
            if isinstance(skyModelFrame, CCDData):
                self.skyModelFrame = skyModelFrame
            else:
                self.skyModelFrame = CCDData.read(skyModelFrame, hdu=0, unit=u.electron,
                                                  hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS',
                                                  key_uncertainty_type='UTYPE')

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=self.skySubtractedFrame, show=False, ext='data', stdWindow=3, title=False, surfacePlot=True, saveToPath=False)

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = self.skySubtractedFrame.header[kw("SEQ_ARM")]
        self.dateObs = self.skySubtractedFrame.header[kw("DATE_OBS")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=self.log,
            settings=self.settings
        ).get(self.arm)

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
            self.filenameTemplate = filenamer(
                log=self.log,
                frame=self.skySubtractedFrame,
                settings=self.settings
            )

        home = expanduser("~")
        self.qcDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}/"
        self.qcDir = self.qcDir.replace("//", "/")
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.qcDir):
            os.makedirs(self.qcDir)

        # OPEN AND UNPACK THE 2D IMAGE MAP
        self.twoDMap = fits.open(twoDMapPath)

        # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        binx = 1
        biny = 1
        try:
            binx = int(self.skySubtractedFrame.header[kw("WIN_BINX")])
            biny = int(self.skySubtractedFrame.header[kw("WIN_BINY")])
        except:
            pass

        xdim = int(self.twoDMap[0].data.shape[1] / binx)
        ydim = int(self.twoDMap[0].data.shape[0] / biny)
        xarray = np.tile(np.arange(0, xdim), ydim)
        yarray = np.repeat(np.arange(0, ydim), xdim)

        self.skySubtractedFrame.data[self.skySubtractedFrame.data == 0] = np.nan

        if binx > 1 or biny > 1:
            from astropy.nddata import block_reduce
            self.twoDMap["WAVELENGTH"].data = block_reduce(self.twoDMap["WAVELENGTH"].data, (biny, binx), func=np.mean)
            self.twoDMap["SLIT"].data = block_reduce(self.twoDMap["SLIT"].data, (biny, binx), func=np.mean)
            self.twoDMap["ORDER"].data = block_reduce(self.twoDMap["ORDER"].data, (biny, binx), func=np.mean)

        self.imageMap = pd.DataFrame.from_dict({
            "x": xarray,
            "y": yarray,
            "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().astype(float),
            "slit_position": self.twoDMap["SLIT"].data.flatten().astype(float),
            "order": self.twoDMap["ORDER"].data.flatten().astype(float),
            "flux": self.skySubtractedFrame.data.flatten().astype(float)
        })
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
        mask = (self.imageMap['wavelength'] == 0) & (self.imageMap['slit_position'] == 0)
        self.imageMap = self.imageMap.loc[~mask]

        # xpd-update-filter-dataframe-column-values

        # FIND THE OBJECT TRACE IN EACH ORDER
        detector = detect_continuum(
            log=self.log,
            traceFrame=self.skySubtractedFrame,
            dispersion_map=self.dispersionMap,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            sofName=self.sofName,
            recipeName=self.recipeName,
            qcTable=self.qc,
            productsTable=self.products,
            locationSetIndex=locationSetIndex
        )
        productPath, self.qc, self.products, orderPolyTable, self.orderPixelTable, orderMetaTable = detector.get()

        # UNPACK THE ORDER TABLE
        orderPolyTable, self.orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=productPath)

        # ORDER CHECK INCASE OF POOR CONTINUUM FITTING
        # GET UNIQUE VALUES IN COLUMN
        self.inst = self.skySubtractedFrame.header[self.kw("INSTRUME")]
        if self.arm.upper() == "VIS" and self.inst.upper() == "SOXS":
            keepOrders = []
            uniqueOrders = self.orderPixelTable['order'].unique()
            for o in uniqueOrders:
                mask = (self.orderPixelTable['order'] == o)
                if self.orderPixelTable.loc[mask][f"{self.axisA}coord_centre"].std() < 20:
                    keepOrders.append(o)
                else:
                    self.log.warning(f"Bad continuum fit to order {o}; this order will not be extracted")
            mask = (self.orderPixelTable['order'].isin(keepOrders))
            self.orderPixelTable = self.orderPixelTable.loc[mask]

        # xpd-update-filter-dataframe-column-values

    def extract(self):
        """*extract the full spectrum order-by-order and return FITS Binary table containing order-merged spectrum*

        **Return:**

        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``mergedSpectumDF`` -- path to the FITS binary table containing the merged spectrum
        """
        self.log.debug('starting the ``extract`` method')

        import matplotlib.pyplot as plt
        import pandas as pd
        from astropy.table import Table
        import copy
        from contextlib import suppress
        from astropy.io import fits
        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser
        from datetime import datetime
        from soxspipe.commonutils.toolkit import read_spectral_format
        from soxspipe.commonutils import dispersion_map_to_pixel_arrays
        import numpy as np
        import scipy.ndimage
        from astropy.stats import sigma_clip

        kw = self.kw
        arm = self.arm

        # GET UNIQUE VALUES IN COLUMN
        uniqueOrders = self.orderPixelTable['order'].unique()
        # print("FIX ME")
        # uniqueOrders = [25]
        extractions = []

        self.log.print("\n# PERFORMING OPTIMAL SOURCE EXTRACTION (Horne Method)\n\n")

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
            log=self.log, settings=self.settings, arm=self.arm, dispersionMap=self.dispersionMap, extended=False, binx=binx, biny=biny)

        zoomFactor = 35
        if self.detectorParams["dispersion-axis"] == "x":
            zoomTuple = (1, zoomFactor)
        else:
            zoomTuple = (zoomFactor, 1)
        wlZoom = scipy.ndimage.zoom(self.twoDMap["WAVELENGTH"].data, zoomTuple, order=0)
        rawFluxZoom = scipy.ndimage.zoom(self.skySubtractedFrame.data, zoomTuple, order=0)
        if self.skyModelFrame:
            skyZoom = scipy.ndimage.zoom(self.skyModelFrame.data, zoomTuple, order=0)
        errorZoom = scipy.ndimage.zoom(self.skySubtractedFrame.uncertainty.array, zoomTuple, order=0)
        bpmZoom = scipy.ndimage.zoom(self.skySubtractedFrame.mask, zoomTuple, order=0)

        def rebin(arr, binx, biny):
            """Rebin 2D array arr to shape new_shape by averaging."""
            new_shape = (arr.shape[0] // binx, arr.shape[1] // biny)
            shape = (new_shape[0], arr.shape[0] // new_shape[0], new_shape[1], arr.shape[1] // new_shape[1])
            return arr.reshape(shape).mean(-1).mean(1)

        def initial_sigma_clipping(
                rawFluxArray,
                bpmArray,
                wlArray):
            """Run some initial sigma clipping to catch more bad pixels
            """

            import numpy as np
            from astropy.stats import sigma_clip

            newBpm = []
            for r, m, w in zip(rawFluxArray, bpmArray, wlArray):
                m[m > 0] = 1
                maskedArray = np.ma.array(r, mask=m)
                # SIGMA-CLIP THE DATA TO REMOVE COSMIC/BAD-PIXELS
                newMask = sigma_clip(
                    maskedArray, sigma_lower=2, sigma_upper=5, maxiters=1, cenfunc='mean', stdfunc="std")
                newBpm.append(newMask.mask)

            return np.array(newBpm)

        # ADD SOME DATA TO THE SLICES
        orderSlices = []
        wlMinMax = []
        # uniqueOrders = [16]
        for order, amin, amax, wlmin, wlmax in zip(orderNums, amins, amaxs, waveLengthMin, waveLengthMax):
            if order in uniqueOrders:
                orderTable = self.orderPixelTable.loc[(self.orderPixelTable['order'] == order) & (self.orderPixelTable[f"{self.axisB}coord"] > amin) & (self.orderPixelTable[f"{self.axisB}coord"] < amax)]

                self.axisAstart = np.round((orderTable[f"{self.axisA}coord_centre"] * zoomFactor)).astype(int) - self.slitHalfLength * zoomFactor
                self.axisAstop = np.round((orderTable[f"{self.axisA}coord_centre"] * zoomFactor)).astype(int) + self.slitHalfLength * zoomFactor
                self.axisBcoord = orderTable[f"{self.axisB}coord"].round().astype(int)
                self.axisAcoords = list(map(lambda x: list(range(x[0], x[1])), zip(self.axisAstart, self.axisAstop)))
                self.axisBcoords = list(map(lambda x: [x] * self.slitHalfLength * 2 * zoomFactor, self.axisBcoord))

                if self.detectorParams["dispersion-axis"] == "x":
                    orderTable["wavelength"] = list(rebin(wlZoom[self.axisBcoords, self.axisAcoords], zoomTuple[0], zoomTuple[1]))
                    orderTable["sliceRawFlux"] = list(rebin(rawFluxZoom[self.axisBcoords, self.axisAcoords], zoomTuple[0], zoomTuple[1]))
                    if self.skyModelFrame:
                        orderTable["sliceSky"] = list(rebin(skyZoom[self.axisBcoords, self.axisAcoords], zoomTuple[0], zoomTuple[1]))
                    else:
                        orderTable["sliceSky"] = list([0] * len(self.axisBcoords))
                    orderTable["sliceError"] = list(rebin(errorZoom[self.axisBcoords, self.axisAcoords], zoomTuple[0], zoomTuple[1]))

                    newBpm = initial_sigma_clipping(rawFluxZoom[self.axisBcoords, self.axisAcoords], bpmZoom[self.axisBcoords, self.axisAcoords], wlZoom[self.axisBcoords, self.axisAcoords])
                    orderTable["bpMask"] = list(rebin(newBpm, zoomTuple[0], zoomTuple[1]))

                else:
                    orderTable["wavelength"] = list(rebin(wlZoom[self.axisAcoords, self.axisBcoords], zoomTuple[1], zoomTuple[0]))
                    orderTable["sliceRawFlux"] = list(rebin(rawFluxZoom[self.axisAcoords, self.axisBcoords], zoomTuple[1], zoomTuple[0]))
                    if self.skyModelFrame:
                        orderTable["sliceSky"] = list(rebin(skyZoom[self.axisAcoords, self.axisBcoords], zoomTuple[1], zoomTuple[0]))
                    else:
                        orderTable["sliceSky"] = list([0] * len(self.axisAcoords))
                    orderTable["sliceError"] = list(rebin(errorZoom[self.axisAcoords, self.axisBcoords], zoomTuple[1], zoomTuple[0]))
                    newBpm = initial_sigma_clipping(rawFluxZoom[self.axisAcoords, self.axisBcoords], bpmZoom[self.axisAcoords, self.axisBcoords], wlZoom[self.axisAcoords, self.axisBcoords])
                    orderTable["bpMask"] = list(rebin(newBpm, zoomTuple[1], zoomTuple[0]))

                orderSlices.append(orderTable)
                wlMinMax.append((wlmin, wlmax))

        from fundamentals import fmultiprocess
        extractions = fmultiprocess(log=self.log, function=extract_single_order,
                                    inputArray=orderSlices, poolSize=False, timeout=300, funclog=self.log, ron=self.ron, slitHalfLength=self.slitHalfLength, clippingSigma=self.clippingSigma, clippingIterationLimit=self.clippingIterationLimit, globalClippingSigma=self.globalClippingSigma, axisA=self.axisA, axisB=self.axisB, turnOffMP=True)

        updatedExtractions = []
        for e, wlTuple in zip(extractions, wlMinMax):
            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            if e is not None:
                mask = ((e['wavelengthMedian'] > wlTuple[0]) & (e['wavelengthMedian'] < wlTuple[1]))
                e = e.loc[mask]
                updatedExtractions.append(e)

        extractions = updatedExtractions

        self.plot_extracted_spectrum_qc(extractions=extractions, uniqueOrders=uniqueOrders)

        # MERGE THE ORDER SPECTRA
        extractedOrdersDF = pd.concat(extractions, ignore_index=True)

        mergedSpectumDF = self.merge_extracted_orders(extractedOrdersDF)

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
            filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_ORDERS{self.noddingSequence}.fits")
            filePath = f"{self.outDir}/{filename}"
            hduList.writeto(filePath, checksum=True, overwrite=True)

            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": "soxs-stare",
                "product_label": f"EXTRACTED_ORDERS_TABLE{self.noddingSequence}",
                "file_name": filename,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"Table of the extracted source in each order",
                "file_path": filePath,
                "label": "PROD"
            }).to_frame().T], ignore_index=True)

            # NOW MERGED SPECTRUM
            filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED{self.noddingSequence}.fits")
            filePath = f"{self.outDir}/{filename}"
            mergedTable = Table.from_pandas(mergedSpectumDF)
            BinTableHDU = fits.table_to_hdu(mergedTable)
            hduList = fits.HDUList([priHDU, BinTableHDU])
            hduList.writeto(filePath, checksum=True, overwrite=True)

            # EXPORTING SPECTRUM IN ASCII FORMAT

            # CHECKING IF WE ARE IN A NODDING SEQUENCE
            if self.noddingSequence or len(self.noddingSequence) > 0:
                pass
            else:
                # SAVE THE MERGED ASTROPY TABLE TO TXT FILE
                # SAVE THE TABLE stackedSpectrum TO DISK IN ASCII FORMAT
                asciiFilepath = filePath.replace(".fits", f".txt")
                asciiFilename = filename.replace(".fits", f".txt")
                mergedTable.write(asciiFilepath, format='ascii', overwrite=True)
                self.products = pd.concat([self.products, pd.Series({
                    "soxspipe_recipe": self.recipeName,
                    "product_label": "EXTRACTED_MERGED_ASCII",
                    "file_name": asciiFilename,
                    "file_type": "TXT",
                    "obs_date_utc": self.dateObs,
                    "reduction_date_utc": utcnow,
                    "product_desc": f"Ascii version of extracted source spectrum",
                    "file_path": filePath,
                    "label": "PROD"
                }).to_frame().T], ignore_index=True)

            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": "soxs-stare",
                "product_label": f"EXTRACTED_MERGED_TABLE{self.noddingSequence}",
                "file_name": filename,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"Table of the extracted, order-merged",
                "file_path": filePath,
                "label": "PROD"
            }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``extract`` method')
        return self.qc, self.products, mergedSpectumDF

    def weighted_average(self, group):
        import numpy as np
        group['FLUX_COUNTS'] = (group['flux_resampled'] * np.abs(group['snr'])).sum() / np.abs(group['snr']).sum()
        # group['wf'] = group['flux_resampled'].mean()

        return group['FLUX_COUNTS']

    def residual_merge(self, group):
        import numpy as np
        # residual = np.abs(group[0]['flux_resampled']) -  np.abs(group[1]['flux_resampled'])
        # if residual <= 0:
        #     group['wf'] = group[0]['flux_resampled']
        # else:
        #     group['wf'] = group[1]['flux_resampled']
        if len(group) > 1:
            if np.abs(group.iloc[0]['flux_resampled']) >= np.abs(group.iloc[1]['flux_resampled']):
                group['FLUX_COUNTS'] = group.iloc[1]['flux_resampled']
            else:
                group['FLUX_COUNTS'] = group.iloc[0]['flux_resampled']
        else:
            group['FLUX_COUNTS'] = group.iloc[0]['flux_resampled']

        return group['FLUX_COUNTS']

    def merge_extracted_orders(
            self,
            extractedOrdersDF):
        """*merge the extracted order spectra in one continuous spectrum*

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders

        **Return:**

        - None
        """

        self.log.debug('starting the ``merge_extracted_orders`` method')

        # THINGS TO TRY
        # - experiment with FluxConserving, Linear and Spline resampling
        # - need to dynamically define the ends of each order (possibly from order centre traces) and don't extract beyong these points
        # - set a S/N threshold, below which the data point is ignored
        # - run some kind of median smoothing to remove obvious spikes (with higher resolution than spectrograph)

        self.log.print(f"\n# MERGING ORDERS INTO SINGLE SPECTRUM")

        import numpy as np
        from astropy.table import Table
        from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler
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

        # ASTROPY HAS RESET LOGGING LEVEL -- FIX
        import logging
        logging.getLogger().setLevel(logging.INFO + 5)

        kw = self.kw
        arm = self.arm

        # PARAMETERS FROM INPUT FILE
        # THIS IS THE STEP SIZE IN NM (0.06 nm IS SIMILAR TO XHSOOTER EXTRACTION)
        if self.arm.upper() == "NIR":
            stepWavelengthOrderMerge = 0.06
        elif self.arm.upper() == "UVB":
            stepWavelengthOrderMerge = 0.02
        elif self.arm.upper() == "VIS":
            stepWavelengthOrderMerge = 0.02

        ratio = 1 / stepWavelengthOrderMerge
        order_list = []

        # A FIX FOR SORTING OF SOXS VIS ORDERSK
        mask = (extractedOrdersDF['order'] == 4)
        extractedOrdersDF.loc[mask, 'order'] = 11
        mask = (extractedOrdersDF['order'] == 3)
        extractedOrdersDF.loc[mask, 'order'] = 13
        mask = (extractedOrdersDF['order'] == 2)
        extractedOrdersDF.loc[mask, 'order'] = 14
        mask = (extractedOrdersDF['order'] == 1)
        extractedOrdersDF.loc[mask, 'order'] = 12
        extractedOrdersDF['order'] = extractedOrdersDF['order'] - 10

        # MORE CAREFUL TREATMENT OF UVB ORDER MERGING
        # if self.arm.upper() in ["UVB", "NIR"]:
        uniqueOrders = np.sort(extractedOrdersDF['order'].unique())

        # self.inst = self.skySubtractedFrame.header[self.kw("INSTRUME")]
        # if self.arm.upper() == "VIS" and self.inst.upper() == "SOXS":
        #     uniqueOrders = [4, 1, 3, 2]

        lastOrderMin = False
        orderJoins = {}
        orderGaps = {}
        for o in uniqueOrders:
            mask = (extractedOrdersDF['order'] == o)
            if lastOrderMin:
                order_join_wl = (extractedOrdersDF.loc[mask]['wavelengthMedian'].max() + lastOrderMin) / 2.
                orderJoins[f'{o-1}{o}'] = order_join_wl
                gap = extractedOrdersDF.loc[mask]['wavelengthMedian'].max() - lastOrderMin
                orderGaps[f'{o-1}{o}'] = gap
                # print(f"ORDER: {o}, JOIN: {order_join_wl}, GAP: {gap}")
            lastOrderMin = extractedOrdersDF.loc[mask]['wavelengthMedian'].min()

        # A FIX FOR THE UV XSHOOTER DATA (FLATS TAKEN WITH DIFFERENT LAMPS)
        # GET UNIQUE VALUES IN COLUMN
        if self.arm.upper() == "UVB":
            orderJoins["2021"] = 363.5

        stepRatio = 20

        for o in uniqueOrders:

            thisKey = f'{o-1}{o}'
            if thisKey in orderJoins.keys():
                mask = (extractedOrdersDF['order'] == o - 1)
                gap = orderGaps[f'{o-1}{o}']
                if gap > stepWavelengthOrderMerge * stepRatio * 2.1:
                    mask = ((extractedOrdersDF['order'] == o - 1) & (extractedOrdersDF['wavelengthMedian'] <= orderJoins[thisKey] - stepWavelengthOrderMerge * stepRatio))
                    extractedOrdersDF = extractedOrdersDF.loc[~mask]

                if f'{o}{o+1}' in orderGaps.keys():
                    gap = orderGaps[f'{o-1}{o}']
                if gap > stepWavelengthOrderMerge * stepRatio * 2.1:
                    mask = (extractedOrdersDF['order'] == o)
                    mask = ((extractedOrdersDF['order'] == o) & (extractedOrdersDF['wavelengthMedian'] > orderJoins[thisKey] + stepWavelengthOrderMerge * stepRatio))
                    extractedOrdersDF = extractedOrdersDF.loc[~mask]

        if self.arm.upper() == "UVB":
            # CLIP DICHROICH REGION
            mask = (extractedOrdersDF['wavelengthMedian'] > 556)
            extractedOrdersDF = extractedOrdersDF.loc[~mask]

        # SORT BY COLUMN NAME
        extractedOrdersDF.sort_values(['wavelengthMedian'],
                                      ascending=[True], inplace=True)

        # DEFINE THE WAVELENGTH ARRAY
        wave_resample_grid = np.arange(float(format(np.min(extractedOrdersDF['wavelengthMedian']) * ratio, '.0f')) / ratio, float(format(np.max(extractedOrdersDF['wavelengthMedian']) * ratio, '.0f')) / ratio, step=stepWavelengthOrderMerge)
        wave_resample_grid = wave_resample_grid * u.nm

        # INTERPOLATE THE ORDER SPECTRUM INTO THIS NEW ARRAY WITH A SINGLE STEP SIZE
        if "PAE" in self.settings and self.settings["PAE"] and True:
            flux_orig = extractedOrdersDF['extractedFluxBoxcarRobust'].values * u.electron
        else:
            flux_orig = extractedOrdersDF['extractedFluxOptimal'].values * u.electron
        # PASS ORIGINAL RAW SPECTRUM AND RESAMPLE
        spectrum_orig = Spectrum1D(flux=flux_orig, spectral_axis=extractedOrdersDF['wavelengthMedian'].values * u.nm, uncertainty=VarianceUncertainty(extractedOrdersDF["varianceSpectrum"].values))

        resampler = FluxConservingResampler()
        flux_resampled = resampler(spectrum_orig, wave_resample_grid)
        # flux_resampled = median_smooth(flux_resampled, width=3)
        merged_orders = pd.DataFrame()
        merged_orders['WAVE'] = flux_resampled.spectral_axis
        merged_orders['FLUX_COUNTS'] = flux_resampled.flux

        self.plot_merged_spectrum_qc(merged_orders)

        return merged_orders

    def plot_extracted_spectrum_qc(
            self,
            uniqueOrders,
            extractions):
        """*plot extracted spectrum QC plot*

        **Key Arguments:**

        - ``uniqueOrders`` -- the unique orders of extraction
        - ``extractions`` -- dataframes hosting order extractions

        **Usage:**

        ```python
        optimalExtractor.plot_extracted_spectrum_qc(uniqueOrders, extractions)
        ```

        """
        self.log.debug('starting the ``plot_extracted_spectrum_qc`` method')

        # DO NOT PLOT IF PRODUCT TABLE HAS NOT BEEN PASSED
        if isinstance(self.products, bool) and not self.products:
            return

        import matplotlib.pyplot as plt
        from datetime import datetime
        import pandas as pd
        from astropy.stats import sigma_clipped_stats

        fig = plt.figure(figsize=(7, 7), constrained_layout=True, dpi=320)
        gs = fig.add_gridspec(1, 1)
        toprow = fig.add_subplot(gs[0:1, :])
        addedLegend = True

        allExtractions = pd.concat(extractions, ignore_index=True)

        mean, median, std = sigma_clipped_stats(allExtractions["extractedFluxBoxcarRobust"], sigma=5., stdfunc="mad_std", cenfunc="median", maxiters=3)

        maxFlux = allExtractions['extractedFluxBoxcarRobust'].max() + std
        # if maxFlux > median + 5 * std:
        #     maxFlux = median + 5 * std

        for df, o in zip(extractions, uniqueOrders):

            extracted_wave_spectrum = df["wavelengthMedian"]
            extracted_spectrum = df["extractedFluxOptimal"]
            extracted_spectrum_nonopt = df["extractedFluxBoxcarRobust"]
            extracted_variance_spectrum = df["varianceSpectrum"]
            extracted_snr = df["snr"]

            try:
                if "PAE" in self.settings and self.settings["PAE"] and True:
                    line = toprow.plot(extracted_wave_spectrum[10:-10], extracted_spectrum_nonopt[10:-10], zorder=2, linewidth=0.1)
                    toprow.text(extracted_wave_spectrum[10:-10].mean(), extracted_spectrum_nonopt[10:-10].mean() + 1.4 * extracted_spectrum_nonopt[10:-10].std(), int(o), fontsize=10, c=line[0].get_color(), verticalalignment='bottom')

                else:
                    if addedLegend:
                        label = "Robust Boxcar Extraction"
                    else:
                        label = None
                    toprow.plot(extracted_wave_spectrum[10:-10], extracted_spectrum_nonopt[10:-10], color="gray", alpha=0.8, zorder=1, label=label, linewidth=0.5)
                    # plt.plot(extracted_wave_spectrum, extracted_spectrum_nonopt, color="gray", alpha=0.1, zorder=1)
                    line = toprow.plot(extracted_wave_spectrum[10:-10], extracted_spectrum[10:-10], zorder=2, linewidth=0.5)
                    if extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std() < maxFlux:
                        toprow.text(extracted_wave_spectrum[10:-10].mean(), extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std(), int(o), fontsize=10, c=line[0].get_color(), verticalalignment='bottom')
                    addedLegend = False
            except:
                self.log.warning(f"Order skipped: {o}")

        # toprow.legend(loc='lower right', bbox_to_anchor=(1, -0.5),
        #               fontsize=8)
        toprow.set_title(
            f"Optimally Extracted Object Spectrum ({self.arm.upper()})", fontsize=11)
        toprow.set_ylabel('flux ($e^{-}$)', fontsize=10)
        toprow.set_xlabel(f'wavelength (nm)', fontsize=10)

        plt.ylim(-200, maxFlux)
        # toprow.set_title(
        #     f"Optimally Extracted Object Spectrum ({self.arm.upper()})", fontsize=11)
        filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_ORDERS_QC_PLOT{self.noddingSequence}.pdf")
        filePath = f"{self.qcDir}/{filename}"
        # plt.tight_layout()
        # plt.show()
        plt.savefig(filePath, dpi='figure', bbox_inches='tight')
        plt.close()
        # plt.show()

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": f"EXTRACTED_ORDERS_QC_PLOT{self.noddingSequence}",
            "file_name": filename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"QC plot of extracted source",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``plot_extracted_spectrum_qc`` method')
        return None

    def plot_merged_spectrum_qc(
            self,
            merged_orders):
        """*plot merged spectrum QC plot*

        **Key Arguments:**

        - ``merged_orders`` -- the dataframe containing the merged order spectrum.

        **Usage:**

        ```python
        optimalExtractor.plot_merged_spectrum_qc(merged_orders)
        ```
        """
        self.log.debug('starting the ``plot_merged_spectrum_qc`` method')

        # DO NOT PLOT IF PRODUCT TABLE HAS NOT BEEN PASSED
        if isinstance(self.products, bool) and not self.products:
            return

        import matplotlib.pyplot as plt
        from datetime import datetime
        import pandas as pd
        from astropy.stats import sigma_clipped_stats

        fig = plt.figure(figsize=(7, 7), constrained_layout=True, dpi=320)
        gs = fig.add_gridspec(1, 1)
        toprow = fig.add_subplot(gs[0:1, :])
        # toprow.legend(loc='lower right', bbox_to_anchor=(1, -0.5),
        #               fontsize=8)
        toprow.set_ylabel('flux ($e^{-}$)', fontsize=10)
        toprow.set_xlabel(f'wavelength (nm)', fontsize=10)
        toprow.set_title(
            f"Optimally Extracted Order-Merged Object Spectrum ({self.arm.upper()})", fontsize=11)

        # for order in extractedOrdersDF['order'].unique():
        #     # LIMIT DATAFRAME TO JUST THIS ORDER
        #     orderDF = extractedOrdersDF.loc[extractedOrdersDF['order'] == order]
        #     plt.plot(orderDF['wavelengthMedian'], orderDF['extractedFluxOptimal'])

        mean, median, std = sigma_clipped_stats(merged_orders['FLUX_COUNTS'], sigma=5.0, stdfunc="mad_std", cenfunc="median", maxiters=3)
        plt.plot(merged_orders['WAVE'], merged_orders['FLUX_COUNTS'], linewidth=0.2, color="#dc322f")
        # sys.exit(0)

        maxFlux = merged_orders['FLUX_COUNTS'].max() + std
        # if maxFlux > median + 5 * std:
        #     maxFlux = median + 5 * std

        plt.ylim(-200, maxFlux)
        try:
            plt.xlim(merged_orders['WAVE'].min().value, merged_orders['WAVE'].max().value)
        except:
            plt.xlim(merged_orders['WAVE'].min(), merged_orders['WAVE'].max())

        filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED_QC_PLOT{self.noddingSequence}.pdf")
        filePath = f"{self.qcDir}/{filename}"
        # plt.show()
        plt.savefig(filePath, dpi='figure', bbox_inches='tight')

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": f"EXTRACTED_MERGED_QC_PLOT{self.noddingSequence}",
            "file_name": filename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"QC plot of extracted order-merged source",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``plot_merged_spectrum_qc`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method


def extract_single_order(crossDispersionSlices, funclog, ron, slitHalfLength, clippingSigma, clippingIterationLimit, globalClippingSigma, axisA, axisB):
    """
    *extract the object spectrum for a single order*

    **Return:**

    - ``crossDispersionSlices`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    """
    #log.debug('starting the ``extract_single_order`` method')

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

    # CREATE THE SLICES AND DROP SLICES WITH ALL NANs (TYPICALLY PIXELS WITH NANs IN 2D IMAGE MAP)
    # sys.stdout.flush()
    # sys.stdout.write("\x1b[1A\x1b[2K")
    # log.print(f"\t## SLICING ORDER INTO CROSS-DISPERSION SLICES - ORDER {order}")

    # REMOVE SLICES WITH ALL NANS
    crossDispersionSlices["sliceRawFlux"] = [np.nan if np.isnan(x).all() else x for x in crossDispersionSlices["sliceRawFlux"]]
    crossDispersionSlices.dropna(axis='index', how='any', subset=["sliceRawFlux"], inplace=True)

    # MASK THE MOST DEVIANT PIXELS IN EACH SLICE
    crossDispersionSlices = crossDispersionSlices.apply(lambda x: create_cross_dispersion_slice(x), axis=1)
    crossDispersionSlices["sliceMask"] = [x.mask for x in crossDispersionSlices["sliceRawFluxMasked"]]

    crossDispersionSlices["sliceRawFluxMaskedSum"] = [x.sum() for x in crossDispersionSlices["sliceRawFluxMasked"]]

    # WEIGHTS ARE NOT YET USED
    crossDispersionSlices["sliceWeights"] = ron + np.abs(crossDispersionSlices["sliceRawFlux"]) / (crossDispersionSlices["sliceRawFluxMaskedSum"].pow(2))

    # NORMALISE THE FLUX
    crossDispersionSlices["sliceFluxNormalised"] = crossDispersionSlices["sliceRawFluxMasked"] / crossDispersionSlices["sliceRawFluxMaskedSum"]

    crossDispersionSlices["sliceFluxNormalisedSum"] = [x.sum() for x in crossDispersionSlices["sliceFluxNormalised"]]

    # REMOVE SLICES WITH FULLY MASKED WAVELENGTH
    crossDispersionSlices["wavelength"] = [np.nan if x.mask.sum() > x.mask.shape[0] / 1.1 else x for x in crossDispersionSlices["wavelengthMask"]]
    crossDispersionSlices.dropna(axis='index', how='any', subset=["wavelength"], inplace=True)

    if not len(crossDispersionSlices.index):
        return None

    # VERTICALLY STACK THE SLICES INTO PSEDUO-RECTIFIED IMAGE
    fluxRawImage = np.vstack(crossDispersionSlices["sliceRawFlux"])
    fluxNormalisedImage = np.vstack(crossDispersionSlices["sliceFluxNormalised"])
    weightImage = np.vstack(crossDispersionSlices["sliceWeights"])
    maskImage = np.vstack(crossDispersionSlices["sliceMask"])
    errorImage = np.vstack(crossDispersionSlices["sliceError"])
    bpMaskImage = np.vstack(crossDispersionSlices["bpMask"])
    wavelengthImage = np.vstack(crossDispersionSlices["wavelength"])

    # PLOT THE RECTIFIED IMAGES
    if False:

        # fig = plt.figure(
        #     num=None,
        #     figsize=(135, 1),
        #     dpi=None,
        #     facecolor=None,
        #     edgecolor=None,
        #     frameon=True)
        # plt.imshow(wavelengthImage.T, interpolation='none', aspect='auto')
        # plt.show()

        # fig = plt.figure(
        #     num=None,
        #     figsize=(135, 1),
        #     dpi=None,
        #     facecolor=None,
        #     edgecolor=None,
        #     frameon=True)
        # plt.imshow(maskImage.T, interpolation='none', aspect='auto')
        # plt.show()

        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True)
        plt.imshow(fluxRawImage.T, interpolation='none', aspect='auto')
        plt.show()

        # fig = plt.figure(
        #     num=None,
        #     figsize=(135, 1),
        #     dpi=None,
        #     facecolor=None,
        #     edgecolor=None,
        #     frameon=True)
        # plt.imshow(fluxNormalisedImage.T, interpolation='none', aspect='auto')
        # plt.show()

    # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
    # ITERATE FIRST PIXEL IN EACH SLICE AND THEN MOVE TO NEXT
    # sys.stdout.flush()
    # sys.stdout.write("\x1b[1A\x1b[2K")
    # log.print(f"\t## FITTING CROSS-SLIT FLUX NORMALISED PROFILES - ORDER {order}")
    for slitPixelIndex in range(0, slitHalfLength * 2):

        iteration = 1
        clipped_count = 1

        fractions = fluxNormalisedImage[:, slitPixelIndex]
        wave_px = crossDispersionSlices[f"{axisB}coord"]
        weights = weightImage[:, slitPixelIndex]
        mask = maskImage[:, slitPixelIndex]

        # fractions WILL STILL CONTAIN BAD-PIXEL/CRHs SO EXCLUDE PIXELS MASKED IN STEP 1 ABOVE
        a = [fractions, wave_px, weights]
        fractions, wave_px, weights = [np.ma.compressed(np.ma.masked_array(
            i, mask)) for i in a]

        startCount = len(fractions)
        while (iteration < clippingIterationLimit) and (clipped_count > 0):

            if len(wave_px):
                coeff = np.polyfit(wave_px, fractions, deg=2)
            else:
                break
            residuals = fractions - np.polyval(coeff, wave_px)

            # REMOVE REMAINING OUTLIERS
            masked_residuals = sigma_clip(residuals, sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
            # REDUCE ARRAYS TO NON-MASKED VALUES
            a = [fractions, wave_px, weights]
            fractions, wave_px, weights = [np.ma.compressed(np.ma.masked_array(
                i, masked_residuals.mask)) for i in a]
            clipped_count = startCount - len(fractions)
            percent = (float(clipped_count) / float(startCount)) * 100.
            # self.log.print(f"\tProfile fitting iteration {iteration}, slice index {slitPixelIndex+1}/{slitHalfLength * 2}. {clipped_count} clipped ({percent:0.2f}%) - ORDER {order}")
            iteration = iteration + 1
            # if iteration > 1:
            #     sys.stdout.flush()
            #     sys.stdout.write("\x1b[1A\x1b[2K")

        # GENERATE THE FINAL FITTING PROFILE FOR THIS SLIT POSITION
        profile = np.polyval(coeff, crossDispersionSlices[f"{axisB}coord"])
        # REMOVE -VE VALUE
        profile[profile < 0] = 0
        crossSlitProfiles.append(profile)

        if False:
            plt.scatter(wave_px, fractions, alpha=0.2)
            plt.plot(wave_px, np.polyval(coeff, wave_px), color='red')
            plt.ylim([-1, 1])
            plt.show()

    # ADD THE NORMALISED PROFILES TO DATAFRAME
    crossSlitProfiles = np.array(crossSlitProfiles)
    transposedProfiles = crossSlitProfiles.T.tolist()
    crossDispersionSlices["sliceFittedProfile"] = [np.array(t) for t in transposedProfiles]

    # sys.stdout.flush()
    # sys.stdout.write("\x1b[1A\x1b[2K")
    # log.print(f"\t## EXTRACTING THE SPECTRUM - ORDER {order}")

    # NORMALISE THE FLUX IN EACH SLICE
    sliceFittedProfileSums = [x.sum() for x in crossDispersionSlices["sliceFittedProfile"]]
    crossDispersionSlices["sliceFittedProfileNormalised"] = crossDispersionSlices["sliceFittedProfile"] / sliceFittedProfileSums

    # TODO: USE DETECTOR GAIN
    gain = 1.0

    # VARIANCE FROM HORNE 86 PAPER
    crossDispersionSlices["sliceVariance"] = ron + np.abs(crossDispersionSlices["sliceRawFlux"] * crossDispersionSlices["sliceFittedProfileNormalised"] + crossDispersionSlices["sliceSky"]) / gain
    # VARIANCE REJECTION NUMBER FROM HORNE 86 PAPER
    sliceRejection = np.power((crossDispersionSlices["sliceRawFlux"] - crossDispersionSlices["sliceRawFluxMaskedSum"] * crossDispersionSlices["sliceFittedProfileNormalised"]), 2) / crossDispersionSlices["sliceVariance"]
    # CREATE A MASK FOR PIXELS WHERE VARIANCE IS TOO HIGH
    mask = np.zeros_like(np.stack(sliceRejection, axis=0))
    mask[np.stack(sliceRejection, axis=0) > globalClippingSigma] = 1
    flipMask = 1 - mask
    goodRowCounts = flipMask.sum(axis=1).astype(int)

    # 1D ARRAY OF GOOD VALUES RETURNED - RESHAPE INTO SAME LENGTH AS DATAFRAME - FOR MULTIPLE CALCUATIONS
    oneDGood = np.ma.masked_array(np.stack(crossDispersionSlices["sliceFittedProfileNormalised"], axis=0), mask).compressed()
    crossDispersionSlices["sliceFittedProfileNormalisedGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])
    oneDGood = np.ma.masked_array(np.stack(crossDispersionSlices["sliceRawFlux"], axis=0), mask).compressed()
    crossDispersionSlices["sliceRawFluxGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])
    oneDGood = np.ma.masked_array(np.stack(crossDispersionSlices["sliceVariance"], axis=0), mask).compressed()
    crossDispersionSlices["sliceVarianceGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])
    oneDGood = np.ma.masked_array(np.stack(crossDispersionSlices["wavelength"], axis=0), mask).compressed()
    crossDispersionSlices["wavelengthGood"] = np.split(oneDGood, np.cumsum(goodRowCounts)[:-1])

    # CALCULATE HORNE 86 NUMERATOR AND DENOMINATOR (EQU ?)
    crossDispersionSlices['horneNumerator'] = crossDispersionSlices["sliceRawFluxGood"] * crossDispersionSlices["sliceFittedProfileNormalisedGood"] / crossDispersionSlices["sliceVarianceGood"]
    crossDispersionSlices['horneNumeratorSum'] = [x.sum() for x in crossDispersionSlices["horneNumerator"]]
    crossDispersionSlices["horneDenominator"] = np.power(crossDispersionSlices["sliceFittedProfileNormalisedGood"], 2) / crossDispersionSlices["sliceVarianceGood"]
    crossDispersionSlices['horneDenominatorSum'] = [x.sum() for x in crossDispersionSlices["horneDenominator"]]
    crossDispersionSlices["wavelengthMedian"] = [np.ma.median(x) for x in crossDispersionSlices["wavelengthMask"]]
    crossDispersionSlices["fudged"] = False

    if False:
        # print("FIX ME")
        crossDispersionSlices['sliceRawFluxGoodSum'] = [x.sum() for x in crossDispersionSlices["sliceRawFluxGood"]]
        crossDispersionSlices['sliceFittedProfileNormalisedGoodSum'] = [x.sum() for x in crossDispersionSlices["sliceFittedProfileNormalisedGood"]]
        crossDispersionSlices["sliceVarianceGoodSum"] = [x.sum() for x in crossDispersionSlices["sliceVarianceGood"]]

    # TODO: IF ALL ABOVE sliceVarianceRejectLimit  ... what?
    mask = (crossDispersionSlices["horneDenominatorSum"] == 0)
    crossDispersionSlices.loc[mask, "horneNumerator"] = crossDispersionSlices.loc[mask, "sliceRawFlux"] * crossDispersionSlices.loc[mask, "sliceFittedProfileNormalised"] / crossDispersionSlices.loc[mask, "sliceVariance"]
    crossDispersionSlices.loc[mask, 'horneNumeratorSum'] = [x.sum() for x in crossDispersionSlices.loc[mask, "horneNumerator"]]
    crossDispersionSlices.loc[mask, "horneDenominator"] = np.power(crossDispersionSlices.loc[mask, "sliceFittedProfileNormalised"], 2) / crossDispersionSlices.loc[mask, "sliceVariance"]
    crossDispersionSlices.loc[mask, 'horneDenominatorSum'] = [x.sum() for x in crossDispersionSlices.loc[mask, "horneDenominator"]]
    crossDispersionSlices.loc[mask, "fudged"] = True

    # CALCULATE THE FINAL EXTRACTED SPECTRA
    crossDispersionSlices["varianceSpectrum"] = 1 / crossDispersionSlices["horneDenominatorSum"]
    crossDispersionSlices["extractedFluxOptimal"] = crossDispersionSlices["horneNumeratorSum"] / crossDispersionSlices["horneDenominatorSum"]
    crossDispersionSlices["extractedFluxBoxcar"] = [x.sum() for x in crossDispersionSlices["sliceRawFlux"]]
    crossDispersionSlices["extractedFluxBoxcarRobust"] = crossDispersionSlices["sliceRawFluxMaskedSum"]
    crossDispersionSlices["snr"] = crossDispersionSlices["extractedFluxOptimal"] / np.power(crossDispersionSlices["varianceSpectrum"], 0.5)

    if False:
        import sqlite3 as sql
        # CONNECT TO THE DATABASE

        # REGISTER SQL CONVERTERS
        from astropy.nddata.nduncertainty import StdDevUncertainty
        sql.register_adapter(StdDevUncertainty, lambda arr: str(arr.array.tolist()))
        sql.register_adapter(list, lambda arr: str(arr))
        sql.register_adapter(np.array, lambda arr: str(arr.tolist()))
        sql.register_adapter(np.ndarray, lambda arr: str(arr.tolist()))
        sql.register_adapter(np.float64, lambda this: this.item())
        sql.register_adapter(np.ma.core.MaskedArray, lambda arr: str(arr.tolist()))

        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser
        home = expanduser("~")
        conn = sql.connect(f"{home}/Desktop/pandas_export.db")
        # SEND TO DATABASE
        crossDispersionSlices.to_sql(f'order_{order}', con=conn,
                                     index=False, if_exists='replace')

    # SORT BY COLUMN NAME
    crossDispersionSlices.sort_values(['wavelengthMedian'],
                                      ascending=[True], inplace=True)

    # REMOVE 0 WAVELENGTH
    crossDispersionSlices = crossDispersionSlices.loc[crossDispersionSlices['wavelengthMedian'] > 0]

    # DETERMINE THE PIXEL SCALE FOR EACH PIXEL - NEEDED FOR ORDER MERGING
    this = (crossDispersionSlices['wavelengthMedian'].values[2:] - crossDispersionSlices['wavelengthMedian'].values[:-2]) / 2
    this = np.insert(this, 0, np.nan)
    this = np.append(this, np.nan)
    crossDispersionSlices['pixelScaleNm'] = this

    crossDispersionSlices['extractedFluxBoxcarRobust'] = crossDispersionSlices['extractedFluxBoxcarRobust'].astype(float)
    crossDispersionSlices.dropna(how="any", subset=["pixelScaleNm", "wavelengthMedian", "extractedFluxOptimal", "snr", "varianceSpectrum", "extractedFluxBoxcarRobust"], inplace=True)
    # CONVERT COLUMN TYPE

    # FILTER DATA FRAME
    # FIRST CREATE THE MASK
    mask = (crossDispersionSlices["fullColumnMask"] == False)
    crossDispersionSlices = crossDispersionSlices.loc[mask]

    #log.debug('completed the ``extract_single_order`` method')

    return crossDispersionSlices[['order', f'{axisA}coord_centre', f'{axisB}coord', 'wavelengthMedian', 'pixelScaleNm', 'varianceSpectrum', 'snr', 'extractedFluxOptimal', 'extractedFluxBoxcar', 'extractedFluxBoxcarRobust']]


def create_cross_dispersion_slice(
        series):
    """This function is used to create a single, 1-pixel wide cross-dispersion slice of object data. When applied to the dataframe, a single slice is created for each discrete pixel position in the dispersion direction
    """

    import numpy as np
    from astropy.stats import sigma_clip

    series['bpMask'][series['bpMask'] > 0] = 1
    maskedArray = np.ma.array(series["sliceRawFlux"], mask=series['bpMask'])

    # SIGMA-CLIP THE DATA TO REMOVE COSMIC/BAD-PIXELS
    series["sliceRawFluxMasked"] = sigma_clip(
        maskedArray, sigma_lower=3, sigma_upper=5, maxiters=1, cenfunc='mean', stdfunc="std")

    # series["sliceRawFluxMasked"].data[series["sliceRawFluxMasked"].mask]=series["sliceRawFluxMasked"].data[~series["sliceRawFluxMasked"].mask].median()
    series["sliceRawFluxMasked"].data[series["sliceRawFluxMasked"].mask] = 0

    series["fullColumnMask"] = False
    if np.ma.count_masked(series["sliceRawFluxMasked"]) > 1:
        series["sliceRawFluxMasked"].mask = True
        series["fullColumnMask"] = True

    # SIGMA-CLIP WAVELENGTH
    series['wavelengthMask'] = series['wavelength'].copy()
    series['wavelengthMask'][series['wavelengthMask'] > 0] = 3
    series['wavelengthMask'][series['wavelengthMask'] < 0.1] = 1
    series['wavelengthMask'][series['wavelengthMask'] > 2] = 0
    maskedArray = np.ma.array(series["wavelength"], mask=series['wavelengthMask'])
    series["wavelengthMask"] = sigma_clip(
        maskedArray, sigma_lower=1, sigma_upper=1, maxiters=3, cenfunc='mean', stdfunc="std")

    return series
