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
from .base_util import base_util

os.environ["TERM"] = "vt100"

# TODO: revisit how the wavelength for each slice is calculated ... take from the continuum, or the central 3-5 pixels?


class horne_extraction(base_util):
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
        

        super(horne_extraction, self).__init__(log, settings, associatedFrame=skySubtractedFrame, dispersionMap=dispersionMap, twoDMapPath=twoDMapPath)

        log.debug("instantiating a new 'horne_extraction' object")
        self.twoDMapPath = twoDMapPath
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
        # self.kw = keyword_lookup(log=self.log, settings=self.settings).get
        # kw = self.kw
        # self.arm = self.skySubtractedFrame.header[kw("SEQ_ARM")]
        # self.dateObs = self.skySubtractedFrame.header[kw("DATE_OBS")]

        # # DETECTOR PARAMETERS LOOKUP OBJECT
        # self.detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

        # GET SKYLINES DATAFRAME
        # self.skylinesDF = get_skylines_dataframe(
        #     self.log, self.settings, self.arm, minBrightnessVIS=1, minBrightnessNIR=0.5
        # )

        # if self.twoDMapPath:
        #     self.mapDF, self.interOrderMaskNDArray = twoD_disp_map_image_to_dataframe(
        #         log=self.log,
        #         slit_length=self.detectorParams["slit_length"],
        #         twoDMapPath=twoDMapPath,
        #         associatedFrame=self.skySubtractedFrame,
        #         kw=self.kw,
        #         dispAxis=self.detectorParams["dispersion-axis"],
        #     )

        #     self.skySubtractedFrame.data[self.interOrderMaskNDArray == 1] = np.nan

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

        # # SET IMAGE ORIENTATION
        # if self.detectorParams["dispersion-axis"] == "x":
        #     self.axisA = "x"
        #     self.axisB = "y"
        # else:
        #     self.axisA = "y"
        #     self.axisB = "x"

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

        # # OPEN AND UNPACK THE 2D IMAGE MAP
        # self.twoDMap = fits.open(twoDMapPath)

        # try:
        #     dpBinx = self.twoDMap[0].header[kw("WIN_BINX")]
        #     dpBiny = self.twoDMap[0].header[kw("WIN_BINY")]
        # except:
        #     dpBinx = 1
        #     dpBiny = 1

        # # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        # binx = 1
        # biny = 1
        # try:
        #     binx = int(self.skySubtractedFrame.header[kw("WIN_BINX")])
        #     biny = int(self.skySubtractedFrame.header[kw("WIN_BINY")])
        # except:
        #     pass

        # ADJUST SLIT HEIGHT IF BINNING IN USE
        if self.binx > 1 or self.biny > 1:
            if self.detectorParams["dispersion-axis"] == "x":
                self.slitHalfLength /= self.binx
            else:
                self.slitHalfLength /= self.biny
            self.slitHalfLength = round(self.slitHalfLength)

        # binxRatio = binx / dpBinx
        # binyRatio = biny / dpBiny

        # xdim = int(self.twoDMap[0].data.shape[1] / binxRatio)
        # ydim = int(self.twoDMap[0].data.shape[0] / binyRatio)
        # xarray = np.tile(np.arange(0, xdim), ydim)
        # yarray = np.repeat(np.arange(0, ydim), xdim)

        # self.skySubtractedFrame.data[self.skySubtractedFrame.data == 0] = np.nan

        # if binxRatio > 1 or binyRatio > 1:
        #     from astropy.nddata import block_reduce

        #     self.twoDMap["WAVELENGTH"].data = block_reduce(
        #         self.twoDMap["WAVELENGTH"].data, (binyRatio, binxRatio), func=np.mean
        #     )
        #     self.twoDMap["SLIT"].data = block_reduce(self.twoDMap["SLIT"].data, (binyRatio, binxRatio), func=np.mean)
        #     self.twoDMap["ORDER"].data = block_reduce(self.twoDMap["ORDER"].data, (binyRatio, binxRatio), func=np.mean)

        # self.imageMap = pd.DataFrame.from_dict(
        #     {
        #         "x": xarray,
        #         "y": yarray,
        #         "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().astype(np.float32),
        #         "slit_position": self.twoDMap["SLIT"].data.flatten().astype(np.float32),
        #         "order": self.twoDMap["ORDER"].data.flatten().astype(np.float32),
        #         "flux": self.skySubtractedFrame.data.flatten().astype(np.float32),
        #     }
        # )
        # self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)

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
        from soxspipe.commonutils.phase3 import write_fits_table_to_disk

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

        # BUILD THE DICTIONARY OF ARRAYS TO ZOOM AND REBIN — KEYS BECOME COLUMN NAMES IN THE ORDER TABLE
        # THIS IS (D-S) IN HORNE 1986, I.E. THE SKY-SUBTRACTED FLUX
        frameArrays = {
            "wavelength": self.twoDMap["WAVELENGTH"].data,
            "fluxRaw": self.skySubtractedFrame.data,
            "slit": self.twoDMap["SLIT"].data,
            "sliceError": self.skySubtractedFrame.uncertainty.array,
        }
        if self.subtractedFrame:
            frameArrays["fluxSky"] = self.subtractedFrame.data

        from soxspipe.commonutils.image_transformer import image_transformer

        # ZOOM AND REBIN ALL ARRAYS ORDER-BY-ORDER, ALSO SIGMA-CLIPPING THE BAD-PIXEL MASK
        transformer = image_transformer(
            log=self.log,
            settings=self.settings,
            orderPixelTable=self.orderPixelTable,
            twoDMapPath=self.twoDMapPath,
            dispersionMap=self.dispersionMap,
            associatedFrame=self.skySubtractedFrame,
            slitHalfLength=self.slitHalfLength,
        )
        transformer.cache_image("fluxRaw", self.skySubtractedFrame.data, associatedMask=self.skySubtractedFrame.mask)
        transformer.cache_image("variance", self.skySubtractedFrame.uncertainty.array ** 2)
        if self.subtractedFrame:
            transformer.cache_image("fluxSky", self.subtractedFrame.data)

        orderSlices = transformer.get_order_slices()
        wlMinMax = transformer.get_order_wavelength_ranges()
        orderRectifiedImages = transformer.get_order_rectified()


        # IF NO SKY FRAME WAS PROVIDED, FILL THE SKY COLUMN WITH ZEROS
        if not self.subtractedFrame:
            for orderTable in orderSlices:
                orderTable["fluxSky"] = [0] * len(orderTable)

        self.log.print(f"\tExtracting {len(orderSlices)} orders\n\n")

        from fundamentals import fmultiprocess

        if self.debug or self.turnOffMP:
            turnOffMP = True
        else:
            turnOffMP = False

        inputArray = [[s,r] for s,r in zip(orderSlices, orderRectifiedImages)]

        extractions = fmultiprocess(
            log=self.log,
            function=extract_single_order,
            inputArray=inputArray,
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
            turnOffMP=True,
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
        extractedOrdersDF = extractedOrdersDF.loc[
            extractedOrdersDF["pixelScaleNm"] < 3
        ]  # FILTER OUT ANY REMAINING BAD PIXELS WITH UNREALISTICALLY LARGE PIXEL SCALE (I.E. WAVELENGTH JUMPS BETWEEN ADJACENT PIXELS)
        if False:
            extractedOrdersDF = self.tune_wavelength_calibration_to_skylines(extractedOrdersDF, arm=arm)
            extractedOrdersDF = self.tune_wavelength_calibration_to_skylines(extractedOrdersDF, arm=arm, byOrder=False)

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

            header[kw("SEQ_ARM")] = arm
            header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
            header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{arm}".upper()

            # DISCRETE ORDERS
            filename = self.filenameTemplate.replace(
                ".fits",
                f"_EXTRACTED_ORDERS{self.noddingSequence}{self.notFlattened}.fits",
            )
            filePath = f"{self.productDir}/{filename}"

            write_fits_table_to_disk(
                log=self.log,
                settings=self.settings,
                header=header,
                tables=[extractedOrdersDF],
                filePath=filePath,
                qc=self.qc,
            )

            utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

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
            write_fits_table_to_disk(
                log=self.log,
                settings=self.settings,
                header=header,
                tables=[mergedTable],
                filePath=filePath,
                qc=self.qc,
            )

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

    def tune_wavelength_calibration_to_skylines(self, extractedOrdersDF, arm, byOrder=True):
        """Skyline-based wavelength calibration tuning helper.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders
        - ``arm`` -- the arm of the instrument
        - ``byOrder`` -- whether to tune calibration order by order, or globally (default: True)

        **Return:**

        - ``extractedOrdersDF`` -- wavelength-corrected dataframe
        """

        self.log.debug("starting the ``tune_wavelength_calibration_to_skylines`` method")

        import numpy as np
        import pandas as pd

        # WAVELENGTH OR B-AXIS (DISPERSION AXIS) COLUMN NAMES
        if True:
            calibrationCol = "wavelengthMean"
            calibrationColSkylines = "wavelength"
            shiftLabel = "wavelength (nm)"
            shiftUnits = "nm"
            tolerance=5
        else:
            calibrationCol = f"{self.axisB}coord"
            shiftLabel = f"{self.axisB} (pixels)"
            calibrationColSkylines = f"fit_{self.axisB}"
            shiftUnits = "px"
            tolerance=15

        uniqueOrders = np.sort(extractedOrdersDF["order"].unique())
        orderShifts = []

        if not byOrder and calibrationCol != "wavelengthMean":
            # CAN'T PIXEL SHIFT ALL ORDERS TOGETHER (ONLY WAVELENGTH SHIFTING MAKES SENSE GLOBALLY)
            return extractedOrdersDF
            

        if not byOrder:
            uniqueOrders = ["all"]

        for o in uniqueOrders:
            if o == "all":
                # DUMMY MASK TO SELECT ALL ORDERS
                mask = extractedOrdersDF["order"] < 100
            else:
                mask = extractedOrdersDF["order"] == o
            orderDF = extractedOrdersDF.loc[mask].sort_values(calibrationCol).copy()
            if orderDF.empty:
                continue

            # PRIME THE SHIFT COLUMN WITH THE ORIGINAL VALUES
            orderDF["wavelength_shifted"] = orderDF["wavelengthMean"]

            for iteration in range(2):
                # EXTRACT NUMERIC ARRAYS FROM ORDER DATAFRAME
                wave, axisBcoord, sky, objectFlux = self._extract_order_arrays(orderDF)
                valid = wave.notna() & sky.notna()
                if not valid.any():
                    continue

                # DETECT PEAKS IN SMOOTHED SKY SPECTRUM
                skyValsOriginal, skyVals, peaks, waveVals, objectVals, axisBVals = self._detect_sky_peaks(wave, sky, objectFlux, axisBcoord, valid)
                wmin, wmax = np.nanmin(waveVals), np.nanmax(waveVals)
                pixelScaleMedian = np.median(orderDF["pixelScaleNm"])

                # MAP CATALOGUE SKYLINES TO PIXEL/WAVELENGTH COORDINATES FOR THIS ORDER
                localSkylines, calibrationSkylines = self._get_local_skylines_for_order(wmin, wmax, o, calibrationColSkylines)

                if calibrationCol == "wavelengthMean":
                    shiftArray = waveVals
                else:
                    shiftArray = axisBVals

                # MATCH DETECTED PEAKS TO ISOLATED CATALOGUE SKYLINES
                matchedSkylinePixels, matchedShifts = self._match_peaks_to_skylines(shiftArray, peaks, calibrationSkylines, tolerance=tolerance)

                # SIGMA-CLIP SHIFT DISTRIBUTION AND COMPUTE MEDIAN
                clippedWave, clippedShifts, medianShift = self._compute_clipped_shift(matchedSkylinePixels, matchedShifts)

                # APPLY MEDIAN SHIFT TO CURRENT ITERATION
                orderDF["wavelength_shifted"] += medianShift

                if iteration == 0:
                    nmShift = medianShift * pixelScaleMedian
                    # print(f"{medianShift:0.2f} {shiftUnits} shift applied to order", o)
                    orderShifts.append(medianShift)
                    if byOrder:
                        self._record_order_shift_qc(o, medianShift)

                # DIAGNOSTIC PLOT SHOWING SKY, SHIFTS, AND OBJECT SPECTRA
                if self.debug:
                    self._plot_skyline_shift_diagnostic(
                        shiftArray, skyValsOriginal, skyVals, peaks, objectVals,
                        localSkylines, calibrationSkylines, matchedSkylinePixels, matchedShifts,
                        clippedWave, clippedShifts, medianShift, o, wmin, wmax,
                        pixelScaleMedian, shiftLabel, iteration, shiftUnits
                    )

        # SHIFT FROM ORDER 1 IS USED AS FALLBACK FOR VIS ORDERS WITH NO SKYLINE MATCHES
        if arm.upper() == "VIS" and byOrder:
            rShift = next((s for o, s in zip(uniqueOrders, orderShifts) if o == 1), 0)
        else:
            rShift = 0

        # APPLY PER-ORDER PIXEL SHIFTS TO WAVELENGTH COLUMN IN FULL DATAFRAME
        for o, s in zip(uniqueOrders, orderShifts):
            if o == "all":
                mask = extractedOrdersDF["order"] < 100
            else:
                mask = extractedOrdersDF["order"] == o
            # PROPAGATE RED-ORDER SHIFT TO VIS ORDERS THAT HAD NO SKYLINE MATCHES
            effectiveShift = rShift if (s == 0 and arm.upper() == "VIS" and o in [2, 3]) else s
            if calibrationCol == "wavelengthMean":
                self.log.print(f"\t\t{effectiveShift:0.2f} {shiftUnits} shift applied to order {o}")
                extractedOrdersDF.loc[mask, "wavelengthMean"] += effectiveShift
            else:
                extractedOrdersDF.loc[mask, "wavelengthMean"] += (
                    effectiveShift * extractedOrdersDF.loc[mask, "pixelScaleNm"]
                )
            

        return extractedOrdersDF

    def _extract_order_arrays(self, orderDF):
        """Extract numeric series for wave, spatial coord, sky, and object flux from an order slice."""
        import pandas as pd
        import numpy as np

        wave = pd.to_numeric(orderDF["wavelength_shifted"], errors="coerce")
        axisBcoord = pd.to_numeric(orderDF[f"{self.axisB}coord"], errors="coerce")
        sky = pd.to_numeric(orderDF["skyFlux"], errors="coerce")
        if "extractedFluxOptimal" in orderDF.columns:
            objectFlux = pd.to_numeric(orderDF["extractedFluxOptimal"], errors="coerce")
        else:
            objectFlux = pd.Series(np.nan, index=orderDF.index)

        return wave, axisBcoord, sky, objectFlux

    def _detect_sky_peaks(self, wave, sky, objectFlux, axisBcoord, valid):
        """Smooth sky spectrum with Savitzky-Golay filter and detect peaks above median."""
        import numpy as np
        from scipy.signal import find_peaks, savgol_filter

        skyValsOriginal = sky[valid].to_numpy()
        skyVals = savgol_filter(sky[valid], window_length=21, polyorder=2)
        peaks, _ = find_peaks(skyVals, height=np.median(skyVals), distance=25)

        waveVals = wave.loc[valid].to_numpy()
        objectVals = objectFlux.loc[valid].to_numpy()
        axisBVals = axisBcoord.loc[valid].to_numpy()

        return skyValsOriginal, skyVals, peaks, waveVals, objectVals, axisBVals

    def _get_local_skylines_for_order(self, wmin, wmax, order, calibrationCol):
        """Fetch catalogue skylines within wavelength range and project to pixel/wavelength coordinates."""
        import pandas as pd
        from soxspipe.commonutils import dispersion_map_to_pixel_arrays

        localSkylinesDF = self.skylinesDF.loc[self.skylinesDF["WAVELENGTH"].between(wmin, wmax)].copy()
        localSkylinesDF.rename(columns={"WAVELENGTH": "wavelength"}, inplace=True)
        localSkylinesDF["order"] = order
        localSkylinesDF["slit_position"] = 0

        if order != "all":
            localSkylinesDF = dispersion_map_to_pixel_arrays(
                log=self.log, dispersionMapPath=self.dispersionMap, orderPixelTable=localSkylinesDF
            )
        localSkylines = pd.to_numeric(localSkylinesDF[calibrationCol], errors="coerce").dropna().to_numpy()

        # ISOLATED FLAG MARKS LINES SUITABLE FOR CENTROIDING
        if "ISOLATED" in self.skylinesDF.columns:
            calibrationSkylinesDF = localSkylinesDF.loc[localSkylinesDF["ISOLATED"] == True]
            calibrationSkylines = pd.to_numeric(calibrationSkylinesDF[calibrationCol], errors="coerce").dropna().to_numpy()
        else:
            calibrationSkylines = []

        return localSkylines, calibrationSkylines

    def _match_peaks_to_skylines(self, shiftArray, peaks, calibrationSkylines, tolerance=15):
        """Match observed sky-peak positions to catalogue skyline positions within a pixel tolerance."""
        import numpy as np

        peakPositions = shiftArray[peaks]
        matchedSkylinePixels = []
        matchedShifts = []

        for ww in calibrationSkylines:
            if len(peakPositions) == 0:
                continue
            idx = np.where(np.abs(peakPositions - ww) < tolerance)[0]
            for i in idx:
                matchedSkylinePixels.append(ww)
                matchedShifts.append(ww - peakPositions[i])

        return matchedSkylinePixels, matchedShifts

    def _compute_clipped_shift(self, matchedSkylinePixels, matchedShifts):
        """Sigma-clip matched shift distribution and return outlier arrays plus median shift."""
        import numpy as np
        from astropy.stats import sigma_clip

        maskedShifts = sigma_clip(matchedShifts, sigma_lower=1.5, sigma_upper=1.5, maxiters=5, cenfunc="mean", stdfunc="std")
        clippedWave = np.asarray(matchedSkylinePixels)[np.asarray(maskedShifts.mask)]
        clippedShifts = np.asarray(maskedShifts.data)[np.asarray(maskedShifts.mask)]
        goodShifts = np.asarray(maskedShifts.data)[~np.asarray(maskedShifts.mask)]

        medianShift = float(np.median(goodShifts)) if len(goodShifts) > 0 else 0.0
        if np.isnan(medianShift):
            medianShift = 0.0

        return clippedWave, clippedShifts, medianShift

    def _record_order_shift_qc(self, order, medianShift):
        """Append a QC entry recording the sky-shift applied to the given order."""
        import pandas as pd
        from datetime import datetime

        utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        self.qc = pd.concat([
            self.qc,
            pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": f"SKY SHIFT O{int(order)}",
                "qc_value": round(medianShift, 3),
                "qc_comment": "Shift applied to wavelength solution based on skyline matches",
                "qc_order": int(order),
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True,
            }).to_frame().T,
        ], ignore_index=True)

    def _plot_skyline_shift_diagnostic(self, shiftArray, skyValsOriginal, skyVals, peaks, objectVals,
                                       localSkylines, calibrationSkylines, matchedSkylinePixels, matchedShifts,
                                       clippedWave, clippedShifts, medianShift, order, wmin, wmax,
                                       pixelScaleMedian, shiftLabel, iteration, shiftUnits):
        """Three-panel diagnostic plot: sky spectrum with skyline markers, shift scatter, object spectrum."""
        import matplotlib
        import matplotlib.pyplot as plt

        matplotlib.use("MacOSX")
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 9), dpi=180, gridspec_kw={"height_ratios": [2, 1, 2]})
        smoothColor = "tab:orange" if iteration == 0 else "tab:green"

        # SKY SPECTRUM WITH DETECTED PEAKS AND CATALOGUE SKYLINE POSITIONS
        ax1.plot(shiftArray, skyValsOriginal, color="black", alpha=0.7, linewidth=0.7, label="skyFlux - original")
        ax1.plot(shiftArray, skyVals, color=smoothColor, linewidth=0.7, label="skyFlux - smoothed")
        for ww in localSkylines:
            ax1.axvline(ww, color="grey", alpha=0.25, linewidth=0.5)
        if len(localSkylines):
            ax1.plot([], [], color="grey", alpha=0.7, linewidth=0.8, label="skyline")
        for ww in calibrationSkylines:
            ax1.axvline(ww, color="blue", alpha=0.5, linewidth=0.6)
        if len(calibrationSkylines):
            ax1.plot([], [], color="blue", alpha=0.7, linewidth=0.8, label="calibration skyline")
        ax1.plot(shiftArray[peaks], skyVals[peaks], "x", color="red", label="skyFlux peaks")
        ax1.set_ylabel("sky flux ($e^{-}$)")
        ax1.legend(loc="best", fontsize=8)

        # MATCHED PEAK SHIFTS WITH SIGMA-CLIPPED OUTLIERS AND MEDIAN INDICATOR
        ax2.axhline(0, color="k", linestyle="--", linewidth=0.7)
        if matchedShifts:
            ax2.scatter(matchedSkylinePixels, matchedShifts, c="tab:green", s=30)
            ax2.scatter(clippedWave, clippedShifts, c="tab:red", s=30, marker="x")
            ax2.axhline(medianShift, color="tab:orange", linestyle=":", linewidth=0.8,
                        label=f"median={medianShift:.3f} {shiftUnits} (catalogue - observed)")
            ax2.legend(loc="best", fontsize=8)
        ax2.set_ylabel(f"shift ({shiftUnits})")
        ax2.set_xlim(ax1.get_xlim())

        # OBJECT SPECTRUM FOR REFERENCE
        ax3.plot(shiftArray, objectVals, color="tab:purple", linewidth=0.2, label="objectFlux")
        ax3.set_xlabel(shiftLabel)
        ax3.set_ylabel("object flux ($e^{-}$)")
        ax3.legend(loc="best", fontsize=8)
        ax3.set_xlim(ax1.get_xlim())

        if order == "all":
            fig.suptitle(f"\nGlobal skyline matching and shifting. WL {wmin:0.2f}-{wmax:0.2f}nm. Mean pixel {pixelScaleMedian:0.2f} nm.")
        else:
            fig.suptitle(f"\nOrder {int(order)} skyline matching and shifting. WL {wmin:0.2f}-{wmax:0.2f}nm. Mean pixel {pixelScaleMedian:0.2f} nm.")
        plt.show()
        plt.close(fig)

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

        utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

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
    inputData,
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

    - ``crossDispersionSlicesDF`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    """
    # log.debug('starting the ``extract_single_order`` method')

    import pandas as pd
    import numpy as np
    from astropy.stats import sigma_clip
    import matplotlib.pyplot as plt

    crossDispersionSlicesDF, orderRectifiedImages = inputData[0], inputData[1]

    # WE ARE BUILDING A SET OF CROSS-SLIT OBJECT PROFILES
    # ALONG THE DISPERSION AXIS
    crossSlitProfiles = []

    # 1) SELECTING THE ORDER FROM THE ORDER PIXEL TABLE - THIS IS THE CONTINUUM OF THE OBJECT
    order = crossDispersionSlicesDF["order"].values[0]

    # SLICE SINGLE IMAGE INTO CROSS-DISPERSION SLICES
    crossDispersionSlicesDF, orderRectifiedImages = generate_masks(crossDispersionSlicesDF=crossDispersionSlicesDF, orderRectifiedImages=orderRectifiedImages)

    # RETURN IF NO SLICES WERE CREATED
    if not len(crossDispersionSlicesDF.index):
        return None

    # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
    # ITERATE FIRST PIXEL IN EACH SLICE AND THEN MOVE TO NEXT
    crossDispersionSlicesDF, orderRectifiedImages = fit_object_profile(
        crossDispersionSlicesDF=crossDispersionSlicesDF,
        orderRectifiedImages=orderRectifiedImages,
        slitHalfLength=slitHalfLength,
        clippingSigma=clippingSigma,
        clippingIterationLimit=clippingIterationLimit,
        hornePolyOrder=hornePolyOrder,
        axisB=axisB,
        order=order,
        debug=debug,
        plt=plt,
    )


    # PLOT THE RECTIFIED IMAGES
    if debug:
        plot_rectified_images(orderRectifiedImages=orderRectifiedImages, order=order)

    extractions = compute_extractions(crossDispersionSlicesDF=crossDispersionSlicesDF, orderRectifiedImages=orderRectifiedImages, order=order)

    return extractions[
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


def compute_extractions(crossDispersionSlicesDF, orderRectifiedImages, order):
    """Compute optimal and boxcar extracted spectra from cross-dispersion slices."""
    import numpy as np

    # CALCULATE HORNE 86 NUMERATOR (EQU 8)
    orderRectifiedImages["horneNumerator"] = np.ma.masked_array(orderRectifiedImages["fluxRaw"] * orderRectifiedImages["objectProfile"] / orderRectifiedImages["variance"], mask=orderRectifiedImages["mask"])
    horneNumeratorSum = orderRectifiedImages["horneNumerator"].sum(axis=1)

    # CALCULATE HORNE 86 DENOMINATOR (EQU 8)
    orderRectifiedImages["horneDenominator"] = np.ma.masked_array(np.power(orderRectifiedImages["objectProfile"],2) / orderRectifiedImages["variance"], mask=orderRectifiedImages["mask"])
    horneDenominatorSum = orderRectifiedImages["horneDenominator"].sum(axis=1)

    orderRectifiedImages["optimalExtraction"] = np.ma.masked_array(orderRectifiedImages["horneNumerator"] / orderRectifiedImages["horneDenominator"], mask=orderRectifiedImages["mask"])

    # plot_rectified_images(orderRectifiedImages=orderRectifiedImages, order=order)

    wavelengthMasked = np.ma.masked_array(orderRectifiedImages["wavelength"], mask=orderRectifiedImages["mask"])
    crossDispersionSlicesDF["wavelengthMean"] = np.ma.mean(wavelengthMasked, axis=1)

    # CALCULATE THE FINAL EXTRACTED SPECTRA
    crossDispersionSlicesDF["varianceSpectrum"] = 1 / horneDenominatorSum
    crossDispersionSlicesDF["extractedFluxOptimal"] = (
        horneNumeratorSum / horneDenominatorSum
    )
    crossDispersionSlicesDF["extractedFluxBoxcar"] = orderRectifiedImages["fluxRaw"].sum(axis=1)
    crossDispersionSlicesDF["skyFlux"] = orderRectifiedImages["fluxSky"].mean(axis=1)
    crossDispersionSlicesDF["extractedFluxBoxcarRobust"] = np.ma.masked_array(orderRectifiedImages["fluxRaw"], mask=orderRectifiedImages["mask"]).sum(axis=1).astype(float)
    crossDispersionSlicesDF["snr"] = crossDispersionSlicesDF["extractedFluxOptimal"] / np.power(
        crossDispersionSlicesDF["varianceSpectrum"], 0.5
    )

    # SORT BY COLUMN NAME
    crossDispersionSlicesDF.sort_values(["wavelengthMean"], ascending=[True], inplace=True)

    # REMOVE 0 WAVELENGTH
    crossDispersionSlicesDF = crossDispersionSlicesDF.loc[crossDispersionSlicesDF["wavelengthMean"] > 0]

    crossDispersionSlicesDF.dropna(
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

    return crossDispersionSlicesDF


def plot_rectified_images(orderRectifiedImages, order):
    """Plot available rectified order image layers for debug inspection."""

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.stats import sigma_clipped_stats

    if not isinstance(orderRectifiedImages, dict) or len(orderRectifiedImages) == 0:
        return None
    
    import matplotlib
    matplotlib.use("MacOSX")

    for key, value in orderRectifiedImages.items():

        
        mean, median, std = sigma_clipped_stats(value, sigma=5.0, stdfunc="std", cenfunc="mean", maxiters=3)
        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True,
        )
        fig.suptitle(f"{key} for Order {order}", fontsize=16)
        if "mask" in key.lower():
            plt.imshow(value.T, interpolation="none", aspect="auto", vmin=0, vmax=1, cmap="viridis")
        else:
            plt.imshow(value.T, interpolation="none", aspect="auto", vmin=mean-2*std, vmax=mean+2*std, cmap="viridis")
        plt.show()

    return None


def generate_masks(crossDispersionSlicesDF, orderRectifiedImages):
    """This function is used to create masks for the cross-dispersion slices and to calculate the pixel scale in wavelength space. 

    **Key Arguments:**

    - ``crossDispersionSlicesDF`` -- the order dataframe
    - ``orderRectifiedImages`` -- the rectified images for the order

    **Return:**

    - ``crossDispersionSlicesDF`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    - ``orderRectifiedImages`` -- the rectified images for the order with updated masks
    """

    import numpy as np
    from astropy.stats import sigma_clip

    # UNPACK RETIFIED ORDER IMAGES
    bpMask = orderRectifiedImages["bpMask"]
    fluxRaw = orderRectifiedImages["fluxRaw"]
    wavelength = orderRectifiedImages["wavelength"]

    # CALCULATE THE PIXEL SCALE BEFORE ANY CLIPPING OCCURS
    crossDispersionSlicesDF["pixelScaleNm"] = np.ma.mean(wavelength, axis=1)
    this = (crossDispersionSlicesDF["pixelScaleNm"].values[2:] - crossDispersionSlicesDF["pixelScaleNm"].values[:-2]) / 2
    this = np.insert(this, 0, np.nan)
    this = np.append(this, np.nan)
    crossDispersionSlicesDF["pixelScaleNm"] = np.abs(this)

    ## REMOVE BAD PIXELS AND COSMIC RAYS FROM THE FLUX ARRAY
    fluxRawMask = _sigma_clip_and_mask(
        fluxRaw, bpMask
    )

    # FAIL SAFE FOR BAD WAVELENGTH VALUES - SOME ODD WAVELENGTHS FROM DISPERSION SOLUTION CAN CAUSE PROBLEMS WITH PROFILE FITTING
    wlMask = wavelength.copy()
    wlMask[wlMask > 0] = 3
    wlMask[wlMask < 0.1] = 1
    wlMask[wlMask > 2] = 0
    wavelengthMasked = np.ma.array(wavelength, mask=wlMask)
    wavelengthMasked = sigma_clip(
        wavelengthMasked,
        sigma_lower=1,
        sigma_upper=1,
        maxiters=3,
        cenfunc="mean",
        stdfunc="std",
        axis=1,
    )

    orderRectifiedImages["mask"] = fluxRawMask | (wlMask.astype(bool))

    # IF THERE IS MORE THAN 1 PIXEL MASKED IN THE CROSS-DISPERSION DIRECTION, THEN FLAG THE ENTIRE COLUMN AS BAD
    fullColumnMaskFlags = np.sum(orderRectifiedImages["mask"], axis=1)
    fullColumnMaskFlags = fullColumnMaskFlags > 1
    orderRectifiedImages["mask"][fullColumnMaskFlags] = True
    crossDispersionSlicesDF["mask"] = [x for x in orderRectifiedImages["mask"]]

    return crossDispersionSlicesDF, orderRectifiedImages


def fit_object_profile(
    crossDispersionSlicesDF,
    orderRectifiedImages,
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

    fluxRawMasked = np.ma.masked_array(orderRectifiedImages["fluxRaw"], mask=orderRectifiedImages["mask"])
    # RETURN THE SUM OF THE ARRAY ELEMENTS OVER THE GIVEN AXIS. MASKED ELEMENTS ARE SET TO 0 INTERNALLY.
    fluxRawMaskedSum = fluxRawMasked.sum(axis=1)    

    ## THIS IS THE NORMALISED FLUX USED FOR FITTING THE DISPERSION PROFILES - THIS IS THE FRACTIONAL FLUX IN HORNE 1986 PAPER
    fluxRawNormalisedMasked = fluxRawMasked / fluxRawMaskedSum[:, np.newaxis]
    
    # DETERMINE LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS
    for slitPixelIndex in range(0, slitHalfLength * 2):

        iteration = 1
        clipped_count = 1

        fractions = fluxRawNormalisedMasked[:, slitPixelIndex]
        dispersionAxisPixels = crossDispersionSlicesDF[f"{axisB}coord"]
        mask = orderRectifiedImages["mask"][:, slitPixelIndex]

        # fractions MAY STILL CONTAIN BAD-PIXEL/CRHs SO DROP PIXELS MASKED IN STEP 1 ABOVE
        a = [fractions, dispersionAxisPixels]
        fractions, dispersionAxisPixels = [np.ma.compressed(np.ma.masked_array(i, mask)) for i in a]

        startCount = len(fractions)
        coeff = []

        while (iteration < clippingIterationLimit) and (clipped_count > 0):
            # FIT A POLYNOMIAL TO THE FRACTIONAL FLUXES
            if len(dispersionAxisPixels):
                coeff = np.polyfit(dispersionAxisPixels, fractions, deg=hornePolyOrder)
            else:
                coeff = []
                break
            residuals = fractions - np.polyval(coeff, dispersionAxisPixels)

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
            a = [fractions, dispersionAxisPixels]
            fractions, dispersionAxisPixels = [np.ma.compressed(np.ma.masked_array(i, masked_residuals.mask)) for i in a]
            clipped_count = startCount - len(fractions)
            percent = (float(clipped_count) / float(startCount)) * 100.0
            # print(f"\tProfile fitting iteration {iteration}, slice index {slitPixelIndex+1}/{slitHalfLength * 2}. {clipped_count} clipped ({percent:0.2f}%) - ORDER {order}")
            iteration = iteration + 1

        # GENERATE THE FINAL FITTING PROFILE FOR THIS SLIT POSITION
        if len(coeff):
            profile = np.polyval(coeff, crossDispersionSlicesDF[f"{axisB}coord"])
            profile[profile < 0] = 0
        else:
            profile = np.zeros_like(crossDispersionSlicesDF[f"{axisB}coord"])
        crossSlitProfiles.append(profile)

        if debug:
            plt.scatter(dispersionAxisPixels, fractions, alpha=0.2)
            plt.plot(crossDispersionSlicesDF[f"{axisB}coord"], profile, color="red")
            plt.title(f"Fitted Profile for Order {order}")
            plt.ylim([-1, 1])
            plt.show()

    crossSlitProfiles = np.array(crossSlitProfiles)
    transposedProfiles = crossSlitProfiles.T.tolist()
    crossDispersionProfile = np.array([np.array(t) for t in transposedProfiles])

    crossDispersionProfileSums = np.array([x.sum() for x in crossDispersionProfile])
    orderRectifiedImages["objectProfile"] = (
        crossDispersionProfile / crossDispersionProfileSums[:, np.newaxis]
    )
    crossDispersionSlicesDF["objectProfile"] = [x for x in orderRectifiedImages["objectProfile"]]

    return crossDispersionSlicesDF, orderRectifiedImages


def _sigma_clip_and_mask(fluxRaw, bpMask):
    import numpy as np
    from astropy.stats import sigma_clip

    # BUILDS A 2D DISTANCE-FROM-CENTER WEIGHTING ARRAY (SAME SHAPE AS fluxRaw) .. ONLY USED FOR SIGMA-CLIPPING THE FLUX ARRAY, NOT FOR WEIGHTING THE PROFILE FITTING
    # EACH ROW GETS IDENTICAL CENTER-DISTANCE WEIGHTS, LATER USED TO BIAS SIGMA-CLIPPING SO OUTLIERS FARTHER FROM THE CENTER ARE PENALIZED MORE.
    sliceIndexArray = np.tile(
        np.abs(np.arange(fluxRaw.shape[1]) - fluxRaw.shape[1] // 2) * 50 + 1,
        (fluxRaw.shape[0], 1),
    )
    bpMask[bpMask > 1] = 1
    fluxRawMasked = np.ma.array(fluxRaw, mask=np.isnan(fluxRaw) | bpMask)
    fluxRawMasked = sigma_clip(
        fluxRawMasked * sliceIndexArray,
        sigma_lower=3000,
        sigma_upper=7,
        maxiters=3,
        cenfunc="mean",
        stdfunc="std",
        axis=1,
    )
    
    ## COMPUTE NUMBER OF NEWLY MASKED PIXELS IN EACH ROW AFTER SIGMA-CLIPPING
    newlyMaskedPixels = np.zeros(fluxRaw.shape[0])
    for i in range(fluxRaw.shape[0]):
        newlyMaskedPixels[i] = np.sum(fluxRawMasked.mask[i]) - np.sum(bpMask[i])
        # print(f"Row {i}: {np.sum(fluxRawMasked.mask[i])} pixels masked, {np.sum(bpMask[i])} pixels were already masked, {newlyMaskedPixels[i]} newly masked")
    #print(f"Total newly masked pixels across all slices: {int(np.sum(newlyMaskedPixels))}")

    return fluxRawMasked.mask
