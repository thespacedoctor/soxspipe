#!/usr/bin/env python
# encoding: utf-8
"""
*perform optimal source extraction using the Horne method (Horne 1986)*

:Author:
    Marco Landoni & David Young

:Date Created:
    May 17, 2023
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


class horne_extraction(object):
    """
    *perform optimal source extraction using the Horne method (Horne 1986)*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``skyModelFrame`` -- sky model frame
        - ``skySubtractedFrame`` -- sky subtracted frame
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
        - ``twoDMapPath`` -- 2D dispersion map image path
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``dispersionMap`` -- the FITS binary table containing dispersion map polynomial
        - ``sofName`` -- the set-of-files filename

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a horne_extraction object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``horne_extraction`` to documentation
        - create a blog post about what ``horne_extraction`` does
    ```

    ```python
    usage code
    ```

    """

    def __init__(
            self,
            log,
            settings,
            skyModelFrame,
            skySubtractedFrame,
            twoDMapPath,
            recipeName=False,
            qcTable=False,
            productsTable=False,
            dispersionMap=False,
            sofName=False

    ):
        import numpy as np
        import pandas as pd
        from astropy.io import fits
        from soxspipe.commonutils import keyword_lookup
        from astropy.nddata import CCDData
        from astropy import units as u
        from os.path import expanduser

        self.log = log
        log.debug("instansiating a new 'horne_extraction' object")
        self.dispersionMap = dispersionMap
        self.twoDMapPath = twoDMapPath
        self.settings = settings
        self.productsTable = productsTable
        self.qcTable = qcTable
        self.recipeName = recipeName
        self.sofName = sofName

        # COLLECT SETTINGS FROM SETTINGS FILE
        self.slitHalfLength = int(self.settings["soxs-stare"]["horne-extraction-slit-length"] / 2)
        self.clippingSigma = self.settings["soxs-stare"]["horne-extraction-profile-clipping-sigma"]
        self.clippingIterationLimit = self.settings["soxs-stare"]["horne-extraction-profile-clipping-iteration-count"]
        self.globalClippingSigma = self.settings["soxs-stare"]["horne-extraction-profile-global-clipping-sigma"]

        # TODO: replace this value with true value from FITS header
        self.ron = 3.0

        # OPEN THE SKY-SUBTRACTED FRAME
        self.skySubtractedFrame = CCDData.read(skySubtractedFrame, hdu=0, unit=u.electron,
                                               hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS',
                                               key_uncertainty_type='UTYPE')

        # OPEN THE SKY-MODEL FRAME
        self.skyModelFrame = CCDData.read(skyModelFrame, hdu=0, unit=u.electron,
                                          hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS',
                                          key_uncertainty_type='UTYPE')

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = self.skySubtractedFrame.header[kw("SEQ_ARM")]
        self.dateObs = self.skySubtractedFrame.header[kw("DATE_OBS")]

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
        xdim = self.twoDMap[0].data.shape[1]
        ydim = self.twoDMap[0].data.shape[0]
        xarray = np.tile(np.arange(0, xdim), ydim)
        yarray = np.repeat(np.arange(0, ydim), xdim)
        self.imageMap = pd.DataFrame.from_dict({
            "x": xarray,
            "y": yarray,
            "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().byteswap().newbyteorder(),
            "slit_position": self.twoDMap["SLIT"].data.flatten().byteswap().newbyteorder(),
            "order": self.twoDMap["ORDER"].data.flatten().byteswap().newbyteorder(),
            "flux": self.skySubtractedFrame.data.flatten().byteswap().newbyteorder()
        })
        self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)

        # FIND THE OBJECT TRACE IN EACH ORDER
        from soxspipe.commonutils import detect_continuum
        detector = detect_continuum(
            log=self.log,
            pinholeFlat=self.skySubtractedFrame,
            dispersion_map=self.dispersionMap,
            settings=self.settings,
            sofName=self.sofName,
            recipeName=self.recipeName,
            qcTable=self.qcTable,
            productsTable=self.productsTable
        )
        productPath, self.qcTable, self.productsTable = detector.get()
        from soxspipe.commonutils.toolkit import unpack_order_table
        # UNPACK THE ORDER TABLE
        orderPolyTable, self.orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=productPath)

    def extract(
            self):
        """*extract the full spectrum order-by-order and return FITS Binary table containing order-merged spectrum*

        **Return:**
            - None
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

        kw = self.kw
        arm = self.arm

        # GET UNIQUE VALUES IN COLUMN
        uniqueOrders = self.orderPixelTable['order'].unique()
        extractions = []

        for order in uniqueOrders:
            crossDispersionSlices = self.extract_single_order(order)
            # from tabulate import tabulate
            # print(tabulate(crossDispersionSlices.tail(10), headers='keys', tablefmt='psql'))
            extractions.append(crossDispersionSlices)

        fig = plt.figure(figsize=(16, 2), constrained_layout=True, dpi=320)
        gs = fig.add_gridspec(1, 1)
        toprow = fig.add_subplot(gs[0:1, :])
        addedLegend = True
        for df, o in zip(extractions, uniqueOrders):
            extracted_wave_spectrum = df["wavelengthMedian"]
            extracted_spectrum = df["extractedFluxOptimal"]
            extracted_spectrum_nonopt = df["extractedFluxBoxcarRobust"]

            try:
                if addedLegend:
                    label = "Robust Boxcar Extraction"
                else:
                    label = None
                toprow.plot(extracted_wave_spectrum[10:-10], extracted_spectrum_nonopt[10:-10], color="gray", alpha=0.8, zorder=1, label=label)
                # plt.plot(extracted_wave_spectrum, extracted_spectrum_nonopt, color="gray", alpha=0.1, zorder=1)
                line = toprow.plot(extracted_wave_spectrum[10:-10], extracted_spectrum[10:-10], zorder=2)

                toprow.text(extracted_wave_spectrum[10:-10].mean(), extracted_spectrum[10:-10].mean() + 5 * extracted_spectrum[10:-10].std(), int(o), fontsize=10, c=line[0].get_color(), verticalalignment='bottom')
                addedLegend = False
            except:
                self.log.warning(f"Order skipped: {o}")

        toprow.legend(loc='lower right', bbox_to_anchor=(1, -0.5),
                      fontsize=8)
        toprow.set_ylabel('flux ($e^{-}$)', fontsize=10)
        toprow.set_xlabel(f'wavelength (nm)', fontsize=10)
        toprow.set_title(
            f"Optimally Extracted Object Spectrum ({arm.upper()})", fontsize=11)

        filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_ORDERS_QC_PLOT.pdf")
        filePath = f"{self.qcDir}/{filename}"
        # plt.tight_layout()
        # plt.show()
        plt.savefig(filePath, dpi='figure')
        plt.close()
        plt.show()

        # MERGE THE ORDER SPECTRA
        extractedOrdersDF = pd.concat(extractions, ignore_index=True)
        mergedSpectum = self.merge_extracted_orders(extractedOrdersDF)

        # CONVERT TO FITS BINARY TABLE
        header = copy.deepcopy(self.skySubtractedFrame.header)
        header.pop(kw("DPR_CATG"))
        header.pop(kw("DPR_TYPE"))
        with suppress(KeyError):
            header.pop(kw("DET_READ_SPEED"))
        with suppress(KeyError):
            header.pop(kw("CONAD"))
        with suppress(KeyError):
            header.pop(kw("GAIN"))
        with suppress(KeyError):
            header.pop(kw("RON"))

        header["HIERARCH " + kw("PRO_TECH")] = header.pop(kw("DPR_TECH"))
        extractedOrdersTable = Table.from_pandas(extractedOrdersDF)
        BinTableHDU = fits.table_to_hdu(extractedOrdersTable)

        header[kw("SEQ_ARM")] = arm
        header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
        header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{arm}".upper()
        priHDU = fits.PrimaryHDU(header=header)

        hduList = fits.HDUList([priHDU, BinTableHDU])

        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"

        filePath = f"{outDir}/{self.filenameTemplate}".replace(".fits", "_EXTRACTED_ORDERS.fits")
        hduList.writeto(filePath, checksum=True, overwrite=True)

        self.merge_extracted_orders(extractedOrdersDF)

        self.log.debug('completed the ``extract`` method')
        return None

    def extract_single_order(self, order):
        """
        *extract the object spectrum for a single order*

        **Return:**
            - ``crossDispersionSlices`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
        """
        self.log.debug('starting the ``extract_single_order`` method')

        import yaml
        import pandas as pd
        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt

        # WE ARE BUILDING A SET OF CROSS-SLIT OBJECT PROFILES
        # ALONG THE DISPERSION AXIS
        crossSlitProfiles = []

        # 1) SELECTING THE ORDER FROM THE ORDER PIXEL TABLE - THIS IS THE CONTINUUM OF THE OBJECT
        crossDispersionSlices = self.orderPixelTable.loc[self.orderPixelTable['order'] == order]

        # ADD SOME DATA TO THE SLICES
        xstart = crossDispersionSlices["xcoord_centre"].astype(int) - self.slitHalfLength
        xstop = crossDispersionSlices["xcoord_centre"].astype(int) + self.slitHalfLength
        ycoord = crossDispersionSlices["ycoord"].astype(int)
        xcoords = list(map(lambda x: list(range(x[0], x[1])), zip(xstart, xstop)))
        ycoords = list(map(lambda x: [x] * self.slitHalfLength * 2, ycoord))
        crossDispersionSlices["wavelength"] = list(self.twoDMap["WAVELENGTH"].data[ycoords, xcoords])
        crossDispersionSlices["sliceRawFlux"] = list(self.skySubtractedFrame.data[ycoords, xcoords])
        crossDispersionSlices["sliceSky"] = list(self.skyModelFrame.data[ycoords, xcoords])
        crossDispersionSlices["sliceError"] = list(self.skySubtractedFrame.uncertainty[ycoords, xcoords])

        # CREATE THE SLICES AND DROP SLICES WITH ALL NANs (TYPICALLY PIXELS WITH NANs IN 2D IMAGE MAP)
        print(f"\n# SLICING ORDER INTO CROSS-DISPERSION SLICES - ORDER {order}")
        crossDispersionSlices = crossDispersionSlices.apply(lambda x: self.create_cross_dispersion_slice(x), axis=1)
        crossDispersionSlices.dropna(axis='index', how='any', subset=["sliceRawFlux"], inplace=True)

        # VERTICALLY STACK THE SLICES INTO PSEDUO-RECTIFIED IMAGE
        fluxNormalisedImage = np.vstack(crossDispersionSlices["sliceFluxNormalised"])
        weightImage = np.vstack(crossDispersionSlices["sliceWeights"])
        maskImage = np.vstack(crossDispersionSlices["sliceMask"])
        errorImage = np.vstack(crossDispersionSlices["sliceError"])

        # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
        # ITERATE FIRST PIXEL IN EACH SLICE AND THEN MOVE TO NEXT
        print(f"\n# FITTING CROSS-SLIT FLUX NORMALISED PROFILES - ORDER {order}")
        for slitPixelIndex in range(0, self.slitHalfLength * 2):

            iteration = 1
            clipped_count = 1

            fractions = fluxNormalisedImage[:, slitPixelIndex]
            wave_px = crossDispersionSlices["ycoord"]
            weights = weightImage[:, slitPixelIndex]
            mask = maskImage[:, slitPixelIndex]

            # fractions WILL STILL CONTAIN BAD-PIXEL/CRHs SO EXCLUDE PIXELS MASKED IN STEP 1 ABOVE
            a = [fractions, wave_px, weights]
            fractions, wave_px, weights = [np.ma.compressed(np.ma.masked_array(
                i, mask)) for i in a]

            startCount = len(fractions)
            while (iteration < self.clippingIterationLimit) and (clipped_count > 0):

                coeff = np.polyfit(wave_px, fractions, deg=2)
                residuals = fractions - np.polyval(coeff, wave_px)

                # REMOVE REMAINING OUTLIERS
                masked_residuals = sigma_clip(residuals, sigma_lower=self.clippingSigma, sigma_upper=self.clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
                # REDUCE ARRAYS TO NON-MASKED VALUES
                a = [fractions, wave_px, weights]
                fractions, wave_px, weights = [np.ma.compressed(np.ma.masked_array(
                    i, masked_residuals.mask)) for i in a]
                clipped_count = startCount - len(fractions)
                percent = (float(clipped_count) / float(startCount)) * 100.
                print(f"\tProfile fitting iteration {iteration}, slice index {slitPixelIndex+1}/{self.slitHalfLength * 2}. {clipped_count} clipped ({percent:0.2f}%) - ORDER {order}")
                iteration = iteration + 1
                # if iteration > 1:
                #     sys.stdout.write("\x1b[1A\x1b[2K")

            # GENERATE THE FINAL FITTING PROFILE FOR THIS SLIT POSITION
            profile = np.polyval(coeff, crossDispersionSlices["ycoord"])
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
        crossDispersionSlices["sliceFittedProfile"] = transposedProfiles

        print(f"\n# EXTRACTING THE SPECTRUM - ORDER {order}")
        crossDispersionSlices = crossDispersionSlices.apply(lambda x: self.extract_spectrum(x), axis=1)

        if True:
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

            conn = sql.connect("/Users/Dave/Desktop/pandas_export.db")
            # SEND TO DATABASE
            crossDispersionSlices.to_sql(f'order_{order}', con=conn,
                                         index=False, if_exists='replace')

        self.log.debug('completed the ``extract_single_order`` method')

        return crossDispersionSlices[['order', 'xcoord_centre', 'ycoord', 'wavelengthMedian', 'extractedFluxOptimal', 'extractedFluxBoxcar', 'extractedFluxBoxcarRobust']]

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

        from tabulate import tabulate
        print(tabulate(extractedOrdersDF.head(10), headers='keys', tablefmt='psql'))

        self.log.debug('completed the ``merge_extracted_orders`` method')
        return None

    def create_cross_dispersion_slice(
            self,
            series):
        """This function is used to create a single, 1-pixel wide cross-dispersion slice of object data. When applied to the dataframe, a single slice is created for each discrete pixel position in the dispersion direction
        """

        import numpy as np
        from astropy.stats import sigma_clip

        # SET SLICES WITH ALL NAN FLUX TO A SINGLE NAN - WILL BE REMOVED FROM THE DATAFRAME LATER
        if np.isnan(series["sliceRawFlux"]).all():
            series["sliceRawFlux"] = np.nan
            return series

        # SIGMA-CLIP THE DATA TO REMOVE COSMIC/BAD-PIXELS
        sliceMasked = sigma_clip(
            series["sliceRawFlux"], sigma_lower=7, sigma_upper=7, maxiters=1, cenfunc='median', stdfunc="mad_std")
        series["sliceMask"] = sliceMasked.mask
        series["sliceRawFluxMaskedSum"] = sliceMasked.sum()

        # WEIGHTS ARE NOT YET USED
        series["sliceWeights"] = self.ron + np.abs(series["sliceRawFlux"]) / (sliceMasked.sum() * sliceMasked.sum())

        # NORMALISE THE FLUX
        series["sliceFluxNormalised"] = series["sliceRawFlux"] / sliceMasked.sum()
        series["sliceFluxNormalisedSum"] = series["sliceFluxNormalised"].sum().item()
        return series

    def extract_spectrum(
            self,
            series):
        """This function is used to optimal extract the spectrum of the object at each slice location (via Horne method)
        """
        import numpy as np

        gain = 1.0

        series["sliceFittedProfileNormalised"] = np.array(series["sliceFittedProfile"]) / np.array(series["sliceFittedProfile"]).sum()
        series["sliceVariance"] = self.ron + np.abs(series["sliceRawFlux"] * series["sliceFittedProfileNormalised"] + series["sliceSky"]) / gain
        sliceRejection = np.power((series["sliceRawFlux"] - series["sliceRawFluxMaskedSum"] * series["sliceFittedProfileNormalised"]), 2) / series["sliceVariance"]
        series["sliceVarianceRejectNumber"] = np.power((series["sliceRawFlux"] - series["sliceRawFluxMaskedSum"] * series["sliceFittedProfileNormalised"]), 2) / series["sliceVariance"]
        series["sliceVarianceRejectLimit"] = self.globalClippingSigma
        mask = np.zeros_like(sliceRejection)
        mask[sliceRejection > self.globalClippingSigma] = 1
        series["sliceFittedProfileNormalisedGood"] = np.ma.masked_array(series["sliceFittedProfileNormalised"], mask)
        series["sliceFittedProfileNormalisedGoodSum"] = series["sliceFittedProfileNormalisedGood"].sum()
        series["sliceRawFluxGood"] = np.ma.masked_array(series["sliceRawFlux"], mask)
        series["sliceVarianceGood"] = np.ma.masked_array(series["sliceVariance"], mask)
        series["horneNumerator"] = series["sliceRawFluxGood"] * series["sliceFittedProfileNormalisedGood"] / series["sliceVarianceGood"]
        series["horneNumeratorSum"] = series["horneNumerator"].sum()
        series["horneDenominator"] = np.power(series["sliceFittedProfileNormalisedGood"], 2) / series["sliceVarianceGood"]
        series["horneDenominatorSum"] = series["horneDenominator"].sum()
        series["wavelengthMedian"] = np.median(series["wavelength"])
        series["fudged"] = False

        # TODO: if all above sliceVarianceRejectLimit  ... what?
        if series["horneDenominatorSum"] == 0 or not isinstance(series["horneDenominatorSum"], float):
            series["horneNumerator"] = series["sliceRawFlux"] * series["sliceFittedProfileNormalised"] / series["sliceVariance"]
            series["horneNumeratorSum"] = series["horneNumerator"].sum()
            series["horneDenominator"] = np.power(series["sliceFittedProfileNormalised"], 2) / series["sliceVariance"]
            series["horneDenominatorSum"] = series["horneDenominator"].sum()
            series["fudged"] = True

        series["extractedFluxOptimal"] = series["horneNumeratorSum"] / series["horneDenominatorSum"]
        series["extractedFluxBoxcar"] = series["sliceRawFlux"].sum()
        series["extractedFluxBoxcarRobust"] = series["sliceRawFluxMaskedSum"]

        return series

    # use the tab-trigger below for new method
    # xt-class-method
