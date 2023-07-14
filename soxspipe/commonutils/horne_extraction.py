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

# TODO: pass in error frame and use instead of recalculating errors later on
# TODO: replace RON value with true value from FITS header
# TODO: move extract code to extract_order and have extract be a loop over all orders
# TODO: remove sliceFluxNormalisedSum -- this was a debug variable
# TODO: include the BPM in create_cross_dispersion_slice (at least)


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
        from astropy.nddata import CCDData
        from astropy import units as u

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

    def extract(self, order):
        """
        *extract the object spectrum*

        **Return:**
            - ``extracted_wave_spectrum`` -- pixel wavelength array
            - ``extracted_spectrum`` -- pixel flux values
            - ``extracted_spectrum_nonopt`` -- the ~boxcar extraction (follows order curvature). Pixel flux values.

        **Usage:**

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed
        ```

        ```python
        usage code
        ```
        """
        self.log.debug('starting the ``extract`` method')

        import yaml
        import pandas as pd
        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt

        # TODO: these can likely be removed once dataframes added
        extracted_spectrum = []
        extracted_spectrum_nonopt = []
        extracted_wave_spectrum = []
        total_weights = []

        # WE ARE BUILDING A SET OF CROSS-SLIT OBJECT PROFILES
        # ALONG THE DISPERSION AXIS
        crossSlitProfiles = []

        # GET UNIQUE VALUES IN COLUMN
        uniquecolNames = self.orderPixelTable['order'].unique()

        # 1) SELECTING THE ORDER FROM THE ORDER PIXEL TABLE - THIS IS THE CONTINUUM OF THE OBJECT
        crossDispersionSlices = self.orderPixelTable.loc[self.orderPixelTable['order'] == order]

        # CREATE THE SLICES AND DROP SLICES WITH ALL NANs (TYPICALLY PIXELS WITH NANs IN 2D IMAGE MAP)
        print(f"\n# SLICING ORDER INTO CROSS-DISPERSION SLICES")
        crossDispersionSlices = crossDispersionSlices.apply(lambda x: self.create_cross_dispersion_slice(x), axis=1)
        crossDispersionSlices.dropna(axis='index', how='any', subset=["sliceRawFlux"], inplace=True)

        # VERTICALLY STACK THE SLICES INTO PSEDUO-RECTIFIED IMAGE
        fluxNormalisedImage = np.vstack(crossDispersionSlices["sliceFluxNormalised"])
        weightImage = np.vstack(crossDispersionSlices["sliceWeights"])
        maskImage = np.vstack(crossDispersionSlices["sliceMask"])
        errorImage = np.vstack(crossDispersionSlices["sliceError"])

        # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
        # ITERATE FIRST PIXEL IN EACH SLICE AND THEN MOVE TO NEXT
        print(f"\n# FITTING CROSS-SLIT FLUX NORMALISED PROFILES")
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
                print(f"\tProfile fitting iteration {iteration}, slice index {slitPixelIndex+1}/{self.slitHalfLength * 2}. {clipped_count} clipped ({percent:0.2f}%).")
                iteration = iteration + 1
                if iteration > 1:
                    sys.stdout.write("\x1b[1A\x1b[2K")

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

        def extract_spectrum(
                series):
            """This function is used to optimal extract the spectrum of the object at each slice location (via Horne method)
            """
            gain = 1.0

            series["sliceFittedProfileNormalised"] = np.array(series["sliceFittedProfile"]) / np.array(series["sliceFittedProfile"]).sum()
            series["sliceVariance"] = self.ron + np.abs(series["sliceRawFlux"] * series["sliceFittedProfileNormalised"] + series["sliceSky"]) / gain
            sliceRejection = np.power((series["sliceRawFlux"] - series["sliceRawFluxMaskedSum"] * series["sliceFittedProfileNormalised"]), 2) / series["sliceVariance"]
            series["sliceVarianceRejectNumber"] = np.power((series["sliceRawFlux"] - series["sliceRawFluxMaskedSum"] * series["sliceFittedProfileNormalised"]), 2) / series["sliceVariance"]
            series["sliceVarianceRejectLimit"] = self.globalClippingSigma
            mask = np.zeros_like(sliceRejection)
            mask[sliceRejection > self.globalClippingSigma] = 1
            series["sliceFittedProfileNormalisedGood"] = np.ma.masked_array(series["sliceFittedProfileNormalised"], mask)
            series["sliceRawFluxGood"] = np.ma.masked_array(series["sliceRawFlux"], mask)
            series["sliceVarianceGood"] = np.ma.masked_array(series["sliceVariance"], mask)
            series["horneNumerator"] = series["sliceRawFluxGood"] * series["sliceFittedProfileNormalisedGood"] / series["sliceVarianceGood"]
            series["horneNumeratorSum"] = series["horneNumerator"].sum()
            series["horneDenominator"] = np.power(series["sliceFittedProfileNormalisedGood"], 2) / series["sliceVarianceGood"]
            series["horneDenominatorSum"] = series["horneDenominator"].sum()
            series["wavelengthMedian"] = np.median(series["wavelength"])
            series["fudged"] = False

            # TODO: if all above sliceVarianceRejectLimit  ... what?
            if series["horneDenominatorSum"] == 0:
                series["horneNumerator"] = series["sliceRawFlux"] * series["sliceFittedProfileNormalised"] / series["sliceVariance"]
                series["horneNumeratorSum"] = series["horneNumerator"].sum()
                series["horneDenominator"] = np.power(series["sliceFittedProfileNormalised"], 2) / series["sliceVariance"]
                series["horneDenominatorSum"] = series["horneDenominator"].sum()
                series["fudged"] = True

            series["extractedFluxOptimal"] = series["horneNumeratorSum"] / series["horneDenominatorSum"]
            series["extractedFluxBoxcar"] = series["sliceRawFlux"].sum()

            return series

        print(f"\n# EXTRACTING THE SPECTRUM")
        crossDispersionSlices = crossDispersionSlices.apply(extract_spectrum, axis=1)

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
            crossDispersionSlices.to_sql('my_export_table', con=conn,
                                         index=False, if_exists='replace')

        extracted_wave_spectrum = crossDispersionSlices["wavelengthMedian"]
        extracted_spectrum = crossDispersionSlices["extractedFluxOptimal"]
        extracted_spectrum_nonopt = crossDispersionSlices["extractedFluxBoxcar"]

        self.log.debug('completed the ``extract`` method')

        return (extracted_wave_spectrum, extracted_spectrum, extracted_spectrum_nonopt)

    def extract_order(
            self,
            order):
        """*extract object spectrum from a single order*

        **Key Arguments:**
            # -

        **Return:**
            - ``order`` -- 

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
        self.log.debug('starting the ``extract_order`` method')

        self.log.debug('completed the ``extract_order`` method')
        return order

    def create_cross_dispersion_slice(
            self,
            series):
        """This function is used to create a single, 1-pixel wide cross-dispersion slice of object data. When applied to the dataframe, a single slice is created for each discrete pixel position in the dispersion direction
        """

        import numpy as np
        from astropy.stats import sigma_clip

        x = series['xcoord_centre']
        y = series['ycoord']
        series["wavelength"] = self.twoDMap["WAVELENGTH"].data[int(y), int(x) - self.slitHalfLength:int(x) + self.slitHalfLength]
        series["sliceRawFlux"] = self.skySubtractedFrame.data[int(y), int(x) - self.slitHalfLength:int(x) + self.slitHalfLength]
        series["sliceSky"] = self.skyModelFrame.data[int(y), int(x) - self.slitHalfLength:int(x) + self.slitHalfLength]
        series["sliceError"] = self.skySubtractedFrame.uncertainty[int(y), int(x) - self.slitHalfLength:int(x) + self.slitHalfLength]

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

    # use the tab-trigger below for new method
    # xt-class-method
