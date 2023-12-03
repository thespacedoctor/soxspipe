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
# TODO: replace 'gain' with true gain from fits headers

class horne_extraction(object):
    """
    *perform optimal source extraction using the Horne method (Horne 1986)*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``skyModelFrame`` -- path to sky model frame
        - ``skySubtractedFrame`` -- path to sky subtracted frame
        - ``twoDMapPath`` -- path to 2D dispersion map image path
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
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
        from os.path import expanduser
        from soxspipe.commonutils import keyword_lookup
        from soxspipe.commonutils import detect_continuum
        from soxspipe.commonutils.toolkit import unpack_order_table

        self.log = log
        log.debug("instansiating a new 'horne_extraction' object")
        self.dispersionMap = dispersionMap
        self.twoDMapPath = twoDMapPath
        self.settings = settings
        self.products = productsTable
        self.qc = qcTable
        self.recipeName = recipeName
        self.sofName = sofName

        home = expanduser("~")
        self.outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"

        # COLLECT SETTINGS FROM SETTINGS FILE
        self.slitHalfLength = int(self.settings["soxs-stare"]["horne-extraction-slit-length"] / 2)
        self.clippingSigma = self.settings["soxs-stare"]["horne-extraction-profile-clipping-sigma"]
        self.clippingIterationLimit = self.settings["soxs-stare"]["horne-extraction-profile-clipping-iteration-count"]
        self.globalClippingSigma = self.settings["soxs-stare"]["horne-extraction-profile-global-clipping-sigma"]

        # TODO: replace this value with true value from FITS header
        self.ron = 3.0

        # OPEN THE SKY-SUBTRACTED FRAME
        if isinstance(skySubtractedFrame, CCDData):
            self.skySubtractedFrame = skySubtractedFrame
        else:
            self.skySubtractedFrame = CCDData.read(skySubtractedFrame, hdu=0, unit=u.electron,
                                                   hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS',
                                                   key_uncertainty_type='UTYPE')

        # OPEN THE SKY-MODEL FRAME
        if isinstance(skyModelFrame, CCDData):
            self.skyModelFrame = skyModelFrame
        else:
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

        self.skySubtractedFrame.data[self.skySubtractedFrame.data == 0] = np.nan

        try:
            self.imageMap = pd.DataFrame.from_dict({
                "x": xarray,
                "y": yarray,
                "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().byteswap().newbyteorder(),
                "slit_position": self.twoDMap["SLIT"].data.flatten().byteswap().newbyteorder(),
                "order": self.twoDMap["ORDER"].data.flatten().byteswap().newbyteorder(),
                "flux": self.skySubtractedFrame.data.flatten()
            })
            self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)
        except:
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
        detector = detect_continuum(
            log=self.log,
            pinholeFlat=self.skySubtractedFrame,
            dispersion_map=self.dispersionMap,
            settings=self.settings,
            sofName=self.sofName,
            recipeName=self.recipeName,
            qcTable=self.qc,
            productsTable=self.products
        )
        productPath, self.qc, self.products = detector.get()

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
        from datetime import datetime

        kw = self.kw
        arm = self.arm

        # GET UNIQUE VALUES IN COLUMN
        uniqueOrders = self.orderPixelTable['order'].unique()
        # uniqueOrders = uniqueOrders[3:5]
        extractions = []

        self.log.print("\n# PERFORMING OPTIMAL SOURCE EXTRACTION (Horne Method)\n\n")

        # ADD SOME DATA TO THE SLICES
        orderSlices = []
        for order in uniqueOrders:
            orderTable = self.orderPixelTable.loc[self.orderPixelTable['order'] == order]
            xstart = orderTable["xcoord_centre"].astype(int) - self.slitHalfLength
            xstop = orderTable["xcoord_centre"].astype(int) + self.slitHalfLength
            ycoord = orderTable["ycoord"].astype(int)
            xcoords = list(map(lambda x: list(range(x[0], x[1])), zip(xstart, xstop)))
            ycoords = list(map(lambda x: [x] * self.slitHalfLength * 2, ycoord))
            orderTable["wavelength"] = list(self.twoDMap["WAVELENGTH"].data[ycoords, xcoords])
            orderTable["sliceRawFlux"] = list(self.skySubtractedFrame.data[ycoords, xcoords])
            orderTable["sliceSky"] = list(self.skyModelFrame.data[ycoords, xcoords])
            orderTable["sliceError"] = list(self.skySubtractedFrame.uncertainty[ycoords, xcoords])
            orderSlices.append(orderTable)

        from fundamentals import fmultiprocess
        extractions = fmultiprocess(log=self.log, function=extract_single_order,
                                    inputArray=orderSlices, poolSize=False, timeout=300, ron=self.ron, slitHalfLength=self.slitHalfLength, clippingSigma=self.clippingSigma, clippingIterationLimit=self.clippingIterationLimit, globalClippingSigma=self.globalClippingSigma)

        fig = plt.figure(figsize=(16, 2), constrained_layout=True, dpi=320)
        gs = fig.add_gridspec(1, 1)
        toprow = fig.add_subplot(gs[0:1, :])
        addedLegend = True
        for df, o in zip(extractions, uniqueOrders):
            extracted_wave_spectrum = df["wavelengthMedian"]
            extracted_spectrum = df["extractedFluxOptimal"]
            extracted_spectrum_nonopt = df["extractedFluxBoxcarRobust"]
            extracted_variance_spectrum = df["varianceSpectrum"]
            extracted_snr = df["snr"]

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
        # plt.show()

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "EXTRACTED_ORDERS_QC_PLOT",
            "file_name": filename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"QC plot of extracted source",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        # MERGE THE ORDER SPECTRA
        extractedOrdersDF = pd.concat(extractions, ignore_index=True)
        mergedSpectumDF = self.merge_extracted_orders(extractedOrdersDF)

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

        # DISCRETE ORDERS
        filename = self.filenameTemplate.replace(".fits", "_EXTRACTED_ORDERS.fits")
        filePath = f"{self.outDir}/{filename}"
        hduList.writeto(filePath, checksum=True, overwrite=True)

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "EXTRACTED_ORDERS_TABLE",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"Table of the extracted source in each order",
            "file_path": filePath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

        # NOW MERGED SPECTRUM
        filename = self.filenameTemplate.replace(".fits", "_EXTRACTED_MERGED.fits")
        filePath = f"{self.outDir}/{filename}"
        mergedTable = Table.from_pandas(mergedSpectumDF)
        BinTableHDU = fits.table_to_hdu(mergedTable)
        hduList = fits.HDUList([priHDU, BinTableHDU])
        hduList.writeto(filePath, checksum=True, overwrite=True)

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "EXTRACTED_MERGED_TABLE",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"Table of the extracted, order-merged",
            "file_path": filePath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``extract`` method')
        return self.qc, self.products

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

        # ASTROPY HAS RESET LOGGING LEVEL -- FIX
        import logging
        logging.getLogger().setLevel(logging.INFO + 5)

        kw = self.kw
        arm = self.arm

        # PARAMETERS FROM INPUT FILE
        # THIS IS THE STEP SIZE IN NM (0.06 nm IS SIMILAR TO XHSOOTER EXTRACTION)
        stepWavelengthOrderMerge = 0.06
        ratio = 1 / stepWavelengthOrderMerge
        order_list = []

        # SORT BY COLUMN NAME
        extractedOrdersDF.sort_values(['wavelengthMedian'],
                                      ascending=[True], inplace=True)

        # DEFINE THE WAVELENGTH ARRAY
        wave_resample_grid = np.arange(float(format(np.min(extractedOrdersDF['wavelengthMedian']) * ratio, '.0f')) / ratio, float(format(np.max(extractedOrdersDF['wavelengthMedian']) * ratio, '.0f')) / ratio, step=stepWavelengthOrderMerge)
        wave_resample_grid = wave_resample_grid * u.nm

        # INTERPOLATE THE ORDER SPECTRUM INTO THIS NEW ARRAY WITH A SINGLE STEP SIZE
        flux_orig = extractedOrdersDF['extractedFluxOptimal'].values * u.electron
        # PASS ORIGINAL RAW SPECTRUM AND RESAMPLE
        spectrum_orig = Spectrum1D(flux=flux_orig, spectral_axis=extractedOrdersDF['wavelengthMedian'].values * u.nm, uncertainty=VarianceUncertainty(extractedOrdersDF["varianceSpectrum"]))
        resampler = FluxConservingResampler()
        flux_resampled = resampler(spectrum_orig, wave_resample_grid)
        merged_orders = pd.DataFrame()
        merged_orders['FLUX_COUNTS'] = flux_resampled.flux
        merged_orders['WAVE'] = flux_resampled.spectral_axis

        fig = plt.figure(figsize=(16, 2), constrained_layout=True, dpi=320)
        gs = fig.add_gridspec(1, 1)
        toprow = fig.add_subplot(gs[0:1, :])
        toprow.legend(loc='lower right', bbox_to_anchor=(1, -0.5),
                      fontsize=8)
        toprow.set_ylabel('flux ($e^{-}$)', fontsize=10)
        toprow.set_xlabel(f'wavelength (nm)', fontsize=10)
        toprow.set_title(
            f"Optimally Extracted Order-Merged Object Spectrum ({arm.upper()})", fontsize=11)

        for order in extractedOrdersDF['order'].unique():
            # LIMIT DATAFRAME TO JUST THIS ORDER
            orderDF = extractedOrdersDF.loc[extractedOrdersDF['order'] == order]
            plt.plot(orderDF['wavelengthMedian'], orderDF['extractedFluxOptimal'])

        plt.plot(merged_orders['WAVE'], merged_orders['FLUX_COUNTS'], linewidth=0.2, color="#002b36")

        filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED_QC_PLOT.pdf")
        filePath = f"{self.qcDir}/{filename}"
        # plt.tight_layout()
        # plt.show()
        plt.savefig(filePath, dpi='figure')

        # plt.show()
        # self.log.debug('completed the ``merge_extracted_orders`` method')

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": "EXTRACTED_MERGED_QC_PLOT",
            "file_name": filename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"QC plot of extracted order-merged source",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        return merged_orders

    # use the tab-trigger below for new method
    # xt-class-method


def extract_single_order(crossDispersionSlices, log, ron, slitHalfLength, clippingSigma, clippingIterationLimit, globalClippingSigma):
    """
    *extract the object spectrum for a single order*

    **Return:**
        - ``crossDispersionSlices`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    """
    log.debug('starting the ``extract_single_order`` method')

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
    sys.stdout.flush()
    sys.stdout.write("\x1b[1A\x1b[2K")
    log.print(f"\t## SLICING ORDER INTO CROSS-DISPERSION SLICES - ORDER {order}")

    # REMOVE SLICES WITH ALL NANS
    crossDispersionSlices["sliceRawFlux"] = [np.nan if np.isnan(x).all() else x for x in crossDispersionSlices["sliceRawFlux"]]
    crossDispersionSlices.dropna(axis='index', how='any', subset=["sliceRawFlux"], inplace=True)

    crossDispersionSlices = crossDispersionSlices.apply(lambda x: create_cross_dispersion_slice(x), axis=1)
    crossDispersionSlices["sliceMask"] = [x.mask for x in crossDispersionSlices["sliceRawFluxMasked"]]
    crossDispersionSlices["sliceRawFluxMaskedSum"] = [x.sum() for x in crossDispersionSlices["sliceRawFluxMasked"]]

    # WEIGHTS ARE NOT YET USED
    crossDispersionSlices["sliceWeights"] = ron + np.abs(crossDispersionSlices["sliceRawFlux"]) / (crossDispersionSlices["sliceRawFluxMaskedSum"].pow(2))

    # NORMALISE THE FLUX
    crossDispersionSlices["sliceFluxNormalised"] = crossDispersionSlices["sliceRawFlux"] / crossDispersionSlices["sliceRawFluxMaskedSum"]
    crossDispersionSlices["sliceFluxNormalisedSum"] = [x.sum() for x in crossDispersionSlices["sliceFluxNormalised"]]

    # VERTICALLY STACK THE SLICES INTO PSEDUO-RECTIFIED IMAGE
    fluxNormalisedImage = np.vstack(crossDispersionSlices["sliceFluxNormalised"])
    weightImage = np.vstack(crossDispersionSlices["sliceWeights"])
    maskImage = np.vstack(crossDispersionSlices["sliceMask"])
    errorImage = np.vstack(crossDispersionSlices["sliceError"])

    # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
    # ITERATE FIRST PIXEL IN EACH SLICE AND THEN MOVE TO NEXT
    sys.stdout.flush()
    sys.stdout.write("\x1b[1A\x1b[2K")
    log.print(f"\t## FITTING CROSS-SLIT FLUX NORMALISED PROFILES - ORDER {order}")
    for slitPixelIndex in range(0, slitHalfLength * 2):

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
        while (iteration < clippingIterationLimit) and (clipped_count > 0):

            coeff = np.polyfit(wave_px, fractions, deg=2)
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
    crossDispersionSlices["sliceFittedProfile"] = [np.array(t) for t in transposedProfiles]

    sys.stdout.flush()
    sys.stdout.write("\x1b[1A\x1b[2K")
    log.print(f"\t## EXTRACTING THE SPECTRUM - ORDER {order}")

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
    crossDispersionSlices["wavelengthMedian"] = [np.median(x) for x in crossDispersionSlices["wavelengthGood"]]
    crossDispersionSlices["fudged"] = False

    # TODO: IF ALL ABOVE sliceVarianceRejectLimit  ... what?
    mask = (crossDispersionSlices["horneDenominatorSum"] == 0)
    crossDispersionSlices.loc[mask, "horneNumerator"] = crossDispersionSlices.loc[mask, "sliceRawFlux"] * crossDispersionSlices.loc[mask, "sliceFittedProfileNormalised"] / crossDispersionSlices.loc[mask, "sliceVariance"]
    crossDispersionSlices.loc[mask, 'horneNumeratorSum'] = [x.sum() for x in crossDispersionSlices.loc[mask, "horneNumerator"]]
    crossDispersionSlices.loc[mask, "horneDenominator"] = np.power(crossDispersionSlices.loc[mask, "sliceFittedProfileNormalised"], 2) / crossDispersionSlices.loc[mask, "sliceVariance"]
    crossDispersionSlices.loc[mask, 'horneDenominatorSum'] = [x.sum() for x in crossDispersionSlices.loc[mask, "horneDenominator"]]
    crossDispersionSlices.loc[mask, "fudged"] = True

    # CALULCATE THE FINAL EXTRACTED SPECTRA
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

    crossDispersionSlices.dropna(how="any", subset=["pixelScaleNm", "wavelengthMedian"], inplace=True)

    log.debug('completed the ``extract_single_order`` method')

    return crossDispersionSlices[['order', 'xcoord_centre', 'ycoord', 'wavelengthMedian', 'pixelScaleNm', 'varianceSpectrum', 'snr', 'extractedFluxOptimal', 'extractedFluxBoxcar', 'extractedFluxBoxcarRobust']]


def create_cross_dispersion_slice(
        series):
    """This function is used to create a single, 1-pixel wide cross-dispersion slice of object data. When applied to the dataframe, a single slice is created for each discrete pixel position in the dispersion direction
    """

    import numpy as np
    from astropy.stats import sigma_clip

    # SIGMA-CLIP THE DATA TO REMOVE COSMIC/BAD-PIXELS
    series["sliceRawFluxMasked"] = sigma_clip(
        series["sliceRawFlux"], sigma_lower=7, sigma_upper=7, maxiters=1, cenfunc='median', stdfunc="mad_std")

    return series
