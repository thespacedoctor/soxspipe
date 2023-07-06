#!/usr/bin/env python
# encoding: utf-8
"""
*perform optimal source extraction using the Horne method (Horne 1986)*

:Author:
    Marco Landoni & David Young

:Date Created:
    May 17, 2023
"""

# Please do import into the methods for speed reasons.

from fundamentals import tools
from builtins import object
import sys
import os
from astropy import units as u
from astropy.nddata import CCDData


os.environ['TERM'] = 'vt100'


class horne_extraction(object):
    """
    *The worker class for the horne_extraction module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``skyModelFrame`` -- sky model frame
        - ``skySubFrame`` -- sky subtracted frame
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
        - ``twoDMapPath`` -- 2D dispersion map image path
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products

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

        self.log = log
        log.debug("instansiating a new 'horne_extraction' object")
        self.dispersionMap = dispersionMap
        self.twoDMapPath = twoDMapPath
        self.settings = settings
        self.productsTable = productsTable
        self.qcTable = qcTable
        self.recipeName = recipeName
        self.sofName = sofName

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
        hdul = fits.open(twoDMapPath)
        # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        xdim = hdul[0].data.shape[1]
        ydim = hdul[0].data.shape[0]
        xarray = np.tile(np.arange(0, xdim), ydim)
        yarray = np.repeat(np.arange(0, ydim), xdim)
        self.imageMap = pd.DataFrame.from_dict({
            "x": xarray,
            "y": yarray,
            "wavelength": hdul["WAVELENGTH"].data.flatten().byteswap().newbyteorder(),
            "slit_position": hdul["SLIT"].data.flatten().byteswap().newbyteorder(),
            "order": hdul["ORDER"].data.flatten().byteswap().newbyteorder(),
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
        # UNPACK THE ORDER TABLE (CENTRE LOCATION ONLY AT THIS STAGE)
        orderPolyTable, self.orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=productPath)

    def extract(self, order):
        """
        *get the horne_extraction object*

        **Return:**
            - ``extracted_wave_spectrum`` -- pixel wavelength array
            - ``extracted_spectrum`` -- pixel flux values

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

        # EXTRACTION STARTING
        i = 0
        extracted_spectrum = []
        # NON-OPTIMAL EXTRACTED SPECTRUM
        extracted_spectrum_nonopt = []
        extracted_wave_spectrum = []
        total_weights = []

        # TODO: FIX NAN = 0.0 ... COULD BE AFFECTING SLICE PIXEL-FLUX FRACTIONS
        ourFrame_nonan = np.nan_to_num(self.skySubtractedFrame.data, nan=0.0, copy=True)
        skyModelFrame_nonan = np.nan_to_num(self.skyModelFrame.data, nan=0.0, copy=True)

        ourFrame_nonan = self.skySubtractedFrame.data.copy()
        skyModelFrame_nonan = self.skyModelFrame.data.copy()

        poly_coeffs = []

        # 1) SELECTING THE ORDER FROM THE 2DMAP TABLE
        # THIS IS THE CONTINUUM OF THE OBJECT
        selection = self.orderPixelTable.loc[self.orderPixelTable['order'] == order]

        # APPLY FUNCTION TO COLUMN VALUES
        def create_slices(
                series):
            xcoord_centre = series['xcoord_centre']
            y = series['ycoord']
            std = np.abs(series['std'])
            series["slice"] = ourFrame_nonan[int(y), int(xcoord_centre) - self.slitHalfLength:int(xcoord_centre) + self.slitHalfLength]
            series["weights"] = self.ron + np.abs(series["slice"]) / (series["slice"].sum() * series["slice"].sum())
            if np.isnan(series["slice"]).all():
                series["slice"] = np.nan
                return series

            # NORMALISE
            series["slice_normalised"] = series["slice"] / series["slice"].sum()
            return series

        selection = selection.apply(create_slices, axis=1)
        # DROP COLUMNS WITH ALL NANs
        selection.dropna(axis='index', how='any', subset=['slice'], inplace=True)
        fractionImage = np.vstack(selection["slice_normalised"])
        weightImage = np.vstack(selection["weights"])

        # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
        # iterate first pixel in each slice and then move to next
        for slitPixelIndex in range(0, self.slitHalfLength * 2):

            # NOW FIT (NOTE: CLIPPED FITTING HERE SHALL BE ASSUMED)
            iteration = 1
            clipped_count = 1
            fractions = fractionImage[:, slitPixelIndex]
            startCount = len(fractions)
            wave_px = selection["ycoord"]
            weights = weightImage[:, slitPixelIndex]
            while (iteration < self.clippingIterationLimit) and (clipped_count > 0):

                coeff = np.polyfit(wave_px, fractions, deg=2)
                residuals = fractions - np.polyval(coeff, wave_px)
                masked_residuals = sigma_clip(residuals, sigma_lower=self.clippingSigma, sigma_upper=self.clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
                # REDUCE ARRAYS TO NON-MASKED VALUES
                a = [fractions, wave_px, weights]
                fractions, wave_px, weights = [np.ma.compressed(np.ma.masked_array(
                    i, masked_residuals.mask)) for i in a]

                # REMOVE OUTLIERS
                clipped_count = startCount - len(fractions)
                percent = (float(clipped_count) / float(startCount)) * 100.

                print(f"Iteration {iteration}, slice index {slitPixelIndex+1}/{self.slitHalfLength * 2}. {clipped_count} clipped ({percent:0.2f}%).")
                iteration = iteration + 1
                if iteration > 1:
                    sys.stdout.write("\x1b[1A\x1b[2K")

            # poly_coeffs = collection of multiple polynominals of fractional fluxes along trace
            poly_coeffs.append(coeff)

            if False:
                plt.scatter(wave_px, fractions, alpha=0.2)
                plt.plot(wave_px, np.polyval(coeff, wave_px), color='red')
                plt.ylim([-1, 1])
                plt.show()

        # 3) OPTIMAL EXTRACTION HERE
        for index, row in selection.iterrows():

            # TODO: this is done above and can be move outside
            xcoord_centre = row['xcoord_centre']
            y = row['ycoord']
            std = np.abs(row['std'])
            slice_data = ourFrame_nonan[int(y), int(xcoord_centre) - self.slitHalfLength:int(xcoord_centre) + self.slitHalfLength]
            slice_data_sky = skyModelFrame_nonan[int(y), int(xcoord_centre) - self.slitHalfLength:int(xcoord_centre) + self.slitHalfLength]
            extracted_val = 0.0
            denominator_val = 0.0
            weigths = 0.0

            # CHECKING BOUNDARIES
            if len(slice_data) == 0:
                continue

            # ENFORCING NORMALISATION OF THE PROFILE
            norm_factor = 0.0
            for js in range(0, self.slitHalfLength * 2):
                coeffs = poly_coeffs[js]
                # norm_factor will be almost 1 but not quite (polynomial fitting introduce fluxuation)
                # y == wave_pix in the code above
                norm_factor += np.max([0, np.polyval(coeffs, y)])

            for js in range(0, self.slitHalfLength * 2):
                coeffs = poly_coeffs[js]
                # ALL Pxl need to sum to 1 ... hence norm_factor
                # RENORMALISING ENFORCES +VE FLUX AND FLUX-CONSERVATION
                Pxl = np.max([0, np.polyval(coeffs, y)]) / norm_factor
                # Vxl is the variance estimate for each pixel

                # TODO: this slice_data will be affected by the same NAN=0.0 issue reported above
                # TODO: pass gain from the detection ...  hardwired to 1.0 here
                gain = 1.0
                # FLUX VALUE = RON + OBJECT (PIXEL FRACTION OF THE DATA SLICE) + SKY
                Vxl = self.ron + np.abs(np.sum(slice_data) * Pxl + slice_data_sky[js]) / gain

                # REJECTING COSMIC / OUTLIERS pixel from the summation along the XD
                if (slice_data[js] - np.sum(slice_data) * Pxl) ** 2 / Vxl < self.globalClippingSigma:
                    extracted_val += slice_data[js] * Pxl / Vxl
                    denominator_val += Pxl ** 2 / Vxl
                    weigths += Pxl ** 2 / Vxl

            # TODO: THIS IS A FUDGE FOR WHEN ALL PIXELS IN A SLICE ARE CLIPPED ... NEED TO RETHINK
            if extracted_val == 0.0:
                for js in range(0, self.slitHalfLength * 2):
                    coeffs = poly_coeffs[js]
                    # ALL Pxl need to sum to 1 ... hence norm_factor
                    # RENORMALISING ENFORCES +VE FLUX AND FLUX-CONSERVATION
                    Pxl = np.max([0, np.polyval(coeffs, y)]) / norm_factor
                    # Vxl is the variance estimate for each pixel

                    # TODO: this slice_data will be affected by the same NAN=0.0 issue reported above
                    # TODO: pass gain from the detection ...  hardwired to 1.0 here
                    gain = 1.0
                    # FLUX VALUE = RON + OBJECT (PIXEL FRACTION OF THE DATA SLICE) + SKY
                    Vxl = self.ron + np.abs(np.sum(slice_data) * Pxl + slice_data_sky[js]) / gain

                    extracted_val += slice_data[js] * Pxl / Vxl
                    denominator_val += Pxl ** 2 / Vxl
                    weigths += Pxl ** 2 / Vxl

            # TODO: use this to generate the variance spectrum (see Horne page 612 equ 9)
            total_weights.append(weigths)

            # NOW ASSIGN THE WAVELENGHT AT THE X_COORD CENTRE POSITION WITH A FIT

            # TODO: wave is the approximate wavelegth of 6 pixel in slice - CAN MOVE THIS OUTSIDE FOR LOOP AND POSSIBLY TAKE THE MEAN OF THE SLICE
            wave = self.imageMap[(self.imageMap['order'] == float(order)) & (self.imageMap['y'] == np.round(y)) & (
                self.imageMap['x'] > xcoord_centre - self.slitHalfLength) & (self.imageMap['x'] < xcoord_centre + self.slitHalfLength)]
            if wave.size != 0:
                # EQU 8 of Horne
                extracted_spectrum.append(extracted_val / denominator_val)
                # BOXCAR
                extracted_spectrum_nonopt.append(np.sum(slice_data))
                coeff_wl = np.polyfit(wave['x'], wave['wavelength'], deg=1)
                extracted_wave_spectrum.append(np.polyval(coeff_wl, xcoord_centre))

        self.log.debug('completed the ``extract`` method')

        return (extracted_wave_spectrum, extracted_spectrum, extracted_spectrum_nonopt)

    # xt-class-method
