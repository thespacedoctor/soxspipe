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

        self.slitLength = int(self.settings["soxs-stare"]["horne-extraction-slit-length"] / 2)
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

        from soxspipe.commonutils import detect_continuum
        import yaml
        import pandas as pd
        import numpy as np
        from astropy.stats import sigma_clip
        import matplotlib.pyplot as plt

        # FIX BELOW
        # TODO: fix the detect_continuum function to work directly on the object frame WIHTOUT having to adjust header keywords
        # TODO: continum fitting need to be done in __init__ (or another method)

        # TODO: PRODUCT: /Users/Dave/soxspipe-unittests/intermediate/xsh/product/soxs-stare/2019.08.22T23.12.18.5011_NIR_OBJECT_STARE_EG_274_SLIT.fits --> /Users/Dave/soxspipe-unittests/intermediate/xsh/product/soxs-stare/2019.08.22T23.12.18.5011_NIR_OBJECT_STARE_EG_274_OBJECT_TRACE.fits
        # TODO: Number of order centre traces found > Number of <orders an object> traces found
        # TODO: ORDER_CENTRES_RES > OBJECT_TRACE_RES
        # TODO: 2019.08.22T23.12.18.5011_NIR_OBJECT_STARE_EG_274_SLIT.fitsORDER_CENTRES_residuals.pdf > 2019.08.22T23.12.18.5011_NIR_OBJECT_STARE_EG_274_OBJECT_TRACE_residuals.pdf
        # TODO: FIX PLOT TITLES
        # TODO: FIX TABLE FITS HEADERS

        # HACK TO GET DETECT CONTINUUM TO RUN
        # self.skySubtractedFrame.header["ESO DPR TYPE"] = "LAMP,FLAT"
        # self.skySubtractedFrame.header["ESO DPR TECH"] = "IMAGE"
        # self.settings["intermediate-data-root"] = "./"

        # NOTE TO MARCO - I HAVE TRICKED THE CODE INTO THINKING THIS IS A LAMP-FLAT PINHOLE FRAME
        # SO detect_continuum OUTPUTS 20190831T001327_NIR_MORDER_CENTRES_residuals.pdf AND 20190831T001327_NIR_ORDER_LOCATIONS.fits
        # BUT THESE ARE REALLY THE OBJECT CONTINUUM NOT ORDER CENTRE TRACES

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

        print(productPath)

        from tabulate import tabulate
        print(tabulate(self.qcTable, headers='keys', tablefmt='psql'))
        print(tabulate(self.productsTable, headers='keys', tablefmt='psql'))

        from soxspipe.commonutils.toolkit import unpack_order_table
        # UNPACK THE ORDER TABLE (CENTRE LOCATION ONLY AT THIS STAGE)
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=productPath)

        # FIX ABOVE ^

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

        poly_coeffs = []

        # 1) SELECTING THE ORDER FROM THE 2DMAP TABLE
        # THIS IS THE CONTINUUM OF THE OBJECT
        selection = orderPixelTable.loc[orderPixelTable['order'] == order]

        # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
        # sl_px = pixel index in each slice
        # iterate first pixel in each slice and then move to next
        for sl_px in range(0, self.slitLength * 2):
            fractions = []
            weigths = []
            wave_px = []
            for index, row in selection.iterrows():
                i = i + 1
                xcoord_centre = row['xcoord_centre']
                y = row['ycoord']
                std = np.abs(row['std'])

                # TODO: could move this outside the for loop and store slices in a lookup
                slice_data = ourFrame_nonan[int(y), int(xcoord_centre) - self.slitLength:int(xcoord_centre) + self.slitLength]
                # for j in range(0, len(slice_data)):

                try:

                    # TODO: inverstigate if average pass ... might retrun inf not nan
                    if np.isnan(slice_data[sl_px]) or np.isnan(slice_data[sl_px] / (np.sum(slice_data))):
                        # print('Here')
                        pass
                    else:
                        # QC PLOTS
                        if False and i > 490 and i < 520:
                            x_axis = np.arange(int(xcoord_centre) - self.slitLength, int(xcoord_centre) + self.slitLength)
                            plt.plot(x_axis, slice_data)
                            # plt.plot(x_axis, py)
                            plt.show()
                        # RATIO OF PIXEL FLUX TO SUM IN SLICE
                        # TODO: FIX NAN = 0.0 ... COULD BE AFFECTING SLICE PIXEL-FLUX FRACTIONS (see above)
                        fractions.append(slice_data[sl_px] / (np.sum(slice_data)))
                        wave_px.append(y)
                        # TODO: COULD USE ERROR MAP AND NOT RON? Current NOT using these weights
                        w = self.ron + np.abs(slice_data[sl_px]) / (np.sum(slice_data) * np.sum(slice_data))
                        weigths.append(w)
                except Exception as e:
                    print(e)
                    pass
            # NOW FIT (NOTE: CLIPPED FITTING HERE SHALL BE ASSUMED)
            iteration = 1
            clipped_count = 1
            startCount = len(fractions)
            fractions = np.array(fractions)
            wave_px = np.array(wave_px)
            weigths = np.array(weigths)
            coeffs = []

            while (iteration < self.clippingIterationLimit) and (clipped_count > 0):

                # TODO: add this deg to settings file?
                coeffs = np.polyfit(wave_px, fractions, deg=2)

                # CALCULATE RESIDUALS
                res = fractions - np.polyval(coeffs, wave_px)

                # CLIP

                # TODO: use MAD instread of STD?
                masked_residuals = sigma_clip(
                    res, sigma=self.clippingSigma, maxiters=1)

                # REMOVE OUTLIERS

                remove_mask = np.where(masked_residuals.mask == True)

                fractions = np.delete(fractions, remove_mask)
                wave_px = np.delete(wave_px, remove_mask)
                weigths = np.delete(weigths, remove_mask)
                iteration = iteration + 1
                clipped_count = startCount - len(fractions)
                if iteration > 1:
                    sys.stdout.write("\x1b[1A\x1b[2K")

            # poly_coeffs = collection of multiple polynominals of fractional fluxes along trace
            poly_coeffs.append(coeffs)
            #plt.scatter(wave_px, fractions, alpha=0.2)
            #plt.plot(wave_px, np.polyval(coeffs, wave_px), color='red')
            #plt.ylim([-1, 1])
            # plt.show()

        # 3) OPTIMAL EXTRACTION HERE
        for index, row in selection.iterrows():

            # TODO: this is done above and can be move outside
            xcoord_centre = row['xcoord_centre']
            y = row['ycoord']
            std = np.abs(row['std'])
            slice_data = ourFrame_nonan[int(y), int(xcoord_centre) - self.slitLength:int(xcoord_centre) + self.slitLength]
            slice_data_sky = skyModelFrame_nonan[int(y), int(xcoord_centre) - self.slitLength:int(xcoord_centre) + self.slitLength]
            extracted_val = 0.0
            denominator_val = 0.0
            weigths = 0.0

            # CHECKING BOUNDARIES
            if len(slice_data) == 0:
                continue

            # ENFORCING NORMALISATION OF THE PROFILE
            norm_factor = 0.0
            for js in range(0, self.slitLength * 2):
                coeffs = poly_coeffs[js]
                # norm_factor will be almost 1 but not quite (polynomial fitting introduce fluxuation)
                # y == wave_pix in the code above
                norm_factor += np.max([0, np.polyval(coeffs, y)])

            for js in range(0, self.slitLength * 2):
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
                for js in range(0, self.slitLength * 2):
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
                self.imageMap['x'] > xcoord_centre - self.slitLength) & (self.imageMap['x'] < xcoord_centre + self.slitLength)]
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
