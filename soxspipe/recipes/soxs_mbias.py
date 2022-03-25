#!/usr/bin/env python
# encoding: utf-8
"""
*The recipe for creating master-bias frames *

:Author:
    David Young & Marco Landoni

:Date Created:
    January 22, 2020
"""
################# GLOBAL IMPORTS ####################
from astropy.stats import sigma_clip, mad_std
from soxspipe.commonutils.toolkit import generic_quality_checks
from datetime import datetime
from soxspipe.commonutils import keyword_lookup
import ccdproc
from astropy import units as u
from astropy.nddata import CCDData
import math
import numpy as np
from ._base_recipe_ import _base_recipe_
from soxspipe.commonutils import set_of_files
from fundamentals import tools
from builtins import object
import sys
import os
from astropy.stats import sigma_clip, mad_std
from scipy.stats import median_absolute_deviation
import matplotlib.pyplot as plt
os.environ['TERM'] = 'vt100'


class soxs_mbias(_base_recipe_):
    """
    *The* `soxs_mbias` *recipe is used to generate a master-bias frame from a set of input raw bias frames. The recipe is used only for the UV-VIS arm as NIR frames have bias (and dark current) removed by subtracting an off-frame of equal expsoure length.*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*

    **Usage**

    ```python
    from soxspipe.recipes import soxs_mbias
    mbiasFrame = soxs_mbias(
        log=log,
        settings=settings,
        inputFrames=fileList
    ).produce_product()
    ```

    ---

    ```eval_rst
    .. todo::
        - add a tutorial about ``soxs_mbias`` to documentation
    ```
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[],
            verbose=False

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_mbias, self).__init__(log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mbias' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-mbias"
        self.recipeSettings = settings[self.recipeName]
        # xt-self-arg-tmpx

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MBIAS - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # print("\n# RAW INPUT BIAS FRAMES - SUMMARY")
        # SORT IMAGE COLLECTION
        self.inputFrames.sort(['mjd-obs'])
        # print(self.inputFrames.summary, "\n")

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])

        return None

    def verify_input_frames(
            self):
        """*verify the input frame match those required by the soxs_mbias recipe*

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        # MIXED INPUT IMAGE TYPES ARE BAD
        if len(imageTypes) > 1:
            imageTypes = " and ".join(imageTypes)
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of %(imageTypes)s" % locals())
        # NON-BIAS INPUT IMAGE TYPES ARE BAD
        elif imageTypes[0] != 'BIAS':
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames not BIAS frames" % locals())

        self.imageType = imageTypes[0]

        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*generate a master bias frame*

        **Return:**
            - ``productPath`` -- the path to the master bias frame
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        combined_bias_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mbias")
        self.dateObs = combined_bias_mean.header[kw("DATE_OBS")]

        self.qc_periodic_pattern_noise(frames=self.inputFrames)

        self.qc_ron(
            frameType="MBIAS",
            frameName="master bias",
            masterFrame=combined_bias_mean
        )
        self.qc_bias_structure(combined_bias_mean)

        # INSPECTING THE THE UNCERTAINTY MAPS
        # print("individual frame data")
        # for a in ccds:
        #     print(a.data[0][0])
        # print("\ncombined frame data")
        # print(combined_bias_mean.data[0][0])
        # print("individual frame error")
        # for a in ccds:
        #     print(a.uncertainty[0][0])
        # print("combined frame error")
        # print(combined_bias_mean.uncertainty[0][0])

        # combined_bias_mean.data = combined_bias_mean.data.astype('float32')
        # combined_bias_mean.uncertainty = combined_bias_mean.uncertainty.astype(
        #     'float32')

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log, frame=combined_bias_mean, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)

        medianFlux = self.qc_median_flux_level(
            frame=combined_bias_mean,
            frameType="MBIAS",
            frameName="master bias"
        )

        # WRITE TO DISK
        productPath = self._write(
            frame=combined_bias_mean,
            filedir=self.intermediateRootPath,
            filename=False,
            overwrite=True
        )
        filename = os.path.basename(productPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.products = self.products.append({
            "soxspipe_recipe": self.recipeName,
            "product_label": "MBIAS",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} Master bias frame",
            "file_path": productPath
        }, ignore_index=True)

        self.report_output()
        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    def qc_bias_structure(
            self,
            combined_bias_mean):
        """*calculate the structure of the bias*

        **Key Arguments:**
            - ``combined_bias_mean`` -- the mbias frame

        **Return:**
            - ``structx`` -- slope of BIAS in X direction
            - ``structx`` -- slope of BIAS in Y direction

        **Usage:**

        ```python
        structx, structy = self.qc_bias_structure(combined_bias_mean)
        ```
        """
        self.log.debug('starting the ``qc_bias_structure`` method')
        plot = False

        collaps_ax1 = np.sum(combined_bias_mean, axis=0)
        collaps_ax2 = np.sum(combined_bias_mean, axis=1)

        x_axis = np.linspace(0, len(collaps_ax1), len(collaps_ax1), dtype=int)
        y_axis = np.linspace(0, len(collaps_ax2), len(collaps_ax2), dtype=int)

        # Fitting with a line and collect the slope
        coeff_ax1 = np.polyfit(x_axis, collaps_ax1, deg=1)
        coeff_ax2 = np.polyfit(y_axis, collaps_ax2, deg=1)

        if plot == True:
            plt.plot(x_axis, collaps_ax1)
            plt.plot(x_axis, np.polyval(coeff_ax1, x_axis))
            plt.xlabel('x-axis')
            plt.ylabel('Summed Pixel Values')
            plt.show()

            plt.plot(y_axis, collaps_ax2)
            plt.plot(y_axis, np.polyval(coeff_ax2, y_axis))
            plt.xlabel('y-axis')
            plt.ylabel('Summed Pixel Values')
            plt.show()

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "STRUCTX",
            "qc_value": coeff_ax1[0],
            "qc_comment": "Slope of BIAS in X direction",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }, ignore_index=True)

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "STRUCTY",
            "qc_value": coeff_ax2[0],
            "qc_comment": "Slope of BIAS in Y direction",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }, ignore_index=True)

        self.log.debug('completed the ``qc_bias_structure`` method')
        return coeff_ax1[0], coeff_ax2[0]

    def qc_periodic_pattern_noise(
            self,
            frames):
        """*calculate the periodic pattern noise based on the raw input bias frames*

        A 2D FFT is applied to each of the raw bias frames and the standard deviation and median absolute deviation calcualted for each result. The maximum std/mad is then added as the ppnmax QC in the master bias frame header.

        **Key Arguments:**
            - ``frames`` -- the raw bias frames (imageFileCollection)

        **Return:**
            - ``ppnmax``

        **Usage:**

        ```python
        self.qc_periodic_pattern_noise(frames=self.inputFrames)
        ```
        """
        self.log.debug('starting the ``qc_periodic_pattern_noise`` method')

        # LIST OF CCDDATA OBJECTS
        ccds = [c for c in frames.ccds(ccd_kwargs={"hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        ratios = []
        for frame in ccds:
            # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
            maskedDataArray = np.ma.array(frame.data, mask=frame.mask)
            dark_image_grey_fourier = np.fft.fftshift(np.fft.fft2(maskedDataArray.filled(np.median(frame.data))))

            # SIGMA-CLIP THE DATA
            masked_dark_image_grey_fourier = sigma_clip(
                dark_image_grey_fourier, sigma_lower=1000, sigma_upper=1000, maxiters=1, cenfunc='mean')
            goodData = np.ma.compressed(masked_dark_image_grey_fourier)

            # frame_mad = median_absolute_deviation(dark_image_grey_fourier, axis=None)
            # frame_std = np.std(dark_image_grey_fourier)
            # print(frame_std, frame_mad, frame_std / frame_mad)

            frame_mad = median_absolute_deviation(goodData, axis=None)
            frame_std = np.std(goodData)

            from soxspipe.commonutils.toolkit import quicklook_image
            quicklook_image(
                log=self.log, CCDObject=abs(masked_dark_image_grey_fourier), show=False, ext=None, stdWindow=0.1)

            ratios.append(frame_std / frame_mad)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        ppnmax = max(ratios)

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "PPNMAX",
            "qc_value": ppnmax,
            "qc_comment": "Max periodic pattern noise ratio in raw bias frames",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }, ignore_index=True)

        self.log.debug('completed the ``qc_periodic_pattern_noise`` method')
        return ppnmax

    # use the tab-trigger below for new method
    # xt-class-method
