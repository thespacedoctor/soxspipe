#!/usr/bin/env python
"""
*The recipe for creating master-bias frames *

Author
: David Young & Marco Landoni

Date Created
: January 22, 2020
"""

################# GLOBAL IMPORTS ####################

import os
import sys

from soxspipe.commonutils.toolkit import generic_quality_checks, utcnow_string

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"


class soxs_mbias(base_recipe):
    """
    *The* `soxs_mbias` *recipe is used to generate a master-bias frame from a set of input raw bias frames. The recipe is used only for the UV-VIS arm as NIR frames have bias (and dark current) removed by subtracting an off-frame of equal exposure length.*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- debug mode. True or False. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.

    **Usage**

    ```python
    from soxspipe.recipes import soxs_mbias
    mbiasFrame = soxs_mbias(
        log=log,
        settings=settings,
        inputFrames=fileList
    ).produce_product()
    ```
    """

    # Initialisation

    def __init__(
        self,
        log,
        settings=False,
        inputFrames=[],
        verbose=False,
        overwrite=False,
        command=False,
        debug=False,
        turnOffMP=False,
    ):
        # INHERIT INITIALISATION FROM  base_recipe
        this = super().__init__(
            log=log,
            settings=settings,
            inputFrames=inputFrames,
            overwrite=overwrite,
            recipeName="soxs-mbias",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        log.debug("instantiating a new 'soxs_mbias' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames,
            ext=self.settings["data-extension"],
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MBIAS - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # SORT IMAGE COLLECTION
        self.inputFrames.sort(["MJD-OBS"])

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        return

    def verify_input_frames(self):
        """*verify the input frame match those required by the soxs_mbias recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        error = False

        # MIXED INPUT IMAGE TYPES ARE BAD
        if len(imageTypes) > 1:
            error = "Input frames are a mix of %(imageTypes)s" % locals()
        # NON-BIAS INPUT IMAGE TYPES ARE BAD
        elif imageTypes[0] != "BIAS":
            error = "Input frames not BIAS frames" % locals()

        if error:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.error("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            self.log.print("")
            raise TypeError(error)

        self.imageType = imageTypes[0]

        self.log.debug("completed the ``verify_input_frames`` method")
        return

    def produce_product(self):
        """*generate a master bias frame*

        **Return:**

        - ``productPath`` -- the path to the master bias frame
        """
        self.log.debug("starting the ``produce_product`` method")

        from soxspipe.commonutils import toolkit

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        (
            combined_bias_mean,
            masterMedianBiasLevel,
            rawRon,
            masterRon,
        ) = self._combine_bias_frames()

        # OPTIMISE: 24%
        self.qc_periodic_pattern_noise(frames=self.inputFrames)

        self.qc_ron(
            frameType="MBIAS",
            frameName="master bias",
            masterFrame=combined_bias_mean,
            rawRon=rawRon,
            masterRon=masterRon,
        )

        self.qc_bias_structure(combined_bias_mean)

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log,
            frame=combined_bias_mean,
            settings=self.settings,
            recipeName=self.recipeName,
            qcTable=self.qc,
        )

        medianFlux = self.qc_median_flux_level(
            frame=combined_bias_mean,
            frameType="MBIAS",
            frameName="master bias",
            medianFlux=masterMedianBiasLevel,
        )

        self.update_fits_keywords(frame=combined_bias_mean)

        # WRITE TO DISK
        toolkit.frame_to_32(combined_bias_mean)
        productPath = self._write(
            frame=combined_bias_mean,
            filedir=self.workspaceRootPath,
            filename=False,
            overwrite=True,
        )
        filename = os.path.basename(productPath)

        utcnow = utcnow_string()

        self.dateObs = combined_bias_mean.header[self.kw("DATE_OBS")]

        self.add_product(
            productLabel="MBIAS",
            fileName=filename,
            filePath=productPath,
            productDesc=f"{self.arm} Master bias frame",
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        qcTable = self.report_output()
        self.clean_up()

        self.log.debug("completed the ``produce_product`` method")
        return productPath, qcTable

    def _combine_bias_frames(self):
        """*stack the raw bias frames into the master-bias frame*

        **Return:**

        - ``combined_bias_mean`` -- the stacked master-bias frame
        - ``masterMedianBiasLevel`` -- the median of the per-frame mean bias levels
        - ``rawRon`` -- the mean of the per-frame read-out noise values
        - ``masterRon`` -- the standard deviation of the stacked noise frame

        **Usage:**

        ```python
        combined_bias_mean, masterMedianBiasLevel, rawRon, masterRon = self._combine_bias_frames()
        ```
        """
        self.log.debug("starting the ``_combine_bias_frames`` method")

        import numpy as np

        # LIST OF CCDDATA OBJECTS
        # OPTIMISE: 9%
        ccds = [
            c
            for c in self.inputFrames.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]

        # OPTIMISE: 33%
        meanBiasLevels, rons, noiseFrames = zip(*[self.subtract_mean_flux_level(c) for c in ccds])
        masterMeanBiasLevel = np.mean(meanBiasLevels)
        masterMedianBiasLevel = np.median(meanBiasLevels)
        rawRon = np.mean(rons)

        # OPTIMISE: 19%
        combined_noise = self.clip_and_stack(
            frames=list(noiseFrames),
            recipe="soxs_mbias",
            ignore_input_masks=False,
            post_stack_clipping=True,
        )

        masterRon = np.std(combined_noise.data)

        # USE COMBINED NOISE MASK AS MBIAS MASK
        combined_noise.data = (
            np.ma.array(combined_noise.data, mask=combined_noise.mask, fill_value=0).filled() + masterMeanBiasLevel
        )
        combined_noise.uncertainty = np.ma.array(
            combined_noise.uncertainty.array,
            mask=combined_noise.mask,
            fill_value=rawRon,
        ).filled()
        combined_bias_mean = combined_noise
        combined_bias_mean.mask = combined_noise.mask

        self.log.debug("completed the ``_combine_bias_frames`` method")
        return combined_bias_mean, masterMedianBiasLevel, rawRon, masterRon

    def qc_bias_structure(self, combined_bias_mean):
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
        self.log.debug("starting the ``qc_bias_structure`` method")

        import numpy as np

        plot = False

        collaps_ax1 = np.nansum(combined_bias_mean, axis=0)
        collaps_ax2 = np.nansum(combined_bias_mean, axis=1)

        x_axis = np.linspace(0, len(collaps_ax1), len(collaps_ax1), dtype=int)
        y_axis = np.linspace(0, len(collaps_ax2), len(collaps_ax2), dtype=int)

        # Fitting with a line and collect the slope
        coeff_ax1 = np.polyfit(x_axis, collaps_ax1, deg=1)
        coeff_ax2 = np.polyfit(y_axis, collaps_ax2, deg=1)

        if plot == True:
            import matplotlib.pyplot as plt

            plt.plot(x_axis, collaps_ax1)
            plt.plot(x_axis, np.polyval(coeff_ax1, x_axis))
            plt.xlabel("x-axis")
            plt.ylabel("Summed Pixel Values")
            plt.show()

            plt.plot(y_axis, collaps_ax2)
            plt.plot(y_axis, np.polyval(coeff_ax2, y_axis))
            plt.xlabel("y-axis")
            plt.ylabel("Summed Pixel Values")
            plt.show()

        # ONE TIMESTAMP DELIBERATELY SHARED ACROSS BOTH QC ROWS BELOW
        utcnow = utcnow_string()

        self.add_qc(
            qcName="STRUCTX",
            qcValue=coeff_ax1[0],
            qcComment="Slope of BIAS in X direction",
            reductionDateUtc=utcnow,
            qcUnit=None,
            toHeader=True,
        )

        self.add_qc(
            qcName="STRUCTY",
            qcValue=coeff_ax2[0],
            qcComment="Slope of BIAS in Y direction",
            reductionDateUtc=utcnow,
            qcUnit=None,
            toHeader=True,
        )

        self.log.debug("completed the ``qc_bias_structure`` method")
        return coeff_ax1[0], coeff_ax2[0]

    def _periodic_noise_ratio(self, frame):
        """*return the periodic-pattern-noise ratio of a single raw bias frame*

        **Key Arguments:**

        - ``frame`` -- a single raw bias frame (CCDData)

        **Return:**

        - ``ratio`` -- the standard deviation of the sigma-clipped 2D FFT divided by its median absolute deviation

        **Usage:**

        ```python
        ratio = self._periodic_noise_ratio(frame)
        ```
        """
        import numpy as np
        from astropy.stats import sigma_clip
        from ccdproc import block_reduce
        from scipy.stats import median_abs_deviation

        from soxspipe.commonutils.toolkit import quicklook_image

        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        maskedDataArray = np.ma.array(frame.data, mask=frame.mask)
        # BIN THE FRAME TO INCREASE SPEED
        maskedDataArray = block_reduce(maskedDataArray, 5, np.mean)
        dark_image_grey_fourier = np.fft.fftshift(np.fft.fft2(maskedDataArray.filled(np.median(frame.data))))

        # SIGMA-CLIP THE DATA
        masked_dark_image_grey_fourier = sigma_clip(
            dark_image_grey_fourier,
            sigma_lower=100,
            sigma_upper=100,
            maxiters=1,
            cenfunc="mean",
        )
        goodData = np.ma.compressed(masked_dark_image_grey_fourier)

        frame_mad = median_abs_deviation(goodData, axis=None)
        frame_std = np.std(goodData)

        quicklook_image(
            log=self.log,
            CCDObject=abs(masked_dark_image_grey_fourier),
            show=False,
            ext=None,
            stdWindow=0.1,
        )
        quicklook_image(
            log=self.log,
            CCDObject=abs(dark_image_grey_fourier),
            show=False,
            ext=None,
            stdWindow=0.1,
        )

        return frame_std / frame_mad

    def qc_periodic_pattern_noise(self, frames):
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
        self.log.debug("starting the ``qc_periodic_pattern_noise`` method")

        # LIST OF CCDDATA OBJECTS
        ccds = [
            c
            for c in frames.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]

        ratios = [self._periodic_noise_ratio(frame) for frame in ccds]

        utcnow = utcnow_string()

        ppnmax = max(ratios)

        self.add_qc(
            qcName="FPN FRACMAX",
            qcValue=ppnmax,
            qcComment="Max periodic pattern noise ratio in raw bias frames",
            reductionDateUtc=utcnow,
            qcUnit=None,
            toHeader=True,
        )

        self.log.debug("completed the ``qc_periodic_pattern_noise`` method")
        return ppnmax

    # use the tab-trigger below for new method
    # xt-class-method
