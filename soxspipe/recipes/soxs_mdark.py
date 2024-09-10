#!/usr/bin/env python
# encoding: utf-8
"""
*The recipe to generate a master dark frame*

Author
: David Young & Marco Landoni

Date Created
: January 27, 2020
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup

from .base_recipe import base_recipe

from fundamentals import tools
from builtins import object
from datetime import datetime
from soxspipe.commonutils.toolkit import generic_quality_checks
import sys
import os
os.environ['TERM'] = 'vt100'


class soxs_mdark(base_recipe):
    """
    *The `soxs_mdark` recipe generates a master-dark frame used to remove flux attributed to the dark current from other frames.*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*

    **Usage**

    ```python
    from soxspipe.recipes import soxs_mdark
    mdarkFrame = soxs_mdark(
        log=log,
        settings=settings,
        inputFrames=fileList,
        verbose=False,
        overwrite=False
    ).produce_product()
    ```
    """

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[],
            verbose=False,
            overwrite=False

    ):
        # INHERIT INITIALISATION FROM  base_recipe
        super(soxs_mdark, self).__init__(log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-mdark")
        self.log = log
        log.debug("instantiating a new 'soxs_mdark' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        # xt-self-arg-tmpx

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        from soxspipe.commonutils.set_of_files import set_of_files
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames,
            ext=self.settings['data-extension']
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MDARK - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # SORT IMAGE COLLECTION
        self.inputFrames.sort(['MJD-OBS'])
        if self.verbose:
            self.log.print("# RAW INPUT FRAMES - SUMMARY")
            self.log.print(self.inputFrames.summary, "\n")

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])

        return None

    def verify_input_frames(
            self):
        """*verify input frame match those required by the soxs_mdark recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        error = False

        # MIXED INPUT IMAGE TYPES ARE BAD
        if not error:
            if len(imageTypes) > 1:
                imageTypes = " and ".join(imageTypes)
                error = "Input frames are a mix of %(imageTypes)s" % locals()
            # NON-BIAS INPUT IMAGE TYPES ARE BAD
            elif imageTypes[0] != 'DARK':
                error = "Input frames not DARK frames" % locals()

        if not error:
            exptimes = self.inputFrames.values(
                keyword=kw("EXPTIME"), unique=True)
            # MIXED INPUT IMAGE TYPES ARE BAD
            if len(exptimes) > 1:
                exptimes = [str(e) for e in exptimes]
                exptimes = " and ".join(exptimes)
                error = "Input frames have differing exposure-times %(exptimes)s" % locals()

        if error:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            self.log.print("")
            raise TypeError(error)

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*generate a master dark frame*

        **Return:**

        - ``productPath`` -- the path to master dark frame
        """
        self.log.debug('starting the ``produce_product`` method')

        import numpy as np
        import pandas as pd

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        # LIST OF CCDDATA OBJECTS
        ccds = [c for c in self.inputFrames.ccds(ccd_kwargs={"hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        meanFluxLevels, rons, noiseFrames = zip(*[self.subtract_mean_flux_level(c) for c in ccds])
        masterMeanFluxLevel = np.mean(meanFluxLevels)
        masterMedianFluxLevel = np.median(meanFluxLevels)
        rawRon = np.mean(rons)

        combined_noise = self.clip_and_stack(
            frames=list(noiseFrames), recipe="soxs_mdark", ignore_input_masks=False, post_stack_clipping=True)

        maskedDataArray = np.ma.array(combined_noise.data, mask=combined_noise.mask)
        masterRon = np.std(maskedDataArray)

        # FILL MASKED PIXELS WITH 0
        combined_noise.data = np.ma.array(combined_noise.data, mask=combined_noise.mask, fill_value=0).filled() + masterMeanFluxLevel
        combined_noise.uncertainty = np.ma.array(combined_noise.uncertainty.array, mask=combined_noise.mask, fill_value=rawRon).filled()
        combined_dark_mean = combined_noise
        combined_dark_mean.mask = combined_noise.mask

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log, frame=combined_dark_mean, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)

        medianFlux = self.qc_median_flux_level(
            frame=combined_dark_mean,
            frameType="MDARK",
            frameName="master dark",
            medianFlux=masterMedianFluxLevel
        )

        self.qc_ron(
            frameType="DARK"
        )

        self.update_fits_keywords(
            frame=combined_dark_mean
        )

        # WRITE TO DISK
        productPath = self._write(
            frame=combined_dark_mean,
            filedir=self.workspaceRootPath,
            filename=False,
            overwrite=True
        )
        filename = os.path.basename(productPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.dateObs = combined_dark_mean.header[kw("DATE_OBS")]

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "MDARK",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} Master dark frame",
            "file_path": productPath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

        self.report_output()
        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method
