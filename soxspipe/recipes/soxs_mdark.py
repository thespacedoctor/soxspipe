#!/usr/bin/env python
# encoding: utf-8
"""
*The recipe to generate a master dark frame*

:Author:
    David Young & Marco Landoni

:Date Created:
    January 27, 2020
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
import ccdproc
from astropy.nddata import CCDData
import numpy as np
from ._base_recipe_ import _base_recipe_
from soxspipe.commonutils import set_of_files
from fundamentals import tools
from builtins import object
from datetime import datetime
import sys
import os
os.environ['TERM'] = 'vt100'


class soxs_mdark(_base_recipe_):
    """
    *The soxs_mdark recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
        - ``verbose`` -- verbose. True or False. Default *False*

    **Usage**

    ```python
    from soxspipe.recipes import soxs_mdark
    mdarkFrame = soxs_mdark(
        log=log,
        settings=settings,
        inputFrames=fileList
    )..produce_product()
    ```

    ---

    ```eval_rst
    .. todo::

        - add a tutorial about ``soxs_mdark`` to documentation
    ```
    """

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[],
            verbose=False

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_mdark, self).__init__(log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mdark' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-mdark"
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MDARK - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # SORT IMAGE COLLECTION
        self.inputFrames.sort(['mjd-obs'])
        if self.verbose:
            print("# RAW INPUT FRAMES - SUMMARY")
            print(self.inputFrames.summary, "\n")

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])

        return None

    def verify_input_frames(
            self):
        """*verify input frame match those required by the soxs_mdark recipe*

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        self._verify_input_frames_basics()

        imageTypes = self.inputFrames.values(
            keyword=kw("DPR_TYPE").lower(), unique=True)
        # MIXED INPUT IMAGE TYPES ARE BAD
        if len(imageTypes) > 1:
            imageTypes = " and ".join(imageTypes)
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of %(imageTypes)s" % locals())
        # NON-BIAS INPUT IMAGE TYPES ARE BAD
        elif imageTypes[0] != 'DARK':
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames not DARK frames" % locals())

        exptimes = self.inputFrames.values(
            keyword=kw("EXPTIME").lower(), unique=True)
        # MIXED INPUT IMAGE TYPES ARE BAD
        if len(exptimes) > 1:
            exptimes = [str(e) for e in exptimes]
            exptimes = " and ".join(exptimes)
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames have differing exposure-times %(exptimes)s" % locals())

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

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        combined_dark_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mdark")

        medianFlux = self.qc_median_flux_level(
            frame=combined_dark_mean,
            frameType="MDARK",
            frameName="master dark"
        )

        self.qc_dark_ron()

        # WRITE TO DISK
        productPath = self._write(
            frame=combined_dark_mean,
            filedir=self.intermediateRootPath,
            filename=False,
            overwrite=True
        )
        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    def qc_dark_ron(
            self):
        """*calculate the read-out-noise on the raw dark frames*

        **Return:**
            - ``rawRon`` -- raw read-out-noise in electrons

        **Usage:**

        ```python
        darkRon = self.qc_dark_ron()
        ```
        """
        self.log.debug('starting the ``qc_dark_ron`` method')

        # LIST OF CCDDATA OBJECTS
        ccds = [c for c in self.inputFrames.ccds(ccd_kwargs={
                                                 "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        # SINGLE FRAME RON
        raw_one = ccds[0]
        raw_two = ccds[1]
        raw_diff = raw_one.subtract(raw_two)

        raw_diff = raw_diff[100:110, 100:110]

        # # PLOT CCDDATA OBJECT
        # import matplotlib.pyplot as plt
        # rotatedImg = np.rot90(raw_one[100:110, 100:110], 1)
        # std = np.std(raw_one[100:110, 100:110])
        # mean = np.mean(raw_one[100:110, 100:110])
        # vmax = mean + 3 * std
        # vmin = mean - 3 * std
        # plt.figure(figsize=(12, 5))
        # plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
        #            cmap='gray', alpha=1)
        # plt.colorbar()
        # plt.xlabel(
        #     "y-axis", fontsize=10)
        # plt.ylabel(
        #     "x-axis", fontsize=10)
        # plt.show()

        # # PLOT CCDDATA OBJECT
        # import matplotlib.pyplot as plt
        # rotatedImg = np.rot90(raw_two[100:110, 100:110], 1)
        # std = np.std(raw_two[100:110, 100:110])
        # mean = np.mean(raw_two[100:110, 100:110])
        # vmax = mean + 3 * std
        # vmin = mean - 3 * std
        # plt.figure(figsize=(12, 5))
        # plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
        #            cmap='gray', alpha=1)
        # plt.colorbar()
        # plt.xlabel(
        #     "y-axis", fontsize=10)
        # plt.ylabel(
        #     "x-axis", fontsize=10)
        # plt.show()

        # # PLOT CCDDATA OBJECT
        # import matplotlib.pyplot as plt
        # rotatedImg = np.rot90(raw_diff, 1)
        # std = np.std(raw_diff)
        # mean = np.mean(raw_diff)
        # vmax = mean + 3 * std
        # vmin = mean - 3 * std
        # plt.figure(figsize=(12, 5))
        # plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
        #            cmap='gray', alpha=1)
        # plt.colorbar()
        # plt.xlabel(
        #     "y-axis", fontsize=10)
        # plt.ylabel(
        #     "x-axis", fontsize=10)
        # plt.show()

        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        raw_diff = np.ma.array(raw_diff.data, mask=raw_diff.mask)
        # raw_diff = raw_two - raw_one
        def imstats(dat): return (dat.min(), dat.max(), dat.mean(), dat.std())
        dmin, dmax, dmean, dstd = imstats(raw_diff)
        print(dmin, dmax, dmean, dstd)

        # from astropy.stats import sigma_clip, mad_std
        # # SIGMA-CLIP THE DATA
        # masked_diff = sigma_clip(
        #     raw_diff, sigma_lower=10, sigma_upper=10, maxiters=5, cenfunc='median', stdfunc=mad_std)
        # dmin, dmax, dmean, dstd = imstats(masked_diff)

        # print(dmin, dmax, dmean, dstd)
        darkRon = dstd
        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "RON DETECTOR",
            "qc_value": darkRon,
            "qc_comment": "RON in single dark",
            "qc_unit": "electrons",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }, ignore_index=True)

        self.log.debug('completed the ``qc_dark_ron`` method')
        return darkRon

    # use the tab-trigger below for new method
    # xt-class-method
