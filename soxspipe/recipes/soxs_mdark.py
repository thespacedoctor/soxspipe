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
from soxspipe.commonutils.toolkit import generic_quality_checks
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
            inputFrames=self.inputFrames,
            ext=self.settings['data-extension']
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
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

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

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log, frame=combined_dark_mean, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)

        medianFlux = self.qc_median_flux_level(
            frame=combined_dark_mean,
            frameType="MDARK",
            frameName="master dark"
        )

        self.qc_ron(
            frameType="DARK"
        )

        # WRITE TO DISK
        productPath = self._write(
            frame=combined_dark_mean,
            filedir=self.intermediateRootPath,
            filename=False,
            overwrite=True
        )
        filename = os.path.basename(productPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.products = self.products.append({
            "soxspipe_recipe": self.recipeName,
            "product_label": "MDARK",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} Master dark frame",
            "file_path": productPath
        }, ignore_index=True)

        self.report_output()
        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method
