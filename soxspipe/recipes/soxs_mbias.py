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
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils import set_of_files
from ._base_recipe_ import _base_recipe_
import numpy as np
from astropy.nddata import CCDData
from astropy import units as u
import ccdproc
from soxspipe.commonutils import keyword_lookup


class soxs_mbias(_base_recipe_):
    """
    *The* `soxs_mbias` *recipe is used to generate a master-bias frame from a set of input raw bias frames. The recipe is used only for the UV-VIS arm as NIR frames have bias (and dark current) removed by subtracting an off-frame of equal expsoure length.*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.

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
            inputFrames=[]

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_mbias, self).__init__(log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mbias' object")
        self.settings = settings
        self.recipeSettings = settings["soxs-mbias"]
        self.inputFrames = inputFrames
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

        print("\n# RAW INPUT BIAS FRAMES - SUMMARY")
        # SORT IMAGE COLLECTION
        self.inputFrames.sort(['mjd-obs'])
        print(self.inputFrames.summary, "\n")

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

        # WRITE TO DISK
        productPath = self._write(
            frame=combined_bias_mean,
            filedir=self.intermediateRootPath,
            filename=False,
            overwrite=True
        )

        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method
