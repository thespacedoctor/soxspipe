#!/usr/bin/env python
# encoding: utf-8
"""
*enhance the wavelength solution achieved with `soxs_disp_solution` by expanding the solution into the spatial dimension (along the slit)*

:Author:
    David Young & Marco Landoni

:Date Created:
    March 17, 2021
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
import ccdproc
from astropy import units as u
from soxspipe.commonutils import keyword_lookup


class soxs_spatial_solution(_base_recipe_):
    """
    *The soxs_spatial_solution recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.   

    See `produce_product` method for usage.

    ```eval_rst
    .. todo::

        - add a tutorial about ``soxs_spatial_solution`` to documentation
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
        super(soxs_spatial_solution, self).__init__(
            log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_spatial_solution' object")
        self.settings = settings
        self.recipeSettings = settings["soxs-spatial-solution"]
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_spatial_solution - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        print("\n# RAW INPUT FRAMES - SUMMARY")
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
        """*verify input frames match those required by the `soxs_spatial_solution` recipe*

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        self._verify_input_frames_basics()

        imageTypes = self.inputFrames.values(
            keyword=kw("DPR_TYPE").lower(), unique=True)
        imageTech = self.inputFrames.values(
            keyword=kw("DPR_TECH").lower(), unique=True)
        imageCat = self.inputFrames.values(
            keyword=kw("DPR_CATG").lower(), unique=True)

        if self.arm == "NIR":
            # WANT ON AND OFF PINHOLE FRAMES
            # MIXED INPUT IMAGE TYPES ARE BAD
            if len(imageTypes) > 1:
                imageTypes = " and ".join(imageTypes)
                print(self.inputFrames.summary)
                raise TypeError(
                    "Input frames are a mix of %(imageTypes)s" % locals())

            if imageTypes[0] != "LAMP,WAVE":
                raise TypeError(
                    "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames for NIR" % locals())

            for i in imageTech:
                if i not in ['ECHELLE,MULTI-PINHOLE', 'IMAGE']:
                    raise TypeError(
                        "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames for NIR" % locals())

        else:
            for i in imageTypes:
                if i not in ["LAMP,WAVE", "BIAS", "DARK"]:
                    raise TypeError(
                        "Input frames for soxspipe spatial_solution need to be LAMP,WAVE and a master-bias and possibly a master dark for UVB/VIS" % locals())

        # LOOK FOR ****
        arm = self.arm
        if arm not in self.supplementaryInput or "DISP_MAP" not in self.supplementaryInput[arm]:
            raise TypeError(
                "Need a DISP_MAP for %(arm)s - none found with the input files" % locals())

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*generate the 2D dispersion map*

        **Return:**
            - ``productPath`` -- the path to the 2D dispersion map

        **Usage**

        ```python
        from soxspipe.recipes import soxs_spatial_solution
        recipe = soxs_spatial_solution(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        disp_map = recipe.produce_product()
        ```
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None

        master_bias = False
        dark = False
        multi_pinhole_image = False

        add_filters = {kw("DPR_CATG"): 'MASTER_BIAS_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS DARK
        add_filters = {kw("DPR_CATG"): 'MASTER_DARK_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # NIR DARK
        add_filters = {kw("DPR_TYPE"): 'LAMP,WAVE',
                       kw("DPR_TECH"): 'IMAGE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # MULTIPINHOLE IMAGE
        add_filters = {kw("DPR_TYPE"): 'LAMP,WAVE',
                       kw("DPR_TECH"): 'ECHELLE,MULTI-PINHOLE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            multi_pinhole_image = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                               hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        self.multiPinholeFrame = self.subtract_calibrations(
            inputFrame=multi_pinhole_image, master_bias=master_bias, dark=dark)

        if self.settings["save-intermediate-products"]:
            fileDir = self.intermediateRootPath
            filepath = self._write(
                self.multiPinholeFrame, fileDir, filename=False, overwrite=True)
            print(f"\nCalibrated multi pinhole frame frame saved to {filepath}\n")

        # GENERATE AN UPDATED DISPERSION MAP
        from soxspipe.commonutils import create_dispersion_map
        mapPath = create_dispersion_map(
            log=self.log,
            settings=self.settings,
            pinholeFrame=self.multiPinholeFrame,
            firstGuessMap=self.supplementaryInput[arm]["DISP_MAP"]
        ).get()

        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return mapPath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
