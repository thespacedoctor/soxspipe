#!/usr/bin/env python
# encoding: utf-8
"""
*Recipe to generate a first approximation of the dispersion solution from single pinhole frames*

:Author:
    David Young & Marco Landoni

:Date Created:
    August 25, 2020
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
import unicodecsv as csv


class soxs_disp_solution(_base_recipe_):
    """
    *generate a first approximation of the dispersion solution from single pinhole frames*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.   

    **Usage**

    ```python
    from soxspipe.recipes import soxs_disp_solution
    disp_map_path = soxs_disp_solution(
        log=log,
        settings=settings,
        inputFrames=sofPath
    ).produce_product()
    ```

    --- 

    ```eval_rst
    .. todo::

        - add a tutorial about ``soxs_disp_solution`` to documentation
    ```
    """

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[]

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_disp_solution, self).__init__(
            log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_disp_solution' object")
        self.settings = settings
        self.recipeSettings = settings["soxs-disp-solution"]
        self.inputFrames = inputFrames

        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_disp_solution - NO MORE, NO LESS.
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
        """*verify input frames match those required by the `soxs_disp_solution` recipe*

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

            if imageTypes[0] != "LAMP,FMTCHK":
                raise TypeError(
                    "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals())

            for i in imageTech:
                if i not in ['ECHELLE,PINHOLE', 'IMAGE']:
                    raise TypeError(
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals())

        else:
            for i in imageTypes:
                if i not in ["LAMP,FMTCHK", "BIAS", "DARK"]:
                    raise TypeError(
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on and a master-bias and possibly a master dark for UVB/VIS" % locals())

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*generate a fisrt guess of the dispersion solution*

        **Return:**
            - ``productPath`` -- the path to the first guess dispersion map
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        # self.inputFrames.summary.pprint_all()

        master_bias = False
        dark = False
        pinhole_image = False

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
        add_filters = {kw("DPR_TYPE"): 'LAMP,FMTCHK', kw("DPR_TECH"): 'IMAGE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        add_filters = {kw("DPR_TYPE"): 'LAMP,FMTCHK',
                       kw("DPR_TECH"): 'ECHELLE,PINHOLE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            pinhole_image = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                         hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        self.pinholeFrame = self.subtract_calibrations(
            inputFrame=pinhole_image, master_bias=master_bias, dark=dark)

        if self.settings["save-intermediate-products"]:
            outDir = self.intermediateRootPath
            filePath = self._write(
                frame=self.pinholeFrame,
                filedir=outDir,
                filename=False,
                overwrite=True
            )
            print(f"\nCalibrated single pinhole frame: {filePath}\n")

        from soxspipe.commonutils import create_dispersion_map
        productPath = create_dispersion_map(
            log=self.log,
            settings=self.settings,
            pinholeFrame=self.pinholeFrame
        ).get()

        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method
