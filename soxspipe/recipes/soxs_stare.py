#!/usr/bin/env python
# encoding: utf-8
"""
*Reduce SOXS data taken in stare mode*

:Author:
    David Young & Marco Landoni

:Date Created:
    February 28, 2022
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
from astropy import units as u
import ccdproc
from astropy.nddata import CCDData
import numpy as np
from ._base_recipe_ import _base_recipe_
from soxspipe.commonutils import set_of_files
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class soxs_stare(_base_recipe_):
    """
    *The soxs_stare recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.   


    See `produce_product` method for usage.

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``soxs_stare`` to documentation
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
        super(soxs_stare, self).__init__(
            log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_stare' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-stare"
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_stare - NO MORE, NO LESS.
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
        """*verify the input frame match those required by the soxs_stare recipe*

        **Return:**
            - ``None``

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if self.arm == "NIR":

            for i in imageTypes:
                if i not in ["OBJECT", "LAMP,FLAT", "DARK"]:
                    raise TypeError(
                        f"Found a {i} file. Input frames for soxspipe stare need to be OBJECT, ***. Can optionally supply a master-flat for NIR.")

            for i in imageTech:
                if i not in ['ECHELLE,SLIT,STARE', "IMAGE", "ECHELLE,SLIT"]:
                    raise TypeError(
                        "Input frames for soxspipe stare need to be ********* lamp on and lamp off frames for NIR" % locals())

        else:
            for i in imageTypes:
                if i not in ["LAMP,FMTCHK", "BIAS", "DARK"]:
                    raise TypeError(
                        "Input frames for soxspipe stare need to be ********* and a master-bias and possibly a master dark for UVB/VIS" % locals())

        # LOOK FOR ****
        arm = self.arm
        # if arm not in self.supplementaryInput or "DISP_MAP" not in self.supplementaryInput[arm]:
        #     raise TypeError(
        #         "Need a **** for %(arm)s - none found with the input files" % locals())

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*The code to generate the product of the soxs_stare recipe*

        **Return:**
            - ``productPath`` -- the path to the final product

        **Usage**

        ```python
        from soxspipe.recipes import soxs_stare
        recipe = soxs_stare(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        stareFrame = recipe.produce_product()
        ```
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None

        # OBJECT FRAMES
        add_filters = {kw("DPR_TYPE"): 'OBJECT',
                       kw("DPR_TECH"): 'ECHELLE,SLIT,STARE'}
        allObjectFrames = []
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            singleFrame = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
            allObjectFrames.append(singleFrame)

        combined_object = self.clip_and_stack(
            frames=allObjectFrames, recipe="soxs_stare")
        self.dateObs = combined_object.header[kw("DATE_OBS")]

        # NIR DARK
        add_filters = {kw("DPR_TYPE"): 'DARK',
                       kw("DPR_TECH"): 'IMAGE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS/NIR FLAT
        add_filters = {kw("DPR_CATG"): 'MASTER_LAMP-FLAT_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_flat = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        combined_object = self.detrend(
            inputFrame=combined_object, master_bias=False, dark=dark, master_flat=master_flat, order_table=False)

        from os.path import expanduser
        home = expanduser("~")
        fileDir = self.settings["intermediate-data-root"].replace("~", home)
        # filename = "/override/filename.fits"
        filepath = self._write(combined_object, fileDir, filename="tmp.fits", overwrite=True)
        print(f"\nxxx frame saved to {filepath}\n")

        sys.exit(0)

        # OBJECT
        # add_filters = {kw("DPR_CATG"): 'OBJECT'}
        # for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
        #     master_flat = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
        #                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # xsoxs-append-to-product-report-table

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log, frame=mflat, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)
        self.qc = spectroscopic_image_quality_checks(
            log=self.log, frame=mflat, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc, orderTablePath=orderTablePath)

        self.clean_up()
        self.report_output()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
