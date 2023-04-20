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
from ._base_recipe_ import _base_recipe_
from soxspipe.commonutils import subtract_sky
from soxspipe.commonutils.toolkit import generic_quality_checks, spectroscopic_image_quality_checks
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
        - ``verbose`` -- verbose. True or False. Default *False*
        - ``overwrite`` -- overwrite the prodcut file if it already exists. Default *False*


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
            verbose=False,
            overwrite=False

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_stare, self).__init__(
            log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-stare")
        self.log = log
        log.debug("instansiating a new 'soxs_stare' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeSettings = settings[self.recipeName]
        # xt-self-arg-tmpx

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        from soxspipe.commonutils.set_of_files import set_of_files
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
        self.inputFrames.sort(['MJD-OBS'])
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

        error = False

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if self.arm == "NIR":
            if not error:
                for i in imageTypes:
                    if i not in ["OBJECT", "LAMP,FLAT", "DARK", "STD,FLUX"]:
                        error = "Found a {i} file. Input frames for soxspipe stare need to be OBJECT, ***. Can optionally supply a master-flat for NIR."

            if not error:
                for i in imageTech:
                    if i not in ['ECHELLE,SLIT,STARE', "IMAGE", "ECHELLE,SLIT", "ECHELLE,MULTI-PINHOLE", "ECHELLE,SLIT,NODDING"]:
                        error = "Input frames for soxspipe stare need to be ********* lamp on and lamp off frames for NIR" % locals()

        else:
            if not error:
                for i in imageTypes:
                    if i not in ["OBJECT", "LAMP,FLAT", "BIAS", "DARK"]:
                        error = "Input frames for soxspipe stare need to be ********* and a master-bias and possibly a master dark for UVB/VIS" % locals()

            if not error:
                for i in [f"MASTER_BIAS_{self.arm}", f"DISP_TAB_{self.arm}"]:
                    if i not in imageCat:
                        error = "Input frames for soxspipe stare need to be ********* and a master-bias and possibly a master dark for UVB/VIS" % locals()

        # LOOK FOR ****
        arm = self.arm
        # if arm not in self.supplementaryInput or "DISP_MAP" not in self.supplementaryInput[arm]:
        #     raise TypeError(
        #         "Need a **** for %(arm)s - none found with the input files" % locals())

        if error:
            sys.stdout.write("\x1b[1A\x1b[2K")
            print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            print(self.inputFrames.summary)
            print()
            raise TypeError(error)

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

        from astropy.nddata import CCDData
        from astropy import units as u

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None
        master_bias = False
        dark = False

        # OBJECT FRAMES
        add_filters = {kw("DPR_TYPE"): 'OBJECT',
                       kw("DPR_TECH"): 'ECHELLE,SLIT,STARE'}
        allObjectFrames = []
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            singleFrame = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
            allObjectFrames.append(singleFrame)

        # FLUX STD FRAMES
        if not len(allObjectFrames):
            add_filters = {kw("DPR_TYPE"): 'STD,FLUX',
                           kw("DPR_TECH"): 'ECHELLE,SLIT,STARE'}
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)

        if not len(allObjectFrames):
            if not len(allObjectFrames):
                add_filters = {kw("DPR_TYPE"): 'STD,FLUX',
                               kw("DPR_TECH"): 'ECHELLE,SLIT,NODDING'}
                allObjectFrames = []
                for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                    singleFrame = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                               hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                    allObjectFrames.append(singleFrame)
                self.log.warning("Processing a NODDING frame with the stare-mode recipe")

        combined_object = self.clip_and_stack(
            frames=allObjectFrames, recipe="soxs_stare", ignore_input_masks=True, post_stack_clipping=False)
        self.dateObs = combined_object.header[kw("DATE_OBS")]

        add_filters = {kw("PRO_CATG"): 'MASTER_BIAS_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS DARK
        add_filters = {kw("PRO_CATG"): 'MASTER_DARK_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # NIR DARK
        add_filters = {kw("DPR_TYPE"): 'DARK',
                       kw("DPR_TECH"): 'IMAGE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS/NIR FLAT
        add_filters = {kw("PRO_CATG"): 'MASTER_FLAT_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_flat = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # FIND THE ORDER TABLE
        filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}"}
        orderTablePath = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE 2D MAP TABLE
        filterDict = {kw("PRO_CATG"): f"DISP_TAB_{arm}"}
        dispMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE 2D MAP IMAGE
        filterDict = {kw("PRO_CATG"): f"DISP_IMAGE_{arm}"}
        twoDMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        if False:
            combined_object = self.detrend(
                inputFrame=combined_object, master_bias=master_bias, dark=dark, master_flat=master_flat, order_table=orderTablePath)
        else:
            combined_object = self.detrend(
                inputFrame=combined_object, master_bias=master_bias, dark=dark, master_flat=master_flat)

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=combined_object, show=False, ext=False, stdWindow=1, title=False, surfacePlot=True, dispMap=dispMap, dispMapImage=twoDMap, settings=self.settings, skylines=True)

        skymodel = subtract_sky(
            log=self.log,
            settings=self.settings,
            objectFrame=combined_object,
            twoDMap=twoDMap,
            qcTable=self.qc,
            productsTable=self.products,
            dispMap=dispMap,
            sofName=self.sofName,
            recipeName=self.recipeName
        )
        skymodelCCDData, skySubtractedCCDData, self.qc, self.products = skymodel.subtract()

        from os.path import expanduser
        home = expanduser("~")
        fileDir = self.settings["workspace-root-dir"].replace("~", home)
        # filename = "/override/filename.fits"
        filepath = self._write(combined_object, fileDir, filename="combined_object_frame.fits", overwrite=True)
        # print(f"\nxxx frame saved to {filepath}\n")

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log, frame=skySubtractedCCDData, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)
        self.qc = spectroscopic_image_quality_checks(
            log=self.log, frame=skySubtractedCCDData, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc, orderTablePath=orderTablePath)

        self.clean_up()
        self.report_output()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
