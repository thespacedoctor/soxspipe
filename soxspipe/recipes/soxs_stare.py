#!/usr/bin/env python
# encoding: utf-8
"""
*Reduce SOXS/Xshooter data taken in stare mode*

Author
: David Young & Marco Landoni

Date Created
: February 28, 2022
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
from .base_recipe import base_recipe
from soxspipe.commonutils import subtract_sky
from soxspipe.commonutils.toolkit import generic_quality_checks, spectroscopic_image_quality_checks
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class soxs_stare(base_recipe):
    """
    *Reduce SOXS/Xshooter data taken in stare mode*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*


    See `produce_product` method for usage.
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
        # INHERIT INITIALISATION FROM  base_recipe
        super(soxs_stare, self).__init__(
            log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-stare")
        self.log = log
        log.debug("instantiating a new 'soxs_stare' object")
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
            inputFrames=self.inputFrames
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_stare - NO MORE, NO LESS.
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

        # GET A TEMPLATE FILENAME USED TO NAME PRODUCTS
        if self.sofName:
            self.filenameTemplate = self.sofName + ".fits"
        else:
            self.filenameTemplate = filenamer(
                log=self.log,
                frame=self.objectFrame,
                settings=self.settings
            )

        return None

    def verify_input_frames(
            self):
        """*verify the input frame match those required by the soxs_stare recipe*

        **Return:**

        - ``None``

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        error = False

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()
        arm = self.arm

        if self.arm == "NIR":
            if not error:
                okList = ["OBJECT", "LAMP,FLAT", "DARK", "STD,FLUX", "STD,TELLURIC", "OBJECT,ASYNC"]
                if "PAE" in self.settings and self.settings["PAE"]:
                    okList.append("FLAT,LAMP")
                for i in imageTypes:
                    if i not in okList:
                        print("1")
                        error = f"Found a {i} file. Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark (MASTER_DARK_{arm}) or off-frame for NIR."

            if not error:
                for i in imageTech:
                    okList = ['ECHELLE,SLIT,STARE', "IMAGE", "ECHELLE,SLIT", "ECHELLE,MULTI-PINHOLE", "ECHELLE,SLIT,NODDING"]
                    if "PAE" in self.settings and self.settings["PAE"]:
                        okList.append("ECHELLE,PINHOLE")
                    if i not in okList:
                        print("2")
                        error = f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark (MASTER_DARK_{arm}) or off-frame for NIR. The sof file is missing a {i} frame."

        else:
            if not error:
                for i in imageTypes:
                    print(i)
                    if i not in ["OBJECT", "LAMP,FLAT", "BIAS", "DARK", "STD,FLUX", "STD,TELLURIC", "OBJECT,ASYNC"]:
                        print("11")
                        error = f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} frame."

            if not error:
                for i in [f"MASTER_BIAS_{self.arm}", f"DISP_TAB_{self.arm}"]:
                    if i not in imageCat:
                        print("22")
                        error = f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} frame."

        # if arm not in self.supplementaryInput or "DISP_MAP" not in self.supplementaryInput[arm]:
        #     raise TypeError(
        #         "Need a **** for %(arm)s - none found with the input files" % locals())

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
        import pandas as pd
        from datetime import datetime

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None
        master_bias = False
        dark = False

        self.subtractSky = True

        # OBJECT FRAMES
        filter_list = [
            {kw("DPR_TYPE"): 'OBJECT', kw("DPR_TECH"): 'ECHELLE,SLIT,STARE'},
            {kw("DPR_TYPE"): 'OBJECT,ASYNC', kw("DPR_TECH"): 'ECHELLE,SLIT,STARE'}
        ]
        allObjectFrames = []
        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)

        # FLUX STD FRAMES
        if not len(allObjectFrames):
            add_filters = {kw("DPR_TYPE"): 'STD,FLUX',
                           kw("DPR_TECH"): 'ECHELLE,SLIT,STARE'}
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)

        # FLUX STD FRAMES
        if not len(allObjectFrames):
            add_filters = {kw("DPR_TYPE"): 'STD,TELLURIC',
                           kw("DPR_TECH"): 'ECHELLE,SLIT,STARE'}
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)

        if not len(allObjectFrames) and "PAE" in self.settings and self.settings["PAE"]:
            add_filters = {kw("DPR_TYPE"): 'LAMP,FLAT',
                           kw("DPR_TECH"): 'ECHELLE,PINHOLE'}
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)
            self.log.warning("Processing a ORDER-TRACE frame with the stare-mode recipe")
            self.subtractSky = False

        if "PAE" in self.settings and self.settings["PAE"]:
            self.subtractSky = False

        if not len(allObjectFrames):
            add_filters = {kw("DPR_TYPE"): 'STD,FLUX',
                           kw("DPR_TECH"): 'ECHELLE,SLIT,NODDING'}
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)
            self.log.warning("Processing a NODDING frame with the stare-mode recipe")

        combined_object = self.clip_and_stack(
            frames=allObjectFrames, recipe="soxs_stare", ignore_input_masks=True, post_stack_clipping=False)
        self.dateObs = combined_object.header[kw("DATE_OBS")]

        add_filters = {kw("PRO_CATG"): 'MASTER_BIAS_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # MASTER DARK
        add_filters = {kw("PRO_CATG"): 'MASTER_DARK_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        if not dark:
            # NIR DARK
            add_filters = {kw("DPR_TYPE"): 'OBJECT',
                           kw("DPR_TECH"): 'IMAGE'}
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                dark = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                    hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        if "PAE" in self.settings and self.settings["PAE"]:
            add_filters = {kw("DPR_TYPE"): 'FLAT,LAMP',
                           kw("DPR_TECH"): 'IMAGE'}
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                dark = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                    hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS/NIR FLAT
        add_filters = {kw("PRO_CATG"): 'MASTER_FLAT_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_flat = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
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
        try:
            if not self.recipeSettings["use_flat"]:
                master_flat = False
        except:
            master_flat = False

        combined_object = self.detrend(
            inputFrame=combined_object, master_bias=master_bias, dark=dark, master_flat=master_flat, order_table=orderTablePath)

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=combined_object, show=False, ext=False, stdWindow=3, title=False, surfacePlot=False, dispMap=dispMap, dispMapImage=twoDMap, settings=self.settings, skylines=False)

        if self.subtractSky:

            skymodel = subtract_sky(
                log=self.log,
                settings=self.settings,
                recipeSettings=self.recipeSettings,
                objectFrame=combined_object,
                twoDMap=twoDMap,
                qcTable=self.qc,
                productsTable=self.products,
                dispMap=dispMap,
                sofName=self.sofName,
                recipeName=self.recipeName
            )
            skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData, self.qc, self.products = skymodel.subtract()

            # WRITE SKY-SUBTRACTON TO DISK
            filename = self.filenameTemplate.replace(".fits", "_SKYSUB.fits")
            productPath = self._write(
                frame=skySubtractedCCDData,
                filedir=self.workspaceRootPath,
                filename=filename,
                overwrite=True,
                maskToZero=True
            )
            filename = os.path.basename(productPath)
            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": "soxs-stare",
                "product_label": "SKY_SUBTRACTED_OBJECT",
                "file_name": filename,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"The sky-subtracted object",
                "file_path": productPath,
                "label": "PROD"
            }).to_frame().T], ignore_index=True)

            # WRITE SKY-MODEL TO DISK
            filename = self.filenameTemplate.replace(".fits", "_SKYMODEL.fits")
            productPath = self._write(
                frame=skymodelCCDData,
                filedir=self.workspaceRootPath,
                filename=filename,
                overwrite=True
            )
            filename = os.path.basename(productPath)
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": "soxs-stare",
                "product_label": "SKY_MODEL",
                "file_name": filename,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"The sky background model",
                "file_path": productPath,
                "label": "PROD"
            }).to_frame().T], ignore_index=True)

            if True:
                # WRITE SKY-MODEL TO DISK
                filename = self.filenameTemplate.replace(".fits", "_SKYSUB_RESIDUALS.fits")
                productPath = self._write(
                    frame=skySubtractedResidualsCCDData,
                    filedir=self.workspaceRootPath,
                    filename=filename,
                    overwrite=True
                )
                filename = os.path.basename(productPath)
                self.products = pd.concat([self.products, pd.Series({
                    "soxspipe_recipe": "soxs-stare",
                    "product_label": "SKY_SUB_RESIDUALS",
                    "file_name": filename,
                    "file_type": "FITS",
                    "obs_date_utc": self.dateObs,
                    "reduction_date_utc": utcnow,
                    "product_desc": f"The sky subtraction residuals",
                    "file_path": productPath,
                    "label": "PROD"
                }).to_frame().T], ignore_index=True)

            # ADD QUALITY CHECKS
            self.qc = generic_quality_checks(
                log=self.log, frame=skySubtractedCCDData, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)
            self.qc = spectroscopic_image_quality_checks(
                log=self.log, frame=skySubtractedCCDData, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc, orderTablePath=orderTablePath)
        else:
            skymodelCCDData = False
            skySubtractedCCDData = combined_object

        from soxspipe.commonutils import horne_extraction
        optimalExtractor = horne_extraction(
            log=self.log,
            skyModelFrame=skymodelCCDData,
            skySubtractedFrame=skySubtractedCCDData,
            twoDMapPath=twoDMap,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            recipeName=self.recipeName,
            qcTable=self.qc,
            productsTable=self.products,
            dispersionMap=dispMap,
            sofName=self.sofName
        )

        self.qc, self.products, mergedSpectumDF = optimalExtractor.extract()

        self.clean_up()
        self.report_output()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
