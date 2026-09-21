#!/usr/bin/env python
"""
*Reduce SOXS/Xshooter data taken in stare mode*

Author
: David Young & Marco Landoni

Date Created
: February 28, 2022
"""

################# GLOBAL IMPORTS ####################
import os
import sys

from soxspipe.commonutils import detector_lookup, subtract_sky
from soxspipe.commonutils.toolkit import (
    generic_quality_checks,
    get_calibrations_path,
    spectroscopic_image_quality_checks,
    utcnow_string,
)

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"


class soxs_stare(base_recipe):
    """
    *Reduce SOXS/Xshooter data taken in stare mode*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- show debug plots. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.



    See `produce_product` method for usage.
    """

    # INITIALISATION

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
        # INHERIT INITIALISATION FROM  BASE_RECIPE
        super().__init__(
            log=log,
            settings=settings,
            inputFrames=inputFrames,
            overwrite=overwrite,
            recipeName="soxs-stare",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        self.log = log
        log.debug("instantiating a new 'soxs_stare' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeSettings = self.get_recipe_settings()
        # xt-self-arg-tmpx

        # INITIAL ACTIONS
        self._collect_input_frames()
        self._verify_and_announce_input_frames()
        self._sort_and_report_input_frames()

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        # GET A TEMPLATE FILENAME USED TO NAME PRODUCTS
        if self.sofName:
            self.filenameTemplate = self.sofName + ".fits"
        else:
            self.filenameTemplate = filenamer(log=self.log, frame=self.objectFrame, settings=self.settings)

        self.generateReponseCurve = False

        return

    def _collect_input_frames(self):
        """*convert the input files to a ccdproc image collection*

        Sets ``self.inputFrames`` and ``self.supplementaryInput``.
        """
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (INPUTFRAMES >
        # IMAGEFILECOLLECTION)
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(log=self.log, settings=self.settings, inputFrames=self.inputFrames)
        self.inputFrames, self.supplementaryInput = sof.get()

        return

    def _verify_and_announce_input_frames(self):
        """*verify the collected frames and report the result to the user*

        Sets ``self.imageType``, through ``verify_input_frames``.
        """
        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_STARE - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        return

    def _sort_and_report_input_frames(self):
        """*sort the image collection by observation date and, when verbose, print it*"""
        # SORT IMAGE COLLECTION
        self.inputFrames.sort(["MJD-OBS"])
        if self.verbose:
            self.log.print("# RAW INPUT FRAMES - SUMMARY")
            self.log.print(self.inputFrames.summary)

        return

    def verify_input_frames(self):
        """*verify the input frame match those required by the soxs_stare recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()
        arm = self.arm

        if self.arm == "NIR":
            error = self._nir_input_frame_error(imageTypes, imageTech, arm)
        else:
            error = self._uvb_vis_input_frame_error(imageTypes, imageCat, arm)

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
        self.log.debug("completed the ``verify_input_frames`` method")
        return

    def _nir_input_frame_error(self, imageTypes, imageTech, arm):
        """*check the NIR input frame types and techniques against those the recipe accepts*

        **Key Arguments:**

        - ``imageTypes`` -- the image types of the input frames
        - ``imageTech`` -- the observing techniques of the input frames
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``error`` -- the error message for the last unaccepted frame, or False when every frame is accepted
        """
        error = False

        if not error:
            okList = [
                "OBJECT",
                "LAMP,FLAT",
                "DARK",
                "STD,FLUX",
                "STD,TELLURIC",
                "OBJECT,ASYNC",
            ]
            if "PAE" in self.settings and self.settings["PAE"]:
                okList.append("FLAT,LAMP")
            for i in imageTypes:
                if i not in okList:
                    error = f"Found a {i} file. Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark (MASTER_DARK_{arm}) or off-frame for NIR."

        if not error:
            for i in imageTech:
                okList = [
                    "ECHELLE,SLIT,STARE",
                    "IMAGE",
                    "ECHELLE,SLIT",
                    "ECHELLE,MULTI-PINHOLE",
                    "ECHELLE,SLIT,NODDING",
                ]
                if "PAE" in self.settings and self.settings["PAE"]:
                    okList.append("ECHELLE,PINHOLE")
                if i not in okList:
                    error = f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark (MASTER_DARK_{arm}) or off-frame for NIR. The sof file is missing a {i} frame."

        return error

    def _uvb_vis_input_frame_error(self, imageTypes, imageCat, arm):
        """*check the UVB/VIS input frame types and catalogues against those the recipe requires*

        **Key Arguments:**

        - ``imageTypes`` -- the image types of the input frames
        - ``imageCat`` -- the product categories of the input frames
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``error`` -- the error message for the last missing or unaccepted frame, or False when every frame is accepted
        """
        error = False

        if not error:
            for i in imageTypes:
                if i not in [
                    "OBJECT",
                    "LAMP,FLAT",
                    "BIAS",
                    "DARK",
                    "STD,FLUX",
                    "STD,TELLURIC",
                    "OBJECT,ASYNC",
                ]:
                    error = f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} frame."

        if not error:
            for i in [f"MASTER_BIAS_{self.arm}", f"DISP_TAB_{self.arm}"]:
                if i not in imageCat:
                    error = f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} frame."

        return error

    def produce_product(self):
        """*The code to generate the product of the soxs_stare recipe*

        **Return:**

        - ``productPath`` -- the path to the final product
        - ``qcTable`` -- the quality control table the recipe reports

        **Usage**

        ```python
        from soxspipe.recipes import soxs_stare
        recipe = soxs_stare(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        productPath, qcTable = recipe.produce_product()
        ```
        """
        self.log.debug("starting the ``produce_product`` method")

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        self.subtractSky = self.recipeSettings["sky-subtraction"]["subtract_sky"]

        allObjectFrames = self._read_stare_object_frames(kw)

        combined_object_notflattened = self.clip_and_stack(
            frames=allObjectFrames,
            recipe="soxs_stare",
            ignore_input_masks=True,
            post_stack_clipping=False,
        )
        self.dateObs = combined_object_notflattened.header[kw("DATE_OBS")]

        master_bias, dark, master_flat = self._read_calibration_frames(kw, arm)
        orderTablePath, dispMap, twoDMap, responseFunctionPath = self._locate_calibration_tables(kw, arm)

        try:
            if not self.recipeSettings["use_flat"]:
                master_flat = False
        except KeyError as e:
            self.log.debug(f"produce_product: `if not self.recipeSettings['use_flat']: mas...` failed, continuing: {e}")
            master_flat = False

        combined_object = self.detrend(
            inputFrame=combined_object_notflattened,
            master_bias=master_bias,
            dark=dark,
            master_flat=master_flat,
            order_table=orderTablePath,
        )

        # INJECT KEYWORDS INTO HEADER
        self.update_fits_keywords(frame=combined_object)

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=combined_object,
            show=False,
            ext=False,
            stdWindow=3,
            title=False,
            surfacePlot=False,
            dispMap=dispMap,
            dispMapImage=twoDMap,
            settings=self.settings,
            skylines=False,
        )

        skymodelCCDData, skySubtractedCCDData, productPath = self._subtract_sky(
            combined_object, twoDMap, dispMap, orderTablePath
        )

        unflattenedSkySubtractedCCDData = self._unflatten_sky_subtracted_frame(
            skySubtractedCCDData, master_flat, combined_object, combined_object_notflattened
        )


        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=combined_object_notflattened, show=True, ext='data', stdWindow=3, title=False, surfacePlot=True)


        from soxspipe.commonutils import horne_extraction

        optimalExtractor = horne_extraction(
            log=self.log,
            skySubtractedFrame=skySubtractedCCDData,
            unflattenedFrame=unflattenedSkySubtractedCCDData,
            subtractedFrame=skymodelCCDData,
            twoDMapPath=twoDMap,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            recipeName=self.recipeName,
            qcTable=self.qc,
            productsTable=self.products,
            dispersionMap=dispMap,
            sofName=self.sofName,
            startNightDate=self.startNightDate,
            debug=self.debug,
            turnOffMP=self.turnOffMP,
        )
        self.qc, self.products, mergedSpectumDF, orderJoins, extractionPath = optimalExtractor.extract()

        # CHECK IF FLUX CALIBRATION IS NEEDED
        filePath_fluxcal = None

        forceFailure = False
        if responseFunctionPath:
            filePath_fluxcal = self._flux_calibrate_spectrum(responseFunctionPath, mergedSpectumDF, combined_object)

        elif self.generateReponseCurve:
            forceFailure = self._generate_response_curve(
                unflattenedSkySubtractedCCDData,
                skymodelCCDData,
                twoDMap,
                dispMap,
                extractionPath,
                orderJoins,
            )

        self._plot_merged_spectrum_qcs(mergedSpectumDF, orderJoins, filePath_fluxcal)

        qcTable = self.report_output()
        self.clean_up(forceFail=forceFailure)

        self.log.debug("completed the ``produce_product`` method")
        return productPath, qcTable

    def _read_stare_object_frames(self, kw):
        """*read the science frames of the first stare input type that has any*

        Flux-standard stare frames set ``self.generateReponseCurve``. PAE mode sets ``self.subtractSky`` to False.

        **Key Arguments:**

        - ``kw`` -- the keyword lookup

        **Return:**

        - ``allObjectFrames`` -- the science frames, as CCDData objects
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        # OBJECT FRAMES
        filter_list = [
            {kw("DPR_TYPE"): "OBJECT", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
            {kw("DPR_TYPE"): "OBJECT,ASYNC", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
        ]
        allObjectFrames = []
        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )
                allObjectFrames.append(singleFrame)

        # FLUX STD FRAMES
        if not len(allObjectFrames):
            add_filters = {
                kw("DPR_TYPE"): "STD,FLUX",
                kw("DPR_TECH"): "ECHELLE,SLIT,STARE",
            }
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )
                allObjectFrames.append(singleFrame)
                self.generateReponseCurve = True

        # FLUX STD FRAMES
        if not len(allObjectFrames):
            add_filters = {
                kw("DPR_TYPE"): "STD,TELLURIC",
                kw("DPR_TECH"): "ECHELLE,SLIT,STARE",
            }
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )
                allObjectFrames.append(singleFrame)

        if not len(allObjectFrames) and "PAE" in self.settings and self.settings["PAE"]:
            add_filters = {
                kw("DPR_TYPE"): "LAMP,FLAT",
                kw("DPR_TECH"): "ECHELLE,PINHOLE",
            }
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )
                allObjectFrames.append(singleFrame)
            self.log.warning("Processing a ORDER-TRACE frame with the stare-mode recipe")
            self.subtractSky = False

        if "PAE" in self.settings and self.settings["PAE"]:
            self.subtractSky = False

        if not len(allObjectFrames):
            add_filters = {
                kw("DPR_TYPE"): "STD,FLUX",
                kw("DPR_TECH"): "ECHELLE,SLIT,NODDING",
            }
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )
                allObjectFrames.append(singleFrame)
            self.log.warning("Processing a NODDING frame with the stare-mode recipe")

        return allObjectFrames

    def _read_calibration_frames(self, kw, arm):
        """*read the master bias, the dark and the master flat*

        **Key Arguments:**

        - ``kw`` -- the keyword lookup
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``master_bias`` -- the master bias, or False when none was supplied
        - ``dark`` -- the master dark, the NIR off-frame or, in PAE mode, the flat-lamp image; False when none
          was supplied
        - ``master_flat`` -- the master flat, or False when none was supplied
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        master_bias = False
        master_flat = False
        dark = False

        add_filters = {kw("PRO_CATG"): "MASTER_BIAS_" + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        # MASTER DARK
        add_filters = {kw("PRO_CATG"): "MASTER_DARK_" + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        if not dark:
            # NIR DARK
            add_filters = {kw("DPR_TYPE"): "OBJECT", kw("DPR_TECH"): "IMAGE"}
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                dark = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        if "PAE" in self.settings and self.settings["PAE"]:
            add_filters = {kw("DPR_TYPE"): "FLAT,LAMP", kw("DPR_TECH"): "IMAGE"}
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                dark = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        # UVB/VIS/NIR FLAT
        add_filters = {kw("PRO_CATG"): "MASTER_FLAT_" + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_flat = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        return master_bias, dark, master_flat

    def _locate_calibration_tables(self, kw, arm):
        """*locate the order table, the dispersion map table and image, and the response table*

        **Key Arguments:**

        - ``kw`` -- the keyword lookup
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``orderTablePath`` -- the path to the order table
        - ``dispMap`` -- the path to the dispersion map table
        - ``twoDMap`` -- the path to the dispersion map image
        - ``responseFunctionPath`` -- the path to the response table, or False when none was supplied
        """
        # FIND THE ORDER TABLE
        filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}"}
        orderTablePath = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE 2D MAP TABLE
        filterDict = {kw("PRO_CATG"): f"DISP_TAB_{arm}"}
        dispMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE 2D MAP IMAGE
        filterDict = {kw("PRO_CATG"): f"DISP_IMAGE_{arm}"}
        twoDMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE RESPONSE FUNCTION, IF PRESENT
        try:
            filterDict = {kw("PRO_CATG"): f"RESP_TAB_{arm}"}
            responseFunctionPath = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        except IndexError as e:
            self.log.debug(f"produce_product: no response function frame for this arm, continuing: {e}")
            responseFunctionPath = False

        return orderTablePath, dispMap, twoDMap, responseFunctionPath

    def _subtract_sky(self, combined_object, twoDMap, dispMap, orderTablePath):
        """*subtract the sky from the detrended frame, when sky subtraction is on, and write the sky products*

        Updates ``self.qc`` and ``self.products``. A sky subtraction that returns no model sets ``self.subtractSky``
        and the recipe setting to False and widens ``horne-extraction-slit-length`` to 30.

        **Key Arguments:**

        - ``combined_object`` -- the detrended object frame
        - ``twoDMap`` -- the path to the dispersion map image
        - ``dispMap`` -- the path to the dispersion map table
        - ``orderTablePath`` -- the path to the order table

        **Return:**

        - ``skymodelCCDData`` -- the sky model, or False when no sky was subtracted
        - ``skySubtractedCCDData`` -- the sky-subtracted frame, or ``combined_object`` when no sky was subtracted
        - ``productPath`` -- the path to the last sky product written, or None when none was written
        """
        productPath = None

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
                recipeName=self.recipeName,
                startNightDate=self.startNightDate,
                debug=self.debug,
            )
            (
                skymodelCCDData,
                skySubtractedCCDData,
                skySubtractedResidualsCCDData,
                self.qc,
                self.products,
            ) = skymodel.subtract()

            if skymodelCCDData is None:
                self.subtractSky = False
                self.recipeSettings["sky-subtraction"]["subtract_sky"] = False
                skymodelCCDData = False
                skySubtractedCCDData = combined_object
                self.recipeSettings["horne-extraction-slit-length"] = 30.0
            else:
                productPath = self._write_sky_products(
                    skySubtractedCCDData, skymodelCCDData, skySubtractedResidualsCCDData
                )

                # ADD QUALITY CHECKS
                self.qc = generic_quality_checks(
                    log=self.log,
                    frame=skySubtractedCCDData,
                    settings=self.settings,
                    recipeName=self.recipeName,
                    qcTable=self.qc,
                )
                self.qc = spectroscopic_image_quality_checks(
                    log=self.log,
                    frame=skySubtractedCCDData,
                    settings=self.settings,
                    recipeName=self.recipeName,
                    qcTable=self.qc,
                    orderTablePath=orderTablePath,
                )
        else:
            skymodelCCDData = False
            skySubtractedCCDData = combined_object

        return skymodelCCDData, skySubtractedCCDData, productPath

    def _write_sky_products(self, skySubtractedCCDData, skymodelCCDData, skySubtractedResidualsCCDData):
        """*write the sky-subtracted frame, the sky model and the residuals, and record each as a product*

        All three product rows share one reduction timestamp. Updates ``self.products``.

        **Key Arguments:**

        - ``skySubtractedCCDData`` -- the sky-subtracted frame
        - ``skymodelCCDData`` -- the sky model
        - ``skySubtractedResidualsCCDData`` -- the sky subtraction residuals

        **Return:**

        - ``productPath`` -- the path to the residuals frame, the last product written
        """
        # WRITE SKY-SUBTRACTON TO DISK
        filename = self.filenameTemplate.replace(".fits", "_SKYSUB.fits")
        productPath = self._write(
            frame=skySubtractedCCDData,
            filedir=self.workspaceRootPath,
            filename=filename,
            overwrite=True,
            maskToZero=True,
        )
        filename = os.path.basename(productPath)
        utcnow = utcnow_string()
        self.add_product(
            recipeName="soxs-stare",
            productLabel="SKY_SUBTRACTED_OBJECT",
            fileName=filename,
            filePath=productPath,
            productDesc="The sky-subtracted object",
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        # WRITE SKY-MODEL TO DISK
        filename = self.filenameTemplate.replace(".fits", "_SKYMODEL.fits")
        productPath = self._write(
            frame=skymodelCCDData,
            filedir=self.workspaceRootPath,
            filename=filename,
            overwrite=True,
        )
        filename = os.path.basename(productPath)
        self.add_product(
            recipeName="soxs-stare",
            productLabel="SKY_MODEL",
            fileName=filename,
            filePath=productPath,
            productDesc="The sky background model",
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        if True:
            # WRITE SKY-MODEL TO DISK
            filename = self.filenameTemplate.replace(".fits", "_SKYSUB_RESIDUALS.fits")
            productPath = self._write(
                frame=skySubtractedResidualsCCDData,
                filedir=self.workspaceRootPath,
                filename=filename,
                overwrite=True,
            )
            filename = os.path.basename(productPath)
            self.add_product(
                recipeName="soxs-stare",
                productLabel="SKY_SUB_RESIDUALS",
                fileName=filename,
                filePath=productPath,
                productDesc="The sky subtraction residuals",
                reductionDateUtc=utcnow,
                fileType="FITS",
                label="PROD",
            )

        return productPath

    def _unflatten_sky_subtracted_frame(
        self, skySubtractedCCDData, master_flat, combined_object, combined_object_notflattened
    ):
        """*return the frame the optimal extraction uses as its unflattened input*

        **Key Arguments:**

        - ``skySubtractedCCDData`` -- the sky-subtracted frame
        - ``master_flat`` -- the master flat, or False when none is used
        - ``combined_object`` -- the detrended object frame, whose uncertainty the unflattened frame takes
        - ``combined_object_notflattened`` -- the stacked frame before detrending

        **Return:**

        - ``unflattenedSkySubtractedCCDData`` -- the sky-subtracted frame multiplied back by the master flat, the
          sky-subtracted frame itself when there is no flat, or the stacked frame when no sky was subtracted
        """
        if self.subtractSky:
            if master_flat:
                unflattenedSkySubtractedCCDData = skySubtractedCCDData.multiply(master_flat)
                unflattenedSkySubtractedCCDData.uncertainty = combined_object.uncertainty.array
                unflattenedSkySubtractedCCDData.header = skySubtractedCCDData.header
            else:
                unflattenedSkySubtractedCCDData = skySubtractedCCDData
        else:
            unflattenedSkySubtractedCCDData = combined_object_notflattened
            unflattenedSkymodelCCDData = False

        return unflattenedSkySubtractedCCDData

    def _flux_calibrate_spectrum(self, responseFunctionPath, mergedSpectumDF, combined_object):
        """*flux calibrate the merged spectrum with the supplied response function*

        Updates ``self.products``.

        **Key Arguments:**

        - ``responseFunctionPath`` -- the path to the response table
        - ``mergedSpectumDF`` -- the merged extracted spectrum
        - ``combined_object`` -- the detrended object frame, whose header supplies the airmass and exposure time

        **Return:**

        - ``filePath_fluxcal`` -- the path to the flux-calibrated spectrum
        """
        import pandas as pd

        from soxspipe.commonutils import flux_calibration

        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)

        detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

        self.log.print("# FLUX CALIBRATING THE SPECTRUM\n")
        fluxCalibrator = flux_calibration(
            log=self.log,
            responseFunction=responseFunctionPath,
            extractedSpectrum=mergedSpectumDF,
            settings=self.settings,
            airmass=combined_object.header.get("HIERARCH ESO TEL AIRM END"),
            exptime=combined_object.header.get("EXPTIME"),
            extinctionPath=calibrationRootPath + "/" + detectorParams["extinction"],
            arm=self.arm,
            header=combined_object.header,
            recipeName=self.recipeName,
            startNightDate=self.startNightDate,
            sofName=self.sofName,
            debug=self.debug,
        )
        filePath_fluxcal, prod = fluxCalibrator.calibrate()
        self.products = pd.concat([self.products, prod], ignore_index=True)
        # self.qc, self.products, calibratedSpectrumDF, calibrationPath = fluxCalibrator.calibrate()
        self.log.print("# FLUX CALIBRATION COMPLETED\n")

        return filePath_fluxcal

    def _generate_response_curve(
        self,
        unflattenedSkySubtractedCCDData,
        skymodelCCDData,
        twoDMap,
        dispMap,
        extractionPath,
        orderJoins,
    ):
        """*extract the unflattened standard-star spectrum and compute the response function from it*

        Updates ``self.qc`` and ``self.products``.

        **Key Arguments:**

        - ``unflattenedSkySubtractedCCDData`` -- the unflattened frame to extract
        - ``skymodelCCDData`` -- the sky model, or False when no sky was subtracted
        - ``twoDMap`` -- the path to the dispersion map image
        - ``dispMap`` -- the path to the dispersion map table
        - ``extractionPath`` -- the path to the flat-fielded extraction of the standard
        - ``orderJoins`` -- the order joins of the flat-fielded extraction

        **Return:**

        - ``forceFailure`` -- the failure the response function reports, which forces the recipe to fail
        """
        from soxspipe.commonutils import horne_extraction

        optimalExtractor = horne_extraction(
            log=self.log,
            skySubtractedFrame=unflattenedSkySubtractedCCDData,
            unflattenedFrame=unflattenedSkySubtractedCCDData,
            subtractedFrame=skymodelCCDData,
            twoDMapPath=twoDMap,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            recipeName=self.recipeName,
            qcTable=self.qc,
            productsTable=self.products,
            dispersionMap=dispMap,
            sofName=self.sofName,
            startNightDate=self.startNightDate,
            debug=self.debug,
            notFlattened=True,
            turnOffMP=self.turnOffMP,
        )
        self.qc, self.products, _, _, extractionPath_notflat = optimalExtractor.extract()

        # GETTING THE RESPONSE
        from soxspipe.commonutils import response_function

        self.log.print("# CALCULATING RESPONSE FUNCTION\n")
        response = response_function(
            log=self.log,
            settings=self.settings,
            recipeName=self.recipeName,
            sofName=self.sofName,
            stdExtractionPath=extractionPath,
            qcTable=self.qc,
            productsTable=self.products,
            startNightDate=self.startNightDate,
            stdNotFlatExtractionPath=extractionPath_notflat,
            orderJoins=orderJoins,
        )
        self.qc, self.products, forceFailure = response.get()

        return forceFailure

    def _plot_merged_spectrum_qcs(self, mergedSpectumDF, orderJoins, filePath_fluxcal):
        """*plot the merged spectrum QC and, when the spectrum was flux calibrated, the calibrated one too*

        Updates ``self.products``.

        **Key Arguments:**

        - ``mergedSpectumDF`` -- the merged extracted spectrum
        - ``orderJoins`` -- the order joins of the extraction
        - ``filePath_fluxcal`` -- the path to the flux-calibrated spectrum, or None when it was not calibrated
        """
        from soxspipe.commonutils.toolkit import plot_merged_spectrum_qc

        self.products, filePath = plot_merged_spectrum_qc(
            merged_orders=mergedSpectumDF,
            products=self.products,
            log=self.log,
            qcDir=self.qcDir,
            filenameTemplate=self.filenameTemplate,
            noddingSequence=None,
            dateObs=self.dateObs,
            arm=self.arm,
            recipeName=self.recipeName,
            orderJoins=orderJoins,
            debug=self.debug,
            fluxCalibrated=False,
            qcTable=self.qc,
            settings=self.settings,
        )

        if filePath_fluxcal:
            from astropy import units as u
            from astropy.table import Table

            fluxcal_spec = Table.read(filePath_fluxcal, format="fits")
            fluxcal_spec["WAVE"] = fluxcal_spec["WAVE"] * u.nm
            fluxcal_spec["FLUX_COUNTS"] = fluxcal_spec["FLUX_CALIBRATED"]  # BACK COMPATIBILITY WITH THE CODE
            fluxcal_spec["SNR"] = mergedSpectumDF["SNR"]
            fluxcal_spec["SKY_COUNTS"] = mergedSpectumDF["SKY_COUNTS"]
            self.products, filePath = plot_merged_spectrum_qc(
                merged_orders=fluxcal_spec,
                products=self.products,
                log=self.log,
                qcDir=self.qcDir,
                filenameTemplate=self.filenameTemplate,
                noddingSequence=None,
                dateObs=self.dateObs,
                arm=self.arm,
                recipeName=self.recipeName,
                orderJoins=orderJoins,
                debug=self.debug,
                fluxCalibrated=True,
                qcTable=self.qc,
                settings=self.settings,
            )

        return

    # USE THE TAB-TRIGGER BELOW FOR NEW METHOD
    # xt-class-method

    # OVERRIDE METHOD ATTRIBUTES
    # method-override-tmpx
