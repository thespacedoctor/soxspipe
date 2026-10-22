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
import re
import sys
from contextlib import closing

from soxspipe.commonutils import detector_lookup, filenamer, subtract_sky
from soxspipe.commonutils.slit_illumination_qc import read_rotator_angle
from soxspipe.commonutils.toolkit import (
    generic_quality_checks,
    get_calibrations_path,
    spectroscopic_image_quality_checks,
    utcnow_string,
)

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"

# THE SLIT OFFSETS FOLLOW THE ROTATOR ANGLE, SO INPUT FRAMES MORE THAN THIS FAR APART IN ABSROT (DEGREES) ARE FLAGGED
ABSROT_SPREAD_WARNING_DEGREES = 10.0

# THE FAILURE OF A STANDARD STAR WITH NO OBJECT TRACE, WHERE ONLY A BOXCAR AT THE SLIT CENTRE COULD BE EXTRACTED
NO_TRACE_STANDARD_MESSAGE = (
    "No object trace was found in the standard-star frame, so no response curve was made. "
    "The spectrum is a robust boxcar extraction; check the pointing and conditions."
)


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
            # NO SET-OF-FILES NAME TO FALL BACK ON, SO NAME PRODUCTS FROM THE
            # OBJECT/STANDARD FRAME INSTEAD (DY-111 REVIEW FOLLOW-UP)
            self.filenameTemplate = self._resolve_filename_template_frame()

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

    def _resolve_filename_template_frame(self):
        """*read the header of the first science-type prepared frame and derive the product filename template*

        Filters the prepared frames with the same OBJECT / OBJECT,ASYNC / STD,FLUX / STD,TELLURIC idiom
        `_read_stare_object_frames` uses to pick out the science frames, so a calibration frame that
        ``verify_input_frames`` also accepts (bias, dark, flat, master calibrations) never drives the
        product name. Falls back to the first prepared frame only when no science-type frame is present.

        **Return:**

        - ``filenameTemplate`` -- the filename template `filenamer` derives from the selected frame's header
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        kw = self.kw
        scienceFilters = [
            {kw("DPR_TYPE"): "OBJECT", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
            {kw("DPR_TYPE"): "OBJECT,ASYNC", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
            {kw("DPR_TYPE"): "STD,FLUX", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
            {kw("DPR_TYPE"): "STD,TELLURIC", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
        ]

        framePath = None
        for filters in scienceFilters:
            matchingPaths = self.inputFrames.files_filtered(include_path=True, **filters)
            if len(matchingPaths):
                framePath = matchingPaths[0]
                break

        if framePath is None:
            # NO SCIENCE-TYPE FRAME PRESENT (E.G. A CALIBRATION-ONLY SET OF
            # FILES) -- FALL BACK TO THE FIRST PREPARED FRAME
            self.log.debug(
                "_resolve_filename_template_frame: no OBJECT/standard frame found, naming from the first "
                "prepared frame instead"
            )
            framePath = self.inputFrames.files_filtered(include_path=True)[0]

        # HEADER-ONLY READ -- `filenamer` NEEDS THE HEADER (AND WCS) ONLY, SO
        # SKIP LOADING THE UNCERTAINTY, MASK AND FLAGS EXTENSIONS
        frame = CCDData.read(framePath, hdu=0, unit=u.electron)

        return filenamer(log=self.log, frame=frame, settings=self.settings)

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
                    error = (
                        f"Found a {i} file. Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), "
                        f"a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an "
                        f"order-location table (ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark "
                        f"(MASTER_DARK_{arm}) or off-frame for NIR."
                    )

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
                    error = (
                        f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map "
                        f"image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table "
                        f"(ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark (MASTER_DARK_{arm}) or "
                        f"off-frame for NIR. The sof file is missing a {i} frame."
                    )

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
                    error = (
                        f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map "
                        f"image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table "
                        f"(ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and "
                        f"optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} "
                        f"frame."
                    )

        if not error:
            for i in [f"MASTER_BIAS_{self.arm}", f"DISP_TAB_{self.arm}"]:
                if i not in imageCat:
                    error = (
                        f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), a dispersion map "
                        f"image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table "
                        f"(ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and "
                        f"optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} "
                        f"frame."
                    )

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

        kw = self.kw

        self.subtractSky = self.recipeSettings["sky-subtraction"]["subtract_sky"]

        allObjectFrames, allObjectNames = self._read_stare_object_frames(kw)
        # THE ATTACHED-ARC TEST MEASURES THE ARC FROM THE FIRST START AND LAST END OF THE BLOCK (DY-1666)
        self.scienceHeaders = [frame.header.copy() for frame in allObjectFrames]

        extractEach = self._should_extract_each_frame(allObjectFrames)
        alsoExtractEachStandard = self._also_extract_each_standard_frame(allObjectFrames)
        self._delete_stale_qc_rows(extractEach, alsoExtractEachStandard)
        if alsoExtractEachStandard:
            productPath, qcTable = self._reduce_stacked_then_each_object_frame(allObjectFrames, allObjectNames)
        elif extractEach:
            productPath, qcTable = self._reduce_each_object_frame(allObjectFrames, allObjectNames)
        else:
            productPath, forceFailure = self._reduce_object_frames(allObjectFrames)
            qcTable = self.report_output()
            self.clean_up(forceFail=forceFailure)

        self.log.debug("completed the ``produce_product`` method")
        return productPath, qcTable

    def _should_extract_each_frame(self, allObjectFrames):
        """*decide whether to extract each object frame on its own instead of stacking them (DY-1611)*

        The `stack_before_extraction` recipe setting (default True) chooses, and the command-line option
        ``--no-stack`` overrides it. A flux standard is always stacked here, because the reducer predicts one
        response curve per SOF and later SOFs read it. With ``--no-stack`` its frames are extracted on their own as
        well (see `_also_extract_each_standard_frame`). Frames that are not stare frames (the PAE flat-lamp and
        nodding-standard fallbacks) are always stacked too.

        **Key Arguments:**

        - ``allObjectFrames`` -- the object frames of the SOF

        **Return:**

        - ``extractEach`` -- True to extract each frame on its own
        """
        isOverridden = self._stacking_override_requested()
        isStacking = self.recipeSettings.get("stack_before_extraction", True) and not isOverridden
        if isStacking or len(allObjectFrames) < 2:
            return False
        if not getattr(self, "readStareFrames", False):
            return False
        if self.generateReponseCurve:
            if not isOverridden:
                self.log.warning(
                    "`stack_before_extraction` is False, but flux-standard frames are always stacked so that the SOF "
                    "gives one response curve."
                )
            return False
        return True

    def _also_extract_each_standard_frame(self, allObjectFrames):
        """*decide whether a stacked flux standard is extracted frame by frame as well (``--no-stack``)*

        The stacked pass always runs, because it writes the SOF-level response curve that later SOFs read. The
        command-line option ``--no-stack`` adds a pass over each frame on its own, with its own response curve.

        **Key Arguments:**

        - ``allObjectFrames`` -- the object frames of the SOF

        **Return:**

        - ``alsoExtractEach`` -- True to follow the stacked pass with a pass over each frame
        """
        return bool(
            self._stacking_override_requested()
            and self.generateReponseCurve
            and getattr(self, "readStareFrames", False)
            and len(allObjectFrames) >= 2
        )

    def _delete_stale_qc_rows(self, extractEach, alsoExtractEachStandard):
        """*delete the QC rows an earlier run of this SOF left that this run will not rewrite*

        A frame's `<sof>_F<N>.sof` QC rows decide the status of `<sof>.sof` when the QC ranges are applied again.
        There are three paths:

        - stacked: the SOF's own rows are rewritten, so the earlier run's frame rows go
        - frame by frame: the frame rows are rewritten and the SOF has none, so its own rows go
        - stacked, then frame by frame (``--no-stack`` on a standard): both are rewritten, but an earlier run with
          more frames left rows beyond this run's frames, so every frame row goes and the SOF's own rows stay

        **Key Arguments:**

        - ``extractEach`` -- True when this run extracts each frame on its own
        - ``alsoExtractEachStandard`` -- True when this run is stacked and then extracts each frame on its own
        """
        if not getattr(self, "conn", None):
            return
        with closing(self.conn.cursor()) as cursor:
            # PLAIN TUPLES, WHATEVER ROW FACTORY THE CONNECTION CARRIES
            cursor.row_factory = None
            if extractEach and not alsoExtractEachStandard:
                staleNames = [f"{self.sofName}.sof"]
            else:
                framePrefix = f"{self.sofName}_F"
                rows = cursor.execute(
                    "select distinct sof_name from quality_control where substr(sof_name, 1, ?) = ?",
                    (len(framePrefix), framePrefix),
                ).fetchall()
                frameName = re.compile(re.escape(framePrefix) + r"\d+\.sof")
                staleNames = [name for (name,) in rows if frameName.fullmatch(name)]
            cursor.executemany("delete from quality_control where sof_name = ?", [(name,) for name in staleNames])

    def _reduce_each_object_frame(self, allObjectFrames, allObjectNames):
        """*reduce and extract each object frame on its own, with its own products and QCs*

        Frame N of the SOF is named ``<sof>_F<N>``. Its QCs are reported under that name, and its product
        header lists only its own raw frame. The recipe passes or fails as a whole: one frame with a failed QC
        fails the SOF.

        **Key Arguments:**

        - ``allObjectFrames`` -- the object frames, in observation order
        - ``allObjectNames`` -- the raw file names of those frames

        **Return:**

        - ``productPath`` -- the path to the last frame's product
        - ``qcTable`` -- the QCs of every frame
        """
        productPath, self.qc, forceFailure = self._extract_frames_on_their_own(allObjectFrames, allObjectNames)
        self.clean_up(forceFail=forceFailure)
        return productPath, self.qc

    def _reduce_stacked_then_each_object_frame(self, allObjectFrames, allObjectNames):
        """*reduce the stacked flux standard, then each of its frames on its own (``--no-stack``)*

        The stacked pass writes the SOF-level response curve and gives the product path. The frame passes write the
        ``<sof>_F<N>`` products, each with its own response curve and QCs. The QCs of both passes are kept, and one
        failed pass fails the SOF. The recipe cleans up once, at the end.

        **Key Arguments:**

        - ``allObjectFrames`` -- the object frames, in observation order
        - ``allObjectNames`` -- the raw file names of those frames

        **Return:**

        - ``productPath`` -- the path to the stacked product
        - ``qcTable`` -- the QCs of the stacked pass and of every frame
        """
        import copy

        import pandas as pd

        # THE RECIPE SETTINGS ARE settings[recipeName] ITSELF AND THE SKY-MODEL FALLBACK EDITS THEM IN PLACE, SO THE
        # STACKED PASS RUNS ON A COPY AND THE ORIGINAL IS RESTORED, UNTOUCHED, FOR THE FRAMES AND THE CALLER
        recipeSettings, subtractSky = self.recipeSettings, self.subtractSky
        try:
            self.recipeSettings = copy.deepcopy(recipeSettings)
            productPath, stackedFailure = self._reduce_object_frames(allObjectFrames)
            stackedQc = self.report_output()
        finally:
            self.recipeSettings, self.subtractSky = recipeSettings, subtractSky
        self.qc, self.products = self._empty_qc_and_product_tables(pd)

        _, frameQc, frameFailure = self._extract_frames_on_their_own(allObjectFrames, allObjectNames)
        self.qc = self._concat_qc_tables([stackedQc, frameQc])
        self.clean_up(forceFail=self._combine_failures(stackedFailure, frameFailure))
        return productPath, self.qc

    def _extract_frames_on_their_own(self, allObjectFrames, allObjectNames):
        """*reduce and extract each object frame on its own, leaving the cleaning up to the caller*

        Frame N of the SOF is named ``<sof>_F<N>`` while it is reduced. The SOF name, filename template, recipe
        settings and sky flag are restored afterwards, whatever happens.

        **Key Arguments:**

        - ``allObjectFrames`` -- the object frames, in observation order
        - ``allObjectNames`` -- the raw file names of those frames

        **Return:**

        - ``productPath`` -- the path to the last frame's product
        - ``qcTable`` -- the QCs of every frame
        - ``forceFailure`` -- the combined failure of the frames: False, True or a message
        """
        import copy

        import pandas as pd

        sofName, filenameTemplate = self.sofName, self.filenameTemplate
        # THE SKY-MODEL FALLBACK EDITS THE SETTINGS IN PLACE, SO EACH FRAME RUNS ON A COPY OF THE ORIGINAL, WHICH IS
        # RESTORED UNTOUCHED
        recipeSettings, subtractSky = self.recipeSettings, self.subtractSky
        templateStem, templateExtension = os.path.splitext(filenameTemplate)
        qcTables = []
        productPath = None
        forceFailure = False
        try:
            for index, (frame, name) in enumerate(zip(allObjectFrames, allObjectNames, strict=True), start=1):
                self.log.print(f"\n# EXTRACTING OBJECT FRAME {index} OF {len(allObjectFrames)} ({name}) ON ITS OWN")
                self.sofName = f"{sofName}_F{index}"
                self.filenameTemplate = f"{templateStem}_F{index}{templateExtension}"
                self.recipeSettings, self.subtractSky = copy.deepcopy(recipeSettings), subtractSky
                productPath, frameFailure = self._reduce_object_frames([frame], rawFrames=[name])
                forceFailure = self._combine_failures(forceFailure, frameFailure)
                qcTables.append(self.report_output())
                self.qc, self.products = self._empty_qc_and_product_tables(pd)
        finally:
            self.sofName, self.filenameTemplate = sofName, filenameTemplate
            self.recipeSettings, self.subtractSky = recipeSettings, subtractSky

        return productPath, self._concat_qc_tables(qcTables), forceFailure

    @staticmethod
    def _combine_failures(current, new):
        """*combine the failure of one frame with the failures of the frames before it*

        A failure message wins over a bare ``True``, and the first message wins over later ones, so the reason a
        frame failed still reaches the database.

        **Key Arguments:**

        - ``current`` -- the combined failure so far: False, True or a message
        - ``new`` -- the failure of the next frame: False, True or a message

        **Return:**

        - ``failure`` -- the combined failure
        """
        if isinstance(current, str):
            return current
        if isinstance(new, str):
            return new
        return current or new

    def _reduce_object_frames(self, allObjectFrames, rawFrames=False):
        """*stack the object frames, then detrend, sky-subtract, extract and flux-calibrate the stack*

        **Key Arguments:**

        - ``allObjectFrames`` -- the object frames to stack; one frame is used as it is
        - ``rawFrames`` -- the raw file names the product header lists. Default *False* lists every raw frame of
          the SOF

        **Return:**

        - ``productPath`` -- the path to the product
        - ``forceFailure`` -- True when the response curve or a slit-edge spill fails the recipe, or the failure message
          when a standard star has no object trace
        """
        arm = self.arm
        kw = self.kw

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
        self.update_fits_keywords(frame=combined_object, rawFrames=rawFrames)

        excludeMask = self._slit_illumination_exclude_mask(combined_object, master_flat, twoDMap, allObjectFrames)

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

        # THE SKY FIT TAKES THE FLAT ONLY WHEN THE FRAME WAS FLAT-FIELDED (DY-1517)
        skyMasterFlat = master_flat if master_flat is not False else None
        skymodelCCDData, skySubtractedCCDData, productPath = self._subtract_sky(
            combined_object, twoDMap, dispMap, orderTablePath, excludeMask, skyMasterFlat
        )

        unflattenedSkySubtractedCCDData = self._unflatten_sky_subtracted_frame(
            skySubtractedCCDData, master_flat, combined_object_notflattened, excludeMask
        )

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=combined_object_notflattened,
            show=self.debug,
            ext="data",
            stdWindow=3,
            title=False,
            surfacePlot=True,
        )

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
            noTraceFallback=True,
            scienceHeaders=self.scienceHeaders,
        )
        self.qc, self.products, mergedSpectumDF, orderJoins, extractionPath = optimalExtractor.extract()
        slitEdgeSpill = bool(optimalExtractor.slitEdgeOrders)
        extractionNotes = list(optimalExtractor.extractionNotes)

        # CHECK IF FLUX CALIBRATION IS NEEDED
        filePath_fluxcal = None

        forceFailure = False
        if responseFunctionPath:
            filePath_fluxcal = self._flux_calibrate_spectrum(
                responseFunctionPath, mergedSpectumDF, combined_object, allObjectFrames
            )

        elif self.generateReponseCurve and optimalExtractor.traceFallback:
            # A BOXCAR AT THE SLIT CENTRE IS NOT A STANDARD-STAR SPECTRUM TO FIT A RESPONSE CURVE TO
            self.log.error(NO_TRACE_STANDARD_MESSAGE)
            forceFailure = NO_TRACE_STANDARD_MESSAGE
            extractionNotes.append(("ERROR", NO_TRACE_STANDARD_MESSAGE))

        elif self.generateReponseCurve:
            forceFailure = (
                self._generate_response_curve(
                    unflattenedSkySubtractedCCDData,
                    skymodelCCDData,
                    twoDMap,
                    dispMap,
                    extractionPath,
                    orderJoins,
                    optimalExtractor.flatPassHandoff,
                    objectFrames=allObjectFrames,
                )
                or slitEdgeSpill
            )
            if slitEdgeSpill:
                self.log.error(
                    "The standard star spills over the slit edge, so this recipe and its response curve are flagged "
                    "as failed and will not be used to flux-calibrate other frames."
                )

        self._plot_merged_spectrum_qcs(
            mergedSpectumDF, orderJoins, filePath_fluxcal, extractionNotes, optimalExtractor.traceFallback
        )

        return productPath, forceFailure

    def _read_stare_object_frames(self, kw):
        """*read the science frames of the first stare input type that has any*

        Flux-standard stare frames set ``self.generateReponseCurve``. PAE mode sets ``self.subtractSky`` to False.
        ``self.readStareFrames`` is True only when the frames are stare object or standard frames, not one of the
        PAE flat-lamp or nodding-standard fallbacks.

        **Key Arguments:**

        - ``kw`` -- the keyword lookup

        **Return:**

        - ``allObjectFrames`` -- the science frames, as CCDData objects
        - ``allObjectNames`` -- the raw file names of those frames, as the product header records them
        """
        self.readStareFrames = False

        # OBJECT FRAMES
        allObjectFrames, allObjectNames = self._read_science_frames(
            [
                {kw("DPR_TYPE"): "OBJECT", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
                {kw("DPR_TYPE"): "OBJECT,ASYNC", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"},
            ]
        )

        # FLUX STD FRAMES
        if not len(allObjectFrames):
            allObjectFrames, allObjectNames = self._read_science_frames(
                [{kw("DPR_TYPE"): "STD,FLUX", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"}]
            )
            if len(allObjectFrames):
                self.generateReponseCurve = True

        # TELLURIC STD FRAMES
        if not len(allObjectFrames):
            allObjectFrames, allObjectNames = self._read_science_frames(
                [{kw("DPR_TYPE"): "STD,TELLURIC", kw("DPR_TECH"): "ECHELLE,SLIT,STARE"}]
            )

        self.readStareFrames = bool(len(allObjectFrames))

        if not len(allObjectFrames) and "PAE" in self.settings and self.settings["PAE"]:
            allObjectFrames, allObjectNames = self._read_science_frames(
                [{kw("DPR_TYPE"): "LAMP,FLAT", kw("DPR_TECH"): "ECHELLE,PINHOLE"}]
            )
            self.log.warning("Processing a ORDER-TRACE frame with the stare-mode recipe")
            self.subtractSky = False

        if "PAE" in self.settings and self.settings["PAE"]:
            self.subtractSky = False

        if not len(allObjectFrames):
            allObjectFrames, allObjectNames = self._read_science_frames(
                [{kw("DPR_TYPE"): "STD,FLUX", kw("DPR_TECH"): "ECHELLE,SLIT,NODDING"}]
            )
            self.log.warning("Processing a NODDING frame with the stare-mode recipe")

        return allObjectFrames, allObjectNames

    def _read_science_frames(self, filterList):
        """*read every input frame that matches one of the filters, in filter order*

        **Key Arguments:**

        - ``filterList`` -- the header filters, as a list of dictionaries

        **Return:**

        - ``frames`` -- the matching frames, as CCDData objects
        - ``names`` -- their raw file names: the prepared file name without its `_pre` suffix
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        frames = []
        names = []
        for filters in filterList:
            for path in self.inputFrames.files_filtered(include_path=True, **filters):
                frames.append(
                    CCDData.read(
                        path,
                        hdu=0,
                        unit=u.electron,
                        hdu_uncertainty="ERRS",
                        hdu_mask="QUAL",
                        hdu_flags="FLAGS",
                        key_uncertainty_type="UTYPE",
                    )
                )
                # SAME NAME AS update_fits_keywords GIVES THE RAW FRAME, SO IT CAN FILTER ON IT
                names.append(os.path.basename(path).replace("_pre", ""))
        return frames, names

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

        # REFUSE A RESP PRODUCT OF AN OLD BASIS NOW, NOT AFTER THE FULL EXTRACTION
        if responseFunctionPath:
            from soxspipe.commonutils.response_basis import require_response_product_basis

            require_response_product_basis(responseFunctionPath)

        return orderTablePath, dispMap, twoDMap, responseFunctionPath

    def _slit_illumination_exclude_mask(self, combined_object, master_flat, twoDMap, allObjectFrames):
        """*measure the slit offset of each order and return the mask of the slit ends to leave out of the sky model*

        Nothing is measured, recorded or plotted unless the ``slit-illumination`` settings block has ``measure-offset``
        set to True. Updates ``self.qc`` and ``self.products``.

        **Key Arguments:**

        - ``combined_object`` -- the detrended object frame
        - ``master_flat`` -- the master flat CCDData, or False when no flat was used
        - ``twoDMap`` -- the path to the dispersion map image
        - ``allObjectFrames`` -- the stacked input frames, whose rotator angles are compared

        **Return:**

        - ``excludeMask`` -- boolean image, True for a pixel to leave out of the sky model, or *None* when the offset
          is not measured or could not be applied
        """
        from astropy.nddata import CCDData

        from soxspipe.commonutils.slit_illumination import measure_slit_illumination
        from soxspipe.commonutils.slit_illumination_mask import exclude_mask_from_result

        illuminationSettings = self._slit_illumination_settings()
        if illuminationSettings is None:
            return None
        binx, biny = self._read_frame_binning(combined_object)
        self._check_rotator_angles(combined_object, allObjectFrames)
        # PARTNERS ARE DELIBERATELY NOT PASSED. VIS ORDERS 2 AND 3 CAN NEVER BE MEASURED (NO SKYLINES BELOW 557.7 NM)
        # AND THE R AND I ORDERS MOVE WITH THE OPPOSITE SIGN, SO A BORROWED OFFSET WOULD BE WRONG
        result = measure_slit_illumination(
            log=self.log,
            frame=combined_object,
            masterFlat=master_flat if isinstance(master_flat, CCDData) else None,
            twoDMapPath=twoDMap,
            settings=self.settings,
            arm=self.arm,
            binx=binx,
            biny=biny,
            illuminationSettings=illuminationSettings,
            partners=None,
        )
        self._record_slit_illumination_qcs(result)
        self._write_slit_illumination_qc_plot(result)
        return exclude_mask_from_result(self.log, result, twoDMap, binx, biny, self.kw)

    def _slit_illumination_settings(self):
        """*the ``slit-illumination`` settings block, or None when the offset is not to be measured*

        A ``measure-offset`` value that is not a boolean is logged as a warning and treated as off.
        """
        block = self.recipeSettings.get("slit-illumination")
        if not isinstance(block, dict) or "measure-offset" not in block:
            return None
        measureOffset = block["measure-offset"]
        if not isinstance(measureOffset, bool):
            self.log.warning(
                f"slit-illumination measure-offset = {measureOffset!r} is not True or False, so the slit offset is "
                "not measured"
            )
            return None
        return block if measureOffset else None

    def _read_frame_binning(self, frame):
        """*the frame's binning from its header: 1 by 1 on the NIR arm, or with a warning if either key is unusable*"""
        if self.arm == "NIR":
            return 1, 1
        try:
            return int(frame.header[self.kw("WIN_BINX")]), int(frame.header[self.kw("WIN_BINY")])
        except (KeyError, TypeError, ValueError) as error:
            self.log.warning(
                f"The frame header has no usable integer binning keywords ({type(error).__name__}: {error}), "
                "so a binning of 1 by 1 is used to read the 2D map."
            )
            return 1, 1

    def _check_rotator_angles(self, combined_object, allObjectFrames):
        """*log the rotator angle of the stacked frame and warn if the stacked input frames differ widely in it*

        The slit offsets depend on the rotator angle, and one offset is measured for the whole stack. A frame with no
        usable ABSROT keyword is left out and never fails the recipe.
        """
        from soxspipe.commonutils.slit_illumination_qc import rotator_angle_spread

        key = self.kw("ABSROT")
        stackedAngle = read_rotator_angle(combined_object.header, key)
        if stackedAngle is None:
            self.log.info("The stacked frame has no usable ABSROT keyword, so its rotator angle is not logged.")
        else:
            self.log.info(f"The ABSROT of the stacked frame is {stackedAngle:.1f} degrees.")
        angles = [read_rotator_angle(frame.header, key) for frame in allObjectFrames]
        spread = rotator_angle_spread([angle for angle in angles if angle is not None])
        if spread > ABSROT_SPREAD_WARNING_DEGREES:
            self.log.warning(
                f"The ABSROT of the stacked input frames spans {spread:.1f} degrees, more than "
                f"{ABSROT_SPREAD_WARNING_DEGREES:.1f} degrees, so one slit offset may not fit every frame."
            )

    def _record_slit_illumination_qcs(self, result):
        """*add the ``SLIT ILLUM`` QCs: an offset for each order with a usable one, how many, and whether applied*

        The QCs are deliberately absent from ``qc-acceptable-ranges``, so flexure never fails a science recipe.
        """
        from soxspipe.commonutils.slit_illumination_mask import usable_offsets
        from soxspipe.commonutils.toolkit import append_qc, order_label

        utcnow = utcnow_string()
        offsets = usable_offsets(result.offsets)
        qcs = [
            (
                f"SLIT ILLUM OFFSET O{order}",
                round(offset, 3),
                f"Slit offset of order {order_label(order, self.arm)}",
                "arcsec",
            )
            for order, offset in offsets.items()
        ]
        qcs.append(("SLIT ILLUM APPLIED", int(result.isApplied), "1 if the slit offsets masked the sky model", None))
        qcs.append(("SLIT ILLUM N ORDERS", len(offsets), "Number of orders with a usable slit offset", None))
        for qcName, qcValue, qcComment, qcUnit in qcs:
            self.qc = append_qc(
                self.qc,
                recipeName=self.recipeName,
                qcName=qcName,
                qcValue=qcValue,
                qcComment=qcComment,
                obsDateUtc=self.dateObs,
                reductionDateUtc=utcnow,
                qcUnit=qcUnit,
                toHeader=True,
            )

    def _write_slit_illumination_qc_plot(self, result):
        """*write the slit-illumination QC plot and add it to the products, if there is a profile to plot*"""
        from soxspipe.commonutils.slit_illumination_qc import plot_slit_illumination_qc

        if result.profiles.empty:
            return
        filename = self.filenameTemplate.replace(".fits", "_SLIT_ILLUM_QC_PLOT.pdf")
        filePath = f"{self.qcDir}/{filename}"
        plot_slit_illumination_qc(result, f"Slit illumination offsets ({self.arm.upper()})", filePath)
        self.add_product(
            productLabel="SLIT_ILLUMINATION_QC_PLOT",
            fileName=filename,
            filePath=filePath,
            productDesc="Skyline slit profile of each order against the master flat, with the fitted slit offset",
            reductionDateUtc=utcnow_string(),
            fileType="PDF",
            label="QC",
        )

    def _subtract_sky(self, combined_object, twoDMap, dispMap, orderTablePath, excludeMask=None, masterFlat=None):
        """*subtract the sky from the detrended frame, when sky subtraction is on, and write the sky products*

        Updates ``self.qc`` and ``self.products``. A sky subtraction that returns no model sets ``self.subtractSky``
        and the recipe setting to False and widens ``horne-extraction-slit-length`` to 30.

        **Key Arguments:**

        - ``combined_object`` -- the detrended object frame
        - ``twoDMap`` -- the path to the dispersion map image
        - ``dispMap`` -- the path to the dispersion map table
        - ``orderTablePath`` -- the path to the order table
        - ``excludeMask`` -- boolean image, True for a pixel to leave out of the sky model. Default *None*
        - ``masterFlat`` -- the master flat the object frame was divided by, so the sky fit can weight by the
          flat-fielded noise and leave out the low-flat pixels. Default *None*, for a frame that was not flat-fielded

        **Return:**

        - ``skymodelCCDData`` -- the sky model, or False when no sky was subtracted
        - ``skySubtractedCCDData`` -- the sky-subtracted frame, or when no sky was subtracted ``combined_object``
          itself (a copy with ``excludeMask`` added to its mask, if there is an ``excludeMask``)
        - ``productPath`` -- the path to the last sky product written, or None when none was written
        """
        import numpy as np

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
                excludeMask=excludeMask,
                masterFlat=masterFlat,
                ron=None if self.readNoise is None else self.readNoise.skyFit,
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
                skySubtractedCCDData = self._flag_excluded_pixels(combined_object, excludeMask)
                self.recipeSettings["horne-extraction-slit-length"] = 30.0
            else:
                productPath = self._write_sky_products(
                    skySubtractedCCDData, skymodelCCDData, skySubtractedResidualsCCDData
                )

                # ADD QUALITY CHECKS. PIXELS THE SKY SUBTRACTION FLAGGED AS MISSING (E.G. BETWEEN THE ORDERS)
                # ARE NOT DETECTOR BAD PIXELS, SO THEY STAY OUT OF THE BAD-PIXEL COUNT (DY-1284)
                skyFlags = np.asarray(skySubtractedCCDData.mask, dtype=bool)
                if combined_object.mask is None:
                    skyOnlyFlags = skyFlags
                else:
                    skyOnlyFlags = skyFlags & ~np.asarray(combined_object.mask, dtype=bool)
                self.qc = generic_quality_checks(
                    log=self.log,
                    frame=skySubtractedCCDData,
                    settings=self.settings,
                    recipeName=self.recipeName,
                    qcTable=self.qc,
                    excludeMask=skyOnlyFlags,
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
            skySubtractedCCDData = self._flag_excluded_pixels(combined_object, excludeMask)

        return skymodelCCDData, skySubtractedCCDData, productPath

    @staticmethod
    def _flag_excluded_pixels(combined_object, excludeMask):
        """*the frame to extract when no sky was subtracted, with the overlap-mask pixels flagged*

        **Key Arguments:**

        - ``combined_object`` -- the detrended object frame, which is not modified
        - ``excludeMask`` -- boolean image, True for a pixel the sky model left out, or *None*

        **Return:**

        - ``frame`` -- ``combined_object`` itself if ``excludeMask`` is None, else a copy whose mask is the union of
          its mask (none counts as nothing masked) and ``excludeMask``
        """
        import numpy as np

        if excludeMask is None:
            return combined_object
        flagged = combined_object.copy()
        existing = np.zeros(excludeMask.shape, dtype=bool) if combined_object.mask is None else combined_object.mask
        flagged.mask = np.asarray(existing, dtype=bool) | np.asarray(excludeMask, dtype=bool)
        return flagged

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
            productDesc="The sky-subtracted object (ERRS excludes the sky-model uncertainty)",
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
            productDesc="The sky background model (ERRS is a copy of the object uncertainty)",
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
        self, skySubtractedCCDData, master_flat, combined_object_notflattened, excludeMask=None
    ):
        """*return the frame the optimal extraction uses as its unflattened input*

        **Key Arguments:**

        - ``skySubtractedCCDData`` -- the sky-subtracted, flat-corrected frame
        - ``master_flat`` -- the master flat, or False when none is used
        - ``combined_object_notflattened`` -- the stacked frame before detrending
        - ``excludeMask`` -- boolean image of the pixels left out of the sky model, or *None*. Default *None*

        **Return:**

        - ``unflattenedSkySubtractedCCDData`` -- the sky-subtracted frame multiplied back by the master flat, the
          sky-subtracted frame itself when there is no flat, or the stacked frame when no sky was subtracted (a
          copy with ``excludeMask`` added to its mask, if there is an ``excludeMask``)

        With a flat, the data is the sky-subtracted data times the flat, the mask is the OR of the sky-subtracted
        and flat masks, and the header is the sky-subtracted frame's. The uncertainty is the sky-subtracted
        (flat-corrected) standard deviation scaled back by the absolute flat. The flat's own uncertainty is not
        propagated a second time, because the flat correction in ``detrend`` already folded it into the
        flat-corrected uncertainty.
        """
        import numpy as np
        from astropy.nddata import StdDevUncertainty

        if self.subtractSky:
            if master_flat:
                unflattenedSkySubtractedCCDData = skySubtractedCCDData.multiply(master_flat)
                # SCALE THE ERROR BACK BY THE FLAT; MULTIPLY WOULD COUNT THE FLAT'S OWN ERROR TWICE
                unflattenedSkySubtractedCCDData.uncertainty = StdDevUncertainty(
                    np.abs(master_flat.data) * skySubtractedCCDData.uncertainty.array
                )
                unflattenedSkySubtractedCCDData.header = skySubtractedCCDData.header
            else:
                unflattenedSkySubtractedCCDData = skySubtractedCCDData
        else:
            unflattenedSkySubtractedCCDData = self._flag_excluded_pixels(combined_object_notflattened, excludeMask)

        return unflattenedSkySubtractedCCDData

    def _flux_calibrate_spectrum(self, responseFunctionPath, mergedSpectumDF, combined_object, objectFrames):
        """*flux calibrate the merged spectrum with the supplied response function*

        Updates ``self.products``.

        **Key Arguments:**

        - ``responseFunctionPath`` -- the path to the response table
        - ``mergedSpectumDF`` -- the merged extracted spectrum
        - ``combined_object`` -- the detrended object frame, whose header supplies the exposure time
        - ``objectFrames`` -- the raw object frames stacked into ``combined_object``; their mean airmass is used

        **Return:**

        - ``filePath_fluxcal`` -- the path to the flux-calibrated spectrum
        """
        import pandas as pd

        from soxspipe.commonutils import flux_calibration
        from soxspipe.commonutils.sequence_airmass import mean_airmass

        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)

        detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

        self.log.print("# FLUX CALIBRATING THE SPECTRUM\n")
        fluxCalibrator = flux_calibration(
            log=self.log,
            responseFunction=responseFunctionPath,
            extractedSpectrum=mergedSpectumDF,
            settings=self.settings,
            airmass=mean_airmass([frame.header for frame in objectFrames], self.kw, self.log),
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
        flatPass=None,
        *,
        objectFrames,
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
        - ``flatPass`` -- the flat-fielded extraction's ``FlatPassHandoff``, whose sky shifts the unflattened
          extraction applies instead of tuning its own (DY-1680). Default *None*
        - ``objectFrames`` -- the raw standard-star frames stacked into the extraction; the response uses their mean
          airmass, as flux calibration does (DY-1713)

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
            scienceHeaders=self.scienceHeaders,
            flatPass=flatPass,
        )
        self.qc, self.products, _, _, extractionPath_notflat = optimalExtractor.extract()

        # GETTING THE RESPONSE
        from soxspipe.commonutils import response_function
        from soxspipe.commonutils.sequence_airmass import mean_airmass

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
            airmass=mean_airmass([frame.header for frame in objectFrames], self.kw, self.log),
        )
        self.qc, self.products, forceFailure = response.get()

        return forceFailure

    def _plot_merged_spectrum_qcs(
        self, mergedSpectumDF, orderJoins, filePath_fluxcal, extractionNotes=(), boxcarExtraction=False
    ):
        """*plot the merged spectrum QC and, when the spectrum was flux calibrated, the calibrated one too*

        Updates ``self.products``.

        **Key Arguments:**

        - ``mergedSpectumDF`` -- the merged extracted spectrum
        - ``orderJoins`` -- the order joins of the extraction
        - ``filePath_fluxcal`` -- the path to the flux-calibrated spectrum, or None when it was not calibrated
        - ``extractionNotes`` -- the ``(level, message)`` warnings and errors of the extraction, to print under
          both plots. Default *()*
        - ``boxcarExtraction`` -- title both plots as a robust boxcar extraction. Default *False*
        """
        from soxspipe.commonutils.extraction_notes import format_extraction_notes
        from soxspipe.commonutils.toolkit import plot_merged_spectrum_qc

        noteLines = format_extraction_notes(extractionNotes)
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
            extractionNotes=noteLines,
            boxcarExtraction=boxcarExtraction,
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
                extractionNotes=noteLines,
                boxcarExtraction=boxcarExtraction,
            )

        return

    # USE THE TAB-TRIGGER BELOW FOR NEW METHOD
    # xt-class-method

    # OVERRIDE METHOD ATTRIBUTES
    # method-override-tmpx
