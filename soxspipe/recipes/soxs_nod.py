#!/usr/bin/env python
"""
*Reduce SOXS/Xshooter data taken in nodding mode*

Author
: David Young & Marco Landoni

Date Created
: February 27, 2024
"""

################# GLOBAL IMPORTS ####################
import os
import sys

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.toolkit import (
    generic_quality_checks,
    get_calibrations_path,
    spectroscopic_image_quality_checks,
    utcnow_string,
)

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"

# TODO: WHEN COMBINING SPECTRA AT THE END, WE USE A SIMPLE SUM. IF WE USE SIGMA-CLIPPING FOLLOWED BY A MEAN
# COMBINE, WE CAN REMOVE CRHS FOR DATA SETS WITH MORE THAN 1 AB CYCLE.


class soxs_nod(base_recipe):
    """
    *Reduce SOXS/Xshooter data taken in nodding mode*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- generate debug plots. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.
    - ``recipeName`` -- the name of the recipe. Default "soxs-nod". This is used to retrieve the recipe settings from the settings dictionary and to name the product file (soxs_offset inherits soxs_nod)

    **Usage**

    ```python
    from soxspipe.recipes import soxs_nod
    productPath, qcTable = soxs_nod(
        log=log,
        settings=settings,
        inputFrames=fileList
    ).produce_product()
    ```
    """

    # EXTRACTION WARNINGS AND ERRORS FOR THE MERGED QC PLOTS.
    # AN IMMUTABLE DEFAULT, REBOUND (NEVER MUTATED) WHEN NOTES ARE ADDED
    extractionNotes = ()

    # THE UVB/VIS MASTER BIAS THAT TAKES THE BIAS PEDESTAL OUT OF THE A AND B ERRORS, OR False (DY-1390)
    masterBias = False

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
        recipeName="soxs-nod",
    ):
        # INHERIT INITIALISATION FROM  base_recipe
        super().__init__(
            log=log,
            settings=settings,
            inputFrames=inputFrames,
            overwrite=overwrite,
            recipeName=recipeName,
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        self.log = log
        log.debug("instansiating a new 'soxs_nod' object")
        self.inputFrames = inputFrames
        self.verbose = verbose

        # INITIAL ACTIONS
        self._collect_input_frames()
        self._verify_and_announce_input_frames()
        self._sort_and_report_input_frames()

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        return

    def _collect_input_frames(self):
        """*convert the input files to a ccdproc image collection*

        Sets ``self.inputFrames`` and ``self.supplementaryInput``.
        """
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (INPUTFRAMES >
        # IMAGEFILECOLLECTION)
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames,
            recipeName=self.recipeName,
            ext=self.settings["data-extension"],
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        return

    def _verify_and_announce_input_frames(self):
        """*verify the collected frames and report the result to the user*

        Sets ``self.imageType``, through ``verify_input_frames``.
        """
        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_NOD - NO MORE, NO LESS.
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
            self.log.print("\n")

        return

    def verify_input_frames(self):
        """*verify the input frame match those required by the soxs_nod recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        error = False

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()
        arm = self.arm

        # AN OFFSET REDUCTION INHERITS THIS METHOD, SO NAME THE RECIPE ACTUALLY RUNNING
        isOffset = "offset" in self.recipeName
        recipe, frames = ("offset", "offset") if isOffset else ("nod", "nodding")

        if not error:
            for i in imageTypes:
                if i not in ["OBJECT", "LAMP,FLAT", "STD,FLUX", "STD,TELLURIC"]:
                    error = (
                        f"Found a {i} file. Input frames for soxspipe {recipe} need to be an object/std {frames} "
                        f"frames, a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table "
                        f"(DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}) and a master-flat "
                        f"(MASTER_FLAT_{arm})."
                    )

        if not error:
            if isOffset:
                error = self._offset_input_frame_tech_error(imageTech, imageTypes, arm)
            else:
                error = self._nod_input_frame_tech_error(imageTech, imageTypes, arm)

        if not error:
            for i in [
                f"DISP_TAB_{self.arm}",
                f"ORDER_TAB_{self.arm}",
                f"DISP_IMAGE_{self.arm}",
            ]:
                if i not in imageCat:
                    error = (
                        f"Input frames for soxspipe {recipe} need to be an object/std {frames} frames, a dispersion "
                        f"map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location "
                        f"table (ORDER_TAB_{arm}) and a master-flat (MASTER_FLAT_{arm}). The sof file is missing a "
                        f"{i} frame."
                    )

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

    def _offset_input_frame_tech_error(self, imageTech, imageTypes, arm):
        """*report the first offset-mode technique the input frames should not carry*

        **Key Arguments:**

        - ``imageTech`` -- the DPR TECH value of each input frame
        - ``imageTypes`` -- the DPR TYPE value of each input frame, in the same order
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``error`` -- the rejection message of the last offending frame, or False when every frame passes
        """
        error = False
        for i, ii in zip(imageTech, imageTypes, strict=False):
            if ii in ["STD,FLUX", "STD,TELLURIC"]:
                pass
            elif i not in [
                "IMAGE",
                "ECHELLE,SLIT",
                "ECHELLE,MULTI-PINHOLE",
                "ECHELLE,SLIT,OFFSET",
                "ECHELLE,SLIT,NODDING",
            ]:
                error = (
                    f"Found a {i} file. Input frames for soxspipe offset need to be an object/std offset "
                    f"frames, a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table "
                    f"(DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}) and a master-flat "
                    f"(MASTER_FLAT_{arm})."
                )

        return error

    def _nod_input_frame_tech_error(self, imageTech, imageTypes, arm):
        """*report the first nodding-mode technique the input frames should not carry*

        **Key Arguments:**

        - ``imageTech`` -- the DPR TECH value of each input frame
        - ``imageTypes`` -- the DPR TYPE value of each input frame, in the same order
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``error`` -- the rejection message of the last offending frame, or False when every frame passes
        """
        error = False
        for i, ii in zip(imageTech, imageTypes, strict=False):
            if ii in ["STD,FLUX", "STD,TELLURIC"]:
                pass
            elif i not in [
                "IMAGE",
                "ECHELLE,SLIT",
                "ECHELLE,MULTI-PINHOLE",
                "ECHELLE,SLIT,NODDING",
            ]:
                error = (
                    f"Found a {i} file. Input frames for soxspipe nod need to be an object/std nodding "
                    f"frames, a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table "
                    f"(DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}) and a master-flat "
                    f"(MASTER_FLAT_{arm})."
                )

        return error

    def produce_product(self):
        """*The code to generate the product of the soxs_nod recipe*

        **Return:**

        - ``productPath`` -- the path to the final product. Always None for this recipe.
        - ``qcTable`` -- the quality control table the recipe reports

        **Usage**

        ```python
        from soxspipe.recipes import soxs_nod
        recipe = soxs_nod(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        productPath, qcTable = recipe.produce_product()
        ```
        """
        self.log.debug("starting the ``produce_product`` method")

        from soxspipe.commonutils.toolkit import quicklook_image

        arm = self.arm
        kw = self.kw

        productPath = None

        allObjectFrames, allFilenames = self._read_nod_object_frames(kw)
        # THE ATTACHED-ARC TEST MEASURES THE ARC FROM THE FIRST START AND LAST END OF THE BLOCK (DY-1666)
        self.scienceHeaders = [frame.header.copy() for frame in allObjectFrames]
        master_flat, orderTablePath, responseFunctionPath = self._read_nod_calibration_inputs(kw, arm)

        quicklook_image(
            log=self.log,
            CCDObject=allObjectFrames[0],
            show=False,
            ext=False,
            stdWindow=3,
            title=False,
            surfacePlot=True,
            saveToPath=False,
        )

        (
            allFrameA,
            allFrameB,
            allFrameAOffsets,
            allFrameBOffsets,
            allFrameANames,
            allFrameBNames,
        ) = self._split_ab_frames(allObjectFrames, allFilenames, kw)
        # THE SPLIT HOLDS THE FRAMES NOW. DROPPING THIS LIST LETS THE ALIGNED FRAMES REPLACE THE ORIGINALS IN MEMORY
        del allObjectFrames
        self._check_ab_balance(allFrameAOffsets, allFrameBOffsets)

        forceFailure = False
        if self._should_align_and_stack():
            allFrameA, allFrameB, shiftAndAdd = self._align_ab_frames(
                allFrameA, allFrameB, allFrameANames, allFrameBNames
            )
            stackedSpectrum, orderJoins, forceFailure = self._reduce_stacked_ab_pair(
                allFrameA,
                allFrameB,
                orderTablePath=orderTablePath,
                master_flat=master_flat,
                shiftAndAdd=shiftAndAdd,
            )
        else:
            stackedSpectrum, orderJoins = self._reduce_each_nodding_cycle(
                allFrameA,
                allFrameB,
                allFrameANames,
                allFrameBNames,
                orderTablePath=orderTablePath,
                master_flat=master_flat,
            )

        # CHECK IF FLUX CALIBRATION IS REQUESTED
        filePath_fluxcal = None
        if responseFunctionPath:
            filePath_fluxcal = self._flux_calibrate_nod_stack(responseFunctionPath, stackedSpectrum, allFrameA)

        self._plot_nod_stacked_spectrum_qcs(stackedSpectrum, orderJoins, filePath_fluxcal)

        qcTable = self.report_output()
        self.clean_up(forceFail=forceFailure)

        self.log.debug("completed the ``produce_product`` method")

        return productPath, qcTable

    def _should_align_and_stack(self):
        """*decide whether to align and stack the frames before extraction, or to extract each cycle on its own*

        The ``stack_before_extraction`` recipe setting decides, and defaults to True. A standard star is always
        aligned and stacked, so that its SOF still yields one response curve.

        **Return:**

        - ``shouldStack`` -- True to align the frames along the slit and stack them before extraction
        """
        isEnabled = bool(self.recipeSettings.get("stack_before_extraction", True))
        isStandard = self.generateReponseCurve or "-std" in self.recipeName
        if isStandard and not isEnabled:
            self.log.warning(
                "stack_before_extraction is off, but a standard star is always aligned and stacked so that the SOF "
                "still yields one response curve."
            )
            return True
        return isEnabled

    def _aligns_each_order(self):
        """*decide whether each order of a frame gets its own slit shift*

        Only SOXS VIS does. SOXS NIR and every X-Shooter arm get one shift per frame.

        **Return:**

        - ``isPerOrder`` -- True for SOXS VIS
        """
        return self.inst.upper() == "SOXS" and self.arm.upper() == "VIS"

    def _align_ab_frames(self, allFrameA, allFrameB, allFrameANames, allFrameBNames):
        """*align the A frames to each other and the B frames to each other along the slit*

        Each frame is traced on its cycle difference with its partner, so both groups are measured before either is
        resampled. The two lists are changed in place: each original frame is replaced by its aligned frame as soon
        as that exists. Updates ``self.qc`` and ``self.products``.

        The position of the object in the reference frame of each group is measured here too, even for a group of one
        frame. Their difference is the slit shift that puts the B-A difference onto the A-B difference.

        **Key Arguments:**

        - ``allFrameA`` -- the A frames
        - ``allFrameB`` -- the B frames
        - ``allFrameANames`` -- the A frames' file names
        - ``allFrameBNames`` -- the B frames' file names

        **Return:**

        - ``allFrameA`` -- the aligned A frames
        - ``allFrameB`` -- the aligned B frames
        - ``shiftAndAdd`` -- the ``shift_and_add_plan`` that combines the A-B and B-A differences, or *None* when the
          object position was not measured at both nod positions
        """
        from soxspipe.commonutils.nod_shift_and_add import build_shift_and_add_plan
        from soxspipe.commonutils.slit_alignment import frame_group

        results = self._align_frame_groups(
            [frame_group(allFrameA, allFrameB, allFrameANames), frame_group(allFrameB, allFrameA, allFrameBNames)],
            measureSingleFrames=True,
        )
        self._record_slit_alignment(
            results,
            obsDateUtc=allFrameA[0].header[self.kw("DATE_OBS")],
            hasShifts=len(allFrameA) > 1,
        )
        return allFrameA, allFrameB, build_shift_and_add_plan(self.log, *results)

    def _align_frame_groups(self, groups, measureSingleFrames=False):
        """*align each group of frames along the slit, to the first of them in time*

        **Key Arguments:**

        - ``groups`` -- the ``frame_group`` objects to align. Their ``frames`` lists are changed in place
        - ``measureSingleFrames`` -- if True, a group of one frame has its trace position measured too. Default *False*

        **Return:**

        - ``results`` -- the ``slit_alignment_result`` of each group
        """
        from soxspipe.commonutils.slit_alignment import align_frame_groups_along_slit

        return align_frame_groups_along_slit(
            log=self.log,
            groups=groups,
            twoDMapPath=self.twoDMap,
            dispersionMap=self.dispMap,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            recipeName=self.recipeName,
            perOrder=self._aligns_each_order(),
            startNightDate=self.startNightDate,
            measureSingleFrames=measureSingleFrames,
        )

    def _record_slit_alignment(self, results, obsDateUtc, hasShifts):
        """*record the fallback count as a QC and the shifts as a QC plot product*

        ``self.dateObs`` is not set until the extractions are stacked, so the observation date is passed in.
        Updates ``self.qc`` and ``self.products``.

        **Key Arguments:**

        - ``results`` -- the ``slit_alignment_result`` of each group aligned
        - ``obsDateUtc`` -- the observation date to record against the QC and the product
        - ``hasShifts`` -- True when the groups held more than one frame, so there are shifts to plot
        """
        import pandas as pd

        from soxspipe.commonutils.slit_alignment_qc import plot_slit_alignment_qc

        reductionDateUtc = utcnow_string()
        self.add_qc(
            qcName="SLIT ALIGN FALLBACK NUM",
            qcValue=sum(result.fallbackCount for result in results),
            qcComment="Number of slit shifts taken from other orders or frames, as no object trace was found",
            reductionDateUtc=reductionDateUtc,
            qcUnit=None,
            toHeader=True,
            obsDateUtc=obsDateUtc,
        )
        if not hasShifts:
            return

        filename = self.sofName + "_SLIT_ALIGN_QC_PLOT.pdf"
        filePath = f"{self.qcDir}/{filename}"
        plot_slit_alignment_qc(
            shifts=pd.concat([result.shifts for result in results], ignore_index=True),
            title=f"Slit alignment ({self.arm.upper()})",
            filePath=filePath,
        )
        self.add_product(
            productLabel="SLIT_ALIGNMENT_QC_PLOT",
            fileName=filename,
            filePath=filePath,
            productDesc="Slit shift of each frame, or order of a frame, applied before stacking",
            reductionDateUtc=reductionDateUtc,
            fileType="PDF",
            label="QC",
            obsDateUtc=obsDateUtc,
        )

    def _read_nod_object_frames(self, kw):
        """*read the nodding science frames of the first frame type that has any*

        Sets ``self.masterHeaderFrame``. Flux-standard input also renames the recipe to its ``-std`` variant,
        re-reads ``self.recipeSettings`` and moves ``self.productDir`` to match.

        **Key Arguments:**

        - ``kw`` -- the keyword lookup

        **Return:**

        - ``allObjectFrames`` -- the science frames, as CCDData objects
        - ``allFilenames`` -- the frames' file names, with the ``_pre`` suffix removed
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        # OBJECT/STANDARD FRAMES
        types = ["OBJECT", "STD,FLUX", "STD,TELLURIC"]
        allObjectFrames, allFilenames = [], []
        self.masterHeaderFrame = False
        for t in types:
            add_filters = {kw("DPR_TYPE"): t, kw("DPR_TECH"): "ECHELLE,SLIT,NODDING"}
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                if t == "STD,FLUX" and "-std" not in self.recipeName:
                    self.recipeName += "-std"
                    self.recipeSettings = self.get_recipe_settings()
                    self.productDir = self.productDir.replace("soxs-nod", "soxs-nod-std")
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
                allFilenames.append(os.path.basename(i).replace("_pre.", "."))
                if not self.masterHeaderFrame:
                    self.masterHeaderFrame = singleFrame.copy()
            if len(allObjectFrames):
                break

        return allObjectFrames, allFilenames

    def _read_nod_calibration_inputs(self, kw, arm):
        """*read the master flat and locate the order table, dispersion maps and response table*

        Sets ``self.dispMap``, ``self.twoDMap`` and ``self.masterBias``.

        **Key Arguments:**

        - ``kw`` -- the keyword lookup
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``master_flat`` -- the master flat, or False when none was supplied
        - ``orderTablePath`` -- the path to the order table
        - ``responseFunctionPath`` -- the path to the response table, or False when none was supplied
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        master_flat = False

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

        self.masterBias = self._read_master_bias(kw, arm)

        # FIND THE ORDER TABLE
        filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}"}
        orderTablePath = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE 2D MAP TABLE
        filterDict = {kw("PRO_CATG"): f"DISP_TAB_{arm}"}
        self.dispMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE 2D MAP IMAGE
        filterDict = {kw("PRO_CATG"): f"DISP_IMAGE_{arm}"}
        self.twoDMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # CHECK IF FLUX CALIBRATION IS REQUESTED
        try:
            filterDict = {kw("PRO_CATG"): f"RESP_TAB_{arm}"}
            responseFunctionPath = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        except IndexError as e:
            self.log.debug(f"produce_product: no response function frame for this arm, continuing: {e}")
            responseFunctionPath = False

        return master_flat, orderTablePath, responseFunctionPath

    def _read_master_bias(self, kw, arm):
        """*read the UVB/VIS master bias the SOF names, so the bias pedestal can leave the A and B errors*

        NIR frames carry no bias pedestal, so a NIR reduction never looks for one.

        **Key Arguments:**

        - ``kw`` -- the keyword lookup
        - ``arm`` -- the arm the frames were taken with

        **Return:**

        - ``masterBias`` -- the master bias, or False on NIR or when the SOF names none
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        if arm.upper() not in ("UVB", "VIS"):
            return False

        masterBias = False
        add_filters = {kw("PRO_CATG"): f"MASTER_BIAS_{arm}"}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            masterBias = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        if masterBias is False:
            self.log.warning(
                f"The sof file names no MASTER_BIAS_{arm} frame, so the {arm} difference-frame errors still count "
                "the bias pedestal as photon noise."
            )

        return masterBias

    def _split_ab_frames(self, allObjectFrames, allFilenames, kw):
        """*split the science frames into A and B frames by the sign of their cumulative nodding offset*

        **Key Arguments:**

        - ``allObjectFrames`` -- the science frames
        - ``allFilenames`` -- the science frames' file names, in the same order
        - ``kw`` -- the keyword lookup

        **Return:**

        - ``allFrameA``, ``allFrameB``, ``allFrameAOffsets``, ``allFrameBOffsets``, ``allFrameANames``,
          ``allFrameBNames`` -- the A and B frames, their offsets and their file names
        """
        # DIVIDING IN A AND B SEQUENCES
        (
            allFrameA,
            allFrameB,
            allFrameAOffsets,
            allFrameBOffsets,
            allFrameANames,
            allFrameBNames,
        ) = (
            [],
            [],
            [],
            [],
            [],
            [],
        )

        # CUMOFF Y IS THE OFFSET IN THE Y DIRECTION OF THE NODDING SEQUENCE. POSITIVE A, NEGATIVE B
        for frame, filename in zip(allObjectFrames, allFilenames, strict=False):
            # offset = frame.header[kw(f"NOD_CUMULATIVE_OFFSET{self.axisA.upper()}")]
            offset = frame.header[kw("NOD_CUMULATIVE_OFFSETY")]
            if offset == 0:
                pass
            if offset > 0:
                allFrameAOffsets.append(offset)
                allFrameA.append(frame)
                allFrameANames.append(filename)
            else:
                allFrameBOffsets.append(offset)
                allFrameB.append(frame)
                allFrameBNames.append(filename)

        return (
            allFrameA,
            allFrameB,
            allFrameAOffsets,
            allFrameBOffsets,
            allFrameANames,
            allFrameBNames,
        )

    def _check_ab_balance(self, allFrameAOffsets, allFrameBOffsets):
        """*refuse an unbalanced or A-less set of nodding frames and announce the cycles to process*

        **Key Arguments:**

        - ``allFrameAOffsets`` -- the cumulative offsets of the A frames
        - ``allFrameBOffsets`` -- the cumulative offsets of the B frames

        **Return:**

        - ``uniqueOffsets`` -- the distinct A-frame offsets
        """
        uniqueOffsets = list(set(allFrameAOffsets))

        if len(allFrameAOffsets) != len(allFrameBOffsets):
            error = f"Found {len(allFrameAOffsets)} A frames and {len(allFrameBOffsets)} B frames. The number of A and B frames must be the same for nodding reductions."
            self.log.error(
                f"Found {len(allFrameAOffsets)} A frames and {len(allFrameBOffsets)} B frames. The number of A and B frames must be the same for nodding reductions."
            )
            raise Exception(error)

        if len(uniqueOffsets) == 0:
            error = "Did not find frames with a positive offset. Please check the `NOD_CUMULATIVE_OFFSETY` header keyword in the providing nodding frames."
            self.log.error(
                "Did not find frames with a positive offset. Please check the `NOD_CUMULATIVE_OFFSETY` header keyword in the providing nodding frames."
            )
            raise Exception(error)

        s = "S" if len(uniqueOffsets) > 1 else ""
        self.log.print(
            f"# PROCESSING {len(allFrameAOffsets)} AB NODDING CYCLES WITH {len(uniqueOffsets)} UNIQUE "
            f"PAIR{s} OF OFFSET LOCATIONS"
        )

        return uniqueOffsets

    def _reduce_each_nodding_cycle(
        self, allFrameA, allFrameB, allFrameANames, allFrameBNames, orderTablePath, master_flat
    ):
        """*extract each AB cycle on its own and stack the extractions*

        Sorts ``allFrameA`` and ``allFrameB`` in place by observation date.

        **Key Arguments:**

        - ``allFrameA`` -- the A frames
        - ``allFrameB`` -- the B frames
        - ``allFrameANames`` -- the A frames' file names
        - ``allFrameBNames`` -- the B frames' file names
        - ``orderTablePath`` -- the path to the order table
        - ``master_flat`` -- the master flat, or False

        **Return:**

        - ``stackedSpectrum`` -- the stacked extraction
        - ``orderJoins`` -- the order joins of the last cycle that found an object trace (*None* if none did)
        """
        from soxspipe.commonutils.toolkit import quicklook_image

        allSpectrumA = []
        allSpectrumB = []
        orderJoins = None
        sequenceCount = 1
        # SORT FRAMEA AND FRAMEB LOOPING AT THEIR MJDOBS KEYWORD IN THE HEADER IN ORDER TO THE CLOSEST A AND B
        # FRAMES IN TIME
        allFrameA.sort(key=lambda x: x.header["MJD-OBS"])
        allFrameB.sort(key=lambda x: x.header["MJD-OBS"])

        for frameA, frameB, frameAName, frameBName in zip(
            allFrameA, allFrameB, allFrameANames, allFrameBNames, strict=False
        ):
            self.log.print(f"Processing AB Nodding Sequence {sequenceCount}")
            if False:
                import matplotlib

                matplotlib.use("MacOSX")
                quicklook_image(
                    log=self.log,
                    CCDObject=frameA,
                    show=True,
                    ext=False,
                    stdWindow=1,
                    title=False,
                    surfacePlot=False,
                    saveToPath=False,
                )
                quicklook_image(
                    log=self.log,
                    CCDObject=frameB,
                    show=False,
                    ext=False,
                    stdWindow=1,
                    title=False,
                    surfacePlot=False,
                    saveToPath=False,
                )
                # SAVE FRAMEA AND FRAMEB TO DISK IN TEMPORARY FILE
                filenameA = self.sofName + f"_A_{sequenceCount}.fits"
                filenameB = self.sofName + f"_B_{sequenceCount}.fits"
                filePathA = f"{self.productDir}/{filenameA}"
                filePathB = f"{self.productDir}/{filenameB}"
                frameA.write(filePathA, overwrite=True, checksum=True)
                frameB.write(filePathB, overwrite=True, checksum=True)

            rawFrames = []
            if "ARCFILE" in frameA.header:
                rawFrames.append(frameA.header["ARCFILE"])
                rawFrames.append(frameB.header["ARCFILE"])
            elif "ORIGFILE" in frameA.header:
                rawFrames.append(frameA.header["ORIGFILE"])
                rawFrames.append(frameB.header["ORIGFILE"])
            else:
                rawFrames.append(frameAName)
                rawFrames.append(frameBName)

            # INJECT KEYWORDS INTO HEADER
            self.update_fits_keywords(frame=frameA, rawFrames=rawFrames)
            self.update_fits_keywords(frame=frameB, rawFrames=rawFrames)

            masterFlat = master_flat if self.recipeSettings["use_flat"] and master_flat else False

            # PROCESSING SINGLE SEQUENCE
            mergedSpectrumDF_A, mergedSpectrumDF_B, cycleJoins = self.process_single_ab_nodding_cycle(
                aFrame=frameA,
                bFrame=frameB,
                locationSetIndex=sequenceCount,
                orderTablePath=orderTablePath,
                masterFlat=masterFlat,
            )
            # A CYCLE WITH NO OBJECT TRACE REPORTS NO ORDER JOINS, SO KEEP THE LATEST ONES AN EARLIER CYCLE FOUND
            if cycleJoins is not None:
                orderJoins = cycleJoins
            if sequenceCount == 1:
                allSpectrumA = mergedSpectrumDF_A
                allSpectrumB = mergedSpectrumDF_B
            else:
                allSpectrumA = self._append_spectrum(allSpectrumA, mergedSpectrumDF_A)
                allSpectrumB = self._append_spectrum(allSpectrumB, mergedSpectrumDF_B)

            sequenceCount += 1
        stackedSpectrum, extractionPath = self.stack_extractions([allSpectrumA, allSpectrumB], orderJoins=orderJoins)

        return stackedSpectrum, orderJoins

    @staticmethod
    def _append_spectrum(spectra, newSpectrum):
        """*append a cycle's spectrum to the spectra gathered so far, skipping cycles that found no object trace*

        **Key Arguments:**

        - ``spectra`` -- the spectra gathered so far (dataframe), or None when no cycle found a trace yet
        - ``newSpectrum`` -- the cycle's spectrum (dataframe), or None when it found no trace

        **Return:**

        - ``spectra`` -- the combined spectra, or None when neither input holds a spectrum
        """
        import pandas as pd

        found = [spectrum for spectrum in (spectra, newSpectrum) if spectrum is not None]
        return pd.concat(found) if found else None

    def _reduce_stacked_ab_pair(self, allFrameA, allFrameB, orderTablePath, master_flat, shiftAndAdd=None):
        """*stack the A and B frames, extract the pair, and build a response curve when requested*

        With a ``shiftAndAdd`` plan, the B-A difference is shifted onto the A-B difference and the two are averaged,
        so each pass makes one extraction. Without one, each pass extracts the A-B and the B-A difference.

        **Key Arguments:**

        - ``allFrameA`` -- the A frames
        - ``allFrameB`` -- the B frames
        - ``orderTablePath`` -- the path to the order table
        - ``master_flat`` -- the master flat, or False
        - ``shiftAndAdd`` -- the ``shift_and_add_plan`` that combines the two differences, or *None*. Default *None*

        **Return:**

        - ``stackedSpectrum`` -- the stacked extraction
        - ``orderJoins`` -- the order joins of the B-A extraction, or of the A-B extraction when the B-A
          extraction found no object trace. With a ``shiftAndAdd`` plan, those of the one extraction
        - ``forceFailure`` -- the response curve's failure flag, or False when no response was built
        """
        forceFailure = False

        # STACKING A AND B SEQUENCES, WHICH ARE ALREADY ALIGNED ALONG THE SLIT
        aFrame = self.clip_and_stack(
            frames=allFrameA,
            recipe="soxs_nod",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )

        bFrame = self.clip_and_stack(
            frames=allFrameB,
            recipe="soxs_nod",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )

        # INJECT KEYWORDS INTO HEADER
        self.update_fits_keywords(frame=aFrame)
        self.update_fits_keywords(frame=bFrame)

        masterFlat = master_flat if self.recipeSettings["use_flat"] and master_flat else False

        mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins = self.process_single_ab_nodding_cycle(
            aFrame=aFrame,
            bFrame=bFrame,
            locationSetIndex=1,
            orderTablePath=orderTablePath,
            masterFlat=masterFlat,
            shiftAndAdd=shiftAndAdd,
        )
        stackedSpectrum, extractionPath = self.stack_extractions(
            self._spectra_to_stack(mergedSpectrumDF_A, mergedSpectrumDF_B), orderJoins=orderJoins
        )

        if self.generateReponseCurve:
            from soxspipe.commonutils import response_function

            mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins = self.process_single_ab_nodding_cycle(
                aFrame=aFrame,
                bFrame=bFrame,
                locationSetIndex=1,
                orderTablePath=orderTablePath,
                notFlattened=True,
                masterFlat=masterFlat,
                shiftAndAdd=shiftAndAdd,
            )
            stackedSpectrum_notflat, extractionPath_notflat = self.stack_extractions(
                self._spectra_to_stack(mergedSpectrumDF_A, mergedSpectrumDF_B),
                notFlattened=True,
                orderJoins=orderJoins,
            )
            # GETTING THE RESPONSE
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

        return stackedSpectrum, orderJoins, forceFailure

    @staticmethod
    def _spectra_to_stack(mergedSpectrumDF_A, mergedSpectrumDF_B):
        """*the extractions of one pass that go into the stack: both nod positions, or the one combined extraction*

        **Key Arguments:**

        - ``mergedSpectrumDF_A`` -- the order merged spectrum of nodding location A, or of the combined frame
        - ``mergedSpectrumDF_B`` -- the order merged spectrum of nodding location B, or False when the pass made one
          extraction

        **Return:**

        - ``spectra`` -- the list of spectra to stack
        """
        if mergedSpectrumDF_B is False:
            return [mergedSpectrumDF_A]
        return [mergedSpectrumDF_A, mergedSpectrumDF_B]

    def _flux_calibrate_nod_stack(self, responseFunctionPath, stackedSpectrum, allFrameA):
        """*flux calibrate the stacked extraction and record its products*

        **Key Arguments:**

        - ``responseFunctionPath`` -- the path to the response table
        - ``stackedSpectrum`` -- the stacked extraction
        - ``allFrameA`` -- the A frames; the first one's header supplies the airmass and exposure time

        **Return:**

        - ``filePath_fluxcal`` -- the path to the flux-calibrated spectrum
        """
        import pandas as pd

        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)
        from soxspipe.commonutils.flux_calibration import flux_calibration

        self.log.print("# PERFORMING FLUX CALIBRATION\n")
        # TODO CHECK IF TAKING THE HEADER OF ONE FRAME IS OK
        fluxCalibrator = flux_calibration(
            log=self.log,
            responseFunction=responseFunctionPath,
            extractedSpectrum=stackedSpectrum,
            settings=self.settings,
            airmass=allFrameA[0].header.get("HIERARCH ESO TEL AIRM END"),
            exptime=allFrameA[0].header.get("EXPTIME"),
            extinctionPath=calibrationRootPath + "/" + self.detectorParams["extinction"],
            arm=self.arm,
            header=allFrameA[0].header,
            recipeName=self.recipeName,
            startNightDate=self.startNightDate,
            sofName=self.sofName,
            debug=self.debug,
        )
        filePath_fluxcal, products = fluxCalibrator.calibrate()
        self.products = pd.concat([self.products, products], ignore_index=True)
        self.log.print("# FLUX CALIBRATION COMPLETED\n")

        return filePath_fluxcal

    def _plot_nod_stacked_spectrum_qcs(self, stackedSpectrum, orderJoins, filePath_fluxcal):
        """*plot the merged-spectrum QC, and the flux-calibrated one when there is one*

        **Key Arguments:**

        - ``stackedSpectrum`` -- the stacked extraction
        - ``orderJoins`` -- the order joins
        - ``filePath_fluxcal`` -- the path to the flux-calibrated spectrum, or None
        """
        from soxspipe.commonutils.toolkit import plot_merged_spectrum_qc

        self.products, filePath = plot_merged_spectrum_qc(
            merged_orders=stackedSpectrum,
            products=self.products,
            log=self.log,
            qcDir=self.qcDir,
            filenameTemplate=self.filenameTemplate,
            noddingSequence=False,
            dateObs=self.dateObs,
            arm=self.arm,
            recipeName=self.recipeName,
            orderJoins=orderJoins,
            debug=self.debug,
            fluxCalibrated=False,
            qcTable=self.qc,
            settings=self.settings,
            extractionNotes=list(self.extractionNotes),
        )

        if filePath_fluxcal:
            from astropy import units as u
            from astropy.table import Table

            fluxcal_spec = Table.read(filePath_fluxcal, format="fits")
            fluxcal_spec["WAVE"] = fluxcal_spec["WAVE"] * u.nm
            fluxcal_spec["FLUX_COUNTS"] = fluxcal_spec["FLUX_CALIBRATED"]  # BACK COMPATIBILITY WITH THE CODE
            # ADD THE SNR COLUMN AND COPY VALUES FROM STACKEDSPECTRUM
            fluxcal_spec["SNR"] = stackedSpectrum["SNR"]

            self.products, filePath = plot_merged_spectrum_qc(
                merged_orders=fluxcal_spec,
                products=self.products,
                log=self.log,
                qcDir=self.qcDir,
                filenameTemplate=self.filenameTemplate,
                noddingSequence=False,
                dateObs=self.dateObs,
                arm=self.arm,
                recipeName=self.recipeName,
                orderJoins=orderJoins,
                debug=self.debug,
                fluxCalibrated=True,
                qcTable=self.qc,
                settings=self.settings,
                extractionNotes=list(self.extractionNotes),
            )

        return

    def process_single_ab_nodding_cycle(
        self,
        aFrame,
        bFrame,
        locationSetIndex,
        orderTablePath,
        notFlattened=False,
        masterFlat=False,
        shiftAndAdd=None,
    ):
        """*process a single AB nodding cycle*

        **Key Arguments:**

        - ``aFrame`` -- the frame taken at the A location. CCDData object.
        - ``bFrame`` -- the frame taken at the B location. CCDDate object.
        - ``locationSetIndex`` -- the cycle index (int), or the ``STACKED_LOCATION_SET`` label (str)
          that a stacked offset pair passes
        - ``orderTablePath`` -- path to the order table
        - ``notFlattened`` -- if True, the extraction is performed on non-flattened data. Default *False*
        - ``masterFlat`` -- path to the master flat frame. Default *False*
        - ``shiftAndAdd`` -- the ``shift_and_add_plan`` of a stacked nod pair, or *None*. With a plan, the B-A
          difference is shifted onto the A-B difference, the two are averaged, and that one frame is extracted,
          tagged as the A sequence. Default *None*, for one extraction of each difference


        **Return:**

        - ``mergedSpectrumDF_A`` -- the order merged spectrum of nodding location A (dataframe), or of the combined
          frame when there is a ``shiftAndAdd`` plan
        - ``mergedSpectrumDF_B`` -- the order merged spectrum of nodding location B (dataframe), or False
          outside nodding mode or when there is a ``shiftAndAdd`` plan
        - ``orderJoins`` -- the order joins of the B-A extraction, or of the A-B extraction when the B-A
          extraction found no object trace (*None* if neither did). An offset run has the A-B extraction only, and a
          ``shiftAndAdd`` run the one extraction of the combined frame

        **Usage:**

        ```python
        mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins = soxs_nod.process_single_ab_nodding_cycle(
            aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath=orderTablePath, masterFlat=masterFlat)
        ```
        """
        self.log.debug("starting the ``process_single_ab_nodding_cycle`` method")

        A_minus_B_notflattened, B_minus_A_notflattened = self._difference_cycle_frames(aFrame, bFrame)

        self._write_cycle_difference_frames(
            A_minus_B_notflattened, B_minus_A_notflattened, locationSetIndex, notFlattened
        )

        self._run_cycle_quality_checks(A_minus_B_notflattened, orderTablePath)

        combinedNotFlattened = self._combine_unflattened_differences(
            A_minus_B_notflattened, B_minus_A_notflattened, shiftAndAdd, locationSetIndex, notFlattened
        )

        theseProducts = self._single_frame_products_table()

        A_minus_B, B_minus_A, aFrame, bFrame = self._detrend_cycle_frames(
            A_minus_B_notflattened,
            B_minus_A_notflattened,
            aFrame,
            bFrame,
            notFlattened=notFlattened,
            masterFlat=masterFlat,
            orderTablePath=orderTablePath,
        )

        combinedFrames = self._combined_cycle_frames(
            A_minus_B, B_minus_A, A_minus_B_notflattened, combinedNotFlattened, shiftAndAdd
        )

        mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins, theseProducts = self._extract_cycle_spectra(
            A_minus_B,
            B_minus_A,
            A_minus_B_notflattened,
            B_minus_A_notflattened,
            aFrame,
            bFrame,
            theseProducts=theseProducts,
            locationSetIndex=locationSetIndex,
            notFlattened=notFlattened,
            combinedFrames=combinedFrames,
        )

        # THE `== True` COMPARISON IS THE OTHER HALF OF THE ASYMMETRY PINNED IN _single_frame_products_table, SO IT IS
        # NOT `if ...:`.
        if self.recipeSettings["save_single_frame_extractions"] == True:  # noqa: E712
            self.products = theseProducts

        self.log.debug("completed the ``process_single_ab_nodding_cycle`` method")
        return mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins

    def _single_frame_products_table(self):
        """*the products table the extractions add their single-frame products to, or False when they add none*

        **Return:**

        - ``theseProducts`` -- ``self.products``, or False when ``save_single_frame_extractions`` is False
        """
        # THIS `== False` AND THE `== True` AT THE SAVE IN process_single_ab_nodding_cycle ARE NOT COMPLEMENTARY,
        # AND THAT IS PINNED BEHAVIOUR.
        # test_save_single_frame_extractions_e712_comparison_controls_the_products_table WALKS SIX SETTING VALUES
        # AND SHOWS THAT A VALUE CAN SEND `False` TO THE EXTRACTOR HERE AND STILL HAVE ITS RESULT DROPPED THERE,
        # OR THE REVERSE. REWRITING EITHER AS `not ...` OR `if ...:` COLLAPSES ONE HALF OF THAT. SIM108 IS
        # SUPPRESSED WITH IT: THE TERNARY FORM PLUS THIS SUPPRESSION GOES PAST 120 CHARACTERS.
        if self.recipeSettings["save_single_frame_extractions"] == False:  # noqa: E712, SIM108
            theseProducts = False
        else:
            theseProducts = self.products
        return theseProducts

    def _difference_cycle_frames(self, aFrame, bFrame):
        """*subtract the two frames of the cycle from each other and reapply their headers*

        The B-A difference belongs to nodding mode only. An offset run never builds one, and never reads
        the False this returns in its place.

        The bias pedestal cancels in the difference, but the error of each frame still counts it as photon
        noise. When ``self.masterBias`` holds a master bias, the A and B errors are first rescaled to the
        photon noise of their bias-subtracted counts, as ``detrend`` does for a stare frame. The data are
        unchanged.

        **Key Arguments:**

        - ``aFrame`` -- the frame taken at the A location
        - ``bFrame`` -- the frame taken at the B location

        **Return:**

        - ``A_minus_B_notflattened`` -- the A-B difference, carrying the A header
        - ``B_minus_A_notflattened`` -- the B-A difference carrying the B header, or False outside nodding mode
        """
        aOperand, bOperand = aFrame, bFrame
        if self.masterBias is not False:
            aOperand = self._remove_bias_pedestal_from_uncertainty(aFrame, self.masterBias)
            bOperand = self._remove_bias_pedestal_from_uncertainty(bFrame, self.masterBias)

        # SUBTRACTING A FROM B
        A_minus_B_notflattened = aOperand.subtract(bOperand)
        B_minus_A_notflattened = False
        if "nod" in self.recipeName:
            B_minus_A_notflattened = bOperand.subtract(aOperand)

        # REAPPLYING HEADERS
        hdr_A = aFrame.header
        hdr_B = bFrame.header
        A_minus_B_notflattened.header = hdr_A
        if "nod" in self.recipeName:
            B_minus_A_notflattened.header = hdr_B

        return A_minus_B_notflattened, B_minus_A_notflattened

    def _combine_cycle_differences(self, abDifference, baDifference, shiftAndAdd):
        """*average the A-B difference with the B-A difference shifted onto it along the slit*

        **Key Arguments:**

        - ``abDifference`` -- the A-B difference
        - ``baDifference`` -- the B-A difference
        - ``shiftAndAdd`` -- the ``shift_and_add_plan`` of the pair

        **Return:**

        - ``combined`` -- the combined frame, carrying the A header
        """
        from soxspipe.commonutils.nod_shift_and_add import combine_nod_differences

        return combine_nod_differences(abDifference, baDifference, shiftAndAdd.orderShifts, shiftAndAdd.mapImages)

    def _combine_unflattened_differences(self, abDifference, baDifference, shiftAndAdd, locationSetIndex, notFlattened):
        """*combine the unflattened differences of a cycle, and in debug mode write the result*

        **Key Arguments:**

        - ``abDifference`` -- the unflattened A-B difference
        - ``baDifference`` -- the unflattened B-A difference
        - ``shiftAndAdd`` -- the ``shift_and_add_plan`` of the pair, or *None*
        - ``locationSetIndex`` -- the cycle index
        - ``notFlattened`` -- if True, the cycle is the unflattened repeat, which writes nothing

        **Return:**

        - ``combined`` -- the combined unflattened frame, or *None* when there is no plan
        """
        if shiftAndAdd is None:
            return None
        combined = self._combine_cycle_differences(abDifference, baDifference, shiftAndAdd)
        if not notFlattened:
            self._write_combined_difference_frame(combined, locationSetIndex)
        return combined

    def _combined_cycle_frames(self, abDifference, baDifference, abUnflattened, combinedUnflattened, shiftAndAdd):
        """*the combined frame to extract and the combined frame to trace, for a cycle with a plan*

        The two are one frame when the differences were not flat-fielded.

        **Key Arguments:**

        - ``abDifference`` -- the A-B difference to extract
        - ``baDifference`` -- the B-A difference to extract
        - ``abUnflattened`` -- the unflattened A-B difference
        - ``combinedUnflattened`` -- the combined unflattened frame, or *None* when there is no plan
        - ``shiftAndAdd`` -- the ``shift_and_add_plan`` of the pair, or *None*

        **Return:**

        - ``combinedFrames`` -- the combined frame and its unflattened counterpart, or *None* when there is no plan
        """
        if shiftAndAdd is None:
            return None
        if abDifference is abUnflattened:
            return combinedUnflattened, combinedUnflattened
        return self._combine_cycle_differences(abDifference, baDifference, shiftAndAdd), combinedUnflattened

    def _write_combined_difference_frame(self, combined, locationSetIndex):
        """*in debug mode, write the combined difference frame and add it to the products*

        Updates ``self.products``.

        **Key Arguments:**

        - ``combined`` -- the combined frame
        - ``locationSetIndex`` -- the cycle index
        """
        if not self.debug:
            return
        filename = self.sofName + f"_AB_COMBINED_{locationSetIndex}.fits"
        filePath = f"{self.productDir}/{filename}"
        combined.write(filePath, overwrite=True, checksum=True)
        self.add_product(
            productLabel="AB_COMBINED",
            fileName=filename,
            filePath=filePath,
            productDesc="The A-B difference averaged with the B-A difference shifted onto it along the slit",
            reductionDateUtc=utcnow_string(),
            fileType="FITS",
            label="PROD",
            obsDateUtc=combined.header[self.kw("DATE_OBS")],
        )

    def _write_cycle_difference_frames(
        self, A_minus_B_notflattened, B_minus_A_notflattened, locationSetIndex, notFlattened
    ):
        """*announce the cycle and write its difference frames to the product directory*

        **Key Arguments:**

        - ``A_minus_B_notflattened`` -- the A-B difference
        - ``B_minus_A_notflattened`` -- the B-A difference, or False outside nodding mode
        - ``locationSetIndex`` -- the cycle index (int), or the ``STACKED_LOCATION_SET`` label (str)
          that a stacked offset pair passes
        - ``notFlattened`` -- if True, the cycle is the unflattened pass used to calculate efficiency
        """
        # WRITE IN A FITS FILE THE A-B AND B-A FRAMES
        extraText = " (not flattened this time - needed to calculate efficiency)" if notFlattened else ""
        if "nod" in self.recipeName:
            self.log.print(f"\n# PROCESSING AB NODDING CYCLE {locationSetIndex} {extraText}")
        else:
            self.log.print(f"\n# PROCESSING ON-OFF OFFSET CYCLE {locationSetIndex} {extraText}")
        if "nod" in self.recipeName:
            filename = self.sofName + f"_AB_{locationSetIndex}.fits"
            filePath = f"{self.productDir}/{filename}"
            A_minus_B_notflattened.write(filePath, overwrite=True, checksum=True)

            filename = self.sofName + f"_BA_{locationSetIndex}.fits"
            filePath = f"{self.productDir}/{filename}"
            B_minus_A_notflattened.write(filePath, overwrite=True, checksum=True)
        else:
            filename = self.sofName + f"_ONOFF_{locationSetIndex}.fits"
            filePath = f"{self.productDir}/{filename}"
            A_minus_B_notflattened.write(filePath, overwrite=True, checksum=True)

        if False:
            from soxspipe.commonutils.toolkit import quicklook_image

            quicklook_image(
                log=self.log,
                CCDObject=A_minus_B_notflattened,
                show=True,
                ext="data",
                stdWindow=3,
                title=False,
                surfacePlot=True,
                saveToPath=False,
            )
            if "nod" in self.recipeName:
                quicklook_image(
                    log=self.log,
                    CCDObject=B_minus_A_notflattened,
                    show=True,
                    ext="data",
                    stdWindow=3,
                    title=False,
                    surfacePlot=True,
                    saveToPath=False,
                )

        return

    def _run_cycle_quality_checks(self, A_minus_B_notflattened, orderTablePath):
        """*run the generic and spectroscopic quality checks over the cycle's A-B difference*

        Sets ``self.qc``.

        **Key Arguments:**

        - ``A_minus_B_notflattened`` -- the A-B difference
        - ``orderTablePath`` -- the path to the order table
        """
        # TODO: ADD THESE CHECKS .... LIKELY FOR EACH AB CYCLE INDEX
        if True:
            self.qc = generic_quality_checks(
                log=self.log,
                frame=A_minus_B_notflattened,
                settings=self.settings,
                recipeName=self.recipeName,
                qcTable=self.qc,
            )
            self.qc = spectroscopic_image_quality_checks(
                log=self.log,
                frame=A_minus_B_notflattened,
                settings=self.settings,
                recipeName=self.recipeName,
                qcTable=self.qc,
                orderTablePath=orderTablePath,
            )

        return

    def _detrend_cycle_frames(
        self,
        A_minus_B_notflattened,
        B_minus_A_notflattened,
        aFrame,
        bFrame,
        notFlattened,
        masterFlat,
        orderTablePath,
    ):
        """*flat-field the cycle's differences and its two input frames, when a master flat was supplied*

        **Key Arguments:**

        - ``A_minus_B_notflattened`` -- the A-B difference
        - ``B_minus_A_notflattened`` -- the B-A difference, or False outside nodding mode
        - ``aFrame`` -- the frame taken at the A location
        - ``bFrame`` -- the frame taken at the B location
        - ``notFlattened`` -- if True, no frame is flat-fielded
        - ``masterFlat`` -- the master flat, or False
        - ``orderTablePath`` -- the path to the order table

        **Return:**

        - ``A_minus_B`` -- the A-B difference, flat-fielded when a master flat was used
        - ``B_minus_A`` -- the B-A difference, or False outside nodding mode
        - ``aFrame`` -- the A frame, flat-fielded when a master flat was used
        - ``bFrame`` -- the B frame, flat-fielded when a master flat was used
        """
        B_minus_A = False

        # THE `== False` COMPARISON IS DELIBERATE.
        # test_not_flattened_e712_comparison_controls_whether_detrend_runs PINS THAT `notFlattened=None`
        # SKIPS DETRENDING, WHICH `not notFlattened` WOULD REVERSE.
        if not isinstance(masterFlat, bool) and notFlattened == False:  # noqa: E712
            A_minus_B = self.detrend(
                inputFrame=A_minus_B_notflattened,
                master_bias=False,
                dark=False,
                master_flat=masterFlat,
                order_table=orderTablePath,
            )
            if "nod" in self.recipeName:
                B_minus_A = self.detrend(
                    inputFrame=B_minus_A_notflattened,
                    master_bias=False,
                    dark=False,
                    master_flat=masterFlat,
                    order_table=orderTablePath,
                )
            bFrame = self.detrend(
                inputFrame=bFrame,
                master_bias=False,
                dark=False,
                master_flat=masterFlat,
                order_table=orderTablePath,
            )
            aFrame = self.detrend(
                inputFrame=aFrame,
                master_bias=False,
                dark=False,
                master_flat=masterFlat,
                order_table=orderTablePath,
            )

        else:
            A_minus_B = A_minus_B_notflattened
            if "nod" in self.recipeName:
                B_minus_A = B_minus_A_notflattened

        return A_minus_B, B_minus_A, aFrame, bFrame

    def _extract_cycle_spectra(
        self,
        A_minus_B,
        B_minus_A,
        A_minus_B_notflattened,
        B_minus_A_notflattened,
        aFrame,
        bFrame,
        theseProducts,
        locationSetIndex,
        notFlattened,
        combinedFrames=None,
    ):
        """*optimally extract the cycle's difference frames*

        Sets ``self.qc``. An offset run extracts the A-B difference only. With ``combinedFrames``, the one combined
        frame is extracted in place of both differences.

        **Key Arguments:**

        - ``A_minus_B`` -- the A-B difference to extract
        - ``B_minus_A`` -- the B-A difference to extract, or False outside nodding mode
        - ``A_minus_B_notflattened`` -- the unflattened A-B difference
        - ``B_minus_A_notflattened`` -- the unflattened B-A difference, or False outside nodding mode
        - ``aFrame`` -- the frame taken at the A location
        - ``bFrame`` -- the frame taken at the B location
        - ``theseProducts`` -- the products table to add the single-frame extractions to, or False
        - ``locationSetIndex`` -- the cycle index (int), or the ``STACKED_LOCATION_SET`` label (str)
          that a stacked offset pair passes
        - ``notFlattened`` -- if True, the extraction is performed on non-flattened data
        - ``combinedFrames`` -- the combined frame and its unflattened counterpart, or *None* to extract the
          differences on their own. Default *None*

        **Return:**

        - ``mergedSpectrumDF_A`` -- the order merged spectrum of nodding location A (dataframe), or of the combined
          frame
        - ``mergedSpectrumDF_B`` -- the order merged spectrum of nodding location B (dataframe), or False
        - ``orderJoins`` -- the order joins of the B-A extraction, or of the A-B extraction when the B-A
          extraction found no object trace (*None* if neither did). An offset run has the A-B extraction only
        - ``theseProducts`` -- the products table the extractions returned
        """
        # EXTRACT THE A MINUS B FRAME, OR THE COMBINED FRAME WHICH CARRIES THE A HEADER TOO
        skySubtractedFrame, unflattenedFrame = combinedFrames or (A_minus_B, A_minus_B_notflattened)
        theseProducts, mergedSpectrumDF_A, orderJoins = self._extract_difference(
            skySubtractedFrame, unflattenedFrame, bFrame, theseProducts, locationSetIndex, notFlattened
        )
        if combinedFrames is not None or "nod" not in self.recipeName:
            return mergedSpectrumDF_A, False, orderJoins, theseProducts

        # EXTRACT THE B MINUS A FRAME
        theseProducts, mergedSpectrumDF_B, orderJoinsB = self._extract_difference(
            B_minus_A, B_minus_A_notflattened, aFrame, theseProducts, locationSetIndex, notFlattened
        )
        # AN EXTRACTION WITH NO OBJECT TRACE REPORTS NO ORDER JOINS, SO IT MUST NOT ERASE THOSE OF THE OTHER
        # POSITION
        if orderJoinsB is not None:
            orderJoins = orderJoinsB

        return mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins, theseProducts

    def _extract_difference(
        self, skySubtractedFrame, unflattenedFrame, subtractedFrame, theseProducts, locationSetIndex, notFlattened
    ):
        """*optimally extract one difference frame*

        Sets ``self.qc`` and ``self.extractionNotes``.

        **Key Arguments:**

        - ``skySubtractedFrame`` -- the difference frame to extract
        - ``unflattenedFrame`` -- the same difference without flat-fielding, on which the object trace is found
        - ``subtractedFrame`` -- the frame that was subtracted, whose flux gives the sky spectrum
        - ``theseProducts`` -- the products table to add the extraction to, or False
        - ``locationSetIndex`` -- the cycle index (int), or the ``STACKED_LOCATION_SET`` label (str)
        - ``notFlattened`` -- if True, the extraction is performed on non-flattened data

        **Return:**

        - ``theseProducts`` -- the products table the extraction returned
        - ``mergedSpectrumDF`` -- the order merged spectrum (dataframe), or *None* when no object trace was found
        - ``orderJoins`` -- the order joins (*None* when no object trace was found)
        """
        from soxspipe.commonutils import horne_extraction

        optimalExtractor = horne_extraction(
            log=self.log,
            skySubtractedFrame=skySubtractedFrame,
            unflattenedFrame=unflattenedFrame,
            subtractedFrame=subtractedFrame,
            twoDMapPath=self.twoDMap,
            settings=self.settings,
            recipeName=self.recipeName,
            recipeSettings=self.recipeSettings,
            qcTable=self.qc,
            productsTable=theseProducts,
            dispersionMap=self.dispMap,
            sofName=self.sofName,
            locationSetIndex=locationSetIndex,
            startNightDate=self.startNightDate,
            debug=self.debug,
            notFlattened=notFlattened,
            turnOffMP=self.turnOffMP,
            scienceHeaders=self.scienceHeaders,
        )
        self.qc, theseProducts, mergedSpectrumDF, orderJoins, _ = optimalExtractor.extract()
        self._collect_extraction_notes(optimalExtractor, notFlattened)
        return theseProducts, mergedSpectrumDF, orderJoins

    def _collect_extraction_notes(self, extractor, notFlattened):
        """*add an extractor's warnings and errors to the notes printed on the merged QC plots*

        The notes are tagged with the extractor's nodding sequence (e.g. ``A1``) and de-duplicated in first-seen
        order. The ``_NOTFLAT`` re-extraction is skipped: its spectrum feeds the response curve and is never
        plotted on the merged QC plots, and its own per-order plot already carries its notes.

        Sets ``self.extractionNotes``.

        **Key Arguments:**

        - ``extractor`` -- the extractor that has just run
        - ``notFlattened`` -- if True, the extraction was performed on non-flattened data and nothing is collected
        """
        from soxspipe.commonutils.extraction_notes import format_extraction_notes

        if notFlattened:
            return

        tag = extractor.noddingSequence.removeprefix("_")
        newNotes = [
            note
            for note in format_extraction_notes(extractor.extractionNotes, tag=tag)
            if note not in self.extractionNotes
        ]
        self.extractionNotes = (*self.extractionNotes, *newNotes)

    def stack_extractions(self, dataFrameList, notFlattened=False, orderJoins=None):
        """*merge individual AB cycles into a master extraction*

        **Key Arguments:**

        - ``dataFrameList`` -- a list of order-merged spectrum dataframes
        - ``notFlattened`` -- if True, the extraction was performed on non-flattened data and the products
          are written with a ``_NOTFLAT`` suffix. Default *False*
        - ``orderJoins`` -- the order joins the extraction reported. Default *None*

        **Return:**

        - ``stackedSpectrum`` -- the combined spectrum, as an astropy table
        - ``filePath`` -- the path the combined spectrum was written to

        **Usage:**

        ```python
        stackedSpectrum, filePath = soxs_nod.stack_extractions(
            [mergedSpectrumDF_A, mergedSpectrumDF_B], orderJoins=orderJoins)
        ```
        """
        self.log.debug("starting the ``stack_extractions`` method")

        import numpy as np
        import pandas as pd
        from astropy import units as u
        from astropy.table import Table
        from specutils import Spectrum1D

        from soxspipe.commonutils.phase3 import write_fits_table_to_disk
        from soxspipe.commonutils.toolkit import (
            add_snr_efficiency_qcs,
        )

        postfix = "_NOTFLAT" if notFlattened else ""

        # AN EXTRACTION THAT FOUND NO OBJECT TRACE RETURNS NO SPECTRUM, AND WITH NONE TO STACK THERE IS NO PRODUCT
        if all(spectrum is None for spectrum in dataFrameList):
            where = (
                "in the offset frame of any cycle"
                if "offset" in self.recipeName
                else "in either nod position of any nodding cycle"
            )
            error = (
                f"No object trace was found {where}, so there is no extracted spectrum to stack. "
                "Please check that the object is in the slit and bright enough to be traced."
            )
            self.log.error(error)
            raise Exception(error)

        # MERGE THE PANDAS DATAFRAMES MERDGED_ORDERS_A AND MERGEDSPECTRUMDF_B INTO A SINGLE DATAFRAME, THEN
        # GROUP BY WAVE AND SUM THE FLUXES

        merged_dataframe = pd.concat(dataFrameList)
        # BEFORE GROUPING, WE NEED TO TRUNCATE THE WAVELENGTH TO THE 4 DIGITS
        merged_dataframe["WAVE"] = merged_dataframe["WAVE"].apply(lambda x: round(float(x.value), 4))
        groupedDataframe = merged_dataframe.groupby(by="WAVE", as_index=False).median()

        self.filenameTemplate = self.sofName + ".fits"

        # PREPARING THE HEADER
        kw = keyword_lookup(log=self.log, settings=self.settings).get

        # SELECTING HEADER A_MINUS_B (IS THIS THE SAME?)
        self.update_fits_keywords(frame=self.masterHeaderFrame)
        header = self.masterHeaderFrame.header

        header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
        header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{self.arm}".upper()

        flux_orig = groupedDataframe["FLUX_COUNTS"].values * u.electron
        # THE RESULT IS UNUSED, BUT BUILDING IT VALIDATES THE GROUPED FLUX AND WAVELENGTH ARRAYS AGAINST EACH
        # OTHER AND RAISES WHEN THEY DISAGREE. DELETING IT WOULD REMOVE THAT CHECK.
        spectrum_orig = Spectrum1D(  # noqa: F841
            flux=flux_orig,
            spectral_axis=groupedDataframe["WAVE"].values * u.nm,
            bin_specification="center",
        )

        groupedDataframe["SNR"] = groupedDataframe["FLUX_COUNTS"].values / np.sqrt(groupedDataframe["VARIANCE"].values)

        # groupedDataframe = calculate_rolling_snr(dataframe=groupedDataframe, flux_column="FLUX_COUNTS", window_size=300)

        # groupedDataframe['signal'] = groupedDataframe['FLUX_COUNTS'].rolling(
        #     window=15, center=True).median().fillna(method='bfill').fillna(method='ffill').values
        # groupedDataframe['normalised_flux'] = groupedDataframe['FLUX_COUNTS'] / \
        #     groupedDataframe['signal']

        # PREPARING THE HDU
        for col, decimals in [
            ("WAVE", 2),
            ("FLUX_COUNTS", 3),
            ("SNR", 2),
            ("FLUX_DENSITY_COUNTS", 3),
        ]:
            # `decimals` IS BOUND IN A DEFAULT ARGUMENT SO THE LAMBDA CANNOT SEE A LATER LOOP VALUE. `.apply`
            # CALLS IT WITHIN THE SAME ITERATION, SO THIS CHANGES NOTHING; IT ONLY MAKES THAT EXPLICIT.
            groupedDataframe[col] = groupedDataframe[col].apply(lambda x, decimals=decimals: round(float(x), decimals))
        stackedSpectrum = Table.from_pandas(groupedDataframe, index=False)

        self.utcnow = utcnow_string()
        self.dateObs = header[kw("DATE_OBS")]

        self.qc = add_snr_efficiency_qcs(
            log=self.log,
            spectrumDF=groupedDataframe,
            qcTable=self.qc,
            orderJoins=orderJoins,
            recipeName=self.recipeName,
            dateObs=self.dateObs,
            arm=self.arm,
        )

        # WRITE PRODUCT TO DISK
        filename = self.filenameTemplate.replace(".fits", "_EXTRACTED_MERGED" + postfix + ".fits")
        filePath = f"{self.productDir}/{filename}"

        write_fits_table_to_disk(
            log=self.log,
            settings=self.settings,
            header=header,
            tables=[stackedSpectrum],
            filePath=filePath,
            qc=self.qc,
        )

        # SAVE THE TABLE STACKEDSPECTRUM TO DISK IN ASCII FORMAT
        asciiFilename = self.filenameTemplate.replace(".fits", "_EXTRACTED_MERGED" + postfix + ".txt")
        asciiFilePath = f"{self.productDir}/{asciiFilename}"
        stackedSpectrum2 = stackedSpectrum.copy()
        stackedSpectrum2["WAVE"] = stackedSpectrum2["WAVE"] * 10
        stackedSpectrum2["WAVE"].format = "{:.2f}"  # CONVERTING TO ANGSTROMS
        stackedSpectrum2.write(asciiFilePath, format="ascii", overwrite=True)

        self.add_product(
            productLabel="EXTRACTED_MERGED_TABLE",
            fileName=filename,
            filePath=filePath,
            productDesc="Table of the extracted source in each order. All nodding cycles combined.",
            reductionDateUtc=self.utcnow,
            fileType="FITS",
            label="PROD",
        )

        self.add_product(
            productLabel="EXTRACTED_MERGED_ASCII",
            fileName=asciiFilename,
            filePath=asciiFilePath,
            productDesc="Ascii version of extracted source spectrum",
            reductionDateUtc=self.utcnow,
            fileType="TXT",
            label="PROD",
        )

        self.log.debug("completed the ``stack_extractions`` method")
        return stackedSpectrum, filePath
