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
        self.settings = settings
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
        uniqueOffsets = self._check_ab_balance(allFrameAOffsets, allFrameBOffsets)

        forceFailure = False
        if len(allFrameAOffsets) > 1 and len(uniqueOffsets) > 1:
            stackedSpectrum, orderJoins = self._reduce_each_nodding_cycle(
                allFrameA,
                allFrameB,
                allFrameANames,
                allFrameBNames,
                orderTablePath=orderTablePath,
                master_flat=master_flat,
            )
        else:
            stackedSpectrum, orderJoins, forceFailure = self._reduce_stacked_ab_pair(
                allFrameA,
                allFrameB,
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

        Sets ``self.dispMap`` and ``self.twoDMap``.

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
        - ``orderJoins`` -- the order joins of the last cycle
        """
        import pandas as pd

        from soxspipe.commonutils.toolkit import quicklook_image

        allSpectrumA = []
        allSpectrumB = []
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
            mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins = self.process_single_ab_nodding_cycle(
                aFrame=frameA,
                bFrame=frameB,
                locationSetIndex=sequenceCount,
                orderTablePath=orderTablePath,
                masterFlat=masterFlat,
            )
            if sequenceCount == 1:
                allSpectrumA = mergedSpectrumDF_A
                allSpectrumB = mergedSpectrumDF_B
            else:
                allSpectrumA = pd.concat([allSpectrumA, mergedSpectrumDF_A])
                allSpectrumB = pd.concat([allSpectrumB, mergedSpectrumDF_B])

            sequenceCount += 1
        stackedSpectrum, extractionPath = self.stack_extractions([allSpectrumA, allSpectrumB], orderJoins=orderJoins)

        return stackedSpectrum, orderJoins

    def _reduce_stacked_ab_pair(self, allFrameA, allFrameB, orderTablePath, master_flat):
        """*stack the A and B frames, extract the pair, and build a response curve when requested*

        **Key Arguments:**

        - ``allFrameA`` -- the A frames
        - ``allFrameB`` -- the B frames
        - ``orderTablePath`` -- the path to the order table
        - ``master_flat`` -- the master flat, or False

        **Return:**

        - ``stackedSpectrum`` -- the stacked extraction
        - ``orderJoins`` -- the order joins of the extraction
        - ``forceFailure`` -- the response curve's failure flag, or False when no response was built
        """
        forceFailure = False

        # STACKING A AND B SEQUENCES - ONLY IF JITTER IS NOT PRESENT
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
        )
        stackedSpectrum, extractionPath = self.stack_extractions(
            [mergedSpectrumDF_A, mergedSpectrumDF_B], orderJoins=orderJoins
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
            )
            stackedSpectrum_notflat, extractionPath_notflat = self.stack_extractions(
                [mergedSpectrumDF_A, mergedSpectrumDF_B],
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
    ):
        """*process a single AB nodding cycle*

        **Key Arguments:**

        - ``aFrame`` -- the frame taken at the A location. CCDData object.
        - ``bFrame`` -- the frame taken at the B location. CCDDate object.
        - ``locationSetIndex`` -- the index of the AB cycle
        - ``orderTablePath`` -- path to the order table
        - ``notFlattened`` -- if True, the extraction is performed on non-flattened data. Default *False*
        - ``masterFlat`` -- path to the master flat frame. Default *False*


        **Return:**

        - ``mergedSpectrumDF_A`` -- the order merged spectrum of nodding location A (dataframe)
        - ``mergedSpectrumDF_B`` -- the order merged spectrum of nodding location B (dataframe), or False
          outside nodding mode
        - ``orderJoins`` -- the order joins of the last extraction

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

        # THIS `== False` AND THE `== True` AT THE SAVE BELOW ARE NOT COMPLEMENTARY, AND THAT IS PINNED
        # BEHAVIOUR. test_save_single_frame_extractions_e712_comparison_controls_the_products_table WALKS SIX
        # SETTING VALUES AND SHOWS THAT A VALUE CAN SEND `False` TO THE EXTRACTOR HERE AND STILL HAVE ITS
        # RESULT DROPPED BELOW, OR THE REVERSE. REWRITING EITHER AS `not ...` OR `if ...:` COLLAPSES ONE HALF
        # OF THAT. SIM108 IS SUPPRESSED WITH IT: THE TERNARY FORM PLUS THIS SUPPRESSION GOES PAST 120
        # CHARACTERS.
        if self.recipeSettings["save_single_frame_extractions"] == False:  # noqa: E712, SIM108
            theseProducts = False
        else:
            theseProducts = self.products

        A_minus_B, B_minus_A, aFrame, bFrame = self._detrend_cycle_frames(
            A_minus_B_notflattened,
            B_minus_A_notflattened,
            aFrame,
            bFrame,
            notFlattened=notFlattened,
            masterFlat=masterFlat,
            orderTablePath=orderTablePath,
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
        )

        # THE `== True` COMPARISON IS THE OTHER HALF OF THE ASYMMETRY PINNED ABOVE, SO IT IS NOT `if ...:`.
        if self.recipeSettings["save_single_frame_extractions"] == True:  # noqa: E712
            self.products = theseProducts

        self.log.debug("completed the ``process_single_ab_nodding_cycle`` method")
        return mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins

    def _difference_cycle_frames(self, aFrame, bFrame):
        """*subtract the two frames of the cycle from each other and reapply their headers*

        The B-A difference belongs to nodding mode only. An offset run never builds one, and never reads
        the False this returns in its place.

        **Key Arguments:**

        - ``aFrame`` -- the frame taken at the A location
        - ``bFrame`` -- the frame taken at the B location

        **Return:**

        - ``A_minus_B_notflattened`` -- the A-B difference, carrying the A header
        - ``B_minus_A_notflattened`` -- the B-A difference carrying the B header, or False outside nodding mode
        """
        # SUBTRACTING A FROM B
        A_minus_B_notflattened = aFrame.subtract(bFrame)
        B_minus_A_notflattened = False
        if "nod" in self.recipeName:
            B_minus_A_notflattened = bFrame.subtract(aFrame)

        # REAPPLYING HEADERS
        hdr_A = aFrame.header
        hdr_B = bFrame.header
        A_minus_B_notflattened.header = hdr_A
        if "nod" in self.recipeName:
            B_minus_A_notflattened.header = hdr_B

        return A_minus_B_notflattened, B_minus_A_notflattened

    def _write_cycle_difference_frames(
        self, A_minus_B_notflattened, B_minus_A_notflattened, locationSetIndex, notFlattened
    ):
        """*announce the cycle and write its difference frames to the product directory*

        **Key Arguments:**

        - ``A_minus_B_notflattened`` -- the A-B difference
        - ``B_minus_A_notflattened`` -- the B-A difference, or False outside nodding mode
        - ``locationSetIndex`` -- the index of the AB cycle
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
    ):
        """*optimally extract the cycle's difference frames*

        Sets ``self.qc``. An offset run extracts the A-B difference only.

        **Key Arguments:**

        - ``A_minus_B`` -- the A-B difference to extract
        - ``B_minus_A`` -- the B-A difference to extract, or False outside nodding mode
        - ``A_minus_B_notflattened`` -- the unflattened A-B difference
        - ``B_minus_A_notflattened`` -- the unflattened B-A difference, or False outside nodding mode
        - ``aFrame`` -- the frame taken at the A location
        - ``bFrame`` -- the frame taken at the B location
        - ``theseProducts`` -- the products table to add the single-frame extractions to, or False
        - ``locationSetIndex`` -- the index of the AB cycle
        - ``notFlattened`` -- if True, the extraction is performed on non-flattened data

        **Return:**

        - ``mergedSpectrumDF_A`` -- the order merged spectrum of nodding location A (dataframe)
        - ``mergedSpectrumDF_B`` -- the order merged spectrum of nodding location B (dataframe), or False
        - ``orderJoins`` -- the order joins of the last extraction
        - ``theseProducts`` -- the products table the extractions returned
        """
        from soxspipe.commonutils import horne_extraction

        # EXTRACT THE A MINUS B FRAME
        optimalExtractor = horne_extraction(
            log=self.log,
            skySubtractedFrame=A_minus_B,
            unflattenedFrame=A_minus_B_notflattened,
            subtractedFrame=bFrame,
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
        )

        self.qc, theseProducts, mergedSpectrumDF_A, orderJoins, extractionFITSPathA = optimalExtractor.extract()

        # EXTRACT THE B MINUS A FRAME
        if "nod" in self.recipeName:
            optimalExtractor = horne_extraction(
                log=self.log,
                skySubtractedFrame=B_minus_A,
                unflattenedFrame=B_minus_A_notflattened,
                subtractedFrame=aFrame,
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
            )
            (
                self.qc,
                theseProducts,
                mergedSpectrumDF_B,
                orderJoins,
                extractionFITSPathB,
            ) = optimalExtractor.extract()
        else:
            mergedSpectrumDF_B = False

        return mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins, theseProducts

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
