#!/usr/bin/env python
"""
*Reduce SOXS/Xshooter data taken in offset mode*

Author
: David Young & Marco Landoni

Date Created
: February 27, 2024
"""

################# GLOBAL IMPORTS ####################
import os
import sys

from soxspipe.commonutils.toolkit import (
    get_calibrations_path,
)

from .soxs_nod import soxs_nod

os.environ["TERM"] = "vt100"

# TODO: WHEN COMBINING SPECTRA AT THE END, WE USE A SIMPLE SUM. IF WE USE SIGMA-CLIPPING FOLLOWED BY A MEAN
# COMBINE, WE CAN REMOVE CRHS FOR DATA SETS WITH MORE THAN 1 AB CYCLE.


class soxs_offset(soxs_nod):
    """
    *Reduce SOXS/Xshooter data taken in offset mode*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- generate debug plots. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.

    **Usage**

    ```python
    from soxspipe.recipes import soxs_offset
    productPath, qcTable = soxs_offset(
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
    ):
        # INHERIT INITIALISATION FROM  BASE_RECIPE
        super().__init__(
            log=log,
            settings=settings,
            inputFrames=inputFrames,
            overwrite=overwrite,
            recipeName="soxs-offset",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        self.log = log
        log.debug("instantiating a new 'soxs_offset' object")
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
        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_OFFSET - NO MORE, NO LESS.
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


    def produce_product(self):
        """*The code to generate the product of the soxs_offset recipe*

        **Return:**

        - ``productPath`` -- the path to the final product
        - ``qcTable`` -- the quality control table the recipe reports

        **Usage**

        ```python
        from soxspipe.recipes import soxs_offset
        recipe = soxs_offset(
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

        allObjectFrames, allFilenames = self._read_offset_object_frames(kw)
        master_flat, orderTablePath, responseFunctionPath = self._read_calibration_inputs(kw, arm)

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
            allFrameON,
            allFrameOFF,
            allFrameONOffsets,
            allFrameOFFOffsets,
            allFrameONNames,
            allFrameOFFNames,
        ) = self._split_on_off_frames(allObjectFrames, allFilenames, kw)
        uniqueOffsets = self._check_on_off_balance(allFrameONOffsets, allFrameOFFOffsets)

        if len(allFrameONOffsets) > 1 and len(uniqueOffsets) > 1:
            stackedSpectrum, productPath, orderJoins, forceFailure = self._reduce_each_offset_cycle(
                allFrameON,
                allFrameOFF,
                allFrameONNames,
                allFrameOFFNames,
                orderTablePath=orderTablePath,
                master_flat=master_flat,
            )
        else:
            stackedSpectrum, productPath, orderJoins, forceFailure = self._reduce_stacked_offset_pair(
                allFrameON,
                allFrameOFF,
                orderTablePath=orderTablePath,
                master_flat=master_flat,
            )

        # CHECK IF FLUX CALIBRATION IS REQUESTED
        filePath_fluxcal = None
        if responseFunctionPath:
            filePath_fluxcal = self._flux_calibrate_stack(responseFunctionPath, stackedSpectrum, allFrameON)

        self._plot_stacked_spectrum_qcs(stackedSpectrum, orderJoins, filePath_fluxcal)
        if filePath_fluxcal:
            productPath = filePath_fluxcal

        qcTable = self.report_output()
        self.clean_up(forceFail=forceFailure)

        self.log.debug("completed the ``produce_product`` method")

        return productPath, qcTable

    def _read_offset_object_frames(self, kw):
        """*read the offset science frames of the first frame type that has any*

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

            add_filters = {kw("DPR_TYPE"): t, kw("DPR_TECH"): "ECHELLE,SLIT,OFFSET"}
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                if t == "STD,FLUX" and "-std" not in self.recipeName:
                    self.recipeName += "-std"
                    self.recipeSettings = self.get_recipe_settings()
                    self.productDir = self.productDir.replace("soxs-offset", "soxs-offset-std")
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

    def _read_calibration_inputs(self, kw, arm):
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

    def _split_on_off_frames(self, allObjectFrames, allFilenames, kw):
        """*split the science frames into ON and OFF frames by the sign of their net offset*

        **Key Arguments:**

        - ``allObjectFrames`` -- the science frames
        - ``allFilenames`` -- the science frames' file names, in the same order
        - ``kw`` -- the keyword lookup

        **Return:**

        - ``allFrameON``, ``allFrameOFF``, ``allFrameONOffsets``, ``allFrameOFFOffsets``, ``allFrameONNames``,
          ``allFrameOFFNames`` -- the ON and OFF frames, their declination offsets and their file names
        """
        # DIVIDING IN A AND B SEQUENCES
        (
            allFrameON,
            allFrameOFF,
            allFrameONOffsets,
            allFrameOFFOffsets,
            allFrameONNames,
            allFrameOFFNames,
        ) = (
            [],
            [],
            [],
            [],
            [],
            [],
        )

        # SPLIT FRAMES INTO ON (NEGATIVE NET OFFSET: OFFSETRA + OFFSETDEC < 0) AND OFF (ZERO OR POSITIVE NET OFFSET)
        for frame, filename in zip(allObjectFrames, allFilenames):

            offsetRA = frame.header[kw("OFFSET_RA")]
            offsetDec = frame.header[kw("OFFSET_DEC")]

            if (offsetRA + offsetDec) < 0:
                allFrameONOffsets.append(offsetDec)
                allFrameON.append(frame)
                allFrameONNames.append(filename)
            else:
                allFrameOFFOffsets.append(offsetDec)
                allFrameOFF.append(frame)
                allFrameOFFNames.append(filename)

        return (
            allFrameON,
            allFrameOFF,
            allFrameONOffsets,
            allFrameOFFOffsets,
            allFrameONNames,
            allFrameOFFNames,
        )

    def _check_on_off_balance(self, allFrameONOffsets, allFrameOFFOffsets):
        """*refuse an unbalanced or ON-less set of offset frames and announce the cycles to process*

        **Key Arguments:**

        - ``allFrameONOffsets`` -- the declination offsets of the ON frames
        - ``allFrameOFFOffsets`` -- the declination offsets of the OFF frames

        **Return:**

        - ``uniqueOffsets`` -- the distinct ON-frame offsets
        """
        uniqueOffsets = list(set(allFrameONOffsets))

        if len(allFrameONOffsets) != len(allFrameOFFOffsets):
            error = f"Found {len(allFrameONOffsets)} ON frames and {len(allFrameOFFOffsets)} OFF frames. The number of ON and OFF frames must be the same for offset reductions."
            self.log.error(
                f"Found {len(allFrameONOffsets)} ON frames and {len(allFrameOFFOffsets)} OFF frames. The number of ON and OFF frames must be the same for offset reductions."
            )
            raise Exception(error)

        if len(uniqueOffsets) == 0:
            error = "Did not find any ON frames (frames with a negative net offset, i.e. offsetRA + offsetDec < 0). Please check the `HIERARCH ESO SEQ FIXOFF` header keywords in the provided offset frames."
            self.log.error(
                "Did not find any ON frames (frames with a negative net offset, i.e. offsetRA + offsetDec < 0). Please check the `HIERARCH ESO SEQ FIXOFF` header keywords in the provided offset frames."
            )
            raise Exception(error)

        if len(uniqueOffsets) > 1:
            s = "S"
        else:
            s = ""
        self.log.print(
            f"# PROCESSING {len(allFrameONOffsets)} ON-OFF OFFSET CYCLES WITH {len(uniqueOffsets)} UNIQUE PAIR{s} OF OFFSET LOCATIONS"
        )

        return uniqueOffsets

    def _reduce_each_offset_cycle(
        self, allFrameON, allFrameOFF, allFrameONNames, allFrameOFFNames, orderTablePath, master_flat
    ):
        """*extract each ON-OFF cycle on its own and stack the extractions*

        Sorts ``allFrameON`` and ``allFrameOFF`` in place by observation date.

        **Key Arguments:**

        - ``allFrameON`` -- the ON frames
        - ``allFrameOFF`` -- the OFF frames
        - ``allFrameONNames`` -- the ON frames' file names
        - ``allFrameOFFNames`` -- the OFF frames' file names
        - ``orderTablePath`` -- the path to the order table
        - ``master_flat`` -- the master flat, or False

        **Return:**

        - ``stackedSpectrum`` -- the stacked extraction
        - ``productPath`` -- the path to the stacked extraction
        - ``orderJoins`` -- the order joins of the last cycle
        - ``forceFailure`` -- True when a response curve was requested, which this path cannot build
        """
        import pandas as pd

        from soxspipe.commonutils.toolkit import quicklook_image

        forceFailure = False

        allSpectrumA = []
        sequenceCount = 1
        # SORT FRAMEON AND FRAMEOFF LOOPING AT THEIR MJDOBS KEYWORD IN THE HEADER IN ORDER TO THE CLOSEST A AND B
        # FRAMES IN TIME
        allFrameON.sort(key=lambda x: x.header["MJD-OBS"])
        allFrameOFF.sort(key=lambda x: x.header["MJD-OBS"])

        for frameON, frameOFF, frameONName, frameOFFName in zip(
            allFrameON, allFrameOFF, allFrameONNames, allFrameOFFNames, strict=False
        ):

            self.log.print(f"Processing ON-OFF Offset Sequence {sequenceCount}")
            if False:
                quicklook_image(
                    log=self.log,
                    CCDObject=frameON,
                    show=True,
                    ext=False,
                    stdWindow=1,
                    title=False,
                    surfacePlot=False,
                    saveToPath=False,
                )
                quicklook_image(
                    log=self.log,
                    CCDObject=frameOFF,
                    show=False,
                    ext=False,
                    stdWindow=1,
                    title=False,
                    surfacePlot=False,
                    saveToPath=False,
                )
                # SAVE FRAMEON AND FRAMEOFF TO DISK IN TEMPORARY FILE
                filenameON = self.sofName + f"_A_{sequenceCount}.fits"
                filenameOFF = self.sofName + f"_B_{sequenceCount}.fits"
                filePathON = f"{self.productDir}/{filenameON}"
                filePathOFF = f"{self.productDir}/{filenameOFF}"
                frameON.write(filePathON, overwrite=True, checksum=True)
                frameOFF.write(filePathOFF, overwrite=True, checksum=True)

            rawFrames = []
            if "ARCFILE" in frameON.header:
                rawFrames.append(frameON.header["ARCFILE"])
                rawFrames.append(frameOFF.header["ARCFILE"])
            elif "ORIGFILE" in frameON.header:
                rawFrames.append(frameON.header["ORIGFILE"])
                rawFrames.append(frameOFF.header["ORIGFILE"])
            else:
                rawFrames.append(frameONName)
                rawFrames.append(frameOFFName)

            # INJECT KEYWORDS INTO HEADER
            self.update_fits_keywords(frame=frameON, rawFrames=rawFrames)
            self.update_fits_keywords(frame=frameOFF, rawFrames=rawFrames)

            if self.recipeSettings["use_flat"] and master_flat:
                masterFlat = master_flat
            else:
                masterFlat = False

            # PROCESSING SINGLE SEQUENCE
            mergedSpectrumDF_A, _, orderJoins = self.process_single_ab_nodding_cycle(
                aFrame=frameON,
                bFrame=frameOFF,
                locationSetIndex=sequenceCount,
                orderTablePath=orderTablePath,
                masterFlat=masterFlat,
            )
            allSpectrumA = mergedSpectrumDF_A if sequenceCount == 1 else pd.concat([allSpectrumA, mergedSpectrumDF_A])

            sequenceCount += 1
        stackedSpectrum, extractionPath = self.stack_extractions(
            [
                allSpectrumA,
            ],
            orderJoins=orderJoins,
        )
        productPath = extractionPath

        if self.generateReponseCurve:
            forceFailure = True

        return stackedSpectrum, productPath, orderJoins, forceFailure

    def _reduce_stacked_offset_pair(self, allFrameON, allFrameOFF, orderTablePath, master_flat):
        """*stack the ON and OFF frames, extract the pair, and build a response curve when requested*

        **Key Arguments:**

        - ``allFrameON`` -- the ON frames
        - ``allFrameOFF`` -- the OFF frames
        - ``orderTablePath`` -- the path to the order table
        - ``master_flat`` -- the master flat, or False

        **Return:**

        - ``stackedSpectrum`` -- the stacked extraction
        - ``productPath`` -- the path to the stacked extraction
        - ``orderJoins`` -- the order joins of the extraction
        - ``forceFailure`` -- the response curve's failure flag, or False when no response was built
        """
        forceFailure = False

        # STACKING A AND B SEQUENCES - ONLY IF JITTER IS NOT PRESENT
        aFrame = self.clip_and_stack(
            frames=allFrameON,
            recipe="soxs_offset",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )

        bFrame = self.clip_and_stack(
            frames=allFrameOFF,
            recipe="soxs_offset",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )

        # INJECT KEYWORDS INTO HEADER
        self.update_fits_keywords(frame=aFrame)
        self.update_fits_keywords(frame=bFrame)

        masterFlat = master_flat if self.recipeSettings["use_flat"] and master_flat else False

        mergedSpectrumDF_A, _, orderJoins = self.process_single_ab_nodding_cycle(
            aFrame=aFrame,
            bFrame=bFrame,
            locationSetIndex=1,
            orderTablePath=orderTablePath,
            masterFlat=masterFlat,
        )
        stackedSpectrum, extractionPath = self.stack_extractions(
            [
                mergedSpectrumDF_A,
            ],
            orderJoins=orderJoins,
        )
        productPath = extractionPath

        if self.generateReponseCurve:
            from soxspipe.commonutils import response_function

            mergedSpectrumDF_A, _, orderJoins = self.process_single_ab_nodding_cycle(
                aFrame=aFrame,
                bFrame=bFrame,
                locationSetIndex=1,
                orderTablePath=orderTablePath,
                notFlattened=True,
                masterFlat=masterFlat,
            )
            stackedSpectrum_notflat, extractionPath_notflat = self.stack_extractions(
                [
                    mergedSpectrumDF_A,
                ],
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

        return stackedSpectrum, productPath, orderJoins, forceFailure

    def _flux_calibrate_stack(self, responseFunctionPath, stackedSpectrum, allFrameON):
        """*flux calibrate the stacked extraction and record its products*

        **Key Arguments:**

        - ``responseFunctionPath`` -- the path to the response table
        - ``stackedSpectrum`` -- the stacked extraction
        - ``allFrameON`` -- the ON frames; the first one's header supplies the airmass and exposure time

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
            airmass=allFrameON[0].header.get("HIERARCH ESO TEL AIRM END"),
            exptime=allFrameON[0].header.get("EXPTIME"),
            extinctionPath=calibrationRootPath + "/" + self.detectorParams["extinction"],
            arm=self.arm,
            header=allFrameON[0].header,
            recipeName=self.recipeName,
            startNightDate=self.startNightDate,
            sofName=self.sofName,
            debug=self.debug,
        )
        filePath_fluxcal, products = fluxCalibrator.calibrate()
        self.products = pd.concat([self.products, products], ignore_index=True)
        self.log.print("# FLUX CALIBRATION COMPLETED\n")

        return filePath_fluxcal

    def _plot_stacked_spectrum_qcs(self, stackedSpectrum, orderJoins, filePath_fluxcal):
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
