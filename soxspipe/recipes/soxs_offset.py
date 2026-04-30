#!/usr/bin/env python
# encoding: utf-8
"""
*Reduce SOXS/Xshooter data taken in offset mode*

Author
: David Young & Marco Landoni

Date Created
: February 27, 2024
"""

################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
from .soxs_nod import soxs_nod
from soxspipe.commonutils.toolkit import (
    add_snr_efficiency_qcs,
    generic_quality_checks,
    spectroscopic_image_quality_checks,
    get_calibrations_path,
)
from fundamentals import tools
from builtins import object
import sys
import os
from soxspipe.commonutils.filenamer import filenamer
from os.path import expanduser

os.environ["TERM"] = "vt100"

# TODO: When combining spectra at the end, we use a simple sum. If we use sigma-clipping followed by a mean combine, we can remove CRHs for data sets with more than 1 AB cycle.


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
    recipe = soxs_offset(
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
        # INHERIT INITIALISATION FROM  base_recipe
        super(soxs_offset, self).__init__(
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
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames,
            recipeName=self.recipeName,
            ext=self.settings["data-extension"],
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_OFFSET - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # SORT IMAGE COLLECTION
        self.inputFrames.sort(["MJD-OBS"])
        if self.verbose:
            self.log.print("# RAW INPUT FRAMES - SUMMARY")
            self.log.print(self.inputFrames.summary)
            self.log.print("\n")

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        return None

    def produce_product(self):
        """*The code to generate the product of the soxs_offset recipe*

        **Return:**

        - ``productPath`` -- the path to the final product

        **Usage**

        ```python
        from soxspipe.recipes import soxs_offset
        recipe = soxs_offset(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        offsetFrame = recipe.produce_product()
        ```
        """
        self.log.debug("starting the ``produce_product`` method")

        from astropy.nddata import CCDData
        from astropy import units as u
        import pandas as pd
        from datetime import datetime
        from soxspipe.commonutils.toolkit import quicklook_image, plot_merged_spectrum_qc

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None
        master_bias = False
        master_flat = False
        dark = False

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

        except:
            responseFunctionPath = False

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

        # DIVIDING IN A AND B SEQUENCES
        allFrameON, allFrameOFF, allFrameONOffsets, allFrameOFFOffsets, allFrameONNames, allFrameOFFNames = (
            [],
            [],
            [],
            [],
            [],
            [],
        )

        # SPLIT FRAMES INTO ON (negative net offset: offsetRA + offsetDec < 0) AND OFF (zero or positive net offset)
        for frame, filename in zip(allObjectFrames, allFilenames):

            offsetRA = frame.header[kw(f"OFFSET_RA")]
            offsetDec = frame.header[kw(f"OFFSET_DEC")]

            if (offsetRA + offsetDec) < 0:
                allFrameONOffsets.append(offsetDec)
                allFrameON.append(frame)
                allFrameONNames.append(filename)
            else:
                allFrameOFFOffsets.append(offsetDec)
                allFrameOFF.append(frame)
                allFrameOFFNames.append(filename)

        uniqueOffsets = list(set(allFrameONOffsets))

        if len(allFrameONOffsets) != len(allFrameOFFOffsets):
            error = f"Found {len(allFrameONOffsets)} ON frames and {len(allFrameOFFOffsets)} OFF frames. The number of ON and OFF frames must be the same for offset reductions."
            self.log.error(
                f"Found {len(allFrameONOffsets)} ON frames and {len(allFrameOFFOffsets)} OFF frames. The number of ON and OFF frames must be the same for offset reductions."
            )
            raise Exception(error)

        if len(uniqueOffsets) == 0:
            error = f"Did not find any ON frames (frames with a negative net offset, i.e. offsetRA + offsetDec < 0). Please check the `HIERARCH ESO SEQ FIXOFF` header keywords in the provided offset frames."
            self.log.error(
                f"Did not find any ON frames (frames with a negative net offset, i.e. offsetRA + offsetDec < 0). Please check the `HIERARCH ESO SEQ FIXOFF` header keywords in the provided offset frames."
            )
            raise Exception(error)

        elif len(uniqueOffsets) > 1:
            s = "S"
        else:
            s = ""
        self.log.print(
            f"# PROCESSING {len(allFrameONOffsets)} ON-OFF OFFSET CYCLES WITH {len(uniqueOffsets)} UNIQUE PAIR{s} OF OFFSET LOCATIONS"
        )

        forceFailure = False
        if len(allFrameONOffsets) > 1 and len(uniqueOffsets) > 1:

            allSpectrumA = []
            allSpectrumB = []
            sequenceCount = 1
            # SORT frameON and frameOFF looping at their MJDOBS keyword in the header in order to the closest A and B frames in time
            allFrameON.sort(key=lambda x: x.header["MJD-OBS"])
            allFrameOFF.sort(key=lambda x: x.header["MJD-OBS"])

            for frameON, frameOFF, frameONName, frameOFFName in zip(
                allFrameON, allFrameOFF, allFrameONNames, allFrameOFFNames
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
                    # Save frameON and frameOFF to disk in temporary file
                    home = expanduser("~")
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
                if sequenceCount == 1:
                    allSpectrumA = mergedSpectrumDF_A
                else:
                    allSpectrumA = pd.concat([allSpectrumA, mergedSpectrumDF_A])

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

        else:

            # STACKING A AND B SEQUENCES - ONLY IF JITTER IS NOT PRESENT
            aFrame = self.clip_and_stack(
                frames=allFrameON, recipe="soxs_offset", ignore_input_masks=False, post_stack_clipping=True
            )

            bFrame = self.clip_and_stack(
                frames=allFrameOFF, recipe="soxs_offset", ignore_input_masks=False, post_stack_clipping=True
            )

            # INJECT KEYWORDS INTO HEADER
            self.update_fits_keywords(frame=aFrame)
            self.update_fits_keywords(frame=bFrame)

            if self.recipeSettings["use_flat"] and master_flat:
                masterFlat = master_flat
            else:
                masterFlat = False

            mergedSpectrumDF_A, _, orderJoins = self.process_single_ab_nodding_cycle(
                aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath=orderTablePath, masterFlat=master_flat
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
                    masterFlat=master_flat,
                )
                stackedSpectrum_notflat, extractionPath_notflat = self.stack_extractions(
                    [
                        mergedSpectrumDF_A,
                    ],
                    notFlattened=True,
                    orderJoins=orderJoins,
                )
                # GETTING THE RESPONSE
                self.log.print(f"# CALCULATING RESPONSE FUNCTION\n")
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

        # CHECK IF FLUX CALIBRATION IS REQUESTED
        filePath_fluxcal = None
        if responseFunctionPath:

            calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)
            from soxspipe.commonutils.flux_calibration import flux_calibration

            self.log.print(f"# PERFORMING FLUX CALIBRATION\n")
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
            self.log.print(f"# FLUX CALIBRATION COMPLETED\n")

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
            from astropy.table import Table
            from astropy.io import fits
            from astropy import units as u

            fluxcal_spec = Table.read(filePath_fluxcal, format="fits")
            fluxcal_spec["WAVE"] = fluxcal_spec["WAVE"] * u.nm
            fluxcal_spec["FLUX_COUNTS"] = fluxcal_spec["FLUX_CALIBRATED"]  # BACK COMPATIBILITY WITH THE CODE
            # ADD THE SNR COLUMN AND COPY VALUES FROM stackedSpectrum
            fluxcal_spec["SNR"] = stackedSpectrum["SNR"]
            productPath = filePath_fluxcal

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

        qcTable = self.report_output()
        self.clean_up(forceFail=forceFailure)

        self.log.debug("completed the ``produce_product`` method")

        

        return productPath, qcTable
