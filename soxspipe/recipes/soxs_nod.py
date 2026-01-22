#!/usr/bin/env python
# encoding: utf-8
"""
*Reduce SOXS/Xshooter data taken in nodding mode*

Author
: David Young & Marco Landoni

Date Created
: February 27, 2024
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
from .base_recipe import base_recipe
from soxspipe.commonutils.toolkit import (
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

    **Usage**

    ```python
    from soxspipe.recipes import soxs_nod
    recipe = soxs_nod(
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
    ):
        # INHERIT INITIALISATION FROM  base_recipe
        super(soxs_nod, self).__init__(
            log=log,
            settings=settings,
            inputFrames=inputFrames,
            overwrite=overwrite,
            recipeName="soxs-nod",
            command=command,
            debug=debug,
            verbose=verbose,
        )
        self.log = log
        log.debug("instansiating a new 'soxs_nod' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeSettings = settings[self.recipeName]

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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_nod - NO MORE, NO LESS.
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

    def verify_input_frames(self):
        """*verify the input frame match those required by the soxs_nod recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        kw = self.kw

        error = False

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()
        arm = self.arm

        if not error:
            for i in imageTypes:
                if i not in ["OBJECT", "LAMP,FLAT", "STD,FLUX", "STD,TELLURIC"]:
                    error = f"Found a {i} file. Input frames for soxspipe nod need to be an object/std nodding frames, a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}) and a master-flat (MASTER_FLAT_{arm})."

        if not error:
            for i in imageTech:
                if i not in ["IMAGE", "ECHELLE,SLIT", "ECHELLE,MULTI-PINHOLE", "ECHELLE,SLIT,NODDING"]:
                    error = f"Found a {i} file. Input frames for soxspipe nod need to be an object/std nodding frames, a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}) and a master-flat (MASTER_FLAT_{arm})."

        if not error:
            for i in [f"DISP_TAB_{self.arm}", f"ORDER_TAB_{self.arm}", f"DISP_IMAGE_{self.arm}"]:
                if i not in imageCat:
                    error = f"Input frames for soxspipe nod need to be an object/std nodding frames, a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}) and a master-flat (MASTER_FLAT_{arm}). The sof file is missing a {i} frame."

        if error:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            self.log.print("")
            raise TypeError(error)

        self.imageType = imageTypes[0]
        self.log.debug("completed the ``verify_input_frames`` method")
        return None

    def produce_product(self):
        """*The code to generate the product of the soxs_nod recipe*

        **Return:**

        - ``productPath`` -- the path to the final product

        **Usage**

        ```python
        from soxspipe.recipes import soxs_nod
        recipe = soxs_nod(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        nodFrame = recipe.produce_product()
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
            add_filters = {kw("DPR_TYPE"): t, kw("DPR_TECH"): "ECHELLE,SLIT,NODDING"}
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
        allFrameA, allFrameB, allFrameAOffsets, allFrameBOffsets, allFrameANames, allFrameBNames = (
            [],
            [],
            [],
            [],
            [],
            [],
        )

        # CUMOFF Y IS THE OFFSET IN THE Y DIRECTION OF THE NODDING SEQUENCE. POSITIVE A, NEGATIVE B
        for frame, filename in zip(allObjectFrames, allFilenames):
            # offset = frame.header[kw(f"NOD_CUMULATIVE_OFFSET{self.axisA.upper()}")]
            offset = frame.header[kw(f"NOD_CUMULATIVE_OFFSETY")]
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

        uniqueOffsets = list(set(allFrameAOffsets))

        if len(allFrameAOffsets) != len(allFrameBOffsets):
            error = f"Found {len(allFrameAOffsets)} A frames and {len(allFrameBOffsets)} B frames. The number of A and B frames must be the same for nodding reductions."
            self.log.error(
                f"Found {len(allFrameAOffsets)} A frames and {len(allFrameBOffsets)} B frames. The number of A and B frames must be the same for nodding reductions."
            )
            raise Exception(error)

        if len(uniqueOffsets) == 0:
            error = f"Did not find frames with a positive offset. Please check the `NOD_CUMULATIVE_OFFSETY` header keyword in the providing nodding frames."
            self.log.error(
                f"Did not find frames with a positive offset. Please check the `NOD_CUMULATIVE_OFFSETY` header keyword in the providing nodding frames."
            )
            raise Exception(error)

        elif len(uniqueOffsets) > 1:
            s = "S"
        else:
            s = ""
        self.log.print(
            f"# PROCESSING {len(allFrameAOffsets)} AB NODDING CYCLES WITH {len(uniqueOffsets)} UNIQUE PAIR{s} OF OFFSET LOCATIONS"
        )

        if len(allFrameAOffsets) > 1 and len(uniqueOffsets) > 1:

            allSpectrumA = []
            allSpectrumB = []
            sequenceCount = 1
            # SORT frameA and frameB looping at their MJDOBS keyword in the header in order to the closest A and B frames in time
            allFrameA.sort(key=lambda x: x.header["MJD-OBS"])
            allFrameB.sort(key=lambda x: x.header["MJD-OBS"])

            for frameA, frameB, frameAName, frameBName in zip(allFrameA, allFrameB, allFrameANames, allFrameBNames):

                self.log.print(f"Processing AB Nodding Sequence {sequenceCount}")
                if False:
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
                    # Save frameA and frameB to disk in temporary file
                    home = expanduser("~")
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

                if self.recipeSettings["use_flat"] and master_flat:
                    masterFlat = master_flat
                else:
                    masterFlat = False

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
            stackedSpectrum, extractionPath = self.stack_extractions([allSpectrumA, allSpectrumB])

        else:

            # STACKING A AND B SEQUENCES - ONLY IF JITTER IS NOT PRESENT
            aFrame = self.clip_and_stack(
                frames=allFrameA, recipe="soxs_nod", ignore_input_masks=False, post_stack_clipping=True
            )

            bFrame = self.clip_and_stack(
                frames=allFrameB, recipe="soxs_nod", ignore_input_masks=False, post_stack_clipping=True
            )

            # INJECT KEYWORDS INTO HEADER
            self.update_fits_keywords(frame=aFrame)
            self.update_fits_keywords(frame=bFrame)

            if self.recipeSettings["use_flat"] and master_flat:
                masterFlat = master_flat
            else:
                masterFlat = False

            mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins = self.process_single_ab_nodding_cycle(
                aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath=orderTablePath, masterFlat=master_flat
            )
            stackedSpectrum, extractionPath = self.stack_extractions([mergedSpectrumDF_A, mergedSpectrumDF_B])

            if self.generateReponseCurve:
                mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins = self.process_single_ab_nodding_cycle(
                    aFrame=aFrame,
                    bFrame=bFrame,
                    locationSetIndex=1,
                    orderTablePath=orderTablePath,
                    notFlattened=True,
                    masterFlat=master_flat,
                )
                stackedSpectrum_notflat, extractionPath_notflat = self.stack_extractions(
                    [mergedSpectrumDF_A, mergedSpectrumDF_B], notFlattened=True
                )

            if self.generateReponseCurve:

                # GETTING THE RESPONSE
                from soxspipe.commonutils import response_function

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
                )
                self.qc, self.products = response.get()

        # CHECK IF FLUX CALIBRATION IS REQUESTED
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
            filePath, products = fluxCalibrator.calibrate()
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
        )

        self.clean_up()
        self.report_output()

        self.log.debug("completed the ``produce_product`` method")

        return productPath

    def process_single_ab_nodding_cycle(
        self, aFrame, bFrame, locationSetIndex, orderTablePath, notFlattened=False, masterFlat=False
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
        - ``mergedSpectrumDF_B`` -- the order merged spectrum of nodding location B (dataframe)

        **Usage:**

        ```python
        mergedSpectrumDF_A, mergedSpectrumDF_B = soxs_nod.process_single_ab_nodding_cycle(
            aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath=orderTablePath, masterFlat=masterFlat)
        ```
        """
        self.log.debug("starting the ``process_single_ab_nodding_cycle`` method")

        from soxspipe.commonutils import horne_extraction
        import pandas as pd

        # SUBTRACTING A FROM B
        A_minus_B_notflattened = aFrame.subtract(bFrame)
        B_minus_A_notflattened = bFrame.subtract(aFrame)

        # REAPPLYING HEADERS
        hdr_A = aFrame.header
        hdr_B = bFrame.header
        A_minus_B_notflattened.header = hdr_A
        B_minus_A_notflattened.header = hdr_B

        # WRITE IN A FITS FILE THE A-B AND B-A FRAMES
        self.log.print(f"\n# PROCESSING AB NODDING CYCLE {locationSetIndex}")
        home = expanduser("~")
        filename = self.sofName + f"_AB_{locationSetIndex}.fits"
        filePath = f"{self.productDir}/{filename}"
        A_minus_B_notflattened.write(filePath, overwrite=True, checksum=True)

        filename = self.sofName + f"_BA_{locationSetIndex}.fits"
        filePath = f"{self.productDir}/{filename}"
        B_minus_A_notflattened.write(filePath, overwrite=True, checksum=True)

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
                frame=B_minus_A_notflattened,
                settings=self.settings,
                recipeName=self.recipeName,
                qcTable=self.qc,
                orderTablePath=orderTablePath,
            )

        if self.recipeSettings["save_single_frame_extractions"] == False:
            theseProducts = False
        else:
            theseProducts = self.products

        if not isinstance(masterFlat, bool) and notFlattened == False:
            A_minus_B = self.detrend(
                inputFrame=A_minus_B_notflattened,
                master_bias=False,
                dark=False,
                master_flat=masterFlat,
                order_table=orderTablePath,
            )
            B_minus_A = self.detrend(
                inputFrame=B_minus_A_notflattened,
                master_bias=False,
                dark=False,
                master_flat=masterFlat,
                order_table=orderTablePath,
            )

        else:
            A_minus_B = A_minus_B_notflattened
            B_minus_A = B_minus_A_notflattened

        # EXTRACT THE A MINUS B FRAME
        optimalExtractor = horne_extraction(
            log=self.log,
            skyModelFrame=False,
            skySubtractedFrame=A_minus_B,
            unflattenedFrame=A_minus_B_notflattened,
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
        )

        self.qc, theseProducts, mergedSpectrumDF_A, orderJoins, extractionFITSPathA = optimalExtractor.extract()

        # EXTRACT THE B MINUS A FRAME
        optimalExtractor = horne_extraction(
            log=self.log,
            skyModelFrame=False,
            skySubtractedFrame=B_minus_A,
            unflattenedFrame=B_minus_A_notflattened,
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
        )
        self.qc, theseProducts, mergedSpectrumDF_B, orderJoins, extractionFITSPathB = optimalExtractor.extract()

        if self.recipeSettings["save_single_frame_extractions"] == True:
            self.products = theseProducts

        self.log.debug("completed the ``process_single_ab_nodding_cycle`` method")
        return mergedSpectrumDF_A, mergedSpectrumDF_B, orderJoins

    def stack_extractions(self, dataFrameList, notFlattened=False):
        """*merge individual AB cycles into a master extraction*

        **Key Arguments:**

        - ``dataFrameList`` -- a list of order-merged spectrum dataframes

        **Return:**

        - ``stackedSpectrum`` -- the combined spectrum in a dataframe

        **Usage:**

        ```python
        stackedSpectrum = soxs_nod.stack_extractions(
            [mergedSpectrumDF_A, mergedSpectrumDF_B])
        ```
        """
        self.log.debug("starting the ``stack_extractions`` method")

        import pandas as pd
        from datetime import datetime
        from astropy.io import fits
        from astropy.table import Table
        import numpy as np
        from soxspipe.commonutils.toolkit import calculate_rolling_snr
        from astropy import units as u
        from specutils import Spectrum1D

        if notFlattened:
            postfix = "_NOTFLAT"
        else:
            postfix = ""

        # MERGE THE PANDAS DATAFRAMES MERDGED_ORDERS_A AND mergedSpectrumDF_B INTO A SINGLE DATAFRAME, THEN GROUP BY WAVE AND SUM THE FLUXES

        merged_dataframe = pd.concat(dataFrameList)
        # BEFORE GROUPING, WE NEED TO TRUNCATE THE WAVELENGTH TO THE 4 DIGITS
        merged_dataframe["WAVE"] = merged_dataframe["WAVE"].apply(lambda x: round(float(x.value), 4))
        groupedDataframe = merged_dataframe.groupby(by="WAVE", as_index=False).median()

        self.filenameTemplate = self.sofName + ".fits"

        # PREPARING THE HEADER
        kw = keyword_lookup(log=self.log, settings=self.settings).get

        # SELECTING HEADER A_minus_B (is this the same?)
        self.update_fits_keywords(frame=self.masterHeaderFrame)
        header = self.masterHeaderFrame.header

        header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
        header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{self.arm}".upper()

        flux_orig = groupedDataframe["FLUX_COUNTS"].values * u.electron
        spectrum_orig = Spectrum1D(
            flux=flux_orig, spectral_axis=groupedDataframe["WAVE"].values * u.nm, bin_specification="center"
        )

        groupedDataframe = calculate_rolling_snr(dataframe=groupedDataframe, flux_column="FLUX_COUNTS", window_size=300)

        # groupedDataframe['signal'] = groupedDataframe['FLUX_COUNTS'].rolling(
        #     window=15, center=True).median().fillna(method='bfill').fillna(method='ffill').values
        # groupedDataframe['normalised_flux'] = groupedDataframe['FLUX_COUNTS'] / \
        #     groupedDataframe['signal']

        # PREPARING THE HDU
        stackedSpectrum = Table.from_pandas(groupedDataframe, index=False)
        BinTableHDU = fits.table_to_hdu(stackedSpectrum)
        priHDU = fits.PrimaryHDU(header=header)
        hduList = fits.HDUList([priHDU, BinTableHDU])

        # WRITE PRODUCT TO DISK
        home = expanduser("~")
        filename = self.filenameTemplate.replace(".fits", "_EXTRACTED_MERGED" + postfix + ".fits")
        filePath = f"{self.productDir}/{filename}"
        hduList.writeto(filePath, checksum=True, overwrite=True)

        # SAVE THE TABLE stackedSpectrum TO DISK IN ASCII FORMAT
        asciiFilename = self.filenameTemplate.replace(".fits", "_EXTRACTED_MERGED" + postfix + ".txt")
        asciiFilePath = f"{self.productDir}/{asciiFilename}"
        stackedSpectrum.write(asciiFilePath, format="ascii", overwrite=True)

        utcnow = datetime.utcnow()
        self.utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        self.dateObs = header[kw("DATE_OBS")]

        self.products = pd.concat(
            [
                self.products,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "product_label": "EXTRACTED_MERGED_TABLE",
                        "file_name": filename,
                        "file_type": "FITS",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": self.utcnow,
                        "product_desc": f"Table of the extracted source in each order. All nodding cycles combined.",
                        "file_path": filePath,
                        "label": "PROD",
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

        self.products = pd.concat(
            [
                self.products,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "product_label": "EXTRACTED_MERGED_ASCII",
                        "file_name": asciiFilename,
                        "file_type": "TXT",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": self.utcnow,
                        "product_desc": f"Ascii version of extracted source spectrum",
                        "file_path": asciiFilePath,
                        "label": "PROD",
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

        self.log.debug("completed the ``stack_extractions`` method")
        return stackedSpectrum, filePath
