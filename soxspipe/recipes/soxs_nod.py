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
from soxspipe.commonutils.toolkit import generic_quality_checks, spectroscopic_image_quality_checks
from fundamentals import tools
from builtins import object
import sys
import os
from soxspipe.commonutils.filenamer import filenamer
from os.path import expanduser

os.environ['TERM'] = 'vt100'

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
            overwrite=False

    ):
        # INHERIT INITIALISATION FROM  base_recipe
        super(soxs_nod, self).__init__(
            log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-nod")
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
            ext=self.settings['data-extension']
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
        self.inputFrames.sort(['MJD-OBS'])
        if self.verbose:
            self.log.print("# RAW INPUT FRAMES - SUMMARY")
            self.log.print(self.inputFrames.summary, "\n")

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])

        return None

    def verify_input_frames(
            self):
        """*verify the input frame match those required by the soxs_nod recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

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
            for i in [f"DISP_TAB_{self.arm}", f'ORDER_TAB_{self.arm}', f'DISP_IMAGE_{self.arm}']:
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
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
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
        self.log.debug('starting the ``produce_product`` method')

        from astropy.nddata import CCDData
        from astropy import units as u
        import pandas as pd
        from datetime import datetime
        from soxspipe.commonutils.toolkit import quicklook_image

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None

        # OBJECT/STANDARD FRAMES
        types = ['OBJECT', 'STD,FLUX', 'STD,TELLURIC']
        allObjectFrames = []
        self.masterHeader = False
        for t in types:
            add_filters = {kw("DPR_TYPE"): t,
                           kw("DPR_TECH"): 'ECHELLE,SLIT,NODDING'}
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)
                if not self.masterHeader:
                    self.masterHeader = singleFrame.header
            if len(allObjectFrames):
                break

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
        self.dispMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # FIND THE 2D MAP IMAGE
        filterDict = {kw("PRO_CATG"): f"DISP_IMAGE_{arm}"}
        self.twoDMap = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # DETREND ALL NODDING IMAGES
        if self.recipeSettings["use_flat"]:
            allObjectFrames[:] = [self.detrend(inputFrame=f, master_bias=False, dark=False, master_flat=master_flat, order_table=orderTablePath) for f in allObjectFrames]

        quicklook_image(
            log=self.log, CCDObject=allObjectFrames[0], show=False, ext='data', stdWindow=3, title=False, surfacePlot=True, saveToPath=False)

        # DIVIDING IN A AND B SEQUENCES
        allFrameA, allFrameB, allFrameAOffsets, allFrameBOffsets = [], [], [], []

        # CUMOFF Y IS THE OFFSET IN THE Y DIRECTION OF THE NODDING SEQUENCE. POSITIVE A, NEGATIVE B
        for frame in allObjectFrames:
            offset = frame.header[kw("NOD_CUMULATIVE_OFFSETY")]
            if offset > 0:
                allFrameAOffsets.append(offset)
                allFrameA.append(frame)
            else:
                allFrameBOffsets.append(offset)
                allFrameB.append(frame)

        uniqueOffsets = list(set(allFrameAOffsets))
        if len(uniqueOffsets) > 1:
            s = "S"
        else:
            s = ""
        self.log.print(f"# PROCESSING {len(allFrameAOffsets)} AB NODDING CYCLES WITH {len(uniqueOffsets)} UNIQUE PAIR{s} OF OFFSET LOCATIONS")

        if len(allFrameAOffsets) > 1 and len(uniqueOffsets) > 1:
            allSpectrumA = []
            allSpectrumB = []
            sequenceCount = 1
            # SORT frameA and frameB looing at their MJDOBS keyword in the header in order to the closest A and B frames in time
            allFrameA.sort(key=lambda x: x.header["MJD-OBS"])
            allFrameB.sort(key=lambda x: x.header["MJD-OBS"])

            for frameA, frameB in zip(allFrameA, allFrameB):
                self.log.print(f"Processing AB Nodding Sequence {sequenceCount}")
                if False:
                    quicklook_image(log=self.log, CCDObject=frameA, show=False, ext='data', stdWindow=1, title=False, surfacePlot=False, saveToPath=False)
                    quicklook_image(log=self.log, CCDObject=frameB, show=False, ext='data', stdWindow=1, title=False, surfacePlot=False, saveToPath=False)
                    # Save frameA and frameB to disk in temporary file
                    home = expanduser("~")
                    filenameA = self.sofName + f"_A_{sequenceCount}.fits"
                    filenameB = self.sofName + f"_B_{sequenceCount}.fits"
                    outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"
                    filePathA = f"{outDir}/{filenameA}"
                    filePathB = f"{outDir}/{filenameB}"
                    frameA.write(filePathA, overwrite=True)
                    frameB.write(filePathB, overwrite=True)

                # PROCESSING SINGLE SEQUENCE
                mergedSpectrumDF_A, mergedSpectrumDF_B = self.process_single_ab_nodding_cycle(aFrame=frameA, bFrame=frameB, locationSetIndex=sequenceCount, orderTablePath=orderTablePath)
                if sequenceCount == 1:
                    allSpectrumA = mergedSpectrumDF_A
                    allSpectrumB = mergedSpectrumDF_B
                else:
                    allSpectrumA = pd.concat([allSpectrumA, mergedSpectrumDF_A])
                    allSpectrumB = pd.concat([allSpectrumB, mergedSpectrumDF_B])

                sequenceCount += 1
            stackedSpectrum = self.stack_extractions([allSpectrumA, allSpectrumB])
            self.plot_stacked_spectrum_qc(stackedSpectrum)
            self.clean_up()
            self.report_output()
            # sys.exit(0)
        else:
            # STACKING A AND B SEQUENCES - ONLY IF JITTER IS NOT PRESENT
            aFrame = self.clip_and_stack(
                frames=allFrameA,
                recipe="soxs_nod",
                ignore_input_masks=False,
                post_stack_clipping=True)

            bFrame = self.clip_and_stack(
                frames=allFrameB,
                recipe="soxs_nod",
                ignore_input_masks=False,
                post_stack_clipping=True)

            mergedSpectrumDF_A, mergedSpectrumDF_B = self.process_single_ab_nodding_cycle(aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath=orderTablePath)
            stackedSpectrum = self.stack_extractions([mergedSpectrumDF_A, mergedSpectrumDF_B])
            self.plot_stacked_spectrum_qc(stackedSpectrum)

            self.clean_up()
            self.report_output()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    def process_single_ab_nodding_cycle(
            self,
            aFrame,
            bFrame,
            locationSetIndex,
            orderTablePath):
        """*process a single AB nodding cycle*

        **Key Arguments:**

        - ``aFrame`` -- the frame taken at the A location. CCDData object.
        - ``bFrame`` -- the frame taken at the B location. CCDDate object.
        - ``locationSetIndex`` -- the index of the AB cycle
        - `orderTablePath` -- path to the order table

        **Return:**

        - ``mergedSpectrumDF_A`` -- the order merged spectrum of nodding location A (dataframe)
        - ``mergedSpectrumDF_B`` -- the order merged spectrum of nodding location B (dataframe)

        **Usage:**

        ```python
        mergedSpectrumDF_A, mergedSpectrumDF_B = soxs_nod.process_single_ab_nodding_cycle(aFrame=aFrame, bFrame=bFrame, locationSetIndex=1)
        ```
        """
        self.log.debug('starting the ``process_single_ab_nodding_cycle`` method')

        from soxspipe.commonutils import horne_extraction
        import pandas as pd

        # SUBTRACTING A FROM B
        A_minus_B = aFrame.subtract(bFrame)
        B_minus_A = bFrame.subtract(aFrame)

        # REAPPLYING HEADERS
        hdr_A = aFrame.header
        hdr_B = bFrame.header
        A_minus_B.header = hdr_A
        B_minus_A.header = hdr_B

        # Write in a fits file the A-B and B-A frames
        home = expanduser("~")
        filename = self.sofName + f"_AB_{locationSetIndex}.fits"
        outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"
        filePath = f"{outDir}/{filename}"
        A_minus_B.write(filePath, overwrite=True)

        filename = self.sofName + f"_BA_{locationSetIndex}.fits"
        B_minus_A.write(filePath, overwrite=True)

        if False:
            from soxspipe.commonutils.toolkit import quicklook_image
            quicklook_image(
                log=self.log, CCDObject=A_minus_B, show=True, ext='data', stdWindow=3, title=False, surfacePlot=True, saveToPath=False)
            quicklook_image(
                log=self.log, CCDObject=B_minus_A, show=True, ext='data', stdWindow=3, title=False, surfacePlot=True, saveToPath=False)

        # TODO: ADD THESE CHECKS .... LIKELY FOR EACH AB CYCLE INDEX
        if True:
            self.qc = generic_quality_checks(
                log=self.log, frame=A_minus_B, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)
            self.qc = spectroscopic_image_quality_checks(
                log=self.log, frame=A_minus_B, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc, orderTablePath=orderTablePath)

        if self.recipeSettings["save_single_frame_extractions"] == False:
            theseProducts = False
        else:
            theseProducts = self.products

        # EXTRACT THE A MINUS B FRAME
        optimalExtractor = horne_extraction(
            log=self.log,
            skyModelFrame=False,
            skySubtractedFrame=A_minus_B,
            twoDMapPath=self.twoDMap,
            settings=self.settings,
            recipeName=self.recipeName,
            recipeSettings=self.recipeSettings,
            qcTable=self.qc,
            productsTable=theseProducts,
            dispersionMap=self.dispMap,
            sofName=self.sofName,
            locationSetIndex=locationSetIndex
        )

        self.qc, theseProducts, mergedSpectrumDF_A = optimalExtractor.extract()

        # EXTRACT THE B MINUS A FRAME
        optimalExtractor = horne_extraction(
            log=self.log,
            skyModelFrame=False,
            skySubtractedFrame=B_minus_A,
            twoDMapPath=self.twoDMap,
            settings=self.settings,
            recipeName=self.recipeName,
            recipeSettings=self.recipeSettings,
            qcTable=self.qc,
            productsTable=theseProducts,
            dispersionMap=self.dispMap,
            sofName=self.sofName,
            locationSetIndex=locationSetIndex
        )
        self.qc, theseProducts, mergedSpectrumDF_B = optimalExtractor.extract()

        if self.recipeSettings["save_single_frame_extractions"] == True:
            self.products = theseProducts

        self.log.debug('completed the ``process_single_ab_nodding_cycle`` method')
        return mergedSpectrumDF_A, mergedSpectrumDF_B

    def stack_extractions(
            self,
            dataFrameList):
        """*merge individual AB cycles into a master extraction*

        **Key Arguments:**

        - ``dataFrameList`` -- a list of order-merged spectrum dataframes

        **Return:**

        - ``stackedSpectrum`` -- the combined spectrum in a dataframe

        **Usage:**

        ```python
        stackedSpectrum = soxs_nod.stack_extractions([mergedSpectrumDF_A, mergedSpectrumDF_B])
        ```
        """
        self.log.debug('starting the ``stack_extractions`` method')

        import pandas as pd
        from datetime import datetime
        from astropy.io import fits
        from astropy.table import Table

        # MERGE THE PANDAS DATAFRAMES MERDGED_ORDERS_A AND mergedSpectrumDF_B INTO A SINGLE DATAFRAME, THEN GROUP BY WAVE AND SUM THE FLUXES

        merged_dataframe = pd.concat(dataFrameList)
        # BEFORE GROUPING, WE NEED TO TRUNCATE THE WAVELENGHT TO THE 4 DIGITS
        merged_dataframe['WAVE'] = merged_dataframe['WAVE'].apply(lambda x: round(float(x.value), 4))
        groupedDataframe = merged_dataframe.groupby(by='WAVE', as_index=False).median()

        self.filenameTemplate = self.sofName + ".fits"

        # PREPARING THE HEADER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get

        # SELECTING HEADER A_minus_B (is this the same?)
        header = self.masterHeader
        header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
        header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{self.arm}".upper()

        # PREPARING THE HDU
        stackedSpectrum = Table.from_pandas(groupedDataframe, index=False)
        BinTableHDU = fits.table_to_hdu(stackedSpectrum)
        priHDU = fits.PrimaryHDU(header=header)
        hduList = fits.HDUList([priHDU, BinTableHDU])

        # WRITE PRODUCT TO DISK
        home = expanduser("~")
        filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED.fits")
        outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"
        filePath = f"{outDir}/{filename}"
        hduList.writeto(filePath, checksum=True, overwrite=True)

        # SAVE THE TABLE stackedSpectrum TO DISK IN ASCII FORMAT
        asciiFilename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED.txt")
        asciiFilePath = f"{outDir}/{asciiFilename}"
        stackedSpectrum.write(asciiFilePath, format='ascii', overwrite=True)

        utcnow = datetime.utcnow()
        self.utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        self.dateObs = self.masterHeader[kw("DATE_OBS")]

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "EXTRACTED_MERGED_TABLE",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": self.utcnow,
            "product_desc": f"Table of the extracted source in each order. All nodding cycles combined.",
            "file_path": filePath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "EXTRACTED_MERGED_ASCII",
            "file_name": asciiFilename,
            "file_type": "TXT",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": self.utcnow,
            "product_desc": f"Ascii version of extracted source spectrum",
            "file_path": filePath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``stack_extractions`` method')
        return stackedSpectrum

    def plot_stacked_spectrum_qc(
            self,
            stackedSpectrum):
        """*plot stacked spectrum as a QC plot*

        **Key Arguments:**

        - ``stackedSpectrum`` -- the combined spectrum in a dataframe

        **Usage:**

        ```python
        soxs_nod.plot_stacked_spectrum_qc(stackedSpectrum)
        ```
        """
        self.log.debug('starting the ``plot_stacked_spectrum_qc`` method')

        import matplotlib.pyplot as plt
        from astropy.stats import sigma_clipped_stats
        import pandas as pd

        fig = plt.figure(figsize=(16, 2), constrained_layout=True, dpi=320)
        gs = fig.add_gridspec(1, 1)
        toprow = fig.add_subplot(gs[0:1, :])
        toprow.legend(loc='lower right', bbox_to_anchor=(1, -0.5),
                      fontsize=8)
        toprow.set_ylabel('flux ($e^{-}$)', fontsize=10)
        toprow.set_xlabel(f'wavelength (nm)', fontsize=10)
        toprow.set_title(
            f"Optimally Extracted Order-Merged Object Spectrum ({self.arm.upper()})", fontsize=11)

        mean, median, std = sigma_clipped_stats(stackedSpectrum['FLUX_COUNTS'], sigma=5.0, stdfunc="mad_std", cenfunc="median", maxiters=3)
        plt.plot(stackedSpectrum['WAVE'], stackedSpectrum['FLUX_COUNTS'], linewidth=0.2, color="#dc322f")
        plt.ylim(-200, median + 20 * std)
        plt.xlim(stackedSpectrum['WAVE'].min(), stackedSpectrum['WAVE'].max())

        filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED_QC_PLOT.pdf")
        filePath = f"{self.qcDir}/{filename}"
        # plt.show()
        plt.savefig(filePath, dpi='figure')

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": "soxs-stare",
            "product_label": f"EXTRACTED_MERGED_QC_PLOT",
            "file_name": filename,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": self.utcnow,
            "product_desc": f"QC plot of extracted order-merged source. All nodding cycles combined.",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``plot_stacked_spectrum_qc`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
