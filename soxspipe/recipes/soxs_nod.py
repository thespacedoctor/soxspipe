#!/usr/bin/env python
# encoding: utf-8
"""
*Reduce SOXS data taken in nodding mode*

:Author:
    David Young & Marco Landoni

:Date Created:
    February 27, 2024
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
from ._base_recipe_ import _base_recipe_
from soxspipe.commonutils.toolkit import generic_quality_checks, spectroscopic_image_quality_checks
from fundamentals import tools
from builtins import object
import sys
import os
from matplotlib import pyplot as plt  
from soxspipe.commonutils.filenamer import filenamer
from os.path import expanduser
from astropy.io import fits
from astropy.table import Table
os.environ['TERM'] = 'vt100'

class soxs_nod(_base_recipe_):
    """
    *The soxs_nod recipe*

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
        - add a tutorial about ``soxs_nod`` to documentation
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
        super(soxs_nod, self).__init__(
            log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-nod")
        self.log = log
        log.debug("instansiating a new 'soxs_nod' object")
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
            inputFrames=self.inputFrames,
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

        **Return:**
            - ``None``

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        error = False

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()
        arm = self.arm

        if self.arm == "NIR":
            if not error:
                for i in imageTypes:
                    if i not in ["OBJECT", "LAMP,FLAT", "DARK", "STD,FLUX"]:
                        error = f"Found a {i} file. Input frames for soxspipe nod need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark (MASTER_DARK_{arm}) or off-frame for NIR."

            if not error:
                for i in imageTech:
                    if i not in ["IMAGE", "ECHELLE,SLIT", "ECHELLE,MULTI-PINHOLE", "ECHELLE,SLIT,NODDING"]:
                        error = f"Input frames for soxspipe nod need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-flat (MASTER_FLAT_{arm}) and master dark (MASTER_DARK_{arm}) or off-frame for NIR. The sof file is missing a {i} frame."

        else:
            if not error:
                for i in imageTypes:
                    if i not in ["OBJECT", "LAMP,FLAT"]:
                        error = f"Input frames for soxspipe nod need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} frame."

            if not error:
                for i in [f"DISP_TAB_{self.arm}"]:
                    if i not in imageCat:
                        error = f"Input frames for soxspipe nod need to be an object frame (OBJECT_{arm}), a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), a master-flat (MASTER_FLAT_{arm}) and optionally a master dark (MASTER_DARK_{arm}) for UVB/VIS. The sof file is missing a {i} frame."

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

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None

        # OBJECT FRAMES
        add_filters = {kw("DPR_TYPE"): 'OBJECT',
                       kw("DPR_TECH"): 'ECHELLE,SLIT,NODDING'}
        allObjectFrames = []
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            singleFrame = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
            allObjectFrames.append(singleFrame)

        # FLUX STD FRAMES
        if not len(allObjectFrames):
            add_filters = {kw("DPR_TYPE"): 'STD,FLUX',
                           kw("DPR_TECH"): 'ECHELLE,SLIT,NODDING'}
            allObjectFrames = []
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                singleFrame = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                           hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
                allObjectFrames.append(singleFrame)

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

        # FLAT-FIELD ALL NOD IMAGES
        if False:
            allObjectFrames[:] = [self.detrend(inputFrame=f, master_bias=False, dark=False, master_flat=master_flat, order_table=orderTablePath) for f in allObjectFrames]
        else:
            allObjectFrames[:] = [self.detrend(inputFrame=f, master_bias=False, dark=False, master_flat=master_flat) for f in allObjectFrames]


            

        #DIVIDING IN A AND B SEQUENCES
        allFrameA = []
        allFrameB = []


        for frame in allObjectFrames:
            if frame.header['HIERARCH ESO SEQ CUMOFF Y'] > 0:   
                allFrameA.append(frame)

            else:
                allFrameB.append(frame)



        #STACKING A AND B SEQUENCES
        masterA = self.clip_and_stack(
                frames= allFrameA, 
                recipe = "soxs_nod", 
                ignore_input_masks =False, 
                post_stack_clipping= True)
        
        masterB = self.clip_and_stack( 
            frames = allFrameB, 
            recipe = "soxs_nod", 
            ignore_input_masks = False,
            post_stack_clipping = True)

        #SUBTRACTING A FROM B
        A_minus_B = masterA.subtract(masterB)
        B_minus_A = masterB.subtract(masterA)

        #ADD TO A_minus_B header of masterA

        hdr = masterB.header
        A_minus_B.header = hdr
        B_minus_A.header = hdr  

        

        #Write the A-B and B-A frames to disk
        A_minus_B_path = "A_minus_B.fits"
        B_minus_A_path =  "B_minus_A.fits"

        A_minus_B.write(A_minus_B_path, overwrite=True)
        B_minus_A.write(B_minus_A_path, overwrite=True)


        # ADD QUALITY CHECKS
        # TODO: ADD THESE CHECK WHEN WE HAVE A FINAL FRAME TO CHECK ... LIKELY A-B FRAME
        if False:
            self.qc = generic_quality_checks(
                log=self.log, frame=mflat, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)
            self.qc = spectroscopic_image_quality_checks(
                log=self.log, frame=mflat, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc, orderTablePath=orderTablePath)
        from soxspipe.commonutils import horne_extraction
        
        #EXTRACT THE A MINUS B FRAME
        
        optimalExtractor = horne_extraction(
             log=self.log,
             skyModelFrame=False,
             skySubtractedFrame=A_minus_B,
             twoDMapPath=twoDMap,
             settings=self.settings,
             recipeName=self.recipeName,
             qcTable=self.qc,
             productsTable=self.products,
             dispersionMap=dispMap,
             sofName=self.sofName
             
         )
        
        self.qc, self.products, merged_orders_a = optimalExtractor.extract(noddingSequence='A')
        #EXTRACT THE B MINUS A FRAME 

        optimalExtractor = horne_extraction(
             log=self.log,
             skyModelFrame=False,
             skySubtractedFrame=B_minus_A,
             twoDMapPath=twoDMap,
             settings=self.settings,
             recipeName=self.recipeName,
             qcTable=self.qc,
             productsTable=self.products,
             dispersionMap=dispMap,
             sofName=self.sofName
         )
        self.qc, self.products, merged_orders_b = optimalExtractor.extract(noddingSequence='B')    

        #MERGE THE PANDAS DATAFRAMES MERDGED_ORDERS_A AND MERGED_ORDERS_B INTO A SINGLE DATAFRAME, THEN GROUP BY WAVE AND SUM THE FLUXES

        merged_dataframe = pd.concat([merged_orders_a, merged_orders_b])

        groupedDataframe = merged_dataframe.groupby(by='WAVE', as_index=False).sum()

        print(groupedDataframe)

        #WRITING THE DATA ON THE FITS FILE
        if self.sofName:
            self.filenameTemplate = self.sofName + ".fits"
        else:
            self.filenameTemplate = filenamer(
                log=self.log,
                frame=self.skySubtractedFrame,
                settings=self.settings
            )
        
        #PREPARING THE HEADER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get

        #SELECTING HEADER A_minus_B (is this the same?)
        header = A_minus_B.header
        header[kw("SEQ_ARM")] = A_minus_B.header[kw("SEQ_ARM")]
        header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
        header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{arm}".upper()
        
        #PREPARING THE HDU
        groupedDataframeTable = Table.from_pandas(groupedDataframe, index=False)
        BinTableHDU = fits.table_to_hdu(groupedDataframeTable)

        priHDU = fits.PrimaryHDU(header=header)
        hduList = fits.HDUList([priHDU, BinTableHDU])

        #WRTIE PRODUCT TO DISK
        home = expanduser("~")
        filename = self.filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED_AB.fits")
        outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"
        filePath = f"{outDir}/{filename}"
        hduList.writeto(filePath, checksum=True, overwrite=True)
        

        self.clean_up()
        self.report_output()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
