#!/usr/bin/env python
# encoding: utf-8
"""
*Recipe to generate a first approximation of the dispersion solution from single pinhole frames*

:Author:
    David Young & Marco Landoni

:Date Created:
    August 25, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils import set_of_files
from ._base_recipe_ import _base_recipe_
import numpy as np
from astropy.nddata import CCDData
from astropy import units as u
import ccdproc
from soxspipe.commonutils import keyword_lookup


class soxs_disp_solution(_base_recipe_):
    """
    *The soxs_disp_solution recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.   

    See `produce_product` method for usage.

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``soxs_disp_solution`` to documentation
    ```
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[]

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_disp_solution, self).__init__(
            log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_disp_solution' object")
        self.settings = settings
        self.inputFrames = inputFrames
        # xt-self-arg-tmpx

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames
        )
        self.inputFrames = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_disp_solution - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        print("\n# RAW INPUT DARK FRAMES - SUMMARY")
        # SORT IMAGE COLLECTION
        self.inputFrames.sort(['mjd-obs'])
        print(self.inputFrames.summary, "\n")

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])

        return None

    def verify_input_frames(
            self):
        """*verify the input frame match those required by the soxs_disp_solution recipe*

        **Return:**
            - ``None``

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        self._verify_input_frames_basics()

        imageTypes = self.inputFrames.values(
            keyword=kw("DPR_TYPE").lower(), unique=True)
        imageTech = self.inputFrames.values(
            keyword=kw("DPR_TECH").lower(), unique=True)
        imageCat = self.inputFrames.values(
            keyword=kw("DPR_CATG").lower(), unique=True)

        if self.arm == "NIR":
            # WANT ON AND OFF PINHOLE FRAMES
            # MIXED INPUT IMAGE TYPES ARE BAD
            if len(imageTypes) > 1:
                imageTypes = " and ".join(imageTypes)
                print(self.inputFrames.summary)
                raise TypeError(
                    "Input frames are a mix of %(imageTypes)s" % locals())

            if imageTypes[0] != "LAMP,FMTCHK":
                raise TypeError(
                    "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals())

            for i in imageTech:
                if i not in ['ECHELLE,PINHOLE', 'IMAGE']:
                    raise TypeError(
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals())

        else:
            for i in imageTypes:
                if i not in ["LAMP,FMTCHK", "BIAS", "DARK"]:
                    raise TypeError(
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on and a master-bias and possibly a master dark for UVB/VIS" % locals())

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*The code to generate the product of the soxs_disp_solution recipe*

        **Return:**
            - ``productPath`` -- the path to the final product

        **Usage**

        ```python
        from soxspipe.recipes import soxs_disp_solution
        recipe = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        disp_solutionFrame = recipe.produce_product()
        ```
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None

        self.inputFrames.summary.pprint_all()

        master_bias = False
        dark = False
        pinhole_image = False

        add_filters = {kw("DPR_CATG"): 'MASTER_BIAS_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        add_filters = {kw("DPR_CATG"): 'MASTER_DARK_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        add_filters = {kw("DPR_TYPE"): 'LAMP,FMTCHK', kw("DPR_TECH"): 'IMAGE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            print(i)
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        add_filters = {kw("DPR_TYPE"): 'LAMP,FMTCHK',
                       kw("DPR_TECH"): 'ECHELLE,PINHOLE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            pinhole_image = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                         hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        pinhole = self.subtract_calibrations(
            inputFrame=pinhole_image, master_bias=master_bias, dark=dark)

        if self.settings["save-intermediate-products"]:
            outDir = self.intermediateRootPath
            filePath = f"{outDir}/single_pinhole_{arm}_calibrated.fits"
            print(f"Calibrated single pinhole frame: {filePath}")
            self._write(pinhole, filePath, overwrite=True)

        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    def subtract_calibrations(
            self,
            inputFrame,
            master_bias=False,
            dark=False):
        """*subtract calibration frames from an input frame*

        **Key Arguments:**
            - ``inputFrame`` -- the input frame to have calibrations subtracted. CCDData object.
            - ``master_bias`` -- the master bias frame to be subtracted. CCDData object. Default *False*.
            - ``dark`` -- a dark frame to be subtracted. CCDData object. Default *False*.

        **Return:**
            - ``calibration_subtracted_frame`` -- the input frame with the calibration frame(s) subtracted. CCDData object.

        **Usage:**

        ```python
        usage code 
        ```

        ---

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - write a command-line tool for this method
            - update package tutorial with command-line tool info if needed
        ```
        """
        self.log.debug('starting the ``subtract_calibrations`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        # VERIFY DATA IS IN ORDER
        if master_bias == False and dark == False:
            raise TypeError(
                "subtract_calibrations method needs a master-bias frame and/or a dark frame to subtract")
        if master_bias == False and dark.header[kw("EXPTIME")] != inputFrame.header[kw("EXPTIME")]:
            raise AttributeError(
                "Dark and science/calibration frame have differing exposure-times. A master-bias frame needs to be supplied to scale the dark frame to same exposure time as input science/calibration frame")
        if master_bias != False and dark != False and dark.header[kw("EXPTIME")] != inputFrame.header[kw("EXPTIME")]:
            raise AttributeError(
                "CODE NEEDS WRITTEN HERE TO SCALE DARK FRAME TO EXPOSURE TIME OF SCIENCE/CALIBRATION FRAME")

        if dark != False and dark.header[kw("EXPTIME")] == inputFrame.header[kw("EXPTIME")]:
            calibration_subtracted_frame = inputFrame.subtract(dark)
            calibration_subtracted_frame.header = inputFrame.header
            try:
                calibration_subtracted_frame.wcs = inputFrame.wcs
            except:
                pass

        if dark == False and master_bias != False:
            calibration_subtracted_frame = inputFrame.subtract(master_bias)
            calibration_subtracted_frame.header = inputFrame.header
            try:
                calibration_subtracted_frame.wcs = inputFrame.wcs
            except:
                pass

        self.log.debug('completed the ``subtract_calibrations`` method')
        return calibration_subtracted_frame

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
