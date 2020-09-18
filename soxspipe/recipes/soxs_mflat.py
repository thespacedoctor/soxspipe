#!/usr/bin/env python
# encoding: utf-8
"""
*generate a single normalised master flat-field frame*

:Author:
    David Young & Marco Landoni

:Date Created:
    September 16, 2020
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
import ccdproc
from soxspipe.commonutils import keyword_lookup
import matplotlib.pyplot as plt
from soxspipe.commonutils.toolkit import unpack_order_table
import numpy.ma as ma
from soxspipe.commonutils.toolkit import quicklook_image


class soxs_mflat(_base_recipe_):
    """
    *The soxs_mflat recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.

    **Usage**

    ```python
    from soxspipe.recipes import soxs_mflat
    recipe = soxs_mflat(
        log=log,
        settings=settings,
        inputFrames=fileList
    )
    mflatFrame = recipe.produce_product()
    ```

    ---

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``soxs_mflat`` to documentation
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
        super(soxs_mflat, self).__init__(
            log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mflat' object")
        self.settings = settings
        self.inputFrames = inputFrames
        # xt-self-arg-tmpx

        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_mflat - NO MORE, NO LESS.
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
        """*verify the input frame match those required by the soxs_mflat recipe*

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

            if imageTypes[0] != "LAMP,FLAT":
                raise TypeError(
                    "Input frames for soxspipe mflat need to be slit flat lamp on and lamp off frames for NIR" % locals())

            for i in imageTech:
                if i not in ['ECHELLE,SLIT', 'IMAGE']:
                    raise TypeError(
                        "Input frames for soxspipe mflat need to be slit flat lamp on and lamp off frames for NIR" % locals())

        else:
            for i in imageTypes:
                if i not in ["LAMP,FLAT", "LAMP,QFLAT", "LAMP,DFLAT", "BIAS", "DARK"]:
                    raise TypeError(
                        "Input frames for soxspipe mflat need to be slit flat lamp frames and a master-bias and possibly a master dark for UVB/VIS" % locals())

        # LOOK FOR ORDER CENTRE TABLE
        arm = self.arm
        if arm not in self.supplementaryInput or "ORDER_CENT" not in self.supplementaryInput[arm]:
            raise TypeError(
                "Need an order centre for %(arm)s - none found with the input files" % locals())

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*The code to generate the product of the soxs_mflat recipe*

        **Return:**
            - ``productPath`` -- the path to the final product
        """
        self.log.debug('starting the ``produce_product`` method')

        productPath = None
        arm = self.arm

        calibratedFlats = self.calibrate_frame_set()

        normalisedFlats = self.normalise_flats(
            calibratedFlats, orderTablePath=self.supplementaryInput[arm]["ORDER_CENT"])

        combined_normalised_flat = self.clip_and_stack(
            frames=normalisedFlats, recipe="soxs_mflat")

        quicklook_image(
            log=self.log, CCDObject=combined_normalised_flat, show=False)

        if 1 == 1:
            from os.path import expanduser
            home = expanduser("~")
            outDir = self.settings["intermediate-data-root"].replace("~", home)
            filePath = f"{outDir}/first_iteration_{arm}_master_flat.fits"
            self._write(combined_normalised_flat, filePath, overwrite=True)

        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    def calibrate_frame_set(
            self):
        """*given all of the input data calibrate the frames by subtracting bias and/or dark*

        **Return:**
            - ``calibratedFlats`` -- the calibrated frames
        """
        self.log.debug('starting the ``calibrate_frame_set`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        # FIND THE BIAS FRAMES
        try:
            filterDict = {kw("DPR_TYPE").lower(): "BIAS"}
            biasCollection = self.inputFrames.filter(**filterDict)
            # LIST OF CCDDATA OBJECTS
            biases = [c for c in biasCollection.ccds(ccd_kwargs={
                "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]
            bias = biases[0]
        except:
            bias = None
            biasCollection = None

        # FIND THE DARK FRAMES
        try:
            filterDict = {kw("DPR_TYPE").lower(): "DARK"}
            darkCollection = self.inputFrames.filter(**filterDict)
        except:
            darkCollection = None

        if not darkCollection:
            try:
                filterDict = {kw("DPR_TYPE").lower(): "LAMP,FLAT",
                              kw("DPR_TECH").lower(): "IMAGE"}
                darkCollection = self.inputFrames.filter(**filterDict)
            except:
                darkCollection = None

        # FIND THE FLAT FRAMES
        try:
            filterDict = {kw("DPR_TYPE").lower(): "LAMP,(D|Q)?FLAT",
                          kw("DPR_TECH").lower(): "ECHELLE,SLIT"}
            flatCollection = self.inputFrames.filter(**filterDict)
        except:
            raise FileNotFoundError(
                "The mflat recipe needs flat-frames as input, none found")
        # LIST OF CCDDATA OBJECTS
        flats = [c for c in flatCollection.ccds(ccd_kwargs={
                                                "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        # IF NO DARK FRAMES EXIST - JUST A MASTER BIAS. SUBTRACT BIAS.
        calibratedFlats = []
        if not darkCollection and bias:
            for flat in flats:
                calibratedFlats.append(self.subtract_calibrations(
                    inputFrame=flat, master_bias=bias, dark=None))

        # IF DARKS EXIST - FIND CLOSEST IN TIME TO FLAT-FRAME. SUBTRACT BIAS
        # AND?OR DARK
        if darkCollection:
            darkMjds = [h[kw("MJDOBS").lower()]
                        for h in darkCollection.headers()]
            darks = [c for c in darkCollection.ccds(ccd_kwargs={
                "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]
            for flat in flats:
                mjd = flat.header[kw("MJDOBS").lower()]
                matchValue, matchIndex = nearest_neighbour(
                    flat.header[kw("MJDOBS").lower()], darkMjds)
                dark = darks[matchIndex]
                calibratedFlats.append(self.subtract_calibrations(
                    inputFrame=flat, master_bias=bias, dark=dark))

        for frame in calibratedFlats:
            quicklook_image(log=self.log, CCDObject=frame, show=False)

        if 1 == 0:
            from os.path import expanduser
            home = expanduser("~")
            outDir = self.settings["intermediate-data-root"].replace("~", home)
            index = 1
            for frame in calibratedFlats:
                filePath = f"{outDir}/{index:02}_flat_{arm}_calibrated.fits"
                index += 1
                self._write(frame, filePath, overwrite=True)

        self.log.debug('completed the ``calibrate_frame_set`` method')
        return calibratedFlats

    def normalise_flats(
            self,
            inputFlats,
            orderTablePath):
        """*determine the median exposure for each flat frame and normalise the flux to that level*

        **Key Arguments:**
            - ``inputFlats`` -- the input flat field frames
            - ``orderTablePath`` -- path to the order table

        **Return:**
            - ``normalisedFrames`` -- the normalised flat-field frames (CCDData array)
        ```
        """
        self.log.debug('starting the ``normalise_flats`` method')

        window = int(self.settings[
            "soxs-mflat"]["centre-median-window"] / 2)

        # UNPACK THE ORDER TABLE

        orderCentres = unpack_order_table(
            log=self.log, orderTablePath=orderTablePath)
        mask = np.ones_like(inputFlats[0].data)
        for xcoords, ycoords in orderCentres:
            for x, y in zip(xcoords, ycoords):
                x = int(x)
                mask[y][x - window:x + window] = 0

        normalisedFrames = []

        # PLOT ONE OF THE MASKED FRAMES TO CHECK
        for frame in [inputFlats[0]]:
            maskedFrame = ma.array(frame.data, mask=mask)
            quicklook_image(log=self.log, CCDObject=maskedFrame, show=False)

        for frame in inputFlats:
            maskedFrame = ma.array(frame.data, mask=mask)
            median = np.ma.median(maskedFrame)
            normalisedFrame = frame.divide(median)
            normalisedFrames.append(normalisedFrame)

        # PLOT ONE OF THE NORMALISED FRAMES TO CHECK
        for frame in normalisedFrames:
            quicklook_image(log=self.log, CCDObject=frame, show=False)

        self.log.debug('completed the ``normalise_flats`` method')
        return normalisedFrames

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx


def nearest_neighbour(singleValue, listOfValues):
    arrayOfValues = np.asarray(listOfValues)
    dist = np.square(arrayOfValues - singleValue)
    minDist = np.amin(dist)
    minIndex = np.where(dist == minDist)[0][0]
    matchValue = listOfValues[minIndex]
    return matchValue, minIndex


# use the tab-trigger below for new function
# xt-def-function
