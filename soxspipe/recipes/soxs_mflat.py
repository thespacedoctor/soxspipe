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
from soxspipe.commonutils import detect_order_edges
import pandas as pd
from soxspipe.commonutils import subtract_background
from os.path import expanduser
from soxspipe.commonutils.filenamer import filenamer
from astropy.stats import sigma_clip, mad_std


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

        - add a tutorial about ``soxs_mflat`` to documentation
    ```
    """

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
        self.recipeSettings = settings["soxs-mflat"]
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY soxs_mflat - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        print("\n# RAW INPUT FRAMES - SUMMARY")
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
        if arm not in self.supplementaryInput or "ORDER_LOCATIONS" not in self.supplementaryInput[arm]:
            raise TypeError(
                "Need an order centre for %(arm)s - none found with the input files" % locals())

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*generate the master flat frames updated order location table (with egde detection)*

        **Return:**
            - ``productPath`` -- the path to the master flat frame
        """
        self.log.debug('starting the ``produce_product`` method')

        productPath = None
        arm = self.arm

        # CALIBRATE THE FRAMES BY SUBTRACTING BIAS AND/OR DARK
        calibratedFlats = self.calibrate_frame_set()

        quicklook_image(
            log=self.log, CCDObject=calibratedFlats[0], show=False)

        # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
        # FLUX TO THAT LEVEL
        normalisedFlats = self.normalise_flats(
            calibratedFlats, orderTablePath=self.supplementaryInput[arm]["ORDER_LOCATIONS"])

        quicklook_image(
            log=self.log, CCDObject=normalisedFlats[0], show=False)

        # STACK THE NORMALISED FLAT FRAMES
        combined_normalised_flat = self.clip_and_stack(
            frames=normalisedFlats, recipe="soxs_mflat")

        quicklook_image(
            log=self.log, CCDObject=combined_normalised_flat, show=False)

        # DIVIDE THROUGH BY FIRST-PASS MASTER FRAME TO REMOVE CROSS-PLANE
        # ILLUMINATION VARIATIONS
        exposureFrames = []
        exposureFrames[:] = [
            n.divide(combined_normalised_flat) for n in calibratedFlats]

        quicklook_image(
            log=self.log, CCDObject=exposureFrames[0], show=False, stdWindow=0.001)

        # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
        # FLUX TO THAT LEVEL (AGAIN!)
        normalisedFlats = self.normalise_flats(
            calibratedFlats, orderTablePath=self.supplementaryInput[arm]["ORDER_LOCATIONS"], exposureFrames=exposureFrames)

        quicklook_image(
            log=self.log, CCDObject=normalisedFlats[0], show=False)

        # STACK THE RE-NORMALISED FLAT FRAMES
        combined_normalised_flat = self.clip_and_stack(
            frames=normalisedFlats, recipe="soxs_mflat")

        quicklook_image(
            log=self.log, CCDObject=combined_normalised_flat, show=True)

        # DETECT THE ORDER EDGES AND UPDATE THE ORDER LOCATIONS TABLE
        edges = detect_order_edges(
            log=self.log,
            flatFrame=combined_normalised_flat,
            orderCentreTable=self.supplementaryInput[arm]["ORDER_LOCATIONS"],
            settings=self.settings,
        )
        orderTablePath = edges.get()

        quicklook_image(
            log=self.log, CCDObject=combined_normalised_flat, show=False)

        background = subtract_background(
            log=self.log,
            frame=combined_normalised_flat,
            orderTable=orderTablePath,
            settings=self.settings
        )
        backgroundFrame, mflat = background.subtract()

        quicklook_image(
            log=self.log, CCDObject=mflat, show=True)

        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        # filePath = f"{outDir}/first_iteration_{arm}_master_flat.fits"

        productPath = self._write(
            mflat, outDir, overwrite=True)

        if 1 == 1:
            filename = filenamer(
                log=self.log,
                frame=mflat,
                settings=self.settings
            )
            filename = filename.replace(".fits", "_background.fits")
            filepath = self._write(
                backgroundFrame, outDir, filename=filename, overwrite=True)
            print(f"\nbackground frame saved to {filepath}\n")

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
        filterDict = {kw("DPR_TYPE").lower(): "BIAS"}
        biasCollection = self.inputFrames.filter(**filterDict)
        # LIST OF CCDDATA OBJECTS
        biases = [c for c in biasCollection.ccds(ccd_kwargs={
            "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        if len(biasCollection.files) == 0:
            bias = None
            biasCollection = None
        else:
            bias = biases[0]

        # FIND THE DARK FRAMES
        filterDict = {kw("DPR_TYPE").lower(): "DARK"}
        darkCollection = self.inputFrames.filter(**filterDict)

        if len(darkCollection.files) == 0:
            filterDict = {kw("DPR_TYPE").lower(): "LAMP,FLAT",
                          kw("DPR_TECH").lower(): "IMAGE"}
            darkCollection = self.inputFrames.filter(**filterDict)
            if len(darkCollection.files) == 0:
                darkCollection = None

        # FIND THE FLAT FRAMES
        filterDict = {kw("DPR_TYPE").lower(): "LAMP,(D|Q)?FLAT",
                      kw("DPR_TECH").lower(): "ECHELLE,SLIT"}
        flatCollection = self.inputFrames.filter(
            regex_match=True, **filterDict)

        if len(flatCollection.files) == 0:
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
        # AND/OR DARK
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
        orderTablePath,
        exposureFrames=None
    ):
        """*determine the median exposure for each flat frame and normalise the flux to that level*

        **Key Arguments:**
            - ``inputFlats`` -- the input flat field frames
            - ``orderTablePath`` -- path to the order table
            - ``exposureFrames`` -- frames where flux represents the frames exposure level. Default None.

        **Return:**
            - ``normalisedFrames`` -- the normalised flat-field frames (CCDData array)
        """
        self.log.debug('starting the ``normalise_flats`` method')

        # DO WE HAVE SEPARATE EXPOSURE FRAMES?
        if not exposureFrames:
            exposureFrames = inputFlats

        window = int(self.settings[
            "soxs-mflat"]["centre-order-window"] / 2)

        # UNPACK THE ORDER TABLE
        orderTableMeta, orderTablePixels = unpack_order_table(
            log=self.log, orderTablePath=orderTablePath)

        mask = np.ones_like(exposureFrames[0].data)

        xcoords = orderTablePixels["xcoord_centre"].values
        ycoords = orderTablePixels["ycoord"].values
        xcoords = xcoords.astype(int)

        # UPDATE THE MASK
        for x, y in zip(xcoords, ycoords):
            mask[y][x - window:x + window] = 0

        # COMBINE MASK WITH THE BAD PIXEL MASK
        mask = (mask == 1) | (inputFlats[0].mask == 1)
        normalisedFrames = []

        # PLOT ONE OF THE MASKED FRAMES TO CHECK
        for frame in [exposureFrames[0]]:
            maskedFrame = ma.array(frame.data, mask=mask)
            quicklook_image(log=self.log, CCDObject=maskedFrame,
                            show=False, ext=None)

        for frame, exp in zip(inputFlats, exposureFrames):
            maskedFrame = ma.array(exp.data, mask=mask)

            # SIGMA-CLIP THE DATA BEFORE CALCULATING MEAN
            maskedFrame = sigma_clip(
                maskedFrame, sigma_lower=2.5, sigma_upper=2.5, maxiters=3, cenfunc='median', stdfunc=mad_std)

            mean = np.ma.mean(maskedFrame)
            normalisedFrame = frame.divide(mean)
            normalisedFrame.header = frame.header
            normalisedFrames.append(normalisedFrame)

        # PLOT ONE OF THE NORMALISED FRAMES TO CHECK
        quicklook_image(
            log=self.log, CCDObject=normalisedFrames[0], show=False)

        self.log.debug('completed the ``normalise_flats`` method')
        return normalisedFrames

    # use the tab-trigger below for new method
    # xt-class-method


def nearest_neighbour(singleValue, listOfValues):
    arrayOfValues = np.asarray(listOfValues)
    dist = np.square(arrayOfValues - singleValue)
    minDist = np.amin(dist)
    minIndex = np.where(dist == minDist)[0][0]
    matchValue = listOfValues[minIndex]
    return matchValue, minIndex
