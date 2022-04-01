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
from soxspipe.commonutils.toolkit import generic_quality_checks, spectroscopic_image_quality_checks
from datetime import datetime
from astropy.stats import sigma_clip, mad_std, sigma_clipped_stats
from soxspipe.commonutils.filenamer import filenamer
from os.path import expanduser
from soxspipe.commonutils import subtract_background
import pandas as pd
from soxspipe.commonutils import detect_order_edges
from soxspipe.commonutils.toolkit import quicklook_image
import numpy.ma as ma
from soxspipe.commonutils.toolkit import unpack_order_table
import matplotlib.pyplot as plt
from soxspipe.commonutils import keyword_lookup
import ccdproc
from astropy.nddata import CCDData
import numpy as np
from ._base_recipe_ import _base_recipe_
from soxspipe.commonutils import set_of_files
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class soxs_mflat(_base_recipe_):
    """
    *The soxs_mflat recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
        - ``verbose`` -- verbose. True or False. Default *False*  

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
            inputFrames=[],
            verbose=False

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_mflat, self).__init__(
            log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mflat' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-mflat"
        self.recipeSettings = settings[self.recipeName]

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

        # SORT IMAGE COLLECTION
        self.inputFrames.sort(['mjd-obs'])
        if self.verbose:
            print("# RAW INPUT FRAMES - SUMMARY")
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
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

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

        # LOOK FOR ORDER TABLE
        arm = self.arm
        if f"ORDER_TAB_{arm}" not in imageCat:
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
        kw = self.kw

        # CALIBRATE THE FRAMES BY SUBTRACTING BIAS AND/OR DARK
        calibratedFlats = self.calibrate_frame_set()

        quicklook_image(
            log=self.log, CCDObject=calibratedFlats[0], show=False)

        # FIND THE ORDER TABLE
        filterDict = {kw("PRO_CATG").lower(): f"ORDER_TAB_{arm}"}
        orderTablePath = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
        # FLUX TO THAT LEVEL
        normalisedFlats = self.normalise_flats(
            calibratedFlats, orderTablePath=orderTablePath)

        quicklook_image(
            log=self.log, CCDObject=normalisedFlats[0], show=False)

        # STACK THE NORMALISED FLAT FRAMES
        combined_normalised_flat = self.clip_and_stack(
            frames=normalisedFlats, recipe="soxs_mflat")

        quicklook_image(
            log=self.log, CCDObject=combined_normalised_flat, show=False)

        # DIVIDE THROUGH BY FIRST-PASS MASTER FRAME TO REMOVE CROSS-PLANE
        # ILLUMINATION VARIATIONS
        print("\n# DIVIDING EACH ORIGINAL FLAT FRAME BY FIRST PASS MASTER FLAT")
        exposureFrames = []
        exposureFrames[:] = [
            n.divide(combined_normalised_flat) for n in calibratedFlats]

        quicklook_image(
            log=self.log, CCDObject=exposureFrames[0], show=False, stdWindow=0.001)

        # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
        # FLUX TO THAT LEVEL (AGAIN!)
        normalisedFlats = self.normalise_flats(
            calibratedFlats, orderTablePath=orderTablePath, exposureFrames=exposureFrames)

        quicklook_image(
            log=self.log, CCDObject=normalisedFlats[0], show=False)

        # STACK THE RE-NORMALISED FLAT FRAMES
        combined_normalised_flat = self.clip_and_stack(
            frames=normalisedFlats, recipe="soxs_mflat")

        quicklook_image(
            log=self.log, CCDObject=combined_normalised_flat, show=False)

        # DETECT THE ORDER EDGES AND UPDATE THE ORDER LOCATIONS TABLE
        edges = detect_order_edges(
            log=self.log,
            flatFrame=combined_normalised_flat,
            orderCentreTable=orderTablePath,
            settings=self.settings,
            qcTable=self.qc,
            productsTable=self.products
        )
        self.products, self.qc = edges.get()
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (self.products['product_label'] == "ORDER_LOC")
        orderTablePath = self.products.loc[mask]["file_path"].values[0]

        self.dateObs = combined_normalised_flat.header[self.kw("DATE_OBS")]

        quicklook_image(
            log=self.log, CCDObject=combined_normalised_flat, show=False)

        mflat = self.mask_low_sens_and_inter_order_to_unity(
            frame=combined_normalised_flat, orderTablePath=orderTablePath)

        # background = subtract_background(
        #     log=self.log,
        #     frame=combined_normalised_flat,
        #     orderTable=orderTablePath,
        #     settings=self.settings
        # )
        # backgroundFrame, mflat = background.subtract()

        quicklook_image(
            log=self.log, CCDObject=mflat, show=False)

        home = expanduser("~")
        outDir = self.settings["intermediate-data-root"].replace("~", home)
        # filePath = f"{outDir}/first_iteration_{arm}_master_flat.fits"

        # WRITE MFLAT TO FILE

        productPath = self._write(
            mflat, outDir, overwrite=True)
        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        basename = os.path.basename(productPath)
        self.products = self.products.append({
            "soxspipe_recipe": self.recipeName,
            "product_label": "MFLAT",
            "file_name": basename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} master spectroscopic flat frame",
            "file_path": productPath
        }, ignore_index=True)

        if 1 == 0:
            filename = filenamer(
                log=self.log,
                frame=mflat,
                settings=self.settings
            )
            filename = filename.replace(".fits", "_background.fits")
            filepath = self._write(
                backgroundFrame, outDir, filename=filename, overwrite=True)
            filepath = os.path.abspath(filepath)
            self.products = self.products.append({
                "soxspipe_recipe": self.recipeName,
                "product_label": "",
                "file_name": filename,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"modelled scatter background light image (removed from master flat)",
                "file_path": backgroundFrame
            }, ignore_index=True)

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log, frame=mflat, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)
        self.qc = spectroscopic_image_quality_checks(
            log=self.log, frame=mflat, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc, orderTablePath=orderTablePath)

        self.clean_up()
        self.report_output()

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
                print("\n# SUBTRACTING MASTER BIAS FROM FRAMES")
                calibratedFlats.append(self.detrend(
                    inputFrame=flat, master_bias=bias, dark=None))

        # IF DARKS EXIST - FIND CLOSEST IN TIME TO FLAT-FRAME. SUBTRACT BIAS
        # AND/OR DARK
        if darkCollection:
            darkMjds = [h[kw("MJDOBS").lower()]
                        for h in darkCollection.headers()]
            darks = [c for c in darkCollection.ccds(ccd_kwargs={
                "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]
            print("\n# SUBTRACTING MASTER DARK/OFF-LAMP FROM FRAMES")
            for flat in flats:
                mjd = flat.header[kw("MJDOBS").lower()]
                matchValue, matchIndex = nearest_neighbour(
                    flat.header[kw("MJDOBS").lower()], darkMjds)
                dark = darks[matchIndex]
                calibratedFlats.append(self.detrend(
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
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
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

        print("\n# NORMALISING FLAT FRAMES TO THEIR MEAN EXPOSURE LEVEL")
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

    def mask_low_sens_and_inter_order_to_unity(
            self,
            frame,
            orderTablePath):
        """*set inter-order pixels to unity and and low-sensitivity pixels to bad-pixel mask*

        **Key Arguments:**
            - ``frame`` -- the frame to work on
            - ``orderTablePath`` -- path to the order table

        **Return:**
            - ``frame`` -- with inter-order pixels set to 1 and BPM updated with low-sensitivity pixels
        """
        self.log.debug(
            'starting the ``mask_low_sens_and_inter_order_to_unity`` method')

        print("\n# CLIPPING LOW-SENSITIVITY PIXELS AND SETTING INTER-ORDER AREA TO UNITY")

        # UNPACK THE ORDER TABLE
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=orderTablePath)

        print(orderTablePath)

        # BAD PIXEL COUNT AT START
        originalBPM = np.copy(frame.mask)

        interOrderMask = np.ones_like(frame.data)
        xcoords_up = orderTablePixels["xcoord_edgeup"].values
        xcoords_low = orderTablePixels["xcoord_edgelow"].values
        ycoords = orderTablePixels["ycoord"].values
        xcoords_up = xcoords_up.astype(int)
        xcoords_low = xcoords_low.astype(int)

        # UPDATE THE MASK TO INCLUDE INTER-ORDER BACKGROUND
        for u, l, y in zip(xcoords_up, xcoords_low, ycoords):
            interOrderMask[y][l:u] = 0

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        interOrderMask = ma.make_mask(interOrderMask)
        frame.mask = (interOrderMask == 1) | (frame.mask == 1)

        # PLOT MASKED FRAMES TO CHECK
        quicklook_image(log=self.log, CCDObject=frame,
                        show=False, ext=None)

        beforeMask = np.copy(frame.mask)

        # SIGMA-CLIP THE LOW-SENSITIVITY PIXELS
        frameClipped = sigma_clip(
            frame, sigma_lower=self.settings["soxs-mflat"]["low-sensitivity-clipping-simga"], sigma_upper=2000, maxiters=5, cenfunc='median', stdfunc=mad_std)

        # PLOT MASKED FRAMES TO CHECK
        quicklook_image(log=self.log, CCDObject=frame,
                        show=False, ext=None)

        lowSensitivityPixelMask = (frameClipped.mask == 1) & (beforeMask != 1)
        lowSensPixelCount = lowSensitivityPixelMask.sum()

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "low sensitivity pixel count",
            "qc_value": lowSensPixelCount,
            "qc_unit": "pixels",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow
        }, ignore_index=True)
        print(f"        {lowSensPixelCount} low-sensitivity pixels added to bad-pixel mask")

        frame.mask = (lowSensitivityPixelMask == 1) | (originalBPM == 1)

        # SET INTRA-ORDER TO 1
        frame.data[interOrderMask] = 1

        # PLOT MASKED FRAMES TO CHECK
        quicklook_image(log=self.log, CCDObject=frame,
                        show=False, ext=None)

        self.log.debug(
            'completed the ``mask_low_sens_and_inter_order_to_unity`` method')
        return frame

    # use the tab-trigger below for new method
    # xt-class-method


def nearest_neighbour(singleValue, listOfValues):
    arrayOfValues = np.asarray(listOfValues)
    dist = np.square(arrayOfValues - singleValue)
    minDist = np.amin(dist)
    minIndex = np.where(dist == minDist)[0][0]
    matchValue = listOfValues[minIndex]
    return matchValue, minIndex
