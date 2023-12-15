#!/usr/bin/env python
# encoding: utf-8
"""
*generate a single normalised master flat-field frame*

:Author:
    David Young & Marco Landoni

:Date Created:
    September 16, 2020
"""
from soxspipe.commonutils.toolkit import generic_quality_checks, spectroscopic_image_quality_checks
from datetime import datetime

from soxspipe.commonutils.filenamer import filenamer
from os.path import expanduser
from soxspipe.commonutils import subtract_background
from soxspipe.commonutils import detect_order_edges
from soxspipe.commonutils.toolkit import quicklook_image
from soxspipe.commonutils.toolkit import unpack_order_table
from soxspipe.commonutils import keyword_lookup
from ._base_recipe_ import _base_recipe_
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
        - ``overwrite`` -- overwrite the prodcut file if it already exists. Default *False*

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
            verbose=False,
            overwrite=False

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_mflat, self).__init__(
            log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-mflat")
        self.log = log
        log.debug("instansiating a new 'soxs_mflat' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeSettings = settings[self.recipeName]

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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY soxs_mflat - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # SORT IMAGE COLLECTION BY MJD
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
        """*verify the input frames match those required by the soxs_mflat recipe*

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception will be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        error = False

        import warnings
        from astropy.utils.exceptions import AstropyWarning
        warnings.simplefilter('ignore', AstropyWarning)

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if self.arm == "NIR":
            # WANT ON AND OFF PINHOLE FRAMES
            # MIXED INPUT IMAGE TYPES ARE BAD
            if not error:
                if len(imageTypes) > 1:
                    imageTypes = " and ".join(imageTypes)
                    error = "Input frames are a mix of %(imageTypes)s" % locals()

            if not error:
                if imageTypes[0] != "LAMP,FLAT":
                    error = "Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR" % locals()

            if not error:
                for i in imageTech:
                    if i not in ['ECHELLE,SLIT', 'IMAGE']:
                        error = "Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR. You have provided {imageTech}" % locals()

            if not error:
                for i in ['ECHELLE,SLIT', 'IMAGE']:
                    if i not in imageTech:
                        error = "Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR. You have provided {imageTech}" % locals()

        else:
            if not error:
                for i in imageTypes:
                    if i not in ["LAMP,FLAT", "LAMP,QFLAT", "LAMP,DFLAT"]:
                        error = "Input frames for soxspipe mflat need to be flat-lamp frames,a master-bias frame, an order-locations tables and possibly a master dark for UVB/VIS" % locals()

            if not error:
                for i in [f"MASTER_BIAS_{self.arm}", f"ORDER_TAB_{self.arm}"]:
                    if i not in imageCat:
                        error = "Input frames for soxspipe mflat need to be flat-lamp frames,a master-bias frame, an order-locations tables and possibly a master dark for UVB/VIS" % locals()

            if not error:
                found = False
                for i in ["LAMP,FLAT", "LAMP,QFLAT", "LAMP,DFLAT"]:
                    if i in imageTypes:
                        found = True
                if not found:
                    error = "Input frames for soxspipe mflat need to be flat-lamp frames,a master-bias frame, an order-locations tables and possibly a master dark for UVB/VIS" % locals()

        # UV-VIS NEEDS BOTH D AND Q-LAMPS
        if not error:
            if "LAMP,QFLAT" in imageTypes or "LAMP,DFLAT" in imageTypes:
                if "LAMP,QFLAT" in imageTypes and "LAMP,DFLAT" in imageTypes:
                    pass
                else:
                    for i in imageTypes:
                        if "LAMP" in i:
                            error = f'Only "{i}" image types found. Please include both D2 and QTH lamp flats'

        # LOOK FOR ORDER TABLE
        if not error:
            arm = self.arm
            if f"ORDER_TAB_{arm}" not in imageCat:
                error = f"Need an order centre for {arm} - none found with the input files"

        # CHECK EXPTIME
        if not error:
            lamps = ["LAMP,FLAT", "LAMP,DFLAT", "LAMP,QFLAT"]
            for l in lamps:
                filterDict = {kw("DPR_TYPE"): l,
                              kw("DPR_TECH"): "ECHELLE,SLIT"}
                flatCollection = self.inputFrames.filter(**filterDict)
                if len(flatCollection.files):
                    exptime = flatCollection.values(
                        keyword=kw("EXPTIME"), unique=True)
                    if len(exptime) > 1:
                        error = f"Input {l} frames for soxspipe mflat need to have a unique exptime"

        if error:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            self.log.print()
            raise TypeError(error)

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

        import pandas as pd

        productPath = None
        arm = self.arm
        kw = self.kw

        # CALIBRATE THE FRAMES BY SUBTRACTING BIAS AND/OR DARK
        calibratedFlats, dcalibratedFlats, qcalibratedFlats = self.calibrate_frame_set()

        allCalibratedFlats = calibratedFlats + dcalibratedFlats + qcalibratedFlats
        quicklook_image(log=self.log, CCDObject=allCalibratedFlats[0], show=False, ext="Data", surfacePlot=True, title="Single bias and/or dark subtracted flat frame")

        calibratedFlatSet = [calibratedFlats, dcalibratedFlats, qcalibratedFlats]
        flatKeywords = ["LAMP,ORDERDEF", 'LAMP,DORDERDEF', 'LAMP,QORDERDEF']
        lampTag = ["", "_DLAMP", "_QLAMP"]
        normalisedFlatSet = []
        self.combinedNormalisedFlatSet = []
        self.masterFlatSet = []
        self.orderTableSet = []
        self.detectionCountSet = []
        medianOrderFluxDFExists = False

        productTable = self.products
        qcTable = self.qc

        for cf, fk, tag in zip(calibratedFlatSet, flatKeywords, lampTag):

            if len(cf) == 0:
                self.orderTableSet.append(None)
                normalisedFlatSet.append(None)
                self.combinedNormalisedFlatSet.append(None)
                self.masterFlatSet.append(None)
                # productTable = self.products.copy()
                # qcTable = self.qc.copy()
                continue

            if tag:
                filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}", kw("OBJECT"): fk}
            else:
                filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}"}
            orderTablePaths = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)
            if len(orderTablePaths) == 1:
                orderTablePath = orderTablePaths[0]
            else:
                self.orderTableSet.append(None)

            # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
            # FLUX TO THAT LEVEL
            normalisedFlats = self.normalise_flats(
                cf, orderTablePath=orderTablePath)

            quicklook_image(log=self.log, CCDObject=normalisedFlats[0], show=False, ext=True, surfacePlot=True, title="Single normalised flat frame")
            # STACK THE NORMALISED FLAT FRAMES
            combined_normalised_flat = self.clip_and_stack(
                frames=normalisedFlats, recipe="soxs_mflat", ignore_input_masks=False, post_stack_clipping=True)
            quicklook_image(log=self.log, CCDObject=combined_normalised_flat, show=False, ext=None, surfacePlot=True, title="Combined normalised flat frames")

            # DIVIDE THROUGH BY FIRST-PASS MASTER FRAME TO REMOVE CROSS-PLANE
            # ILLUMINATION VARIATIONS
            self.log.print("\n# DIVIDING EACH ORIGINAL FLAT FRAME BY FIRST PASS MASTER FLAT")
            exposureFrames = []
            exposureFrames[:] = [
                n.divide(combined_normalised_flat) for n in cf]
            quicklook_image(log=self.log, CCDObject=exposureFrames[0], show=False, ext=None, surfacePlot=True, title="Single exposure map of flat frame", inst=self.inst)

            # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
            # FLUX TO THAT LEVEL (AGAIN!)
            normalisedFlats = self.normalise_flats(
                cf, orderTablePath=orderTablePath, exposureFrames=exposureFrames)

            quicklook_image(log=self.log, CCDObject=normalisedFlats[0], show=False, ext=None, surfacePlot=True, title="Single re-normalised flat frame")

            # STACK THE RE-NORMALISED FLAT FRAMES
            combined_normalised_flat = self.clip_and_stack(
                frames=normalisedFlats, recipe="soxs_mflat", ignore_input_masks=False, post_stack_clipping=True)

            quicklook_image(log=self.log, CCDObject=combined_normalised_flat, show=False, ext=None, surfacePlot=True, title="Recombined normalised flat frames")

            self.combinedNormalisedFlatSet.append(combined_normalised_flat.copy())

            # DETECT THE ORDER EDGES AND UPDATE THE ORDER LOCATIONS TABLE
            edges = detect_order_edges(
                log=self.log,
                flatFrame=combined_normalised_flat,
                orderCentreTable=orderTablePath,
                settings=self.settings,
                qcTable=qcTable,
                productsTable=productTable,
                tag=tag,
                sofName=self.sofName,
                binx=self.binx,
                biny=self.biny
            )
            productTable, qcTable, orderDetectionCounts = edges.get()

            if tag:
                # NEED TO TRY AND RENAME BOTH ORDER AND COUNT COLUMNS FOR PANDAS 1.X and 2.X
                orderDetectionCounts.rename(columns={"order": tag}, inplace=True)
                orderDetectionCounts.rename(columns={"count": tag}, inplace=True)
                orderDetectionCounts.index.names = ['order']

            self.detectionCountSet.append(orderDetectionCounts)

            mask = (productTable['product_label'] == f"ORDER_LOC{tag}")
            orderTablePath = productTable.loc[mask]["file_path"].values[0]

            self.orderTableSet.append(orderTablePath)

            self.dateObs = combined_normalised_flat.header[self.kw("DATE_OBS")]

            # UNPACK THE ORDER TABLE
            orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
                log=self.log, orderTablePath=orderTablePath)

            if tag in ("_DLAMP", "_QLAMP"):
                zero = True
            else:
                zero = False

            mflat, medianOrderFluxDF = self.mask_low_sens_and_inter_order_to_unity(
                frame=combined_normalised_flat, orderTablePath=orderTablePath, returnMedianOrderFlux=True, zero=zero)

            self.masterFlatSet.append(mflat)

            if tag:
                medianOrderFluxDF.rename(columns={"medianFlux": tag}, inplace=True)

            if tag and not medianOrderFluxDFExists:
                medianOrderFluxDFExists = True
                medianOrderFluxDFFirst = medianOrderFluxDF.copy()
            elif tag:
                medianOrderFluxDF = pd.merge(medianOrderFluxDFFirst, medianOrderFluxDF)

            # background = subtract_background(
            #     log=self.log,
            #     frame=combined_normalised_flat,
            #     orderTable=orderTablePath,
            #     settings=self.settings
            # )
            # backgroundFrame, mflat = background.subtract()

        # UV-STITCHING
        if len(self.detectionCountSet) > 1:
            mflat = self.stitch_uv_mflats(medianOrderFluxDF)
        else:
            self.products = productTable
            self.qc = qcTable

        quicklook_image(log=self.log, CCDObject=combined_normalised_flat, show=False, ext=None, surfacePlot=True, title="Final master flat frame")

        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home)
        # filePath = f"{outDir}/first_iteration_{arm}_master_flat.fits"

        self.update_fits_keywords(
            frame=mflat
        )

        # WRITE MFLAT TO FILE
        productPath = self._write(
            mflat, outDir, overwrite=True)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        basename = os.path.basename(productPath)
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": f"MFLAT",
            "file_name": basename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} master spectroscopic flat frame",
            "file_path": productPath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

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
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": "",
                "file_name": filename,
                "file_type": "FITS",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"modelled scatter background light image (removed from master flat)",
                "file_path": backgroundFrame,
                "label": "PROD"
            }).to_frame().T], ignore_index=True)

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
        filterDict = {kw("PRO_CATG"): f"MASTER_BIAS_{self.arm.upper()}"}
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
        filterDict = {kw("PRO_CATG"): f"MASTER_DARK_{self.arm.upper()}"}
        darkCollection = self.inputFrames.filter(**filterDict)

        if len(darkCollection.files) == 0:
            filterDict = {kw("DPR_TYPE"): "LAMP,FLAT",
                          kw("DPR_TECH"): "IMAGE"}
            darkCollection = self.inputFrames.filter(**filterDict)
            if len(darkCollection.files) == 0:
                darkCollection = None

        # FIND THE FLAT FRAMES
        filterDict = {kw("DPR_TYPE"): "LAMP,FLAT",
                      kw("DPR_TECH"): "ECHELLE,SLIT"}
        flatCollection = self.inputFrames.filter(**filterDict)

        filterDict = {kw("DPR_TYPE"): "LAMP,DFLAT",
                      kw("DPR_TECH"): "ECHELLE,SLIT"}
        dflatCollection = self.inputFrames.filter(**filterDict)
        filterDict = {kw("DPR_TYPE"): "LAMP,QFLAT",
                      kw("DPR_TECH"): "ECHELLE,SLIT"}
        qflatCollection = self.inputFrames.filter(**filterDict)

        if len(flatCollection.files) == 0 and len(dflatCollection.files) == 0 and len(qflatCollection.files) == 0:
            raise FileNotFoundError(
                "The mflat recipe needs flat-frames as input, none found")

        # LIST OF CCDDATA OBJECTS
        flats = [c for c in flatCollection.ccds(ccd_kwargs={
                                                "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]
        dflats = [c for c in dflatCollection.ccds(ccd_kwargs={
            "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]
        qflats = [c for c in qflatCollection.ccds(ccd_kwargs={
            "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        # IF NO DARK FRAMES EXIST - JUST A MASTER BIAS. SUBTRACT BIAS.
        calibratedFlats = []
        dcalibratedFlats = []
        qcalibratedFlats = []
        if not darkCollection and bias:
            self.log.print("\n# SUBTRACTING MASTER BIAS FROM FRAMES")
            for flat in flats:
                calibratedFlats.append(self.detrend(
                    inputFrame=flat, master_bias=bias, dark=None))
            for flat in dflats:
                dcalibratedFlats.append(self.detrend(
                    inputFrame=flat, master_bias=bias, dark=None))
            for flat in qflats:
                qcalibratedFlats.append(self.detrend(
                    inputFrame=flat, master_bias=bias, dark=None))

        # IF DARKS EXIST - FIND CLOSEST IN TIME TO FLAT-FRAME. SUBTRACT BIAS
        # AND/OR DARK
        if darkCollection:
            darkMjds = [h[kw("MJDOBS")]
                        for h in darkCollection.headers()]
            darks = [c for c in darkCollection.ccds(ccd_kwargs={
                "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]
            self.log.print("\n# SUBTRACTING MASTER DARK/OFF-LAMP FROM FRAMES")
            for flat in flats:
                mjd = flat.header[kw("MJDOBS")]
                matchValue, matchIndex = nearest_neighbour(
                    flat.header[kw("MJDOBS")], darkMjds)
                dark = darks[matchIndex]
                calibratedFlats.append(self.detrend(
                    inputFrame=flat, master_bias=bias, dark=dark))

        for frame in calibratedFlats:
            quicklook_image(log=self.log, CCDObject=frame, show=False)

        if 1 == 0:
            from os.path import expanduser
            home = expanduser("~")
            outDir = self.settings["workspace-root-dir"].replace("~", home)
            index = 1
            for frame in calibratedFlats:
                filePath = f"{outDir}/{index:02}_flat_{arm}_calibrated.fits"
                index += 1
                self._write(frame, filePath, overwrite=True)

        self.log.debug('completed the ``calibrate_frame_set`` method')
        return calibratedFlats, dcalibratedFlats, qcalibratedFlats

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

        import numpy.ma as ma
        import numpy as np
        from astropy.stats import sigma_clip
        kw = self.kw

        # DO WE HAVE SEPARATE EXPOSURE FRAMES?
        if not exposureFrames:
            exposureFrames = inputFlats

        try:
            self.binx = exposureFrames[0].header[kw("WIN_BINX")]
            self.biny = exposureFrames[0].header[kw("WIN_BINY")]
        except:
            if self.arm.lower() == "nir":
                self.binx = 1
                self.biny = 1

        window = int(self.settings[
            "soxs-mflat"]["centre-order-window"] / 2)
        if self.axisA == "x" and self.binx > 1:
            window = int(window / self.binx)
            self.recipeSettings["slice-length-for-edge-detection"] = int(self.recipeSettings["slice-length-for-edge-detection"] / self.binx)
        elif self.axisA == "y" and self.biny > 1:
            window = int(window / self.biny)
            self.recipeSettings["slice-length-for-edge-detection"] = int(self.recipeSettings["slice-length-for-edge-detection"] / self.biny)

        # UNPACK THE ORDER TABLE
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=orderTablePath, binx=self.binx, biny=self.biny)

        mask = np.ones_like(exposureFrames[0].data)

        axisAcoords = orderTablePixels[f"{self.axisA}coord_centre"].values
        axisBcoords = orderTablePixels[f"{self.axisB}coord"].values
        axisAcoords = axisAcoords.astype(int)

        # UPDATE THE MASK
        if self.axisA == "x":
            for x, y in zip(axisAcoords, axisBcoords):
                mask[y][x - window:x + window] = 0
        else:
            for y, x in zip(axisAcoords, axisBcoords):
                mask[y][x - window:x + window] = 0

        # COMBINE MASK WITH THE BAD PIXEL MASK
        mask = (mask == 1) | (inputFlats[0].mask == 1)
        normalisedFrames = []

        # PLOT ONE OF THE MASKED FRAMES TO CHECK
        for frame in [exposureFrames[0]]:
            maskedFrame = ma.array(frame.data, mask=mask)
            quicklook_image(log=self.log, CCDObject=maskedFrame, show=False, ext=None, surfacePlot=True, title="Single masked flat frame")

        self.log.print("\n# NORMALISING FLAT FRAMES TO THEIR MEAN EXPOSURE LEVEL")
        for frame, exp in zip(inputFlats, exposureFrames):
            maskedFrame = ma.array(exp.data, mask=mask)

            # SIGMA-CLIP THE DATA BEFORE CALCULATING MEAN
            maskedFrame = sigma_clip(
                maskedFrame, sigma_lower=2.5, sigma_upper=2.5, maxiters=3, cenfunc='median', stdfunc='mad_std')

            mean = np.ma.mean(maskedFrame)
            normalisedFrame = frame.divide(mean)
            normalisedFrame.header = frame.header
            normalisedFrames.append(normalisedFrame)

        # PLOT ONE OF THE NORMALISED FRAMES TO CHECK
        quicklook_image(
            log=self.log, CCDObject=normalisedFrames[0], show=False, ext=None, surfacePlot=True, title="Single normalised flat frame")

        self.log.debug('completed the ``normalise_flats`` method')
        return normalisedFrames

    def mask_low_sens_and_inter_order_to_unity(
            self,
            frame,
            orderTablePath,
            returnMedianOrderFlux=False,
            zero=False):
        """*set inter-order pixels to unity and and low-sensitivity pixels to bad-pixel mask*

        **Key Arguments:**
            - ``frame`` -- the frame to work on
            - ``orderTablePath`` -- path to the order table
            - ``returnMedianOrderFlux`` -- return a table of the median order fluxes. Default *False*.
            - ``zero`` -- set inter-order pixels to zero instead of one

        **Return:**
            - ``frame`` -- with inter-order pixels set to 1 and BPM updated with low-sensitivity pixels
            - ``medianOrderFluxDF`` -- data-frame of the median order fluxes (if ``returnMedianOrderFlux`` is True)
        """
        self.log.debug(
            'starting the ``mask_low_sens_and_inter_order_to_unity`` method')

        import pandas as pd
        import numpy.ma as ma
        import numpy as np
        from astropy.stats import sigma_clip

        self.log.print("\n# CLIPPING LOW-SENSITIVITY PIXELS AND SETTING INTER-ORDER AREA TO UNITY")

        # UNPACK THE ORDER TABLE
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=orderTablePath, binx=self.binx, biny=self.biny, prebinned=True)

        # BAD PIXEL COUNT AT START
        originalBPM = np.copy(frame.mask)

        # SET MASK TO ZEROS OR ONES
        if zero:
            interOrderMask = np.zeros_like(frame.data)
        else:
            interOrderMask = np.ones_like(frame.data)

        orders = orderTablePixels["order"].values
        axisAcoords_up = orderTablePixels[f"{self.axisA}coord_edgeup"].values.round().astype(int)
        axisAcoords_low = orderTablePixels[f"{self.axisA}coord_edgelow"].values.round().astype(int)
        axisBcoords = orderTablePixels[f"{self.axisB}coord"].values
        uniqueOrders = orderTablePixels['order'].unique()

        # CALCULATE AND RETURN MEDIAN FLUXES FOR ORDERS
        if returnMedianOrderFlux:
            orderFluxes = {}
            bAxisMiddles = {}
            for o in orders:
                orderFluxes[o] = []
            for o in uniqueOrders:
                filteredDf = orderTablePixels.loc[(orderTablePixels["order"] == o)]
                bAxisMiddles[o] = int(filteredDf[f"{self.axisB}coord"].mean())
            medianFlux = []

        # UPDATE THE MASK TO ALLOW INTER-ORDER PIXELS
        for u, l, b, o in zip(axisAcoords_up, axisAcoords_low, axisBcoords, orders):
            if l < 0:
                l = 0
            if u < 0:
                u = 0
            if self.axisA == "x":
                interOrderMask[b, l:u] = 0
                if returnMedianOrderFlux and b > bAxisMiddles[o] - 3 and b < bAxisMiddles[o] + 3:
                    orderFluxes[o] = np.append(orderFluxes[o], frame.data[b, l:u])
            else:
                interOrderMask[l:u, b] = 0
                if returnMedianOrderFlux and b > bAxisMiddles[o] - 3 and b < bAxisMiddles[o] + 3:
                    orderFluxes[o] = np.append(orderFluxes[o], frame.data[b, l:u])

        # GET UNIQUE VALUES IN COLUMN

        if returnMedianOrderFlux:
            for o in uniqueOrders:
                medianFlux.append(np.median(orderFluxes[o]))

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        interOrderMask = ma.make_mask(interOrderMask)
        frame.mask = (interOrderMask == 1) | (frame.mask == 1)

        # PLOT MASKED FRAMES TO CHECK
        quicklook_image(log=self.log, CCDObject=frame, show=False, ext=None, surfacePlot=True, title="Masking inter-order pixels")

        beforeMask = np.copy(frame.mask)

        # SIGMA-CLIP THE LOW-SENSITIVITY PIXELS
        frameClipped = sigma_clip(
            frame, sigma_lower=self.settings["soxs-mflat"]["low-sensitivity-clipping-sigma"], sigma_upper=2000, maxiters=5, cenfunc='median', stdfunc='mad_std')

        lowSensitivityPixelMask = (frameClipped.mask == 1) & (beforeMask != 1)
        lowSensPixelCount = lowSensitivityPixelMask.sum()

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "NLOWSENS",
            "qc_value": lowSensPixelCount,
            "qc_comment": "Number of low-sensitivity pixels found in master flat",
            "qc_unit": "pixels",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow
        }).to_frame().T], ignore_index=True)
        self.log.print(f"        {lowSensPixelCount} low-sensitivity pixels added to bad-pixel mask")

        frame.mask = (lowSensitivityPixelMask == 1) | (originalBPM == 1)

        # SET INTRA-ORDER TO 1
        frame.data[interOrderMask] = 1

        # PLOT MASKED FRAMES TO CHECK
        quicklook_image(log=self.log, CCDObject=frame, show=False, ext=None, surfacePlot=True, title="Low Sensitivity pixels masked and inter-order pixel set to 1")

        self.log.debug(
            'completed the ``mask_low_sens_and_inter_order_to_unity`` method')

        if returnMedianOrderFlux:
            medianOrderFluxDF = {
                "order": uniqueOrders,
                "medianFlux": medianFlux
            }
            medianOrderFluxDF = pd.DataFrame(medianOrderFluxDF)
            return frame, medianOrderFluxDF

        return frame

    def stitch_uv_mflats(
            self,
            medianOrderFluxDF):
        """*return a master UV-VIS flat frame after determining the best order to slice and stitch the UV-VIS D-Lamp and QTH-Lamp flat frames*

        **Key Arguments:**
            - ``medianOrderFluxDF`` -- data frame containing median order fluxes for D and QTH frames

        **Return:**
            - ``stitchedFlat`` -- the stitch D and QTH-Lamp master flat frame

        **Usage:**

        ```python
        mflat = self.stitch_uv_mflats(medianOrderFluxDF)
        ```
        """
        self.log.debug('starting the ``stitch_uv_mflats`` method')

        import pandas as pd

        kw = self.kw

        # MERGE DETECTION DATAFRAMES - LISTING ALL ORDER EGDE LOCATIONS FOUND FOR EACH ORDER IN D AND QTH LAMP FRAMES
        detectionCount = pd.merge(self.detectionCountSet[0], self.detectionCountSet[1], left_index=True, right_index=True)

        detectionCount["best_frame"] = detectionCount.idxmax(axis=1)

        mask = (detectionCount['best_frame'] == "_QLAMP")
        orderFlip = detectionCount.loc[mask]
        # THIS IS THE ORDER THAT WE USE TO RESCALE ONE FLAT TO ANOTHER
        orderFlip = orderFlip.index.max() - 1
        filteredDf = medianOrderFluxDF.loc[(medianOrderFluxDF["order"] == orderFlip)]
        filteredDf["scale"] = filteredDf["_DLAMP"] / filteredDf["_QLAMP"]
        DQscale = filteredDf["scale"].values[0]
        # SCALE D FRAME TO QTH FRAME
        dmflatScaled = self.masterFlatSet[1].divide(DQscale)

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=dmflatScaled, show=False, ext=False, stdWindow=3, surfacePlot=True, title="D Flat scaled to Q-Flat")

        # UNPACK THE ORDER TABLE
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTableSet[1], extend=3000, binx=self.binx, biny=self.biny)

        # FIND THE LINE USED TO SLICE AND STITCH THE 2 FRAMES TOGETHER
        filteredDf = orderTablePixels.loc[(orderTablePixels["order"] == orderFlip)]
        axisAStitchCoords = filteredDf[f"{self.axisA}coord_edgeup"].values.astype(int) + 4
        axisBStitchCoords = filteredDf[f"{self.axisB}coord"].values
        stitchedFlat = self.masterFlatSet[2].copy()

        # STITCH FLAT FRAMES AND COMBINED NORMALISED FRAMES (NEEDED FOR BEST ORDER EDGE DETECTION) TOGETHER
        if self.axisA == "x":
            for x, y in zip(axisAStitchCoords, axisBStitchCoords):
                if y < stitchedFlat.data.shape[0] and x < stitchedFlat.data.shape[1]:
                    stitchedFlat.data[y, :x] = dmflatScaled.data[y, :x]
                    stitchedFlat.mask[y, :x] = dmflatScaled.mask[y, :x]
                    stitchedFlat.uncertainty.array[y, :x] = dmflatScaled.uncertainty.array[y, :x]
        else:
            for y, x in zip(axisAStitchCoords, axisBStitchCoords):
                stitchedFlat.data[y, x:] = dmflatScaled.data[y, x:]
                stitchedFlat.mask[y, x:] = dmflatScaled.mask[y, x:]
                stitchedFlat.uncertainty.array[y, x:] = dmflatScaled.uncertainty.array[y, x:]

        stitchedFlat.header[kw("DPR_TYPE")] = stitchedFlat.header[kw("DPR_TYPE")].replace(",D", ",").replace(",Q", ",")

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=stitchedFlat, show=False, ext=False, stdWindow=3, title=False, surfacePlot=True)

        # DETECT THE ORDER EDGES AND UPDATE THE ORDER LOCATIONS TABLE
        edges = detect_order_edges(
            log=self.log,
            flatFrame=stitchedFlat,
            orderCentreTable=self.orderTableSet[2],
            settings=self.settings,
            qcTable=self.qc,
            productsTable=self.products,
            tag="",
            sofName=self.sofName,
            binx=self.binx,
            biny=self.biny
        )
        self.products, self.qc, orderDetectionCounts = edges.get()
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (self.products['product_label'] == f"ORDER_LOC")
        orderTablePath = self.products.loc[mask]["file_path"].values[0]

        stitchedFlat = self.mask_low_sens_and_inter_order_to_unity(
            frame=stitchedFlat, orderTablePath=orderTablePath)

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=stitchedFlat, show=False, ext=False, stdWindow=3, title=False, surfacePlot=True)

        self.log.debug('completed the ``stitch_uv_mflats`` method')
        return stitchedFlat

    # use the tab-trigger below for new method
    # xt-class-method


def nearest_neighbour(singleValue, listOfValues):
    import numpy as np
    arrayOfValues = np.asarray(listOfValues)
    dist = np.square(arrayOfValues - singleValue)
    minDist = np.amin(dist)
    minIndex = np.where(dist == minDist)[0][0]
    matchValue = listOfValues[minIndex]
    return matchValue, minIndex
