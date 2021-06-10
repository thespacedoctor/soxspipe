#!/usr/bin/env python
# encoding: utf-8
"""
*The base recipe class which all other recipes inherit*

:Author:
    David Young & Marco Landoni

:Date Created:
    January 22, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
import numpy as np
from astropy.nddata import CCDData
from astropy import units as u
from astropy.stats import mad_std
import ccdproc
from astropy.nddata.nduncertainty import StdDevUncertainty
from ccdproc import Combiner
from soxspipe.commonutils import set_of_files
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils import detector_lookup
from datetime import datetime
from soxspipe.commonutils import filenamer
import shutil


class _base_recipe_(object):
    """
    The base recipe class which all other recipes inherit

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary

    **Usage**

    To use this base recipe to create a new `soxspipe` recipe, have a look at the code for one of the simpler recipes (e.g. `soxs_mbias`) - copy and modify the code.
    """

    def __init__(
            self,
            log,
            settings=False
    ):
        self.log = log
        log.debug("instansiating a new '__init__' object")
        self.settings = settings
        self.intermediateRootPath = self._absolute_path(
            settings["intermediate-data-root"])
        self.reducedRootPath = self._absolute_path(
            settings["reduced-data-root"])
        self.calibrationRootPath = self._absolute_path(
            settings["calibration-data-root"])

        # SET LATER WHEN VERIFYING FRAMES
        self.arm = None
        self.detectorParams = None

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get

        return None

    def _prepare_single_frame(
            self,
            frame,
            save=False):
        """*prepare a single raw frame by converting pixel data from ADU to electrons and adding mask and uncertainty extensions*

        **Key Arguments:**
            - ``frame`` -- the path to the frame to prepare, of a CCDData object

        **Return:**
            - ``frame`` -- the prepared frame with mask and uncertainty extensions (CCDData object)

        ```eval_rst
        .. todo::

            - write a command-line tool for this method
        ```
        """
        self.log.debug('starting the ``_prepare_single_frame`` method')

        kw = self.kw
        dp = self.detectorParams

        # STORE FILEPATH FOR LATER USE
        filepath = frame

        # CONVERT FILEPATH TO CCDDATA OBJECT
        if isinstance(frame, str):
            # CONVERT RELATIVE TO ABSOLUTE PATHS
            frame = self._absolute_path(frame)
            # OPEN THE RAW FRAME - MASK AND UNCERT TO BE POPULATED LATER
            frame = CCDData.read(frame, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                 hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # CHECK THE NUMBER OF EXTENSIONS IS ONLY 1 AND "SXSPRE" DOES NOT
        # EXIST. i.e. THIS IS A RAW UNTOUCHED FRAME
        if len(frame.to_hdu()) > 1 or "SXSPRE" in frame.header:
            return filepath

        # MANIPULATE XSH DATA
        frame = self.xsh2soxs(frame)
        frame = self._trim_frame(frame)

        # CORRECT FOR GAIN - CONVERT DATA FROM ADU TO ELECTRONS
        frame = ccdproc.gain_correct(frame, dp["gain"])

        # GENERATE UNCERTAINTY MAP AS EXTENSION
        if frame.header[kw("DPR_TYPE")] == "BIAS":
            # ERROR IS ONLY FROM READNOISE FOR BIAS FRAMES
            errorMap = np.ones_like(frame.data) * dp["ron"]
            # errorMap = StdDevUncertainty(errorMap)
            frame.uncertainty = errorMap
        else:
            # GENERATE UNCERTAINTY MAP AS EXTENSION
            frame = ccdproc.create_deviation(
                frame, readnoise=dp["ron"])

        # FIND THE APPROPRIATE BAD-PIXEL BITMAP AND APPEND AS 'FLAG' EXTENSION
        # NOTE FLAGS NOTE YET SUPPORTED BY CCDPROC THIS THIS WON'T GET SAVED OUT
        # AS AN EXTENSION
        arm = self.arm
        if kw('WIN_BINX') in frame.header:
            binx = int(frame.header[kw('WIN_BINX')])
            biny = int(frame.header[kw('WIN_BINY')])
        else:
            binx = 1
            biny = 1

        bitMapPath = self.calibrationRootPath + "/" + dp["bad-pixel map"][f"{binx}x{biny}"]

        if not os.path.exists(bitMapPath):
            message = "the path to the bitMapPath %s does not exist on this machine" % (
                bitMapPath,)
            self.log.critical(message)
            raise IOError(message)
        bitMap = CCDData.read(bitMapPath, hdu=0, unit=u.dimensionless_unscaled)

        # BIAS FRAMES HAVE NO 'FLUX', JUST READNOISE, SO ADD AN EMPTY BAD-PIXEL
        # MAP
        if frame.header[kw("DPR_TYPE")] == "BIAS":
            bitMap.data = np.zeros_like(bitMap.data)

        # print(bitMap.data.shape)
        # print(frame.data.shape)

        frame.flags = bitMap.data

        # FLATTEN BAD-PIXEL BITMAP TO BOOLEAN FALSE (GOOD) OR TRUE (BAD) AND
        # APPEND AS 'UNCERT' EXTENSION
        boolMask = bitMap.data.astype(bool).data
        try:
            # FAILS IN PYTHON 2.7 AS BOOLMASK IS A BUFFER - NEED TO CONVERT TO
            # 2D ARRAY
            boolMask.shape

        except:
            arr = np.frombuffer(boolMask, dtype=np.uint8)
            arr.shape = (frame.data.shape)
            boolMask = arr

        frame.mask = boolMask

        if save:
            outDir = self.intermediateRootPath
        else:
            outDir = self.intermediateRootPath + "/tmp"

        # INJECT THE PRE KEYWORD
        utcnow = datetime.utcnow()
        frame.header["SXSPRE"] = (utcnow.strftime(
            "%Y-%m-%dT%H:%M:%S.%f"), "UTC timestamp")

        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        # CONVERT CCDData TO FITS HDU (INCLUDING HEADER) AND SAVE WITH PRE TAG
        # PREPENDED TO FILENAME
        basename = os.path.basename(filepath)
        filenameNoExtension = os.path.splitext(basename)[0]
        extension = os.path.splitext(basename)[1]
        filePath = outDir + "/" + \
            filenameNoExtension + "_pre" + extension

        # SAVE TO DISK
        self._write(
            frame=frame,
            filedir=outDir,
            filename=filenameNoExtension + "_pre" + extension,
            overwrite=True
        )

        self.log.debug('completed the ``_prepare_single_frame`` method')
        return filePath

    def _absolute_path(
            self,
            path):
        """*convert paths from home directories to absolute paths*

        **Key Arguments:**
            - ``path`` -- path possibly relative to home directory

        **Return:**
            - ``absolutePath`` -- absolute path

        **Usage**

        ```python
        myPath = self._absolute_path(myPath)
        ```
        """
        self.log.debug('starting the ``_absolute_path`` method')

        from os.path import expanduser
        home = expanduser("~")
        if path[0] == "~":
            path = home + "/" + path[1:]

        self.log.debug('completed the ``_absolute_path`` method')
        return path.replace("//", "/")

    def prepare_frames(
            self,
            save=False):
        """*prepare raw frames by converting pixel data from ADU to electrons and adding mask and uncertainty extensions*

        **Key Arguments:**
            - ``save`` -- save out the prepared frame to the intermediate products directory. Default False.

        **Return:**
            - ``preframes`` -- the new image collection containing the prepared frames

        **Usage**

        Usually called within a recipe class once the input frames have been selected and verified (see `soxs_mbias` code for example):

        ```python
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])
        ```
        """
        self.log.debug('starting the ``prepare_frames`` method')

        kw = self.kw

        filepaths = self.inputFrames.files_filtered(include_path=True)

        frameCount = len(filepaths)
        print("# PREPARING %(frameCount)s RAW FRAMES - TRIMMING OVERSCAN, CONVERTING TO ELECTRON COUNTS, GENERATING UNCERTAINTY MAPS AND APPENDING DEFAULT BAD-PIXEL MASK" % locals())
        preframes = []
        preframes[:] = [self._prepare_single_frame(
            frame=frame, save=save) for frame in filepaths]
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=preframes
        )
        preframes, supplementaryInput = sof.get()
        preframes.sort([kw('MJDOBS').lower()])

        print("# PREPARED FRAMES - SUMMARY")
        print(preframes.summary)

        self.log.debug('completed the ``prepare_frames`` method')
        return preframes

    def _verify_input_frames_basics(
            self):
        """*the basic verifications that needs done for all recipes*

        **Return:**
            - None

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``_verify_input_frames_basics`` method')

        kw = self.kw

        # CHECK WE ACTUALLY HAVE IMAGES
        if not len(self.inputFrames.files_filtered(include_path=True)):
            raise FileNotFoundError(
                "No image frames where passed to the recipe")

        arm = self.inputFrames.values(
            keyword=kw("SEQ_ARM").lower(), unique=True)
        # MIXED INPUT ARMS ARE BAD
        if len(arm) > 1:
            arms = " and ".join(arms)
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of %(imageTypes)s" % locals())
        else:
            self.arm = arm[0]

        # CREATE DETECTOR LOOKUP DICTIONARY - SOME VALUES CAN BE OVERWRITTEN
        # WITH WHAT IS FOUND HERE IN FITS HEADERS
        self.detectorParams = detector_lookup(
            log=self.log,
            settings=self.settings
        ).get(self.arm)

        # MIXED BINNING IS BAD
        cdelt1 = self.inputFrames.values(
            keyword=kw("CDELT1").lower(), unique=True)
        cdelt2 = self.inputFrames.values(
            keyword=kw("CDELT2").lower(), unique=True)

        if len(cdelt1) > 1 or len(cdelt2) > 1:
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of binnings" % locals())

        if cdelt1[0] and cdelt2[0]:
            self.detectorParams["binning"] = [int(cdelt2[0]), int(cdelt1[0])]

        # MIXED READOUT SPEEDS IS BAD
        readSpeed = self.inputFrames.values(
            keyword=kw("DET_READ_SPEED").lower(), unique=True)
        if len(readSpeed) > 1:
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of readout speeds" % locals())

        # MIXED GAIN SPEEDS IS BAD
        # HIERARCH ESO DET OUT1 CONAD - Electrons/ADU
        # CONAD IS REALLY GAIN AND HAS UNIT OF Electrons/ADU
        gain = self.inputFrames.values(
            keyword=kw("CONAD").lower(), unique=True)
        if len(gain) > 1:
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of gain" % locals())
        if gain[0]:
            # UVB & VIS
            self.detectorParams["gain"] = gain[0] * u.electron / u.adu
        else:
            # NIR
            self.detectorParams["gain"] = self.detectorParams[
                "gain"] * u.electron / u.adu

        # HIERARCH ESO DET OUT1 RON - Readout noise in electrons
        ron = self.inputFrames.values(
            keyword=kw("RON").lower(), unique=True)

        # MIXED NOISE
        if len(ron) > 1:
            print(self.inputFrames.summary)
            raise TypeError("Input frames are a mix of readnoise" % locals())
        if ron[0]:
            # UVB & VIS
            self.detectorParams["ron"] = ron[0] * u.electron
        else:
            # NIR
            self.detectorParams["ron"] = self.detectorParams[
                "ron"] * u.electron

        self.log.debug('completed the ``_verify_input_frames_basics`` method')
        return None

    def clean_up(
            self):
        """*remove intermediate files once recipe is complete*

        **Usage**

        ```python
        recipe.clean_up()
        ```
        """
        self.log.debug('starting the ``clean_up`` method')

        outDir = self.intermediateRootPath + "/tmp"

        try:
            shutil.rmtree(outDir)
        except:
            pass

        self.log.debug('completed the ``clean_up`` method')
        return None

    def xsh2soxs(
            self,
            frame):
        """*perform some massaging of the xshooter data so it more closely resembles soxs data -  this function can be removed once code is production ready*

        **Key Arguments:**
            - ``frame`` -- the CCDDate frame to manipulate

        **Return:**
            - ``frame`` -- the manipulated soxspipe-ready frame

        **Usage:**

        ```python
        frame = self.xsh2soxs(frame)
        ```
        """
        self.log.debug('starting the ``xsh2soxs`` method')

        kw = self.kw
        dp = self.detectorParams

        # NP ROTATION OF ARRAYS IS IN COUNTER-CLOCKWISE DIRECTION
        rotationIndex = int(dp["clockwise-rotation"] / 90.)

        if self.settings["instrument"] == "xsh" and rotationIndex > 0:
            frame.data = np.rot90(frame.data, rotationIndex)

        self.log.debug('completed the ``xsh2soxs`` method')
        return frame

    def _trim_frame(
            self,
            frame):
        """*return frame with pre-scan and overscan regions removed*

        **Key Arguments:**
            - ``frame`` -- the CCDData frame to be trimmed
        """
        self.log.debug('starting the ``_trim_frame`` method')

        kw = self.kw
        arm = self.arm
        dp = self.detectorParams

        rs, re, cs, ce = dp["science-pixels"]["rows"]["start"], dp["science-pixels"]["rows"][
            "end"], dp["science-pixels"]["columns"]["start"], dp["science-pixels"]["columns"]["end"]

        binning = dp["binning"]
        if binning[0] > 1:
            rs = int(rs / binning[0])
            re = int(re / binning[0])
        if binning[1] > 1:
            cs = int(cs / binning[0])
            ce = int(ce / binning[0])

        trimmed_frame = ccdproc.trim_image(frame[rs: re, cs: ce])

        self.log.debug('completed the ``_trim_frame`` method')
        return trimmed_frame

    def _write(
            self,
            frame,
            filedir,
            filename=False,
            overwrite=True):
        """*write frame to disk at the specified location*

        **Key Arguments:**
            - ``frame`` -- the frame to save to disk (CCDData object)
            - ``filedir`` -- the location to save the frame
            - ``filename`` -- the filename to save the file as. Default: **False** (standardised filename generated in code)
            - ``overwrite`` -- if a file exists at the filepath then choose to overwrite the file. Default: True

        **Usage:**

        Use within a recipe like so:

        ```python
        self._write(frame, filePath)
        ```
        """
        self.log.debug('starting the ``write`` method')

        if not filename:

            filename = filenamer(
                log=self.log,
                frame=frame,
                settings=self.settings
            )

        filepath = filedir + "/" + filename

        HDUList = frame.to_hdu(
            hdu_mask='QUAL', hdu_uncertainty='ERRS', hdu_flags=None)
        HDUList[0].name = "FLUX"
        HDUList.writeto(filepath, output_verify='exception',
                        overwrite=overwrite, checksum=True)

        self.log.debug('completed the ``write`` method')
        return filepath

    def clip_and_stack(
            self,
            frames,
            recipe):
        """*mean combine input frames after sigma-clipping outlying pixels using a median value with median absolute deviation (mad) as the deviation function*

        **Key Arguments:**
            - ``frames`` -- an ImageFileCollection of the frames to stack or a list of CCDData objects
            - ``recipe`` -- the name of recipe needed to read the correct settings from the yaml files

        **Return:**
            - ``combined_frame`` -- the combined master frame (with updated bad-pixel and uncertainty maps)

        **Usage:**

        This snippet can be used within the recipe code to combine individual (using bias frames as an example):

        ```python
        combined_bias_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mbias")
        ```

        ---

        ```eval_rst
        .. todo::

            - revisit error propagation when combining frames: https://github.com/thespacedoctor/soxspipe/issues/42
        ```
        """
        self.log.debug('starting the ``clip_and_stack`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams
        imageType = self.imageType

        # ALLOW FOR UNDERSCORE AND HYPHENS
        recipe = recipe.replace("soxs_", "soxs-")

        # UNPACK SETTINGS
        clipping_lower_sigma = self.settings[
            recipe]["clipping-lower-simga"]
        clipping_upper_sigma = self.settings[
            recipe]["clipping-upper-simga"]
        clipping_iteration_count = self.settings[
            recipe]["clipping-iteration-count"]

        # LIST OF CCDDATA OBJECTS NEEDED BY COMBINER OBJECT
        # ccds = [c for c in self.inputFrames.ccds()]
        if not isinstance(frames, list):
            ccds = [c for c in frames.ccds(ccd_kwargs={"hdu_uncertainty": 'ERRS',
                                                       "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]
        else:
            ccds = frames

        # COMBINE MASKS AND THEN RESET
        combinedMask = ccds[0].mask
        for c in ccds:
            combinedMask = c.mask | combinedMask
            c.mask[:, :] = False

        # COMBINER OBJECT WILL FIRST GENERATE MASKS FOR INDIVIDUAL IMAGES VIA
        # CLIPPING AND THEN COMBINE THE IMAGES WITH THE METHOD SELECTED. PIXEL
        # MASKED IN ALL INDIVIDUAL IMAGES ARE MASK IN THE FINAL COMBINED IMAGE
        combiner = Combiner(ccds)

        print(f"\n# SIGMA-CLIPPING PIXEL WITH OUTLYING VALUES IN INDIVIDUAL {imageType} FRAMES")
        # PRINT SOME INFO FOR USER
        badCount = combinedMask.sum()
        totalPixels = np.size(combinedMask)
        percent = (float(badCount) / float(totalPixels)) * 100.
        print(f"The basic bad-pixel mask for the {arm} detector {imageType} frames contains {badCount} pixels ({percent:0.2}% of all pixels)")

        # GENERATE A MASK FOR EACH OF THE INDIVIDUAL INOUT FRAMES - USING
        # MEDIAN WITH MEDIAN ABSOLUTE DEVIATION (MAD) AS THE DEVIATION FUNCTION
        old_n_masked = -1
        # THIS IS THE SUM OF BAD-PIXELS IN ALL INDIVIDUAL FRAME MASKS
        new_n_masked = combiner.data_arr.mask.sum()
        iteration = 1
        while (new_n_masked > old_n_masked and iteration <= clipping_iteration_count):
            combiner.sigma_clipping(
                low_thresh=clipping_lower_sigma, high_thresh=clipping_upper_sigma, func=np.ma.median, dev_func=mad_std)
            old_n_masked = new_n_masked
            # RECOUNT BAD-PIXELS NOW CLIPPING HAS RUN
            new_n_masked = combiner.data_arr.mask.sum()
            diff = new_n_masked - old_n_masked
            extra = ""
            if diff == 0:
                extra = " - we're done"
            print("    Clipping iteration %(iteration)s finds %(diff)s more rogue pixels in the set of input frames%(extra)s" % locals())
            iteration += 1

        # GENERATE THE COMBINED MEDIAN
        print("\n# MEAN COMBINING FRAMES - WITH UPDATED BAD-PIXEL MASKS")
        combined_frame = combiner.average_combine()

        # RECOMBINE THE COMBINED MASK FROM ABOVE
        combined_frame.mask = combined_frame.mask | combinedMask

        # MASSIVE FUDGE - NEED TO CORRECTLY WRITE THE HEADER FOR COMBINED
        # IMAGES
        combined_frame.header = ccds[0].header
        try:
            combined_frame.wcs = ccds[0].wcs
        except:
            pass
        combined_frame.header[
            kw("DPR_CATG")] = "MASTER_%(imageType)s_%(arm)s" % locals()

        # CALCULATE NEW PIXELS ADDED TO MASK
        newBadCount = combined_frame.mask.sum()
        diff = newBadCount - badCount
        print("%(diff)s new pixels made it into the combined bad-pixel map" % locals())

        self.log.debug('completed the ``clip_and_stack`` method')
        return combined_frame

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

        Within a soxspipe recipe use `subtract_calibrations` like so:

        ```python
        myCalibratedFrame = self.subtract_calibrations(
            inputFrame=inputFrameCCDObject, master_bias=masterBiasCCDObject, dark=darkCCDObject)
        ```

        ---

        ```eval_rst
        .. todo::

            - code needs written to scale dark frame to exposure time of science/calibration frame
        ```
        """
        self.log.debug('starting the ``subtract_calibrations`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        if master_bias == None:
            master_bias = False
        if dark == None:
            dark = False

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

        # DARK WITH MATCHING EXPOSURE TIME
        tolerence = 0.5
        if dark != False and (int(dark.header[kw("EXPTIME")]) < int(inputFrame.header[kw("EXPTIME")]) + tolerence) and (int(dark.header[kw("EXPTIME")]) > int(inputFrame.header[kw("EXPTIME")]) - tolerence):
            calibration_subtracted_frame = inputFrame.subtract(dark)
            calibration_subtracted_frame.header = inputFrame.header
            try:
                calibration_subtracted_frame.wcs = inputFrame.wcs
            except:
                pass

        # ONLY A MASTER BIAS FRAME, NO DARK
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
