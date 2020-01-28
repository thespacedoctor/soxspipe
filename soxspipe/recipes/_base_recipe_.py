#!/usr/bin/env python
# encoding: utf-8
"""
*The base recipe class which all other recipes inherit*

:Author:
    David Young

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
import ccdproc
from soxspipe.commonutils import set_of_files


class _base_recipe_():
    """
    The base recipe class which all other recipes inherit

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary

    **Usage:**

    ```python
    usage code 
    ```

    ---

    ```eval_rst
    .. todo::

        - add usage info, including a template to create a new recipe
        - create a sublime snippet for usage
    ```
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,

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

        # xt-self-arg-tmpx

        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # Initial Actions

        return None

    def prepare_single_frame(
            self,
            frame,
            save=False):
        """*prepare a single raw frame by converting to electron counts and adding mask and uncertainty extensions*

        **Key Arguments:**
            - ``frame`` -- the path to the frame to prepare, of a CCDData object

        **Return:**
            - ``frame`` -- the prepared frame with mask and uncertainty extensions (CCDData object) 

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

        ## Preparing the Raw SOXS Frames

        ### Trim Overscan

        The first thing we need to do is trim off the overscan area of the image.

        ### ADU to Electrons

        The first thing to do is to convert the pixel data from ADU to electron counts by multiplying each pixel value in the raw frame by the gain (`HIERARCH ESO DET OUT1 CONAD` keyword value).

        $$\rm{electron\ count} = \rm{adu\ count} \times \rm{gain}$$

        ### Generating an Uncertainty Map

        Next we need to generate the uncertainty map for the raw image and add this as the 'UNCERT' extension to the image.

        For each pixel the uncertainty is calculated as:

        $$\rm{error} = \sqrt{\rm{readnoise}^2+\rm{electron\ count}}$$

        ### Bitmap Extension

        The appropriate bitmap extension is selected and simply added as the 'FLAG' extension of the frame.

        ### Bad Pixel Mask

        The selected bitmap is converted to a boolean mask, with values >0 becoming TRUE to indicate these pixels need to be masks. All other values are set to FALSE.

        Finally the prepared frames are saved out into the intermediate frames location with the prefix `pre_`.

        Viewing the image in DS9 (using the command `ds9 -multiframe pre_filename.fits` to show all extensions as tiled frames) we can see the 'DATA', 'UNCERT' and 'MASK' extensions are now all present.

        ![](https://live.staticflickr.com/65535/49434143241_580e904616_o.png)
        """
        self.log.debug('starting the ``prepare_single_frame`` method')

        # CONVERT FILEPATH TO CCDDATA OBJECT
        filepath = frame
        if isinstance(frame, str):

            if frame[0] == "~":
                from os.path import expanduser
                home = expanduser("~")
                frame = home + "/" + frame[1:]

            # OPEN THE RAW FRAME
            frame = CCDData.read(frame, hdu=0, unit='electron', hdu_uncertainty='UNCERT',
                                 hdu_mask='MASK', hdu_flags='BITMAP', key_uncertainty_type='UTYPE')

        # CHECK THE NUMBER OF EXTENSIONS IS ONLY 1
        if len(frame.to_hdu()) > 1:
            raise TypeError("%(frame)s is not a raw frame" % locals())

        # TRIM OVERSCAN REGION
        if frame.header['HIERARCH ESO SEQ ARM'] != "NIR":
            frame = ccdproc.trim_image(frame[:, :2048])
        else:
            frame = ccdproc.trim_image(frame[:1056, :2040])

        # HIERARCH ESO DET OUT1 CONAD - Electrons/ADU
        if ("HIERARCH ESO DET OUT1 CONAD" not in frame.header) and ("HIERARCH ESO DET CHIP GAIN" not in frame.header):
            raise AttributeError(
                "'HIERARCH ESO DET OUT1 CONAD/HIERARCH ESO DET CHIP GAIN' keyword not found in %(frame)s" % locals())
        if "HIERARCH ESO DET OUT1 CONAD" in frame.header:
            gain = frame.header["HIERARCH ESO DET OUT1 CONAD"]
        else:
            gain = frame.header["HIERARCH ESO DET CHIP GAIN"
                                ]

        # HIERARCH ESO DET OUT1 RON - Readout noise in electrons
        if ("HIERARCH ESO DET OUT1 RON" not in frame.header) and ("HIERARCH ESO DET CHIP RON" not in frame.header):
            raise AttributeError(
                "'HIERARCH ESO DET OUT1 RON/HIERARCH ESO DET CHIP RON' keyword not found in %(frame)s" % locals())
        if "HIERARCH ESO DET OUT1 RON" in frame.header:
            ron = frame.header["HIERARCH ESO DET OUT1 RON"]
        else:
            ron = frame.header["HIERARCH ESO DET CHIP RON"]

        # CONVERT ADU TO ELECTRONS
        frame.data = frame.data * gain

        # GENERATE UNCERTAINTY MAP AS EXTENSION
        frame = ccdproc.create_deviation(
            frame, readnoise=ron * u.electron)

        # FIND THE APPROPRIATE BAD-PIXEL BITMAP AND APPEND AS 'FLAG' EXTENSION
        # NOTE FLAGS NOTE YET SUPPORTED BY CCDPROC THIS THIS WON'T GET SAVED OUT
        # AS AN EXTENSION
        arm = frame.header["HIERARCH ESO SEQ ARM"]
        if "NIR" in arm:
            bitMapPath = self.calibrationRootPath + \
                "/cal/BP_MAP_RP_%(arm)s.fits" % locals()
        else:
            binx = int(frame.header["HIERARCH ESO DET WIN1 BINX"])
            biny = int(frame.header["HIERARCH ESO DET WIN1 BINY"])
            bitMapPath = self.calibrationRootPath + \
                "/cal/BP_MAP_RP_%(arm)s_%(binx)sx%(biny)s.fits" % locals()

        bitMap = CCDData.read(bitMapPath, hdu=0, unit='electron')
        if "NIR" in arm:
            bitMap.data = np.rot90(bitMap.data, -1)
        frame.flags = bitMap.data

        # FLATTEN BAD-PIXEL BITMAP TO BOOLEAN FALSE (GOOD) OR TRUE (BAD) AND
        # APPEND AS 'UNCERT' EXTENSION
        boolMask = bitMap.data.astype(bool)
        frame.mask = boolMask.data

        if save:
            outDir = self.intermediateRootPath
        else:
            outDir = self.intermediateRootPath + "/tmp"

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
        frame.write(filePath, overwrite=True)

        self.log.debug('completed the ``prepare_single_frame`` method')
        return filePath

    def _absolute_path(
            self,
            path):
        """*convert paths from home directories to absolute paths*

        **Key Arguments:**
            - ``path`` -- path possibly relative to home directory
            -

        **Return:**
            - ``absolutePath`` -- absolute path

        **Usage:**

        ```python
        myPath = self._absolute_path(myPath)
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
        """*prepare all frames in the input data*

        See `prepare_single_frame` for details.

        **Key Arguments:**
            - ``save`` -- save out the prepared frame to the intermediate products directory. Default False.

        **Return:**
            - ``preframes`` -- the new image collection containing the prepared frames

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
        self.log.debug('starting the ``prepare_frames`` method')

        filepaths = self.inputFrames.files_filtered(include_path=True)

        frameCount = len(filepaths)
        print("# PREPARING %(frameCount)s RAW FRAMES - CONVERTING TO ELECTRON COUNTS, GENERATING UNCERTAINTY MAPS AND APPENDING DEFAULT BAD-PIXEL MASK" % locals())
        preframes = []
        preframes[:] = [self.prepare_single_frame(
            frame=frame, save=save) for frame in filepaths]
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=preframes
        )
        preframes = sof.get()

        self.inputFrames.sort(['mjd-obs'])

        self.log.debug('completed the ``prepare_frames`` method')
        return preframes

    # use the tab-trigger below for new method
    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
