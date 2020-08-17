#!/usr/bin/env python
# encoding: utf-8
"""
*The recipe for creating master-bias frames *

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
from soxspipe.commonutils import set_of_files
from ._base_recipe_ import _base_recipe_
from astropy.stats import mad_std
import numpy as np
from astropy.nddata import CCDData
from astropy import units as u
import ccdproc
from ccdproc import Combiner
from soxspipe.commonutils import keyword_lookup


class soxs_mbias(_base_recipe_):
    """
    *The `soxs_mbias` recipe is used to generate a master-bias frame from a set of input raw bias frames. The recipe is used only for the UV-VIS arm as NIR frames have the bias (and dark current) removed by subtracting an off-frame of equal expsoure length.*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths. Default []

    **Usage**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    See `produce_product` method for usage.


    ```eval_rst
    .. todo::

        - add a tutorial about ``soxs_mbias`` to documentation
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
        super(soxs_mbias, self).__init__(log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mbias' object")
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MBIAS - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        print("\n# RAW INPUT BIAS FRAMES - SUMMARY")
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
        """*verify the input frame match those required by the soxs_mbias recipe*

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
        # MIXED INPUT IMAGE TYPES ARE BAD
        if len(imageTypes) > 1:
            imageTypes = " and ".join(imageTypes)
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of %(imageTypes)s" % locals())
        # NON-BIAS INPUT IMAGE TYPES ARE BAD
        elif imageTypes[0] != 'BIAS':
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames not BIAS frames" % locals())

        self.imageType = imageTypes[0]

        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*The code to generate the product of the soxs_mbias recipe*

        **Return:**
            - ``productPath`` -- the path to the final product

        **Usage**

        ```python
        from soxspipe.recipes import soxs_mbias
        recipe = soxs_mbias(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        mbiasFrame = recipe.produce_product()
        ```
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        combined_bias_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mbias")

        x = combined_bias_mean.header[kw("WIN_BINX")]
        y = combined_bias_mean.header[kw("WIN_BINY")]
        productPath = self.intermediateRootPath + \
            "/master_bias_%(arm)s_%(x)sx%(y)s.fits" % locals()

        # INSPECTING THE THE UNCERTAINTY MAPS
        # print("individual frame data")
        # for a in ccds:
        #     print(a.data[0][0])
        # print("\ncombined frame data")
        # print(combined_bias_mean.data[0][0])
        # print("individual frame error")
        # for a in ccds:
        #     print(a.uncertainty[0][0])
        # print("combined frame error")
        # print(combined_bias_mean.uncertainty[0][0])

        # combined_bias_mean.data = combined_bias_mean.data.astype('float32')
        # combined_bias_mean.uncertainty = combined_bias_mean.uncertainty.astype(
        #     'float32')

        # WRITE TO DISK
        self.write(combined_bias_mean, productPath, overwrite=True)
        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    def clip_and_stack(
            self,
            frames,
            recipe):
        """*mean combine input frames after sigma-clipping outlying pixels using a median value with median absolute deviation (mad) as the deviation function*

        **Key Arguments:**
            - ``frames`` -- an ImageFileCollection of the framers to stack
            - ``recipe`` -- the name of recipe needed to read the correct settings from the yaml files

        **Return:**
            - ``combined_frame`` -- the combined master frame (with updated bad-pixel and uncertainty maps)

        ```eval_rst
        .. todo::

        - revisit error propagation when combining frames: https://github.com/thespacedoctor/soxspipe/issues/42
        ```

        **Usage:**

        This snippet can be used within the recipe code to combine individual (using bias frames as an example):

        ```python
        combined_bias_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mbias")
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
        ccds = [c for c in self.inputFrames.ccds(ccd_kwargs={"hdu_uncertainty": 'ERRS',
                                                             "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        # COMBINER OBJECT WILL FIRST GENERATE MASKS FOR INDIVIDUAL IMAGES VIA
        # CLIPPING AND THEN COMBINE THE IMAGES WITH THE METHOD SELECTED. PIXEL
        # MASKED IN ALL INDIVIDUAL IMAGES ARE MASK IN THE FINAL COMBINED IMAGE
        combiner = Combiner(ccds)

        print(f"\n# SIGMA-CLIPPING PIXEL WITH OUTLYING VALUES IN INDIVIDUAL {imageType} FRAMES")
        # PRINT SOME INFO FOR USER
        badCount = ccds[0].mask.sum()
        totalPixels = np.size(ccds[0].mask)
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

        # MASSIVE FUDGE - NEED TO CORRECTLY WRITE THE HEADER FOR COMBINED
        # IMAGES
        combined_frame.header = ccds[0].header
        combined_frame.header[
            kw("DPR_CATG")] = "MASTER_%(imageType)s_%(arm)s" % locals()

        # CALCULATE NEW PIXELS ADDED TO MASK
        newBadCount = combined_frame.mask.sum()
        diff = newBadCount - badCount
        print("%(diff)s new pixels made it into the combined bad-pixel map" % locals())

        self.log.debug('completed the ``clip_and_stack`` method')
        return combined_frame

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
