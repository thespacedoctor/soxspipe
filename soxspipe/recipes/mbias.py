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
import ccdproc
from ccdproc import Combiner


class mbias(_base_recipe_):
    """
    *The mbias recipe*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths. Default []

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    See `produce_product` method for usage.


    ```eval_rst
    .. todo::

        - add a tutorial about ``mbias`` to documentation
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
        super(mbias, self).__init__(log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'mbias' object")
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY MBIAS - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        print("# RAW INPUT BIAS FRAMES - SUMMARY")
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
        """*verify the input frame match those required by the mbias recipe*

        **Return:**
            - ``None``

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        imageTypes = self.inputFrames.values(
            keyword='eso dpr type', unique=True)

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

        arms = self.inputFrames.values(
            keyword='eso seq arm', unique=True)
        # MIXED INPUT ARMS ARE BAD
        if len(arms) > 1:
            arms = " and ".join(arms)
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of %(imageTypes)s" % locals())

        cdelt1 = self.inputFrames.values(
            keyword='cdelt1', unique=True)
        cdelt2 = self.inputFrames.values(
            keyword='cdelt2', unique=True)
        # MIXED BINNING IS BAD
        if len(cdelt1) > 1 or len(cdelt2) > 1:
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of binnings" % locals())

        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*The code to generate the product of the mbias recipe*

        **Return:**
            - ``productPath`` -- the path to the final product

        **Usage:**

        ```python
        from soxspipe.recipes import mbias
        recipe = mbias(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        mbiasFrame = recipe.produce_product()
        ```
        """
        self.log.debug('starting the ``produce_product`` method')

        # IMAGECOLLECTION FILEPATHS
        filepaths = self.inputFrames.files_filtered(include_path=True)
        ccds = []
        firstInputheader = None
        for f in filepaths:
            ccd = CCDData.read(f, hdu=0, unit='electron', hdu_uncertainty='UNCERT',
                               hdu_mask='MASK', hdu_flags='BITMAP', key_uncertainty_type='UTYPE')
            if not firstInputheader:
                firstInputheader = ccd.header
            ccds.append(ccd)

        print("# CLIPPING PIXELS WITH EXTREME VALUES IN INDIVIDUAL FRAMES")
        # PRINT SOME INFO FOR USER
        arm = ccd.header['eso seq arm'].upper()
        badCount = ccd.mask.sum()
        totalPixels = np.size(ccd.mask)
        percent = (float(badCount) / float(totalPixels)) * 100.
        print("    The basic bad-pixel mask for the %(arm)s detector contains %(badCount)s pixels (%(percent)0.2f%% of all pixels)" % locals())

        # CREATE A COMBINER OBJECT
        # https://ccdproc.readthedocs.io/en/latest/api/ccdproc.Combiner.html#ccdproc.Combiner
        combiner = Combiner(ccds)

        # GENERATE A MASK FOR EACH OF THE INDIVIDUAL INOUT FRAMES - USING
        # MEDIAN WITH MEDIAN ABSOLUTE DEVIATION (MAD) AS THE DEVIATION FUNCTION
        old_n_masked = 0
        new_n_masked = combiner.data_arr.mask.sum()
        # print("The basic bad-pixel mask contains %(new_n_masked)s pixels" % locals())
        iteration = 1
        while (new_n_masked > old_n_masked):
            combiner.sigma_clipping(
                low_thresh=5, high_thresh=5, func=np.ma.median, dev_func=mad_std)
            old_n_masked = new_n_masked
            new_n_masked = combiner.data_arr.mask.sum()
            diff = new_n_masked - old_n_masked
            extra = ""
            if diff == 0:
                extra = " - we're done"
            print("    Clipping iteration %(iteration)s finds %(diff)s more rogue pixels in the set of input frames%(extra)s" % locals())
            iteration += 1

        # GENERATE THE COMBINED MEDIAN
        print("# MEDIAN COMBINING INPUT FRAMES - USING UPDATED BAD-PIXEL MASK")
        combined_bias_median = combiner.median_combine()
        combined_bias_median.header = firstInputheader
        combined_bias_median.header[
            "HIERARCH ESO PRO CATG"] = "MASTER_BIAS_%(arm)s" % locals()

        # CALCULATE NEW PIXELS ADDED TO MASK
        newBadCount = ccd.mask.sum()
        diff = newBadCount - badCount
        print("%(diff)s new pixels made it into the propagated bad-pixel map" % locals())

        x = combined_bias_median.header['eso det win1 binx']
        y = combined_bias_median.header['eso det win1 biny']
        productPath = self.intermediateRootPath + \
            "/master_bias_%(arm)s_%(x)sx%(y)s.fits" % locals()

        combined_bias_median.write(
            productPath,
            overwrite=True)

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
