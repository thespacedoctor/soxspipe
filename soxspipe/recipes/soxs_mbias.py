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
from astropy.stats import sigma_clip, mad_std
from soxspipe.commonutils.toolkit import generic_quality_checks
from datetime import datetime
from soxspipe.commonutils import keyword_lookup
import ccdproc
from astropy import units as u
from astropy.nddata import CCDData
import math
import numpy as np
from ._base_recipe_ import _base_recipe_
from soxspipe.commonutils import set_of_files
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class soxs_mbias(_base_recipe_):
    """
    *The* `soxs_mbias` *recipe is used to generate a master-bias frame from a set of input raw bias frames. The recipe is used only for the UV-VIS arm as NIR frames have bias (and dark current) removed by subtracting an off-frame of equal expsoure length.*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*

    **Usage**

    ```python
    from soxspipe.recipes import soxs_mbias
    mbiasFrame = soxs_mbias(
        log=log,
        settings=settings,
        inputFrames=fileList
    ).produce_product()
    ```

    ---

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
            inputFrames=[],
            verbose=False

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_mbias, self).__init__(log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mbias' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-mbias"
        self.recipeSettings = settings[self.recipeName]
        # xt-self-arg-tmpx

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MBIAS - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.write("\x1b[1A\x1b[2K")
        print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # print("\n# RAW INPUT BIAS FRAMES - SUMMARY")
        # SORT IMAGE COLLECTION
        self.inputFrames.sort(['mjd-obs'])
        # print(self.inputFrames.summary, "\n")

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])

        return None

    def verify_input_frames(
            self):
        """*verify the input frame match those required by the soxs_mbias recipe*

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
        """*generate a master bias frame*

        **Return:**
            - ``productPath`` -- the path to the master bias frame
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        combined_bias_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mbias")
        self.dateObs = combined_bias_mean.header[kw("DATE_OBS")]

        rawRon, mbiasRon = self.calc_ron(combined_bias_mean)

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

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log, frame=combined_bias_mean, settings=self.settings, recipeName=self.recipeName, qcTable=self.qc)

        # DETERMINE MEDIAN BIAS LEVEL
        maskedDataArray = np.ma.array(
            combined_bias_mean.data, mask=combined_bias_mean.mask)
        medianBias = np.ma.median(maskedDataArray)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "MBIAS MEDIAN",
            "qc_value": medianBias,
            "qc_comment": "Median level of master bias",
            "qc_unit": "electrons",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }, ignore_index=True)

        # WRITE TO DISK
        productPath = self._write(
            frame=combined_bias_mean,
            filedir=self.intermediateRootPath,
            filename=False,
            overwrite=True
        )
        filename = os.path.basename(productPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.products = self.products.append({
            "soxspipe_recipe": self.recipeName,
            "product_label": "MBIAS",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} Master bias frame",
            "file_path": productPath
        }, ignore_index=True)

        self.report_output()
        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    def calc_ron(
            self,
            combined_bias_mean):
        """*calculate the read-out-noise on the raw frames*

        **Return:**
            - ``rawRon`` -- raw read-out-noise in electrons
            - ``masterRon`` -- combined read-out-noise in mbias

        **Usage:**

        ```python
        rawRon, mbiasRon = self.calc_ron(combined_bias_mean)
        ```
        """
        self.log.debug('starting the ``calc_ron`` method')

        # LIST OF CCDDATA OBJECTS
        ccds = [c for c in self.inputFrames.ccds(ccd_kwargs={
                                                 "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        # SINGLE FRAME RON
        raw_one = ccds[0].data
        raw_two = ccds[1].data
        raw_diff = raw_two - raw_one
        def imstats(dat): return (dat.min(), dat.max(), dat.mean(), dat.std())
        dmin, dmax, dmean, dstd = imstats(raw_diff)
        rawRon = dstd

        # PREDICTED MASTER NOISE
        predictedMasterRon = rawRon/math.sqrt(len(ccds))

        dmin, dmax, dmean, dstd = imstats(combined_bias_mean.data)
        masterRon = dstd

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "RON MASTER",
            "qc_value": masterRon,
            "qc_comment": "Combined RON in MBIAS",
            "qc_unit": "electrons",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }, ignore_index=True)

        self.qc = self.qc.append({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "RON DETECTOR",
            "qc_value": rawRon,
            "qc_comment": "RON in single bias",
            "qc_unit": "electrons",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }, ignore_index=True)

        self.log.debug('completed the ``calc_ron`` method')
        return rawRon, masterRon

    # use the tab-trigger below for new method
    # xt-class-method
