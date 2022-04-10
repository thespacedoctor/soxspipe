#!/usr/bin/env python
# encoding: utf-8
"""
*enhance the wavelength solution achieved with `soxs_disp_solution` by expanding the solution into the spatial dimension (along the slit)*

:Author:
    David Young & Marco Landoni

:Date Created:
    March 17, 2021
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
from astropy import units as u
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


class soxs_spatial_solution(_base_recipe_):
    """
    *The soxs_spatial_solution recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths
        - ``verbose`` -- verbose. True or False. Default *False*

    See `produce_product` method for usage.

    ```eval_rst
    .. todo::

        - add a tutorial about ``soxs_spatial_solution`` to documentation
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
        super(soxs_spatial_solution, self).__init__(
            log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_spatial_solution' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-spatial-solution"
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_spatial_solution - NO MORE, NO LESS.
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
        """*verify input frames match those required by the `soxs_spatial_solution` recipe*

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if self.arm == "NIR":
            # WANT ON AND OFF PINHOLE FRAMES
            for i in imageTypes:
                if i not in ["LAMP,WAVE", "LAMP,FLAT"]:
                    raise TypeError(
                        f"Found a {i} file. Input frames for soxspipe spatial_solution need to be LAMP,WAVE. Can optionally supply a master-flat for NIR.")

            for i in imageTech:
                if i not in ['ECHELLE,MULTI-PINHOLE', 'IMAGE', 'ECHELLE,SLIT', 'ECHELLE,PINHOLE']:
                    raise TypeError(
                        f"Found a {i} file. Input frames for soxspipe spatial_solution need to be LAMP,WAVE. Can optionally supply a master-flat for NIR.")

            if "LAMP,WAVE" not in imageTypes:
                raise TypeError(
                    "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames for NIR" % locals())

            if "ECHELLE,MULTI-PINHOLE" not in imageTech:
                raise TypeError(
                    "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames for NIR" % locals())

        else:
            for i in imageTypes:
                if i not in ["LAMP,WAVE", "BIAS", "DARK", "LAMP,FLAT"]:
                    raise TypeError(
                        f"Found a {i} file. Input frames for soxspipe spatial_solution need to be LAMP,WAVE and a master-bias. Can optionally supply a master-flat and/or master-dark for UVB/VIS.")

        # LOOK FOR DISP MAP
        arm = self.arm
        if f"DISP_TAB_{arm}" not in imageCat:
            raise TypeError(
                "Need a first guess dispersion map for %(arm)s - none found with the input files" % locals())

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
        """*generate the 2D dispersion map*

        **Return:**
            - ``productPath`` -- the path to the 2D dispersion map

        **Usage**

        ```python
        from soxspipe.recipes import soxs_spatial_solution
        recipe = soxs_spatial_solution(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        disp_map = recipe.produce_product()
        ```
        """
        self.log.debug('starting the ``produce_product`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None

        master_bias = False
        dark = False
        master_flat = False
        multi_pinhole_image = False
        order_table = False

        add_filters = {kw("DPR_CATG"): 'MASTER_BIAS_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS DARK
        add_filters = {kw("DPR_CATG"): 'MASTER_DARK_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # NIR DARK
        add_filters = {kw("DPR_TYPE"): 'LAMP,WAVE',
                       kw("DPR_TECH"): 'IMAGE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS/NIR FLAT
        add_filters = {kw("DPR_CATG"): 'MASTER_LAMP-FLAT_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_flat = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # MULTIPINHOLE IMAGE
        add_filters = {kw("DPR_TYPE"): 'LAMP,WAVE',
                       kw("DPR_TECH"): 'ECHELLE,MULTI-PINHOLE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            multi_pinhole_image = CCDData.read(i, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                               hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # FIND THE ORDER TABLE
        filterDict = {kw("PRO_CATG").lower(): f"ORDER_TAB_{arm}"}
        order_table = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        add_filters = {kw("PRO_CATG"): f"DISP_TAB_{arm}".upper()}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            disp_map_table = i

        self.multiPinholeFrame = self.detrend(
            inputFrame=multi_pinhole_image, master_bias=master_bias, dark=dark, master_flat=master_flat, order_table=order_table)

        if self.settings["save-intermediate-products"]:
            fileDir = self.intermediateRootPath
            filepath = self._write(
                self.multiPinholeFrame, fileDir, filename=False, overwrite=True)
            print(f"\nCalibrated multi pinhole frame frame saved to {filepath}\n")

        # GENERATE AN UPDATED DISPERSION MAP
        from soxspipe.commonutils import create_dispersion_map
        mapPath, mapImagePath = create_dispersion_map(
            log=self.log,
            settings=self.settings,
            pinholeFrame=self.multiPinholeFrame,
            firstGuessMap=disp_map_table,
            orderTable=order_table,
            qcTable=False,
            productsTable=False
        ).get()

        from datetime import datetime
        filename = os.path.basename(mapPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.products = self.products.append({
            "soxspipe_recipe": self.recipeName,
            "product_label": "SPAT_SOL",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} full dispersion-spatial solution",
            "file_path": productPath
        }, ignore_index=True)

        filename = os.path.basename(mapImagePath)
        self.products = self.products.append({
            "soxspipe_recipe": self.recipeName,
            "product_label": "2D_MAP",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} 2D detector map of wavelength, slit position and order",
            "file_path": productPath
        }, ignore_index=True)

        self.report_output()

        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return mapPath, mapImagePath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
