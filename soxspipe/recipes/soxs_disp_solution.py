#!/usr/bin/env python
# encoding: utf-8
"""
*Recipe to generate a first approximation of the dispersion solution from single pinhole frames*

Author
: David Young & Marco Landoni

Date Created
: August 25, 2020
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import keyword_lookup
from .base_recipe import base_recipe

from fundamentals import tools
from builtins import object
import sys
import os
from soxspipe.commonutils import create_dispersion_map

os.environ['TERM'] = 'vt100'


class soxs_disp_solution(base_recipe):
    """
    *generate a first approximation of the dispersion solution from single pinhole frames*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``polyOrders`` -- the orders of the x-y polynomials used to fit the dispersion solution. Overrides parameters found in the yaml settings file. e.g 345400 is order_x=3, order_y=4 ,wavelength_x=5 ,wavelength_y=4. Default *False*.

    **Usage**

    ```python
    from soxspipe.recipes import soxs_disp_solution
    disp_map_path = soxs_disp_solution(
        log=log,
        settings=settings,
        inputFrames=sofPath
    ).produce_product()
    ```
    """

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[],
            verbose=False,
            overwrite=False,
            polyOrders=False

    ):
        # INHERIT INITIALISATION FROM  base_recipe
        super(soxs_disp_solution, self).__init__(
            log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-disp-solution")
        self.log = log
        log.debug("instantiating a new 'soxs_disp_solution' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-disp-solution"
        self.polyOrders = polyOrders

        if self.polyOrders:
            try:
                self.polyOrders = int(self.polyOrders)
            except:
                pass
            if not isinstance(self.polyOrders, int):
                raise TypeError("THE poly VALUE NEEDS TO BE A 4 DIGIT INTEGER")

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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_disp_solution - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # SORT IMAGE COLLECTION
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
        """*verify input frames match those required by the `soxs_disp_solution` recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()
        error = False

        if False:
            print(imageTypes)
            print(imageTech)
            print(imageCat)
            print(imageTypes[0])

        if self.arm == "NIR":
            # WANT ON AND OFF PINHOLE FRAMES
            # MIXED INPUT IMAGE TYPES ARE BAD
            if not error:
                if len(imageTypes) > 1:
                    imageTypes = " and ".join(imageTypes)
                    imageTypes = " and ".join(imageTypes)
                    error = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals()

            if not error:
                # FIX ME!
                if imageTypes[0] not in ["LAMP,FMTCHK", 'LAMP,WAVE', 'WAVE,LAMP']:
                    error = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals()

            if not error:
                for i in imageTech:
                    if i not in ['ECHELLE,PINHOLE', 'IMAGE']:
                        error = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals()

            if not error:
                for i in ['ECHELLE,PINHOLE', 'IMAGE']:
                    if i not in imageTech:
                        error = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR" % locals()

        else:
            if not error:
                for i in imageTypes:
                    # FIX ME!
                    if i not in ["LAMP,FMTCHK", 'LAMP,WAVE', 'WAVE,LAMP']:
                        error = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and a master-bias and possibly a master dark for UVB/VIS" % locals()

            if not error:
                for i in ['ECHELLE,PINHOLE']:
                    if i not in imageTech:
                        error = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and a master-bias and possibly a master dark for UVB/VIS"

            if not error:
                for i in [f"MASTER_BIAS_{self.arm}"]:
                    if i not in imageCat:
                        error = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and a master-bias and possibly a master dark for UVB/VIS"

        if error:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            self.log.print("")
            raise TypeError(error)

        self.imageType = imageTypes[0]
        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*generate a fisrt guess of the dispersion solution*

        **Return:**

        - ``productPath`` -- the path to the first guess dispersion map
        """
        self.log.debug('starting the ``produce_product`` method')

        from astropy.nddata import CCDData
        from astropy import units as u
        import pandas as pd
        from datetime import datetime

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        # self.inputFrames.summary.pprint_all()

        master_bias = False
        dark = False
        pinhole_image = False

        add_filters = {kw("PRO_CATG"): 'MASTER_BIAS_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                       hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # UVB/VIS DARK
        add_filters = {kw("PRO_CATG"): 'MASTER_DARK_' + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # NIR DARK
        if self.inst.lower() == "soxs":
            add_filters = {kw("DPR_TYPE"): 'WAVE,LAMP', kw("DPR_TECH"): 'IMAGE'}
        else:
            add_filters = {kw("DPR_TYPE"): 'LAMP,FMTCHK', kw("DPR_TECH"): 'IMAGE'}
        # from tabulate import tabulate
        # self.log.print(tabulate(self.inputFrames.summary, headers='keys', tablefmt='psql'))
        # self.log.print(self.inputFrames.files_filtered(include_path=True, **add_filters))
        # sys.exit(0)

        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        if self.inst.lower() == "soxs":
            filter_list = [
                {kw("DPR_TYPE"): 'LAMP,WAVE', kw("DPR_TECH"): 'ECHELLE,PINHOLE'},
                {kw("DPR_TYPE"): 'WAVE,LAMP', kw("DPR_TECH"): 'ECHELLE,PINHOLE'}
            ]
        else:
            filter_list = [{kw("DPR_TYPE"): 'LAMP,FMTCHK',
                            kw("DPR_TECH"): 'ECHELLE,PINHOLE'}]

        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                pinhole_image = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                             hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        self.pinholeFrame = self.detrend(
            inputFrame=pinhole_image, master_bias=master_bias, dark=dark)

        if self.settings["save-intermediate-products"]:
            outDir = self.workspaceRootPath
            filePath = self._write(
                frame=self.pinholeFrame,
                filedir=outDir,
                filename=False,
                overwrite=True,
                product=False
            )
            self.log.print(f"\nCalibrated single pinhole frame: {filePath}\n")

        if self.settings["tune-pipeline"]:
            from itertools import product
            digits = [2, 3, 4, 5, 6]
            perm = product(digits, repeat=4)
            try:
                os.remove("residuals.txt")
            except:
                pass

            # GET THE LINE DETECTION LIST BEFORE JUMPING TO PERMUTATIONS
            mapPath, mapImagePath, res_plots, qcTable, productsTable, lineDetectionTable = create_dispersion_map(
                log=self.log,
                settings=self.settings,
                recipeSettings=self.recipeSettings,
                pinholeFrame=self.pinholeFrame,
                qcTable=self.qc,
                productsTable=self.products,
                sofName=self.sofName,
            ).get()

            # CHANGE MPL BACKEND OR WE HAVE ISSUES WITH MULTIPROCESSING
            import matplotlib.pyplot as plt
            plt.switch_backend('Agg')
            from fundamentals import fmultiprocess

            permList = list(perm)

            # DEFINE AN INPUT ARRAY

            print("TUNING SOXSPIPE\n")

            results = fmultiprocess(log=self.log, function=parameterTuning,
                                    inputArray=permList, poolSize=100, timeout=3600, recipeSettings=self.recipeSettings, settings=self.settings, pinholeFrame=self.pinholeFrame, qc=self.qc, products=self.products, sofName=self.sofName, lineDetectionTable=lineDetectionTable, turnOffMP=False, mute=True, progressBar=True)
            productPath = None

        else:
            if self.polyOrders:
                self.polyOrders = str(self.polyOrders)
                self.polyOrders = [int(digit) for digit in str(self.polyOrders)]
                self.recipeSettings["order-deg"] = self.polyOrders[:2]
                self.recipeSettings["wavelength-deg"] = self.polyOrders[2:4]

            productPath, mapImagePath, res_plots, qcTable, productsTable, lineDetectionTable = create_dispersion_map(
                log=self.log,
                settings=self.settings,
                recipeSettings=self.recipeSettings,
                pinholeFrame=self.pinholeFrame,
                qcTable=self.qc,
                productsTable=self.products,
                sofName=self.sofName
            ).get()

            filename = os.path.basename(productPath)

            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

            self.products = pd.concat([self.products, productsTable])
            self.qc = pd.concat([self.qc, qcTable])

            self.dateObs = self.pinholeFrame.header[kw("DATE_OBS")]

            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": "DISP_MAP",
                "file_name": filename,
                "file_type": "FITS Table",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"{self.arm} first pass dispersion solution",
                "file_path": productPath,
                "label": "PROD"
            }).to_frame().T], ignore_index=True)

            self.report_output()
            self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method


def parameterTuning(p, log, recipeSettings, settings, pinholeFrame, qc, products, sofName, lineDetectionTable):
    """*tuning the spatial solution*        
    """

    recipeSettings["order-deg"] = list(p[:2])
    recipeSettings["wavelength-deg"] = list(p[2:4])

    from soxspipe.commonutils import create_dispersion_map
    this = create_dispersion_map(
        log=log,
        settings=settings,
        recipeSettings=recipeSettings,
        pinholeFrame=pinholeFrame,
        qcTable=qc,
        productsTable=products,
        sofName=sofName,
        create2DMap=False,
        lineDetectionTable=lineDetectionTable
    )
    productPath, mapImagePath, res_plots, qcTable, productsTable, lineDetectionTable = this.get()

    return None
