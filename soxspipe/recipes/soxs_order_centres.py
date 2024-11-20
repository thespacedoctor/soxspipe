#!/usr/bin/env python
# encoding: utf-8
"""
*further constrain the first guess locations of the order centres derived in `soxs_disp_solution`*

Author
: David Young & Marco Landoni

Date Created
: September  8, 2020
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils import detect_continuum
from soxspipe.commonutils import keyword_lookup
from .base_recipe import base_recipe
from fundamentals import tools
import sys
import os
os.environ['TERM'] = 'vt100'


class soxs_order_centres(base_recipe):
    """
    *further constrain the first guess locations of the order centres derived in `soxs_disp_solution`*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``polyOrders`` -- the orders of the x-y polynomials used to fit the dispersion solution. Overrides parameters found in the yaml settings file. e.g 345400 is order_x=3, order_y=4 ,wavelength_x=5 ,wavelength_y=4. Default *False*. 

    **Usage**

    ```python
    from soxspipe.recipes import soxs_order_centres
    order_table = soxs_order_centres(
        log=log,
        settings=settings,
        inputFrames=a["inputFrames"]
    ).produce_product()
    ```
    """
    # Initialisation

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
        super(soxs_order_centres, self).__init__(
            log=log, settings=settings, inputFrames=inputFrames, overwrite=overwrite, recipeName="soxs-order-centre")
        self.log = log
        log.debug("instantiating a new 'soxs_order_centres' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.polyOrders = polyOrders

        if self.polyOrders:
            try:
                self.polyOrders = int(self.polyOrders)
            except:
                pass
            if not isinstance(self.polyOrders, int):
                raise TypeError("THE poly VALUE NEEDS TO BE A 2 DIGIT INTEGER")

        # INITIAL ACTIONS
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

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_order_centres - NO MORE, NO LESS.
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
        """*verify input frames match those required by the soxs_order_centres recipe*

        **Return:**

        - ``None``

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        kw = self.kw

        error = False

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if self.arm == "NIR":
            # WANT ON AND OFF PINHOLE FRAMES
            # MIXED INPUT IMAGE TYPES ARE BAD
            if not error:
                if len(imageTypes) > 1:
                    imageTypes = " and ".join(imageTypes)
                    erorr = "Input frames are a mix of %(imageTypes)s" % locals()

            if not error:
                if self.inst == "SOXS":
                    good = 'FLAT'
                else:
                    good = "LAMP,ORDERDEF"
                if good not in imageTypes[0]:
                    error = "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off frames and a first-guess dispersion solution table for NIR" % locals()

            if not error:
                for i in imageTech:
                    if i not in ['ECHELLE,PINHOLE', 'IMAGE']:
                        error = "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off frames a first-guess dispersion solution table for NIR" % locals()

            if not error:
                for i in [f"DISP_TAB_{self.arm}"]:
                    if i not in imageCat:
                        error = "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off frames a first-guess dispersion solution table for NIR" % locals()

        else:
            if not error:
                if self.inst == "SOXS":
                    goodList = ['FLAT,LAMP', 'LAMP,DFLAT', 'LAMP,FLAT']
                else:
                    goodList = ["LAMP,ORDERDEF", 'LAMP,DORDERDEF', 'LAMP,QORDERDEF']
                for i in imageTypes:
                    if i not in goodList:
                        error = "Input frames for soxspipe order_centres need to be single pinhole flat-lamp, a master-bias frame, a first-guess dispersion solution table and possibly a master dark for UVB/VIS. Found {i}" % locals()

            if not error:
                for i in [f"MASTER_BIAS_{self.arm}", f"DISP_TAB_{self.arm}"]:
                    if i not in imageCat:
                        error = "Input frames for soxspipe order_centres need to be single pinhole flat-lamp, a master-bias frame, a first-guess dispersion solution table and possibly a master dark for UVB/VIS." % locals()

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
        """*generate the order-table with polynomal fits of order-centres*

        **Return:**

        - ``productPath`` -- the path to the order-table
        """
        self.log.debug('starting the ``produce_product`` method')

        from astropy.nddata import CCDData
        from astropy import units as u
        import pandas as pd
        from datetime import datetime

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        productPath = None

        master_bias = False
        dark = False
        orderDef_image = False

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
        if self.inst == "SOXS":
            add_filters = {kw("DPR_TYPE"): 'FLAT,LAMP',
                           kw("DPR_TECH"): 'IMAGE'}
        else:
            add_filters = {kw("DPR_TYPE"): 'LAMP,ORDERDEF',
                           kw("DPR_TECH"): 'IMAGE'}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        if self.inst == "SOXS":
            filter_list = [
                {kw("DPR_TYPE"): 'FLAT,LAMP', kw("DPR_TECH"): 'ECHELLE,PINHOLE'},
                {kw("DPR_TYPE"): 'LAMP,FLAT', kw("DPR_TECH"): 'ECHELLE,PINHOLE'},
                {kw("DPR_TYPE"): 'LAMP,DFLAT', kw("DPR_TECH"): 'ECHELLE,PINHOLE'},
                # KEYWORD SCREW-UP DURING PAE MEANT WE HAD TO ADD BELOW WITH ECHELLE,SLIT ... SHOULD REMOVE THIS EVENTUALLY
                {kw("DPR_TYPE"): 'FLAT,LAMP', kw("DPR_TECH"): 'ECHELLE,SLIT'},
                {kw("DPR_TYPE"): 'LAMP,FLAT', kw("DPR_TECH"): 'ECHELLE,SLIT'},
                {kw("DPR_TYPE"): 'LAMP,DFLAT', kw("DPR_TECH"): 'ECHELLE,SLIT'},
            ]
        else:
            # UVB XSHOOTER - CHECK FOR D2 LAMP FIRST AND IF NOT FOUND USE THE QTH LAMP
            filter_list = [
                {kw("DPR_TYPE"): 'LAMP,ORDERDEF', kw("DPR_TECH"): 'ECHELLE,PINHOLE'},
                {kw("DPR_TYPE"): 'LAMP,QORDERDEF', kw("DPR_TECH"): 'ECHELLE,PINHOLE'},
                {kw("DPR_TYPE"): 'LAMP,DORDERDEF', kw("DPR_TECH"): 'ECHELLE,PINHOLE'}
            ]

        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                orderDef_image = CCDData.read(i, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                              hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        add_filters = {kw("PRO_CATG"): f"DISP_TAB_{arm}".upper()}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            disp_map_table = i

        self.orderFrame = self.detrend(
            inputFrame=orderDef_image, master_bias=master_bias, dark=dark)

        if self.settings["save-intermediate-products"]:
            fileDir = self.workspaceRootPath
            filepath = self._write(
                self.orderFrame, fileDir, filename=False, overwrite=True, product=False)
            self.log.print(f"\nCalibrated single pinhole frame frame saved to {filepath}\n")

        # FIND THE APPROPRIATE PREDICTED LINE-LIST
        if arm != "NIR" and kw('WIN_BINX') in self.orderFrame.header:
            binx = int(self.orderFrame.header[kw('WIN_BINX')])
            biny = int(self.orderFrame.header[kw('WIN_BINY')])
        else:
            binx = 1
            biny = 1

        if self.settings["tune-pipeline"]:
            from itertools import product
            digits = [2, 3, 4, 5, 6]
            perm = product(digits, repeat=2)
            try:
                os.remove("residuals.txt")
            except:
                pass

            # DETECT THE CONTINUUM OF ORDERE CENTRES - RETURN ORDER TABLE FILE PATH
            # self.log.print("\n# DETECTING ORDER CENTRE CONTINUUM\n")
            detector = detect_continuum(
                log=self.log,
                traceFrame=self.orderFrame,
                dispersion_map=disp_map_table,
                settings=self.settings,
                recipeSettings=self.recipeSettings,
                recipeName="soxs-order-centre",
                qcTable=self.qc,
                productsTable=self.products,
                sofName=self.sofName,
                binx=binx,
                biny=biny
            )
            orderPixelTable = detector.sample_trace()

            print("\n\nTUNING SOXSPIPE\n")

            # CHANGE MPL BACKEND OR WE HAVE ISSUES WITH MULTIPROCESSING
            import matplotlib.pyplot as plt
            plt.switch_backend('Agg')
            # DEFINE AN INPUT ARRAY
            from fundamentals import fmultiprocess
            # NOTE TO SELF: if having issue with multiprocessing stalling, try and import required modules into the mthod/function running this fmultiprocess function instead of at the module level
            results = fmultiprocess(log=self.log, function=parameterTuning,
                                    inputArray=list(perm), poolSize=False, timeout=360000, recipeSettings=self.recipeSettings, settings=self.settings, orderFrame=self.orderFrame, disp_map_table=disp_map_table, orderPixelTable=orderPixelTable, qc=self.qc, products=self.products, sofName=self.sofName, binx=binx, biny=biny, turnOffMP=False, mute=True, progressBar=True)
            return None

        else:
            if self.polyOrders:
                self.polyOrders = str(self.polyOrders)
                self.polyOrders = [int(digit) for digit in str(self.polyOrders)]
                self.recipeSettings["detect-continuum"]["order-deg"] = self.polyOrders[0]
                self.recipeSettings["detect-continuum"]["disp-axis-deg"] = self.polyOrders[1]

            # DETECT THE CONTINUUM OF ORDERE CENTRES - RETURN ORDER TABLE FILE PATH
            # self.log.print("\n# DETECTING ORDER CENTRE CONTINUUM\n")
            detector = detect_continuum(
                log=self.log,
                traceFrame=self.orderFrame,
                dispersion_map=disp_map_table,
                settings=self.settings,
                recipeSettings=self.recipeSettings,
                recipeName="soxs-order-centre",
                qcTable=self.qc,
                productsTable=self.products,
                sofName=self.sofName,
                binx=binx,
                biny=biny
            )
            productPath, qcTable, productsTable, orderPolyTable, orderPixelTable, orderMetaTable = detector.get()

        self.products = pd.concat([self.products, productsTable])
        self.qc = pd.concat([self.qc, qcTable])

        filename = os.path.basename(productPath)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
        self.dateObs = self.orderFrame.header[kw("DATE_OBS")]

        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "ORDER_CENTRES",
            "file_name": filename,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} order centre traces",
            "file_path": productPath,
            "label": "PROD"
        }).to_frame().T], ignore_index=True)

        self.report_output()
        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath


def parameterTuning(p, log, recipeSettings, settings, orderFrame, disp_map_table, orderPixelTable, qc, products, sofName, binx, biny):
    """*tuning the spatial solution*        
    """

    recipeSettings["detect-continuum"]["order-deg"] = p[0]
    recipeSettings["detect-continuum"]["disp-axis-deg"] = p[1]

    from soxspipe.commonutils import create_dispersion_map
    detector = detect_continuum(
        log=log,
        traceFrame=orderFrame,
        dispersion_map=disp_map_table,
        settings=settings,
        recipeSettings=recipeSettings,
        recipeName="soxs-order-centre",
        qcTable=qc,
        productsTable=products,
        sofName=sofName,
        binx=binx,
        biny=biny,
        orderPixelTable=orderPixelTable
    )
    try:
        productPath, qcTable, productsTable, orderPolyTable, orderPixelTable, orderMetaTable = detector.get()
    except:
        pass

    return None

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
