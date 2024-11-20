#!/usr/bin/env python
# encoding: utf-8
"""
*The base recipe class which all other recipes inherit*

Author
: David Young & Marco Landoni

Date Created
: January 22, 2020
"""
################# GLOBAL IMPORTS ####################

from soxspipe.commonutils import filenamer
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils import subtract_background
from fundamentals import tools
from builtins import object
import sys
import os

os.environ['TERM'] = 'vt100'


class base_recipe(object):
    """
    The base recipe class which all other recipes inherit

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``recipeName`` -- name of the recipe as it appears in the settings dictionary. Default *False*

    **Usage**

    To use this base recipe to create a new `soxspipe` recipe, have a look at the code for one of the simpler recipes (e.g. `soxs_mbias`) - copy and modify the code.
    """

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=False,
            verbose=False,
            overwrite=False,
            recipeName=False
    ):
        import yaml
        import pandas as pd
        from soxspipe.commonutils import toolkit
        import sqlite3 as sql

        log.debug("instantiating a new '__init__' object")
        self.recipeName = recipeName
        self.settings = settings
        self.workspaceRootPath = self._absolute_path(
            settings["workspace-root-dir"])

        # CHECK IF PRODUCT ALREADY EXISTS
        if inputFrames and not isinstance(inputFrames, list) and inputFrames.split(".")[-1].lower() == "sof":
            self.sofName = os.path.basename(inputFrames).replace(".sof", "")
            self.productPath = toolkit.predict_product_path(inputFrames, self.recipeName)
            if os.path.exists(self.productPath) and not overwrite:
                print(f"The product of this recipe already exists at '{self.productPath}'. To overwrite this product, rerun the pipeline command with the overwrite flag (-x).")
                raise FileExistsError
            self.log = toolkit.add_recipe_logger(log, self.productPath)
        else:
            self.sofName = False
            self.productPath = False
            self.log = log

        from soxspipe.commonutils.toolkit import get_calibrations_path
        self.calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)

        self.verbose = verbose
        # SET LATER WHEN VERIFYING FRAMES
        self.arm = None
        self.detectorParams = None
        self.dateObs = None

        # FIND THE CURRENT SESSION
        from os.path import expanduser
        home = expanduser("~")
        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=self.log,
            rootDir=self.settings["workspace-root-dir"].replace("~", home)
        )
        self.currentSession, allSessions = do.session_list(silent=True)

        # INITIATE A DB CONNECTION
        self.conn = None
        if self.currentSession and self.sofName:
            self.sessionDb = self.settings["workspace-root-dir"].replace("~", home) + "/soxspipe.db"

            def dict_factory(cursor, row):
                d = {}
                for idx, col in enumerate(cursor.description):
                    d[col[0]] = row[idx]
                return d
            self.conn = sql.connect(self.sessionDb, isolation_level=None)
            self.conn.row_factory = dict_factory

        # SET RECIPE TO 'FAIL' AND SWITCH TO 'PASS' ONLY IF RECIPE COMPLETES
        if self.conn:
            c = self.conn.cursor()
            sqlQuery = f"update product_frames set status_{self.currentSession} = 'fail' where sof = '{self.sofName}.sof'"
            c.execute(sqlQuery)
            c.close()

        # COLLECT ADVANCED SETTINGS IF AVAILABLE
        parentDirectory = os.path.dirname(__file__)
        advs = parentDirectory + "/advanced_settings.yaml"
        level = 0
        exists = False
        count = 1
        while not exists and len(advs) and count < 10:
            count += 1
            level -= 1
            exists = os.path.exists(advs)
            if not exists:
                advs = "/".join(parentDirectory.split("/")
                                [:level]) + "/advanced_settings.yaml"
        if not exists:
            advs = {}
        else:
            with open(advs, 'r') as stream:
                advs = yaml.safe_load(stream)
        # MERGE ADVANCED SETTINGS AND USER SETTINGS (USER SETTINGS OVERRIDE)
        self.settings = {**advs, **self.settings}

        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        self.qc = pd.DataFrame({
            "soxspipe_recipe": [],
            "qc_name": [],
            "qc_value": [],
            "qc_unit": [],
            "qc_comment": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "to_header": []
        })
        self.products = pd.DataFrame({
            "soxspipe_recipe": [],
            "product_label": [],
            "file_name": [],
            "file_type": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "file_path": [],
            "label": []
        })

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get

        self.qcDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}/"
        self.productDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}/"

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

        :::{todo}
            - write a command-line tool for this method
        :::
        """
        self.log.debug('starting the ``_prepare_single_frame`` method')

        from astropy.nddata import CCDData
        import ccdproc
        from astropy import units as u
        import numpy as np
        import logging
        import warnings
        from datetime import datetime

        warnings.filterwarnings(
            action='ignore'
        )
        logging.captureWarnings(True)

        kw = self.kw
        dp = self.detectorParams

        # STORE FILEPATH FOR LATER USE
        filepath = frame

        # CONVERT FILEPATH TO CCDDATA OBJECT
        if isinstance(frame, str):
            # CONVERT RELATIVE TO ABSOLUTE PATHS
            frame = self._absolute_path(frame)
            # OPEN THE RAW FRAME - MASK AND UNCERT TO BE POPULATED LATER
            try:
                frame = CCDData.read(frame, hdu=0, unit=u.adu, hdu_uncertainty='ERRS',
                                     hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
            except TypeError as e:
                self.log.info(f"{filepath} is a FITS Binary Table")
                return filepath

        # CHECK THE NUMBER OF EXTENSIONS IS ONLY 1 AND "SXSPRE" DOES NOT
        # EXIST. i.e. THIS IS A RAW UNTOUCHED FRAME
        if len(frame.to_hdu()) > 1 or "SXSPRE" in frame.header:
            return filepath

        # MANIPULATE XSH DATA
        frame = self.xsh2soxs(frame)
        frame = self._trim_frame(frame)

        # FUDGE BAD-PIXEL MAP CREATION
        # bpmData = np.random.rand(frame.data.shape[0], frame.data.shape[1])
        # bpmData[bpmData < 0.995] = 0
        # bpmData[bpmData > 0] = 1
        # bpmData = CCDData(bpmData, unit="adu")
        # from soxspipe.commonutils.toolkit import quicklook_image
        # quicklook_image(
        #     log=self.log, CCDObject=bpmData, show=False, ext='data', stdWindow=3)
        # HDUList = bpmData.to_hdu()
        # HDUList[0].name = "BPM"
        # HDUList.writeto("/Users/Dave/Desktop/tmp.fits", output_verify='exception',
        #                 overwrite=True, checksum=True)

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
                frame, readnoise=dp["ron"], disregard_nan=True)

        # FIND THE APPROPRIATE BAD-PIXEL BITMAP AND APPEND AS 'FLAG' EXTENSION
        # NOTE FLAGS NOT YET SUPPORTED BY CCDPROC THIS THIS WON'T GET SAVED OUT
        # AS AN EXTENSION
        arm = self.arm
        if arm != "NIR" and kw('WIN_BINX') in frame.header:
            binx = int(frame.header[kw('WIN_BINX')])
            biny = int(frame.header[kw('WIN_BINY')])
        else:
            binx = 1
            biny = 1

        bitMapPath = self.calibrationRootPath + "/" + dp["bad-pixel map"][f"{binx}x{biny}"]

        if not os.path.exists(bitMapPath):
            message = "the path to the bitMapPath %s does not exist on this machine" % (
                bitMapPath,)

            if True:
                # CREATE A DUMMY BAD-PIXEL MAP
                import numpy as np
                from astropy.nddata import CCDData
                frame = CCDData(np.full_like(frame.data, 0), unit="adu")
                # WRITE CCDDATA OBJECT TO FILE
                HDUList = frame.to_hdu()
                HDUList.writeto(bitMapPath, output_verify='exception',
                                overwrite=True, checksum=True)

            self.log.critical(message)
            raise IOError(message)

        bitMap = CCDData.read(bitMapPath, hdu=0, unit=u.dimensionless_unscaled)

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
            outDir = self.workspaceRootPath
        else:
            outDir = self.workspaceRootPath + "/tmp"

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
            overwrite=True,
            product=False
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

        from os.path import expanduser
        home = expanduser("~")
        if path[0] == "~":
            path = home + "/" + path[1:]

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

        from soxspipe.commonutils.set_of_files import set_of_files
        import numpy as np

        kw = self.kw

        filepaths = self.inputFrames.files_filtered(include_path=True)

        frameCount = len(filepaths)

        self.log.print("\n# PREPARING %(frameCount)s RAW FRAMES - TRIMMING OVERSCAN, CONVERTING TO ELECTRON COUNTS, GENERATING UNCERTAINTY MAPS AND APPENDING DEFAULT BAD-PIXEL MASK" % locals())
        preframes = []
        preframes[:] = [self._prepare_single_frame(
            frame=frame, save=save) for frame in filepaths]

        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=preframes,
            recipeName=self.recipeName,
            verbose=self.verbose
        )
        preframes, supplementaryInput = sof.get()
        preframes.sort([kw('MJDOBS')])

        self.log.print("# PREPARED FRAMES - SUMMARY")

        slitname = kw(f"SLIT_{self.arm}".upper())
        try:
            preframes.summary["SLIT"] = preframes.summary[slitname]
        except:
            pass

        preframes.summary["LAMP"] = "------------"
        columns = preframes.summary.colnames

        for i in range(7):
            thisLamp = kw(f"LAMP{i+1}")
            # FIRST FIND THE NAME OF THE LAMP
            newLamp = preframes.summary[thisLamp][np.where(preframes.summary[thisLamp].filled(999) != 999)]
            if len(newLamp):
                newLamp = newLamp[0]
                newLamp = newLamp.replace("_Lamp", "")
                newLamp = newLamp.replace("Argo", "Ar").replace("Neon", "Ne").replace("Merc", "Hg").replace("Xeno", "Xe")

                updatedList = list(preframes.summary["LAMP"][np.where(preframes.summary[thisLamp].filled(999) != 999)].data)
                updatedList[:] = [u.replace("-", "") + newLamp for u in updatedList]
                preframes.summary["LAMP"][np.where(preframes.summary[thisLamp].filled(999) != 999)] = updatedList
            columns.remove(thisLamp)

        preframes.summary["LAMP"][np.where(preframes.summary["LAMP"] == "------------")] = "--"

        try:
            columns.remove(kw("SLIT_NIR"))
        except:
            pass
        try:
            columns.remove(kw("SLIT_VIS"))
        except:
            pass
        try:
            columns.remove(kw("SLIT_UVB"))
        except:
            pass

        if "filename" in columns:
            columns.remove("file")
            columns.remove("filename")
            columns = ["filename"] + columns

        # MAKE A COPY TO RENAME COLUMNS
        this = preframes.summary.copy()

        newColumns = []
        for c in columns:
            newC = c
            if "ESO SEQ " in c:
                newC = c.replace("ESO SEQ ", "")
                this.rename_column(c, newC)
            if "ESO DPR " in c:
                newC = c.replace("ESO DPR ", "")
                this.rename_column(c, newC)
            newColumns.append(newC)

        self.log.print(this[newColumns])
        self.log.print("\n")

        # SORT RECIPE AND ARM SETTINGS
        self.recipeSettings = self.get_recipe_settings()

        self.log.debug('completed the ``prepare_frames`` method')
        return preframes

    def _verify_input_frames_basics(
            self):
        """*the basic verifications that needs done for all recipes*

        **Return:**

        - None

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug('starting the ``_verify_input_frames_basics`` method')

        from astropy import units as u
        from contextlib import suppress

        kw = self.kw

        # CHECK WE ACTUALLY HAVE IMAGES
        if not len(self.inputFrames.files_filtered(include_path=True)):
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            raise FileNotFoundError(
                "No image frames where passed to the recipe")

        arm = self.inputFrames.values(
            keyword=kw("SEQ_ARM"), unique=True)

        inst = self.inputFrames.values(kw("INSTRUME"), unique=True)
        with suppress(ValueError):
            inst.remove(None)
        self.inst = inst[0]

        # MIXED INPUT ARMS ARE BAD
        if None in arm:
            arm.remove(None)
        if len(arm) > 1:
            arms = " and ".join(arm)
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
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

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

        # MIXED BINNING IS BAD
        if self.arm == "NIR":
            # NIR ARRAY NEVER BINNED
            cdelt1 = [1]
            cdelt2 = [1]
        else:
            cdelt1 = self.inputFrames.values(
                keyword=kw("CDELT1"), unique=True)
            cdelt2 = self.inputFrames.values(
                keyword=kw("CDELT2"), unique=True)
            try:
                cdelt1.remove(None)
                cdelt2.remove(None)
            except:
                pass

        if len(cdelt1) > 1 or len(cdelt2) > 1:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            raise TypeError(
                "Input frames are a mix of binnings" % locals())

        if cdelt1[0] and cdelt2[0]:
            self.detectorParams["binning"] = [int(cdelt2[0]), int(cdelt1[0])]

        # MIXED READOUT SPEEDS IS BAD
        readSpeed = self.inputFrames.values(
            keyword=kw("DET_READ_SPEED"), unique=True)

        with suppress(ValueError):
            readSpeed.remove(None)

        if len(readSpeed) > 1:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            raise TypeError(
                f"Input frames are a mix of readout speeds. {readSpeed}" % locals())

        # MIXED GAIN SPEEDS IS BAD
        # HIERARCH ESO DET OUT1 CONAD - Electrons/ADU
        # CONAD IS REALLY GAIN AND HAS UNIT OF Electrons/ADU
        if self.settings["instrument"] == "xsh":
            gain = self.inputFrames.values(
                keyword=kw("CONAD"), unique=True)
        else:
            gain = self.inputFrames.values(
                keyword=kw("GAIN"), unique=True)

        with suppress(ValueError):
            gain.remove(None)

        if len(gain) > 1:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of gain" % locals())
        if len(gain) and gain[0]:
            # UVB & VIS
            self.detectorParams["gain"] = gain[0] * u.electron / u.adu
        else:
            # NIR
            self.detectorParams["gain"] = self.detectorParams[
                "gain"] * u.electron / u.adu

        # CONVERT TO DATAFRAME AND FILTER TO CHECK SLIT WIDTHS
        filteredDf = self.inputFrames.summary.to_pandas()
        matchList = ["LAMP,FLAT", "LAMP,QFLAT", "LAMP,DFLAT", "OBJECT", "STD,FLUX", f"MASTER_FLAT_{self.arm}"]
        mask = (filteredDf[kw("DPR_TYPE")].isin(matchList)) | (filteredDf[kw("PRO_CATG")].isin(matchList))
        filteredDf = filteredDf.loc[mask]

        # MIXED SLIT-WIDTH IS BAD
        slitWidth = list(filteredDf[kw(f"SLIT_{self.arm}")].unique())
        with suppress(ValueError):
            slitWidth.remove(None)

        # ITEMS TO REMOVE
        removeItems = ["pin", "blind"]
        for i in slitWidth:
            for j in removeItems:
                if j in str(i).lower():
                    slitWidth.remove(i)

        slitWidth = [str(s).replace("JH", "") for s in slitWidth]
        slitWidth = list(set(slitWidth))

        if len(slitWidth) > 1:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            raise TypeError(
                f"Input frames are a mix of slit-width ({slitWidth})" % locals())

        # HIERARCH ESO DET OUT1 RON - Readout noise in electrons
        ron = self.inputFrames.values(
            keyword=kw("RON"), unique=True)
        with suppress(ValueError):
            ron.remove(None)

        # MIXED NOISE
        if len(ron) > 1:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            raise TypeError(f"Input frames are a mix of readnoise. {ron}" % locals())
        if len(ron) and ron[0]:
            # UVB & VIS
            self.detectorParams["ron"] = ron[0] * u.electron
        else:
            # NIR
            self.detectorParams["ron"] = self.detectorParams[
                "ron"] * u.electron

        imageTypes = self.inputFrames.values(
            keyword=kw("DPR_TYPE"), unique=True) + self.inputFrames.values(
            keyword=kw("PRO_TYPE"), unique=True)
        imageTech = self.inputFrames.values(
            keyword=kw("DPR_TECH"), unique=True) + self.inputFrames.values(
            keyword=kw("PRO_TECH"), unique=True)
        imageCat = self.inputFrames.values(
            keyword=kw("DPR_CATG"), unique=True) + self.inputFrames.values(
            keyword=kw("PRO_CATG"), unique=True)

        def clean_list(myList):
            myList = list(set(myList))
            try:
                myList.remove(None)
            except:
                pass
            try:
                myList.remove("REDUCED")
            except:
                pass

            return myList

        imageTypes = clean_list(imageTypes)
        imageTech = clean_list(imageTech)
        imageCat = clean_list(imageCat)

        self.log.debug('completed the ``_verify_input_frames_basics`` method')
        return imageTypes, imageTech, imageCat

    def clean_up(
            self):
        """*update product status in DB and remove intermediate files once recipe is complete*

        **Usage**

        ```python
        recipe.clean_up()
        ```
        """
        self.log.debug('starting the ``clean_up`` method')

        import shutil

        # SET RECIPE PRODUCTS TO 'PASS'
        if self.conn:
            c = self.conn.cursor()
            sqlQuery = f"update product_frames set status_{self.currentSession} = 'pass' where sof = '{self.sofName}.sof'"
            c.execute(sqlQuery)
            c.close()

        outDir = self.workspaceRootPath + "/tmp"

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
        import numpy as np

        kw = self.kw
        dp = self.detectorParams

        # NP ROTATION OF ARRAYS IS IN COUNTER-CLOCKWISE DIRECTION
        rotationIndex = int(dp["clockwise-rotation"] / 90.)

        if rotationIndex > 0:
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

        import ccdproc

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
            overwrite=True,
            product=True,
            maskToZero=False):
        """*write frame to disk at the specified location*

        **Key Arguments:**

        - ``frame`` -- the frame to save to disk (CCDData object)
        - ``filedir`` -- the location to save the frame
        - ``filename`` -- the filename to save the file as. Default: **False** (standardised filename generated in code)
        - ``overwrite`` -- if a file exists at the filepath then choose to overwrite the file. Default: True
        - ``product`` -- is this a recipe product?
        - ``maskToZero`` -- set masked pixels to zero before writing to file?

        **Usage:**

        Use within a recipe like so:

        ```python
        self._write(frame, filePath)
        ```
        """
        self.log.debug('starting the ``write`` method')

        kw = self.kw

        # WRITE QCs TO HEADERS
        for n, v, c, h in zip(self.qc["qc_name"].values, self.qc["qc_value"].values, self.qc["qc_comment"].values, self.qc["to_header"].values):
            if h:
                frame.header[f"ESO QC {n}".upper()] = (v, c)

        # NEATLY SORT KEYWORDS
        keywords = [k for k in frame.header if len(k)]
        values = [frame.header[k] for k in frame.header if len(k)]
        comments = [frame.header.comments[k] for k in frame.header if len(k)]
        keywords, values, comments = zip(
            *sorted(zip(keywords, values, comments)))
        if "COMMENT" not in keywords and "HISTORY" not in keywords:
            frame.header.clear()
            for k, v, c in zip(keywords, values, comments):
                if k == "COMMENT":
                    frame.header[k] = v
                else:
                    frame.header[k] = (v, c)

        if not filename and self.sofName:
            filename = self.sofName + ".fits"
        if not filename:

            filename = filenamer(
                log=self.log,
                frame=frame,
                settings=self.settings
            )

        if product:
            removeKw = ["DPR_TECH", "DPR_CATG", "DPR_TYPE"]
            for k in removeKw:
                try:
                    frame.header.pop(kw(k))
                except:
                    pass

            filedir += f"/product/{self.recipeName}/"
            filedir = filedir.replace("//", "/")
            # Recursively create missing directories
            if not os.path.exists(filedir):
                os.makedirs(filedir)

        filepath = filedir + "/" + filename

        # SET BAD-PIXELS TO 0 IN DATA FRAME
        if maskToZero:
            self.log.print(f"\nSetting {frame.mask.sum()} bad-pixels to a value of 0 while saving '{filename}'.")
            frame.data[frame.mask] = 1

        HDUList = frame.to_hdu(
            hdu_mask='QUAL', hdu_uncertainty='ERRS', hdu_flags=None)
        HDUList[0].name = "FLUX"
        if product:
            HDUList.writeto(filepath, output_verify='fix+warn',
                            overwrite=overwrite, checksum=True)
        else:
            HDUList.writeto(filepath,
                            overwrite=overwrite, checksum=False)

        filepath = os.path.abspath(filepath)

        self.log.debug('completed the ``write`` method')
        return filepath

    def clip_and_stack(
            self,
            frames,
            recipe,
            ignore_input_masks=False,
            post_stack_clipping=True):
        """*mean combine input frames after sigma-clipping outlying pixels using a median value with median absolute deviation (mad) as the deviation function*

        **Key Arguments:**

        - ``frames`` -- an ImageFileCollection of the frames to stack or a list of CCDData objects
        - ``recipe`` -- the name of recipe needed to read the correct settings from the yaml files
        - ``ignore_input_masks`` -- ignore the input masks during clip and stacking?
        - ``post_stack_clipping`` -- allow cross-plane clipping on combined frame. Clipping settings in setting file. Default *True*.

        **Return:**

        - ``combined_frame`` -- the combined master frame (with updated bad-pixel and uncertainty maps)

        **Usage:**

        This snippet can be used within the recipe code to combine individual (using bias frames as an example):

        ```python
        combined_bias_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mbias", ignore_input_masks=False, post_stack_clipping=True)
        ```
        """
        self.log.debug('starting the ``clip_and_stack`` method')

        from astropy.stats import sigma_clip, mad_std
        from soxspipe.commonutils.combiner import Combiner
        from astropy import units as u
        import numpy as np

        if len(frames) == 1:
            self.log.info("Only 1 frame was sent to the clip and stack method. Returning the frame with no further processing.")
            return frames[0]

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams
        imageType = self.imageType

        # ALLOW FOR UNDERSCORE AND HYPHENS
        recipe = recipe.replace("soxs_", "soxs-")

        # UNPACK SETTINGS
        stacked_clipping_sigma = self.recipeSettings["stacked-clipping-sigma"]
        stacked_clipping_iterations = self.recipeSettings["stacked-clipping-iterations"]

        # LIST OF CCDDATA OBJECTS NEEDED BY COMBINER OBJECT
        if not isinstance(frames, list):
            ccds = [c for c in frames.ccds(ccd_kwargs={"hdu_uncertainty": 'ERRS',
                                                       "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE', "unit": u.electron})]
        else:
            ccds = frames

        imageType = ccds[0].header[kw("DPR_TYPE")].replace(",", "-")
        imageTech = ccds[0].header[kw("DPR_TECH")].replace(",", "-")
        imageCat = ccds[0].header[kw("DPR_CATG")].replace(",", "-")

        self.log.print(f"\n# MEAN COMBINING {len(ccds)} {arm} {imageCat} {imageTech} {imageType} FRAMES")

        # COMBINE MASKS AND THEN RESET
        combinedMask = ccds[0].mask
        for c in ccds:
            combinedMask = c.mask | combinedMask
            if ignore_input_masks:
                c.mask[:, :] = False

        # COMBINER OBJECT WILL FIRST GENERATE MASKS FOR INDIVIDUAL IMAGES VIA
        # CLIPPING AND THEN COMBINE THE IMAGES WITH THE METHOD SELECTED. PIXEL
        # MASKED IN ALL INDIVIDUAL IMAGES ARE MASK IN THE FINAL COMBINED IMAGE
        combiner = Combiner(ccds)

        # self.log.print(f"\n# SIGMA-CLIPPING PIXEL WITH OUTLYING VALUES IN INDIVIDUAL {imageType} FRAMES")
        # PRINT SOME INFO FOR USER
        badCount = combinedMask.sum()
        totalPixels = np.size(combinedMask)
        percent = (float(badCount) / float(totalPixels)) * 100.
        if imageType != "BIAS":
            self.log.print(f"\tThe basic bad-pixel mask for the {arm} detector {imageType} frames contains {badCount} pixels ({percent:0.2}% of all pixels)")

        # GENERATE A MASK FOR EACH OF THE INDIVIDUAL INPUT FRAMES - USING
        # MEDIAN WITH MEDIAN ABSOLUTE DEVIATION (MAD) AS THE DEVIATION FUNCTION
        old_n_masked = -1
        # THIS IS THE SUM OF BAD-PIXELS IN ALL INDIVIDUAL FRAME MASKS
        new_n_masked = combiner.data_arr.mask.sum()

        # SIGMA CLIPPING OVERWRITES ORIGINAL MASKS - COPY HERE TO READD
        # preclipped_masks = np.copy(combiner.data_arr.mask)
        totalPixels = np.size(combinedMask)

        combiner.data_arr.mask = sigma_clip(combiner.data_arr.data,
                                            sigma_lower=stacked_clipping_sigma,
                                            sigma_upper=stacked_clipping_sigma,
                                            axis=0,
                                            copy=False,
                                            maxiters=stacked_clipping_iterations,
                                            cenfunc='median',
                                            stdfunc='mad_std',
                                            masked=True).mask
        old_n_masked = new_n_masked
        # RECOUNT BAD-PIXELS NOW CLIPPING HAS RUN
        new_n_masked = combiner.data_arr.mask.sum()
        diff = new_n_masked - old_n_masked
        if self.verbose:
            percent = 100 * combiner.data_arr.mask[0].sum() / totalPixels
            self.log.print(f"\tClipping found {diff} more rogue pixels in the set of all input frames (~{percent:0.2}% per-frame)")

        # GENERATE THE COMBINED MEAN
        # self.log.print("\n# MEAN COMBINING FRAMES - WITH UPDATED BAD-PIXEL MASKS")
        combined_frame = combiner.average_combine()

        # RECOMBINE THE COMBINED MASK FROM ABOVE
        combined_frame.mask = combined_frame.mask | combinedMask

        # INVIDUAL UPDATED MASKS (POST CLIPPING)
        new_individual_masks = combiner.data_arr.mask
        masked_values = new_individual_masks.sum(axis=0)

        # A HACK TO THE COMBINER OBJECT TO COMBINE ERROR MAPS EXACTLY AS DATA WAS COMBINED
        for i, ccd in enumerate(ccds):
            combiner.data_arr.data[i] = ccd.uncertainty.array
        combined_uncertainty = combiner.average_combine()
        combined_frame.uncertainty = combined_uncertainty.data / (np.sqrt(len(new_individual_masks) - masked_values))

        # MASSIVE FUDGE - NEED TO CORRECTLY WRITE THE HEADER FOR COMBINED
        # IMAGES
        combined_frame.header = ccds[0].header
        try:
            combined_frame.wcs = ccds[0].wcs
        except:
            pass

        # CALCULATE NEW PIXELS ADDED TO MASK
        if imageType != "BIAS":
            newBadCount = combined_frame.mask.sum()
            diff = newBadCount - badCount
            totalPixels = np.size(combinedMask)
            percent = (float(newBadCount) / float(totalPixels)) * 100.
            self.log.print(f"\t{diff} new pixels made it into the combined bad-pixel map (bad pixels now account for {percent:0.2f}% of all pixels)")

        self.log.debug('completed the ``clip_and_stack`` method')
        return combined_frame

    def detrend(
            self,
            inputFrame,
            master_bias=False,
            dark=False,
            master_flat=False,
            order_table=False):
        """*subtract calibration frames from an input frame*

        **Key Arguments:**

        - ``inputFrame`` -- the input frame to have calibrations subtracted. CCDData object.
        - ``master_bias`` -- the master bias frame to be subtracted. CCDData object. Default *False*.
        - ``dark`` -- a dark frame to be subtracted. CCDData object. Default *False*.
        - ``master_flat`` -- divided input frame by this master flat frame. CCDData object. Default *False*.
        - ``order_table`` -- order table with order edges defined. Used to subtract scattered light background from frames. Default *False*.

        **Return:**

        - ``calibration_subtracted_frame`` -- the input frame with the calibration frame(s) subtracted. CCDData object.

        **Usage:**

        Within a soxspipe recipe use `detrend` like so:

        ```python
        myCalibratedFrame = self.detrend(
            inputFrame=inputFrameCCDObject, master_bias=masterBiasCCDObject, dark=darkCCDObject)
        ```
        """
        self.log.debug('starting the ``detrend`` method')

        import ccdproc
        from astropy import units as u
        import copy
        from astropy.io import fits
        import pandas as pd
        from datetime import datetime
        from os.path import expanduser

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        if master_bias == None:
            master_bias = False
        if dark == None:
            dark = False

        # VERIFY DATA IS IN ORDER
        if master_bias == False and dark == False and master_flat == False:
            raise TypeError(
                "detrend method needs at least a master-bias frame, a dark frame or a master flat frame")
        if master_bias == False and dark != False and dark.header[kw("EXPTIME")] != inputFrame.header[kw("EXPTIME")]:
            self.log.warning("Dark and science/calibration frame have differing exposure-times.")
        if master_bias != False and dark != False and dark.header[kw("EXPTIME")] != inputFrame.header[kw("EXPTIME")]:
            raise AttributeError(
                "CODE NEEDS WRITTEN HERE TO SCALE DARK FRAME TO EXPOSURE TIME OF SCIENCE/CALIBRATION FRAME")

        processedFrame = inputFrame.copy()

        if master_bias != False:
            processedFrame = ccdproc.subtract_bias(processedFrame, master_bias)

        # DARK WITH MATCHING EXPOSURE TIME
        tolerence = 0.5
        if dark != False and (int(dark.header[kw("EXPTIME")]) < int(processedFrame.header[kw("EXPTIME")]) + tolerence) and (int(dark.header[kw("EXPTIME")]) > int(processedFrame.header[kw("EXPTIME")]) - tolerence):
            processedFrame = ccdproc.subtract_bias(processedFrame, dark)
        elif dark != False:
            self.log.print(f"Scaling the dark to the exposure time of {inputFrame.header[kw('EXPTIME')]}s")
            processedFrame = ccdproc.subtract_dark(processedFrame, dark, exposure_time=kw("EXPTIME"), exposure_unit=u.second, scale=True)

        doSubtraction = True
        if "subtract_background" in self.recipeSettings and not self.recipeSettings["subtract_background"]:
            doSubtraction = False

        if order_table != False and doSubtraction:

            background = subtract_background(
                log=self.log,
                frame=processedFrame,
                sofName=self.sofName,
                recipeName=self.recipeName,
                orderTable=order_table,
                settings=self.settings,
                productsTable=self.products,
                qcTable=self.qc
            )
            backgroundFrame, processedFrame, self.products = background.subtract()

            from soxspipe.commonutils.toolkit import quicklook_image
            quicklook_image(
                log=self.log, CCDObject=backgroundFrame, show=False, ext='data', stdWindow=3, title="Background Light", surfacePlot=True)

            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

            # WRITE FITS FRAME OF BACKGROUND IMAGE ... PDF BEING GENERATED INSTEAD
            if False:
                # DETERMINE WHERE TO WRITE THE FILE
                home = expanduser("~")

                if self.currentSession:
                    outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/sessions/{self.currentSession}/qc/{self.recipeName}"
                else:
                    outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}"
                outDir = outDir.replace("//", "/")
                # RECURSIVELY CREATE MISSING DIRECTORIES
                if not os.path.exists(outDir):
                    os.makedirs(outDir)

                # GET THE EXTENSION (WITH DOT PREFIX)
                filename = self.sofName + "_BKGROUND.fits"
                filepath = f"{outDir}/{filename}"
                header = copy.deepcopy(inputFrame.header)
                primary_hdu = fits.PrimaryHDU(backgroundFrame.data, header=header)
                hdul = fits.HDUList([primary_hdu])
                hdul.writeto(filepath, output_verify='exception',
                             overwrite=True, checksum=True)

                self.products = pd.concat([self.products, pd.Series({
                    "soxspipe_recipe": self.recipeName,
                    "product_label": "BKGROUND",
                    "file_name": filename,
                    "file_type": "FITS",
                    "obs_date_utc": self.dateObs,
                    "reduction_date_utc": utcnow,
                    "product_desc": f"Fitted intra-order image background",
                    "file_path": filepath,
                    "label": "QC"
                }).to_frame().T], ignore_index=True)

        if master_flat != False:
            processedFrame = ccdproc.flat_correct(processedFrame, master_flat)

        self.log.debug('completed the ``detrend`` method')
        return processedFrame

    def report_output(
            self,
            rformat="stdout"):
        """*a method to report QC values alongside intermediate and final products*

        **Key Arguments:**

        - ``rformat`` -- the format to outout reports as. Default *stdout*. [stdout|....]

        **Usage:**

        ```python
        self.report_output(rformat="stdout")
        ```
        """
        self.log.debug('starting the ``report_output`` method')

        from tabulate import tabulate

        if not self.verbose:
            # REMOVE COLUMN FROM DATA FRAME
            self.products.drop(columns=['file_path'], inplace=True)

        # REMOVE DUPLICATE ENTRIES IN COLUMN 'qc_name' AND KEEP THE LAST ENTRY
        self.qc = self.qc.drop_duplicates(subset=['qc_name'], keep='last')
        # SORT BY COLUMN NAME
        self.qc.sort_values(['qc_name'], inplace=True)
        columns = list(self.qc.columns)
        columns.remove("to_header")
        columns.remove("obs_date_utc")
        columns.remove("reduction_date_utc")
        columns.remove("soxspipe_recipe")

        # SORT BY COLUMN NAME
        self.products.sort_values(['label'],
                                  ascending=[True], inplace=True)

        self.products.drop_duplicates(inplace=True)

        columns2 = list(self.products.columns)
        columns2.remove("reduction_date_utc")
        columns2.remove("soxspipe_recipe")

        try:
            soxspipe_recipe = self.qc["soxspipe_recipe"].values[0].upper()
        except:
            soxspipe_recipe = self.recipeName.upper()

        if rformat == "stdout":
            self.log.print(f"\n# {soxspipe_recipe} QC METRICS")
            self.log.print(tabulate(self.qc[columns], headers='keys', tablefmt='psql', showindex=False, stralign="right"))
            self.log.print(f"\n# {soxspipe_recipe} RECIPE PRODUCTS & QC OUTPUTS")
            self.log.print(tabulate(self.products[columns2], headers='keys', tablefmt='psql', showindex=False, stralign="right"))

        self.log.debug('completed the ``report_output`` method')
        return None

    def qc_ron(
            self,
            frameType=False,
            frameName=False,
            masterFrame=False,
            rawRon=False,
            masterRon=False):
        """*calculate the read-out-noise from bias/dark frames*

        **Key Arguments:**

        - ``frameType`` -- the type of the frame for reporting QC values. Default *False*
        - ``frameName`` -- the name of the frame in human readable words. Default *False*
        - ``masterFrame`` -- the master frame (only makes sense to measure RON on master bias). Default *False*
        - ``rawRon`` -- if serendipitously calculated elsewhere don't recalculate. Default *False*
        - ``masterRon`` -- if serendipitously calculated elsewhere don't recalculate. Default *False*

        **Return:**

        - ``rawRon`` -- raw read-out-noise in electrons
        - ``masterRon`` -- combined read-out-noise in mbias

        **Usage:**

        ```python
        rawRon, mbiasRon = self.qc_ron(
            frameType="MBIAS",
            frameName="master bias",
            masterFrame=masterFrame
        )
        ```
        """
        self.log.debug('starting the ``qc_bias_ron`` method')

        from astropy.stats import sigma_clip
        import numpy as np
        import pandas as pd
        import math
        from datetime import datetime

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        if not rawRon and len(self.inputFrames.files) > 1:
            # LIST OF RAW CCDDATA OBJECTS
            ccds = [c for c in self.inputFrames.ccds(ccd_kwargs={
                "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

            # SINGLE FRAME RON
            raw_one = ccds[0]
            raw_two = ccds[1]
            raw_diff = raw_one.subtract(raw_two)

            # SIGMA-CLIP THE DATA (AT HIGH LEVEL)
            masked_diff = sigma_clip(
                raw_diff, sigma_lower=10, sigma_upper=10, maxiters=2, cenfunc='median', stdfunc='mad_std')
            combinedMask = raw_diff.mask | masked_diff.mask

            # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
            raw_diff = np.ma.array(raw_diff.data, mask=combinedMask)

            def imstats(dat): return (dat.min(), dat.max(), dat.mean(), dat.std())
            dmin, dmax, dmean, dstd = imstats(raw_diff)

            # ACCOUNT FOR EXTRA NOISE ADDED FROM SUBTRACTING FRAMES
            rawRon = dstd / math.sqrt(2)

        if rawRon:
            singleFrameType = frameType
            if frameType[0] == "M":
                singleFrameType = frameType[1:]

            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "RON DETECTOR",
                "qc_value": rawRon,
                "qc_comment": f"[e-] RON in single {singleFrameType}",
                "qc_unit": "electrons",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

        if masterFrame and not masterRon:

            # PREDICTED MASTER NOISE
            # predictedMasterRon = rawRon / math.sqrt(len(ccds))

            # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
            tmp = np.ma.array(masterFrame.data, mask=combinedMask)

            dmin, dmax, dmean, dstd = imstats(tmp)
            masterRon = dstd

        elif masterRon:
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "RON MASTER",
                "qc_value": masterRon,
                "qc_comment": f"[e-] Combined RON in {frameType}",
                "qc_unit": "electrons",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
        else:
            masterRon = None

        self.log.debug('completed the ``qc_bias_ron`` method')
        return rawRon, masterRon

    def qc_median_flux_level(
            self,
            frame,
            frameType="MBIAS",
            frameName="master bias",
            medianFlux=False):
        """*calculate the median flux level in the frame, excluding masked pixels*

        **Key Arguments:**

        - ``frame`` -- the frame (CCDData object) to determine the median level.
        - ``frameType`` -- the type of the frame for reporting QC values Default "MBIAS"
        - ``frameName`` -- the name of the frame in human readable words. Default "master bias"
        - ``medianFlux`` -- if serendipitously calculated elsewhere don't recalculate. Default *False*

        **Return:**

        - ``medianFlux`` -- median flux level in electrons

        **Usage:**

        ```python
        medianFlux = self.qc_median_flux_level(
            frame=myFrame,
            frameType="MBIAS",
            frameName="master bias")
        ```
        """
        self.log.debug('starting the ``qc_median_flux_level`` method')

        import numpy as np
        import pandas as pd
        from datetime import datetime

        if not medianFlux:
            # DETERMINE MEDIAN BIAS LEVEL
            maskedDataArray = np.ma.array(
                frame.data, mask=frame.mask)
            medianFlux = np.ma.median(maskedDataArray)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": f"{frameType} MEDIAN".upper(),
            "qc_value": medianFlux,
            "qc_comment": f"[e-] Median flux level of {frameName}",
            "qc_unit": "electrons",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``qc_median_flux_level`` method')
        return medianFlux

    def subtract_mean_flux_level(
            self,
            rawFrame):
        """*iteratively median sigma-clip raw bias data frames before calculating and removing the mean bias level*

        **Key Arguments:**

        - ``rawFrame`` -- the raw bias frame

        **Return:**

        - `meanFluxLevel` -- the frame mean bias level
        - `fluxStd` -- the standard deviation of the flux distribution (RON)
        - `noiseFrame` -- the raw bias frame with mean bias level removed

        **Usage:**

        ```python
        meanFluxLevel, fluxStd, noiseFrame = self.subtract_mean_flux_level(rawFrame)
        ```
        """
        self.log.debug('starting the ``subtract_mean_flux_level`` method')

        from astropy.stats import sigma_clip, mad_std
        import numpy as np

        # UNPACK SETTINGS
        clipping_lower_sigma = self.recipeSettings["frame-clipping-sigma"]
        clipping_upper_sigma = clipping_lower_sigma
        clipping_iteration_count = self.recipeSettings["frame-clipping-iterations"]

        maskedFrame = sigma_clip(
            rawFrame, sigma=clipping_lower_sigma, maxiters=clipping_iteration_count, cenfunc='median', stdfunc='mad_std')

        # DETERMINE MEDIAN BIAS LEVEL
        maskedDataArray = np.ma.array(
            maskedFrame.data, mask=maskedFrame.mask)
        meanFluxLevel = np.ma.mean(maskedDataArray)
        fluxStd = np.ma.std(maskedDataArray)
        rawFrame.data -= meanFluxLevel

        self.log.debug('completed the ``subtract_mean_flux_level`` method')
        return (meanFluxLevel, fluxStd, rawFrame)

    def update_fits_keywords(
            self,
            frame):
        """*update fits keywords to comply with ESO Phase 3 standards*

        **Key Arguments:**

        - ``frame`` -- the frame to update

        **Return:**

        - None

        **Usage:**

        ```python
        usage code
        ```

        :::{todo}
            - add usage info
            - create a sublime snippet for usage
            - write a command-line tool for this method
            - update package tutorial with command-line tool info if needed
        :::
        """
        self.log.debug('starting the ``update_fits_keywords`` method')

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams
        imageType = self.imageType
        if "FLAT" in imageType:
            imageType = "FLAT"

        frame.header[kw("SEQ_ARM").upper()] = arm
        frame.header[kw("PRO_TYPE").upper()] = "REDUCED"

        # PROD CATG
        if imageType in ["BIAS", "DARK", "FLAT"]:
            frame.header[
                kw("PRO_CATG")] = f"MASTER_{imageType}_{arm}".replace("QLAMP", "LAMP").replace("DLAMP", "LAMP")
            frame.header[
                kw("PRO_TECH")] = "IMAGE"

        self.log.debug('completed the ``update_fits_keywords`` method')
        return None

    def get_recipe_settings(
            self):
        """*get the recipe and arm specific settings*

        **Return:**

        - ``recipeSettings`` -- the recipe specific settings

        **Usage:**

        ```python
        usage code
        ```
        """
        self.log.debug('starting the ``get_recipe_settings`` method')

        recipeSettings = False
        if self.recipeName:
            recipeSettings = self.settings[self.recipeName]
        if recipeSettings and self.arm and self.arm.lower() in recipeSettings:
            for k, v in recipeSettings[self.arm.lower()].items():
                recipeSettings[k] = v

        self.log.debug('completed the ``get_recipe_settings`` method')
        return recipeSettings

    # use the tab-trigger below for new method
    # xt-class-method
