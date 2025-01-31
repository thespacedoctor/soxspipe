from __future__ import print_function
from builtins import str
import os
import unittest
import shutil
import unittest
import yaml
from soxspipe.utKit import utKit
from fundamentals import tools
from os.path import expanduser
from soxspipe.commonutils.filenamer import filenamer
from os.path import expanduser
import pandas as pd
home = expanduser("~")

packageDirectory = utKit("").get_project_root()
settingsFile = packageDirectory + "/test_settings_xsh.yaml"
# settingsFile = home + \
#     "/git_repos/_misc_/settings/soxspipe/test_settings.yaml"

su = tools(
    arguments={"settingsFile": settingsFile},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName=None,
    defaultSettingsFile=False
)
arguments, settings, log, dbConn = su.setup()

# SETUP PATHS TO COMMON DIRECTORIES FOR TEST DATA
moduleDirectory = os.path.dirname(__file__)
pathToInputDir = moduleDirectory + "/input/"
pathToOutputDir = moduleDirectory + "/output/"

try:
    shutil.rmtree(pathToOutputDir)
except:
    pass
# COPY INPUT TO OUTPUT DIR
shutil.copytree(pathToInputDir, pathToOutputDir)

# Recursively create missing directories
if not os.path.exists(pathToOutputDir):
    os.makedirs(pathToOutputDir)


# xt-setup-unit-testing-files-and-folders
# xt-utkit-refresh-database

class test_subtract_sky(unittest.TestCase):

    import pytest

    @pytest.mark.full
    def test_xsh_2D_image_to_DF_function(self):

        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_single.fits"
        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_multi.fits"
        twoDMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/20190830T184348_NIR_2D_MAP_IMAGE.fits"

        from soxspipe.commonutils import keyword_lookup
        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=log,
            settings=settings
        ).get

        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
        mapDF, interOrderMask = twoD_disp_map_image_to_dataframe(log=log, twoDMapPath=twoDMap, slit_length=11, kw=kw)

        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
        from astropy.nddata import CCDData
        from astropy import units as u
        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser
        home = expanduser("~")
        if objectPath[0] == "~":
            objectPath = objectPath.replace("~", home)

        objectFrame = CCDData.read(objectPath, hdu=0, unit=u.electron, hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
        mapDF, interOrderMask = twoD_disp_map_image_to_dataframe(log=log, twoDMapPath=twoDMap, associatedFrame=objectFrame, slit_length=11, kw=kw)

        from tabulate import tabulate
        print(tabulate(mapDF.head(100), headers='keys', tablefmt='psql'))

    @pytest.mark.full
    def test_xsh_subtract_sky_function(self):

        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_single.fits"
        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_multi.fits"
        twoDMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/20190830T184348_NIR_2D_MAP_IMAGE.fits"
        dispMap = "~/xshooter-pipeline-dTata/unittest_data/xsh/xshooter-subtract-sky/20190830T184348_NIR_2D_MAP.fits"

        # UNIT READ FROM BUNIT KEYWORD OF FITS FILE UNLESS EXPLICITLY SUPPLIED
        # CCDDATA BEHAVES LIKE A NUMPY (MASKED IF MASK SET) ARRAY
        # UNMASKED DATA ACCESSED VIA frame.data
        from astropy.nddata import CCDData
        from astropy import units as u
        from os.path import expanduser
        home = expanduser("~")
        if objectPath[0] == "~":
            objectPath = objectPath.replace("~", home)

        objectFrame = CCDData.read(objectPath, hdu=0, unit=u.electron, hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
        # objectFrame.data.byteswap().newbyteorder() ## REMOVE IF BELOW .astype(float) AND .astype(bool) WORK
        # objectFrame.mask.byteswap().newbyteorder()
        objectFrame.data.astype(float)
        objectFrame.mask.astype(bool)

        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        qc = pd.DataFrame({
            "soxspipe_recipe": [],
            "qc_name": [],
            "qc_value": [],
            "qc_unit": [],
            "qc_comment": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "to_header": []
        })
        products = pd.DataFrame({
            "soxspipe_recipe": [],
            "product_label": [],
            "file_name": [],
            "file_type": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "file_path": [],
            "label": []
        })

        from soxspipe.commonutils import subtract_sky
        this = subtract_sky(
            log=log,
            settings=settings,
            recipeSettings=settings["soxs-stare"],
            objectFrame=objectFrame,
            twoDMap=twoDMap,
            qcTable=qc,
            productsTable=products,
            dispMap=dispMap
        )
        skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData, qcTable, productsTable = this.subtract()

        from tabulate import tabulate
        print(tabulate(productsTable, headers='keys', tablefmt='psql'))

    @pytest.mark.full
    def test_soxs_subtract_sky_function_exception(self):

        from soxspipe.commonutils import subtract_sky
        try:
            this = subtract_sky(
                log=log,
                settings=settings,
                fakeKey="break the code"
            )
            this.get()
            assert False
        except Exception as e:
            assert True
            print(str(e))

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
