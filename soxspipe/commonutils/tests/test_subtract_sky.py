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

    def test_2D_image_to_DF_function(self):

        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_single.fits"
        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_multi.fits"
        twoDMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/20190830T184348_NIR_2D_MAP_IMAGE.fits"

        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
        mapDF = twoD_disp_map_image_to_dataframe(log=log, twoDMapPath=twoDMap)

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
        mapDF = twoD_disp_map_image_to_dataframe(log=log, twoDMapPath=twoDMap, assosiatedFrame=objectFrame)

        from tabulate import tabulate
        print(tabulate(mapDF.head(100), headers='keys', tablefmt='psql'))

    def test_subtract_sky_function(self):

        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_multi.fits"
        objectPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/stare_mode_cal_single.fits"
        twoDMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/20190830T184348_NIR_2D_MAP_IMAGE.fits"
        dispMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-subtract-sky/20190830T184348_NIR_2D_MAP.fits"

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

        from soxspipe.commonutils import subtract_sky
        this = subtract_sky(
            log=log,
            settings=settings,
            objectFrame=objectFrame,
            twoDMap=twoDMap,
            dispMap=dispMap
        )
        this.get()

    def test_subtract_sky_function_exception(self):

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
