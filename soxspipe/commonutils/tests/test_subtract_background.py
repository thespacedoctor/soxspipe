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
from astropy.nddata import CCDData
from astropy import units as u
packageDirectory = utKit("").get_project_root()
settingsFile = packageDirectory + "/test_settings.yaml"
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

class test_subtract_background(unittest.TestCase):

    def test_subtract_background_function(self):

        flatPath = "~/xshooter-pipeline-data/unittest_data/xshooter-subtract-background/first_iteration_NIR_master_flat.fits"
        orderTable = "~/xshooter-pipeline-data/unittest_data/xshooter-subtract-background/20170819T132045_NIR_ORDER_LOCATIONS.csv"
        home = expanduser("~")
        flatPath = flatPath.replace("~", home)
        orderTable = orderTable.replace("~", home)
        flatFrame = CCDData.read(flatPath, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                 hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        from soxspipe.commonutils.subtract_background import subtract_background
        background = subtract_background(
            log=log,
            frame=flatFrame,
            orderTable=orderTable,
            settings=settings
        )
        backgroundSubtractedFrame = background.subtract()

    def test_subtract_background_function_exception(self):

        from soxspipe.commonutils.subtract_background import subtract_background
        try:
            this = subtract_background(
                log=log,
                settings=settings,
                fakeKey="break the code"
            )
            this.subtract()
            assert False
        except Exception as e:
            assert True
            print(str(e))

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
