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

class test_flux_calibration(unittest.TestCase):

    def test_flux_calibration_function(self):

        responseFunction = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-fluxcal/nir/PLACEHOLDER.fits"
        extractedSpectrum = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-fluxcal/nir/PLACEHOLDER.fits"

        from soxspipe.commonutils import flux_calibration
        fluxcal = flux_calibration(
            log=log,
            settings=settings,
            responseFunction=responseFunction,
            extractedSpectrum=extractedSpectrum
        )
        fluxcal.calibrate()

    def test_flux_calibration_function_exception(self):

        from soxspipe.commonutils import flux_calibration
        try:
            this = flux_calibration(
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
