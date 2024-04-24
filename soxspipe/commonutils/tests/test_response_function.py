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
# pip install -e . (in the repo dir)

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

class test_response_function(unittest.TestCase):

    import pytest

    def test_xsh_response_function_uvb_function(self):
        # print(packageDirectory)
        stdExtractionPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-reponse-function/uvb/2019.08.23T23.24.12.925_UVB_1X2_SLOW_STARE_XSHOOTER_EG_274_EXTRACTED_MERGED.fits"
        from soxspipe.commonutils import response_function
        response = response_function(
            log=log,
            settings=settings,
            stdExtractionPath=stdExtractionPath
        )
        response.get()

    def test_xsh_response_function_vis_function(self):
        # print(packageDirectory)
        stdExtractionPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-reponse-function/vis/2019.08.23T23.10.18.163_VIS_1X2_SLOW_STARE_XSHOOTER_EG_274_EXTRACTED_MERGED.fits"
        from soxspipe.commonutils import response_function
        response = response_function(
            log=log,
            settings=settings,
            stdExtractionPath=stdExtractionPath
        )
        response.get()

    def test_xsh_response_function_nir_function(self):
        # print(packageDirectory)
        stdExtractionPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-reponse-function/nir/2019.08.22T23.12.18.5011_NIR_STARE_205PT0_EG_274_EXTRACTED_MERGED.fits"
        from soxspipe.commonutils import response_function
        response = response_function(
            log=log,
            settings=settings,
            stdExtractionPath=stdExtractionPath
        )
        response.get()

    def test_soxs_response_function_function_exception(self):

        from soxspipe.commonutils import response_function
        try:
            this = response_function(
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


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(test_response_function("test_xsh_response_function_function"))
    runner = unittest.TextTestRunner()
    runner.run(suite)

    # unittest.main()

    # x-class-to-test-named-worker-function
