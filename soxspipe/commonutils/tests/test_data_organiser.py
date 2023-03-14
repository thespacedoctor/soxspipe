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
pathToInputDir = home + "/xshooter-pipeline-data/unittest_data/xsh/data-organiser/"
pathToOutputDir = home + "/xshooter-pipeline-data/unittest_data/xsh/data-organiser-output/"

# try:
#     shutil.rmtree(pathToOutputDir)
# except:
#     pass
# # COPY INPUT TO OUTPUT DIR
# shutil.copytree(pathToInputDir, pathToOutputDir)

# # Recursively create missing directories
# if not os.path.exists(pathToOutputDir):
#     os.makedirs(pathToOutputDir)


# xt-setup-unit-testing-files-and-folders
# xt-utkit-refresh-database

class test_data_organiser(unittest.TestCase):

    def test_data_organiser_function(self):

        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=log,
            settings=settings,
            rootDir=pathToOutputDir + "01_EG274"
        )
        do.sync_raw_frames()
        do.populate_product_frames_db_table()

    def test_data_organiser_function_exception(self):

        from soxspipe.commonutils import data_organiser
        try:
            this = data_organiser(
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
