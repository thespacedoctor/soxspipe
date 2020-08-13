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


class test_detector_lookup(unittest.TestCase):

    def test_detector_lookup_function(self):

        from soxspipe.commonutils import detector_lookup
        this = detector_lookup(
            log=log,
            settings=settings
        )
        detectorDict = this._select_dictionary()
        print(this._select_dictionary())

    def test_detector_lookup_get_function(self):

        from soxspipe.commonutils import detector_lookup
        detectDict = detector_lookup(
            log=log,
            settings=settings
        ).get("NIR")
        print(detectDict)

        detectDict = detector_lookup(
            log=log,
            settings=settings
        ).get("UVB")
        print(detectDict)

        detectDict = detector_lookup(
            log=log,
            settings=settings
        ).get("VIS")
        print(detectDict)

        # HOW ABOUT LOWERCASE?
        detectDict = detector_lookup(
            log=log,
            settings=settings
        ).get("nir")
        print(detectDict)

    def test_detector_lookup_function_wrong_arm(self):

        from soxspipe.commonutils import detector_lookup
        try:
            detectDict = detector_lookup(
                log=log,
                settings=settings
            ).get("RUBBISH")
            assert False
        except Exception as e:
            assert True
            print(str(e))

    def test_detector_lookup_function_exception(self):

        from soxspipe.commonutils import detector_lookup
        try:
            this = detector_lookup(
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
