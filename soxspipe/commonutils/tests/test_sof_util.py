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

packageDirectory = utKit("", dbConn=False).get_project_root()
settingsFile = packageDirectory + "/test_settings.yaml"
# settingsFile = home + "/.config/soxspipe/soxspipe.yaml"
su = tools(
    arguments={"settingsFile": settingsFile},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName=None,
    defaultSettingsFile=False
)
arguments, settings, log, dbConn = su.setup()

# SETUP AND TEARDOWN FIXTURE FUNCTIONS FOR THE ENTIRE MODULE
moduleDirectory = os.path.dirname(__file__)
utKit = utKit(moduleDirectory)
log, dbConn, pathToInputDir, pathToOutputDir = utKit.setupModule()
utKit.tearDownModule()

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


class test_sof_util(unittest.TestCase):

    def test_sof_util_function(self):

        # utKit.refresh_database() # reset database to database found in
        # soxspipe/test/input
        directory = settings["test-data-root"] + "/xshooter-bias/vis"
        other_output = settings[
            "reduced-data-root"].replace("reduced", "other_output")

        sofPath = other_output + "/test.sof"
        from soxspipe.commonutils import sof_util
        sof = sof_util(
            log=log,
            settings=settings
        )
        sofFile = sof.generate_sof_file_from_directory(
            directory=directory, sofPath=sofPath)
        print("sof file written to %(sofPath)s" % locals())

    def test_sof_util_function_exception(self):

        from soxspipe.commonutils import sof_util
        try:
            this = sof_util(
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
