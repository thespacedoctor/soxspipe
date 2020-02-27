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
# settingsFile = home +
# "/.config/soxspipe.commonutils/soxspipe.commonutils.yaml"
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


class test_keyword_lookup(unittest.TestCase):

    def test_keyword_dictionary_selection_function(self):

        from soxspipe.commonutils import keyword_lookup
        this = keyword_lookup(
            log=log,
            settings=settings
        )
        kw = this._select_dictionary()
        print(kw["DET_NDITSKIP"])
        print(kw["MJDOBS"])
        print(kw["NAXIS1"])

    def test_keyword_lookup_function(self):

        from soxspipe.commonutils import keyword_lookup
        kw = keyword_lookup(
            log=log,
            settings=settings
        ).get
        if settings["instrument"] == "xsh":
            assert kw("DET_NDITSKIP") == "ESO DET NDITSKIP"
            assert kw("PROV", 14) == "PROV14"
            assert kw("PROV", 1) == "PROV01"

        # AND HOW ABOUT LISTS
        if settings["instrument"] == "xsh":
            assert kw(["PROV", "DET_NDITSKIP"])[1] == 'ESO DET NDITSKIP'

    def test_keyword_lookup_function_exception(self):

        from soxspipe.commonutils import keyword_lookup
        try:
            this = keyword_lookup(
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
