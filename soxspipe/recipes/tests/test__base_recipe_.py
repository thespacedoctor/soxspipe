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
# settingsFile = home + "/.config/soxspipe.recipes/soxspipe.recipes.yaml"
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


class test__base_recipe_(unittest.TestCase):

    def test__base_recipe__function(self):

        framePath = settings["test-data-root"] + \
            "/xshooter-mbias/uvb/XSHOO.2019-07-03T10:40:24.434.fits"
        interMediatePath = settings["intermediate-data-root"]
        from soxspipe.recipes import _base_recipe_
        recipe = _base_recipe_(
            log=log,
            settings=settings
        )

        from soxspipe.commonutils import detector_lookup
        recipe.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get("UVB")
        recipe.detectorParams["gain"] = 1.75
        recipe.detectorParams["ron"] = 4.5
        recipe.arm = "UVB"

        preFrame = recipe._prepare_single_frame(frame=framePath)

        # NOW TRY SAVING
        preFrame = recipe._prepare_single_frame(frame=framePath, save=settings[
            "save-intermediate-products"])

        # NOW CLEAN UP
        recipe.clean_up()

    def test__base_recipe__function_exception(self):

        from soxspipe.recipes import _base_recipe_
        try:
            this = _base_recipe_(
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
