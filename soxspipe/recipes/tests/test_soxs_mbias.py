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


class test_mbias(unittest.TestCase):

    def test_mbias_from_directory_function(self):
        directory = settings["test-data-root"] + \
            "/xshooter-mbias/uvb/1x1/fast_read"

        from soxspipe.recipes import soxs_mbias
        this = soxs_mbias(
            log=log,
            settings=settings,
            inputFrames=directory
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    def test_mbias_from_sof_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mbias/sof/bias_uvb_1x1.sof"

        # utKit.refresh_database() # reset database to database found in
        # soxspipe.recipes/test/input
        from soxspipe.recipes import soxs_mbias
        this = soxs_mbias(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    def test_mbias_from_list_of_fits_function(self):
        directory = settings["test-data-root"] + \
            "/xshooter-mbias/uvb/1x1/slow_read"
        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser
        home = expanduser("~")
        if directory[0] == "~":
            directory = directory.replace("~", home)

        fileList = []
        for d in os.listdir(directory):
            filename = os.path.join(directory, d)
            if os.path.isfile(filename) and ".fits" in d:
                fileList.append(filename)

        # utKit.refresh_database() # reset database to database found in
        # soxspipe.recipes/test/input
        from soxspipe.recipes import soxs_mbias
        this = soxs_mbias(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    def test_produce_product_function(self):
        directory = settings["test-data-root"] + \
            "/xshooter-mbias/uvb/1x1/fast_read"
        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser
        home = expanduser("~")
        if directory[0] == "~":
            directory = directory.replace("~", home)

        fileList = []
        for d in os.listdir(directory):
            filename = os.path.join(directory, d)
            if os.path.isfile(filename) and ".fits" in d:
                fileList.append(filename)

        from soxspipe.recipes import soxs_mbias
        this = soxs_mbias(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    def test_mbias_mixed_image_type_exception(self):

        directory = settings["test-data-root"] + "/xshooter-lingain/vis"
        try:
            from soxspipe.recipes import soxs_mbias
            this = soxs_mbias(
                log=log,
                settings=settings,
                inputFrames=directory
            )
            this.get()
            assert False
        except Exception as e:
            assert True
            print(str(e))

    def test_mbias_wrong_image_type_exception(self):

        directory = settings["test-data-root"] + "/xshooter-mdark/vis"
        try:
            from soxspipe.recipes import soxs_mbias
            this = soxs_mbias(
                log=log,
                settings=settings,
                inputFrames=directory
            )
            this.get()
            assert False
        except Exception as e:
            assert True
            print(str(e))

    def test_mbias_mixed_binning_exception(self):

        directory = settings["test-data-root"] + "/xshooter-mbias/vis"
        try:
            from soxspipe.recipes import soxs_mbias
            this = soxs_mbias(
                log=log,
                settings=settings,
                inputFrames=directory
            )
            this.get()
            assert False
        except Exception as e:
            assert True
            print(str(e))

    def test_mbias_function_exception(self):

        from soxspipe.recipes import soxs_mbias
        try:
            this = soxs_mbias(
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
