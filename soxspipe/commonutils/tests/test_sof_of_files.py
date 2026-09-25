from __future__ import print_function
from builtins import str
import os
import unittest
import shutil
import tempfile
import unittest
import yaml
from astropy.io import fits
from soxspipe.utKit import utKit
from fundamentals import tools
from os.path import expanduser

home = expanduser("~")

packageDirectory = utKit("", dbConn=False).get_project_root()
settingsFile = packageDirectory + "/test_settings_xsh.yaml"
# settingsFile = home + "/.config/soxspipe/soxspipe.yaml"
su = tools(
    arguments={"settingsFile": settingsFile},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName=None,
    defaultSettingsFile=False,
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


class test_set_of_files(unittest.TestCase):

    def test_sof_supplementary_files_are_excluded_from_image_collection(self):
        from soxspipe.commonutils.set_of_files import set_of_files

        with tempfile.TemporaryDirectory() as temp_directory:
            frame_directory = os.path.join(temp_directory, "frames")
            supplementary_directory = os.path.join(temp_directory, "supplementary")
            os.makedirs(frame_directory)
            os.makedirs(supplementary_directory)

            frame_path = os.path.join(frame_directory, "a.fits")
            supplementary_path = os.path.join(
                supplementary_directory, "VIS_ORDER_LOCATIONS.csv"
            )
            sof_path = os.path.join(temp_directory, "input.sof")
            fits.PrimaryHDU().writeto(frame_path)
            with open(supplementary_path, "w", encoding="utf-8") as stream:
                stream.write("order,centre\n")
            with open(sof_path, "w", encoding="utf-8") as stream:
                stream.write(f"{frame_path} BIAS_VIS\n{supplementary_path}\n")

            sof_reader = set_of_files.__new__(set_of_files)
            sof_reader.log = log
            sof_reader.inputFrames = sof_path
            sof_reader.currentSession = None
            sof_reader.ext = 0
            sof_reader.keys = ["file"]

            collection, supplementary = sof_reader.get()

            self.assertEqual(collection.location, frame_directory)
            self.assertEqual(collection.files, ["a.fits"])
            self.assertEqual(list(collection.summary["file"]), ["a.fits"])
            self.assertEqual(
                supplementary["VIS"]["ORDER_LOCATIONS"], supplementary_path
            )

    def test01_xsh_set_of_files_function(self):
        directory = settings["test-data-root"] + "/xshooter-mbias/vis"
        sofPath = "~/xshooter-pipeline-data/unittest_data/test.sof"
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(log=log, settings=settings)
        sofFile = sof._generate_sof_file_from_directory(
            directory=directory, sofPath=sofPath
        )
        print("sof file written to %(sofPath)s" % locals())

    def test_xsh_sof_to_collection_from_directory_function(self):
        directory = settings["test-data-root"] + "/xshooter-mbias/vis"
        sofPath = "~/xshooter-pipeline-data/unittest_data/test.sof"
        from soxspipe.commonutils.set_of_files import set_of_files

        print(directory)
        sof = set_of_files(log=log, settings=settings, inputFrames=directory)
        sofFile, supplementaryInput = sof.get()
        print(sofFile.summary)

    def test_sof_to_collection_from_sof_file_function(self):
        directory = settings["test-data-root"] + "/xshooter-mbias/vis"
        sofPath = "~/xshooter-pipeline-data/unittest_data/test.sof"
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(log=log, settings=settings, inputFrames=sofPath)
        sofFile, supplementaryInput = sof.get()
        print(sofFile.summary)

    def test_xsh_sof_to_collection_from_list_of_fits_files_function(self):
        directory = settings["test-data-root"] + "/xshooter-mbias/vis"
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

        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(log=log, settings=settings, inputFrames=fileList)
        sofFile, supplementaryInput = sof.get()
        print(sofFile.summary)

    def test_xsh_validate_sof_frames_function(self):
        directory = settings["test-data-root"] + "/xshooter-mbias/vis"
        sofPath = "~/xshooter-pipeline-data/unittest_data/test.sof"
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(log=log, settings=settings)
        sofFile = sof._generate_sof_file_from_directory(
            directory=directory, sofPath=sofPath
        )

        print("sof file written to %(sofPath)s" % locals())

    def test_soxs_set_of_files_function_exception(self):

        from soxspipe.commonutils.set_of_files import set_of_files

        try:
            this = set_of_files(log=log, settings=settings, fakeKey="break the code")
            this.get()
            assert False
        except Exception as e:
            assert True
            print(str(e))

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
