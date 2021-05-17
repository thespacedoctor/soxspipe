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
# xt-utkit-refresh-database

class test_dispersion_map_to_pixel_arrays(unittest.TestCase):

    def test_dispersion_map_to_pixel_arrays_function(self):

        dispersionMapPath = "~/xshooter-pipeline-data/unittest_data/detect_continuum/single_pinhole_NIR_disp_map.csv"
        from soxspipe.commonutils import dispersion_map_to_pixel_arrays
        import pandas as pd
        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {
            "order": [11, 11, 11, 11, 11],
            "wavelength": [2000., 2100., 2200., 2300., 2400.],
            "slit_position": [0, 0, 0, 0, 0]
        }
        orderPixelTable = pd.DataFrame(myDict)
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=log,
            dispersionMapPath=dispersionMapPath,
            orderPixelTable=orderPixelTable
        )

    def test_dispersion_map_to_pixel_arrays_function_exception(self):

        from soxspipe.commonutils import dispersion_map_to_pixel_arrays
        try:
            this = dispersion_map_to_pixel_arrays(
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
