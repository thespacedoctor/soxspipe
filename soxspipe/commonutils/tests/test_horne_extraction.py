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

class test_horne_extraction(unittest.TestCase):

    import pytest

    @pytest.mark.full
    def test_horne_extraction_function(self):

        # XRAY BINARY (FAINT)
        skyModelFrame = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/2019.08.31T00.13.27.1305_NIR_STARE_300PT0_SAX_J1808.43658_SKYMODEL.fits"
        skySubtractedFrame = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/2019.08.31T00.13.27.1305_NIR_STARE_300PT0_SAX_J1808.43658_SKYSUB.fits"
        twoDMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/2019.08.30T18.43.48.7597_NIR_SPAT_SOL_0PT6651_IMAGE.fits"
        dispMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/2019.08.30T18.43.48.7597_NIR_SPAT_SOL_0PT6651.fits"

        # STANDARD STAR (BRIGHT)
        skyModelFrame = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/eg274/2019.08.22T23.12.18.5011_NIR_STARE_205PT0_EG_274_SKYMODEL.fits"
        skySubtractedFrame = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/eg274/2019.08.22T23.12.18.5011_NIR_STARE_205PT0_EG_274_SKYSUB.fits"
        twoDMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/eg274/2019.08.23T15.54.49.8571_NIR_SPAT_SOL_0PT6651_IMAGE.fits"
        dispMap = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-horne-extraction/nir/eg274/2019.08.23T15.54.49.8571_NIR_SPAT_SOL_0PT6651.fits"

        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import figure

        qc = pd.DataFrame({
            "soxspipe_recipe": [],
            "qc_name": [],
            "qc_value": [],
            "qc_unit": [],
            "qc_comment": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "to_header": []
        })
        products = pd.DataFrame({
            "soxspipe_recipe": [],
            "product_label": [],
            "file_name": [],
            "file_type": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "file_path": [],
            "label": []
        })

        from os.path import expanduser
        home = expanduser("~")
        twoDMap = twoDMap.replace("~", home)
        skyModelFrame = skyModelFrame.replace("~", home)
        skySubtractedFrame = skySubtractedFrame.replace("~", home)
        dispMap = dispMap.replace("~", home)

        from soxspipe.commonutils import horne_extraction
        optimalExtractor = horne_extraction(
            log=log,
            skyModelFrame=skyModelFrame,
            skySubtractedFrame=skySubtractedFrame,
            twoDMapPath=twoDMap,
            settings=settings,
            recipeName="soxs-stare",
            qcTable=qc,
            productsTable=products,
            dispersionMap=dispMap,
            sofName="2019.08.22T23.12.18.5011_NIR_STARE_205PT0_EG_274"
        )
        #(w, f) = optimalExtractor.extract(15,3, 10, 2.0)

        #figure(figsize=(13, 6))

        #plt.ylim(-100, 300)
        #plt.plot(w, f)

        # plt.show()
        e_w_l = []
        e_s = []
        en_s = []

        for order in range(11, 12):
            extracted_wave_spectrum, extracted_spectrum, nonopt = optimalExtractor.extract(order)
            if len(extracted_wave_spectrum):
                e_w_l.append(extracted_wave_spectrum)
                e_s.append(extracted_spectrum)
                en_s.append(nonopt)
                #plt.plot(extracted_wave_spectrum, extracted_spectrum)

        plt.clf()
        figure(figsize=(13, 6))

        breakpoint()

        for (w, s) in zip(e_w_l, e_s):
            if len(w) and len(s):
                try:
                    plt.plot(w, s)
                except:
                    breakpoint()
        plt.show()

    @pytest.mark.full
    def test_horne_extraction_function_exception(self):

        from soxspipe.commonutils import horne_extraction
        try:
            this = horne_extraction(
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
