from __future__ import print_function
from astropy import units as u
from astropy.nddata import CCDData
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
#     "/git_repos/_misc_/settings/soxspipe/test_settings_xsh.yaml"

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


class test_detect_order_edges(unittest.TestCase):

    import pytest

    @pytest.mark.full
    def test_xsh_detect_vis_order_edges_function(self):

        flatPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-detect-order-edges/first_iteration_VIS_master_flat.fits"
        orderCentreTable = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-detect-order-edges/20170818T172920_VIS_1X1_FAST_ORDER_LOCATIONS.fits"
        home = expanduser("~")
        flatPath = flatPath.replace("~", home)
        orderCentreTable = orderCentreTable.replace("~", home)
        flatFrame = CCDData.read(flatPath, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                 hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        import pandas as pd
        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
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

        from soxspipe.commonutils import detect_order_edges
        edges = detect_order_edges(
            log=log,
            flatFrame=flatFrame,
            orderCentreTable=orderCentreTable,
            settings=settings,
            recipeSettings=settings["soxs-mflat"],
            qcTable=qc,
            productsTable=products,
        )
        edges.get()

    @pytest.mark.full
    def test_soxs_detect_order_edges_function_exception(self):

        from soxspipe.commonutils import detect_order_edges
        try:
            this = detect_order_edges(
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
