from __future__ import print_function
import unittest
import pandas as pd
import numpy as np
from unittest.mock import patch
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from soxspipe.utKit import utKit
from fundamentals import tools

packageDirectory = utKit("").get_project_root()
settingsFile = packageDirectory + "/test_settings_xsh.yaml"
su = tools(
    arguments={"settingsFile": settingsFile},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName=None,
    defaultSettingsFile=False,
)
arguments, settings, log, dbConn = su.setup()


class test_toolkit(unittest.TestCase):

    def test_spectroscopic_image_quality_checks_defaults_missing_binning_headers_for_non_nir(self):
        self._run_spectroscopic_quality_check_with_header_binning(expected_binx=1, expected_biny=1)

    def test_spectroscopic_image_quality_checks_defaults_when_only_one_binning_header_exists(self):
        self._run_spectroscopic_quality_check_with_header_binning(winx=2, expected_binx=2, expected_biny=1)

    def _run_spectroscopic_quality_check_with_header_binning(self, expected_binx, expected_biny, winx=None, winy=None):
        from soxspipe.commonutils import keyword_lookup
        from soxspipe.commonutils import toolkit

        kw = keyword_lookup(log=log, settings=settings).get
        header = fits.Header()
        header[kw("SEQ_ARM")] = "UVB"
        header[kw("DATE_OBS")] = "2024-01-01T00:00:00"
        header[kw("INSTRUME")] = "XSHOOTER"
        if winx is not None:
            header[kw("WIN_BINX")] = winx
        if winy is not None:
            header[kw("WIN_BINY")] = winy

        frame = CCDData(
            np.ones((4, 4), dtype=float),
            unit=u.electron,
            header=header,
            mask=np.zeros((4, 4), dtype=bool),
        )

        observed = {}

        def fake_unpack_order_table(log, orderTablePath, binx, biny, prebinned):
            observed["binx"] = binx
            observed["biny"] = biny
            observed["prebinned"] = prebinned
            orderTablePixels = pd.DataFrame(
                {
                    "xcoord_edgeup": [2],
                    "xcoord_edgelow": [0],
                    "ycoord": [0],
                    "ycoord_edgeup": [2],
                    "ycoord_edgelow": [0],
                    "xcoord": [0],
                }
            )
            return pd.DataFrame(), orderTablePixels, pd.DataFrame()

        with patch.object(toolkit, "unpack_order_table", side_effect=fake_unpack_order_table), patch.object(
            toolkit, "quicklook_image", return_value=None
        ):
            qcTable = toolkit.spectroscopic_image_quality_checks(
                log=log,
                frame=frame,
                orderTablePath="dummy-order-table.fits",
                settings=settings,
                recipeName="soxs-mflat",
                qcTable=pd.DataFrame(),
            )

        assert observed == {"binx": expected_binx, "biny": expected_biny, "prebinned": True}
        assert {"INNER ORDER PIX MEAN", "INNER ORDER PIX SUM"}.issubset(set(qcTable["qc_name"]))
