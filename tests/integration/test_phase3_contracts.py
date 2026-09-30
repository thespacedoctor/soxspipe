"""Synthetic integration tests for Phase 3 FITS products."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.phase3 import write_fits_table_to_disk
from tests.factories import instrument_header, pipeline_settings

pytestmark = pytest.mark.integration


def test_phase3_table_writer_scrubs_header_and_writes_qc(
    tmp_path: Path,
    log: object,
) -> None:
    header = instrument_header()
    header["ARCFILE"] = "raw.fits"
    header["RADECSYS"] = "ICRS"
    header["HIERARCH ESO TEL TARG EQUINOX"] = 2000.0
    qc = pd.DataFrame(
        [{"qc_name": "RON", "qc_value": 3.2, "qc_comment": "Read noise", "to_header": True}]
    )
    outputPath = tmp_path / "phase3.fits"

    write_fits_table_to_disk(
        log=log,
        settings=pipeline_settings(tmp_path),
        header=header,
        tables=[Table({"order": [10, 11]}), pd.DataFrame({"flux": [1.0, 2.0]})],
        filePath=outputPath,
        qc=qc,
    )

    with fits.open(outputPath, checksum=True) as hdus:
        assert len(hdus) == 3
        assert hdus[0].header["ESO QC RON"] == 3.2
        assert "ESO DPR TYPE" not in hdus[0].header
        assert "ARCFILE" not in hdus[0].header
        assert hdus[0].header["RADESYS"] == "ICRS"
        assert hdus[0].header["EQUINOX"] == 2000.0
        assert list(hdus[1].data["order"]) == [10, 11]
        assert list(hdus[2].data["flux"]) == [1.0, 2.0]


def test_phase3_table_writer_omits_disabled_and_nan_qc_rows(
    tmp_path: Path,
    log: object,
) -> None:
    qc = pd.DataFrame(
        [
            {
                "qc_name": "DISABLED",
                "qc_value": 1.0,
                "qc_comment": "Not for the header",
                "to_header": False,
            },
            {
                "qc_name": "MISSING",
                "qc_value": np.nan,
                "qc_comment": "Missing value",
                "to_header": True,
            },
        ]
    )
    outputPath = tmp_path / "phase3.fits"

    write_fits_table_to_disk(
        log=log,
        settings=pipeline_settings(tmp_path),
        header=instrument_header(),
        tables=[Table({"order": [10]})],
        filePath=outputPath,
        qc=qc,
    )

    with fits.open(outputPath) as hdus:
        assert "ESO QC DISABLED" not in hdus[0].header
        assert "ESO QC MISSING" not in hdus[0].header
