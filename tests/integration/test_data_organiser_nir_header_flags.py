"""Integration contract: `prepare` turns NIR header cards into the simulation and ignore flags (DY-2394)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from astropy.io import fits

from tests.factories import harvestable_raw_fits, workspace_organiser

pytestmark = pytest.mark.integration


def _nir_dark(path: Path, *, mode: str, temperature: float, slitWheelSimulates: bool, slit: str) -> None:
    """Write a NIR dark whose header holds the cards the DY-2394 rules read."""
    harvestable_raw_fits(path)
    with fits.open(path, mode="update") as hdus:
        header = hdus[0].header
        header["ESO SEQ ARM"] = "NIR"
        header["ESO DET3 CAM NAME"] = "NIR"
        header["ESO DPR TYPE"] = "DARK"
        header["ESO DPR TECH"] = "IMAGE"
        header["ESO DET SOFW MODE"] = mode
        header["ESO INS TEMP217 VAL"] = temperature
        header["ESO INS NISE SWSIM"] = slitWheelSimulates
        header["ESO INS NISE NAME"] = slit


def test_prepare_flags_simulated_and_warm_nir_darks_and_keeps_a_good_dark_behind_an_undefined_slit(
    tmp_path, log, monkeypatch
) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    rootPath = Path(organiser.rootDir)
    _nir_dark(rootPath / "good.fits", mode="NORMAL", temperature=45.0, slitWheelSimulates=True, slit="UNDEFINED")
    _nir_dark(rootPath / "hwsim.fits", mode="HW-SIM", temperature=48.5, slitWheelSimulates=False, slit="BLANK")
    _nir_dark(rootPath / "warm.fits", mode="NORMAL", temperature=250.0, slitWheelSimulates=False, slit="BLANK")
    monkeypatch.chdir(rootPath)

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    frames = pd.read_sql('SELECT file, simulation, ignore, "nir temp k" FROM raw_frames ORDER BY file', organiser.conn)
    flags = {row.file: (int(row.simulation), int(row.ignore)) for row in frames.itertuples()}
    assert flags == {"good.fits": (0, 0), "hwsim.fits": (1, 1), "warm.fits": (0, 1)}
