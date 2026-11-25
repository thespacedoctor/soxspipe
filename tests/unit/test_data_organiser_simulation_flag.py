"""Unit contracts for the `simulation` flag of the data organiser (DY-2394)."""

from __future__ import annotations

import pandas as pd
import pytest

from tests.factories import workspace_organiser

pytestmark = pytest.mark.unit

MISSING = -99.99


def _frames(*rows: dict) -> pd.DataFrame:
    """Return frames with the columns the simulation flag reads, filled as the header harvest fills them."""
    defaults = {
        "eso seq arm": "NIR",
        "eso dpr type": "DARK",
        "eso ins nise swsim": MISSING,
        "eso ins vise swsim": MISSING,
        "junk": MISSING,
        "eso det sofw mode": "NORMAL",
        "simulation": "--",
    }
    return pd.DataFrame([{**defaults, **row} for row in rows])


def _soxs_organiser(tmp_path, log):
    """Return a SOXS organiser whose keyword lookup is set up."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser._select_instrument(inst="SOXS")
    return organiser


def _simulation(tmp_path, log, **row) -> int:
    """Return the simulation flag the organiser gives one frame."""
    organiser = _soxs_organiser(tmp_path, log)
    return int(organiser._add_simulation_flag(_frames(row)).loc[0, "simulation"])


@pytest.mark.parametrize("mode", ["HW-SIM", "LCU-SIM"])
def test_nir_dark_is_simulated_when_the_detector_software_mode_is_not_normal(tmp_path, log, mode) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso det sofw mode": mode})

    # ASSERT
    assert flag == 1


def test_nir_dark_is_not_simulated_when_only_the_slit_wheel_simulates(tmp_path, log) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso ins nise swsim": "T", "eso det sofw mode": "NORMAL"})

    # ASSERT
    assert flag == 0


def test_nir_bias_is_not_simulated_when_only_the_slit_wheel_simulates(tmp_path, log) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso dpr type": "BIAS", "eso ins nise swsim": "T"})

    # ASSERT
    assert flag == 0


@pytest.mark.parametrize("dprType", ["LAMP,FLAT", "LAMP,WAVE", "LAMP,ORDERDEF", "OBJECT", "STD,FLUX"])
def test_nir_frame_that_needs_the_slit_is_simulated_when_the_slit_wheel_simulates(tmp_path, log, dprType) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso dpr type": dprType, "eso ins nise swsim": "T"})

    # ASSERT
    assert flag == 1


def test_nir_frame_that_needs_the_slit_is_simulated_when_the_detector_simulates(tmp_path, log) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso dpr type": "LAMP,WAVE", "eso det sofw mode": "HW-SIM"})

    # ASSERT
    assert flag == 1


@pytest.mark.parametrize("mode", [MISSING, "--"])
def test_nir_frame_without_a_detector_software_mode_is_not_simulated(tmp_path, log, mode) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso det sofw mode": mode})

    # ASSERT
    assert flag == 0


def test_nir_frame_with_a_normal_detector_and_no_slit_wheel_simulation_is_not_simulated(tmp_path, log) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso dpr type": "LAMP,FLAT"})

    # ASSERT
    assert flag == 0


def test_vis_dark_keeps_its_slit_wheel_simulation_flag(tmp_path, log) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso seq arm": "VIS", "eso ins vise swsim": "T"})

    # ASSERT
    assert flag == 1


def test_vis_dark_is_not_simulated_by_the_nir_detector_software_mode(tmp_path, log) -> None:
    # ARRANGE, ACT
    flag = _simulation(tmp_path, log, **{"eso seq arm": "VIS", "eso det sofw mode": "HW-SIM"})

    # ASSERT
    assert flag == 0


def test_uvb_frame_is_never_simulated(tmp_path, log) -> None:
    # ARRANGE, ACT
    flag = _simulation(
        tmp_path, log, **{"eso seq arm": "UVB", "eso det sofw mode": "HW-SIM", "eso ins nise swsim": "T"}
    )

    # ASSERT
    assert flag == 0


def test_every_frame_of_a_mixed_table_gets_its_own_flag(tmp_path, log) -> None:
    # ARRANGE
    organiser = _soxs_organiser(tmp_path, log)
    frames = _frames(
        {"eso det sofw mode": "HW-SIM"},
        {"eso ins nise swsim": "T"},
        {"eso seq arm": "VIS", "eso ins vise swsim": "T"},
        {"eso dpr type": "OBJECT", "eso ins nise swsim": "T"},
    )

    # ACT
    flags = organiser._add_simulation_flag(frames)["simulation"].astype(int).tolist()

    # ASSERT
    assert flags == [1, 0, 1, 1]
