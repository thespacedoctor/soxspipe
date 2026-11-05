"""Contracts for a `--vlt` workspace, whose raw and reduced folders sit outside the workspace (DY-2113)."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest
from astropy.io import fits

from soxspipe.commonutils import data_organiser
from soxspipe.commonutils.data_organiser import _UnsafePathError
from tests.factories import harvestable_raw_fits

pytestmark = pytest.mark.integration

NIGHT = "2024-01-02"
BIAS_SET_SIZE = 3
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
# RUNS THE REAL CLI WITH THE VLT FOLDERS OF THE TEST, BECAUSE THE PACKAGED SETTINGS NAME `/data/raw` AND `/data/reduced`
CLI_DRIVER = """
import importlib
import sys

organiserModule = importlib.import_module("soxspipe.commonutils.data_organiser")
vltFolders = (sys.argv[1], sys.argv[2])
organiserModule._vlt_data_folders = lambda: vltFolders
sys.argv = ["soxspipe"] + sys.argv[3:]

from soxspipe.cl_utils import main

main()
"""


@dataclass(frozen=True)
class VltMachine:
    """The folders of a VLT workstation: `vltRaw` is itself a symlink, as `/data/raw` is on `wsxpl`."""

    workspace: Path
    vltRaw: Path
    vltReduced: Path
    diskRaw: Path
    diskReduced: Path
    driver: Path
    home: Path


def _bias_frame(destination: Path, exposure: int) -> None:
    harvestable_raw_fits(destination)
    with fits.open(destination, mode="update") as hdus:
        header = hdus[0].header
        header["ESO TPL NEXP"] = BIAS_SET_SIZE
        header["ESO TPL EXPNO"] = exposure + 1
        header["ESO TPL START"] = f"{NIGHT}T03:00:00"
        header["DATE-OBS"] = f"{NIGHT}T03:04:0{exposure}.000"
        header["MJD-OBS"] = 60311.75 + exposure * 1e-5
        # THE PACKAGE SHIPS THE 1X1 BAD-PIXEL MAP; ANY OTHER BINNING WOULD WRITE A NEW MAP INTO THE PACKAGE DIRECTORY
        header["ESO DET BINY"] = 1
        header["ESO INS TEMP301 VAL"] = 10.0
        header["ESO INS TEMP104 VAL"] = -100.0


@pytest.fixture(name="vltMachine")
def vlt_machine(tmp_path: Path) -> VltMachine:
    """Build a workstation whose raw night folder holds one complete bias set."""
    diskData = tmp_path / "diskb" / "data"
    diskRaw = diskData / "raw"
    diskReduced = diskData / "reduced"
    (diskRaw / NIGHT).mkdir(parents=True)
    diskReduced.mkdir()
    dataPath = tmp_path / "data"
    dataPath.mkdir()
    vltRaw = dataPath / "raw"
    vltReduced = dataPath / "reduced"
    vltRaw.symlink_to(diskRaw, target_is_directory=True)
    vltReduced.symlink_to(diskReduced, target_is_directory=True)
    workspace = tmp_path / "diska" / "home" / "pipeline" / "soxspipe_reductions"
    workspace.mkdir(parents=True)
    home = tmp_path / "user-home"
    (home / ".config" / "soxspipe").mkdir(parents=True)
    shutil.copy(REPOSITORY_ROOT / "soxspipe" / "soxs_default_settings.yaml", home / ".config/soxspipe/soxspipe.yaml")
    driver = tmp_path / "cli_driver.py"
    driver.write_text(CLI_DRIVER, encoding="utf-8")
    for exposure in range(BIAS_SET_SIZE):
        _bias_frame(diskRaw / NIGHT / f"bias{exposure}.fits", exposure)
    return VltMachine(workspace, vltRaw, vltReduced, diskRaw, diskReduced, driver, home)


def _soxspipe(machine: VltMachine, *arguments: str) -> subprocess.CompletedProcess[str]:
    """Run one `soxspipe` command in the workspace, in its own process, as a user would."""
    environment = {**os.environ, "PYTHONPATH": str(REPOSITORY_ROOT), "HOME": str(machine.home)}
    return subprocess.run(  # noqa: S603
        [sys.executable, str(machine.driver), str(machine.vltRaw), str(machine.vltReduced), *arguments],
        cwd=machine.workspace,
        env=environment,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )


def _assert_passed(result: subprocess.CompletedProcess[str]) -> None:
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"


def test_every_command_works_after_the_first_prep_of_a_vlt_workspace(vltMachine: VltMachine) -> None:
    """Run the command sequence that failed on `wsxpl`, each in a new process."""
    firstPrep = _soxspipe(vltMachine, "prep", ".", "--vlt")
    _assert_passed(firstPrep)
    assert "PREPARING" in firstPrep.stdout

    sessionList = _soxspipe(vltMachine, "session", "ls")
    _assert_passed(sessionList)
    assert "base" in sessionList.stdout

    _assert_passed(_soxspipe(vltMachine, "prep", "."))
    _assert_passed(_soxspipe(vltMachine, "prep", ".", "--vlt"))
    _assert_passed(_soxspipe(vltMachine, "prep", ".", "--refresh"))
    reduction = _soxspipe(vltMachine, "reduce", "all", ".")
    _assert_passed(reduction)

    assert (vltMachine.workspace / "raw").is_symlink()
    assert sorted(path.name for path in (vltMachine.diskRaw / NIGHT).iterdir()) == [
        f"bias{exposure}.fits" for exposure in range(BIAS_SET_SIZE)
    ]
    # THE SYNTHETIC FRAMES ARE TOO SMALL FOR THE DETECTOR MASK, SO THE RECIPE STOPS AFTER IT READS THEM THROUGH THE LINK
    assert list(vltMachine.diskReduced.rglob("soxs-mbias/*.log")), reduction.stdout

    sofName = next((vltMachine.workspace / "sof").iterdir()).name
    _assert_passed(_soxspipe(vltMachine, "raw", "sof", sofName))
    assert (vltMachine.workspace / "exported" / "bias0.fits").is_file()


@pytest.fixture(name="vltOrganiserModule")
def vlt_organiser_module(vltMachine: VltMachine, monkeypatch: pytest.MonkeyPatch) -> object:
    """Point the organiser module at the VLT folders of the test, from inside the workspace."""
    organiserModule = sys.modules[data_organiser.__module__]
    monkeypatch.setattr(
        organiserModule, "_vlt_data_folders", lambda: (str(vltMachine.vltRaw), str(vltMachine.vltReduced))
    )
    monkeypatch.chdir(vltMachine.workspace)
    return organiserModule


@pytest.mark.usefixtures("vltOrganiserModule")
def test_vlt_raw_link_is_accepted_without_the_vlt_flag(vltMachine: VltMachine, log: object) -> None:
    """A later command has no `--vlt` flag, so the constructor must accept the link on its own."""
    (vltMachine.workspace / "raw").symlink_to(vltMachine.vltRaw, target_is_directory=True)

    organiser = data_organiser(log=log, rootDir=".", dbConnect=False)

    assert organiser.rawDir == str(vltMachine.diskRaw.resolve())


@pytest.mark.usefixtures("vltOrganiserModule")
def test_vlt_flag_replaces_the_link_and_leaves_the_vlt_folder_in_place(vltMachine: VltMachine, log: object) -> None:
    """A repeated `--vlt` acts on the workspace link, never on the folder the link leads to."""
    (vltMachine.workspace / "raw").symlink_to(vltMachine.vltRaw, target_is_directory=True)

    data_organiser(log=log, rootDir=".", vlt=True, dbConnect=False)

    assert (vltMachine.workspace / "raw").is_symlink()
    assert vltMachine.diskRaw.is_dir()


@pytest.mark.usefixtures("vltOrganiserModule")
def test_raw_link_to_a_folder_below_the_vlt_raw_folder_is_refused(vltMachine: VltMachine, log: object) -> None:
    """Only the configured folder itself is allowed, not everything below it."""
    (vltMachine.workspace / "raw").symlink_to(vltMachine.diskRaw / NIGHT, target_is_directory=True)

    with pytest.raises(_UnsafePathError, match="raw directory"):
        data_organiser(log=log, rootDir=".", dbConnect=False)


@pytest.mark.usefixtures("vltOrganiserModule")
def test_raw_link_to_another_outside_folder_is_refused(vltMachine: VltMachine, log: object, tmp_path: Path) -> None:
    """Configured VLT folders do not widen the check to other locations."""
    outsidePath = tmp_path / "outside"
    outsidePath.mkdir()
    (vltMachine.workspace / "raw").symlink_to(outsidePath, target_is_directory=True)

    with pytest.raises(_UnsafePathError, match="raw directory"):
        data_organiser(log=log, rootDir=".", dbConnect=False)


@pytest.mark.usefixtures("vltOrganiserModule")
@pytest.mark.parametrize("linkTarget", ["vltRaw", "vltReduced"])
def test_sessions_link_to_a_vlt_folder_is_refused(vltMachine: VltMachine, log: object, linkTarget: str) -> None:
    """The VLT folders are allowed for `raw` only; `sessions` must stay in the workspace."""
    (vltMachine.workspace / "sessions").symlink_to(getattr(vltMachine, linkTarget), target_is_directory=True)

    with pytest.raises(_UnsafePathError, match="sessions directory"):
        data_organiser(log=log, rootDir=".", dbConnect=False)


@pytest.mark.usefixtures("vltOrganiserModule")
def test_raw_link_loop_is_refused(vltMachine: VltMachine, log: object) -> None:
    """A link loop is an unsafe path, not a crash."""
    rawLink = vltMachine.workspace / "raw"
    rawLink.symlink_to(rawLink)

    with pytest.raises(_UnsafePathError, match="raw directory"):
        data_organiser(log=log, rootDir=".", dbConnect=False)


@pytest.mark.usefixtures("vltOrganiserModule")
def test_database_filepath_outside_workspace_and_vlt_raw_is_refused(
    vltMachine: VltMachine, log: object, tmp_path: Path
) -> None:
    """Allowing the VLT raw folder does not allow other paths held in the database."""
    (vltMachine.workspace / "raw").symlink_to(vltMachine.vltRaw, target_is_directory=True)
    organiser = data_organiser(log=log, rootDir=".", dbConnect=False)

    insidePath = organiser._validate_database_filepath(Path("raw") / NIGHT / "bias0.fits")

    assert insidePath == (vltMachine.diskRaw / NIGHT / "bias0.fits").resolve()
    with pytest.raises(_UnsafePathError, match="database filepath"):
        organiser._validate_database_filepath(tmp_path / "outside.fits")


def test_prepare_indexes_a_new_vlt_frame_with_an_absolute_workspace_path(
    vltMachine: VltMachine, vltOrganiserModule: object, log: object
) -> None:
    """The `watch` daemon passes an absolute workspace path, which is not part of a path in the VLT raw folder."""
    data_organiser(log=log, rootDir=str(vltMachine.workspace), vlt=True).prepare(report=False)
    _bias_frame(vltMachine.diskRaw / NIGHT / "late.fits", 0)

    organiser = data_organiser(log=log, rootDir=str(vltMachine.workspace))
    organiser.prepare(report=False)

    indexed = {row[0] for row in organiser.conn.execute("SELECT file FROM raw_frames")}
    assert "late.fits" in indexed
    assert (vltMachine.diskRaw / NIGHT / "late.fits").is_file()
    assert not (vltMachine.workspace / "late.fits").exists()


@pytest.mark.usefixtures("vltOrganiserModule")
def test_superseded_frame_in_the_vlt_raw_folder_is_not_deleted(vltMachine: VltMachine, log: object) -> None:
    """The VLT raw folder is not owned by the workspace, so `prep` deletes nothing there."""
    (vltMachine.workspace / "raw").symlink_to(vltMachine.vltRaw, target_is_directory=True)
    twinPath = vltMachine.diskRaw / NIGHT / "twin.fits"
    twinPath.write_bytes(b"plain")
    (vltMachine.diskRaw / NIGHT / "twin.fits.Z").write_bytes(b"compressed")
    organiser = data_organiser(log=log, rootDir=".", dbConnect=False)

    organiser._delete_superseded_frames()

    assert twinPath.exists()
