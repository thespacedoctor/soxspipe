"""A recipe run end to end from a list of frames and from a directory, with no `.sof` file (DY-90).

`soxs_mbias` runs for real over synthetic raw bias frames written by the
factories. Only the workspace session lookup is replaced, so the run leaves no
trace outside `tmp_path`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml
from astropy.io import fits
from astropy.time import Time

import soxspipe.commonutils as commonutils
from soxspipe.recipes import soxs_mbias
from tests.factories import raw_fits

pytestmark = pytest.mark.integration

# RAW VIS FRAME SIZE: THE SCIENCE REGION OF THE DETECTOR PLUS ITS PRE-SCAN COLUMNS
RAW_SHAPE = (4096, 858)
FRAME_TIMES = ("2024-01-02T20:05:00", "2024-01-02T20:00:00")
EXPECTED_NIGHT = "2024-01-02"
EXPECTED_SOF_NAME = f"soxs-mbias_VIS_{EXPECTED_NIGHT}.sof"


class _IsolatedOrganiser:
    """Session lookup boundary with no workspace session."""

    def __init__(self, **_: object) -> None:
        return None

    def session_list(self, *, silent: bool) -> tuple[bool, list[str]]:
        return False, []

    def close(self) -> None:
        return None


def _settings(workspace: Path) -> dict[str, Any]:
    """Return the packaged advanced and SOXS default settings rooted in a temporary workspace."""
    import soxspipe

    packageDirectory = Path(soxspipe.__file__).parent
    settings: dict[str, Any] = {}
    for fileName in ("advanced_settings.yaml", "soxs_default_settings.yaml"):
        settings.update(yaml.safe_load((packageDirectory / fileName).read_text(encoding="utf-8")))
    return {**settings, "workspace-root-dir": str(workspace)}


def _bias_frames(rawDirectory: Path) -> list[str]:
    """Write two raw VIS bias frames, and return their paths."""
    rawDirectory.mkdir(parents=True, exist_ok=True)
    paths = []
    for seed, isoTime in enumerate(FRAME_TIMES, start=1):
        path = raw_fits(
            rawDirectory / f"bias_{seed}.fits",
            shape=RAW_SHAPE,
            seed=seed,
            headerOverrides={
                "MJD-OBS": float(Time(isoTime, scale="utc").mjd),
                "SEQ_ARM": "VIS",
                "DATE_OBS": isoTime,
                "ESO DET BINX": 1,
                "ESO DET BINY": 1,
                "ESO INS TEMP104 VAL": 150.0,
                "ESO INS TEMP301 VAL": 10.0,
            },
        )
        paths.append(str(path))
    return paths


def test_master_bias_runs_to_completion_from_a_frame_list_without_a_sof_file(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A frame list produces a master bias, QC rows and products under the start of night."""
    # ARRANGE
    monkeypatch.setattr(commonutils, "data_organiser", _IsolatedOrganiser)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    framePaths = _bias_frames(tmp_path / "raw")

    # ACT
    productPath, qcTable = soxs_mbias(
        log=log,
        settings=_settings(workspace),
        inputFrames=framePaths,
        turnOffMP=True,
    ).produce_product()

    # ASSERT
    expectedDirectory = workspace / "reduced" / EXPECTED_NIGHT / "soxs-mbias"
    assert Path(productPath).parent == expectedDirectory
    assert Path(productPath).exists()
    with fits.open(productPath) as hdus:
        assert hdus[0].data.shape == (4096, 808)
    assert set(qcTable["sof_name"]) == {EXPECTED_SOF_NAME}
    assert (workspace / "qc" / EXPECTED_NIGHT / "soxs-mbias").is_dir()
