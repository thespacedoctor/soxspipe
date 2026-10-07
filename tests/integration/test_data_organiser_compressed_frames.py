"""Integration contracts for preparing a workspace of LZW-compressed `.fits.Z` raw frames (DY-310)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.set_of_files import set_of_files
from tests.factories import (
    harvestable_raw_fits,
    lzw_compressed_fits,
    pipeline_settings,
    workspace_organiser,
)

pytestmark = pytest.mark.integration

LZW_MAGIC = b"\x1f\x9d"
NIGHT = "2024-01-02"
# THE raw_frames_valid VIEW ADMITS A BIAS SET ONLY WHEN IT HOLDS ALL `ESO TPL NEXP` EXPOSURES
BIAS_SET_SIZE = 3
MBIAS_SOF = "20240102T030400_VIS_1X2_2_MBIAS_SOXS.sof"


def _compressed_frame(tmp_path: Path, destination: Path, *, biasExposure: int | None = None) -> Path:
    """Write a harvestable bias frame as ``destination`` (a `.fits.Z` path), leaving no plain copy behind.

    With ``biasExposure`` set, the frame is that exposure of a complete bias template, so `prep` groups it.
    """
    plainPath = harvestable_raw_fits(tmp_path / f"plain-{destination.name}.fits")
    if biasExposure is not None:
        with fits.open(plainPath, mode="update") as hdus:
            header = hdus[0].header
            header["ESO TPL NEXP"] = BIAS_SET_SIZE
            header["ESO TPL EXPNO"] = biasExposure + 1
            header["ESO TPL START"] = f"{NIGHT}T03:00:00"
            header["DATE-OBS"] = f"{NIGHT}T03:04:0{biasExposure}.000"
            header["MJD-OBS"] = 60311.75 + biasExposure * 1e-5
    compressedPath = lzw_compressed_fits(plainPath, destination)
    plainPath.unlink()
    return compressedPath


def _sof_settings(tmp_path: Path) -> dict[str, object]:
    summaryKeys = ["DPR_TYPE", "SEQ_ARM"]
    return pipeline_settings(
        tmp_path,
        overrides={"summary-keys": {"default": summaryKeys, "verbose": summaryKeys, "nodding_extras": []}},
    )


def _indexed_frames(organiser) -> list[dict[str, str]]:
    return pd.read_sql("SELECT file, filepath FROM raw_frames ORDER BY file", organiser.conn).to_dict("records")


@pytest.fixture
def organiser(tmp_path, log, monkeypatch):
    organiser = workspace_organiser(tmp_path, log=log)
    monkeypatch.chdir(organiser.rootDir)
    yield organiser
    if getattr(organiser, "conn", None) is not None:
        organiser.conn.close()


def test_prepare_indexes_a_root_frame_and_moves_it_into_raw_still_compressed(tmp_path, organiser) -> None:
    rootPath = Path(organiser.rootDir)
    _compressed_frame(tmp_path, rootPath / "bias.fits.Z")

    organiser.prepare(report=False)

    assert _indexed_frames(organiser) == [{"file": "bias.fits.Z", "filepath": f"./raw/{NIGHT}/bias.fits.Z"}]
    movedPath = rootPath / "raw" / NIGHT / "bias.fits.Z"
    assert movedPath.read_bytes()[:2] == LZW_MAGIC
    assert list(rootPath.rglob("bias.fits")) == []


def test_prepare_indexes_a_frame_already_nested_in_the_raw_tree(tmp_path, organiser) -> None:
    rootPath = Path(organiser.rootDir)
    nestedPath = rootPath / "raw" / NIGHT
    nestedPath.mkdir(parents=True)
    _compressed_frame(tmp_path, nestedPath / "bias.fits.Z")

    organiser.prepare(report=False)

    assert _indexed_frames(organiser) == [{"file": "bias.fits.Z", "filepath": f"./raw/{NIGHT}/bias.fits.Z"}]
    assert (nestedPath / "bias.fits.Z").read_bytes()[:2] == LZW_MAGIC
    assert not (rootPath / "bias.fits.Z").exists()


def test_prepare_writes_compressed_members_into_a_sof_a_recipe_can_read(tmp_path, organiser, log) -> None:
    rootPath = Path(organiser.rootDir)
    for exposure in range(BIAS_SET_SIZE):
        _compressed_frame(tmp_path, rootPath / f"bias{exposure}.fits.Z", biasExposure=exposure)

    organiser.prepare(report=False)

    sofPath = rootPath / "sessions" / "base" / "sof" / MBIAS_SOF
    assert sofPath.read_text().splitlines() == [
        f"./raw/{NIGHT}/bias{exposure}.fits.Z  BIAS_VIS" for exposure in range(BIAS_SET_SIZE)
    ]
    collection, _ = set_of_files(
        log=log, settings=_sof_settings(tmp_path), inputFrames=str(sofPath), verbose=False
    ).get()
    assert list(collection.summary["file"]) == [f"bias{exposure}.fits.Z" for exposure in range(BIAS_SET_SIZE)]
    assert list(rootPath.rglob("*.fits")) == []


def test_prepare_keeps_the_compressed_frame_and_deletes_its_uncompressed_twin_at_the_root(
    tmp_path, organiser
) -> None:
    rootPath = Path(organiser.rootDir)
    _compressed_frame(tmp_path, rootPath / "bias.fits.Z")
    harvestable_raw_fits(rootPath / "bias.fits")

    organiser.prepare(report=False)

    assert _indexed_frames(organiser) == [{"file": "bias.fits.Z", "filepath": f"./raw/{NIGHT}/bias.fits.Z"}]
    assert list(rootPath.rglob("bias.fits")) == []


def test_prepare_deletes_an_uncompressed_twin_found_in_the_raw_tree(tmp_path, organiser) -> None:
    rootPath = Path(organiser.rootDir)
    nestedPath = rootPath / "raw" / NIGHT
    nestedPath.mkdir(parents=True)
    _compressed_frame(tmp_path, nestedPath / "bias.fits.Z")
    harvestable_raw_fits(nestedPath / "bias.fits")

    organiser.prepare(report=False)

    assert _indexed_frames(organiser) == [{"file": "bias.fits.Z", "filepath": f"./raw/{NIGHT}/bias.fits.Z"}]
    assert list(rootPath.rglob("bias.fits")) == []


def test_a_later_prepare_deletes_an_uncompressed_copy_of_an_indexed_compressed_frame(tmp_path, organiser) -> None:
    rootPath = Path(organiser.rootDir)
    _compressed_frame(tmp_path, rootPath / "bias.fits.Z")
    organiser.prepare(report=False)

    harvestable_raw_fits(rootPath / "bias.fits")
    organiser.prepare(report=False)

    assert _indexed_frames(organiser) == [{"file": "bias.fits.Z", "filepath": f"./raw/{NIGHT}/bias.fits.Z"}]
    assert list(rootPath.rglob("bias.fits")) == []
