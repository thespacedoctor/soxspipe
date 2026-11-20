"""A SOXS NIR recipe takes its read noise from the sample count of its own frames (DY-2341).

The raw frames carry `ESO DET NDSAMPLES`. The read noise is looked up from that count in the NIR table of the SOXS
detector parameters, whether or not the SOF names a master dark. The master dark's `ESO QC RAW RON` sets nothing. The
stare sky fit adds the master dark's noise, looked up at the dark's own sample count and divided by the number of
frames the dark combined.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.set_of_files import ImageFileCollection
from soxspipe.recipes import (
    base_recipe,
    soxs_disp_solution,
    soxs_mdark,
    soxs_mflat,
    soxs_nod,
    soxs_offset,
    soxs_order_centres,
    soxs_spatial_solution,
    soxs_stare,
    soxs_straighten,
)

pytestmark = pytest.mark.unit

SOXS_SETTINGS = {"instrument": "soxs"}
TABLE = {2: 17.50, 4: 16.94, 7: 14.62, 10: 13.14, 39: 8.95, 77: 7.62, 128: 7.55}
FALLBACK_RON = 7.55
SCIENCE_NDSAMPLES = 128
DARK_NDSAMPLES = 10
DARK_NCOMBINED = 4
# A QC RAW RON FAR FROM EVERY TABLE VALUE, SO A RECIPE THAT READS IT SHOWS
DARK_QC_RAW_RON = 99.0
NIR_RAW = {"ESO DET CHIP RON": 0.0}
RECIPE_CLASSES = [
    soxs_disp_solution,
    soxs_mdark,
    soxs_mflat,
    soxs_nod,
    soxs_offset,
    soxs_order_centres,
    soxs_spatial_solution,
    soxs_stare,
    soxs_straighten,
]


def _write_frame(directory: Path, name: str, cards: dict[str, object]) -> None:
    header = fits.Header()
    header["INSTRUME"] = "SOXS"
    header["HIERARCH ESO SEQ ARM"] = "NIR"
    for key, value in cards.items():
        header[f"HIERARCH {key}"] = value
    directory.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=np.zeros((2, 2), dtype=np.float32), header=header).writeto(directory / name)


def _raw(ndsamples: int | None = SCIENCE_NDSAMPLES) -> dict[str, object]:
    cards: dict[str, object] = dict(NIR_RAW)
    if ndsamples is not None:
        cards["ESO DET NDSAMPLES"] = ndsamples
    return cards


def _dark(**overrides: object) -> dict[str, object]:
    cards: dict[str, object] = {
        "ESO PRO CATG": "MASTER_DARK_NIR",
        "ESO DET NDSAMPLES": DARK_NDSAMPLES,
        "ESO PRO DATANCOM": DARK_NCOMBINED,
        "ESO QC RAW RON": DARK_QC_RAW_RON,
    }
    cards.update(overrides)
    return {key: value for key, value in cards.items() if value is not None}


def _write_scene(
    directory: Path, rawCards: dict[str, object], darkCards: dict[str, object] | None, nRaw: int = 2
) -> None:
    for index in range(1, nRaw + 1):
        _write_frame(directory, f"raw_{index}.fits", rawCards)
    if darkCards is not None:
        _write_frame(directory, "master_dark.fits", darkCards)


def _recipe(
    directory: Path, log: Any, recipeClass: type[base_recipe] = soxs_stare, workspace: Path | None = None
) -> base_recipe:
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    recipe = recipeClass.__new__(recipeClass)
    recipe.log = log
    recipe.kw = kw
    recipe.arm = "NIR"
    recipe.settings = SOXS_SETTINGS
    recipe.inputFrames = ImageFileCollection(
        location=str(directory),
        keywords=[kw("SEQ_ARM"), kw("PRO_CATG"), kw("RON"), kw("CHIP_RON"), kw("DET_NDSAMPLES"), "file"],
    )
    recipe.detectorParams = detector_lookup(log=log, settings=SOXS_SETTINGS).get("NIR").copy()
    recipe.workspaceRootPath = str(workspace if workspace is not None else directory / "no_workspace")
    recipe._report_verification_error = lambda **_: None
    return recipe


def _verify_read_noise(recipe: base_recipe) -> base_recipe:
    recipe._verify_single_readnoise()
    recipe._apply_measured_read_noise()
    return recipe


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def _ron(recipe: base_recipe) -> float:
    return float(recipe.detectorParams["ron"].value)


@pytest.mark.parametrize("hasDark", [False, True], ids=["no master dark", "master dark"])
@pytest.mark.parametrize("recipeClass", RECIPE_CLASSES, ids=lambda recipeClass: recipeClass.__name__)
def test_every_soxs_nir_recipe_takes_the_read_noise_from_its_frames_ndsamples(
    tmp_path: Path, log: Any, recipeClass: type[base_recipe], hasDark: bool
) -> None:
    # ARRANGE: 10 SAMPLES GIVE 13.14 e-; THE DARK'S OWN QC RAW RON AND THE 128-SAMPLE FALLBACK DIFFER FROM IT
    _write_scene(tmp_path, _raw(DARK_NDSAMPLES), _dark(**{"ESO DET NDSAMPLES": 128}) if hasDark else None)
    recipe = _recipe(tmp_path, log, recipeClass)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert _ron(recipe) == pytest.approx(TABLE[10])
    assert recipe.detectorParams["ron"].unit == "electron"
    assert _warnings(log) == []


def test_a_stare_recipe_interpolates_the_read_noise_at_an_unmeasured_sample_count(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_scene(tmp_path, _raw(97), None)
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert 7.55 < _ron(recipe) < 7.62


def test_the_stare_sky_fit_adds_the_master_darks_noise_at_the_darks_own_sample_count(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE SCIENCE FRAMES HAVE 128 SAMPLES (7.55 e-); THE MASTER DARK HAS 10 (13.14 e-) AND COMBINES 4 DARKS
    _write_scene(tmp_path, _raw(SCIENCE_NDSAMPLES), _dark())
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    expected = math.sqrt(TABLE[128] ** 2 + TABLE[10] ** 2 / DARK_NCOMBINED)
    assert recipe.readNoise.skyFit == pytest.approx(expected)
    assert recipe.readNoise.singleFrame == pytest.approx(TABLE[128])
    assert recipe.readNoise.nCombined == DARK_NCOMBINED
    assert recipe.readNoise.source == "NDSAMPLES table"
    assert _ron(recipe) == pytest.approx(TABLE[128])


def test_the_stare_sky_fit_of_a_frame_without_a_master_dark_is_the_frames_own_read_noise(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE
    _write_scene(tmp_path, _raw(DARK_NDSAMPLES), None)
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.readNoise.skyFit == pytest.approx(TABLE[10])
    assert recipe.readNoise.nCombined is None
    assert recipe.readNoise.source == "NDSAMPLES table"


@pytest.mark.parametrize("recipeClass", [soxs_nod, soxs_offset], ids=["nod", "offset"])
def test_the_nod_sky_fit_is_the_frames_read_noise_whatever_master_dark_the_sof_names(
    tmp_path: Path, log: Any, recipeClass: type[base_recipe]
) -> None:
    # ARRANGE: A-B REMOVES THE DARK, SO TWO MASTER DARKS RAISE NO WARNING AND ADD NOTHING
    _write_scene(tmp_path, _raw(SCIENCE_NDSAMPLES), _dark())
    _write_frame(tmp_path, "master_dark_b.fits", _dark())
    recipe = _recipe(tmp_path, log, recipeClass)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.readNoise.skyFit == pytest.approx(TABLE[128])
    assert recipe.readNoise.nCombined is None
    assert _warnings(log) == []


def test_a_master_dark_without_the_ndsamples_card_takes_it_from_the_first_raw_frame_it_lists(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE: THE LISTED RAW FRAME SITS IN A DATED FOLDER OF THE WORKSPACE RAW DIRECTORY
    workspace = tmp_path / "workspace"
    _write_frame(workspace / "raw" / "2026-01-30", "SOXS.2026-01-31T19:16:38.123.fits", _raw(DARK_NDSAMPLES))
    sofDirectory = tmp_path / "sof"
    _write_scene(
        sofDirectory,
        _raw(SCIENCE_NDSAMPLES),
        _dark(**{"ESO DET NDSAMPLES": None, "ESO PRO REC1 RAW1 NAME": "SOXS.2026-01-31T19:16:38.123.fits"}),
    )
    recipe = _recipe(sofDirectory, log, soxs_stare, workspace)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    expected = math.sqrt(TABLE[128] ** 2 + TABLE[10] ** 2 / DARK_NCOMBINED)
    assert recipe.readNoise.skyFit == pytest.approx(expected)
    assert _warnings(log) == []


def test_a_master_dark_whose_sample_count_cannot_be_found_uses_the_fallback_noise_with_a_warning(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE: NO CARD ON THE DARK AND NO RAW FRAME LISTED
    _write_scene(tmp_path, _raw(SCIENCE_NDSAMPLES), _dark(**{"ESO DET NDSAMPLES": None}))
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    expected = math.sqrt(TABLE[128] ** 2 + FALLBACK_RON**2 / DARK_NCOMBINED)
    assert recipe.readNoise.skyFit == pytest.approx(expected)
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "master_dark.fits" in warnings[0]
    assert "NDSAMPLES" in warnings[0]


def test_a_master_dark_that_does_not_say_how_many_frames_it_combined_is_counted_as_one(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE
    _write_scene(tmp_path, _raw(SCIENCE_NDSAMPLES), _dark(**{"ESO PRO DATANCOM": None}))
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.readNoise.nCombined == 1
    assert recipe.readNoise.skyFit == pytest.approx(math.sqrt(TABLE[128] ** 2 + TABLE[10] ** 2))
    assert any("combined" in message for message in _warnings(log))


def test_an_unreadable_master_dark_header_leaves_the_frames_read_noise_and_warns(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    _write_scene(tmp_path, _raw(SCIENCE_NDSAMPLES), _dark())
    recipe = _recipe(tmp_path, log, soxs_stare)

    def unreadable(*_: object, **__: object) -> None:
        raise OSError("truncated file")

    monkeypatch.setattr(fits, "getheader", unreadable)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert _ron(recipe) == pytest.approx(TABLE[128])
    assert recipe.readNoise.skyFit == pytest.approx(TABLE[128])
    assert any("master dark" in message and "could not be read" in message for message in _warnings(log))


def test_several_master_darks_use_the_last_by_name_for_the_sky_fit_and_warn(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE LATER NAME HAS 2 SAMPLES (17.50 e-) AND COMBINES 1 DARK
    _write_scene(tmp_path, _raw(SCIENCE_NDSAMPLES), None)
    _write_frame(tmp_path, "master_dark_b.fits", _dark(**{"ESO DET NDSAMPLES": 2, "ESO PRO DATANCOM": 1}))
    _write_frame(tmp_path, "master_dark_a.fits", _dark())
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.readNoise.skyFit == pytest.approx(math.sqrt(TABLE[128] ** 2 + TABLE[2] ** 2))
    assert any("2 master darks" in message and "master_dark_b.fits" in message for message in _warnings(log))


@pytest.mark.parametrize("recipeClass", [soxs_stare, soxs_nod], ids=["stare", "nod"])
def test_a_raw_frame_without_the_ndsamples_card_uses_the_fallback_noise_and_warns_with_its_name(
    tmp_path: Path, log: Any, recipeClass: type[base_recipe]
) -> None:
    # ARRANGE
    _write_scene(tmp_path, _raw(None), None, nRaw=1)
    recipe = _recipe(tmp_path, log, recipeClass)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert _ron(recipe) == pytest.approx(FALLBACK_RON)
    assert recipe.readNoise.source == "detector default"
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "raw_1.fits" in warnings[0]


def test_raw_frames_with_different_sample_counts_use_the_larger_read_noise_and_warn(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_frame(tmp_path, "raw_1.fits", _raw(SCIENCE_NDSAMPLES))
    _write_frame(tmp_path, "raw_2.fits", _raw(DARK_NDSAMPLES))
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert _ron(recipe) == pytest.approx(TABLE[10])
    assert any("sample counts" in message for message in _warnings(log))


def test_a_master_dark_far_from_the_table_value_raises_no_nominal_disagreement_message(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE: THE DY-1581 15% MESSAGE COMPARED THE DARK'S QC RAW RON WITH THE NOMINAL VALUE
    _write_scene(tmp_path, _raw(SCIENCE_NDSAMPLES), _dark())
    recipe = _recipe(tmp_path, log, soxs_stare)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert not any("nominal" in message for _, message in log.messages)


@pytest.mark.parametrize("recipeClass", [soxs_stare, soxs_nod], ids=["stare", "nod"])
def test_an_xsh_nir_recipe_ignores_the_ndsamples_table(
    tmp_path: Path, log: Any, recipeClass: type[base_recipe]
) -> None:
    # ARRANGE: X-SHOOTER FRAMES CARRY NO TABLE; THE DARK STILL SETS THE STARE READ NOISE AND THE NOD KEEPS THE NOMINAL
    xshSettings = {"instrument": "xsh"}
    _write_scene(tmp_path, _raw(DARK_NDSAMPLES), _dark(**{"ESO QC RAW RON": 9.0}))
    kw = keyword_lookup(log=log, settings=xshSettings).get
    recipe = recipeClass.__new__(recipeClass)
    recipe.log = log
    recipe.kw = kw
    recipe.arm = "NIR"
    recipe.settings = xshSettings
    recipe.inputFrames = ImageFileCollection(
        location=str(tmp_path), keywords=[kw("SEQ_ARM"), kw("PRO_CATG"), kw("RON"), kw("CHIP_RON"), "file"]
    )
    recipe.detectorParams = {"ron": 8.0}
    recipe._report_verification_error = lambda **_: None

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    expected = 9.0 if recipeClass is soxs_stare else 8.0
    assert _ron(recipe) == pytest.approx(expected)


@pytest.mark.parametrize("isVerbose", [False, True], ids=["default keys", "verbose keys"])
def test_the_shipped_frame_summary_carries_the_ndsamples_column(tmp_path: Path, log: Any, isVerbose: bool) -> None:
    # ARRANGE: THE SHIPPED SUMMARY KEYS, AS A RECIPE READS THEM
    import yaml

    import soxspipe
    from soxspipe.commonutils.set_of_files import set_of_files

    advancedPath = Path(soxspipe.__file__).parent / "advanced_settings.yaml"
    settings = {**yaml.safe_load(advancedPath.read_text()), **SOXS_SETTINGS, "workspace-root-dir": str(tmp_path)}
    _write_frame(tmp_path / "sof", "raw_1.fits", _raw(DARK_NDSAMPLES))

    # ACT
    frames, _ = set_of_files(
        log=log, settings=settings, inputFrames=str(tmp_path / "sof"), ext=0, verbose=isVerbose
    ).get()

    # ASSERT
    assert list(frames.summary["ESO DET NDSAMPLES"]) == [DARK_NDSAMPLES]


def test_an_unreadable_raw_frame_listed_by_the_master_dark_warns_and_uses_the_fallback_noise(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE: THE LISTED RAW FRAME IS ON DISK BUT NOT A FITS FILE
    workspace = tmp_path / "workspace"
    rawDirectory = workspace / "raw" / "2026-01-30"
    rawDirectory.mkdir(parents=True)
    (rawDirectory / "SOXS.2026-01-31T19:16:38.123.fits").write_text("not a fits file")
    sofDirectory = tmp_path / "sof"
    _write_scene(
        sofDirectory,
        _raw(SCIENCE_NDSAMPLES),
        _dark(**{"ESO DET NDSAMPLES": None, "ESO PRO REC1 RAW1 NAME": "SOXS.2026-01-31T19:16:38.123.fits"}),
    )
    recipe = _recipe(sofDirectory, log, soxs_stare, workspace)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    expected = math.sqrt(TABLE[128] ** 2 + FALLBACK_RON**2 / DARK_NCOMBINED)
    assert recipe.readNoise.skyFit == pytest.approx(expected)
    warnings = _warnings(log)
    assert any("could not be read" in message for message in warnings)
    assert any("master_dark.fits" in message and "NDSAMPLES" in message for message in warnings)
