"""Characterization of `base_recipe._verify_input_frames_basics`.

These tests pin what the method does today, before it is split into smaller
methods. The method is 232 lines long and most of its uncovered lines are the
rejection branches: every "input frames are a mix of X" guard, the NIR paths
that read the detector parameter file instead of the FITS header, and the one
exemption that lets a SOXS NIR nodding recipe carry two slit widths.

The method mutates `self` rather than returning its results, so each test
asserts the attributes it set as well as the exception it raised.

The mixed-arm guard used to interpolate `%(imageTypes)s` against `locals()`
before `imageTypes` existed, so it raised `KeyError` instead of the
`TypeError` it read as raising (DY-89). It now raises `TypeError` naming the
arms that differ.
"""

from __future__ import annotations

import importlib
from typing import Any

import pytest
from astropy.table import Table

from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

RECIPE_MODULE = importlib.import_module("soxspipe.recipes.base_recipe")

# THE COLUMNS `_verify_input_frames_basics` READS OUT OF THE INVENTORY
# SUMMARY. A COLUMN MISSING HERE FAILS THE METHOD WITH A `KeyError` RATHER
# THAN THE GUARD UNDER TEST.
UNIFORM_VIS_SUMMARY = {
    "SEQ_ARM": ["VIS"],
    "INSTRUME": ["SOXS"],
    "PRO_CATG": ["SCIENCE_VIS"],
    "CDELT1": [1],
    "CDELT2": [1],
    "DET_READ_SPEED": ["fast"],
    "CONAD": [1.1],
    "GAIN": [1.0],
    "DPR_TYPE": ["OBJECT"],
    "SLIT_VIS": ["SLIT1.0"],
    "RON": [3.5],
    "PRO_TYPE": [None],
    "DPR_TECH": ["ECHELLE,SLIT"],
    "PRO_TECH": [None],
    "DPR_CATG": ["SCIENCE"],
}


class FrameInventory:
    """Inventory exposing the two `ImageFileCollection` seams the method uses."""

    def __init__(self, summary: Table, *, files: list[str] | None = None) -> None:
        self.summary = summary
        self._files = ["science.fits"] if files is None else files

    def files_filtered(self, *, include_path: bool) -> list[str]:
        assert include_path is True
        return self._files

    def values(self, keyword: str, *, unique: bool) -> list[Any]:
        assert unique is True
        # `ImageFileCollection.values` DE-DUPLICATES WHILE KEEPING THE ORDER
        # OF FIRST APPEARANCE, WHICH IS WHAT THE GUARDS COUNT.
        seen: list[Any] = []
        for value in list(self.summary[keyword]):
            cleaned = None if value is None else value
            if cleaned not in seen:
                seen.append(cleaned)
        return seen


def _detector_lookup(
    *,
    dispersionAxis: str = "x",
    gain: float = 2.0,
    ron: float = 4.0,
) -> type:
    """Return a detector lookup class serving one fixed parameter set."""

    class DetectorLookup:
        """Deterministic replacement for the resource-backed lookup."""

        def __init__(self, **_: object) -> None:
            return None

        def get(self, arm: str) -> dict[str, object]:
            return {"dispersion-axis": dispersionAxis, "gain": gain, "ron": ron}

    return DetectorLookup


# THE BANNER EVERY REJECTION PRINTS BEFORE IT RAISES. A REFACTOR THAT DROPS OR
# REORDERS IT CHANGES WHAT THE USER SEES, SO EACH REJECTION TEST ASSERTS IT.
ERROR_BANNER = "# VERIFYING INPUT FRAMES - **ERROR**\n"


def _printed(log: Any) -> list[str]:
    """Return everything the recipe printed, in order.

    The recording logger stores what it was handed as a string, so a printed
    summary table is compared as its rendered text.
    """
    return [message for level, message in log.messages if level == "print"]


def _recipe(
    log: Any,
    summaryColumns: dict[str, list[Any]],
    *,
    instrument: str = "soxs",
    recipeName: str = "soxs-stare",
    files: list[str] | None = None,
) -> base_recipe:
    """Build a recipe carrying only the attributes the verification reads."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.kw = lambda name: name
    recipe.settings = {"instrument": instrument}
    recipe.recipeName = recipeName
    recipe.get_recipe_settings = lambda: {"frame-clipping-sigma": 3.0}
    recipe.inputFrames = FrameInventory(Table(summaryColumns), files=files)
    return recipe


def test_empty_input_inventory_is_rejected_before_any_header_is_read(
    log: Any,
) -> None:
    """An inventory with no files fails with `FileNotFoundError`."""
    # ARRANGE
    recipe = _recipe(log, UNIFORM_VIS_SUMMARY, files=[])

    # ACT / ASSERT
    with pytest.raises(FileNotFoundError, match="No image frames where passed"):
        recipe._verify_input_frames_basics()

    assert _printed(log) == [ERROR_BANNER]


def test_mixed_arms_raise_type_error_naming_the_arms(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The mixed-arm guard fails with `TypeError` naming the arms that differ."""
    # ARRANGE
    summary = {**UNIFORM_VIS_SUMMARY, "SEQ_ARM": ["VIS", "NIR"]}
    summary = {key: (value * 2)[:2] for key, value in summary.items()}
    summary["SEQ_ARM"] = ["VIS", "NIR"]
    recipe = _recipe(log, summary)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT / ASSERT
    with pytest.raises(TypeError, match="mix of arms: VIS and NIR"):
        recipe._verify_input_frames_basics()


def test_a_dispersion_axis_of_y_swaps_the_recorded_image_axes(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A `y` dispersion axis sets `axisA` to `y` and `axisB` to `x`."""
    # ARRANGE
    recipe = _recipe(log, UNIFORM_VIS_SUMMARY)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup(dispersionAxis="y"))

    # ACT
    recipe._verify_input_frames_basics()

    # ASSERT
    assert recipe.axisA == "y"
    assert recipe.axisB == "x"


def test_mixed_binning_is_rejected(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Two different `CDELT1` values fail with `TypeError`."""
    # ARRANGE
    summary = {key: value * 2 for key, value in UNIFORM_VIS_SUMMARY.items()}
    summary["CDELT1"] = [1, 2]
    recipe = _recipe(log, summary)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT / ASSERT
    with pytest.raises(TypeError, match="mix of binnings"):
        recipe._verify_input_frames_basics()

    assert _printed(log) == [ERROR_BANNER]


def test_mixed_readout_speeds_are_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two different readout speeds fail with `TypeError` naming the values."""
    # ARRANGE
    summary = {key: value * 2 for key, value in UNIFORM_VIS_SUMMARY.items()}
    summary["DET_READ_SPEED"] = ["fast", "slow"]
    recipe = _recipe(log, summary)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT / ASSERT
    with pytest.raises(TypeError, match="mix of readout speeds"):
        recipe._verify_input_frames_basics()

    printed = _printed(log)
    assert printed[0] == ERROR_BANNER
    assert printed[1] == str(recipe.inputFrames.summary)
    assert printed[2] == "\n\n"
    assert len(printed) == 3


def test_mixed_gain_is_rejected(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Two different gains fail with `TypeError` naming both values."""
    # ARRANGE
    summary = {key: value * 2 for key, value in UNIFORM_VIS_SUMMARY.items()}
    summary["CONAD"] = [1.1, 2.2]
    summary["GAIN"] = [1.0, 2.0]
    recipe = _recipe(log, summary)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT / ASSERT
    with pytest.raises(TypeError, match="mix of gain"):
        recipe._verify_input_frames_basics()

    assert _printed(log) == [ERROR_BANNER, str(recipe.inputFrames.summary)]


def test_xsh_reads_its_gain_from_the_conad_keyword_alone(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Under `instrument = xsh` the gain comes from `CONAD`, not `max(CONAD, GAIN)`."""
    # ARRANGE
    summary = {**UNIFORM_VIS_SUMMARY, "CONAD": [1.1], "GAIN": [9.9]}
    recipe = _recipe(log, summary, instrument="xsh")
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT
    recipe._verify_input_frames_basics()

    # ASSERT
    assert recipe.detectorParams["gain"].value == pytest.approx(1.1)


def test_nir_frames_take_gain_and_readnoise_from_the_detector_parameters(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A NIR inventory is never binned and reads its gain and readnoise from file."""
    # ARRANGE
    summary = {
        **UNIFORM_VIS_SUMMARY,
        "SEQ_ARM": ["NIR"],
        "PRO_CATG": ["SCIENCE_NIR"],
        "CONAD": [None],
        "GAIN": [None],
        "RON": [None],
        "SLIT_NIR": ["SLIT1.0"],
    }
    del summary["SLIT_VIS"]
    recipe = _recipe(log, summary, recipeName="soxs-nod")
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup(gain=2.5, ron=7.0))

    # ACT
    recipe._verify_input_frames_basics()

    # ASSERT
    assert recipe.detectorParams["binning"] == [1, 1]
    assert recipe.detectorParams["gain"].value == pytest.approx(2.5)
    assert recipe.detectorParams["ron"].value == pytest.approx(7.0)
    assert "Gain is being read from the detector parameter file (not the FITS header)" in "".join(
        str(message) for message in log.messages
    )


def test_mixed_slit_widths_are_rejected_outside_the_nodding_exemption(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two slit widths fail for a VIS staring recipe."""
    # ARRANGE
    summary = {key: value * 2 for key, value in UNIFORM_VIS_SUMMARY.items()}
    summary["SLIT_VIS"] = ["SLIT1.0", "SLIT5.0"]
    summary["DPR_TYPE"] = ["OBJECT", "OBJECT"]
    recipe = _recipe(log, summary)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT / ASSERT
    with pytest.raises(TypeError, match="mix of slit-width"):
        recipe._verify_input_frames_basics()

    assert _printed(log) == [ERROR_BANNER, str(recipe.inputFrames.summary)]


def test_a_soxs_nir_nodding_recipe_accepts_the_five_and_one_point_five_slit_pair(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`SLIT5.0` with `SLIT1.5` passes for SOXS NIR nodding, staring and offset."""
    # ARRANGE
    summary = {key: value * 2 for key, value in UNIFORM_VIS_SUMMARY.items()}
    summary["SEQ_ARM"] = ["NIR", "NIR"]
    summary["PRO_CATG"] = ["SCIENCE_NIR", "SCIENCE_NIR"]
    summary["SLIT_NIR"] = ["SLIT5.0", "SLIT1.5"]
    summary["DPR_TYPE"] = ["OBJECT", "OBJECT"]
    del summary["SLIT_VIS"]
    recipe = _recipe(log, summary, recipeName="soxs-nod-std")
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT
    imageTypes, imageTech, imageCat = recipe._verify_input_frames_basics()

    # ASSERT
    assert recipe.arm == "NIR"
    assert set(imageTypes) == {"OBJECT"}
    assert set(imageTech) == {"ECHELLE,SLIT"}
    assert set(imageCat) == {"SCIENCE", "SCIENCE_NIR"}


def test_mixed_readnoise_is_rejected(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Two different `RON` values fail with `TypeError` naming the values."""
    # ARRANGE
    summary = {key: value * 2 for key, value in UNIFORM_VIS_SUMMARY.items()}
    summary["RON"] = [3.5, 4.5]
    recipe = _recipe(log, summary)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT / ASSERT
    with pytest.raises(TypeError, match="mix of readnoise"):
        recipe._verify_input_frames_basics()

    assert _printed(log) == [ERROR_BANNER, str(recipe.inputFrames.summary)]


def test_the_returned_classifications_drop_none_and_the_reduced_marker(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`REDUCED` and `None` are stripped from the three returned lists."""
    # ARRANGE
    summary = {
        **UNIFORM_VIS_SUMMARY,
        "PRO_TYPE": ["REDUCED"],
        "PRO_TECH": ["ECHELLE,SLIT"],
    }
    recipe = _recipe(log, summary)
    monkeypatch.setattr(RECIPE_MODULE, "detector_lookup", _detector_lookup())

    # ACT
    imageTypes, imageTech, imageCat = recipe._verify_input_frames_basics()

    # ASSERT
    assert set(imageTypes) == {"OBJECT"}
    assert set(imageTech) == {"ECHELLE,SLIT"}
    assert set(imageCat) == {"SCIENCE", "SCIENCE_VIS"}
