"""Characterization of `soxs_mflat.verify_input_frames`.

These tests pin what the recipe's own verification logic does today, before a
later commit splits `soxs_mflat.py`'s functions into smaller methods. Every scenario
below calls the real `verify_input_frames` directly on a bare recipe built
with `soxs_mflat.__new__`, with the inherited `base_recipe.
_verify_input_frames_basics` stubbed to return the classifications the
scenario needs and to set `self.arm`/`self.inst`, exactly as the real
collaborator does.

Several of the pinned messages come from a `"..." % locals()` pattern with no
`%` conversion specifiers in the string, an F507/UP031 lint site the tests
guard against being "fixed" into a different rendered string during the
split. Others pin a bug where a validation loop keeps overwriting its error
variable, so only the *last* offending value survives in the message -- pinned
as found, not as intended.
"""

from __future__ import annotations

from typing import Any

import pytest

from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_mflat import soxs_mflat

pytestmark = pytest.mark.unit

# THE MESSAGE THE NIR BRANCH RAISES WHEN NEITHER LAMP-ON TYPE IS PRESENT.
NIR_LAMP_MISSING_MESSAGE = "Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR"

# THE MESSAGE EVERY UVB/VIS VALIDATION FAILURE SHARES, REGARDLESS OF WHICH OF
# THE THREE CHECKS TRIPPED IT.
UVB_VIS_SHARED_MESSAGE = (
    "Input frames for soxspipe mflat need to be flat-lamp frames,"
    "a master-bias frame, an order-locations tables and possibly a master dark for UVB/VIS"
)


class _SubCollection:
    """A single filtered slice of the recipe's input-frame inventory."""

    def __init__(self, files: list[str], exptimes: list[float]) -> None:
        self.files = files
        self._exptimes = exptimes

    def values(self, keyword: str, unique: bool) -> list[float]:
        assert unique is True
        return list(dict.fromkeys(self._exptimes))


class _VerificationInventory:
    """Inventory double exposing exactly the surface `verify_input_frames` reads.

    `filter` is keyed by the exact filter dict a caller passes (order
    independent), so a scenario can either serve a specific lamp's exptime
    collection or leave every filter unmatched, which behaves like an
    ordinary empty `ImageFileCollection.filter(...)` result.
    """

    def __init__(
        self,
        *,
        matches: dict[tuple[tuple[str, str], ...], dict[str, Any]] | None = None,
        summary: str = "FAKE SUMMARY TABLE",
    ) -> None:
        self.summary = summary
        self._matches = {} if matches is None else matches
        self.filterCalls: list[dict[str, str]] = []

    def filter(self, **filters: str) -> _SubCollection:
        self.filterCalls.append(dict(filters))
        matched = self._matches.get(tuple(sorted(filters.items())))
        if matched is None:
            return _SubCollection(files=[], exptimes=[])
        return _SubCollection(files=matched["files"], exptimes=matched.get("exptimes", []))


def _stub_basics(
    monkeypatch: pytest.MonkeyPatch,
    *,
    imageTypes: list[str],
    imageTech: list[str],
    imageCategories: list[str],
    arm: str = "VIS",
    inst: str = "SOXS",
) -> None:
    """Replace the inherited basic verification with the classifications a scenario needs."""

    def fake_basics(self: soxs_mflat) -> tuple[list[str], list[str], list[str]]:
        self.arm = arm
        self.inst = inst
        return list(imageTypes), list(imageTech), list(imageCategories)

    monkeypatch.setattr(base_recipe, "_verify_input_frames_basics", fake_basics)


def _bare_recipe(log: Any, inventory: _VerificationInventory) -> soxs_mflat:
    """Return a recipe carrying only what `verify_input_frames` reads."""
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.inputFrames = inventory
    return recipe


def test_nir_without_a_lamp_flat_type_raises_the_locals_formatted_message(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Neither `LAMP,FLAT` nor `FLAT,LAMP` in the NIR image types raises the exact pinned text."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["DARK"],
        imageTech=["ECHELLE,SLIT", "IMAGE"],
        imageCategories=["ORDER_TAB_NIR"],
        arm="NIR",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == NIR_LAMP_MISSING_MESSAGE


def test_nir_unrecognised_tech_reports_only_the_last_offending_value(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With two unrecognised techs, the loop overwrites the message, so only the last survives."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["BADONE", "BADTWO"],
        imageCategories=["ORDER_TAB_NIR"],
        arm="NIR",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == (
        "Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR. You have provided BADTWO"
    )


def test_nir_missing_tech_reports_only_the_last_missing_value(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With both required techs missing, the loop overwrites the message, so `IMAGE` wins last."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=[],
        imageCategories=["ORDER_TAB_NIR"],
        arm="NIR",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == (
        "Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR. "
        "You have are missing TECH=IMAGE"
    )


@pytest.mark.parametrize(
    "imageTypes",
    [
        pytest.param(["LAMP,FLAT", "DARK"], id="dark_mix_of_length_two"),
        pytest.param(["LAMP,FLAT", "OTHER"], id="non_dark_mix_of_length_two"),
    ],
)
def test_nir_mixed_image_types_of_length_two_pass_verification_either_way(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageTypes: list[str],
) -> None:
    """The mixed-type rejection is commented out in the source, so both mixes pass silently.

    Pinned as found, not as intended: a length-two mix that is not `DARK`
    should, by the surrounding comment's own logic, be rejected, but the
    `error = "Input frames are a mix of %(imageTypes)s" % locals()` line is
    dead code, so neither mix ever raises.
    """
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=imageTypes,
        imageTech=["ECHELLE,SLIT", "IMAGE"],
        imageCategories=["ORDER_TAB_NIR"],
        arm="NIR",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == imageTypes[0]


def test_uvb_vis_rejects_a_non_lamp_image_type(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A UVB/VIS image type outside the five recognised lamp/dome types raises the shared message."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["DARK"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == UVB_VIS_SHARED_MESSAGE


def test_uvb_vis_rejects_a_missing_master_bias_or_order_table_category(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A UVB/VIS input missing `MASTER_BIAS_<arm>`/`ORDER_TAB_<arm>` raises the shared message."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=[],
        arm="VIS",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == UVB_VIS_SHARED_MESSAGE


def test_uvb_vis_rejects_an_empty_image_type_list_via_the_not_found_branch(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An empty `imageTypes` list is the only way to reach the "no lamp found" branch.

    The first per-type loop is vacuously true over an empty list, and the
    category check passes, so only the `found` guard trips -- raising the
    same shared message text as the other two UVB/VIS checks.
    """
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=[],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == UVB_VIS_SHARED_MESSAGE


def test_x_shooter_rejects_lamp_types_missing_the_other_dq_lamp_reporting_the_last_lamp_type(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """X-Shooter with only one of D2/QTH lamp flats fails, naming the last `LAMP` type found."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,QFLAT", "LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
        inst="XSH",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == 'Only "LAMP,FLAT" image types found. Please include both D2 and QTH lamp flats'


def test_x_shooter_with_both_dq_lamp_flats_passes_the_dq_check(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """X-Shooter carrying both D2 and QTH lamp flats passes the D/Q-lamp completeness check."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,QFLAT", "LAMP,DFLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
        inst="XSH",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "LAMP,QFLAT"


def test_nir_without_an_order_table_category_reports_the_arm_specific_message(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A NIR input with valid types/techs but no `ORDER_TAB_NIR` category fails the order-table check."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT", "IMAGE"],
        imageCategories=[],
        arm="NIR",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == "Need an order centre for NIR - none found with the input files"


def test_a_lamp_type_with_two_distinct_exptimes_raises_the_exptime_message(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A lamp collection reporting two distinct unique exposure times fails the exptime check."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    inventory = _VerificationInventory(
        matches={
            (("DPR_TECH", "ECHELLE,SLIT"), ("DPR_TYPE", "LAMP,FLAT")): {
                "files": ["a.fits", "b.fits"],
                "exptimes": [10.0, 20.0],
            },
        }
    )
    recipe = _bare_recipe(log, inventory)

    # ACT / ASSERT
    with pytest.raises(TypeError) as excinfo:
        recipe.verify_input_frames()
    assert str(excinfo.value) == "Input LAMP,FLAT frames for soxspipe mflat need to have a unique exptime"


def test_the_exptime_check_filters_each_lamp_type_in_order_on_a_passing_call(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A passing call still filters all four lamp types, in the fixed source order, for exptime."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    inventory = _VerificationInventory()
    recipe = _bare_recipe(log, inventory)

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert inventory.filterCalls == [
        {"DPR_TYPE": "LAMP,FLAT", "DPR_TECH": "ECHELLE,SLIT"},
        {"DPR_TYPE": "LAMP,DFLAT", "DPR_TECH": "ECHELLE,SLIT"},
        {"DPR_TYPE": "LAMP,QFLAT", "DPR_TECH": "ECHELLE,SLIT"},
        {"DPR_TYPE": "DOME,FLAT", "DPR_TECH": "ECHELLE,SLIT"},
    ]


def test_a_failed_verification_prints_the_error_banner_then_the_summary_then_a_blank_line(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """On failure, the recipe prints an error banner, the inventory summary, then an empty string."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["DARK"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    recipe = _bare_recipe(log, _VerificationInventory(summary="FAKE SUMMARY TABLE"))

    # ACT
    with pytest.raises(TypeError):
        recipe.verify_input_frames()

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed[-3:] == ["# VERIFYING INPUT FRAMES - **ERROR**\n", "FAKE SUMMARY TABLE", ""]


def test_a_passing_uvb_vis_call_records_the_first_image_type(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A fully passing UVB/VIS verification sets `imageType` to the first declared image type."""
    # ARRANGE
    _stub_basics(
        monkeypatch,
        imageTypes=["DOME,FLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    recipe = _bare_recipe(log, _VerificationInventory())

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "DOME,FLAT"
