"""Doubles for exercising `soxs_mflat.verify_input_frames` without a workspace."""

from __future__ import annotations

from typing import Any

import pytest

from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_mflat import soxs_mflat


class SubCollection:
    """A single filtered slice of the recipe's input-frame inventory."""

    def __init__(self, files: list[str], exptimes: list[float]) -> None:
        self.files = files
        self._exptimes = exptimes
        self.uniqueRequests: list[bool] = []

    def values(self, keyword: str, unique: bool) -> list[float]:
        self.uniqueRequests.append(unique)
        return list(dict.fromkeys(self._exptimes))


class VerificationInventory:
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
        self.subCollections: list[SubCollection] = []

    def filter(self, **filters: str) -> SubCollection:
        self.filterCalls.append(dict(filters))
        matched = self._matches.get(tuple(sorted(filters.items())))
        if matched is None:
            subCollection = SubCollection(files=[], exptimes=[])
        else:
            subCollection = SubCollection(files=matched["files"], exptimes=matched.get("exptimes", []))
        self.subCollections.append(subCollection)
        return subCollection


def stub_mflat_basics(
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


def bare_mflat_recipe(log: Any, inventory: VerificationInventory) -> soxs_mflat:
    """Return a recipe carrying only what `verify_input_frames` reads."""
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.inputFrames = inventory
    return recipe
