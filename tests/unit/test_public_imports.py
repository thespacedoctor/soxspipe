"""Characterization tests for supported import paths."""

import pytest

pytestmark = pytest.mark.unit


def test_commonutils_public_imports_remain_available() -> None:
    from soxspipe.commonutils import detector_lookup, filenamer, keyword_lookup

    assert callable(detector_lookup)
    assert callable(filenamer)
    assert callable(keyword_lookup)


def test_recipe_public_imports_include_every_recipe() -> None:
    from soxspipe import recipes

    expectedNames = {
        "base_recipe",
        "soxs_disp_solution",
        "soxs_mbias",
        "soxs_mdark",
        "soxs_mflat",
        "soxs_nod",
        "soxs_offset",
        "soxs_order_centres",
        "soxs_spatial_solution",
        "soxs_stare",
        "soxs_straighten",
    }

    assert expectedNames <= set(vars(recipes))


def test_utkit_compatibility_import_remains_available() -> None:
    from soxspipe.utKit import utKit

    assert callable(utKit)
