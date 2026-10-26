"""Truth table for the run-serially rule held by each class that owns it (DY-1607)."""

from __future__ import annotations

import pytest

from soxspipe.commonutils.create_dispersion_map import create_dispersion_map
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_disp_solution import soxs_disp_solution
from soxspipe.recipes.soxs_order_centres import soxs_order_centres
from soxspipe.recipes.soxs_spatial_solution import soxs_spatial_solution

pytestmark = pytest.mark.unit

RULE_OWNERS = [base_recipe, create_dispersion_map, soxs_disp_solution, soxs_order_centres, soxs_spatial_solution]


@pytest.mark.parametrize("ownerClass", RULE_OWNERS, ids=lambda c: c.__name__)
@pytest.mark.parametrize(
    ("turnOffMP", "debug", "expectedRunSerially"),
    [(False, False, False), (True, False, True), (False, True, True), (True, True, True)],
)
def test_run_serially_when_multiprocessing_is_off_or_debug_is_on(
    ownerClass: type, turnOffMP: bool, debug: bool, expectedRunSerially: bool
) -> None:
    owner = object.__new__(ownerClass)
    owner.turnOffMP = turnOffMP
    owner.debug = debug

    assert owner.runSerially is expectedRunSerially


@pytest.mark.parametrize("ownerClass", RULE_OWNERS, ids=lambda c: c.__name__)
def test_run_serially_is_read_only(ownerClass: type) -> None:
    owner = object.__new__(ownerClass)
    owner.turnOffMP = False
    owner.debug = False

    with pytest.raises(AttributeError):
        owner.runSerially = True
