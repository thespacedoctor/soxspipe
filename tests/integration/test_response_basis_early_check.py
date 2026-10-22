"""The stare, nod and offset recipes refuse an old-basis RESP product before any extraction runs (DY-1719).

`flux_calibration.load_response_factors` refuses a RESP product without the per-unbinned-pixel basis card, but it
runs only after the full extraction. Each observing recipe locates the RESP product long before that, so it checks
the card there and fails at once with the same message.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pytest

from soxspipe.commonutils.response_basis import RESPONSE_BASIS, require_response_basis
from tests.factories import response_product_fits
from tests.integration.test_observing_recipe_orchestration import _route
from tests.integration.test_soxs_nod_reduction_characterization import SINGLE_PAIR as NOD_PAIR
from tests.integration.test_soxs_nod_reduction_characterization import _nod_recipe
from tests.integration.test_soxs_nod_reduction_characterization import _patch as _patch_nod
from tests.integration.test_soxs_offset_reduction_characterization import SINGLE_PAIR as OFFSET_PAIR
from tests.integration.test_soxs_offset_reduction_characterization import _offset_recipe
from tests.integration.test_soxs_offset_reduction_characterization import _patch as _patch_offset
from tests.integration.test_soxs_stare_reduction_characterization import (
    _patch_stare,
    _spy_response_and_plot,
    _stare_recipe,
)

pytestmark = pytest.mark.integration

# A RESP PRODUCT FROM BEFORE THE BASIS CARD EXISTED, AND ONE THAT STATES A BASIS THIS RELEASE DOES NOT USE
OLD_BASES = [None, "PER_BINNED_PIXEL"]


def _old_response(tmpPath: Path, basis: str | None) -> tuple[Path, dict[frozenset[tuple[str, object]], list[str]]]:
    """Write an old-basis RESP product and return its path and the SOF route to it."""
    responsePath = response_product_fits(tmpPath / "RESP_TAB_VIS.fits", basis=basis)
    return responsePath, {_route(PRO_CATG="RESP_TAB_VIS"): [str(responsePath)]}


def _expected_message(responsePath: Path, basis: str | None) -> str:
    """Return the message `require_response_basis` gives for this product."""
    header = {} if basis is None else {"ESO PRO RESP BASIS": basis}
    with pytest.raises(ValueError) as caught:
        require_response_basis(header, str(responsePath))
    return str(caught.value)


def _record_extractor_construction(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, object]]:
    """Replace `horne_extraction` with a recorder and return the list of its constructions."""
    constructed: list[dict[str, object]] = []

    class RecordingExtractor:
        def __init__(self, **kwargs: object) -> None:
            constructed.append(kwargs)
            raise AssertionError("the extraction was constructed")

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "horne_extraction", RecordingExtractor)
    return constructed


@pytest.mark.parametrize("basis", OLD_BASES)
def test_stare_refuses_an_old_basis_response_before_the_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, basis: str | None
) -> None:
    # ARRANGE
    responsePath, routes = _old_response(tmp_path, basis)
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)
    constructed = _record_extractor_construction(monkeypatch)

    # ACT
    with pytest.raises(ValueError) as caught:
        recipe.produce_product()

    # ASSERT: THE SAME MESSAGE AS THE FLUX-CALIBRATION GUARD, AND NO EXTRACTION OR DETRENDING WAS STARTED
    assert str(caught.value) == _expected_message(responsePath, basis)
    assert constructed == []
    assert captured["extractor"] == []
    assert captured["detrend"] == []


@pytest.mark.parametrize("basis", OLD_BASES)
def test_nod_refuses_an_old_basis_response_before_the_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, basis: str | None
) -> None:
    # ARRANGE
    responsePath, routes = _old_response(tmp_path, basis)
    recipe, _ = _nod_recipe(log, tmp_path, NOD_PAIR, extraRoutes=routes)
    calls, captured = _patch_nod(recipe, monkeypatch, tmp_path)
    constructed = _record_extractor_construction(monkeypatch)

    # ACT
    with pytest.raises(ValueError) as caught:
        recipe.produce_product()

    # ASSERT
    assert str(caught.value) == _expected_message(responsePath, basis)
    assert constructed == []
    assert captured["extract_cycle"] == []
    assert "extract_cycle" not in calls
    assert "stack_extractions" not in calls


@pytest.mark.parametrize("basis", OLD_BASES)
def test_offset_refuses_an_old_basis_response_before_the_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, basis: str | None
) -> None:
    # ARRANGE
    responsePath, routes = _old_response(tmp_path, basis)
    recipe, _ = _offset_recipe(log, tmp_path, OFFSET_PAIR, extraRoutes=routes)
    calls, captured = _patch_offset(recipe, monkeypatch, tmp_path)
    constructed = _record_extractor_construction(monkeypatch)

    # ACT
    with pytest.raises(ValueError) as caught:
        recipe.produce_product()

    # ASSERT
    assert str(caught.value) == _expected_message(responsePath, basis)
    assert constructed == []
    assert captured["extract_cycle"] == []
    assert "extract_cycle" not in calls
    assert "stack_extractions" not in calls


def test_stare_extracts_when_the_response_has_the_current_basis(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The control case: a current-basis RESP product passes the check and the extraction runs."""
    # ARRANGE
    responsePath = response_product_fits(tmp_path / "RESP_TAB_VIS.fits", basis=RESPONSE_BASIS)
    routes = {_route(PRO_CATG="RESP_TAB_VIS"): [str(responsePath)]}
    recipe = _stare_recipe(log, tmp_path, extraRoutes=routes)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)
    calls = _spy_response_and_plot(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["extractor"]) == 1
    assert calls["fluxCal"][0][0] == str(responsePath)
