"""The DATANCOM keyword alias resolves to the same ESO card for SOXS and X-Shooter."""

from __future__ import annotations

from typing import Any

import pytest

from soxspipe.commonutils.keyword_lookup import keyword_lookup

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("instrument", ["soxs", "xsh"])
def test_datancom_alias_resolves_to_the_eso_pro_card(log: Any, instrument: str) -> None:
    # ARRANGE
    lookup = keyword_lookup(log=log, settings={"instrument": instrument})

    # ACT
    card = lookup.get("DATANCOM")

    # ASSERT
    assert card == "ESO PRO DATANCOM"
