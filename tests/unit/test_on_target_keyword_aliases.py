"""The keyword aliases that the on-target frame list reads resolve for SOXS and X-Shooter (DY-2356)."""

from __future__ import annotations

from typing import Any

import pytest

from soxspipe.commonutils.keyword_lookup import keyword_lookup

pytestmark = pytest.mark.unit

ALIASES = {
    "soxs": {
        "ARCFILE": "ARCFILE",
        "ORIGFILE": "ORIGFILE",
        "OBS_PROG_ID": "ESO OBS PROG ID",
        "OBS_ID": "ESO OBS ID",
        "TPL_START": "ESO TPL START",
        "DET_DIT": "ESO DET SEQ1 DIT",
        "DET_NDIT": "ESO DET NDIT",
        "MJDOBS": "MJD-OBS",
    },
    "xsh": {
        "ARCFILE": "ARCFILE",
        "ORIGFILE": "ORIGFILE",
        "OBS_PROG_ID": "ESO OBS PROG ID",
        "OBS_ID": "ESO OBS ID",
        "TPL_START": "ESO TPL START",
        "DET_DIT": "ESO DET DIT",
        "DET_NDIT": "ESO DET NDIT",
        "MJDOBS": "MJD-OBS",
    },
}


@pytest.mark.parametrize("instrument", sorted(ALIASES))
def test_the_frame_list_aliases_resolve_to_the_raw_header_cards(log: Any, instrument: str) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings={"instrument": instrument}).get

    # ACT
    resolved = {alias: kw(alias) for alias in ALIASES[instrument]}

    # ASSERT
    assert resolved == ALIASES[instrument]
