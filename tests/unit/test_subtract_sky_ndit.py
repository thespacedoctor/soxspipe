"""``subtract_sky`` reads the frame's NDIT for its sky-fit noise model (DY-2323)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tests.unit.test_subtract_sky_default_path_golden import ARMS, run_default_subtractor

pytestmark = pytest.mark.unit

NDIT = 4
NDIT_CARD = "ESO DET NDIT"


@pytest.mark.parametrize("arm", ARMS)
def test_the_subtractor_ndit_is_one_when_the_frame_has_no_ndit_card(log: Any, tmp_path: Path, arm: str) -> None:
    # ACT
    subtractor, _ = run_default_subtractor(log, arm, tmp_path)

    # ASSERT
    assert subtractor.ndit == 1


def test_the_subtractor_reads_the_ndit_card_of_the_frame(log: Any, tmp_path: Path) -> None:
    # ACT
    subtractor, _ = run_default_subtractor(log, "NIR", tmp_path, {NDIT_CARD: NDIT})

    # ASSERT
    assert subtractor.ndit == NDIT


@pytest.mark.parametrize("arm", ARMS)
def test_the_products_of_an_ndit_1_frame_are_bit_identical_to_those_of_a_frame_without_the_card(
    log: Any, tmp_path: Path, arm: str
) -> None:
    # ARRANGE
    withoutCard = run_default_subtractor(log, arm, tmp_path / "none")[1]

    # ACT
    withCard = run_default_subtractor(log, arm, tmp_path / "one", {NDIT_CARD: 1})[1]

    # ASSERT
    for name, values in withoutCard.items():
        np.testing.assert_array_equal(withCard[name], values, err_msg=name)
