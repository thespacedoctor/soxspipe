"""The nod recipe blocks of the default settings ship ``reduce-as-stare: False`` and document it (DY-2194)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
CASES = [
    ("soxs_default_settings.yaml", "soxs-nod"),
    ("soxs_default_settings.yaml", "soxs-nod-std"),
    ("xsh_default_settings.yaml", "soxs-nod"),
]


@pytest.mark.parametrize(("fileName", "recipe"), CASES)
def test_each_nod_block_reduces_in_nod_mode_by_default(fileName: str, recipe: str) -> None:
    # ACT
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    # ASSERT
    assert settings[recipe]["reduce-as-stare"] is False


@pytest.mark.parametrize(("fileName", "recipe"), CASES)
def test_each_reduce_as_stare_key_has_an_uppercase_comment_above_it(fileName: str, recipe: str) -> None:
    # ARRANGE
    lines = (PACKAGE_DIR / fileName).read_text().splitlines()
    blockStart = lines.index(f"{recipe}:")

    # ACT
    keyIndex = next(i for i in range(blockStart, len(lines)) if lines[i].strip().startswith("reduce-as-stare:"))
    commentLines = []
    for line in reversed(lines[:keyIndex]):
        if not line.strip().startswith("#"):
            break
        commentLines.append(line.strip())
    comment = " ".join(reversed(commentLines))

    # ASSERT
    assert comment
    assert comment == comment.upper()
    assert "STARE" in comment


def test_the_xshooter_defaults_gain_no_standard_nod_block() -> None:
    # ACT
    settings = yaml.safe_load((PACKAGE_DIR / "xsh_default_settings.yaml").read_text())

    # ASSERT
    assert "soxs-nod-std" not in settings
