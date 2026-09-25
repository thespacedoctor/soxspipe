"""Characterization tests pinning the current behaviour of `set_of_files.py`.

These tests describe what the module does today, including behaviour that
looks like a defect (DY-73, DY-74, DY-75, and others noted inline). They must
all pass unchanged against the current implementation -- they are not
RED/GREEN tests, they exist to freeze behaviour ahead of a refactor.
"""

from __future__ import annotations

from pathlib import Path

import astropy.io.fits as astropy_fits
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.set_of_files import (
    ImageFileCollection,
    _join_fits_summaries_in_input_order,
    _supplementary_path_from_sof_line,
)
from tests.factories import raw_fits

pytestmark = pytest.mark.unit


def _settings(
    tmp_path: Path,
    *,
    default: list[str] | None = None,
    verbose: list[str] | None = None,
    nodding_extras: list[str] | None = None,
) -> dict[str, object]:
    """Return a minimal settings mapping with overridable summary-keys lists."""
    from tests.factories import pipeline_settings

    return pipeline_settings(
        tmp_path,
        overrides={
            "summary-keys": {
                "default": [] if default is None else default,
                "verbose": [] if verbose is None else verbose,
                "nodding_extras": [] if nodding_extras is None else nodding_extras,
            }
        },
    )


# ---------------------------------------------------------------------------
# MODULE-LEVEL HELPERS
# ---------------------------------------------------------------------------


def test_supplementary_path_returns_an_existing_path_unchanged(
    tmp_path: Path,
) -> None:
    existingPath = tmp_path / "existing.csv"
    existingPath.write_text("data", encoding="utf-8")

    result = _supplementary_path_from_sof_line(str(existingPath), str(tmp_path))

    assert result == str(existingPath)


def test_supplementary_path_strips_a_recognised_arm_tag(tmp_path: Path) -> None:
    line = f"{tmp_path / 'missing.csv'} ORDER_TAB_VIS"

    result = _supplementary_path_from_sof_line(line, str(tmp_path))

    assert result == str(tmp_path / "missing.csv")


def test_supplementary_path_keeps_a_lowercase_tag(tmp_path: Path) -> None:
    line = f"{tmp_path / 'missing.csv'} order_tab_vis"

    result = _supplementary_path_from_sof_line(line, str(tmp_path))

    assert result == line


def test_supplementary_path_keeps_a_tag_with_an_unrecognised_suffix(
    tmp_path: Path,
) -> None:
    line = f"{tmp_path / 'missing.csv'} ORDER_TAB_XYZ"

    result = _supplementary_path_from_sof_line(line, str(tmp_path))

    assert result == line


def test_supplementary_path_expands_a_leading_tilde(tmp_path: Path) -> None:
    result = _supplementary_path_from_sof_line("~/foo.csv", "/fake/home")

    assert result == "/fake/home/foo.csv"


def test_join_fits_summaries_preserves_extension_row_order() -> None:
    # THE EXTENSION TABLE'S ROW ORDER WINS, EVEN THOUGH THE PRIMARY TABLE IS
    # ORDERED DIFFERENTLY
    primarySummary = Table({"file": ["b", "a"], "primaryValue": [1, 2]})
    extensionSummary = Table({"file": ["a", "b"], "extensionValue": [10, 20]})

    joined = _join_fits_summaries_in_input_order(primarySummary, extensionSummary)

    assert list(joined["file"]) == ["a", "b"]
    assert list(joined["primaryValue"]) == [2, 1]
    assert list(joined["extensionValue"]) == [10, 20]
    assert joined.colnames == ["file", "primaryValue", "extensionValue"]


# ---------------------------------------------------------------------------
# ImageFileCollection._dict_from_fits_header
# ---------------------------------------------------------------------------


def test_dict_from_fits_header_uses_full_paths_when_files_span_two_directories(
    tmp_path: Path,
) -> None:
    firstDir = tmp_path / "dirA"
    secondDir = tmp_path / "dirB"
    firstDir.mkdir()
    secondDir.mkdir()
    firstPath = raw_fits(firstDir / "one.fits", seed=1)
    secondPath = raw_fits(secondDir / "two.fits", seed=2)

    collection = ImageFileCollection(
        filenames=[str(firstPath), str(secondPath)],
        location=None,
        keywords=["DPR_TYPE"],
    )

    # NO COMMON LOCATION: THE FILE COLUMN KEEPS THE FULL PATHS, IN INPUT ORDER
    assert list(collection.summary["file"]) == [str(firstPath), str(secondPath)]


def test_dict_from_fits_header_keeps_first_value_and_warns_on_duplicate_keyword(
    tmp_path: Path,
) -> None:
    filePath = tmp_path / "one.fits"
    header = fits.Header()
    header["DPR_TYPE"] = "BIAS"
    header.append(fits.Card("DUPKEY", "first"), useblanks=False, bottom=True)
    header.append(fits.Card("DUPKEY", "second"), useblanks=False, bottom=True)
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=header).writeto(filePath)

    with pytest.warns(UserWarning) as recordedWarnings:
        collection = ImageFileCollection(
            filenames=[str(filePath)],
            location=str(tmp_path),
            keywords=["DPR_TYPE", "DUPKEY"],
        )

    assert list(collection.summary["DUPKEY"]) == ["first"]
    assert str(recordedWarnings[0].message) == (
        f'Header from file "{filePath}" contains multiple entries for '
        f'"DUPKEY", the pair "DUPKEY=second" will be ignored.'
    )


def test_dict_from_fits_header_skips_a_blank_keyword(tmp_path: Path) -> None:
    filePath = tmp_path / "one.fits"
    header = fits.Header()
    header["DPR_TYPE"] = "BIAS"
    # A NON-EMPTY BLANK-KEYWORD CARD ROUND-TRIPS THROUGH WRITETO/GETHEADER
    # WITH AN EMPTY-STRING KEYWORD, WHICH MUST BE SKIPPED. (A TRULY EMPTY
    # `fits.Card("", "")` IS SILENTLY DROPPED BY ASTROPY ON WRITE, SO IT
    # CANNOT BE USED TO EXERCISE THIS BRANCH.)
    header.add_blank("free text")
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=header).writeto(filePath)

    collection = ImageFileCollection(
        filenames=[str(filePath)], location=str(tmp_path), keywords=["DPR_TYPE"]
    )

    assert list(collection.summary["DPR_TYPE"]) == ["BIAS"]


def test_dict_from_fits_header_marks_a_keyword_missing_from_a_later_file(
    tmp_path: Path,
) -> None:
    firstPath = tmp_path / "one.fits"
    secondPath = tmp_path / "two.fits"
    firstHeader = fits.Header({"DPR_TYPE": "BIAS", "EXTRA": 5})
    secondHeader = fits.Header({"DPR_TYPE": "BIAS"})
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=firstHeader).writeto(firstPath)
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=secondHeader).writeto(secondPath)

    collection = ImageFileCollection(
        filenames=[str(firstPath), str(secondPath)],
        location=str(tmp_path),
        keywords=["DPR_TYPE", "EXTRA"],
    )

    assert collection.summary["EXTRA"][0] == 5
    assert collection.summary["EXTRA"].mask.tolist() == [False, True]


def test_dict_from_fits_header_rejects_a_header_containing_a_file_keyword(
    tmp_path: Path,
) -> None:
    filePath = tmp_path / "one.fits"
    header = fits.Header({"FILE": "oops"})
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=header).writeto(filePath)

    # HEADER `in` MEMBERSHIP IS CASE-INSENSITIVE, SO EVEN AN UPPERCASE "FILE"
    # KEYWORD TRIPS THE `assert "file" not in h` GUARD
    with pytest.raises(AssertionError):
        ImageFileCollection(
            filenames=[str(filePath)], location=str(tmp_path), keywords=["FILE"]
        )


def test_dict_from_fits_header_real_duplicate_comment_card_warns_instead_of_joining(
    tmp_path: Path,
) -> None:
    """Pin that genuine duplicate COMMENT cards do NOT get joined.

    `astropy.io.fits.Header` always normalises the COMMENT/HISTORY keyword to
    uppercase, even when the card is built from a lowercase string, so the
    module's `k in ["comment", "history"]` check (comparing against a
    lowercase literal) never matches a key read from a real header. A second
    COMMENT card therefore falls through to the generic duplicate-keyword
    warning branch instead of being merged, and only the first COMMENT value
    survives.
    """
    filePath = tmp_path / "one.fits"
    header = fits.Header()
    header["DPR_TYPE"] = "BIAS"
    header.add_comment("first comment")
    header.add_comment("second comment")
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=header).writeto(filePath)

    with pytest.warns(UserWarning) as recordedWarnings:
        collection = ImageFileCollection(
            filenames=[str(filePath)],
            location=str(tmp_path),
            keywords=["DPR_TYPE", "COMMENT"],
        )

    assert list(collection.summary["COMMENT"]) == ["first comment"]
    assert 'contains multiple entries for "COMMENT"' in str(
        recordedWarnings[0].message
    )


def test_dict_from_fits_header_joins_comment_and_history_when_keys_are_lowercase(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pin the comment/history join logic in isolation.

    This is the only way to exercise it: a genuine `fits.getheader()` result
    never yields a lowercase "comment"/"history" key (see the test above), so
    `fits.getheader` is monkeypatched to return a header-like stand-in that
    does, to prove the join code is correct even though it is unreachable
    through the real FITS reading path.
    """

    class _LowercaseCommentHeader:
        def __init__(self) -> None:
            self._items = [
                ("DPR_TYPE", "BIAS"),
                ("comment", "first comment"),
                ("comment", "second comment"),
                ("history", "only history"),
            ]

        def __contains__(self, key: str) -> bool:
            return any(k.lower() == key.lower() for k, _ in self._items)

        def items(self) -> list[tuple[str, str]]:
            return self._items

    filePath = tmp_path / "one.fits"
    filePath.touch()
    monkeypatch.setattr(
        astropy_fits, "getheader", lambda *_args, **_kwargs: _LowercaseCommentHeader()
    )

    collection = ImageFileCollection(
        filenames=[str(filePath)],
        location=str(tmp_path),
        keywords=["DPR_TYPE", "comment", "history"],
    )

    assert list(collection.summary["comment"]) == ["first comment,second comment"]
    assert list(collection.summary["history"]) == ["only history"]
