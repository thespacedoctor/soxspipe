#!/usr/bin/env python
"""*fail when a changelog bullet is repeated within one release block*

A release block starts at each `## ` heading. The lines before the first `## ` heading
form a block of their own, which is where unreleased entries live. A bullet is a line
that starts with `* ` in the first column, the only list style `CHANGES.md` uses; `- `
and `+ ` bullets, indented bullets and fenced code blocks are not recognised. Two
bullets match when their text is equal once trailing whitespace is stripped.

The same bullet in two different release blocks is allowed, since a fix can legitimately
be noted again in a later release. Only a repeat inside one block fails.

**Usage:**

```bash
python tools/check_changelog.py
python tools/check_changelog.py CHANGES.md
```

Each duplicate is printed with its block heading and the line numbers of both copies.
Exit status is 0 when the file is clean, 1 when duplicates are found, and 2 when the
file could not be read.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import NamedTuple

DEFAULT_CHANGELOG = Path("CHANGES.md")

# THE HEADING THAT STARTS A RELEASE BLOCK
RELEASE_HEADING_PREFIX = "## "

# THE MARKER THAT STARTS A BULLET
BULLET_PREFIX = "* "

# THE LABEL FOR THE LINES BEFORE THE FIRST RELEASE HEADING, WHERE UNRELEASED ENTRIES LIVE
UNRELEASED_HEADING = "(unreleased, before the first ## heading)"


class Duplicate(NamedTuple):
    """*one bullet repeated within a single release block*"""

    heading: str
    text: str
    firstLine: int
    secondLine: int

    def __str__(self) -> str:
        return f"{self.heading}: lines {self.firstLine} and {self.secondLine}: {self.text}"


def find_duplicates(text: str) -> list[Duplicate]:
    """*find every bullet repeated within one release block*

    A bullet that appears more than twice is reported once per repeat, each against
    the line of its first occurrence.

    **Key Arguments:**

    - ``text`` -- the changelog text

    **Return:**

    - ``duplicates`` -- the repeated bullets, in file order
    """
    duplicates = []
    heading = UNRELEASED_HEADING
    firstSeen: dict[str, int] = {}

    for lineNumber, line in enumerate(text.splitlines(), start=1):
        if line.startswith(RELEASE_HEADING_PREFIX):
            heading = line.rstrip()
            firstSeen = {}
        elif line.startswith(BULLET_PREFIX):
            bullet = line.rstrip()
            if bullet in firstSeen:
                duplicates.append(Duplicate(heading, bullet, firstSeen[bullet], lineNumber))
            else:
                firstSeen[bullet] = lineNumber

    return duplicates


def main(argv: list[str] | None = None) -> int:
    """*run the checker from the command line*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments, excluding the program name. Default *None*, i.e. `sys.argv[1:]`

    **Return:**

    - ``status`` -- 0 when the file is clean, 1 when duplicates were found, 2 when the file could not be read
    """
    parser = argparse.ArgumentParser(description="Reject a changelog bullet repeated within one release block.")
    parser.add_argument(
        "path", nargs="?", type=Path, default=DEFAULT_CHANGELOG, help="the changelog to check. Default CHANGES.md"
    )
    arguments = parser.parse_args(argv)

    try:
        # utf-8-sig DROPS A LEADING BOM, WHICH WOULD OTHERWISE HIDE A ## HEADING ON LINE 1
        text = arguments.path.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError) as error:
        print(f"{arguments.path}: could not read the changelog: {error}", file=sys.stderr)
        return 2

    duplicates = find_duplicates(text)
    for duplicate in duplicates:
        print(duplicate)

    return 1 if duplicates else 0


if __name__ == "__main__":
    sys.exit(main())
