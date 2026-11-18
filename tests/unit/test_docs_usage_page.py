"""*tests that the tracked docs usage page matches what the docs build generates*

The docs build rewrites `docs/source/usage.md` from the `cl_utils` docstring. The tracked file must
equal that output, or every docs build leaves the work tree dirty.
"""

from __future__ import annotations

import ast
import os
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_SOURCE = REPO_ROOT / "docs" / "source"
USAGE_FUNCTION_NAME = "updateUsageMd"


def _generate_usage_page(configDirectory: Path) -> str:
    """*run the docs build's `updateUsageMd()` against a scratch directory*

    The function is taken from `conf.py` without importing the module, because importing it also runs the
    rest of the Sphinx configuration and rewrites the tracked file.

    **Key Arguments:**

    - ``configDirectory`` -- scratch directory that `updateUsageMd()` writes `usage.md` into

    **Return:**

    - ``usagePage`` -- the text `updateUsageMd()` wrote
    """
    confSource = (DOCS_SOURCE / "conf.py").read_text(encoding="utf-8")
    functionNodes = [
        node
        for node in ast.parse(confSource).body
        if isinstance(node, ast.FunctionDef) and node.name == USAGE_FUNCTION_NAME
    ]
    assert len(functionNodes) == 1, f"conf.py must define {USAGE_FUNCTION_NAME}() exactly once"

    namespace = {"__file__": str(configDirectory / "conf.py"), "os": os}
    exec(compile(ast.Module(body=functionNodes, type_ignores=[]), "conf.py", "exec"), namespace)  # noqa: S102
    namespace[USAGE_FUNCTION_NAME]()

    usagePath = configDirectory / "usage.md"
    assert usagePath.exists(), f"{USAGE_FUNCTION_NAME}() wrote no usage.md; check the cl_utils docstring"
    return usagePath.read_text(encoding="utf-8")


def test_tracked_usage_page_matches_the_generated_page(tmp_path: Path) -> None:
    # ARRANGE
    trackedPage = (DOCS_SOURCE / "usage.md").read_text(encoding="utf-8")

    # ACT
    generatedPage = _generate_usage_page(tmp_path)

    # ASSERT
    assert trackedPage == generatedPage, (
        "docs/source/usage.md is out of date; run updateUsageMd() (or build the docs) and commit usage.md"
    )
