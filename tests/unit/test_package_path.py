"""Package-location contract for bundled resources."""

from __future__ import annotations

from pathlib import Path

import pytest

from soxspipe.commonutils.getpackagepath import getpackagepath

pytestmark = pytest.mark.unit


def test_getpackagepath_returns_the_installed_package_root() -> None:
    """Resource callers receive the directory that contains ``commonutils``."""
    packagePath = Path(getpackagepath()).resolve()

    assert packagePath.name == "soxspipe"
    assert (packagePath / "commonutils").is_dir()
    assert (packagePath / "resources").is_dir()
