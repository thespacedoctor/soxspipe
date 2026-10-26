"""*guard for the load-bearing import order in the package ``__init__`` files*

``soxspipe.commonutils`` imports ``keyword_lookup`` and ``detector_lookup`` before everything else on purpose. If an
import sorter reorders the block, a submodule that runs ``from soxspipe.commonutils import keyword_lookup`` while the
package is half loaded binds the submodule instead of the class, and later calls fail with
``TypeError: 'module' object is not callable``. Each entry point is checked in a fresh interpreter because the
failure depends on which package is imported first.
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
SUBPROCESS_TIMEOUT_SECONDS = 120

CHECK_SCRIPT = textwrap.dedent(
    """
    import importlib
    import inspect
    import sys
    import types
    from pathlib import Path

    LOOKUP_NAMES = ("keyword_lookup", "detector_lookup")

    importlib.import_module(sys.argv[2])
    import soxspipe
    import soxspipe.commonutils as commonutils

    repoRoot = Path(sys.argv[1])
    packageFile = Path(soxspipe.__file__).resolve()
    isUnderRepo = any(parent.samefile(repoRoot) for parent in packageFile.parents)
    assert isUnderRepo, f"imported {packageFile}, not the package under {repoRoot}"

    for name in LOOKUP_NAMES:
        member = getattr(commonutils, name)
        assert inspect.isclass(member), f"soxspipe.commonutils.{name} is {member!r}, not a class"

    # A SUBMODULE THAT IMPORTED THE NAME TOO EARLY HOLDS THE SUBMODULE, WHICH THE PACKAGE ATTRIBUTE LATER HIDES
    misbound = sorted(
        f"{moduleName}.{name}"
        for moduleName, module in list(sys.modules.items())
        if moduleName.startswith("soxspipe") and module is not None
        for name in LOOKUP_NAMES
        if isinstance(vars(module).get(name), types.ModuleType)
    )
    assert not misbound, f"names bound to a module instead of a class: {misbound}"
    """
)


def _run_in_fresh_interpreter(firstModule: str) -> subprocess.CompletedProcess[str]:
    """*run the class check in a new interpreter that imports the repository's package*

    **Key Arguments:**

    - ``firstModule`` -- the module to import before anything else

    **Return:**

    - ``completed`` -- the finished process, with captured output
    """
    environment = dict(os.environ)
    existingPath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = os.pathsep.join([str(REPO_ROOT), existingPath] if existingPath else [str(REPO_ROOT)])

    return subprocess.run(  # noqa: S603
        [sys.executable, "-c", CHECK_SCRIPT, str(REPO_ROOT), firstModule],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT_SECONDS,
        check=False,
    )


@pytest.mark.parametrize("firstModule", ["soxspipe", "soxspipe.commonutils", "soxspipe.recipes"])
def test_lookup_names_are_classes_whichever_package_is_imported_first(firstModule):
    # ARRANGE: THE PACKAGE A FRESH INTERPRETER IMPORTS FIRST COMES FROM THE PARAMETRISATION

    # ACT
    completed = _run_in_fresh_interpreter(firstModule)

    # ASSERT
    assert completed.returncode == 0, f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
