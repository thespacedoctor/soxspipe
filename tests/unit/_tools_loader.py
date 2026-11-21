"""*load the developer scripts in tools/ by path, since tools/ is not an installed package*

A script imports its shared helpers by bare name (`import measurement_common`), which works when it runs from the
command line because its own folder is first on `sys.path`. The loader keeps every module under the prefixed name
`soxspipe_tools.<name>` in `sys.modules`, so no bare name stays behind after the test run starts. While it executes a
script, it points the bare names of the already-loaded helpers at those same module objects, so the script's own imports
resolve to them and the tests and the script share one copy of each helper.
"""

from __future__ import annotations

import importlib.util
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType

TOOLS_DIR = Path(__file__).resolve().parents[2] / "tools"
PREFIX = "soxspipe_tools."
SHARED_MODULES = ("measurement_common", "line_fitting", "spec_res_frames")


@contextmanager
def _bare_names() -> Iterator[None]:
    """*point the bare names of the loaded helpers at the prefixed modules, and restore the old entries on exit*"""
    saved = {name: sys.modules.get(name) for name in SHARED_MODULES}
    for name in SHARED_MODULES:
        if PREFIX + name in sys.modules:
            sys.modules[name] = sys.modules[PREFIX + name]
    try:
        yield
    finally:
        for name, module in saved.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


def load_tool(name: str) -> ModuleType:
    """*import a script or shared helper of tools/ by path, once*

    **Key Arguments:**

    - ``name`` -- the file name without ``.py``

    **Return:**

    - ``module`` -- the imported module. A module that is already loaded is returned, not run again. The shared helpers
      that come before ``name`` in ``SHARED_MODULES`` are loaded first.

    **Usage:**

    ```python
    step_tool = load_tool("measure_log_wavelength_step")
    ```
    """
    key = PREFIX + name
    if key in sys.modules:
        return sys.modules[key]
    # A HELPER NEEDS ONLY THE HELPERS BEFORE IT IN THE TUPLE, A SCRIPT NEEDS THEM ALL
    helperCount = SHARED_MODULES.index(name) if name in SHARED_MODULES else len(SHARED_MODULES)
    for helper in SHARED_MODULES[:helperCount]:
        load_tool(helper)
    path = TOOLS_DIR / f"{name}.py"
    specification = importlib.util.spec_from_file_location(key, path)
    if specification is None or specification.loader is None:
        raise AssertionError(f"{name} is not importable from {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[key] = module
    try:
        with _bare_names():
            specification.loader.exec_module(module)
    except BaseException:
        # LIKE IMPORTLIB, A MODULE THAT FAILED TO RUN IS NOT LEFT REGISTERED
        sys.modules.pop(key, None)
        raise
    return module
