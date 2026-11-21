"""*load the developer scripts in tools/ by path, since tools/ is not an installed package*

The scripts import their shared helpers as plain top-level modules (`import measurement_common`), which works when a
script runs from the command line because its own folder is first on `sys.path`. The loader registers those helpers in
`sys.modules` before it loads a script, so a script loaded by path finds the same helper module the tests hold.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

TOOLS_DIR = Path(__file__).resolve().parents[2] / "tools"
SHARED_MODULES = ("measurement_common", "line_fitting", "spec_res_frames")


def _load(name: str) -> ModuleType:
    """*import ``tools/<name>.py`` and register it in ``sys.modules`` under its plain name*"""
    path = TOOLS_DIR / f"{name}.py"
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise AssertionError(f"{name} is not importable from {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def load_tool(name: str) -> ModuleType:
    """*import a script or shared helper of tools/ by path*

    **Key Arguments:**

    - ``name`` -- the file name without ``.py``

    **Return:**

    - ``module`` -- the imported module. The shared helpers are loaded once and reused.

    **Usage:**

    ```python
    step_tool = load_tool("measure_log_wavelength_step")
    ```
    """
    for shared in SHARED_MODULES:
        if shared not in sys.modules and (TOOLS_DIR / f"{shared}.py").exists():
            _load(shared)
    if name in SHARED_MODULES:
        return sys.modules[name]
    return _load(name)
