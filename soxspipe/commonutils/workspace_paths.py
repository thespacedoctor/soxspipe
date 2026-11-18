#!/usr/bin/env python
"""
*Join a path under a workspace root*

``data_organiser``, ``set_of_files`` and ``toolkit`` all hold workspace-relative paths such as
``./raw/2024-01-01/a.fits``. This module gives them one rule for turning such a path into a path that does not
depend on the current working directory.

This module imports nothing from soxspipe, so any module can import it without creating an import cycle.

Author
: David Young

Date Created
: November 18, 2026
"""

import os

CURRENT_DIRECTORY_ROOTS = (".", "./", "")


def path_under_root(path: str, rootDir: str) -> str:
    """*return a path joined under the workspace root*

    **Key Arguments:**

    - ``path`` -- a file path, normally relative to the workspace
    - ``rootDir`` -- the workspace root. A leading ``~`` is expanded. ``"."``, ``"./"`` and ``""`` mean the current
      working directory

    **Return:**

    - ``rootedPath`` -- ``path`` unchanged when it is absolute or ``rootDir`` is the current working directory,
      else the tidied ``path`` joined to ``rootDir``

    **Usage:**

    ```python
    from soxspipe.commonutils.workspace_paths import path_under_root
    path_under_root("./raw/2024-01-01/a.fits", "/data/workspace")
    # "/data/workspace/raw/2024-01-01/a.fits"
    ```
    """
    if rootDir in CURRENT_DIRECTORY_ROOTS or os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(os.path.expanduser(rootDir), path))
