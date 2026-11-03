#!/usr/bin/env python
"""
*make, track and remove the scratch directories that recipes keep in* ``<workspace>/tmp/``

Every recipe makes one ``mkdtemp`` directory in ``<workspace>/tmp/`` and writes its prepared frames there.
This module records those directories, and it holds the one deletion policy for them: a directory goes only
when its resolved parent is the scratch root, and a failure to remove one is logged and never raised.

The module imports nothing from the package, so ``base_recipe`` and ``reducer`` can both use it without
disturbing the package's load-bearing import order.

Author
: David Young

Date Created
: November 3, 2026
"""

import os
import shutil
import tempfile
import threading

# THE WORKSPACE ROOT THAT `advanced_settings.yaml` SETS WHEN A SETTINGS FILE HAS NO `workspace-root-dir`
DEFAULT_WORKSPACE_ROOT = "./"

# `mkdtemp` NAMES ITS DIRECTORIES WITH THIS PREFIX, SO THE SWEEP LEAVES EVERYTHING ELSE IN `tmp/` ALONE
SCRATCH_PREFIX = tempfile.template

# THE SCRATCH DIRECTORIES EACH THREAD HAS CREATED, IN ORDER. IT IS PER THREAD SO THAT A CALLER THAT CLEANS UP
# AFTER A FAILED RECIPE CAN NEVER REACH THE LIVE SCRATCH OF A RECIPE RUNNING IN ANOTHER THREAD
_REGISTRY = threading.local()


def absolute_path(path: str) -> str:
    """*convert a path from a home directory to an absolute path*

    **Key Arguments:**

    - ``path`` -- path possibly relative to the home directory

    **Return:**

    - ``absolutePath`` -- the path with a leading ``~`` expanded
    """
    home = os.path.expanduser("~")
    if path[0] == "~":
        path = home + "/" + path[1:]

    return path.replace("//", "/")


def scratch_root(settings: dict) -> str:
    """*return the directory in which recipes make their scratch directories*

    **Key Arguments:**

    - ``settings`` -- the settings dictionary. A missing ``workspace-root-dir`` reads as the shipped default

    **Return:**

    - ``scratchRoot`` -- ``<workspace-root-dir>/tmp`` as an absolute path, so a later change of working
      directory cannot move it
    """
    workspaceRoot = absolute_path(settings.get("workspace-root-dir", DEFAULT_WORKSPACE_ROOT))
    return os.path.abspath(os.path.join(workspaceRoot, "tmp"))


def created_scratch_directories() -> list[str]:
    """*return the list of scratch directories the calling thread has created and not yet removed*

    **Return:**

    - ``directories`` -- the live list. A caller can mark its length and later take the entries added after it.
    """
    if not hasattr(_REGISTRY, "directories"):
        _REGISTRY.directories = []
    return _REGISTRY.directories


def register_scratch_directory(path: str) -> None:
    """*record a scratch directory the calling thread has just created*

    **Key Arguments:**

    - ``path`` -- the new scratch directory
    """
    created_scratch_directories().append(path)


def unregister_scratch_directory(path: str) -> None:
    """*forget a scratch directory the calling thread has removed; one never recorded is ignored*

    **Key Arguments:**

    - ``path`` -- the scratch directory
    """
    directories = created_scratch_directories()
    if path in directories:
        directories.remove(path)


def remove_scratch_tree(log, path: str | None, requiredParent: str | None = None) -> None:
    """*remove a scratch directory tree, logging a failure instead of raising it*

    A missing directory is not a failure. The function is safe to call from a ``finally`` block, because
    it never replaces the exception that is already on its way out, whatever the type of its own failure.

    **Key Arguments:**

    - ``log`` -- the logger to warn on
    - ``path`` -- the directory to remove. An empty value or *None* does nothing.
    - ``requiredParent`` -- when given, ``path`` is removed only if its resolved parent directory is this
      resolved directory; a symbolic link is judged by where it leads. Otherwise a warning is logged and
      nothing is removed. Default *None* (no check).
    """
    if not path:
        return
    try:
        if requiredParent is not None and os.path.dirname(os.path.realpath(path)) != os.path.realpath(requiredParent):
            log.warning(f"Not removing `{path}`: it is not directly inside the scratch directory `{requiredParent}`.")
            return
        shutil.rmtree(path)
    except FileNotFoundError:
        return
    except Exception as e:
        log.warning(f"Could not remove the scratch directory `{path}`: {e}", exc_info=True)


def remove_recipe_scratch(log, settings: dict, scratchDirectories: list[str]) -> None:
    """*remove the scratch directories a recipe run created, each only if it is inside* ``<workspace>/tmp``

    Never raises, so it is safe in a ``finally`` block.

    **Key Arguments:**

    - ``log`` -- the logger to warn on
    - ``settings`` -- the settings dictionary the recipe ran with. ``workspace-root-dir`` names the workspace,
      read the way ``base_recipe`` reads it
    - ``scratchDirectories`` -- the directories the run created
    """
    if not scratchDirectories:
        return
    try:
        scratchPath = scratch_root(settings)
    except Exception as e:
        log.warning(f"Could not find the scratch root, so {scratchDirectories} stay: {e}", exc_info=True)
        return
    for path in scratchDirectories:
        remove_scratch_tree(log, path, requiredParent=scratchPath)


def clear_workspace_scratch(log, workspaceDirectory: str) -> None:
    """*remove the* ``mkdtemp`` *directories in* ``<workspace>/tmp/``, *where recipes keep their prepared frames*

    Only one soxspipe process runs in a workspace at a time (DY-1129), so every such directory is left
    over from a recipe that was killed or failed, or from a bulk run that has finished. ``tmp/`` itself
    and anything in it that ``mkdtemp`` did not make stay in place.

    Nothing is removed, and a warning is logged, when the workspace path is empty, is not a directory, or
    is the filesystem root or the home directory. A ``tmp`` that is a symbolic link is deliberately not
    swept (a warning only), because the sweep must not reach outside the workspace.

    **Key Arguments:**

    - ``log`` -- the logger to warn on
    - ``workspaceDirectory`` -- path to the root of the workspace
    """
    if not workspaceDirectory or not os.path.isdir(workspaceDirectory):
        log.warning(f"Not sweeping scratch: `{workspaceDirectory}` is not a workspace directory.")
        return
    if os.path.realpath(workspaceDirectory) in (os.path.sep, os.path.realpath(os.path.expanduser("~"))):
        log.warning(f"Not sweeping scratch: `{workspaceDirectory}` is the filesystem root or the home directory.")
        return
    scratchPath = os.path.join(workspaceDirectory, "tmp")
    if os.path.islink(scratchPath):
        log.warning(f"Not sweeping scratch: `{scratchPath}` is a symbolic link.")
        return
    if not os.path.isdir(scratchPath):
        return
    try:
        names = sorted(os.listdir(scratchPath))
    except OSError as e:
        log.warning(f"Could not list the scratch directory `{scratchPath}`: {e}")
        return
    for name in names:
        child = os.path.join(scratchPath, name)
        if name.startswith(SCRATCH_PREFIX) and (os.path.islink(child) or os.path.isdir(child)):
            remove_scratch_tree(log, child, requiredParent=scratchPath)
