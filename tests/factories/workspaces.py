"""Factories for isolated data-organiser workspaces."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from soxspipe.commonutils import data_organiser


@dataclass(frozen=True)
class WorkspaceLayout:
    """Paths belonging to one disposable reduction workspace."""

    root: Path
    raw: Path
    misc: Path
    sessions: Path


def workspace_layout(destination: Path) -> WorkspaceLayout:
    """Create and return a fresh, isolated workspace directory layout."""
    rootPath = destination / "workspace"
    rawPath = rootPath / "raw"
    miscPath = rootPath / "misc"
    sessionsPath = rootPath / "sessions"
    for path in (rawPath, miscPath, sessionsPath):
        path.mkdir(parents=True, exist_ok=True)
    return WorkspaceLayout(
        root=rootPath,
        raw=rawPath,
        misc=miscPath,
        sessions=sessionsPath,
    )


def workspace_organiser(destination: Path, *, log: Any) -> data_organiser:
    """Return an organizer rooted in a fresh workspace without opening SQLite."""
    layout = workspace_layout(destination)
    return data_organiser(log=log, rootDir=str(layout.root), dbConnect=False)


def session_workspace(destination: Path, sessionId: str) -> Path:
    """Create a workspace whose current reduction session is `sessionId` and return its root."""
    rootPath = destination / "workspace"
    (rootPath / "sessions" / sessionId).mkdir(parents=True)
    (rootPath / "sessions" / ".sessionid").write_text(sessionId, encoding="utf-8")
    return rootPath
