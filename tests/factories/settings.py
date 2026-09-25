"""Factories for isolated pipeline settings."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def pipeline_settings(
    workspacePath: Path,
    *,
    instrument: str = "soxs",
    overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a fresh minimal settings mapping rooted in a temporary workspace."""
    baseSettings: dict[str, Any] = {
        "instrument": instrument,
        "workspace-root-dir": str(workspacePath),
        "save-intermediate-products": False,
        "calibration-root-dir": str(workspacePath / "calibrations"),
    }
    return {**baseSettings, **(overrides or {})}
