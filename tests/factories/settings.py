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


def logging_settings(logPath: Path) -> dict[str, Any]:
    """Return a `logging settings` block shaped like the shipped default settings.

    The file handler uses the shipped ``w+`` mode, so a consumer that shares the log file must override it.
    """
    return {
        "logging settings": {
            "formatters": {
                "file_style": {"format": "%(levelname)s: %(message)s", "datefmt": "%H:%M:%S"},
                "console_style": {"format": "%(levelname)s: %(message)s", "datefmt": "%H:%M:%S"},
            },
            "handlers": {
                "console": {
                    "class": "logging.StreamHandler",
                    "level": "WARNING",
                    "formatter": "console_style",
                    "stream": "ext://sys.stdout",
                },
                "file": {
                    "class": "logging.handlers.GroupWriteRotatingFileHandler",
                    "level": "INFO",
                    "formatter": "file_style",
                    "filename": str(logPath),
                    "mode": "w+",
                    "maxBytes": 102400,
                    "backupCount": 1,
                    "delay": True,
                },
            },
            "root": {"level": "INFO", "handlers": ["file", "console"]},
        }
    }
