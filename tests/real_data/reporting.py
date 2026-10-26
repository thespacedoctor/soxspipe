"""Baseline reporting shared by the real-data acceptance checks."""

from __future__ import annotations


def report(label: str, value) -> None:
    """Print a measured baseline quantity so the CI log records what this run produced.

    The bands in each check are only as good as the spread they were recorded from, so every
    run leaves its own numbers in the log for the next re-record.
    """
    print(f"BASELINE {label} = {value!r}", flush=True)
