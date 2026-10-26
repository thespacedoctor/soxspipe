"""The chain-wide QC-verdict check of the real-data gate (DY-908).

A recipe can write every product file and still fail QC, so a product-level baseline
alone can pass on a QC-failed reduction. This check reads the verdict the recipes
recorded in the workspace's `quality_control` table instead.
"""

from __future__ import annotations

import sqlite3
from dataclasses import dataclass


@dataclass(frozen=True)
class AllowedQcFailure:
    """One QC failure the gate expects, and why it is accepted."""

    recipe: str
    sofName: str
    qcName: str
    reason: str


def qc_verdict_problems(connection: sqlite3.Connection, allowedFailures: tuple[AllowedQcFailure, ...]) -> list[str]:
    """Return one message per QC row whose verdict is `fail` and that is not allowed.

    An empty list means the check passes. A missing or empty table is a problem too, so a
    reduction that recorded no verdict at all cannot pass as one with no failures.
    """
    hasTable = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'quality_control'"
    ).fetchone()
    if hasTable is None:
        return ["the workspace database has no quality_control table"]
    (rowCount,) = connection.execute("SELECT COUNT(*) FROM quality_control").fetchone()
    if rowCount == 0:
        return ["the quality_control table holds no rows, so no QC verdict was recorded"]

    failingRows = connection.execute(
        "SELECT soxspipe_recipe, sof_name, qc_name, qc_value, qc_value_min, qc_value_max "
        "FROM quality_control WHERE qc_flag = 'fail' "
        "ORDER BY soxspipe_recipe, sof_name, qc_name"
    ).fetchall()
    allowedKeys = {(allowed.recipe, allowed.sofName, allowed.qcName) for allowed in allowedFailures}
    return [
        f"{recipe} {sofName}: {qcName} = {qcValue} failed QC (acceptable range {qcMin} to {qcMax})"
        for recipe, sofName, qcName, qcValue, qcMin, qcMax in failingRows
        if (recipe, sofName, qcName) not in allowedKeys
    ]
