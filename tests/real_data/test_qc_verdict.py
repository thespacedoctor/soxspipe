"""Acceptance check that no recipe in the real-data gate's reduction chain failed QC (DY-908).

The baseline checks read only the science products, so a calibration recipe that writes
its products but fails QC would pass them. This check reads every verdict in the
workspace's `quality_control` table. It reads the QC verdict, never the session status,
because a SOF whose QC failed can still end with session status `pass` (DY-272).

The verdicts are final when this runs. A recipe sets `qc_flag` from its QC-acceptable
ranges as it runs. The data organiser recomputes the flags from the same settings only
in `prep` and when it restores session statuses, and the gate runs neither after the
last reduction.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from tests.real_data.qc_verdict import AllowedQcFailure, qc_verdict_problems
from tests.real_data.reporting import report

pytestmark = pytest.mark.real_data

# EACH ENTRY NAMES THE RECIPE, SOF AND QC NAME OF ONE EXPECTED FAILURE, WITH THE REASON
# IT IS ACCEPTED. KEEP IT EMPTY UNLESS A FAILURE IS UNDERSTOOD AND TRACKED
ALLOWED_QC_FAILURES: tuple[AllowedQcFailure, ...] = ()


def test_no_recipe_in_the_reduction_chain_failed_qc(reduced_workspace: Path) -> None:
    # READ-ONLY, SO A MISSING DATABASE RAISES INSTEAD OF BEING CREATED EMPTY
    connection = sqlite3.connect(f"{(reduced_workspace / 'soxspipe.db').as_uri()}?mode=ro", uri=True)
    try:
        problems = qc_verdict_problems(connection, allowedFailures=ALLOWED_QC_FAILURES)
    finally:
        connection.close()
    report("qc verdict problems", len(problems))
    assert not problems, "QC failures recorded during the gate's reduction:\n" + "\n".join(problems)
