#!/usr/bin/env python
"""
*Quality-control measures for the spline fit of the instrument response*

:Author:
    David Young
"""

import os
from typing import Any

os.environ["TERM"] = "vt100"


def measure_response_qcs(fit: Any, bands: Any) -> dict:
    """*measure how well the spline describes the raw response*

    **Key Arguments:**

    - ``fit`` -- the ``ResponseFit``
    - ``bands`` -- list of ``[lower, upper]`` wavelength bands (nm) for the self-calibration residual

    **Return:**

    - ``measures`` -- dictionary with ``clipFrac`` (fraction of the non-wing bins clipped by the sigma clip; telluric
      wing bins are dropped by design and are left out of both the numerator and the denominator, and a fit with only
      wing bins gives 0.0), ``maxGap`` (nm, widest stretch inside the data range without a kept bin) and
      ``selfCal`` (list of ``(band, median(model / ratio) - 1)`` for every band that holds unmasked samples).
      ``model / ratio`` is the self-calibrated flux over the catalogue flux of each sample, so a positive residual
      means the calibrated standard is too bright. It is linear in the extracted counts, so negative low-signal
      counts do not bias the median, and the median ignores narrow telluric lines.
    """
    import numpy as np

    nonWing = ~fit.wingMask
    nonWingCount = np.count_nonzero(nonWing)
    clipFrac = float(np.count_nonzero(~fit.keptMask & nonWing) / nonWingCount) if nonWingCount else 0.0

    finiteWave = fit.wave[np.isfinite(fit.wave)]
    edges = np.concatenate([[finiteWave.min()], np.sort(fit.binWave[fit.keptMask]), [finiteWave.max()]])
    maxGap = float(np.max(np.diff(edges)))

    selfCal = []
    for lower, upper in bands:
        inBand = fit.mask & (fit.wave > lower) & (fit.wave < upper)
        if not inBand.any():
            continue
        residual = np.median(fit.evaluate(fit.wave[inBand]) / fit.ratio[inBand]) - 1.0
        selfCal.append(([lower, upper], float(residual)))
    return {"clipFrac": clipFrac, "maxGap": maxGap, "selfCal": selfCal}


def response_qc_rows(qcTable: Any, measures: dict, recipeName: str, dateObs: str) -> Any:
    """*append the response-fit measures to the QC table*

    **Key Arguments:**

    - ``qcTable`` -- the QC table to extend
    - ``measures`` -- the dictionary from ``measure_response_qcs``
    - ``recipeName`` -- name of the recipe to record on each row
    - ``dateObs`` -- observation date (ISO format)

    **Return:**

    - ``qcTable`` -- a new QC table with the measures added. The input table is not changed.
    """
    import pandas as pd

    from soxspipe.commonutils.toolkit import utcnow_string

    entries = [
        ("RESP BIN CLIP FRAC", measures["clipFrac"], "Fraction of response bins clipped from the spline fit"),
        ("RESP MAX GAP", measures["maxGap"], "Widest wavelength stretch (nm) inside the data without a kept bin"),
    ]
    for (lower, upper), value in measures["selfCal"]:
        entries.append(
            (
                f"RESP SELFCAL RES {lower}-{upper}",
                value,
                f"Median self-calibrated over catalogue flux minus one in {lower}-{upper} nm",
            )
        )

    utcNow = utcnow_string()
    rows = pd.DataFrame(
        [
            {
                "soxspipe_recipe": recipeName,
                "qc_name": name,
                "qc_value": float(f"{value:0.4f}"),
                "qc_comment": comment,
                "qc_unit": None,
                "obs_date_utc": dateObs,
                "reduction_date_utc": utcNow,
                "to_header": True,
            }
            for name, value, comment in entries
        ]
    )
    if qcTable.empty and len(qcTable.columns):
        # CONCATENATING ONTO AN EMPTY TABLE IS DEPRECATED IN PANDAS; KEEP ITS COLUMN SCHEMA INSTEAD
        return rows.reindex(columns=qcTable.columns)
    return pd.concat([qcTable, rows], ignore_index=True)


def response_qc_warnings(measures: dict, fitSettings: Any) -> list[str]:
    """*list the measures that lie outside the tolerances in the settings*

    **Key Arguments:**

    - ``measures`` -- the dictionary from ``measure_response_qcs``
    - ``fitSettings`` -- the arm settings holding ``clip_frac_max``, ``max_gap_nm_max`` and ``selfcal_res_max``

    **Return:**

    - ``messages`` -- one message per out-of-tolerance measure, empty when all are in tolerance
    """
    messages = []
    if measures["clipFrac"] > fitSettings["clip_frac_max"]:
        messages.append(
            f"{measures['clipFrac']:.0%} of the response bins were clipped (limit {fitSettings['clip_frac_max']:.0%})"
        )
    if measures["maxGap"] > fitSettings["max_gap_nm_max"]:
        messages.append(f"response data gap of {measures['maxGap']:.0f} nm exceeds {fitSettings['max_gap_nm_max']} nm")
    for (lower, upper), value in measures["selfCal"]:
        if abs(value) > fitSettings["selfcal_res_max"]:
            messages.append(
                f"response self-calibration residual {value:+.1%} in {lower}-{upper} nm exceeds "
                f"{fitSettings['selfcal_res_max']:.0%}"
            )
    return messages
