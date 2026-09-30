#!/usr/bin/env python
"""
*Name the calibrations that stop a science set from being reduced*

Author
: David Young

Date Created
: September 28, 2026
"""

from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

CALIBRATION_DESCRIPTIONS = {
    "mbias": "master bias",
    "mdark": "master dark",
    "disp_solution": "dispersion solution",
    "order_centres": "order centres",
    "mflat": "master flat",
    "spat_solution": "spatial solution",
    "std_flux": "flux standard",
}

UNKNOWN_MISSING = "unknown"

# WHICH RECIPE(S) PRODUCE EACH CALIBRATION TYPE, FOR TRACING WHY ONE IS MISSING
CALIBRATION_RECIPES = {
    "mbias": ("mbias",),
    "mdark": ("mdark",),
    "disp_solution": ("disp_solution",),
    "order_centres": ("order_centres",),
    "mflat": ("mflat",),
    "spat_solution": ("spat_solution",),
    "std_flux": ("stare_std", "nod_std", "offset_std"),
}

FAILED_QC = "failed QC"
FAILED_RUN = "failed run"
NOT_YET_REDUCED = "not yet reduced"
NOT_OBSERVED = "not observed"
NO_MATCH = "no match"

# `run_recipe_bulk` (soxspipe/commonutils/reducer.py) WRITES THIS EXACT PREFIX ONLY WHEN A RECIPE COMPLETED
# BUT ONE OF ITS QC VALUES WAS OUT OF RANGE; ANY OTHER `status = 'fail'` ROW (AN EXCEPTION OR CRASH) CARRIES
# THE RAISED EXCEPTION'S OWN TEXT INSTEAD, SO THIS PREFIX IS HOW classify_missing_reason TELLS THE TWO APART
QC_FAILURE_MESSAGE_PREFIX = "QC values outside of acceptable limits:"

# CALIBRATION TYPES WHOSE `cal_<type>` VIEW HARD-GATES A CANDIDATE UPSTREAM PRODUCT ON READOUT SPEED
# (BOTH ARMS) AND, ON VIS ONLY, ON BINNING AND GAIN TOO (SEE soxs_sof_map.yaml-DRIVEN VIEWS IN THE DB:
# cal_mbias/cal_mflat REQUIRE `upstream.rospeed = downstream.rospeed`, AND ON VIS ALSO
# `upstream.binning = downstream.binning AND upstream.gain = downstream.gain`; NIR IS EXEMPT FROM THE
# BINNING/GAIN CHECK. EVERY OTHER TYPE'S VIEW RANKS BY BINNING/SLIT/ROSPEED PROXIMITY BUT DOES NOT
# REQUIRE THEM TO MATCH, SO classify_missing_reason LEAVES THOSE TYPES ARM-ONLY TO STAY FAITHFUL TO THE
# REAL VIEW, RATHER THAN GUESSING A TIGHTER MATCH THAT THE VIEW ITSELF DOES NOT APPLY
ROSPEED_GATED_TYPES = frozenset({"mbias", "mflat"})
VIS_BINNING_GATED_TYPES = frozenset({"mbias", "mflat"})


def required_calibrations(sofMapLookup, recipe, arm):
    """*list the calibration types a recipe needs on one arm, in sof-map order*

    **Key Arguments:**

    - ``sofMapLookup`` -- the parsed `sof_map.yaml`, or `None` if it was never loaded
    - ``recipe`` -- the recipe name to look up
    - ``arm`` -- the instrument arm the calibrations are keyed by

    **Return:**

    - ``calibrationTypes`` -- the type names, or an empty list if the recipe, arm or `calibrations` entry is absent
    """
    for filters in (sofMapLookup or {}).values():
        if filters.get("recipe") != recipe:
            continue
        calibrations = filters.get("calibrations")
        if not isinstance(calibrations, dict):
            return []
        return list(calibrations.get(arm) or [])
    return []


def find_missing_calibrations(conn, incompleteSets, sofMapLookup):
    """*find which required calibrations each incomplete science set lacks*

    A calibration is missing when its `cal_<type>` view has no row for the
    set's sof with a passing (or null) upstream status. This is the same test
    `build_sof_files` uses to mark a set complete.

    **Key Arguments:**

    - ``conn`` -- the open SQLite connection holding the `cal_<type>` views
    - ``incompleteSets`` -- dataframe with `sof`, `recipe` and `eso seq arm` columns (not modified)
    - ``sofMapLookup`` -- the instrument's parsed `sof_map.yaml`, or `None`

    **Return:**

    - ``missingBySof`` -- dict mapping every sof to a dict of `{missing type: position in the sof map}`
    """
    uniqueSets = incompleteSets[["sof", "recipe", "eso seq arm"]].drop_duplicates()
    requiredBySof = {
        sof: required_calibrations(sofMapLookup, recipe, arm) for sof, recipe, arm in uniqueSets.itertuples(index=False)
    }

    sofsByType = {}
    for sof, calibrationTypes in requiredBySof.items():
        for calibrationType in calibrationTypes:
            sofsByType.setdefault(calibrationType, []).append(sof)

    missingBySof = {sof: {} for sof in requiredBySof}
    for calibrationType, sofs in sofsByType.items():
        calTable = validate_sql_identifier(f"cal_{calibrationType}", "calibration table name")
        # IDENTIFIER VALIDATED ABOVE AND NO VALUES ARE INTERPOLATED
        sqlQuery = (
            f"SELECT DISTINCT sof FROM {calTable} "  # noqa: S608
            "WHERE upstream_status = 'pass' OR upstream_status IS NULL"
        )
        presentSofs = {row[0] for row in conn.execute(sqlQuery)}
        for sof in sofs:
            if sof not in presentSofs:
                missingBySof[sof][calibrationType] = requiredBySof[sof].index(calibrationType)
    return missingBySof


def describe_missing(positionsByType, reasonsByType=None):
    """*join missing calibration types into one human-readable string*

    **Key Arguments:**

    - ``positionsByType`` -- dict of `{calibration type: position in the sof map}`
    - ``reasonsByType`` -- optional dict of `{calibration type: reason string}`. When a type has a reason,
      it is appended to that type's label in parentheses.

    **Return:**

    - ``description`` -- comma-separated labels in sof-map order, or `unknown` if the dict is empty
    """
    if not positionsByType:
        return UNKNOWN_MISSING
    reasonsByType = reasonsByType or {}
    ordered = sorted(positionsByType, key=lambda calibrationType: (positionsByType[calibrationType], calibrationType))
    labels = []
    for calibrationType in ordered:
        label = CALIBRATION_DESCRIPTIONS.get(calibrationType, calibrationType)
        reason = reasonsByType.get(calibrationType)
        if reason:
            label = f"{label} ({reason})"
        labels.append(label)
    return ", ".join(labels)


def calibration_match_clause(calibrationType, arm, binning, rospeed, gain):
    """*build the extra `WHERE` clause and params that mirror one `cal_<type>` view's hard-match rules*

    Only `mbias` and `mflat` gate on readout speed (both arms) and, on VIS, binning and gain too; see
    `ROSPEED_GATED_TYPES`/`VIS_BINNING_GATED_TYPES`. A `None` context value (not known for this group)
    is never gated on, since a false exclusion is worse than the coarser arm-only match it replaces.

    **Return:**

    - ``clause``, ``params`` -- an extra SQL fragment (starting with a leading space, or `""`) and its
      bound parameters, to append after the existing `recipe IN (...) AND "eso seq arm" = ?` clause
    """
    clause = ""
    params = []
    if calibrationType in ROSPEED_GATED_TYPES and rospeed is not None:
        clause += " AND rospeed = ?"
        params.append(rospeed)
    if calibrationType in VIS_BINNING_GATED_TYPES and arm == "VIS":
        if binning is not None:
            clause += " AND binning = ?"
            params.append(binning)
        if gain is not None:
            clause += " AND gain = ?"
            params.append(gain)
    return clause, params


def classify_missing_reason(conn, calibrationType, arm, binning=None, rospeed=None, gain=None):
    """*work out why one calibration type is missing for one arm*

    Checked in order: a producing recipe's product failed QC; a producing recipe's product failed for
    some other reason (an exception, not a QC threshold); a producing recipe's set has not finished
    reducing yet; no raw frames for a producing recipe exist at all on this arm; otherwise frames exist
    and presumably passed, but no `cal_<type>` row matched (binning, readout speed or slit mismatch).

    For `mbias`/`mflat`, ``binning``, ``rospeed`` and ``gain`` narrow every check to the same hard-match
    rules their `cal_<type>` view applies (see `calibration_match_clause`), so an unrelated failed or
    unfinished run in a different binning/readout mode on the same arm is not misattributed to this
    calibration. Every other type's view does not hard-gate on these, so they are ignored for it.

    **Key Arguments:**

    - ``conn`` -- the open SQLite connection holding `failed_products` and `raw_frame_sets`
    - ``calibrationType`` -- the calibration type to classify
    - ``arm`` -- the instrument arm to check
    - ``binning`` -- the science set's binning, or `None` if not known
    - ``rospeed`` -- the science set's readout speed, or `None` if not known
    - ``gain`` -- the science set's gain, or `None` if not known

    **Return:**

    - ``reason`` -- one of `FAILED_QC`, `FAILED_RUN`, `NOT_YET_REDUCED`, `NO_MATCH`, `NOT_OBSERVED`
    """
    recipes = CALIBRATION_RECIPES.get(calibrationType, (calibrationType,))
    placeholders = ", ".join("?" for _ in recipes)
    matchClause, matchParams = calibration_match_clause(calibrationType, arm, binning, rospeed, gain)

    sqlQuery = (
        f"SELECT error_message FROM failed_products WHERE recipe IN ({placeholders}) "  # noqa: S608
        f'AND "eso seq arm" = ?{matchClause}'
    )
    # RECIPE NAMES COME FROM CALIBRATION_RECIPES (A MODULE CONSTANT), NOT USER INPUT; ALL VALUES ARE BOUND PARAMETERS
    failures = conn.execute(sqlQuery, (*recipes, arm, *matchParams)).fetchall()
    if failures:
        isQcFailure = any((message or "").startswith(QC_FAILURE_MESSAGE_PREFIX) for (message,) in failures)
        return FAILED_QC if isQcFailure else FAILED_RUN

    sqlQuery = (
        f"SELECT 1 FROM raw_frame_sets WHERE recipe IN ({placeholders}) "  # noqa: S608
        f'AND "eso seq arm" = ? AND complete = 0{matchClause} LIMIT 1'
    )
    if conn.execute(sqlQuery, (*recipes, arm, *matchParams)).fetchone():
        return NOT_YET_REDUCED

    sqlQuery = (
        f"SELECT 1 FROM raw_frame_sets WHERE recipe IN ({placeholders}) "  # noqa: S608
        f'AND "eso seq arm" = ?{matchClause} LIMIT 1'
    )
    if conn.execute(sqlQuery, (*recipes, arm, *matchParams)).fetchone():
        return NO_MATCH

    return NOT_OBSERVED


def missing_reasons(conn, calibrationTypes, arm, binning=None, rospeed=None, gain=None):
    """*classify why each of several calibration types is missing, for one arm*

    **Key Arguments:**

    - ``conn`` -- the open SQLite connection holding `failed_products` and `raw_frame_sets`
    - ``calibrationTypes`` -- iterable of calibration types to classify
    - ``arm`` -- the instrument arm to check
    - ``binning`` -- the science set's binning, or `None` if not known
    - ``rospeed`` -- the science set's readout speed, or `None` if not known
    - ``gain`` -- the science set's gain, or `None` if not known

    **Return:**

    - ``reasonsByType`` -- dict of `{calibration type: reason string}`
    """
    return {
        calibrationType: classify_missing_reason(conn, calibrationType, arm, binning, rospeed, gain)
        for calibrationType in calibrationTypes
    }
