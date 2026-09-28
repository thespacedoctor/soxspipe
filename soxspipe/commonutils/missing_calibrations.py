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


def describe_missing(positionsByType):
    """*join missing calibration types into one human-readable string*

    **Key Arguments:**

    - ``positionsByType`` -- dict of `{calibration type: position in the sof map}`

    **Return:**

    - ``description`` -- comma-separated labels in sof-map order, or `unknown` if the dict is empty
    """
    if not positionsByType:
        return UNKNOWN_MISSING
    ordered = sorted(positionsByType, key=lambda calibrationType: (positionsByType[calibrationType], calibrationType))
    return ", ".join(CALIBRATION_DESCRIPTIONS.get(calibrationType, calibrationType) for calibrationType in ordered)
