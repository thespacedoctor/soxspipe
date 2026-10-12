#!/usr/bin/env python
"""
*Match each incomplete product SOF with its nearest upstream calibration products*

The `cal_<type>` views in `soxspipe.db` read the `calibration_matches` table that this module fills. Each
calibration type has one `CalibrationRule`: the hard matches an upstream product must pass, and how many of the
nearest passing products (by `|Δ mjd-obs|`) a downstream SOF keeps.

Author
: David Young

Date Created
: October 12, 2026
"""

from __future__ import annotations

import contextlib
import re
import sqlite3
from dataclasses import dataclass

import numpy as np
import pandas as pd

PASSING_STATUS = "pass"
NIR_ARM = "NIR"

# THE `product_frames` COLUMNS THE MATCHER READS, FOR BOTH UPSTREAM PRODUCTS AND DOWNSTREAM SOFS
PRODUCT_COLUMNS = [
    "file",
    "sof",
    "eso pro catg",
    "filepath",
    "recipe",
    "recipe_order",
    "status",
    "complete",
    "eso seq arm",
    "mjd-obs",
    "exptime",
    "binning",
    "gain",
    "rospeed",
    "slit",
]
_SELECT_COLUMNS = ", ".join(f'"{column}"' for column in PRODUCT_COLUMNS)


@dataclass(frozen=True)
class CalibrationRule:
    """*the matching rule for one calibration type*

    Every rule also requires the same arm, an upstream `recipe_order` lower than the downstream one, an upstream
    product with `complete = 1` and status `pass` or NULL, and a downstream SOF with `complete < 1`. As in SQL, a
    NULL value never matches anything.

    **Key Arguments:**

    - ``calType`` -- the calibration type name (`cal_<calType>` is its view)
    - ``upstreamRecipes`` -- the recipes whose products can serve as this calibration
    - ``nearestCount`` -- how many of the nearest upstream product rows each downstream SOF keeps
    - ``matchRospeed`` -- upstream readout speed must equal the downstream one
    - ``matchBinningGainUnlessNir`` -- upstream binning and gain must equal the downstream ones, except on NIR
    - ``matchSlit`` -- upstream slit must equal the downstream one
    - ``waivedForDownstreamRecipes`` -- downstream recipes for which `matchRospeed` and `matchSlit` are waived
    - ``nirSlitFallback`` -- `(upstream slit, downstream slit)`: on NIR this upstream slit also serves this
      downstream slit
    - ``upstreamArms`` -- when set, only upstream products from these arms are used
    - ``upstreamSlit`` -- when set, only upstream products with this slit are used
    - ``upstreamTagLike`` -- when set, only upstream products whose `eso pro catg` is SQL `LIKE` this pattern
    - ``downstreamRecipes`` -- when set, only downstream SOFs of these recipes get a match
    - ``excludedDownstreamRecipes`` -- downstream SOFs of these recipes never get a match
    - ``rankByExptime`` -- rank first by `|Δ exptime|`, then by `|Δ mjd-obs|`
    - ``postFilterRecipes``, ``postFilterTagLike`` -- after ranking, downstream SOFs of these recipes keep only
      the matches whose tag is SQL `LIKE` this pattern
    """

    calType: str
    upstreamRecipes: tuple[str, ...]
    nearestCount: int = 1
    matchRospeed: bool = False
    matchBinningGainUnlessNir: bool = False
    matchSlit: bool = False
    waivedForDownstreamRecipes: tuple[str, ...] = ()
    nirSlitFallback: tuple[str, str] | None = None
    upstreamArms: tuple[str, ...] | None = None
    upstreamSlit: str | None = None
    upstreamTagLike: str | None = None
    downstreamRecipes: tuple[str, ...] | None = None
    excludedDownstreamRecipes: tuple[str, ...] = ()
    rankByExptime: bool = False
    postFilterRecipes: tuple[str, ...] = ()
    postFilterTagLike: str | None = None


STD_RECIPES = ("stare_std", "nod_std", "offset_std")

CALIBRATION_RULES = {
    rule.calType: rule
    for rule in (
        CalibrationRule(
            calType="mbias",
            upstreamRecipes=("mbias",),
            matchRospeed=True,
            matchBinningGainUnlessNir=True,
        ),
        CalibrationRule(
            calType="mdark",
            upstreamRecipes=("mdark",),
            matchRospeed=True,
            upstreamArms=("NIR",),
            rankByExptime=True,
        ),
        CalibrationRule(
            calType="disp_solution",
            upstreamRecipes=("disp_solution",),
            downstreamRecipes=("order_centres", "spat_solution"),
        ),
        CalibrationRule(
            calType="order_centres",
            upstreamRecipes=("order_centres",),
            downstreamRecipes=("mflat",),
        ),
        CalibrationRule(
            calType="mflat",
            upstreamRecipes=("mflat",),
            nearestCount=2,
            matchRospeed=True,
            matchBinningGainUnlessNir=True,
            matchSlit=True,
            waivedForDownstreamRecipes=("spat_solution",),
            nirSlitFallback=("SLIT1.5", "SLIT5.0"),
            postFilterRecipes=("spat_solution",),
            postFilterTagLike="ORDER%",
        ),
        CalibrationRule(
            calType="spat_solution",
            upstreamRecipes=("spat_solution",),
            nearestCount=2,
        ),
        CalibrationRule(
            calType="std_flux",
            upstreamRecipes=STD_RECIPES,
            upstreamSlit="SLIT5.0",
            upstreamTagLike="RESP_TAB%",
            excludedDownstreamRecipes=STD_RECIPES,
        ),
    )
}

MATCH_COLUMNS = [
    "sof",
    "file",
    "tag",
    "filepath",
    "complete",
    "upstream_status",
    "upstream_recipe",
    "upstream_recipe_order",
    "downstream_recipe",
    "downstream_recipe_order",
    "mjd_diff",
    "rank",
]


def match_calibrations(upstream, downstream, rule):
    """*pair each incomplete downstream SOF with its nearest upstream products under one rule*

    **Key Arguments:**

    - ``upstream`` -- `product_frames` rows that may serve as the calibration (not modified)
    - ``downstream`` -- one row per SOF, as `z_unique_product_sof_files` gives it (not modified)
    - ``rule`` -- the `CalibrationRule` to apply

    **Return:**

    - ``matches`` -- dataframe with `MATCH_COLUMNS`, one row per matched upstream product. Ties at equal
      distance are broken by upstream `mjd-obs`, then `filepath`
    """
    upstream = _eligible_upstream(upstream, rule)
    downstream = _eligible_downstream(downstream, rule)
    if upstream.empty or downstream.empty:
        return pd.DataFrame(columns=MATCH_COLUMNS)

    # A NULL `filepath` SORTS FIRST, AS IN SQLITE, AND CANNOT BREAK THE STRING COMPARISON
    upstream = upstream.assign(_pathKey=upstream["filepath"].fillna(""))
    upstream = upstream.sort_values(["mjd-obs", "_pathKey"], kind="mergesort").reset_index(drop=True)

    rows = []
    for _, sofs in downstream.groupby(_class_columns(rule), dropna=False, sort=False):
        mask = _hard_match_mask(upstream, sofs.iloc[0], rule)
        candidates = upstream[mask]
        if candidates.empty:
            continue
        for _, sof in sofs.iterrows():
            rows.extend(_nearest_rows(candidates, sof, rule))

    matches = pd.DataFrame(rows, columns=MATCH_COLUMNS)
    return _post_filter(matches, rule)


def _eligible_upstream(upstream, rule):
    """*keep the upstream products the rule can ever use*"""
    keep = (
        upstream["recipe"].isin(rule.upstreamRecipes)
        & (upstream["complete"] == 1)
        & (upstream["status"].isna() | (upstream["status"] == PASSING_STATUS))
        & upstream["mjd-obs"].notna()
    )
    if rule.upstreamArms is not None:
        keep &= upstream["eso seq arm"].isin(rule.upstreamArms)
    if rule.upstreamSlit is not None:
        keep &= upstream["slit"] == rule.upstreamSlit
    if rule.upstreamTagLike is not None:
        keep &= sql_like(upstream["eso pro catg"], rule.upstreamTagLike)
    if rule.rankByExptime:
        keep &= upstream["exptime"].notna()
    return upstream[keep]


def _eligible_downstream(downstream, rule):
    """*keep the incomplete downstream SOFs the rule can match*"""
    keep = (downstream["complete"] < 1) & downstream["mjd-obs"].notna()
    if rule.downstreamRecipes is not None:
        keep &= downstream["recipe"].isin(rule.downstreamRecipes)
    if rule.excludedDownstreamRecipes:
        keep &= downstream["recipe"].notna() & ~downstream["recipe"].isin(rule.excludedDownstreamRecipes)
    return downstream[keep]


def _class_columns(rule):
    """*the downstream columns that decide which upstream products pass the rule's hard matches*"""
    columns = ["eso seq arm", "recipe", "recipe_order"]
    if rule.matchRospeed:
        columns.append("rospeed")
    if rule.matchBinningGainUnlessNir:
        columns += ["binning", "gain"]
    if rule.matchSlit:
        columns.append("slit")
    return columns


def _equals(series, value):
    """*SQL equality against one value: a NULL on either side never matches*"""
    if pd.isna(value):
        return pd.Series(False, index=series.index)
    return series.notna() & (series == value)


def _hard_match_mask(upstream, sof, rule):
    """*the upstream products that pass every hard match for one downstream match class*"""
    arm = sof["eso seq arm"]
    recipe = sof["recipe"]
    mask = _equals(upstream["eso seq arm"], arm)
    if pd.isna(sof["recipe_order"]):
        return mask & False
    mask &= upstream["recipe_order"] < sof["recipe_order"]

    isWaived = recipe in rule.waivedForDownstreamRecipes
    if rule.matchRospeed and not isWaived:
        mask &= _equals(upstream["rospeed"], sof["rospeed"])
    if rule.matchBinningGainUnlessNir and arm != NIR_ARM:
        mask &= _equals(upstream["binning"], sof["binning"]) & _equals(upstream["gain"], sof["gain"])
    if rule.matchSlit and not isWaived:
        slitMask = _equals(upstream["slit"], sof["slit"])
        if rule.nirSlitFallback is not None and arm == NIR_ARM and sof["slit"] == rule.nirSlitFallback[1]:
            slitMask |= _equals(upstream["slit"], rule.nirSlitFallback[0])
        mask &= slitMask
    return mask


def _nearest_rows(candidates, sof, rule):
    """*the `nearestCount` nearest candidate rows for one downstream SOF, as dicts keyed by `MATCH_COLUMNS`*"""
    if rule.rankByExptime and pd.notna(sof["exptime"]):
        expDiff = (candidates["exptime"] - sof["exptime"]).abs()
        candidates = candidates[expDiff == expDiff.min()]

    mjds = candidates["mjd-obs"].to_numpy()
    target = sof["mjd-obs"]
    nearest = _nearest_indices(mjds, target, rule.nearestCount)
    chosen = candidates.iloc[nearest]
    distances = np.abs(chosen["mjd-obs"].to_numpy() - target)
    order = np.lexsort((chosen["_pathKey"].to_numpy(), chosen["mjd-obs"].to_numpy(), distances))
    order = order[: rule.nearestCount]

    rows = []
    for rank, position in enumerate(order, start=1):
        match = chosen.iloc[position]
        rows.append(
            {
                "sof": sof["sof"],
                "file": match["file"],
                "tag": match["eso pro catg"],
                "filepath": match["filepath"],
                "complete": sof["complete"],
                "upstream_status": match["status"],
                "upstream_recipe": match["recipe"],
                "upstream_recipe_order": match["recipe_order"],
                "downstream_recipe": sof["recipe"],
                "downstream_recipe_order": sof["recipe_order"],
                "mjd_diff": distances[position],
                "rank": rank,
            }
        )
    return rows


def _nearest_indices(mjds, target, count):
    """*indices of every value in sorted `mjds` at most as far from `target` as its `count`-th nearest*

    The window holds the `count` nearest values plus any value tied with the last of them, so the caller's
    tie-break sees every tied candidate.
    """
    size = len(mjds)
    insertAt = int(np.searchsorted(mjds, target))
    low = max(insertAt - count, 0)
    high = min(insertAt + count, size)
    distances = np.abs(mjds[low:high] - target)
    kth = np.sort(distances)[min(count, len(distances)) - 1]
    while low > 0 and abs(mjds[low - 1] - target) <= kth:
        low -= 1
    while high < size and abs(mjds[high] - target) <= kth:
        high += 1
    return np.arange(low, high)


def _post_filter(matches, rule):
    """*drop ranked matches the rule's post-filter rejects*"""
    if not rule.postFilterRecipes or matches.empty:
        return matches
    isFiltered = matches["downstream_recipe"].isin(rule.postFilterRecipes)
    keep = ~isFiltered | sql_like(matches["tag"], rule.postFilterTagLike)
    return matches[keep].reset_index(drop=True)


def sql_like(series, pattern):
    """*SQLite `LIKE` on a string series: `%` and `_` wildcards, ASCII case-insensitive, NULL never matches*"""
    regex = "".join(".*" if char == "%" else "." if char == "_" else re.escape(char) for char in pattern)
    compiled = re.compile(regex, re.ASCII | re.IGNORECASE | re.DOTALL)
    return series.map(lambda value: isinstance(value, str) and compiled.fullmatch(value) is not None).astype(bool)


def refresh_calibration_matches(conn, calTypes=None, *, recipe=None, arm=None):
    """*recompute the `calibration_matches` rows that the `cal_<type>` views read*

    The rows for each type are replaced in one savepoint, so a reader never sees a half-written type. On a
    connection already inside a transaction, committing stays with the caller.

    **Key Arguments:**

    - ``conn`` -- the open SQLite connection to a workspace `soxspipe.db`
    - ``calTypes`` -- the calibration types to refresh, or one type name; all of `CALIBRATION_RULES` when `None`
    - ``recipe``, ``arm`` -- give both to refresh only the downstream SOFs of this recipe on this arm; rows for
      every other SOF are left as they are

    **Raises:**

    - ``ValueError`` -- for a calibration type with no rule, or when only one of ``recipe`` and ``arm`` is given
    """
    if calTypes is None:
        calTypes = list(CALIBRATION_RULES)
    elif isinstance(calTypes, str):
        calTypes = [calTypes]
    else:
        calTypes = list(calTypes)
    unknown = [calType for calType in calTypes if calType not in CALIBRATION_RULES]
    if unknown:
        raise ValueError(f"no calibration matching rule for: {', '.join(map(str, unknown))}")
    if (recipe is None) != (arm is None):
        raise ValueError("give both `recipe` and `arm` to scope a calibration-match refresh, or neither")

    downstream = _read_downstream(conn, recipe, arm)
    matchesByType = {
        calType: match_calibrations(_read_upstream(conn, CALIBRATION_RULES[calType]), downstream, rule)
        for calType, rule in ((calType, CALIBRATION_RULES[calType]) for calType in calTypes)
    }
    _replace_matches(conn, matchesByType, recipe, arm)


def _read_downstream(conn, recipe, arm):
    """*one row per incomplete SOF, optionally only those of one recipe on one arm*"""
    sqlQuery = f"SELECT {_SELECT_COLUMNS} FROM z_unique_product_sof_files WHERE complete < 1"  # noqa: S608
    params = ()
    if recipe is not None:
        sqlQuery += ' AND recipe = ? AND "eso seq arm" = ?'
        params = (recipe, arm)
    return pd.read_sql(sqlQuery, con=conn, params=params)


def _read_upstream(conn, rule):
    """*the complete products of the rule's upstream recipes*"""
    placeholders = ", ".join("?" for _ in rule.upstreamRecipes)
    sqlQuery = (
        f"SELECT {_SELECT_COLUMNS} FROM product_frames "  # noqa: S608
        f"WHERE complete = 1 AND recipe IN ({placeholders})"
    )
    return pd.read_sql(sqlQuery, con=conn, params=rule.upstreamRecipes)


def _replace_matches(conn, matchesByType, recipe, arm):
    """*swap in the new rows for each type inside one savepoint*"""
    insertQuery = (
        f"INSERT INTO calibration_matches (cal_type, {', '.join(MATCH_COLUMNS)}) "  # noqa: S608
        f"VALUES (?, {', '.join('?' for _ in MATCH_COLUMNS)})"
    )
    conn.execute("SAVEPOINT refresh_calibration_matches")
    try:
        for calType, matches in matchesByType.items():
            if recipe is None:
                conn.execute("DELETE FROM calibration_matches WHERE cal_type = ?", (calType,))
            else:
                conn.execute(
                    "DELETE FROM calibration_matches WHERE cal_type = ? AND sof IN "
                    '(SELECT sof FROM product_frames WHERE recipe = ? AND "eso seq arm" = ?)',
                    (calType, recipe, arm),
                )
            values = matches.astype(object).where(matches.notna(), None)
            conn.executemany(insertQuery, [(calType, *row) for row in values.itertuples(index=False)])
    except BaseException:
        # SQLITE MAY ALREADY HAVE ROLLED THE SAVEPOINT BACK (E.G. DISK FULL); KEEP THE ORIGINAL ERROR EITHER WAY
        with contextlib.suppress(sqlite3.Error):
            conn.execute("ROLLBACK TO refresh_calibration_matches")
            conn.execute("RELEASE refresh_calibration_matches")
        raise
    conn.execute("RELEASE refresh_calibration_matches")


def has_calibration_matches_table(conn):
    """*whether the database has the `calibration_matches` table (workspaces made before DY-1513 do not)*"""
    sqlQuery = "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'calibration_matches'"
    return conn.execute(sqlQuery).fetchone() is not None
