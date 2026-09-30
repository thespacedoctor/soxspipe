#!/usr/bin/env python
"""
*Find the sets a late-arriving raw frame belongs to, reset them so they are reduced again, and keep their SOF names*

Author
: David Young

Date Created
: September 30, 2026
"""

from __future__ import annotations

import logging
import os
import sqlite3
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from pathlib import Path

from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

# NO IMPORTS FROM OTHER SOXSPIPE MODULES EXCEPT sql_identifiers: THE PACKAGE `__init__` IMPORT ORDER IS
# LOAD-BEARING, AND THIS MODULE MUST STAY SAFE TO IMPORT FROM ANYWHERE.

# THESE TECHNIQUES ARE GROUPED ONE EXPOSURE AT A TIME, SO A LATE FRAME NEVER JOINS AN EXISTING SET
PER_EXPOSURE_TECHS = ("ECHELLE,SLIT,STARE", "ECHELLE,PINHOLE", "ECHELLE,MULTI-PINHOLE")

# PRODUCTS WITH NO FILE OF THEIR OWN CARRY THIS PLACEHOLDER NAME
PLACEHOLDER_PRODUCT_FILE = "XXXX"

# THE MOST VALUES BOUND TO ONE STATEMENT (SQLITE BUILDS OLDER THAN 3.32 ALLOW ONLY 999)
MAX_BOUND_VALUES = 500

# NIR IMAGE FRAMES THAT ARE NOT DARKS ARE "OFF" FRAMES, WHICH THE GROUPING LEAVES OUT (SEE get_raw_frames_and_groups)
# A NULL ARM OR TECHNIQUE IS "--" IN THE GROUPING, SO IT IS "--" HERE TOO AND THE FRAME STAYS IN
_NOT_AN_OFF_FRAME = (
    "NOT (COALESCE(r.`eso dpr tech`, '--') = 'IMAGE' AND COALESCE(r.`eso seq arm`, '--') = 'NIR' "
    "AND COALESCE(r.`eso dpr type`, '--') != 'DARK')"
)


def _placeholders(count: int) -> str:
    return ", ".join("?" for _ in range(count))


def _chunks(values: Sequence[str]) -> Iterator[Sequence[str]]:
    for start in range(0, len(values), MAX_BOUND_VALUES):
        yield values[start : start + MAX_BOUND_VALUES]


def _late_frame_exists(connection: sqlite3.Connection) -> bool:
    """*check whether any unprocessed frame could be a late member of an existing set*"""
    sqlQuery = (
        "SELECT 1 FROM raw_frames_valid WHERE processed = 0 AND `eso tpl start` IS NOT NULL "  # noqa: S608
        f"AND COALESCE(`eso dpr tech`, '--') NOT IN ({_placeholders(len(PER_EXPOSURE_TECHS))}) LIMIT 1;"
    )
    return connection.execute(sqlQuery, PER_EXPOSURE_TECHS).fetchone() is not None


def find_rejoined_sofs(
    connection: sqlite3.Connection, sofMapTable: str, groupingKeys: Sequence[str]
) -> dict[str, tuple[str, frozenset[str]]]:
    """*find the SOFs that a frame added after they were grouped now belongs to*

    A late frame is unprocessed and in no SOF. It joins a SOF when it matches that SOF's frames on every grouping
    key and on `eso tpl start`. A frame with no `eso tpl start` never joins, and neither do frames of a technique
    that is grouped one exposure at a time.

    **Key Arguments:**

    - ``connection`` -- open connection to the workspace database
    - ``sofMapTable`` -- the current session's `sof_map_<id>` table (checked again by `validate_sql_identifier`)
    - ``groupingKeys`` -- the `raw_frames` columns the frames of one SOF share; they never enter the SQL text

    **Return:**

    - ``rejoined`` -- SOF name to (recipe, frozenset of the raw frame filepaths the SOF held), for every SOF that
      gained a late frame
    """
    import pandas as pd

    sofMapTable = validate_sql_identifier(sofMapTable, "sof map table name")
    if not _late_frame_exists(connection):
        return {}

    techPlaceholders = _placeholders(len(PER_EXPOSURE_TECHS))
    # TABLE NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier ABOVE
    sqlQuery = f"""SELECT r.*, s.sof AS member_sof
        FROM raw_frames_valid r
        LEFT JOIN (SELECT DISTINCT sof, filepath FROM {sofMapTable}) s ON s.filepath = r.filepath
        WHERE COALESCE(r.`eso dpr tech`, '--') NOT IN ({techPlaceholders})
        AND {_NOT_AN_OFF_FRAME}
        AND r.`eso tpl start` IN (
            SELECT `eso tpl start` FROM raw_frames_valid
            WHERE processed = 0 AND `eso tpl start` IS NOT NULL
        );"""  # noqa: S608
    frames = pd.read_sql(sqlQuery, con=connection, params=PER_EXPOSURE_TECHS)

    membersBySof: dict[str, set[str]] = {}
    for _, group in frames.groupby(list(groupingKeys) + ["eso tpl start"], dropna=False):
        isLate = (group["processed"] == 0) & group["member_sof"].isna()
        if not isLate.any():
            continue
        for sof, filepath in group.loc[group["member_sof"].notna(), ["member_sof", "filepath"]].itertuples(index=False):
            membersBySof.setdefault(sof, set()).add(filepath)

    recipes = _recipes_of(connection, list(membersBySof))
    return {sof: (recipes[sof], frozenset(members)) for sof, members in sorted(membersBySof.items()) if sof in recipes}


def _recipes_of(connection: sqlite3.Connection, sofs: Sequence[str]) -> dict[str, str]:
    """*read the recipe that produces each SOF from `product_frames`*"""
    recipes: dict[str, str] = {}
    for chunk in _chunks(sofs):
        # `_placeholders` IS ALWAYS LITERAL `?` MARKS; THE VALUES ARE BOUND
        sqlQuery = f"SELECT DISTINCT sof, recipe FROM product_frames WHERE sof IN ({_placeholders(len(chunk))});"  # noqa: S608
        for sof, recipe in connection.execute(sqlQuery, tuple(chunk)):
            recipes.setdefault(sof, recipe)
    return recipes


def status_columns(connection: sqlite3.Connection) -> list[str]:
    """*list `status` and every per-session `status_<id>` column of `product_frames`*

    **Key Arguments:**

    - ``connection`` -- open connection to the workspace database

    **Return:**

    - ``columns`` -- the validated column names
    """
    names = [row[1] for row in connection.execute("PRAGMA table_info(product_frames);")]
    return [
        validate_sql_identifier(name, "status column")
        for name in names
        if name == "status" or name.startswith("status_")
    ]


def reset_rejoined_sofs(
    connection: sqlite3.Connection,
    sofMapTable: str,
    rejoined: Mapping[str, tuple[str, frozenset[str]]],
    workspaceRoot: str | os.PathLike,
    validateOwnedPath: Callable[[Path, Path, str], Path],
    log: logging.Logger,
) -> list[Path]:
    """*queue rejoined SOFs to be grouped and reduced again*

    Marks each SOF's products incomplete with every status and `set_first_file` cleared, frees its raw frames to
    be grouped again, and drops its rows from the SOF map and the raw frame sets. The `product_frames` rows
    themselves are kept. Nothing is deleted from disk: the caller commits, then passes the returned paths to
    `delete_stale_files`, so a failed commit never leaves a product deleted that the database still expects.

    **Key Arguments:**

    - ``connection`` -- open connection to the workspace database
    - ``sofMapTable`` -- the current session's `sof_map_<id>` table
    - ``rejoined`` -- the result of `find_rejoined_sofs`
    - ``workspaceRoot`` -- the workspace root; no file outside it is deleted
    - ``validateOwnedPath`` -- `data_organiser._validate_owned_path`, passed in to avoid a circular import
    - ``log`` -- logger

    **Return:**

    - ``stalePaths`` -- the validated product and `_ERROR.log` paths to delete once the changes are committed
    """
    sofMapTable = validate_sql_identifier(sofMapTable, "sof map table name")
    sofs = sorted(rejoined)
    if not sofs:
        return []

    _clear_product_state(connection, sofs)
    filepaths = sorted({filepath for _, members in rejoined.values() for filepath in members})
    for chunk in _chunks(filepaths):
        # `_placeholders` IS ALWAYS LITERAL `?` MARKS; THE VALUES ARE BOUND
        connection.execute(
            f"UPDATE raw_frames SET processed = 0 WHERE filepath IN ({_placeholders(len(chunk))});",  # noqa: S608
            tuple(chunk),
        )
    for chunk in _chunks(sofs):
        # TABLE NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier ABOVE
        for table in (sofMapTable, "raw_frame_sets"):
            connection.execute(
                f"DELETE FROM {table} WHERE sof IN ({_placeholders(len(chunk))});",  # noqa: S608
                tuple(chunk),
            )
    return stale_product_paths(connection, sofs, workspaceRoot, validateOwnedPath, log)


def _clear_product_state(connection: sqlite3.Connection, sofs: Sequence[str]) -> None:
    """*mark the products of ``sofs`` incomplete and clear `set_first_file`, `status` and every `status_<id>` column*"""
    assignments = ", ".join(
        ["complete = 0", "set_first_file = NULL"] + [f"{column} = NULL" for column in status_columns(connection)]
    )
    for chunk in _chunks(sofs):
        # COLUMN NAMES CANNOT BE BOUND; EACH IS CHECKED BY validate_sql_identifier IN status_columns
        connection.execute(
            f"UPDATE product_frames SET {assignments} WHERE sof IN ({_placeholders(len(chunk))});",  # noqa: S608
            tuple(chunk),
        )


def stale_product_paths(
    connection: sqlite3.Connection,
    sofs: Sequence[str],
    workspaceRoot: str | os.PathLike,
    validateOwnedPath: Callable[[Path, Path, str], Path],
    log: logging.Logger,
) -> list[Path]:
    """*list the product file and `_ERROR.log` of each SOF, so a new reduction does not skip over them*

    Never raises: a path outside the workspace `reduced` directory, or a symlink, is left out and logged as a
    warning.

    **Key Arguments:**

    - ``connection`` -- open connection to the workspace database
    - ``sofs`` -- the SOF names whose products are stale
    - ``workspaceRoot`` -- the workspace root; only paths inside its `reduced` directory are returned
    - ``validateOwnedPath`` -- `data_organiser._validate_owned_path`
    - ``log`` -- logger

    **Return:**

    - ``stalePaths`` -- the validated paths; some may not exist
    """
    rootPath = Path(workspaceRoot)
    # `reduced` IS A SYMLINK TO THE CURRENT SESSION'S DIRECTORY; `validateOwnedPath` COMPARES RESOLVED PATHS
    reducedPath = rootPath / "reduced"
    stalePaths: list[Path] = []
    for chunk in _chunks(sofs):
        # `_placeholders` IS ALWAYS LITERAL `?` MARKS; THE VALUES ARE BOUND
        sqlQuery = (
            f"SELECT filepath FROM product_frames WHERE sof IN ({_placeholders(len(chunk))}) "  # noqa: S608
            "AND filepath IS NOT NULL AND file != ?;"
        )
        for (filepath,) in connection.execute(sqlQuery, (*chunk, PLACEHOLDER_PRODUCT_FILE)).fetchall():
            productPath = rootPath / filepath
            errorLogPath = productPath.with_name(productPath.stem + "_ERROR.log")
            for candidate in (productPath, errorLogPath):
                if candidate.is_symlink():
                    log.warning(f"late_frame_rejoin: `{candidate}` is a symlink and was not removed")
                    continue
                try:
                    stalePaths.append(validateOwnedPath(candidate, reducedPath, "stale product path"))
                except (OSError, ValueError) as error:
                    log.warning(f"late_frame_rejoin: `{candidate}` is not a safe path and was not removed: {error}")
    return stalePaths


def delete_stale_files(paths: Iterable[Path], log: logging.Logger) -> None:
    """*delete the stale product files found by `reset_rejoined_sofs`*

    Never raises: a file that is already gone is fine, and any other failure is logged as a warning.

    **Key Arguments:**

    - ``paths`` -- the validated paths to delete
    - ``log`` -- logger
    """
    for path in paths:
        try:
            path.unlink(missing_ok=True)
        except OSError as error:
            log.warning(f"late_frame_rejoin: stale file `{path}` could not be removed: {error}")


def pin_rejoined_sof_names(rawGroups, rejoined: Mapping[str, tuple[str, frozenset[str]]]):
    """*give each regrouped set the SOF name it had before its late frame arrived*

    A group takes an old name only when its recipe matches and it holds every frame the old SOF held.

    **Key Arguments:**

    - ``rawGroups`` -- the grouped raw frames, with `sof`, `recipe` and `filepaths` columns; never changed
    - ``rejoined`` -- the result of `find_rejoined_sofs`

    **Return:**

    - ``pinnedGroups`` -- a new dataframe, with the old SOF names where they apply
    """
    pinnedGroups = rawGroups.copy()
    if not len(pinnedGroups.index) or not rejoined:
        return pinnedGroups

    def pinned_name(group) -> str:
        groupFilepaths = set(group["filepaths"])
        for sof, (recipe, members) in rejoined.items():
            if group["recipe"] == recipe and groupFilepaths >= members:
                return sof
        return group["sof"]

    pinnedGroups["sof"] = pinnedGroups.apply(pinned_name, axis=1)
    return pinnedGroups


def format_rejoin_summary(sofs: Iterable[str]) -> str:
    """*return the one-line message naming the SOFs that were regrouped*

    **Key Arguments:**

    - ``sofs`` -- the SOF names that gained a late frame

    **Return:**

    - ``summary`` -- the message
    """
    names = sorted(sofs)
    return f"{len(names)} SOF(S) GAINED LATE-ARRIVING FRAMES AND WILL BE REDUCED AGAIN: {', '.join(names)}"
