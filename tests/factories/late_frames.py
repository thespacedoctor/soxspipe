"""Raw-frame factories and a loader for workspaces that receive frames after their first `prep`."""

from __future__ import annotations

import sqlite3

import pandas as pd

from .dataframes import raw_frame_table


def insert_frames(connection: sqlite3.Connection, frames: pd.DataFrame) -> None:
    """Insert dataframe rows into `raw_frames`, keeping only the columns the table has."""
    columns = {row[1] for row in connection.execute("PRAGMA table_info(raw_frames);")}
    known = frames.loc[:, frames.columns.isin(columns)]
    known.to_sql("raw_frames", con=connection, index=False, if_exists="append")


def pinhole_frame(biases: pd.DataFrame, name: str, dprType: str, tplStart: str, obsTime: str) -> pd.DataFrame:
    """Return one VIS pinhole exposure that is its own template run."""
    frame = biases.iloc[[0]].copy()
    frame["file"] = f"{name}.fits"
    frame["filepath"] = f"./raw/2024-01-01/{name}.fits"
    frame["eso dpr type"] = dprType
    frame["eso dpr tech"] = "ECHELLE,PINHOLE"
    frame["eso tpl nexp"] = 1
    frame["eso tpl name"] = "SOXS_orderdef" if dprType == "LAMP,FLAT" else "SOXS_arc"
    frame["exptime"] = 10.0
    frame["date-obs"] = f"{obsTime}.000"
    frame["mjd-obs"] = frame["mjd-obs"] + 0.04
    frame["eso tpl start"] = tplStart
    frame["set_first_file"] = f"{name}.fits"
    return frame


def calibration_chain_frames() -> pd.DataFrame:
    """Return two biases, a pinhole arc and an order-centre flat: everything a slit-flat SOF needs upstream."""
    biases = raw_frame_table()
    arc = pinhole_frame(biases, "arc-1", "LAMP,WAVE", "2024-01-02T04:04:00", "2024-01-02T04:04:05")
    orderCentreFlat = pinhole_frame(biases, "ocflat-1", "LAMP,FLAT", "2024-01-02T05:04:00", "2024-01-02T05:04:05")
    return pd.concat([biases, arc, orderCentreFlat])


def spatial_solution_arc() -> pd.DataFrame:
    """Return a multi-pinhole arc: the input of the spatial-solution recipe."""
    frame = pinhole_frame(raw_frame_table(), "spat-1", "LAMP,WAVE", "2024-01-02T07:04:00", "2024-01-02T07:04:05")
    frame["eso dpr tech"] = "ECHELLE,MULTI-PINHOLE"
    return frame


def flux_standard_stare(count: int = 2) -> pd.DataFrame:
    """Return one OB of flux-standard stare exposures: the input of the `stare_std` recipe."""
    template = raw_frame_table().iloc[0].to_dict()
    rows = []
    for index in range(count):
        name = f"std-{index + 1}.fits"
        rows.append(
            {
                **template,
                "eso dpr catg": "SCIENCE",
                "eso dpr tech": "ECHELLE,SLIT,STARE",
                "eso dpr type": "STD,FLUX",
                "eso obs id": "201",
                # A STANDARD STAR THAT IS NOT IN THE STATIC LIBRARY IS FLAGGED TO IGNORE
                "eso obs name": "GD71",
                "eso obs targ name": "GD71",
                "exptime": 60.0,
                "slit": "1.0",
                "object": "STD",
                "template": "SOXS_slt_cal_Stare",
                "eso tpl name": "SOXS_slt_cal_Stare",
                "eso tpl nexp": count,
                "eso tpl expno": index + 1,
                "eso tpl start": "2024-01-02T08:00:00",
                "file": name,
                "filepath": f"./raw/2024-01-01/{name}",
                "date-obs": f"2024-01-02T08:01:{index:02d}.000",
                "mjd-obs": 60311.34 + index * 0.0001,
                "set_first_file": "std-1.fits",
            }
        )
    return pd.DataFrame(rows)
