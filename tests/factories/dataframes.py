"""Factories for deterministic pipeline tables."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd


def order_table() -> pd.DataFrame:
    """Return a fresh minimal order-coordinate table."""
    return pd.DataFrame(
        {
            "order": [10.0, 11.0, 12.0],
            "wavelength": [500.0, 510.0, 520.0],
            "slit_position": [-1.0, 0.0, 1.0],
            "x": [4.0, 8.0, 12.0],
            "y": [5.0, 9.0, 13.0],
        }
    )


def dispersion_table() -> pd.DataFrame:
    """Return a fresh minimal dispersion-coefficient table."""
    return pd.DataFrame(
        {
            "axis": ["x", "y"],
            "order_deg": [1, 1],
            "wavelength_deg": [1, 1],
            "slit_deg": [0, 0],
            "c000": [1.0, 2.0],
            "c010": [0.1, 0.2],
            "c100": [0.01, 0.02],
            "c110": [0.001, 0.002],
        }
    )


def qc_table() -> pd.DataFrame:
    """Return a fresh QC table with the production column contract."""
    return pd.DataFrame(
        [
            {
                "soxspipe_recipe": "soxs-mbias",
                "qc_name": "RON",
                "qc_value": 3.2,
                "qc_unit": "electron",
                "qc_order": None,
                "qc_comment": "Synthetic read noise",
                "obs_date_utc": "2024-01-02T03:04:05.678",
                "reduction_date_utc": "2024-01-02T04:05:06.789",
                "to_header": True,
            }
        ]
    )


def qc_row(
    *,
    recipeName: str,
    name: str,
    value: object,
    unit: str,
    comment: str,
    toHeader: bool = True,
) -> pd.DataFrame:
    """Return one deterministic QC row with the complete production contract."""
    return qc_table().assign(
        soxspipe_recipe=recipeName,
        qc_name=name,
        qc_value=value,
        qc_unit=unit,
        qc_comment=comment,
        to_header=toHeader,
    )


def product_table() -> pd.DataFrame:
    """Return a fresh product table with deterministic metadata."""
    return pd.DataFrame(
        [
            {
                "soxspipe_recipe": "soxs-mbias",
                "product_label": "MASTER_BIAS_VIS",
                "file_name": "MASTER_BIAS_VIS.fits",
                "file_type": "FITS",
                "obs_date_utc": "2024-01-02T03:04:05.678",
                "reduction_date_utc": "2024-01-02T04:05:06.789",
                "product_desc": "Synthetic master bias",
                "file_path": "products/MASTER_BIAS_VIS.fits",
                "label": "PROD",
            }
        ]
    )


def raw_group_table() -> pd.DataFrame:
    """Return one complete synthetic raw-frame group for organizer contracts."""
    return pd.DataFrame(
        [
            {
                "eso seq arm": "VIS",
                "eso dpr catg": "CALIB",
                "eso dpr tech": "IMAGE",
                "eso dpr type": "BIAS",
                "night start date": "2024-01-01",
                "binning": "1x1",
                "rospeed": "fast",
                "lamp": "--",
                "slit": "--",
                "exptime": 0.0,
                "instrume": "SOXS",
                "recipe": "mbias",
                "recipe_order": 0,
                "sof": "20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof",
                "counts": 3,
                "complete": 1,
                "date-obs": "20240102T030405",
                "filepaths": [
                    "./raw/2024-01-01/bias-1.fits",
                    "./raw/2024-01-01/bias-2.fits",
                ],
            }
        ]
    )


def raw_frame_table() -> pd.DataFrame:
    """Return two bias frames with the complete organizer grouping schema."""
    sharedValues: dict[str, object] = {
        "eso seq arm": "VIS",
        "eso dpr catg": "CALIB",
        "eso dpr tech": "IMAGE",
        "eso dpr type": "BIAS",
        "eso pro catg": "--",
        "eso pro tech": "--",
        "eso pro type": "--",
        "eso obs id": "42",
        "eso obs name": "Synthetic bias",
        "exptime": 0.0,
        "binning": "1x1",
        "rospeed": "fast",
        "slit": "--",
        "slitmask": "--",
        "lamp": "--",
        "gain": 1.2,
        "night start date": "2024-01-01",
        "night start mjd": 60310,
        "mjd-date": "2024-01-02",
        "object": "BIAS",
        "template": "SOXS_cal_bias",
        "instrume": "SOXS",
        "absrot": 12.5,
        "eso tpl name": "SOXS_cal_bias",
        "eso tpl nexp": 2,
        "eso tpl expno": 1,
        "filter": "--",
        "ra": 0.0,
        "dec": 0.0,
        "simulation": 0,
        "eso tel parang end": 0.0,
        "eso tel parang start": 0.0,
        "eso tel az": 0.0,
        "eso tel alt": 0.0,
        "eso tel ambi fwhm end": 0.0,
        "eso tel ambi fwhm start": 0.0,
        "eso tel airm end": 1.0,
        "eso tel airm start": 1.0,
        "eso obs targ name": "BIAS",
        "eso obs tplno": 1,
        "eso obs ntpl": 1,
        "eso tpl start": "2024-01-02T03:04:00",
        "eso obs start": "2024-01-02T03:04:00",
        "nir temp k": 150.0,
        "vis temp c": 10.0,
        "cp temp c": 11.0,
        "afc1 pos1": 0.0,
        "afc1 pos2": 0.0,
        "afc2 pos1": 0.0,
        "afc2 pos2": 0.0,
        "set_first_file": "bias-1.fits",
        "processed": 0,
    }
    rows = []
    for index, second in enumerate((5, 6), start=1):
        rows.append(
            {
                **sharedValues,
                "file": f"bias-{index}.fits",
                "filepath": f"./raw/2024-01-01/bias-{index}.fits",
                "date-obs": f"2024-01-02T03:04:0{second}.000",
                "mjd-obs": 60311.1278 + (index - 1) * 0.0001,
                "eso tpl expno": index,
            }
        )
    return pd.DataFrame(rows)


FLAT_RUN_DEFAULT_TPL_START = "2024-01-02T06:04:00"
SECONDS_PER_DAY = 86400.0


def raw_flat_run_table(
    expnos: Sequence[int],
    tplStart: str = FLAT_RUN_DEFAULT_TPL_START,
    *,
    namePrefix: str = "flat",
    slit: str = "1.0",
    nightOffset: int = 0,
    nexp: int = 10,
) -> pd.DataFrame:
    """Return the slit-flat exposures ``expnos`` of one template run, with the complete organizer schema.

    Exposure ``n`` is taken ``10 + n`` seconds after ``tplStart``, so a run started at 06:04:00 has exposure 1 at
    06:04:11. Passing only some of the exposures models a run that has not fully arrived yet.
    """
    template = raw_frame_table()
    startTime = pd.Timestamp(tplStart)
    rows = []
    for expno in expnos:
        obsTime = startTime + pd.Timedelta(seconds=10 + expno)
        secondsIntoDay = obsTime.hour * 3600 + obsTime.minute * 60 + obsTime.second
        row = template.iloc[0].to_dict()
        row.update(
            {
                "file": f"{namePrefix}-{expno}.fits",
                "filepath": f"./raw/2024-01-01/{namePrefix}-{expno}.fits",
                "eso dpr type": "LAMP,FLAT",
                "eso dpr tech": "ECHELLE,SLIT",
                "eso tpl name": "SOXS_flat",
                "template": "SOXS_flat",
                "object": "FLAT",
                "eso tpl nexp": nexp,
                "eso tpl expno": expno,
                "exptime": 10.0,
                "slit": slit,
                "date-obs": obsTime.strftime("%Y-%m-%dT%H:%M:%S.000"),
                "mjd-obs": 60311 + nightOffset + secondsIntoDay / SECONDS_PER_DAY,
                "eso tpl start": tplStart,
                "night start date": (pd.Timestamp("2024-01-01") + pd.Timedelta(days=nightOffset)).strftime("%Y-%m-%d"),
                "night start mjd": 60310 + nightOffset,
                "set_first_file": f"{namePrefix}-{expnos[0]}.fits",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)
