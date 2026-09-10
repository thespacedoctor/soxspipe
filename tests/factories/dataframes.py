"""Factories for deterministic pipeline tables."""

from __future__ import annotations

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
