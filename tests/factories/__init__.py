"""Fresh deterministic test-data factories."""

from .dataframes import (
    dispersion_table,
    order_table,
    product_table,
    qc_row,
    qc_table,
    raw_flat_run_table,
    raw_frame_table,
    raw_group_table,
)
from .files import (
    dispersion_map_fits,
    harvestable_raw_fits,
    lzw_compress,
    lzw_compressed_fits,
    order_table_fits,
    prepared_fits,
    raw_fits,
    sof_file,
)
from .frames import instrument_header, synthetic_ccd
from .settings import pipeline_settings
from .signals import synthetic_signal
from .workspaces import WorkspaceLayout, workspace_layout, workspace_organiser

__all__ = [
    "WorkspaceLayout",
    "dispersion_map_fits",
    "dispersion_table",
    "harvestable_raw_fits",
    "instrument_header",
    "lzw_compress",
    "lzw_compressed_fits",
    "order_table",
    "order_table_fits",
    "pipeline_settings",
    "prepared_fits",
    "product_table",
    "qc_row",
    "qc_table",
    "raw_fits",
    "raw_flat_run_table",
    "raw_frame_table",
    "raw_group_table",
    "sof_file",
    "synthetic_ccd",
    "synthetic_signal",
    "workspace_layout",
    "workspace_organiser",
]
