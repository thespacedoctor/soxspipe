"""Fresh deterministic test-data factories."""

from .dataframes import (
    dispersion_table,
    merged_spectrum,
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
    response_product_fits,
    sof_file,
)
from .frames import StubInputFrames, instrument_header, noise_frame_pair, synthetic_ccd
from .merge_extractors import soxs_merge_extractor
from .late_frames import (
    calibration_chain_frames,
    flux_standard_stare,
    insert_frames,
    pinhole_frame,
    spatial_solution_arc,
)
from .settings import logging_settings, pipeline_settings
from .signals import synthetic_signal
from .slit_illumination import illumination_scene, synthetic_illumination_scene
from .workspaces import WorkspaceLayout, session_workspace, workspace_layout, workspace_organiser

__all__ = [
    "calibration_chain_frames",
    "flux_standard_stare",
    "insert_frames",
    "pinhole_frame",
    "spatial_solution_arc",
    "illumination_scene",
    "StubInputFrames",
    "WorkspaceLayout",
    "dispersion_map_fits",
    "dispersion_table",
    "harvestable_raw_fits",
    "instrument_header",
    "logging_settings",
    "lzw_compress",
    "lzw_compressed_fits",
    "merged_spectrum",
    "noise_frame_pair",
    "order_table",
    "order_table_fits",
    "pipeline_settings",
    "prepared_fits",
    "product_table",
    "soxs_merge_extractor",
    "qc_row",
    "qc_table",
    "raw_fits",
    "raw_flat_run_table",
    "raw_frame_table",
    "raw_group_table",
    "response_product_fits",
    "session_workspace",
    "sof_file",
    "synthetic_ccd",
    "synthetic_illumination_scene",
    "synthetic_signal",
    "workspace_layout",
    "workspace_organiser",
]
