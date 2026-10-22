"""The order location tables carry an `ORDER_TRACES` extension of sampled trace coordinates (DY-1726)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.toolkit import unpack_order_table
from tests.factories import instrument_header, pipeline_settings

pytestmark = pytest.mark.unit

TRACES_EXTNAME = "ORDER_TRACES"
CENTRE_POLYNOMIALS = {
    "degorder_cent": 1,
    "degy_cent": 1,
    "cent_00": 1.0,
    "cent_01": 2.0,
    "cent_10": 3.0,
    "cent_11": 0.5,
    "std_00": 2.0,
    "std_01": 0.0,
    "std_10": 0.0,
    "std_11": 0.0,
}
EDGE_POLYNOMIALS = {
    "degorder_edgeup": 1,
    "degy_edgeup": 1,
    "edgeup_00": 11.0,
    "edgeup_01": 2.0,
    "edgeup_10": 3.0,
    "edgeup_11": 0.5,
    "degorder_edgelow": 1,
    "degy_edgelow": 1,
    "edgelow_00": -9.0,
    "edgelow_01": 2.0,
    "edgelow_10": 3.0,
    "edgelow_11": 0.5,
}


def _writer(tmp_path: Path, log: object, *, recipeName: str, binx: int = 1, biny: int = 1) -> detect_continuum:
    """A detector holding only the state `write_order_table_to_file` reads."""
    settings = pipeline_settings(tmp_path)
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.settings = settings
    detector.kw = keyword_lookup(log=log, settings=settings).get
    detector.arm = "VIS"
    detector.inst = "SOXS"
    detector.recipeName = recipeName
    detector.sofName = "synthetic_sof"
    detector.lampTag = False
    detector.notFlattened = ""
    detector.noddingSequence = ""
    detector.productDir = str(tmp_path)
    detector.qc = None
    detector.binx = binx
    detector.biny = biny
    return detector


def _write(detector: detect_continuum, polynomials: dict, metadata: pd.DataFrame) -> Path:
    frame = SimpleNamespace(header=instrument_header(arm="VIS"))
    return Path(
        detector.write_order_table_to_file(
            frame=frame, orderPolyTable=pd.DataFrame([polynomials]), orderMetaTable=metadata
        )
    )


def _metadata(ymax: float = 6.0) -> pd.DataFrame:
    return pd.DataFrame({"order": [10, 11], "ymin": [0.0, 2.0], "ymax": [ymax, ymax + 2.0]})


def _traces(path: Path) -> pd.DataFrame:
    return Table.read(path, format="fits", hdu=3).to_pandas()


def test_order_centre_table_gains_a_third_extension_named_order_traces(tmp_path: Path, log: object) -> None:
    # ARRANGE
    detector = _writer(tmp_path, log, recipeName="soxs-order-centre")

    # ACT
    path = _write(detector, CENTRE_POLYNOMIALS, _metadata())

    # ASSERT
    with fits.open(path) as hdus:
        assert len(hdus) == 4
        assert hdus[3].header["EXTNAME"] == TRACES_EXTNAME
    assert list(_traces(path).columns) == ["ycoord", "order", "xcoord_centre"]


def test_order_centre_traces_match_the_unpacked_order_table(tmp_path: Path, log: object) -> None:
    # ARRANGE
    detector = _writer(tmp_path, log, recipeName="soxs-order-centre")
    path = _write(detector, CENTRE_POLYNOMIALS, _metadata())

    # ACT
    traces = _traces(path)
    _, expected, _ = unpack_order_table(log=log, orderTablePath=str(path), pixelDelta=1)

    # ASSERT
    assert len(traces) == len(expected) > 0
    for column in ["order", "ycoord", "xcoord_centre"]:
        np.testing.assert_allclose(traces[column], expected[column], rtol=1e-12, atol=1e-12)


def test_mflat_traces_add_both_edges_to_the_centre(tmp_path: Path, log: object) -> None:
    # ARRANGE
    detector = _writer(tmp_path, log, recipeName="soxs-mflat")

    # ACT
    path = _write(detector, {**CENTRE_POLYNOMIALS, **EDGE_POLYNOMIALS}, _metadata())

    # ASSERT
    traces = _traces(path)
    _, expected, _ = unpack_order_table(log=log, orderTablePath=str(path), pixelDelta=1)
    assert list(traces.columns) == ["ycoord", "order", "xcoord_centre", "xcoord_edgeup", "xcoord_edgelow"]
    for column in traces.columns:
        np.testing.assert_allclose(traces[column], expected[column], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(traces["xcoord_edgeup"] - traces["xcoord_centre"], 10.0)


def test_mflat_traces_are_in_unbinned_pixels_when_the_flat_was_binned(tmp_path: Path, log: object) -> None:
    # ARRANGE: THE FLAT IS BINNED 2 X 2 SO ITS METADATA LIMITS ARE IN BINNED ROWS
    binning = 2
    detector = _writer(tmp_path, log, recipeName="soxs-mflat", binx=binning, biny=binning)
    path = _write(detector, {**CENTRE_POLYNOMIALS, **EDGE_POLYNOMIALS}, _metadata(ymax=3.0))

    # ACT
    traces = _traces(path)
    _, binned, _ = unpack_order_table(
        log=log, orderTablePath=str(path), pixelDelta=1, binx=binning, biny=binning, prebinned=True
    )

    # ASSERT: EVERY UNBINNED ROW IS WRITTEN, AND THE BINNED UNPACK IS THE EVERY-SECOND-ROW SUBSET
    order10 = traces.loc[traces["order"] == 10]
    assert list(order10["ycoord"]) == [0, 1, 2, 3, 4, 5]
    onBinnedGrid = traces.loc[traces["ycoord"] % binning == 0]
    assert len(onBinnedGrid) == len(binned)
    np.testing.assert_allclose(onBinnedGrid["ycoord"] / binning, binned["ycoord"])
    for column in ["xcoord_centre", "xcoord_edgeup", "xcoord_edgelow"]:
        np.testing.assert_allclose(onBinnedGrid[column] / binning, binned[column], rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("recipeName", ["soxs-stare", "soxs-nod", "soxs-offset"])
def test_object_trace_tables_do_not_get_the_order_traces_extension(
    tmp_path: Path, log: object, recipeName: str
) -> None:
    # ARRANGE
    detector = _writer(tmp_path, log, recipeName=recipeName)

    # ACT
    path = _write(detector, CENTRE_POLYNOMIALS, _metadata())

    # ASSERT
    with fits.open(path) as hdus:
        assert len(hdus) == 3
        assert TRACES_EXTNAME not in [hdu.name for hdu in hdus]


def test_a_two_extension_order_table_still_unpacks(tmp_path: Path, log: object) -> None:
    # ARRANGE
    detector = _writer(tmp_path, log, recipeName="soxs-order-centre")
    path = _write(detector, CENTRE_POLYNOMIALS, _metadata())
    with fits.open(path) as hdus:
        hdus[:3].writeto(tmp_path / "legacy.fits")

    # ACT
    _, legacy, _ = unpack_order_table(log=log, orderTablePath=str(tmp_path / "legacy.fits"), pixelDelta=1)
    _, current, _ = unpack_order_table(log=log, orderTablePath=str(path), pixelDelta=1)

    # ASSERT
    pd.testing.assert_frame_equal(legacy, current)
