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


def _write(
    detector: detect_continuum, polynomials: dict, metadata: pd.DataFrame, frameBin: tuple[int, int] = (1, 1)
) -> Path:
    """Write the order table from a frame carrying the given absolute binning."""
    header = instrument_header(arm="VIS")
    header["ESO DET BINX"], header["ESO DET BINY"] = frameBin
    frame = SimpleNamespace(header=header)
    return Path(
        detector.write_order_table_to_file(
            frame=frame, orderPolyTable=pd.DataFrame([polynomials]), orderMetaTable=metadata
        )
    )


def _metadata(ymax: float = 6.0, axis: str = "y") -> pd.DataFrame:
    return pd.DataFrame({"order": [10, 11], f"{axis}min": [0.0, 2.0], f"{axis}max": [ymax, ymax + 2.0]})


def _nir(polynomials: dict) -> dict:
    """Re-key a polynomial dictionary so the traces are solved for y and sampled along x."""
    return {k.replace("degy", "degx"): v for k, v in polynomials.items()}


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
    path = _write(detector, {**CENTRE_POLYNOMIALS, **EDGE_POLYNOMIALS}, _metadata(ymax=3.0), (binning, binning))

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
    # HAND-COMPUTED AT y = 3, order 10: 1 + 2 * 3 + 3 * 10 + 0.5 * 10 * 3 = 52, AND THE EDGES SIT 10 EITHER SIDE
    row = order10.loc[order10["ycoord"] == 3].iloc[0]
    assert (row["xcoord_centre"], row["xcoord_edgeup"], row["xcoord_edgelow"]) == pytest.approx((52.0, 62.0, 42.0))
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


def test_mflat_traces_are_unbinned_when_the_order_centre_table_was_binned_too(tmp_path: Path, log: object) -> None:
    # ARRANGE: THE FLAT AND THE ORDER-CENTRE FRAME ARE BOTH BINNED 2 X 2, SO THE MFLAT BINNING RATIO IS 1
    detector = _writer(tmp_path, log, recipeName="soxs-mflat", binx=1, biny=1)
    path = _write(detector, {**CENTRE_POLYNOMIALS, **EDGE_POLYNOMIALS}, _metadata(ymax=3.0), (2, 2))

    # ACT
    traces = _traces(path)

    # ASSERT: THE COEFFICIENTS ARE IN 2 X 2 BINNED PIXELS, THE LIMITS ARE IN BINNED ROWS, THE OUTPUT IN UNBINNED PIXELS
    order10 = traces.loc[traces["order"] == 10]
    assert list(order10["ycoord"]) == [0, 1, 2, 3, 4, 5]
    # AT UNBINNED y = 3 THE BINNED ROW IS 1.5: 1 + 2 * 1.5 + 30 + 0.5 * 10 * 1.5 = 41.5, WHICH IS 83 UNBINNED
    row = order10.loc[order10["ycoord"] == 3].iloc[0]
    assert row["xcoord_centre"] == pytest.approx(83.0)
    # THE UPPER EDGE IS 11 + 2 * 1.5 + 30 + 0.5 * 10 * 1.5 = 51.5 BINNED
    assert row["xcoord_edgeup"] == pytest.approx(103.0)


def test_order_centre_traces_are_unbinned_when_the_order_centre_frame_was_binned(tmp_path: Path, log: object) -> None:
    # ARRANGE: THE ORDER-CENTRE FRAME IS BINNED 2 X 2, SO ITS COEFFICIENTS AND LIMITS ARE IN BINNED PIXELS
    detector = _writer(tmp_path, log, recipeName="soxs-order-centre", binx=2, biny=2)
    path = _write(detector, CENTRE_POLYNOMIALS, _metadata(ymax=3.0), (2, 2))

    # ACT
    traces = _traces(path)

    # ASSERT
    order10 = traces.loc[traces["order"] == 10]
    assert list(order10["ycoord"]) == [0, 1, 2, 3, 4, 5]
    row = order10.loc[order10["ycoord"] == 3].iloc[0]
    assert row["xcoord_centre"] == pytest.approx(83.0)


def test_nir_traces_name_their_columns_after_the_swapped_axes(tmp_path: Path, log: object) -> None:
    # ARRANGE: NIR ORDERS RUN ALONG x, SO THE TRACES ARE SAMPLED ALONG x AND SOLVED FOR y
    detector = _writer(tmp_path, log, recipeName="soxs-mflat")
    polynomials = _nir({**CENTRE_POLYNOMIALS, **EDGE_POLYNOMIALS})
    path = _write(detector, polynomials, _metadata(axis="x"))

    # ACT
    traces = _traces(path)
    _, expected, _ = unpack_order_table(log=log, orderTablePath=str(path), pixelDelta=1)

    # ASSERT
    assert list(traces.columns) == ["xcoord", "order", "ycoord_centre", "ycoord_edgeup", "ycoord_edgelow"]
    for column in traces.columns:
        np.testing.assert_allclose(traces[column], expected[column], rtol=1e-12, atol=1e-12)


def test_nir_mflat_traces_scale_each_axis_by_its_own_binning(tmp_path: Path, log: object) -> None:
    # ARRANGE: AN MFLAT BINNED 2 ALONG x AND 1 ALONG y, WITH THE SAME BINNING AS THE ORDER-CENTRE FRAME'S UNBINNED SCALE
    detector = _writer(tmp_path, log, recipeName="soxs-mflat", binx=2, biny=1)
    polynomials = _nir({**CENTRE_POLYNOMIALS, **EDGE_POLYNOMIALS})
    path = _write(detector, polynomials, _metadata(ymax=3.0, axis="x"), (2, 1))

    # ACT
    traces = _traces(path)
    _, binned, _ = unpack_order_table(log=log, orderTablePath=str(path), pixelDelta=1, binx=2, biny=1, prebinned=True)

    # ASSERT: x IS SAMPLED AT EVERY UNBINNED COLUMN AND THE BINNED UNPACK IS EVERY SECOND ONE
    order10 = traces.loc[traces["order"] == 10]
    assert list(order10["xcoord"]) == [0, 1, 2, 3, 4, 5]
    onBinnedGrid = traces.loc[traces["xcoord"] % 2 == 0]
    np.testing.assert_allclose(onBinnedGrid["xcoord"] / 2, binned["xcoord"])
    for column in ["ycoord_centre", "ycoord_edgeup", "ycoord_edgelow"]:
        np.testing.assert_allclose(onBinnedGrid[column], binned[column], rtol=1e-12, atol=1e-12)
    # HAND-COMPUTED AT x = 3, order 10
    row = order10.loc[order10["xcoord"] == 3].iloc[0]
    assert row["ycoord_centre"] == pytest.approx(52.0)
