"""Synthetic contracts for master-flat normalization, masking, and stitching."""

from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes.soxs_mflat import soxs_mflat

pytestmark = pytest.mark.unit
mflatModule = importlib.import_module("soxspipe.recipes.soxs_mflat")


def _frame(data: np.ndarray, *, dprType: str = "LAMP,D,Q") -> CCDData:
    """Return a mutable CCD frame with deterministic mask and uncertainty."""
    frame = CCDData(data, unit="electron")
    frame.mask = np.zeros(data.shape, dtype=bool)
    frame.uncertainty = StdDevUncertainty(np.ones(data.shape))
    frame.header["WIN_BINX"] = 1
    frame.header["WIN_BINY"] = 1
    frame.header["DPR_TYPE"] = dprType
    return frame


def _recipe(log: object) -> soxs_mflat:
    """Return the small recipe state required by master-flat helper seams."""
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.kw = lambda key: key
    recipe.arm = "VIS"
    recipe.axisA = "x"
    recipe.axisB = "y"
    recipe.binx = 1
    recipe.biny = 1
    recipe.recipeName = "soxs-mflat"
    recipe.recipeSettings = {
        "centre-order-window": 2,
        "low-sensitivity-clipping-sigma": 2,
        "scale-d2-to-qth": True,
    }
    recipe.qc = pd.DataFrame()
    recipe.products = pd.DataFrame()
    recipe.dateObs = "2024-01-01T00:00:00"
    recipe.debug = False
    return recipe


def _centre_pixels() -> pd.DataFrame:
    """Return a five-row vertical order-centre trace."""
    return pd.DataFrame({"xcoord_centre": [2] * 5, "ycoord": range(5)})


class _Collection:
    """Small in-memory image collection with the recipe's public surface."""

    def __init__(self, files: list[str], frames: list[CCDData]) -> None:
        self.files = files
        self._frames = frames

    def ccds(self, *, ccd_kwargs: dict[str, str]) -> list[CCDData]:
        assert ccd_kwargs["hdu_uncertainty"] == "ERRS"
        return self._frames

    def headers(self) -> list[dict[str, float]]:
        return [{"MJDOBS": frame.header["MJDOBS"]} for frame in self._frames]


class _Collections:
    """Return collections keyed by the exact filter contract passed by mflat."""

    def __init__(self, collections: dict[tuple[tuple[str, str], ...], _Collection]) -> None:
        self._collections = collections

    def filter(self, **filters: str) -> _Collection:
        return self._collections.get(tuple(sorted(filters.items())), _Collection([], []))


def _filters(**filters: str) -> tuple[tuple[str, str], ...]:
    """Make an order-independent collection lookup key."""
    return tuple(sorted(filters.items()))


def test_calibrate_frame_set_bias_corrects_each_soxs_lamp_collection(
    log: object,
) -> None:
    """D, Q, and dome flats retain their collection and normalized filenames."""
    recipe = _recipe(log)
    recipe.inst = "SOXS"
    recipe.detectorParams = {}
    bias = _frame(np.full((3, 3), 2.0))
    dflat = _frame(np.full((3, 3), 10.0))
    qflat = _frame(np.full((3, 3), 20.0))
    domeflat = _frame(np.full((3, 3), 30.0))
    for index, frame in enumerate([bias, dflat, qflat, domeflat]):
        frame.header["MJDOBS"] = float(index)
    recipe.inputFrames = _Collections(
        {
            _filters(PRO_CATG="MASTER_BIAS_VIS"): _Collection(["bias_pre.fits"], [bias]),
            _filters(LAMP2="Deut_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection(
                ["dflat_pre.fits"], [dflat]
            ),
            _filters(LAMP1="Qth_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection(
                ["qflat_pre.fits"], [qflat]
            ),
            _filters(DPR_TYPE="DOME,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection(
                ["dome_pre.fits"], [domeflat]
            ),
        }
    )
    detrendCalls: list[dict[str, Any]] = []

    def detrend(**kwargs: Any) -> CCDData:
        detrendCalls.append(kwargs)
        return kwargs["inputFrame"].copy()

    recipe.detrend = detrend

    calibrated, dcalibrated, qcalibrated, domecalibrated = recipe.calibrate_frame_set()

    assert calibrated == []
    assert [frame.data[0, 0] for frame in dcalibrated] == [10.0]
    assert [frame.data[0, 0] for frame in qcalibrated] == [20.0]
    assert [frame.data[0, 0] for frame in domecalibrated] == [30.0]
    assert all(call["master_bias"] is bias and call["dark"] is None for call in detrendCalls)
    assert recipe.dFlatFiles == ["dflat.fits"]
    assert recipe.qFlatFiles == ["qflat.fits"]
    assert recipe.domeFlatFiles == ["dome.fits"]


def test_calibrate_frame_set_requires_at_least_one_flat_collection(log: object) -> None:
    """The recipe fails clearly rather than attempting to make an empty product."""
    recipe = _recipe(log)
    recipe.inst = "SOXS"
    recipe.detectorParams = {}
    recipe.inputFrames = _Collections({})

    with pytest.raises(FileNotFoundError, match="needs flat-frames"):
        recipe.calibrate_frame_set()


def test_calibrate_frame_set_selects_the_nearest_master_dark_for_nir_flats(
    log: object,
) -> None:
    """NIR flats use the master dark nearest to their observation time."""
    recipe = _recipe(log)
    recipe.arm = "NIR"
    recipe.inst = "SOXS"
    recipe.detectorParams = {}
    darkEarly = _frame(np.full((3, 3), 1.0))
    darkLate = _frame(np.full((3, 3), 2.0))
    flat = _frame(np.full((3, 3), 10.0))
    darkEarly.header["MJDOBS"] = 10.0
    darkLate.header["MJDOBS"] = 20.0
    flat.header["MJDOBS"] = 18.0
    recipe.inputFrames = _Collections(
        {
            _filters(PRO_CATG="MASTER_DARK_NIR"): _Collection(
                ["dark-a.fits", "dark-b.fits"], [darkEarly, darkLate]
            ),
            _filters(DPR_TYPE="FLAT,LAMP", DPR_TECH="ECHELLE,SLIT"): _Collection(
                ["flat_pre.fits"], [flat]
            ),
        }
    )
    detrendCalls: list[dict[str, Any]] = []
    def detrend(**kwargs: Any) -> CCDData:
        detrendCalls.append(kwargs)
        return kwargs["inputFrame"].subtract(kwargs["dark"])

    recipe.detrend = detrend

    calibrated, dcalibrated, qcalibrated, domecalibrated = recipe.calibrate_frame_set()

    assert len(calibrated) == 1
    assert dcalibrated == qcalibrated == domecalibrated == []
    np.testing.assert_allclose(calibrated[0].data, 8.0)
    assert len(detrendCalls) == 1
    assert all(call["master_bias"] is None for call in detrendCalls)
    assert all(call["dark"] is darkLate for call in detrendCalls)
    assert recipe.calibratedFlatFiles == ["flat.fits"]


def test_find_uvb_overlap_order_scales_at_the_last_qth_dominated_order(
    log: object,
) -> None:
    """UVB lamp stitching chooses the order immediately after QTH dominance."""
    recipe = _recipe(log)
    recipe.arm = "UVB"

    class OrderTables:
        """Serve the D2 and QTH order-definition files requested by the recipe."""

        def filter(self, **filters: str) -> object:
            objectName = filters["OBJECT"]
            path = "d2-orders.fits" if objectName == "LAMP,DORDERDEF" else "qth-orders.fits"
            return SimpleNamespace(files_filtered=lambda **_: [path])

    recipe.inputFrames = OrderTables()
    calls: list[str] = []

    def normalise(flats: list[CCDData], *, orderTablePath: str) -> tuple[list[CCDData], pd.DataFrame]:
        calls.append(orderTablePath)
        if orderTablePath == "d2-orders.fits":
            return flats, pd.DataFrame({"order": [10, 20], "90_perc": [20.0, 30.0]})
        return flats, pd.DataFrame({"order": [10, 20], "90_perc": [100.0, 120.0]})

    recipe.normalise_flats = normalise

    overlapOrder = recipe.find_uvb_overlap_order_and_scale([], [])

    assert overlapOrder == 21
    assert calls == ["d2-orders.fits", "qth-orders.fits"]


def test_normalise_flats_scales_each_frame_to_its_order_centre_mean(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    recipe = _recipe(log)
    orderTable = tmp_path / "orders.fits"
    fits.PrimaryHDU().writeto(orderTable)
    monkeypatch.setattr(
        mflatModule,
        "unpack_order_table",
        lambda **kwargs: (None, _centre_pixels(), None),
    )
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frames = [_frame(np.full((5, 5), 10.0)), _frame(np.full((5, 5), 20.0))]

    normalised = recipe.normalise_flats(frames, str(orderTable))

    for frame, expectedUncertainty in zip(normalised, [0.1, 0.05]):
        np.testing.assert_allclose(frame.data, 1.0)
        np.testing.assert_allclose(frame.uncertainty.array, expectedUncertainty)
    np.testing.assert_allclose(frames[0].data, 10.0)
    np.testing.assert_allclose(frames[1].data, 20.0)
    assert "ORDEXP10" in recipe.qc["qc_name"].tolist()
    assert (recipe.binRatioX, recipe.binRatioY) == (1, 1)


def test_mask_low_sensitivity_pixels_returns_order_medians(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _recipe(log)
    pixels = pd.DataFrame(
        {
            "order": [10],
            "xcoord_edgeup": [4],
            "xcoord_edgelow": [1],
            "ycoord": [2],
        }
    )
    monkeypatch.setattr(
        mflatModule,
        "unpack_order_table",
        lambda **kwargs: (None, pixels, None),
    )
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.full((5, 5), 10.0))
    frame.data[2, 2] = 0.0

    masked, medianFluxes = recipe.mask_low_sens_pixels(
        frame,
        "orders.fits",
        returnMedianOrderFlux=True,
        writeQC=False,
    )

    assert bool(masked.mask[2, 2])
    assert not bool(masked.mask[2, 1])
    assert medianFluxes.to_dict("records") == [{"order": 10, "medianFlux": 10.0}]


def test_stitch_uv_mflats_scales_d_lamp_and_uses_the_selected_order_edge(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _recipe(log)
    recipe.settings = {}
    recipe.sofName = "UVB"
    recipe.startNightDate = "2024-01-01"
    recipe.binRatioX = 1
    recipe.binRatioY = 1
    recipe.orderTableSet = ["unused.fits", "centres.fits"]
    recipe.masterFlatSet = [
        _frame(np.zeros((4, 6))),
        _frame(np.full((4, 6), 2.0)),
        _frame(np.full((4, 6), 10.0)),
    ]
    edgePixels = pd.DataFrame(
        {"order": [11], "xcoord_edgeup": [1], "ycoord": [1]}
    )
    monkeypatch.setattr(
        mflatModule,
        "unpack_order_table",
        lambda **kwargs: (None, edgePixels, None),
    )
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)

    class FakeEdges:
        """Return the one generated order-location product used by stitching."""

        def __init__(self, **kwargs: object) -> None:
            self.products = kwargs["productsTable"]
            self.qc = kwargs["qcTable"]

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
            products = pd.DataFrame(
                {"product_label": ["ORDER_LOC"], "file_path": ["edges.fits"]}
            )
            return products, self.qc, {"detected": 1}

    monkeypatch.setattr(mflatModule, "detect_order_edges", FakeEdges)
    captured: dict[str, object] = {}
    recipe.mask_low_sens_pixels = lambda frame, orderTablePath: (
        captured.update({"orderTablePath": orderTablePath}) or frame
    )
    orderFluxes = pd.DataFrame(
        {"order": [10, 11], "_QLAMP": [10.0, 5.0], "_DLAMP": [9.0, 10.0]}
    )

    stitched = recipe.stitch_uv_mflats(orderFluxes, "original.fits")

    np.testing.assert_allclose(stitched.data[1, :5], 2.0)
    np.testing.assert_allclose(stitched.data[1, 5:], 10.0)
    assert stitched.header["DPR_TYPE"] == "LAMP,,"
    assert captured == {"orderTablePath": "edges.fits"}
