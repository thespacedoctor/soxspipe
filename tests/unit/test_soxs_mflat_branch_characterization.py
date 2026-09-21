"""Characterization of `soxs_mflat.produce_product` and `calibrate_frame_set` branches.

These tests pin what `produce_product` and `calibrate_frame_set` do today,
before a later commit splits `soxs_mflat.py`'s functions into smaller methods, including
three defects: an unbound local reused across loop iterations, a stale
order-table path leaking from one lamp into the next, and an instrument
branch that only ever binds `domeflatCollection` for SOXS. Each defect is
pinned with a docstring noting it is pinned as found, not as intended -- the
split must preserve this behaviour, not quietly fix it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes.soxs_mflat import soxs_mflat
from tests.factories import prepared_fits, synthetic_ccd
from tests.unit.test_soxs_mflat_characterization import (
    _new_bare_recipe,
    _stub_shared_collaborators,
    soxs_mflat_module,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# shared helpers -- produce_product
# ---------------------------------------------------------------------------


def _median_flux_df() -> pd.DataFrame:
    """Return a fresh, single-order median-flux table.

    A fresh instance is required on every call: the production code renames
    columns on this frame `inplace=True`, so a shared instance reused across
    lamps would carry a stale column name into the next lamp's iteration.
    """
    return pd.DataFrame({"order": [10], "medianFlux": [100.0]})


class _EmptyOrderTableCollection:
    """Order-table lookup double that never finds a match."""

    def filter(self, **_: object) -> _EmptyOrderTableCollection:
        return self

    def files_filtered(self, *, include_path: bool) -> list[str]:
        assert include_path is True
        return []


class _RecordingOrderTableCollection:
    """Order-table lookup double recording every filter dict it receives."""

    def __init__(self, orderTablePath: str, calls: list[dict[str, str]]) -> None:
        self._orderTablePath = orderTablePath
        self._calls = calls

    def filter(self, **filters: object) -> _RecordingOrderTableCollection:
        self._calls.append(dict(filters))
        return self

    def files_filtered(self, *, include_path: bool) -> list[str]:
        assert include_path is True
        return [self._orderTablePath]


class _FilesFilteredResult:
    """A single filtered order-table lookup result."""

    def __init__(self, files: list[str]) -> None:
        self._files = files

    def files_filtered(self, *, include_path: bool) -> list[str]:
        assert include_path is True
        return list(self._files)


class _TaggedOrderTableCollection:
    """Order-table lookup double keyed by the exact filter dict a lamp passes."""

    def __init__(self, filesByFilter: dict[tuple[tuple[str, str], ...], list[str]]) -> None:
        self._filesByFilter = filesByFilter

    def filter(self, **filters: object) -> _FilesFilteredResult:
        key = tuple(sorted(filters.items()))
        return _FilesFilteredResult(self._filesByFilter.get(key, []))


def _make_tagged_fake_edges(recipe: soxs_mflat, pathsByTag: dict[str, Path]) -> type:
    """Return a `detect_order_edges` double whose product path depends on the lamp tag."""

    class FakeEdges:
        def __init__(self, **kwargs: object) -> None:
            self._products = kwargs["productsTable"]
            self._qc = kwargs["qcTable"]
            self._tag = kwargs["tag"]

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            path = pathsByTag[self._tag]
            row = {
                "soxspipe_recipe": recipe.recipeName,
                "product_label": f"ORDER_LOC{self._tag}",
                "product_desc": "table of coefficients from polynomial fits to order locations",
                "file_name": Path(path).name,
                "file_type": "FITS",
                "obs_date_utc": "2024-01-02T03:04:05.678",
                "reduction_date_utc": "2024-01-02T04:05:06",
                "file_path": str(path),
                "label": "PROD",
            }
            products = pd.concat([self._products, pd.DataFrame([row])], ignore_index=True)
            return products, self._qc.copy(deep=True), pd.DataFrame({"order": [10], "count": [5]})

    return FakeEdges


def test_x_shooter_order_table_filter_carries_object_key_only_for_a_tagged_lamp(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The untagged lamp filters on `PRO_CATG` alone; a tagged X-Shooter lamp adds `OBJECT`."""
    # ARRANGE
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=901)
    flatFrame = synthetic_ccd(seed=902, prepared=True)
    dFrame = synthetic_ccd(seed=903, prepared=True)
    stitched = synthetic_ccd(seed=904, prepared=True)

    recipe = _new_bare_recipe(
        log,
        tmp_path,
        orderPath,
        subtractBackground=False,
        calibratedFlatFiles=["flat.fits"],
        dFlatFiles=["dflat.fits"],
        qFlatFiles=[],
        domeFlatFiles=[],
    )
    recipe.inst = "XSH"
    filterCalls: list[dict[str, str]] = []
    recipe.inputFrames = _RecordingOrderTableCollection(str(orderPath), filterCalls)

    monkeypatch.setattr(recipe, "calibrate_frame_set", lambda: ([flatFrame], [dFrame], [], []))
    monkeypatch.setattr(recipe, "normalise_flats", lambda cf, orderTablePath, **k: [cf[0].copy()])
    monkeypatch.setattr(recipe, "clip_and_stack", lambda **k: k["frames"][0].copy())
    monkeypatch.setattr(recipe, "mask_low_sens_pixels", lambda **k: (k["frame"].copy(), _median_flux_df()))
    monkeypatch.setattr(recipe, "_write", lambda *a, **k: str(tmp_path / "MFLAT.fits"))
    monkeypatch.setattr(recipe, "stitch_uv_mflats", lambda *a, **k: stitched)
    _stub_shared_collaborators(recipe, monkeypatch, orderPath)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert filterCalls[0] == {"PRO_CATG": "ORDER_TAB_VIS"}
    assert filterCalls[1] == {"PRO_CATG": "ORDER_TAB_VIS", "OBJECT": "LAMP,DORDERDEF"}


def test_missing_order_table_on_the_first_lamp_raises_unbound_local_error(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When the first lamp's order table is not found, `orderTablePath` is never bound.

    The `else` branch that handles "no order table found" only appends
    `None` to `self.orderTableSet`; it never assigns the local
    `orderTablePath` the next line reads. On the very first lamp this name
    has never been bound at all, so the call raises `UnboundLocalError`
    rather than a domain-meaningful error. Pinned as found, not as intended.
    """
    # ARRANGE
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=910)
    flatFrame = synthetic_ccd(seed=911, prepared=True)

    recipe = _new_bare_recipe(
        log,
        tmp_path,
        orderPath,
        subtractBackground=False,
        calibratedFlatFiles=["flat.fits"],
        dFlatFiles=[],
        qFlatFiles=[],
        domeFlatFiles=[],
    )
    recipe.inputFrames = _EmptyOrderTableCollection()
    monkeypatch.setattr(recipe, "calibrate_frame_set", lambda: ([flatFrame], [], [], []))
    _stub_shared_collaborators(recipe, monkeypatch, orderPath)

    # ACT / ASSERT
    with pytest.raises(UnboundLocalError):
        recipe.produce_product()
    assert recipe.orderTableSet == [None]


def test_a_missing_dlamp_order_table_reuses_the_plain_lamps_final_order_loc_path(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A missing D-lamp order table leaves `orderTablePath` stale from the plain lamp.

    `orderTablePath` is reassigned twice per successful lamp: once from the
    order-table lookup at the top of the loop, and again near the bottom from
    the just-detected `ORDER_LOC<tag>` row in `self.products`. When a lamp's
    lookup finds nothing, only `None` is appended to `self.orderTableSet`;
    the `orderTablePath` local itself is left untouched, so the *next* use of
    it -- the D-lamp's own call to `normalise_flats` -- silently reads the
    plain lamp's final detected order-location product path instead of its
    own. `self.orderTableSet` ends up with an extra `None` in it as a result,
    which shifts every later entry out of alignment with its lamp. Pinned as
    found, not as intended.
    """
    # ARRANGE
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=920)
    untaggedSource = tmp_path / "untagged_source.fits"
    qlampSource = tmp_path / "qlamp_source.fits"
    untaggedSource.write_bytes(b"")
    qlampSource.write_bytes(b"")

    plainFrame = synthetic_ccd(seed=921, prepared=True)
    dFrame = synthetic_ccd(seed=922, prepared=True)
    qFrame = synthetic_ccd(seed=923, prepared=True)
    stitched = synthetic_ccd(seed=924, prepared=True)

    recipe = _new_bare_recipe(
        log,
        tmp_path,
        orderPath,
        subtractBackground=False,
        calibratedFlatFiles=["flat.fits"],
        dFlatFiles=["dflat.fits"],
        qFlatFiles=["qflat.fits"],
        domeFlatFiles=[],
    )
    recipe.inst = "XSH"
    recipe.inputFrames = _TaggedOrderTableCollection(
        {
            (("PRO_CATG", "ORDER_TAB_VIS"),): [str(untaggedSource)],
            (("OBJECT", "LAMP,DORDERDEF"), ("PRO_CATG", "ORDER_TAB_VIS")): [],
            (("OBJECT", "LAMP,QORDERDEF"), ("PRO_CATG", "ORDER_TAB_VIS")): [str(qlampSource)],
        }
    )

    pathsByTag = {
        "": tmp_path / "ORDER_LOC.fits",
        "_DLAMP": tmp_path / "ORDER_LOC_DLAMP.fits",
        "_QLAMP": tmp_path / "ORDER_LOC_QLAMP.fits",
    }
    receivedOrderTablePaths: dict[str, str] = {}

    def fake_normalise(cf: list[CCDData], orderTablePath: str, lamp: str = "", **_: object) -> list[CCDData]:
        receivedOrderTablePaths[lamp] = orderTablePath
        return [cf[0].copy()]

    monkeypatch.setattr(recipe, "calibrate_frame_set", lambda: ([plainFrame], [dFrame], [qFrame], []))
    monkeypatch.setattr(recipe, "normalise_flats", fake_normalise)
    monkeypatch.setattr(recipe, "clip_and_stack", lambda **k: k["frames"][0].copy())
    monkeypatch.setattr(recipe, "mask_low_sens_pixels", lambda **k: (k["frame"].copy(), _median_flux_df()))
    monkeypatch.setattr(recipe, "_write", lambda *a, **k: str(tmp_path / "MFLAT.fits"))
    monkeypatch.setattr(recipe, "stitch_uv_mflats", lambda *a, **k: stitched)
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **k: None)
    monkeypatch.setattr(recipe, "report_output", lambda: recipe.qc)
    monkeypatch.setattr(recipe, "clean_up", lambda: None)
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **k: None)
    monkeypatch.setattr(
        soxs_mflat_module,
        "unpack_order_table",
        lambda **k: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(soxs_mflat_module, "generic_quality_checks", lambda **k: k["qcTable"])
    monkeypatch.setattr(soxs_mflat_module, "spectroscopic_image_quality_checks", lambda **k: k["qcTable"])
    monkeypatch.setattr(soxs_mflat_module, "detect_order_edges", _make_tagged_fake_edges(recipe, pathsByTag))

    # ACT
    recipe.produce_product()

    # ASSERT
    assert receivedOrderTablePaths[""] == str(untaggedSource)
    # THE D-LAMP RECEIVES THE PLAIN LAMP'S *FINAL* DETECTED ORDER_LOC PATH,
    # NOT THE ORIGINAL SOURCE ORDER TABLE PATH -- THE STALE REUSE.
    assert receivedOrderTablePaths["_DLAMP"] == str(pathsByTag[""])
    assert receivedOrderTablePaths["_QLAMP"] == str(qlampSource)
    assert recipe.orderTableSet == [
        str(pathsByTag[""]),
        None,
        str(pathsByTag["_DLAMP"]),
        str(pathsByTag["_QLAMP"]),
        None,
    ]


# ---------------------------------------------------------------------------
# shared helpers -- calibrate_frame_set
# ---------------------------------------------------------------------------


def _calib_frame(data: np.ndarray, *, dprType: str = "LAMP,D,Q") -> CCDData:
    """Return a mutable CCD frame with deterministic mask and uncertainty."""
    frame = CCDData(data, unit="electron")
    frame.mask = np.zeros(data.shape, dtype=bool)
    frame.uncertainty = StdDevUncertainty(np.ones(data.shape))
    frame.header["WIN_BINX"] = 1
    frame.header["WIN_BINY"] = 1
    frame.header["DPR_TYPE"] = dprType
    return frame


def _calib_recipe(log: Any, *, arm: str = "VIS", inst: str = "SOXS") -> soxs_mflat:
    """Return the small recipe state required by `calibrate_frame_set`."""
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.kw = lambda key: key
    recipe.arm = arm
    recipe.inst = inst
    recipe.detectorParams = {}
    return recipe


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
    """Return collections keyed by the exact filter contract a caller passes, recording every call."""

    def __init__(self, collections: dict[tuple[tuple[str, str], ...], _Collection]) -> None:
        self._collections = collections
        self.calls: list[dict[str, str]] = []

    def filter(self, **filters: str) -> _Collection:
        self.calls.append(dict(filters))
        return self._collections.get(tuple(sorted(filters.items())), _Collection([], []))


def _filters(**filters: str) -> tuple[tuple[str, str], ...]:
    """Make an order-independent collection lookup key."""
    return tuple(sorted(filters.items()))


def test_x_shooter_calibrate_frame_set_calibrates_lamp_flats_and_reports_no_dome_flats(
    log: Any,
) -> None:
    """X-Shooter has no dome flats, so the dome collection is empty and the lamp flats calibrate."""
    # ARRANGE
    recipe = _calib_recipe(log, arm="VIS", inst="XSH")
    bias = _calib_frame(np.full((3, 3), 1.0))
    bias.header["MJDOBS"] = 0.0
    flat = _calib_frame(np.full((3, 3), 5.0))
    flat.header["MJDOBS"] = 1.0
    inputFrames = _Collections(
        {
            _filters(PRO_CATG="MASTER_BIAS_VIS"): _Collection(["bias_pre.fits"], [bias]),
            _filters(PRO_CATG="MASTER_DARK_VIS"): _Collection([], []),
            _filters(DPR_TYPE="LAMP,FLAT", DPR_TECH="IMAGE"): _Collection([], []),
            _filters(DPR_TYPE="DARK", DPR_TECH="IMAGE"): _Collection([], []),
            _filters(DPR_TYPE="LAMP,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection(["flat_pre.fits"], [flat]),
            _filters(DPR_TYPE="LAMP,DFLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="LAMP,QFLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
        }
    )
    recipe.inputFrames = inputFrames
    recipe.detrend = lambda **k: k["inputFrame"].copy()

    # ACT
    calibrated, dcalibrated, qcalibrated, domecalibrated = recipe.calibrate_frame_set()

    # ASSERT
    assert len(calibrated) == 1
    assert dcalibrated == qcalibrated == domecalibrated == []
    assert recipe.calibratedFlatFiles == ["flat.fits"]
    assert recipe.domeFlatFiles == []

    # THE NO-MASTER-DARK, NON-SOXS FALLBACK FILTERS ON `LAMP,FLAT`/`IMAGE`.
    assert _filters(DPR_TYPE="LAMP,FLAT", DPR_TECH="IMAGE") in [_filters(**call) for call in inputFrames.calls]


def test_x_shooter_calibrate_frame_set_raises_file_not_found_when_no_flats_are_given(
    log: Any,
) -> None:
    """With no flat frames of any lamp, an X-Shooter call reports the missing input, not a name error."""
    # ARRANGE
    recipe = _calib_recipe(log, arm="VIS", inst="XSH")
    recipe.inputFrames = _Collections(
        {
            _filters(PRO_CATG="MASTER_BIAS_VIS"): _Collection([], []),
            _filters(PRO_CATG="MASTER_DARK_VIS"): _Collection([], []),
            _filters(DPR_TYPE="LAMP,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="LAMP,DFLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="LAMP,QFLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
        }
    )
    recipe.detrend = lambda **k: k["inputFrame"].copy()

    # ACT / ASSERT
    with pytest.raises(FileNotFoundError, match="needs flat-frames as input"):
        recipe.calibrate_frame_set()


def test_soxs_falls_back_to_the_nearest_raw_dark_when_no_master_or_off_lamp_dark_exists(
    log: Any,
) -> None:
    """With no master dark and no off-lamp dark, SOXS falls back to raw `DARK`/`IMAGE` frames."""
    # ARRANGE
    recipe = _calib_recipe(log, arm="NIR", inst="SOXS")
    darkEarly = _calib_frame(np.full((3, 3), 2.0))
    darkEarly.header["MJDOBS"] = 5.0
    darkLate = _calib_frame(np.full((3, 3), 3.0))
    darkLate.header["MJDOBS"] = 15.0
    flat = _calib_frame(np.full((3, 3), 10.0))
    flat.header["MJDOBS"] = 14.0
    recipe.inputFrames = _Collections(
        {
            _filters(DPR_TYPE="DARK", DPR_TECH="IMAGE"): _Collection(
                ["dark-a.fits", "dark-b.fits"], [darkEarly, darkLate]
            ),
            _filters(DPR_TYPE="LAMP,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection(["flat_pre.fits"], [flat]),
            _filters(DPR_TYPE="FLAT,LAMP", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(LAMP2="Deut_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(LAMP1="Qth_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="DOME,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
        }
    )
    detrendCalls: list[dict[str, Any]] = []

    def detrend(**kwargs: Any) -> CCDData:
        detrendCalls.append(kwargs)
        return kwargs["inputFrame"].subtract(kwargs["dark"])

    recipe.detrend = detrend

    # ACT
    calibrated, dcalibrated, qcalibrated, domecalibrated = recipe.calibrate_frame_set()

    # ASSERT
    assert dcalibrated == qcalibrated == domecalibrated == []
    assert len(calibrated) == 1
    np.testing.assert_allclose(calibrated[0].data, 7.0)
    assert len(detrendCalls) == 1
    assert detrendCalls[0]["master_bias"] is None
    assert detrendCalls[0]["dark"] is darkLate


def test_soxs_falls_back_to_bias_only_subtraction_when_no_dark_exists_at_all(
    log: Any,
) -> None:
    """With no master, off-lamp, or raw dark, `darkCollection` is `None` and only bias is subtracted."""
    # ARRANGE
    recipe = _calib_recipe(log, arm="NIR", inst="SOXS")
    bias = _calib_frame(np.full((3, 3), 1.0))
    bias.header["MJDOBS"] = 0.0
    flat = _calib_frame(np.full((3, 3), 10.0))
    flat.header["MJDOBS"] = 14.0
    recipe.inputFrames = _Collections(
        {
            _filters(PRO_CATG="MASTER_BIAS_NIR"): _Collection(["bias_pre.fits"], [bias]),
            _filters(DPR_TYPE="LAMP,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection(["flat_pre.fits"], [flat]),
            _filters(DPR_TYPE="FLAT,LAMP", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(LAMP2="Deut_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(LAMP1="Qth_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="DOME,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
        }
    )
    detrendCalls: list[dict[str, Any]] = []

    def detrend(**kwargs: Any) -> CCDData:
        detrendCalls.append(kwargs)
        return kwargs["inputFrame"].copy()

    recipe.detrend = detrend

    # ACT
    calibrated, dcalibrated, qcalibrated, domecalibrated = recipe.calibrate_frame_set()

    # ASSERT
    assert dcalibrated == qcalibrated == domecalibrated == []
    assert len(calibrated) == 1
    assert len(detrendCalls) == 1
    assert detrendCalls[0]["master_bias"] is bias
    assert detrendCalls[0]["dark"] is None


def test_nir_flat_lookup_falls_back_to_the_second_filter_when_the_first_finds_nothing(
    log: Any,
) -> None:
    """`FLAT,LAMP` is tried first; when it finds nothing, `LAMP,FLAT` is used instead."""
    # ARRANGE
    recipe = _calib_recipe(log, arm="NIR", inst="SOXS")
    bias = _calib_frame(np.full((3, 3), 1.0))
    bias.header["MJDOBS"] = 0.0
    flat = _calib_frame(np.full((3, 3), 10.0))
    flat.header["MJDOBS"] = 1.0
    recipe.inputFrames = _Collections(
        {
            _filters(PRO_CATG="MASTER_BIAS_NIR"): _Collection(["bias_pre.fits"], [bias]),
            _filters(DPR_TYPE="FLAT,LAMP", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="LAMP,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection(["flat_pre.fits"], [flat]),
            _filters(LAMP2="Deut_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(LAMP1="Qth_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="DOME,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
        }
    )
    recipe.detrend = lambda **k: k["inputFrame"].copy()

    # ACT
    calibrated, _dcalibrated, _qcalibrated, _domecalibrated = recipe.calibrate_frame_set()

    # ASSERT
    assert len(calibrated) == 1
    assert recipe.calibratedFlatFiles == ["flat.fits"]


def test_soxs_uvb_vis_bias_subtracts_dome_flats_into_the_fourth_return_list(
    log: Any,
) -> None:
    """With a bias, no dark, and dome flats, dome flats land in the fourth return value."""
    # ARRANGE
    recipe = _calib_recipe(log, arm="VIS", inst="SOXS")
    bias = _calib_frame(np.full((3, 3), 1.0))
    bias.header["MJDOBS"] = 0.0
    domeflat = _calib_frame(np.full((3, 3), 30.0))
    domeflat.header["MJDOBS"] = 1.0
    recipe.inputFrames = _Collections(
        {
            _filters(PRO_CATG="MASTER_BIAS_VIS"): _Collection(["bias_pre.fits"], [bias]),
            _filters(LAMP2="Deut_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(LAMP1="Qth_Lamp", DPR_TECH="ECHELLE,SLIT"): _Collection([], []),
            _filters(DPR_TYPE="DOME,FLAT", DPR_TECH="ECHELLE,SLIT"): _Collection(["dome_pre.fits"], [domeflat]),
        }
    )
    detrendCalls: list[dict[str, Any]] = []

    def detrend(**kwargs: Any) -> CCDData:
        detrendCalls.append(kwargs)
        return kwargs["inputFrame"].copy()

    recipe.detrend = detrend

    # ACT
    calibrated, dcalibrated, qcalibrated, domecalibrated = recipe.calibrate_frame_set()

    # ASSERT
    assert calibrated == []
    assert dcalibrated == []
    assert qcalibrated == []
    assert len(domecalibrated) == 1
    assert detrendCalls[0]["master_bias"] is bias
    assert detrendCalls[0]["dark"] is None
    assert recipe.domeFlatFiles == ["dome.fits"]
