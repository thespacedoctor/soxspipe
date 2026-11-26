"""The error plane of a UVB/VIS science frame uses the read noise of its own prescan (DY-2437).

The master bias read noise (`QC RAW RON`) is measured on the morning biases, and the night frames can have more.
A science recipe measures the prescan of each raw science frame before the frame is trimmed, builds that frame's
uncertainty with it, and keeps it in the frame's header as `QC SCI RON`. Everything later reads the value from the
frame it works on: the bias-pedestal correction, the sky fit and the sky-line signal to noise of the extraction. A
stack of frames carries the root mean square of its frames. Calibration frames, the NIR arm and X-Shooter are
unchanged.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.read_noise import SCIENCE_RON_CARD, ReadNoise
from soxspipe.recipes import base_recipe
from tests.factories import instrument_header, pipeline_settings
from tests.factories.prepared_frames import (
    KEYWORDS,
    TEST_CRYOSTAT_TEMP,
    TEST_GAIN,
    TEST_NIR_DETECTOR_TEMP,
    TEST_RON,
    TEST_VIS_DETECTOR_TEMP,
)
from tests.factories.prepared_frames import prepare_frame_recipe as _recipe

pytestmark = pytest.mark.unit

ROWS = 200
PRESCAN_COLUMNS = 30
SCIENCE_COLUMNS = 40
SKY_ADU = 100.0
PRESCAN_SIGMA_ADU = 5.0
# THE PRESCAN IS READ AS COLUMNS 4 TO 28 OF EVERY ROW
PRESCAN_WINDOW = {"rows": {"start": 0, "end": ROWS}, "columns": {"start": 4, "end": 28}}
MASTER_BIAS_RON = TEST_RON
PRESCAN_RON = PRESCAN_SIGMA_ADU * TEST_GAIN
TOLERANCE = 0.04
CARD_KEY = f"HIERARCH {SCIENCE_RON_CARD}"


def _raw_science_frame(
    tmpPath: Path, *, prescanSigma: float = PRESCAN_SIGMA_ADU, name: str = "raw.fits", dprType: str = "OBJECT,ASYNC"
) -> str:
    """A raw VIS frame in ADU: a flat science region beside a noisy prescan. The default DPR type is a real one."""
    rng = np.random.default_rng(7)
    data = np.full((ROWS, PRESCAN_COLUMNS + SCIENCE_COLUMNS), SKY_ADU, dtype=np.float32)
    data[:, 4:28] = rng.normal(590.0, prescanSigma, (ROWS, 24)) if prescanSigma else 590.0
    frame = CCDData(data, unit=u.adu)
    frame.header[KEYWORDS["DPR_TYPE"]] = dprType
    frame.header[KEYWORDS["OBS_NAME"]] = "synthetic"
    frame.header[KEYWORDS["CP_TEMP_C"]] = TEST_CRYOSTAT_TEMP
    frame.header[KEYWORDS["NIR_TEMP_K"]] = TEST_NIR_DETECTOR_TEMP
    frame.header[KEYWORDS["VIS_TEMP_C"]] = TEST_VIS_DETECTOR_TEMP
    path = tmpPath / name
    frame.write(str(path), overwrite=True)
    return str(path)


def _prescan_recipe(log: Any, tmpPath: Path, *, recipeName: str = "soxs-stare", arm: str = "VIS", window: Any = None):
    """A recipe for the 200 x 40 science pixels beside the 30-column prescan strip."""
    calibrations = tmpPath / "calibrations"
    calibrations.mkdir(parents=True, exist_ok=True)
    CCDData(np.zeros((ROWS, SCIENCE_COLUMNS), dtype=np.int16), unit=u.dimensionless_unscaled).write(
        str(calibrations / "bpm.fits"), overwrite=True
    )
    recipe = _recipe(log, tmpPath, recipeName=recipeName, arm=arm, recipeSettings={"use_flat": False})
    recipe.dateObs = None
    recipe.detectorParams["science-pixels"] = {
        "rows": {"start": 0, "end": ROWS},
        "columns": {"start": PRESCAN_COLUMNS, "end": PRESCAN_COLUMNS + SCIENCE_COLUMNS},
    }
    recipe.detectorParams["prescan-pixels"] = PRESCAN_WINDOW if window is None else window
    return recipe


def _prepared(recipe: Any, rawPath: str) -> CCDData:
    preparedPath = recipe._prepare_single_frame(rawPath)
    return CCDData.read(preparedPath, hdu=0, unit=u.electron, hdu_uncertainty="ERRS")


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def _expected_error(ron: float) -> float:
    return float(np.sqrt(SKY_ADU * TEST_GAIN + ron**2))


def test_a_science_frame_error_plane_uses_the_read_noise_of_its_prescan(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    errors = prepared.uncertainty.array
    assert errors.shape == (ROWS, SCIENCE_COLUMNS)
    assert errors[10, 10] == pytest.approx(_expected_error(PRESCAN_RON), rel=TOLERANCE)
    assert errors[10, 10] > _expected_error(MASTER_BIAS_RON)


def test_the_prepared_frame_keeps_its_read_noise_in_its_header(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert prepared.header[CARD_KEY] == pytest.approx(PRESCAN_RON, rel=TOLERANCE)


def test_a_standard_star_recipe_also_uses_the_prescan(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path, recipeName="soxs-nod-std")

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path, dprType="STD,FLUX"))

    # ASSERT
    assert prepared.uncertainty.array[10, 10] == pytest.approx(_expected_error(PRESCAN_RON), rel=TOLERANCE)


@pytest.mark.parametrize("dprType", ["BIAS", "DARK", "LAMP,FLAT"])
def test_a_calibration_frame_in_a_science_sof_keeps_the_master_bias_read_noise(
    log: Any, tmp_path: Path, dprType: str
) -> None:
    # ARRANGE: THE CALIBRATION FRAME HAS A PRESCAN SCATTER OF 8 ADU, FAR FROM THE 5 ADU OF A SCIENCE FRAME
    recipe = _prescan_recipe(log, tmp_path)
    rawPath = _raw_science_frame(tmp_path, prescanSigma=8.0, dprType=dprType)

    # ACT
    prepared = _prepared(recipe, rawPath)

    # ASSERT
    expected = MASTER_BIAS_RON if dprType == "BIAS" else _expected_error(MASTER_BIAS_RON)
    assert prepared.uncertainty.array[10, 10] == pytest.approx(expected, rel=1e-6)
    assert CARD_KEY not in prepared.header
    assert _warnings(log) == []


def test_a_calibration_recipe_keeps_the_master_bias_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path, recipeName="soxs-mflat")

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert prepared.uncertainty.array[10, 10] == pytest.approx(_expected_error(MASTER_BIAS_RON), rel=1e-6)
    assert CARD_KEY not in prepared.header


def test_a_nir_science_recipe_keeps_its_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path, arm="NIR")

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert prepared.uncertainty.array[10, 10] == pytest.approx(_expected_error(MASTER_BIAS_RON), rel=1e-6)
    assert CARD_KEY not in prepared.header
    assert _warnings(log) == []


def test_a_detector_with_no_prescan_window_keeps_the_master_bias_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    del recipe.detectorParams["prescan-pixels"]

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert prepared.uncertainty.array[10, 10] == pytest.approx(_expected_error(MASTER_BIAS_RON), rel=1e-6)
    assert _warnings(log) == []


def test_a_frame_with_no_prescan_falls_back_to_the_master_bias_with_a_warning(log: Any, tmp_path: Path) -> None:
    # ARRANGE: THE PRESCAN WINDOW LIES OUTSIDE THE FRAME
    window = {"rows": {"start": 0, "end": ROWS}, "columns": {"start": 200, "end": 240}}
    recipe = _prescan_recipe(log, tmp_path, window=window)

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path, name="raw_a.fits"))

    # ASSERT
    assert prepared.uncertainty.array[10, 10] == pytest.approx(_expected_error(MASTER_BIAS_RON), rel=1e-6)
    assert CARD_KEY not in prepared.header
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "raw_a.fits" in warnings[0]
    assert "outside the frame" in warnings[0]


def test_a_prescan_with_no_scatter_falls_back_to_the_master_bias_with_a_warning(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)

    # ACT
    prepared = _prepared(recipe, _raw_science_frame(tmp_path, prescanSigma=0.0))

    # ASSERT
    assert prepared.uncertainty.array[10, 10] == pytest.approx(_expected_error(MASTER_BIAS_RON), rel=1e-6)
    assert "not a positive number" in _warnings(log)[0]


def test_each_frame_keeps_its_own_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    quiet = _raw_science_frame(tmp_path, prescanSigma=2.0, name="quiet.fits")
    noisy = _raw_science_frame(tmp_path, prescanSigma=8.0, name="noisy.fits")

    # ACT
    quietFrame = _prepared(recipe, quiet)
    noisyFrame = _prepared(recipe, noisy)

    # ASSERT
    assert quietFrame.uncertainty.array[10, 10] == pytest.approx(_expected_error(2.0 * TEST_GAIN), rel=TOLERANCE)
    assert noisyFrame.uncertainty.array[10, 10] == pytest.approx(_expected_error(8.0 * TEST_GAIN), rel=TOLERANCE)
    assert quietFrame.header[CARD_KEY] < noisyFrame.header[CARD_KEY]


def test_a_frame_that_falls_back_after_a_measured_one_takes_the_master_bias_noise_not_the_measured_one(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE: THE SECOND FRAME'S PRESCAN CANNOT BE MEASURED. THE RECIPE'S READ NOISE MUST STILL BE THE MASTER BIAS'S
    recipe = _prescan_recipe(log, tmp_path)
    _prepared(recipe, _raw_science_frame(tmp_path, prescanSigma=8.0, name="measured.fits"))
    recipe.detectorParams["prescan-pixels"] = {"rows": {"start": 0, "end": ROWS}, "columns": {"start": 200, "end": 240}}

    # ACT
    second = _prepared(recipe, _raw_science_frame(tmp_path, name="second.fits"))

    # ASSERT
    assert second.uncertainty.array[10, 10] == pytest.approx(_expected_error(MASTER_BIAS_RON), rel=1e-6)
    assert recipe.detectorParams["ron"].value == MASTER_BIAS_RON


def _header_with_ron(ron: float | None) -> fits.Header:
    header = fits.Header()
    if ron is not None:
        header[CARD_KEY] = ron
    return header


def test_a_frame_read_noise_is_its_header_value_else_the_detector_value(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    withCard = SimpleNamespace(header=_header_with_ron(5.1))
    withoutCard = SimpleNamespace(header=_header_with_ron(None))

    # ACT / ASSERT
    assert recipe._frame_read_noise(withCard) == pytest.approx(5.1)
    assert recipe._frame_read_noise(withoutCard) == pytest.approx(MASTER_BIAS_RON)


def _detrend_recipe(log: Any) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: keyword
    recipe.detectorParams = {"ron": 3.8 * u.electron}
    recipe.recipeSettings = {"subtract_background": False}
    recipe.darkDetrendWarningIssued2 = False
    recipe.recipeName = "soxs-stare"
    recipe.settings = {}
    return recipe


def _prepared_science_frame(counts: float, ron: float | None) -> CCDData:
    """A prepared frame of uniform counts whose error plane is sqrt(counts + ron^2)."""
    shape = (4, 4)
    header = _header_with_ron(ron)
    header["DPR_TYPE"] = "OBJECT"
    return CCDData(
        np.full(shape, counts, dtype=np.float32),
        unit=u.electron,
        meta=header,
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(
            np.full(shape, np.sqrt(counts + (3.8 if ron is None else ron) ** 2)), unit=u.electron
        ),
    )


def test_the_bias_pedestal_correction_uses_the_read_noise_of_each_frame(log: Any) -> None:
    # ARRANGE: TWO FRAMES OF THE SAME COUNTS AND A 60 e- BIAS PEDESTAL, WITH 4 AND 8 e- OF READ NOISE
    recipe = _detrend_recipe(log)
    bias = CCDData(
        np.full((4, 4), 60.0, dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.zeros((4, 4)), unit=u.electron),
    )

    # ACT
    quiet = recipe.detrend(_prepared_science_frame(100.0, 4.0), master_bias=bias)
    noisy = recipe.detrend(_prepared_science_frame(100.0, 8.0), master_bias=bias)

    # ASSERT: THE PHOTON TERM IS RE-DERIVED FROM THE BIAS-SUBTRACTED 40 e-, WITH THE FRAME'S OWN READ NOISE
    assert quiet.uncertainty.array[0, 0] == pytest.approx(np.sqrt(40.0 + 4.0**2), rel=1e-5)
    assert noisy.uncertainty.array[0, 0] == pytest.approx(np.sqrt(40.0 + 8.0**2), rel=1e-5)


def test_a_frame_with_no_prescan_read_noise_is_corrected_with_the_detector_value(log: Any) -> None:
    # ARRANGE
    recipe = _detrend_recipe(log)
    bias = CCDData(
        np.full((4, 4), 60.0, dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.zeros((4, 4)), unit=u.electron),
    )

    # ACT
    corrected = recipe.detrend(_prepared_science_frame(100.0, None), master_bias=bias)

    # ASSERT
    assert corrected.uncertainty.array[0, 0] == pytest.approx(np.sqrt(40.0 + 3.8**2), rel=1e-5)


def _stacking_recipe(tmpPath: Path, log: Any, monkeypatch: pytest.MonkeyPatch) -> base_recipe:
    recipe = object.__new__(base_recipe)
    recipe.log = log
    recipe.settings = pipeline_settings(tmpPath)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.arm = "VIS"
    recipe.detectorParams = {"ron": 3.8 * u.electron}
    recipe.imageType = "OBJECT"
    recipe.recipeSettings = {
        "stacked-clipping-sigma": 3.0,
        "stacked-clipping-iterations": 2,
        "frame-clipping-sigma": 3.0,
        "frame-clipping-iterations": 2,
    }
    recipe.verbose = False
    recipe.debug = False
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    return recipe


def _stackable_frame(sigma: float, ron: float | None) -> CCDData:
    shape = (8, 8)
    header = instrument_header()
    if ron is not None:
        header[CARD_KEY] = ron
    return CCDData(
        np.full(shape, 10.0, dtype=np.float32),
        unit=u.electron,
        meta=header,
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(shape, sigma), unit=u.electron),
    )


def test_a_stack_carries_the_root_mean_square_read_noise_of_its_frames(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _stacking_recipe(tmp_path, log, monkeypatch)
    frames = [_stackable_frame(4.0, 4.0), _stackable_frame(8.0, 8.0)]

    # ACT
    stacked = recipe.clip_and_stack(frames, recipe="soxs_stare")

    # ASSERT: THE FRAMES KEEP THEIR OWN VALUES
    assert stacked.header[CARD_KEY] == pytest.approx(float(np.sqrt((16.0 + 64.0) / 2.0)))
    assert frames[0].header[CARD_KEY] == 4.0
    assert frames[1].header[CARD_KEY] == 8.0


def test_a_stack_counts_a_frame_with_no_prescan_read_noise_at_the_detector_value(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _stacking_recipe(tmp_path, log, monkeypatch)
    frames = [_stackable_frame(4.0, 6.0), _stackable_frame(4.0, None)]

    # ACT
    stacked = recipe.clip_and_stack(frames, recipe="soxs_stare")

    # ASSERT
    assert stacked.header[CARD_KEY] == pytest.approx(float(np.sqrt((36.0 + 3.8**2) / 2.0)))
    assert CARD_KEY not in frames[1].header


def test_a_stack_of_frames_with_no_prescan_read_noise_carries_none(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _stacking_recipe(tmp_path, log, monkeypatch)

    # ACT
    stacked = recipe.clip_and_stack([_stackable_frame(4.0, None), _stackable_frame(4.0, None)], recipe="soxs_stare")

    # ASSERT
    assert CARD_KEY not in stacked.header


def test_the_stack_error_of_frames_with_different_errors_is_their_mean_over_the_square_root_of_n(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _stacking_recipe(tmp_path, log, monkeypatch)
    frames = [_stackable_frame(4.0, 4.0), _stackable_frame(8.0, 8.0)]

    # ACT
    stacked = recipe.clip_and_stack(frames, recipe="soxs_stare")

    # ASSERT: (4 + 8) / 2 / sqrt(2). THE QUADRATURE SUM, sqrt(4^2 + 8^2) / 2, IS 4.47: 5 % HIGHER FOR THIS EXTREME PAIR.
    # THE MEAN IS KEPT, AS EVERY STACKED CALIBRATION PRODUCT USES IT
    assert stacked.uncertainty.array[0, 0] == pytest.approx(6.0 / np.sqrt(2.0), rel=1e-5)


def test_the_sky_fit_read_noise_adds_the_master_bias_noise_to_that_of_the_frame(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.readNoise = ReadNoise(3.83, 4, 3.83 * np.sqrt(1.25), "master bias")
    frame = SimpleNamespace(header=_header_with_ron(5.1))

    # ACT / ASSERT
    assert recipe._sky_fit_read_noise(frame) == pytest.approx(float(np.sqrt(5.1**2 + 3.83**2 / 4)))


def test_the_sky_fit_read_noise_of_a_frame_with_no_prescan_value_is_that_of_the_master_bias(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.readNoise = ReadNoise(3.83, 4, 3.83 * np.sqrt(1.25), "master bias")
    frame = SimpleNamespace(header=_header_with_ron(None))

    # ACT / ASSERT
    assert recipe._sky_fit_read_noise(frame) == pytest.approx(3.83 * np.sqrt(1.25))


def test_the_sky_fit_read_noise_without_a_master_bias_is_none_or_the_frame_value(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.readNoise = None

    # ACT / ASSERT
    assert recipe._sky_fit_read_noise(SimpleNamespace(header=_header_with_ron(None))) is None
    assert recipe._sky_fit_read_noise(SimpleNamespace(header=_header_with_ron(5.1))) == pytest.approx(5.1)


def test_the_extraction_read_noise_is_that_of_the_frame_it_extracts(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)

    # ACT / ASSERT
    assert recipe._extraction_read_noise(SimpleNamespace(header=_header_with_ron(5.1))) == pytest.approx(5.1)
    assert recipe._extraction_read_noise(SimpleNamespace(header=_header_with_ron(None))) == MASTER_BIAS_RON
    assert recipe._extraction_read_noise() == MASTER_BIAS_RON


def test_the_science_read_noise_of_a_frame_is_recorded_as_a_qc_that_goes_to_the_product_header(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]

    # ACT
    recipe._record_science_read_noise_qc(SimpleNamespace(header=_header_with_ron(5.1)))

    # ASSERT
    row = recipe.qc[recipe.qc["qc_name"] == "SCI RON"].iloc[0]
    assert row["qc_value"] == pytest.approx(5.1)
    assert row["qc_unit"] == "electrons"
    assert bool(row["to_header"]) is True


def test_a_later_frame_replaces_the_science_read_noise_qc_of_an_earlier_one(log: Any, tmp_path: Path) -> None:
    # ARRANGE: EACH PASS OF A PER-FRAME REDUCTION RECORDS THE VALUE OF ITS OWN FRAME
    recipe = _prescan_recipe(log, tmp_path)
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]

    # ACT
    recipe._record_science_read_noise_qc(SimpleNamespace(header=_header_with_ron(5.1)))
    recipe._record_science_read_noise_qc(SimpleNamespace(header=_header_with_ron(4.1)))

    # ASSERT
    assert list(recipe.qc["qc_value"]) == [4.1]


def test_a_frame_with_no_prescan_read_noise_records_no_qc(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]

    # ACT
    recipe._record_science_read_noise_qc(SimpleNamespace(header=_header_with_ron(None)))

    # ASSERT
    assert len(recipe.qc) == 0


def test_updating_the_keywords_of_a_frame_records_no_science_read_noise_qc(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE NODDING RECIPES CALL IT FOR THE A FRAME, THE B FRAME AND THE FIRST RAW FRAME IN TURN, SO IT MUST NOT
    # RECORD THE QC: THE LAST CALL WOULD OVERWRITE THE VALUE OF THE FRAME THE PRODUCT IS BUILT FROM
    recipe = _prescan_recipe(log, tmp_path)
    recipe.kw = lambda keyword: keyword
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]
    recipe.imageType = "OBJECT"
    recipe.rawFrames = SimpleNamespace(
        to_pandas=lambda: pd.DataFrame({"filename": ["a_pre.fits"], "TYPE": ["OBJECT"], "ARM": ["VIS"]})
    )
    for name in ("_stamp_raw_frame_records", "_stamp_calibration_frame_records", "_stamp_recipe_parameter_records"):
        monkeypatch.setattr(recipe, name, lambda *args, **kwargs: None)

    # ACT
    recipe.update_fits_keywords(frame=SimpleNamespace(header=_header_with_ron(5.1)))

    # ASSERT
    assert len(recipe.qc) == 0


def test_a_stare_pass_records_the_science_read_noise_of_the_stack_it_extracts(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    from soxspipe.recipes.soxs_stare import soxs_stare

    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.recipeName = "soxs-stare"
    recipe.dateObs = None
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]
    stack = SimpleNamespace(header={**_header_with_ron(5.1), "DATE_OBS": "2026-01-15T03:14:22"})

    class StopError(Exception):
        pass

    def stop(*args: object, **kwargs: object) -> None:
        raise StopError

    monkeypatch.setattr(recipe, "detrend", stop)

    # ACT
    with pytest.raises(StopError):
        recipe._extract_stack(stack, [stack], False, (None,) * 7)

    # ASSERT
    assert list(recipe.qc["qc_value"]) == [5.1]


def _nod_recipe(log: Any, recipeName: str = "soxs-nod") -> Any:
    from soxspipe.recipes.soxs_nod import soxs_nod

    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.arm = "VIS"
    recipe.recipeName = recipeName
    recipe.masterBias = False
    recipe.detectorParams = {"ron": 3.8 * u.electron}
    recipe.dateObs = None
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]
    return recipe


def _nod_frame(counts: float, ron: float | None) -> CCDData:
    shape = (4, 4)
    return CCDData(
        np.full(shape, counts, dtype=np.float32),
        unit=u.electron,
        meta=_header_with_ron(ron),
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(shape, 2.0), unit=u.electron),
    )


@pytest.mark.parametrize("recipeName", ["soxs-nod", "soxs-offset"])
def test_a_difference_frame_carries_the_quadrature_sum_of_the_read_noise_of_its_two_frames(
    log: Any, recipeName: str
) -> None:
    # ARRANGE
    recipe = _nod_recipe(log, recipeName)
    aFrame, bFrame = _nod_frame(100.0, 3.0), _nod_frame(90.0, 4.0)

    # ACT
    difference, reverse = recipe._difference_cycle_frames(aFrame, bFrame)

    # ASSERT: THE FRAMES KEEP THEIR OWN VALUES
    assert difference.header[CARD_KEY] == pytest.approx(5.0)
    assert aFrame.header[CARD_KEY] == 3.0
    assert bFrame.header[CARD_KEY] == 4.0
    if recipeName == "soxs-nod":
        assert reverse.header[CARD_KEY] == pytest.approx(5.0)


def test_a_difference_counts_a_frame_with_no_prescan_read_noise_at_the_detector_value(log: Any) -> None:
    # ARRANGE
    recipe = _nod_recipe(log)

    # ACT
    difference, _ = recipe._difference_cycle_frames(_nod_frame(100.0, 3.0), _nod_frame(90.0, None))

    # ASSERT
    assert difference.header[CARD_KEY] == pytest.approx(float(np.sqrt(3.0**2 + 3.8**2)))


def test_a_difference_of_frames_with_no_prescan_read_noise_carries_none(log: Any) -> None:
    # ARRANGE
    recipe = _nod_recipe(log)

    # ACT
    difference, _ = recipe._difference_cycle_frames(_nod_frame(100.0, None), _nod_frame(90.0, None))

    # ASSERT
    assert CARD_KEY not in difference.header


def test_a_nodding_cycle_records_the_read_noise_of_its_difference(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    recipe = _nod_recipe(log)
    recipe.recipeSettings = {"save_single_frame_extractions": False}
    seen: list[list[float]] = []

    def extract(*args: object, **kwargs: object) -> tuple[None, None, None, bool]:
        seen.append(list(recipe.qc["qc_value"]))
        return None, None, None, False

    for name in ("_write_cycle_difference_frames", "_run_cycle_quality_checks", "_check_host_contamination"):
        monkeypatch.setattr(recipe, name, lambda *args, **kwargs: None)
    monkeypatch.setattr(recipe, "_combine_unflattened_differences", lambda *args: None)
    monkeypatch.setattr(recipe, "_detrend_cycle_frames", lambda ab, ba, a, b, **kwargs: (ab, ba, a, b), raising=True)
    monkeypatch.setattr(recipe, "_combined_cycle_frames", lambda *args: None)
    monkeypatch.setattr(recipe, "_extract_cycle_spectra", extract)

    # ACT: TWO CYCLES. THE PRODUCTS OF EACH ARE WRITTEN WITH THE VALUE OF ITS OWN DIFFERENCE
    recipe.process_single_ab_nodding_cycle(_nod_frame(100.0, 3.0), _nod_frame(90.0, 4.0), 1, "order.fits")
    recipe.process_single_ab_nodding_cycle(_nod_frame(100.0, 6.0), _nod_frame(90.0, 8.0), 2, "order.fits")

    # ASSERT
    assert seen == [[pytest.approx(5.0)], [pytest.approx(10.0)]]
    assert list(recipe.qc["qc_name"]) == ["SCI RON"]


def test_the_stacked_spectrum_of_the_cycles_carries_the_root_mean_square_of_the_cycle_values(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: TWO CYCLES. THE FIRST RAW FRAME, WHICH THE SPECTRUM HEADER IS BUILT FROM, HAS 3 e-
    recipe = _nod_recipe(log)
    seenAtStack: list[list[float]] = []

    def process(aFrame: Any, bFrame: Any, **kwargs: Any) -> tuple[None, None, None]:
        recipe._record_science_read_noise_qc(recipe._difference_cycle_frames(aFrame, bFrame)[0])
        return None, None, None

    def stack(spectra: Any, orderJoins: Any = None) -> tuple[None, None]:
        seenAtStack.append(list(recipe.qc["qc_value"]))
        return None, None

    recipe.recipeSettings = {"use_flat": False}
    monkeypatch.setattr(recipe, "process_single_ab_nodding_cycle", process)
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda *args, **kwargs: None)
    monkeypatch.setattr(recipe, "_leave_untraced_frames_out", lambda *args, **kwargs: None)
    monkeypatch.setattr(recipe, "_record_host_contamination", lambda *args, **kwargs: None)
    monkeypatch.setattr(recipe, "stack_extractions", stack)
    framesA = [_nod_frame(100.0, 3.0), _nod_frame(100.0, 6.0)]
    framesB = [_nod_frame(90.0, 4.0), _nod_frame(90.0, 8.0)]
    for index, frame in enumerate(framesA + framesB):
        frame.header["MJD-OBS"] = 60000.0 + index

    # ACT
    recipe._reduce_each_nodding_cycle(framesA, framesB, ["a1", "a2"], ["b1", "b2"], "order.fits", False)

    # ASSERT: sqrt((5^2 + 10^2) / 2), NOT THE LAST CYCLE'S 10
    assert seenAtStack == [[pytest.approx(float(np.sqrt((25.0 + 100.0) / 2.0)))]]


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (5.1, 5.1),
        (np.float32(5.0), 5.0),
        (4, 4.0),
        ("5.1", 5.1),
        ("abc", None),
        (None, None),
        (float("nan"), None),
        (float("inf"), None),
    ],
)
def test_the_read_noise_card_takes_any_real_number_and_nothing_else(value: object, expected: float | None) -> None:
    # ARRANGE
    # A PLAIN MAPPING, SINCE A FITS HEADER REFUSES NaN AND INFINITY
    frame = SimpleNamespace(header={} if value is None else {SCIENCE_RON_CARD: value})

    # ACT
    result = base_recipe._prescan_read_noise_of(frame)

    # ASSERT
    assert result == pytest.approx(expected) if expected is not None else result is None


@pytest.mark.parametrize("value", [float("inf"), -1.0, 0.0, True, False])
def test_the_read_noise_card_ignores_infinity_zero_negative_and_boolean_values(value: object) -> None:
    # ARRANGE
    header = {SCIENCE_RON_CARD: value}

    # ACT / ASSERT
    assert base_recipe._prescan_read_noise_of(SimpleNamespace(header=header)) is None


def test_the_prescan_read_noise_card_survives_detrending(log: Any) -> None:
    # ARRANGE: THE CARD IS READ FROM THE HEADER OF A DETRENDED FRAME, SO THE HEADER MUST STAY AN astropy FITS HEADER
    recipe = _detrend_recipe(log)
    bias = CCDData(
        np.full((4, 4), 60.0, dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.zeros((4, 4)), unit=u.electron),
    )

    # ACT
    detrended = recipe.detrend(_prepared_science_frame(100.0, 5.1), master_bias=bias)

    # ASSERT
    assert isinstance(detrended.header, fits.Header)
    assert detrended.header[CARD_KEY] == pytest.approx(5.1)


def test_the_stacked_spectrum_of_the_offset_cycles_carries_the_root_mean_square_of_the_cycle_values(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    from soxspipe.recipes.soxs_offset import soxs_offset

    recipe = soxs_offset.__new__(soxs_offset)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.arm = "VIS"
    recipe.recipeName = "soxs-offset"
    recipe.masterBias = False
    recipe.detectorParams = {"ron": 3.8 * u.electron}
    recipe.dateObs = None
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]
    recipe.recipeSettings = {"use_flat": False}
    seenAtStack: list[list[float]] = []

    def process(aFrame: Any, bFrame: Any, **kwargs: Any) -> tuple[None, None, None]:
        recipe._record_science_read_noise_qc(recipe._difference_cycle_frames(aFrame, bFrame)[0])
        return None, None, None

    def stack(spectra: Any, orderJoins: Any = None) -> tuple[None, None]:
        seenAtStack.append(list(recipe.qc["qc_value"]))
        return None, None

    monkeypatch.setattr(recipe, "process_single_ab_nodding_cycle", process)
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda *args, **kwargs: None)
    monkeypatch.setattr(recipe, "_leave_untraced_frames_out", lambda *args, **kwargs: None)
    monkeypatch.setattr(recipe, "stack_extractions", stack)

    # ACT
    recipe._reduce_each_offset_cycle(
        [_nod_frame(100.0, 3.0), _nod_frame(100.0, 6.0)],
        [_nod_frame(90.0, 4.0), _nod_frame(90.0, 8.0)],
        ["on1", "on2"],
        ["off1", "off2"],
        "order.fits",
        False,
    )

    # ASSERT
    assert seenAtStack == [[pytest.approx(float(np.sqrt((25.0 + 100.0) / 2.0)))]]


def _stare_mode_nod_recipe(log: Any, monkeypatch: pytest.MonkeyPatch) -> tuple[Any, list[list[float]]]:
    """A stare recipe that combines nod positions, and the SCI RON QC it holds when the spectra are stacked."""
    from soxspipe.recipes.soxs_stare import soxs_stare

    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.recipeName = "soxs-stare"
    recipe.dateObs = None
    recipe.qc = recipe._empty_qc_and_product_tables(pd)[0]
    seenAtStack: list[list[float]] = []

    def stack(spectra: Any, orderJoins: Any = None, notFlattened: bool = False) -> tuple[None, None]:
        seenAtStack.append(list(recipe.qc["qc_value"]))
        return None, None

    monkeypatch.setattr(recipe, "stack_extractions", stack)
    monkeypatch.setattr(recipe, "_combined_notes", lambda *args: [])
    monkeypatch.setattr(recipe, "_combine_not_flattened", lambda *args: None)
    return recipe, seenAtStack


def _position_extraction(ron: float | None, *, traceFallback: bool = False) -> Any:
    from soxspipe.recipes.nod_stare_stacking import stack_extraction

    return stack_extraction(
        productPath=None,
        combinedObject=SimpleNamespace(header=_header_with_ron(ron)),
        mergedSpectrum="spectrum",
        orderJoins=None,
        extractionPath=False,
        responseInputs=None,
        notFlatPath=None,
        notFlatSpectrum=None,
        traceFallback=traceFallback,
        slitEdgeSpill=False,
        extractionNotes=[],
    )


def test_the_combined_nod_spectrum_in_stare_mode_carries_the_root_mean_square_of_its_positions(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE LAST POSITION (B) HAS 8 e- AND THE FIRST (A) HAS 4 e-
    recipe, seenAtStack = _stare_mode_nod_recipe(log, monkeypatch)
    recipe._record_science_read_noise_qc(SimpleNamespace(header=_header_with_ron(8.0)))

    # ACT
    recipe._combine_extractions({"A": _position_extraction(4.0), "B": _position_extraction(8.0)})

    # ASSERT
    assert seenAtStack == [[pytest.approx(float(np.sqrt((16.0 + 64.0) / 2.0)))]]


def test_the_combined_nod_spectrum_leaves_out_a_position_with_no_object_trace(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, seenAtStack = _stare_mode_nod_recipe(log, monkeypatch)

    # ACT
    recipe._combine_extractions({"A": _position_extraction(4.0), "B": _position_extraction(8.0, traceFallback=True)})

    # ASSERT
    assert seenAtStack == [[pytest.approx(4.0)]]


def test_the_combined_nod_spectrum_skips_a_position_with_no_prescan_read_noise(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, seenAtStack = _stare_mode_nod_recipe(log, monkeypatch)

    # ACT
    recipe._combine_extractions({"A": _position_extraction(None), "B": _position_extraction(8.0)})

    # ASSERT
    assert seenAtStack == [[pytest.approx(8.0)]]


def test_the_combined_nod_spectrum_keeps_the_qc_when_no_position_has_a_prescan_read_noise(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, seenAtStack = _stare_mode_nod_recipe(log, monkeypatch)
    recipe._record_science_read_noise_qc(SimpleNamespace(header=_header_with_ron(5.1)))

    # ACT
    recipe._combine_extractions({"A": _position_extraction(None), "B": _position_extraction(None)})

    # ASSERT
    assert seenAtStack == [[pytest.approx(5.1)]]
