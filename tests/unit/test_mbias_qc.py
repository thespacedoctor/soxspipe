"""Characterization tests for master-bias quality-control calculations."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes.soxs_mbias import soxs_mbias
from tests.factories import product_table

pytestmark = pytest.mark.unit

# THE PACKAGE RE-EXPORTS THE CLASS UNDER THE MODULE NAME, SO FETCH THE MODULE EXPLICITLY
mbiasModule = import_module("soxspipe.recipes.soxs_mbias")


class _Frames:
    """Expose the image-collection surface needed by periodic-noise QC."""

    def __init__(self, frames: list[CCDData]) -> None:
        self._frames = frames

    def ccds(self, *, ccd_kwargs: dict[str, str]) -> list[CCDData]:
        assert ccd_kwargs["hdu_uncertainty"] == "ERRS"
        assert ccd_kwargs["hdu_mask"] == "QUAL"
        return self._frames


def _recipe(log: object) -> soxs_mbias:
    """Return the minimal mutable state consumed by the QC methods."""
    recipe = soxs_mbias.__new__(soxs_mbias)
    recipe.log = log
    recipe.qc = pd.DataFrame()
    recipe.recipeName = "soxs-mbias"
    recipe.dateObs = "2024-01-02T03:04:05"
    return recipe


def test_qc_bias_structure_records_axis_slopes(log: object) -> None:
    """Collapsed linear bias gradients become named STRUCT QC measurements."""
    recipe = _recipe(log)
    yPixels, xPixels = np.indices((4, 5))
    frame = 2.0 * xPixels + 3.0 * yPixels

    structX, structY = recipe.qc_bias_structure(frame)

    expectedX = np.polyfit(
        np.linspace(0, frame.shape[1], frame.shape[1], dtype=int),
        np.nansum(frame, axis=0),
        deg=1,
    )[0]
    expectedY = np.polyfit(
        np.linspace(0, frame.shape[0], frame.shape[0], dtype=int),
        np.nansum(frame, axis=1),
        deg=1,
    )[0]
    assert structX == pytest.approx(expectedX)
    assert structY == pytest.approx(expectedY)
    assert recipe.qc["qc_name"].tolist() == ["STRUCTX", "STRUCTY"]
    assert recipe.qc["qc_value"].tolist() == pytest.approx([structX, structY])
    assert recipe.qc["to_header"].tolist() == [True, True]


def test_qc_periodic_pattern_noise_records_maximum_frame_ratio(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """FFT-derived periodic-noise QC is finite and retains its header contract."""
    recipe = _recipe(log)
    random = np.random.default_rng(42)
    first = CCDData(random.normal(size=(20, 20)), unit="electron")
    second = CCDData(random.normal(size=(20, 20)), unit="electron")
    first.mask = np.zeros(first.shape, dtype=bool)
    second.mask = np.zeros(second.shape, dtype=bool)
    quicklookCalls: list[object] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.quicklook_image",
        lambda **kwargs: quicklookCalls.append(kwargs["CCDObject"]),
    )

    periodicNoise = recipe.qc_periodic_pattern_noise(_Frames([first, second]))

    assert np.isfinite(periodicNoise)
    assert periodicNoise > 0
    assert len(quicklookCalls) == 4
    assert recipe.qc["qc_name"].tolist() == ["FPN FRACMAX"]
    assert recipe.qc.loc[0, "qc_value"] == pytest.approx(periodicNoise)
    # NUMPY BOOL, NOT PYTHON BOOL, NOW THE COLUMN CARRIES A REAL BOOL DTYPE
    assert bool(recipe.qc.loc[0, "to_header"]) is True


def _summary(**overrides: object) -> dict[str, object]:
    """Return a bias-distribution summary for a symmetric raw frame and a tight master."""
    random = np.random.default_rng(3)
    arguments = {
        "rawPixels": random.normal(loc=100.0, scale=2.0, size=5000),
        "frameRon": 2.0,
        "masterPixels": random.normal(loc=100.0, scale=0.5, size=5000),
        "rawRon": 2.5,
        "masterRon": 0.5,
        "frameCount": 7,
        "clipSigma": 3,
        "clipIterations": 4,
    }
    arguments.update(overrides)
    return mbiasModule.bias_distribution_summary(**arguments)


def test_bias_distribution_summary_bins_span_five_sigma_around_raw_mean() -> None:
    """The histogram range is the raw mean plus or minus five frame RONs."""
    # ARRANGE
    rawPixels = np.array([98.0, 100.0, 102.0])

    # ACT
    summary = _summary(rawPixels=rawPixels, frameRon=2.0)

    # ASSERT
    binEdges = summary["binEdges"]
    assert len(binEdges) == mbiasModule.BIAS_HISTOGRAM_BINS + 1
    assert binEdges[0] == pytest.approx(100.0 - 5 * 2.0)
    assert binEdges[-1] == pytest.approx(100.0 + 5 * 2.0)


def test_bias_distribution_summary_reports_distinct_mean_and_median_for_skewed_pixels() -> None:
    """A skewed pixel population gives different mean and median levels."""
    # ARRANGE
    rawPixels = np.array([100.0, 100.0, 100.0, 100.0, 110.0])
    masterPixels = np.array([50.0, 60.0, 60.0, 60.0, 60.0])

    # ACT
    summary = _summary(rawPixels=rawPixels, masterPixels=masterPixels)

    # ASSERT
    assert summary["rawMean"] == pytest.approx(102.0)
    assert summary["rawMedian"] == pytest.approx(100.0)
    assert summary["masterMean"] == pytest.approx(58.0)
    assert summary["masterMedian"] == pytest.approx(60.0)
    annotation = "\n".join(summary["annotation"])
    assert "102.00" in annotation
    assert "58.00" in annotation


def test_bias_distribution_summary_labels_three_ron_values_distinctly() -> None:
    """The plotted-frame RON, the mean raw RON and the master RON each get their own line."""
    # ARRANGE
    expectedFrameRon = "2.13"
    expectedRawRon = "3.57"
    expectedMasterRon = "0.41"

    # ACT
    summary = _summary(frameRon=2.134, rawRon=3.568, masterRon=0.413, frameCount=7, clipSigma=3, clipIterations=4)

    # ASSERT
    lines = summary["annotation"]
    frameLines = [line for line in lines if "plotted raw frame" in line]
    rawLines = [line for line in lines if "RAW RON (QC)" in line]
    masterLines = [line for line in lines if "MASTER RON (QC)" in line]
    assert len(frameLines) == len(rawLines) == len(masterLines) == 1
    assert expectedFrameRon in frameLines[0]
    assert "RAW RON" not in frameLines[0]
    assert "3" in frameLines[0]
    assert "4 iterations" in frameLines[0]
    assert expectedRawRon in rawLines[0]
    assert "7" in rawLines[0]
    assert expectedMasterRon in masterLines[0]


def test_unclipped_raw_pixels_excludes_both_masks_and_restores_bias_level() -> None:
    """Clipped and bad pixels are dropped and the subtracted mean level is added back."""
    # ARRANGE
    noiseFrame = CCDData(np.array([[-1.0, 0.0], [1.0, 2.0]]), unit=u.electron)
    noiseFrame.mask = np.array([[True, False], [False, False]])
    badPixelMask = np.array([[False, False], [True, False]])

    # ACT
    pixels = mbiasModule._unclipped_raw_pixels(noiseFrame, badPixelMask, 100.0)

    # ASSERT
    assert pixels.tolist() == [100.0, 102.0]


def test_bias_distribution_summary_falls_back_to_data_range_when_frame_ron_is_zero() -> None:
    """A zero frame RON still gives finite, increasing bin edges that cover every pixel."""
    # ARRANGE
    rawPixels = np.array([99.0, 100.0, 101.0])
    masterPixels = np.array([100.0, 100.5])

    # ACT
    summary = _summary(rawPixels=rawPixels, masterPixels=masterPixels, frameRon=0.0)

    # ASSERT
    binEdges = summary["binEdges"]
    assert np.all(np.isfinite(binEdges))
    assert np.all(np.diff(binEdges) > 0)
    assert binEdges[0] == pytest.approx(99.0)
    assert binEdges[-1] == pytest.approx(101.0)


def test_bias_distribution_summary_falls_back_to_data_range_when_frame_ron_is_nan() -> None:
    """A non-finite frame RON still gives finite, increasing bin edges."""
    # ARRANGE
    rawPixels = np.array([100.0, 100.0])
    masterPixels = np.array([100.0])

    # ACT
    summary = _summary(rawPixels=rawPixels, masterPixels=masterPixels, frameRon=float("nan"))

    # ASSERT
    binEdges = summary["binEdges"]
    assert np.all(np.isfinite(binEdges))
    assert np.all(np.diff(binEdges) > 0)
    assert binEdges[0] == pytest.approx(100.0 - mbiasModule.BIAS_HISTOGRAM_FLAT_HALF_WIDTH)
    assert binEdges[-1] == pytest.approx(100.0 + mbiasModule.BIAS_HISTOGRAM_FLAT_HALF_WIDTH)


@pytest.mark.parametrize("emptySample", ["rawPixels", "masterPixels"])
def test_bias_distribution_summary_rejects_a_sample_with_no_pixels(emptySample: str) -> None:
    """A sample with no surviving pixels fails with a clear message instead of plotting NaN."""
    # ARRANGE
    overrides = {emptySample: np.array([])}

    # ACT / ASSERT
    with pytest.raises(ValueError, match="no finite, unmasked pixels"):
        _summary(**overrides)


def test_finite_unmasked_pixels_drops_masked_nan_and_infinite_pixels() -> None:
    """Masked pixels from every mask, NaN and infinite values are all left out."""
    # ARRANGE
    data = np.array([[1.0, np.nan], [np.inf, 4.0], [5.0, 6.0]])
    firstMask = np.array([[False, False], [False, False], [True, False]])
    secondMask = np.array([[False, False], [False, True], [False, False]])

    # ACT
    pixels = mbiasModule._finite_unmasked_pixels(data, firstMask, None, secondMask)

    # ASSERT
    assert pixels.tolist() == [1.0, 6.0]


def test_unclipped_raw_pixels_drops_non_finite_pixels() -> None:
    """A NaN or infinite raw pixel is not passed to the plot statistics."""
    # ARRANGE
    noiseFrame = CCDData(np.array([[np.nan, 0.0], [np.inf, 2.0]]), unit=u.electron)

    # ACT
    pixels = mbiasModule._unclipped_raw_pixels(noiseFrame, None, 10.0)

    # ASSERT
    assert pixels.tolist() == [10.0, 12.0]


def test_unclipped_raw_pixels_keeps_everything_when_no_masks_exist() -> None:
    """A missing clip mask and a missing bad-pixel mask exclude nothing."""
    # ARRANGE
    noiseFrame = CCDData(np.array([[-1.0, 0.0], [1.0, 2.0]]), unit=u.electron)

    # ACT
    pixels = mbiasModule._unclipped_raw_pixels(noiseFrame, None, 10.0)

    # ASSERT
    assert pixels.tolist() == [9.0, 10.0, 11.0, 12.0]


def _combine_fixture(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[soxs_mbias, list[CCDData], list[np.ndarray]]:
    """Return a recipe whose collaborators mimic the in-place mask replacement of the real ones."""
    first = CCDData(np.arange(100.0, 109.0).reshape(3, 3), unit=u.electron)
    second = CCDData(np.arange(200.0, 209.0).reshape(3, 3), unit=u.electron)
    first.mask = np.zeros((3, 3), dtype=bool)
    first.mask[0, 0] = True
    originalData = [first.data.copy(), second.data.copy()]
    levels = {id(first): 100.0, id(second): 200.0}
    rons = {id(first): 2.0, id(second): 4.0}
    clipMask = np.zeros((3, 3), dtype=bool)
    clipMask[1, 1] = True

    def fake_subtract(frame: CCDData) -> tuple[float, float, CCDData]:
        frame.mask = clipMask.copy()
        frame.data = frame.data - levels[id(frame)]
        return levels[id(frame)], rons[id(frame)], frame

    def fake_stack(**kwargs: object) -> CCDData:
        return CCDData(
            np.zeros((3, 3)),
            unit=u.electron,
            mask=np.zeros((3, 3), dtype=bool),
            uncertainty=StdDevUncertainty(np.ones((3, 3)), unit=u.electron),
        )

    recipe = _recipe(log)
    recipe.inputFrames = _Frames([first, second])
    monkeypatch.setattr(recipe, "subtract_mean_flux_level", fake_subtract)
    monkeypatch.setattr(recipe, "clip_and_stack", fake_stack)
    return recipe, [first, second], originalData


def test_combine_bias_frames_samples_first_frame_without_clipped_or_bad_pixels(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The sample comes from frame 0 and drops pixels bad before the mask was replaced."""
    # ARRANGE
    recipe, _frames, originalData = _combine_fixture(log, monkeypatch)
    expectedPixels = originalData[0].ravel()[[1, 2, 3, 5, 6, 7, 8]]

    # ACT
    result = recipe._combine_bias_frames()

    # ASSERT
    assert len(result) == 5
    rawFrameSample = result[4]
    assert rawFrameSample["pixels"].tolist() == expectedPixels.tolist()
    assert rawFrameSample["pixels"].max() < 200.0


def test_combine_bias_frames_keeps_mean_raw_ron_as_third_value(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The sample RON is frame 0's while the returned raw RON stays the mean of all frames."""
    # ARRANGE
    recipe, _frames, _originalData = _combine_fixture(log, monkeypatch)

    # ACT
    result = recipe._combine_bias_frames()

    # ASSERT
    assert result[2] == pytest.approx(3.0)
    assert result[4]["frameRon"] == pytest.approx(2.0)
    assert result[4]["frameCount"] == 2


def _plot_recipe(log: object, tmp_path: Path) -> soxs_mbias:
    """Return a recipe holding the state the QC plot method reads."""
    recipe = _recipe(log)
    recipe.arm = "VIS"
    recipe.qcDir = str(tmp_path)
    recipe.products = product_table().iloc[0:0].copy()
    recipe.recipeSettings = {"frame-clipping-sigma": 3, "frame-clipping-iterations": 4}
    return recipe


def _plot_inputs() -> tuple[dict[str, object], CCDData]:
    """Return a raw-frame sample and a master frame with a masked pixel."""
    random = np.random.default_rng(5)
    rawFrameSample = {
        "pixels": random.normal(loc=100.0, scale=2.0, size=4000),
        "frameRon": 2.0,
        "frameCount": 5,
    }
    masterFrame = CCDData(random.normal(loc=100.0, scale=0.5, size=(20, 20)), unit=u.electron)
    masterFrame.mask = np.zeros((20, 20), dtype=bool)
    masterFrame.mask[0, 0] = True
    return rawFrameSample, masterFrame


def test_plot_bias_distribution_qc_writes_pdf_in_qc_directory(log: object, tmp_path: Path) -> None:
    """The plot is a PDF named after the master bias, inside the QC directory."""
    # ARRANGE
    recipe = _plot_recipe(log, tmp_path)
    rawFrameSample, masterFrame = _plot_inputs()

    # ACT
    filePath = recipe.plot_bias_distribution_qc(rawFrameSample, masterFrame, 2.5, 0.5, "/products/MASTER_BIAS_VIS.fits")

    # ASSERT
    assert filePath == str(tmp_path / "MASTER_BIAS_VIS_BIAS_DISTRIBUTION_QC_PLOT.pdf")
    with open(filePath, "rb") as pdfFile:
        assert pdfFile.read(4) == b"%PDF"


def test_plot_bias_distribution_qc_registers_qc_pdf_product(log: object, tmp_path: Path) -> None:
    """The plot is recorded as a QC-labelled PDF product row of the mbias recipe."""
    # ARRANGE
    recipe = _plot_recipe(log, tmp_path)
    rawFrameSample, masterFrame = _plot_inputs()

    # ACT
    filePath = recipe.plot_bias_distribution_qc(rawFrameSample, masterFrame, 2.5, 0.5, "/products/MASTER_BIAS_VIS.fits")

    # ASSERT
    assert len(recipe.products) == 1
    row = recipe.products.iloc[0]
    assert row["label"] == "QC"
    assert row["file_type"] == "PDF"
    assert row["product_label"] == "BIAS_DISTRIBUTION_QC_PLOT"
    assert row["soxspipe_recipe"] == "soxs-mbias"
    assert row["file_name"] == "MASTER_BIAS_VIS_BIAS_DISTRIBUTION_QC_PLOT.pdf"
    assert row["file_path"] == filePath
    assert row["obs_date_utc"] == "2024-01-02T03:04:05"


def test_plot_bias_distribution_qc_leaves_no_open_figures(log: object, tmp_path: Path) -> None:
    """The figure is closed after saving so batch runs do not leak memory."""
    # ARRANGE
    recipe = _plot_recipe(log, tmp_path)
    rawFrameSample, masterFrame = _plot_inputs()
    openFigures = plt.get_fignums()

    # ACT
    recipe.plot_bias_distribution_qc(rawFrameSample, masterFrame, 2.5, 0.5, "/products/MASTER_BIAS_VIS.fits")

    # ASSERT
    assert plt.get_fignums() == openFigures


def test_plot_bias_distribution_qc_closes_figure_when_saving_fails(log: object, tmp_path: Path) -> None:
    """A failed save still closes the figure and records no product."""
    # ARRANGE
    recipe = _plot_recipe(log, tmp_path)
    recipe.qcDir = str(tmp_path / "missing")
    rawFrameSample, masterFrame = _plot_inputs()
    openFigures = plt.get_fignums()

    # ACT
    with pytest.raises(FileNotFoundError):
        recipe.plot_bias_distribution_qc(rawFrameSample, masterFrame, 2.5, 0.5, "/products/MASTER_BIAS_VIS.fits")

    # ASSERT
    assert plt.get_fignums() == openFigures
    assert recipe.products.empty


def test_plot_bias_distribution_qc_accepts_master_frame_without_mask(log: object, tmp_path: Path) -> None:
    """A master frame with no mask plots every pixel."""
    # ARRANGE
    recipe = _plot_recipe(log, tmp_path)
    rawFrameSample, masterFrame = _plot_inputs()
    masterFrame.mask = None

    # ACT
    filePath = recipe.plot_bias_distribution_qc(rawFrameSample, masterFrame, 2.5, 0.5, "/products/MASTER_BIAS_VIS.fits")

    # ASSERT
    assert Path(filePath).stat().st_size > 0


def test_plot_bias_distribution_qc_leaves_out_unmasked_nan_master_pixels(log: object, tmp_path: Path) -> None:
    """An unmasked NaN in the master bias does not turn the plotted master statistics into NaN."""
    # ARRANGE
    recipe = _plot_recipe(log, tmp_path)
    rawFrameSample, masterFrame = _plot_inputs()
    masterFrame.data[1, 1] = np.nan
    summaries = []
    realSummary = mbiasModule.bias_distribution_summary

    def recording_summary(**kwargs: object) -> dict[str, object]:
        summary = realSummary(**kwargs)
        summaries.append(summary)
        return summary

    # ACT
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(mbiasModule, "bias_distribution_summary", recording_summary)
        recipe.plot_bias_distribution_qc(rawFrameSample, masterFrame, 2.5, 0.5, "/products/MASTER_BIAS_VIS.fits")

    # ASSERT
    assert len(summaries) == 1
    assert np.isfinite(summaries[0]["masterMean"])
    assert np.isfinite(summaries[0]["masterMedian"])
