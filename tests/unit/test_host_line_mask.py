"""Finding host emission lines that fill the slit of a stare frame, so the sky fit can leave them out (DY-2242).

A host line that fills the slit is in the pixels the sky model treats as sky. The detector looks for narrow emission
that is present along most of the slit and that changes along it, which a sky line, a star's wings and a cosmic ray
do not.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.host_line_mask import _slit_cells, find_host_line_windows, host_line_pixel_mask
from soxspipe.commonutils.nod_host_contamination import host_feature

pytestmark = pytest.mark.unit

N_LINES = 1200
N_SLIT = 50
FIRST_NM = 600.0
STEP_NM = 0.04
SLIT_STEP_ARCSEC = 0.2
# THE SLIT IS TILTED AGAINST THE DETECTOR LINES, SO THE PIXELS OF EACH SLIT BIN SAMPLE A STEEP LINE AT DIFFERENT PHASES
TILT_LINES_PER_ARCSEC = 0.2
SKY_LEVEL = 30.0
# A HOST LINE 0.05 NM WIDE (SIGMA) THAT IS STRONGEST AT THE SLIT CENTRE, AND MOVES 0.04 NM PER ARCSEC ALONG THE SLIT
HOST_NM = 620.0
HOST_SIGMA_NM = 0.05
HOST_PEAK = 150.0
HOST_SLIT_SIGMA_ARCSEC = 3.0
HOST_ROTATION_NM_PER_ARCSEC = 0.04
# A HOST LINE THAT PASSES A LIMIT OF 5 SIGMA AND FAILS A LIMIT OF 8
FAINT_HOST_PEAK = 45.0
# A HOST LINE AT ABOUT 5.4 SIGMA, WHICH PASSES A LIMIT OF 5 AND FAILS 5.5
MARGINAL_HOST_PEAK = 40.0
SKY_LINE_NM = 640.0
SKY_LINE_PEAK = 800.0
SKY_LINE_SIGMA_NM = 0.05


def _order_frame(seed: int = 7) -> pd.DataFrame:
    """A flat-sky order table with its dispersion along y, as the sky fit receives it, before any feature is added."""
    rng = np.random.default_rng(seed)
    x, y = np.meshgrid(np.arange(N_SLIT), np.arange(N_LINES), indexing="ij")
    frame = pd.DataFrame(
        {
            "x": x.ravel(),
            "y": y.ravel(),
            "order": 1.0,
            "slit_position": (x.ravel() - N_SLIT / 2) * SLIT_STEP_ARCSEC,
        }
    )
    frame["wavelength"] = FIRST_NM + STEP_NM * (frame["y"] + TILT_LINES_PER_ARCSEC * frame["slit_position"])
    frame["error"] = np.sqrt(SKY_LEVEL)
    frame["flux"] = rng.normal(SKY_LEVEL, np.sqrt(SKY_LEVEL), len(frame))
    for flag in ("flagged_object_clipped", "flagged_bad_pixel_clipped", "flagged_edge_clipped"):
        frame[flag] = False
    return frame


def _add_host_line(frame: pd.DataFrame, centreNm: float = HOST_NM, peak: float = HOST_PEAK) -> pd.DataFrame:
    """Add a host line that moves along the slit and is stronger near the slit centre."""
    slit = frame["slit_position"]
    centre = centreNm + HOST_ROTATION_NM_PER_ARCSEC * slit
    lineProfile = np.exp(-0.5 * ((frame["wavelength"] - centre) / HOST_SIGMA_NM) ** 2)
    slitProfile = np.exp(-0.5 * (slit / HOST_SLIT_SIGMA_ARCSEC) ** 2)
    return frame.assign(flux=frame["flux"] + peak * lineProfile * slitProfile)


def _add_sky_line(
    frame: pd.DataFrame,
    slitGradient: float = 0.0,
    centreNm: float = SKY_LINE_NM,
    sigmaNm: float = SKY_LINE_SIGMA_NM,
) -> pd.DataFrame:
    """Add a sky line at one wavelength that changes along the slit by ``slitGradient`` (a fraction per arcsec)."""
    lineProfile = np.exp(-0.5 * ((frame["wavelength"] - centreNm) / sigmaNm) ** 2)
    slitScale = 1.0 + slitGradient * frame["slit_position"]
    return frame.assign(flux=frame["flux"] + SKY_LINE_PEAK * slitScale * lineProfile)


def _windows(frame: pd.DataFrame, **kwargs) -> tuple[host_feature, ...]:
    return find_host_line_windows(frame, **kwargs)


def test_finds_one_window_around_a_host_line_that_fills_the_slit() -> None:
    frame = _add_host_line(_order_frame())

    windows = _windows(frame)

    assert len(windows) == 1
    assert windows[0].wavelengthMinNm < HOST_NM < windows[0].wavelengthMaxNm
    assert windows[0].order == 1


def test_window_covers_the_line_at_every_slit_position() -> None:
    frame = _add_host_line(_order_frame())
    slitEnd = frame["slit_position"].abs().max()

    (window,) = _windows(frame)

    assert window.wavelengthMinNm <= HOST_NM - HOST_ROTATION_NM_PER_ARCSEC * slitEnd * 0.5
    assert window.wavelengthMaxNm >= HOST_NM + HOST_ROTATION_NM_PER_ARCSEC * slitEnd * 0.5


def test_finds_nothing_in_a_sky_with_no_lines() -> None:
    assert _windows(_order_frame()) == ()


def test_finds_nothing_for_a_sky_line_that_is_constant_along_the_slit() -> None:
    frame = _add_sky_line(_order_frame())

    assert _windows(frame) == ()


@pytest.mark.parametrize("phase", range(8))
def test_finds_nothing_for_a_bright_unresolved_sky_line_at_any_phase_of_the_wavelength_bins(phase: int) -> None:
    # THE TILT SAMPLES A LINE ONE DETECTOR LINE WIDE DIFFERENTLY IN EACH SLIT BIN, WHEREVER IT SITS IN THE BIN
    frame = _add_sky_line(_order_frame(), centreNm=SKY_LINE_NM + phase * STEP_NM / 4, sigmaNm=STEP_NM)

    assert _windows(frame) == ()


@pytest.mark.parametrize("phase", range(8))
def test_finds_nothing_for_a_sky_line_whose_core_the_clip_flagged_as_object(phase: int) -> None:
    # THE ROLLING CLIP OF THE SKY FIT FLAGS THE CORE OF A BRIGHT SKY LINE AS OBJECT AT EVERY SLIT POSITION
    centreNm = SKY_LINE_NM + phase * STEP_NM / 4
    frame = _add_sky_line(_order_frame(), centreNm=centreNm, sigmaNm=STEP_NM)
    frame["flagged_object_clipped"] = (frame["wavelength"] - centreNm).abs() < 1.5 * STEP_NM

    assert _windows(frame) == ()


def test_finds_nothing_for_a_sky_line_with_a_small_gradient_along_the_slit() -> None:
    frame = _add_sky_line(_order_frame(), slitGradient=0.02)

    assert _windows(frame) == ()


def test_finds_nothing_for_a_catalogue_sky_line_whose_wings_vary_along_the_slit() -> None:
    frame = _add_host_line(_order_frame(), centreNm=SKY_LINE_NM)

    windows = _windows(frame, skylineWavelengthsNm=[SKY_LINE_NM])

    assert windows == ()


def test_finds_a_host_line_away_from_the_catalogue_sky_lines() -> None:
    frame = _add_host_line(_order_frame())

    windows = _windows(frame, skylineWavelengthsNm=[SKY_LINE_NM, HOST_NM + 5.0])

    assert len(windows) == 1


def test_finds_nothing_for_a_point_source_with_wings_outside_the_object_mask() -> None:
    frame = _order_frame()
    slit = frame["slit_position"]
    # A SMOOTH STAR PROFILE WHOSE CORE IS FLAGGED AS OBJECT AND WHOSE WINGS ARE NOT
    wings = 200.0 * np.exp(-0.5 * (slit / 1.5) ** 2) * (1.0 + 0.3 * np.sin(frame["wavelength"] / 3.0))
    frame = frame.assign(flux=frame["flux"] + wings, flagged_object_clipped=slit.abs() < 1.0)

    assert _windows(frame) == ()


def test_finds_nothing_for_emission_only_in_the_pixels_flagged_as_object() -> None:
    frame = _add_host_line(_order_frame())
    frame["flagged_object_clipped"] = True

    assert _windows(frame) == ()


def test_finds_nothing_for_a_cosmic_ray_hit_in_one_slit_bin() -> None:
    frame = _order_frame()
    isHit = (frame["y"].between(400, 401)) & (frame["slit_position"].between(1.0, 1.9))
    frame.loc[isHit, "flux"] += 500.0

    assert _windows(frame) == ()


def test_finds_nothing_for_telluric_absorption_on_an_extended_continuum() -> None:
    frame = _order_frame()
    slit = frame["slit_position"]
    continuum = 100.0 * np.exp(-0.5 * (slit / 3.0) ** 2)
    absorption = 1.0 - 0.6 * np.exp(-0.5 * ((frame["wavelength"] - HOST_NM) / HOST_SIGMA_NM) ** 2)
    frame = frame.assign(flux=frame["flux"] + continuum * absorption)

    assert _windows(frame) == ()


def test_padding_widens_the_window_by_the_given_number_of_detector_lines() -> None:
    frame = _add_host_line(_order_frame())

    (bare,) = _windows(frame, paddingLines=0)
    (padded,) = _windows(frame, paddingLines=6)

    assert padded.wavelengthMinNm == pytest.approx(bare.wavelengthMinNm - 6 * STEP_NM)
    assert padded.wavelengthMaxNm == pytest.approx(bare.wavelengthMaxNm + 6 * STEP_NM)


def test_finds_the_window_when_the_dispersion_runs_along_x() -> None:
    frame = _add_host_line(_order_frame()).rename(columns={"x": "y", "y": "x"})

    windows = find_host_line_windows(frame)

    assert len(windows) == 1
    assert windows[0].wavelengthMinNm < HOST_NM < windows[0].wavelengthMaxNm


def test_a_higher_sigma_limit_drops_a_faint_host_line() -> None:
    frame = _add_host_line(_order_frame(), peak=FAINT_HOST_PEAK)

    assert len(_windows(frame, sigmaLimit=5.0)) == 1
    assert _windows(frame, sigmaLimit=8.0) == ()


def test_pixel_mask_marks_the_wavelengths_inside_each_window() -> None:
    windows = [host_feature(1, 10.0, 11.0, 9.0), host_feature(1, 20.0, 21.0, 9.0)]
    wavelength = np.array([9.9, 10.0, 10.5, 11.0, 11.1, 20.5, 25.0])

    isMasked = host_line_pixel_mask(wavelength, windows)

    assert isMasked.tolist() == [False, True, True, True, False, True, False]


def test_pixel_mask_is_empty_when_there_are_no_windows() -> None:
    assert not host_line_pixel_mask(np.array([1.0, 2.0]), []).any()


def test_finds_nothing_when_every_pixel_is_flagged_as_a_bad_pixel() -> None:
    frame = _add_host_line(_order_frame())
    frame["flagged_bad_pixel_clipped"] = True

    assert _windows(frame) == ()


def test_finds_nothing_in_an_order_that_is_one_detector_line_long() -> None:
    frame = _add_host_line(_order_frame())
    frame = frame[frame["y"] == 0]

    assert _windows(frame) == ()


def test_finds_nothing_in_an_order_whose_wavelength_does_not_change() -> None:
    frame = _add_host_line(_order_frame()).assign(wavelength=HOST_NM)

    assert _windows(frame) == ()


def test_finds_nothing_in_an_order_with_no_usable_pixel_outside_the_object_trace() -> None:
    frame = _add_host_line(_order_frame())
    frame["flagged_object_clipped"] = frame["slit_position"] < 4.0

    assert _windows(frame) == ()


def test_finds_nothing_and_warns_of_nothing_when_every_slit_position_is_not_a_number() -> None:
    frame = _add_host_line(_order_frame()).assign(slit_position=np.nan)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        windows = _windows(frame)

    assert windows == ()


def test_still_finds_the_line_when_some_slit_positions_are_not_a_number() -> None:
    frame = _add_host_line(_order_frame())
    frame.loc[frame.index[::97], "slit_position"] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        windows = _windows(frame)

    assert len(windows) == 1
    assert windows[0].wavelengthMinNm < HOST_NM < windows[0].wavelengthMaxNm


@pytest.mark.parametrize("tilt", [0.0, 0.2, 5.0])
def test_the_line_step_is_the_median_wavelength_increment_of_one_pixel_whatever_the_slit_tilt(tilt: float) -> None:
    # A STRONG TILT ADDS ITS SPREAD ACROSS THE SLIT TO THE WAVELENGTH RANGE OF THE ORDER, WHICH BIASED THE OLD ESTIMATE
    frame = _order_frame()
    frame["wavelength"] = FIRST_NM + STEP_NM * (frame["y"] + tilt * frame["slit_position"])

    cells = _slit_cells(frame)

    assert cells is not None
    assert cells[3] == pytest.approx(STEP_NM, rel=1e-6)


def test_the_line_step_follows_the_dispersion_axis_when_it_is_x() -> None:
    frame = _order_frame().rename(columns={"x": "y", "y": "x"})

    cells = _slit_cells(frame)

    assert cells is not None
    assert cells[3] == pytest.approx(STEP_NM, rel=1e-6)


def test_the_padding_is_in_true_detector_lines_for_a_strongly_tilted_order() -> None:
    frame = _order_frame()
    frame["wavelength"] = FIRST_NM + STEP_NM * (frame["y"] + 5.0 * frame["slit_position"])
    frame = _add_host_line(frame)

    (bare,) = _windows(frame, paddingLines=0)
    (padded,) = _windows(frame, paddingLines=6)

    assert padded.wavelengthMinNm == pytest.approx(bare.wavelengthMinNm - 6 * STEP_NM, rel=1e-6)


def test_the_default_sigma_limit_is_above_a_marginal_five_sigma_line() -> None:
    # A MARGINAL LINE AT ABOUT 5.4 SIGMA CAME AND WENT WITH THE PHASE OF THE WAVELENGTH BINS ON A NIR STANDARD, SO THE
    # DEFAULT LIMIT IS 5.5 AND A LINE BELOW IT IS NOT FOUND UNLESS THE LIMIT IS LOWERED
    frame = _add_host_line(_order_frame(), peak=MARGINAL_HOST_PEAK)

    assert _windows(frame) == ()
    assert len(_windows(frame, sigmaLimit=5.0)) == 1
