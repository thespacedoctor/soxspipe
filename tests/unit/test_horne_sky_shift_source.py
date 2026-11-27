"""Contracts for the ``SKY SHIFT O<n> SRC`` cards: where each order's applied sky shift came from (DY-2461)."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.phase3 import fit_card_comment
from soxspipe.commonutils.sky_line_snr import combine_line_snrs
from soxspipe.commonutils.ug_wavelength_transfer import MIN_SKY_LINE_SNR
from tests.unit._horne_nir_scene import ORDERS as NIR_ORDERS
from tests.unit._horne_nir_scene import build_extractor as build_nir_extractor
from tests.unit._horne_nir_scene import build_orders as build_nir_orders
from tests.unit._horne_vis_scene import (
    CATALOGUE,
    I_LINE_OFFSET_NM,
    KEYWORDS,
    PIXEL_NM,
    R_LINE_OFFSET_NM,
    build_extractor,
    night_header,
    qc_by_name,
    vis_parameters,
)

pytestmark = pytest.mark.unit

SOURCE_CARD = re.compile(r"^SKY SHIFT O(\d+) SRC$")
SHIFT_CARD = re.compile(r"^SKY SHIFT O(\d+)$")
BRIGHT_HEIGHT = 20.0
FAINT_HEIGHT = 2.0
# THE SCENE SKY IS 1 e- BETWEEN LINES; WITH THE 3.8 e- VIS READ NOISE A LINE OF HEIGHT h HAS SNR h / SQRT(1 + 3.8^2),
# AND THE NEAREST 0.05 nm SAMPLE OF AN i LINE (0.04 nm OFF ITS CATALOGUE WAVELENGTH) IS 0.01 nm OFF ITS PEAK
PHOTON_NOISE = np.sqrt(1.0 + 3.8**2)
SAMPLED_PEAK_FRACTION = np.exp(-0.5 * (0.01 / 0.08) ** 2)
# AN i HEIGHT WHOSE COMBINED SNR (BELOW) IS 5.38: THE QC ROUNDS IT TO THE 5.4 THRESHOLD, BUT THE GUARD REFUSES IT
NEAR_THRESHOLD_HEIGHT = 5.38 * PHOTON_NOISE / (SAMPLED_PEAK_FRACTION * np.sqrt(len(CATALOGUE)))
REFUSED_SNR_RANGE = (5.35, MIN_SKY_LINE_SNR)


def _sky(wavelength: np.ndarray, centres: np.ndarray, height: float) -> np.ndarray:
    return 1.0 + sum(height * np.exp(-0.5 * ((wavelength - c) / 0.08) ** 2) for c in centres)


def _vis_scene(rHeight: float = BRIGHT_HEIGHT, iHeight: float = BRIGHT_HEIGHT) -> pd.DataFrame:
    wavelength = np.arange(500.0, 520.0, PIXEL_NM)
    frames = []
    for order, skyFlux in [
        (1, _sky(wavelength, CATALOGUE + R_LINE_OFFSET_NM, rHeight)),
        (2, np.ones_like(wavelength)),
        (3, np.ones_like(wavelength)),
        (4, _sky(wavelength, CATALOGUE + I_LINE_OFFSET_NM, iHeight)),
    ]:
        frames.append(
            pd.DataFrame({"order": order, "wavelengthMean": wavelength, "pixelScaleNm": PIXEL_NM, "skyFlux": skyFlux})
        )
    return pd.concat(frames, ignore_index=True)


def _without(detectorKey: str) -> dict:
    return {key: value for key, value in vis_parameters().items() if key != detectorKey}


def _sources(extractor: horne_extraction) -> dict[int, str]:
    names = extractor.qc["qc_name"]
    return {
        int(SOURCE_CARD.match(name).group(1)): value
        for name, value in zip(names, extractor.qc["qc_value"], strict=True)
        if SOURCE_CARD.match(name)
    }


def _shift_orders(extractor: horne_extraction) -> set[int]:
    return {int(SHIFT_CARD.match(name).group(1)) for name in extractor.qc["qc_name"] if SHIFT_CARD.match(name)}


def _tune_vis(extractor: horne_extraction, scene: pd.DataFrame) -> pd.DataFrame:
    return extractor.tune_wavelength_calibration_to_skylines(scene.copy(), "VIS")


def _tune_nir(extractor: horne_extraction, scene: pd.DataFrame) -> pd.DataFrame:
    return extractor._tune_or_skip(scene.copy(), "NIR")


def _i_snr(extractor: horne_extraction, scene: pd.DataFrame) -> float:
    orderDF = extractor._order_slice(scene, 4)
    return float(combine_line_snrs(extractor._order_sky_line_snrs(orderDF, 4)))


# VIS i


def test_vis_i_is_measured_when_its_guards_pass(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene())

    # ASSERT
    assert _sources(extractor)[4] == "measured"
    comment = qc_by_name(extractor, "SKY SHIFT O4 SRC")["qc_comment"]
    # THE SIGMA CLIP MAY DROP A LINE, SO ONLY THE PHRASE IS FIXED
    assert re.fullmatch(r"\d+ sky lines kept", comment)


def test_vis_i_has_no_shift_when_fewer_than_two_sky_lines_are_kept(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    extractor._get_local_skylines_for_order = lambda wmin, wmax, order, col: (
        (CATALOGUE[:1], CATALOGUE[:1]) if order == 4 else (CATALOGUE, CATALOGUE)
    )

    # ACT
    _tune_vis(extractor, _vis_scene())

    # ASSERT
    assert qc_by_name(extractor, "SKY SHIFT O4")["qc_value"] == 0.0
    source = qc_by_name(extractor, "SKY SHIFT O4 SRC")
    assert source["qc_value"] == "none"
    assert source["qc_comment"] == "1 sky lines kept (min 2)"


def test_vis_i_has_no_shift_when_the_sky_line_snr_is_below_the_threshold(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene(iHeight=FAINT_HEIGHT))

    # ASSERT
    assert qc_by_name(extractor, "SKY SHIFT O4")["qc_value"] == 0.0
    source = qc_by_name(extractor, "SKY SHIFT O4 SRC")
    assert source["qc_value"] == "none"
    assert source["qc_comment"] == "sky-line SNR below threshold"


def test_vis_i_refused_just_below_the_threshold_is_none_although_its_snr_card_reads_the_threshold(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _vis_scene(iHeight=NEAR_THRESHOLD_HEIGHT)
    extractor._get_local_skylines_for_order = lambda *args: (CATALOGUE, CATALOGUE)

    # ACT
    _tune_vis(extractor, scene)

    # ASSERT
    # THE UNROUNDED SNR IS REFUSED, BUT THE SNR CARD ROUNDS TO 1 DECIMAL AND SO READS THE THRESHOLD
    assert REFUSED_SNR_RANGE[0] <= _i_snr(extractor, scene) < REFUSED_SNR_RANGE[1]
    assert qc_by_name(extractor, "SKY LINE SNR O4")["qc_value"] == MIN_SKY_LINE_SNR
    assert qc_by_name(extractor, "SKY SHIFT O4")["qc_value"] == 0.0
    assert qc_by_name(extractor, "SKY SHIFT O4 SRC")["qc_value"] == "none"


# VIS r


def test_vis_r_measured_alone_is_measured(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene(iHeight=FAINT_HEIGHT))

    # ASSERT
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "measured"
    assert _sources(extractor)[1] == "measured"


def test_vis_r_predicted_from_i_is_model(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene(rHeight=FAINT_HEIGHT))

    # ASSERT
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "predicted"
    source = qc_by_name(extractor, "SKY SHIFT O1 SRC")
    assert source["qc_value"] == "model"
    assert source["qc_comment"] == "predicted from i and ABSROT"


def test_vis_r_combined_with_the_i_prediction_is_model(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene())

    # ASSERT
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "combined"
    source = qc_by_name(extractor, "SKY SHIFT O1 SRC")
    assert source["qc_value"] == "model"
    assert source["qc_comment"] == "combined with i prediction"


def test_vis_r_refused_is_none(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene(rHeight=FAINT_HEIGHT, iHeight=FAINT_HEIGHT))

    # ASSERT
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "none"
    assert qc_by_name(extractor, "SKY SHIFT O1")["qc_value"] == 0.0
    source = qc_by_name(extractor, "SKY SHIFT O1 SRC")
    assert source["qc_value"] == "none"
    assert source["qc_comment"] == "sky-line SNR below threshold"


def test_vis_r_without_r_from_i_coefficients_is_measured_although_no_r_source_card_is_written(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams=_without("r-from-i-transfer"))

    # ACT
    _tune_vis(extractor, _vis_scene())

    # ASSERT
    assert "R SHIFT SOURCE" not in set(extractor.qc["qc_name"])
    assert _sources(extractor)[1] == "measured"


def test_vis_r_without_r_from_i_coefficients_is_none_when_refused(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams=_without("r-from-i-transfer"))

    # ACT
    _tune_vis(extractor, _vis_scene(rHeight=FAINT_HEIGHT))

    # ASSERT
    assert "R SHIFT SOURCE" not in set(extractor.qc["qc_name"])
    assert _sources(extractor)[1] == "none"


# VIS u and g


def test_vis_u_and_g_from_the_transfer_model_are_model(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene())

    # ASSERT
    assert qc_by_name(extractor, "UG SHIFT SOURCE")["qc_value"] == "model"
    assert _sources(extractor)[2] == "model"
    assert _sources(extractor)[3] == "model"
    assert qc_by_name(extractor, "SKY SHIFT O2 SRC")["qc_comment"] == "predicted from r/i and ABSROT"


def test_vis_u_and_g_are_none_when_the_transfer_is_skipped(log: object, tmp_path: Path) -> None:
    # ARRANGE
    placeholder = {KEYWORDS["ABSROT"]: 0.0, KEYWORDS["PAR_ANG_START"]: 0.0}
    extractor = build_extractor(log, tmp_path, placeholder)

    # ACT
    _tune_vis(extractor, _vis_scene())

    # ASSERT
    assert qc_by_name(extractor, "UG SHIFT SOURCE")["qc_value"] == "none"
    assert _sources(extractor)[2] == "none"
    assert _sources(extractor)[3] == "none"
    assert qc_by_name(extractor, "SKY SHIFT O3 SRC")["qc_comment"] == "u/g transfer model not applied"


def test_vis_u_and_g_are_none_when_the_arm_has_no_transfer_coefficients(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams=_without("ug-wavelength-transfer"))

    # ACT
    _tune_vis(extractor, _vis_scene())

    # ASSERT
    assert "UG SHIFT SOURCE" not in set(extractor.qc["qc_name"])
    assert _sources(extractor)[2] == "none"
    assert _sources(extractor)[3] == "none"
    assert qc_by_name(extractor, "SKY SHIFT O2 SRC")["qc_comment"] == "no u/g transfer model"


# NIR


def test_nir_global_fit_orders_with_a_correction_are_measured(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_nir_extractor(log, tmp_path)

    # ACT
    _tune_nir(extractor, build_nir_orders())

    # ASSERT
    assert _sources(extractor) == dict.fromkeys(NIR_ORDERS, "measured")
    assert qc_by_name(extractor, f"SKY SHIFT O{NIR_ORDERS[0]} SRC")["qc_comment"] == "global fit, 10 sky lines"


def test_nir_global_fit_order_without_a_correction_is_none(log: object, tmp_path: Path) -> None:
    # ARRANGE
    unsigned = NIR_ORDERS[-1]
    extractor = build_nir_extractor(log, tmp_path)
    realDispersions = extractor._signed_order_dispersions
    extractor._signed_order_dispersions = lambda df, orders: {
        o: d for o, d in realDispersions(df, orders).items() if o != unsigned
    }

    # ACT
    _tune_nir(extractor, build_nir_orders())

    # ASSERT
    assert qc_by_name(extractor, f"SKY SHIFT O{unsigned}")["qc_value"] == 0.0
    assert _sources(extractor)[unsigned] == "none"
    assert qc_by_name(extractor, f"SKY SHIFT O{unsigned} SRC")["qc_comment"] == "no global-fit correction"
    assert _sources(extractor)[NIR_ORDERS[0]] == "measured"


def test_nir_global_fit_own_seed_order_is_measured(log: object, tmp_path: Path) -> None:
    # ARRANGE: LIKE THE REAL O10, WITH 3 SKY LINES AND AN ARC-SOLUTION OFFSET OF 8 px, THE ORDER KEEPS ITS OWN SHIFT
    offsetOrder = 12
    extractor = build_nir_extractor(log, tmp_path)
    scene = build_nir_orders(offsets={offsetOrder: 8.0}, lineCounts={offsetOrder: 3})

    # ACT
    _tune_nir(extractor, scene)

    # ASSERT
    assert extractor.skyCorrections[offsetOrder].source == "own-seed"
    assert _sources(extractor)[offsetOrder] == "measured"
    assert "own-seed" in qc_by_name(extractor, f"SKY SHIFT O{offsetOrder} SRC")["qc_comment"]


def test_nir_per_order_fallback_order_without_a_kept_sky_line_is_none_although_it_has_a_px_card(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE: ONE ORDER HOLDS LINES, SO THE GLOBAL FIT HAS TOO FEW AND EVERY ORDER IS TUNED ALONE
    extractor = build_nir_extractor(log, tmp_path)

    # ACT
    _tune_nir(extractor, build_nir_orders(linesOnly=NIR_ORDERS[:1]))

    # ASSERT
    assert qc_by_name(extractor, "SKY TUNING")["qc_value"] == "skylines-per-order"
    assert qc_by_name(extractor, f"SKY SHIFT O{NIR_ORDERS[1]} PX")["qc_value"] == 0.0
    sources = _sources(extractor)
    assert sources[NIR_ORDERS[0]] == "measured"
    assert {sources[o] for o in NIR_ORDERS[1:]} == {"none"}
    assert qc_by_name(extractor, f"SKY SHIFT O{NIR_ORDERS[1]} SRC")["qc_comment"] == "no sky line matched"


def test_the_generic_per_order_path_writes_source_cards_too(log: object, tmp_path: Path) -> None:
    # ARRANGE: X-SHOOTER NIR TUNES EVERY ORDER ALONE, WITH NO GLOBAL FIT
    extractor = build_nir_extractor(log, tmp_path, inst="XSH")
    lineless = NIR_ORDERS[-1]

    # ACT
    _tune_nir(extractor, build_nir_orders(linelessOrders=(lineless,)))

    # ASSERT
    sources = _sources(extractor)
    assert sources[lineless] == "none"
    assert {sources[o] for o in NIR_ORDERS[:-1]} == {"measured"}


# EVERY SHIFT CARD HAS A SOURCE CARD


def _run_vis(log: object, tmp_path: Path) -> horne_extraction:
    extractor = build_extractor(log, tmp_path, night_header())
    _tune_vis(extractor, _vis_scene(rHeight=FAINT_HEIGHT))
    return extractor


def _run_nir_global(log: object, tmp_path: Path) -> horne_extraction:
    extractor = build_nir_extractor(log, tmp_path)
    _tune_nir(extractor, build_nir_orders())
    return extractor


def _run_nir_fallback(log: object, tmp_path: Path) -> horne_extraction:
    extractor = build_nir_extractor(log, tmp_path)
    _tune_nir(extractor, build_nir_orders(linesOnly=NIR_ORDERS[:1]))
    return extractor


TUNING_RUNS = pytest.mark.parametrize(
    "tune", [_run_vis, _run_nir_global, _run_nir_fallback], ids=["vis", "nir-global", "nir-per-order"]
)


@TUNING_RUNS
def test_every_order_with_a_shift_card_has_a_source_card_and_the_reverse(tune: object, log: object, tmp_path: Path):
    # ARRANGE & ACT
    extractor = tune(log, tmp_path)

    # ASSERT
    assert _shift_orders(extractor)
    assert _shift_orders(extractor) == set(_sources(extractor))
    sourceRows = extractor.qc.loc[extractor.qc["qc_name"].str.match(SOURCE_CARD)]
    assert set(sourceRows["qc_value"]) <= {"measured", "model", "none"}
    assert sourceRows["to_header"].all()
    assert (
        sourceRows["qc_order"].astype(int).to_numpy()
        == [int(SOURCE_CARD.match(n).group(1)) for n in sourceRows["qc_name"]]
    ).all()
    assert sourceRows["qc_unit"].isna().all()


@TUNING_RUNS
def test_every_source_comment_fits_its_fits_card(tune: object, log: object, tmp_path: Path) -> None:
    # ARRANGE & ACT
    extractor = tune(log, tmp_path)

    # ASSERT
    rows = extractor.qc.loc[extractor.qc["qc_name"].str.match(SOURCE_CARD)]
    assert len(rows)
    for name, value, comment in zip(rows["qc_name"], rows["qc_value"], rows["qc_comment"], strict=True):
        assert fit_card_comment(f"ESO QC {name}", value, comment) == comment, name


def test_the_classification_is_kept_beside_the_shifts_for_the_error_cards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    _tune_vis(extractor, _vis_scene(rHeight=FAINT_HEIGHT))

    # ASSERT
    assert extractor.skyShiftSources == {1: "model", 2: "model", 3: "model", 4: "measured"}
    assert extractor.skyShiftSources == _sources(extractor)


def test_the_nir_global_fit_keeps_its_classification_beside_the_corrections(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_nir_extractor(log, tmp_path)

    # ACT
    _tune_nir(extractor, build_nir_orders())

    # ASSERT
    assert extractor.skyShiftSources == dict.fromkeys(NIR_ORDERS, "measured")


def test_tuning_all_orders_together_writes_no_shift_and_no_source_cards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    extractor.tune_wavelength_calibration_to_skylines(_vis_scene(), "VIS", byOrder=False)

    # ASSERT
    assert not [name for name in extractor.qc.get("qc_name", []) if name.startswith("SKY SHIFT")]
    assert extractor.skyShiftSources is None
