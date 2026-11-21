"""The ON/OFF classification of offset frames, per template run (DY-2354)."""

from __future__ import annotations

from typing import Any

import pytest

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.offset_frames import OffsetClassificationError, classify_offset_frames
from tests.factories.offset_frames import fixed_sky_run, generic_offset_run

pytestmark = pytest.mark.unit


class FakeFrame:
    """A frame that carries only a header."""

    def __init__(self, header: dict[str, object]) -> None:
        self.header = header


def _classify(log: Any, headers: list[dict[str, object]], names: list[str] | None = None) -> tuple[Any, ...]:
    """Classify frames built from ``headers`` with the real SOXS keyword aliases."""
    frames = [FakeFrame(header) for header in headers]
    fileNames = names or [f"frame_{index}.fits" for index in range(len(headers))]
    kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get
    return frames, classify_offset_frames(log, frames, fileNames, kw)


# THE CD-44 12736 CYCLE OF 2025-05-14: OBJECT, SKY, SKY, OBJECT. THE FIRST STEP IS ONLY JITTER
OSSO_STEPS = [(0.0, -0.8), (30.0, 30.8), (0.0, 0.1), (-30.0, -30.1)]


def test_a_fixed_sky_cycle_puts_the_object_positions_in_on_and_the_sky_positions_in_off(log: Any) -> None:
    # ARRANGE
    headers = fixed_sky_run(OSSO_STEPS, jitter=2.0)
    names = ["obj_1", "sky_1", "sky_2", "obj_2"]

    # ACT
    frames, (frameOn, frameOff, namesOn, namesOff) = _classify(log, headers, names)

    # ASSERT
    assert frameOn == [frames[0], frames[3]]
    assert frameOff == [frames[1], frames[2]]
    assert namesOn == ["obj_1", "obj_2"]
    assert namesOff == ["sky_1", "sky_2"]


def test_a_sky_position_whose_offsets_sum_to_a_negative_number_is_still_off(log: Any) -> None:
    # ARRANGE
    headers = fixed_sky_run([(0.0, -0.4), (-20.0, 10.4), (0.0, 0.1), (20.0, -10.5)])

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, headers)

    # ASSERT
    assert frameOn == [frames[0], frames[3]]
    assert frameOff == [frames[1], frames[2]]


@pytest.mark.parametrize("jitterStep", [(0.0, -0.4), (0.0, 0.4), (-0.3, 0.2), (0.3, 0.0)])
def test_a_jitter_only_step_at_the_object_position_is_on_whatever_its_sign(
    log: Any, jitterStep: tuple[float, float]
) -> None:
    # ARRANGE
    headers = fixed_sky_run([jitterStep, (20.0, 20.0), (0.0, 0.1), (-20.0, -20.1)], jitter=1.0)

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, headers)

    # ASSERT
    assert frameOn == [frames[0], frames[3]]
    assert frameOff == [frames[1], frames[2]]


def test_the_first_exposure_of_the_run_is_the_one_with_the_lowest_exposure_number_whatever_the_input_order(
    log: Any,
) -> None:
    # ARRANGE
    headers = fixed_sky_run(OSSO_STEPS, firstExpno=2)
    shuffled = [headers[3], headers[1], headers[0], headers[2]]

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, shuffled)

    # ASSERT
    assert frameOn == [frames[2], frames[0]]
    assert frameOff == [frames[1], frames[3]]


def _steps_to_reach(positions: list[float]) -> list[tuple[float, float]]:
    """Return the ``SEQ FIXOFF`` steps (due north) that take the telescope through ``positions`` in turn."""
    previous = 0.0
    steps = []
    for position in positions:
        steps.append((0.0, position - previous))
        previous = position
    return steps


@pytest.mark.parametrize(
    ("jitter", "distance", "positions", "onIndices"),
    [
        # 2.9 ARCSEC IS WITHIN THE JITTER WIDTH (2.0) PLUS 1 ARCSEC
        (2.0, 2.9, [0.0, 2.9, 20.0, 20.0, 20.0, 0.0], [0, 1, 5]),
        # 2.1 ARCSEC IS BEYOND THE JITTER WIDTH (1.0) PLUS 1 ARCSEC
        (1.0, 2.1, [0.0, 2.1, 20.0, 20.0, 0.0, 0.0], [0, 4, 5]),
    ],
)
def test_a_frame_is_on_within_the_jitter_width_plus_one_arcsec_of_the_first_exposure(
    log: Any, jitter: float, distance: float, positions: list[float], onIndices: list[int]
) -> None:
    # ARRANGE
    headers = fixed_sky_run(_steps_to_reach(positions), jitter=jitter)

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, headers)

    # ASSERT
    assert distance == positions[1]
    assert frameOn == [frames[index] for index in onIndices]
    assert frameOff == [frames[index] for index in range(len(frames)) if index not in onIndices]


def test_a_frame_whose_fixoff_sum_disagrees_with_its_pointing_fails_the_sof_and_names_the_frame(log: Any) -> None:
    # ARRANGE
    headers = fixed_sky_run(OSSO_STEPS)
    headers[1]["HIERARCH ESO SEQ FIXOFF DEC"] += 5.0
    names = ["obj_1.fits", "sky_1.fits", "sky_2.fits", "obj_2.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match="sky_1.fits"):
        _classify(log, headers, names)

    assert any(level == "error" and "sky_1.fits" in message for level, message in log.messages)


@pytest.mark.parametrize(
    "missingCard",
    [
        "ESO TPL START",
        "ESO TPL EXPNO",
        "ESO SEQ JITTER WIDTH",
        "HIERARCH ESO SEQ FIXOFF RA",
        "HIERARCH ESO SEQ FIXOFF DEC",
        "RA",
        "DEC",
    ],
)
def test_a_frame_missing_a_card_the_rule_reads_fails_the_sof_and_names_the_frame_and_the_card(
    log: Any, missingCard: str
) -> None:
    # ARRANGE
    headers = fixed_sky_run(OSSO_STEPS)
    del headers[2][missingCard]
    names = ["obj_1.fits", "sky_1.fits", "sky_2.fits", "obj_2.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"sky_2\.fits.*" + missingCard):
        _classify(log, headers, names)


def test_each_template_run_gets_its_own_first_exposure_as_the_reference(log: Any) -> None:
    # ARRANGE
    first = fixed_sky_run([(0.0, 0.0), (20.0, 20.0)], start="2025-05-14T04:00:00", acquisition=(-9.6, -3.9))
    second = fixed_sky_run(
        [(0.0, 0.0), (20.0, 20.0)], start="2025-05-14T05:00:00", acquisition=(50.0, 40.0), firstMjd=60001.0
    )

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, first + second)

    # ASSERT
    assert frameOn == [frames[0], frames[2]]
    assert frameOff == [frames[1], frames[3]]


def test_on_and_off_frames_are_paired_within_each_run_even_when_the_runs_overlap_in_time(log: Any) -> None:
    # ARRANGE
    runA = fixed_sky_run([(0.0, 0.0), (20.0, 20.0)], start="2025-05-14T04:00:00")
    runB = fixed_sky_run([(0.0, 0.0), (20.0, 20.0)], start="2025-05-14T05:00:00", acquisition=(50.0, 40.0))
    for header, mjd in zip(runA + runB, [1.0, 10.0, 5.0, 6.0], strict=True):
        header["MJD-OBS"] = mjd

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, runB + runA)

    # ASSERT
    pairedRuns = [
        (on.header["ESO TPL START"], off.header["ESO TPL START"]) for on, off in zip(frameOn, frameOff, strict=True)
    ]
    assert pairedRuns == [("2025-05-14T04:00:00",) * 2, ("2025-05-14T05:00:00",) * 2]


def test_a_run_with_more_on_than_off_frames_fails_the_sof_even_when_the_sof_totals_balance(log: Any) -> None:
    # ARRANGE
    threeOnOneOff = fixed_sky_run([(0.0, 0.0), (0.0, 0.1), (20.0, 20.0), (-20.0, -20.0)], start="2025-05-14T04:00:00")
    oneOnThreeOff = fixed_sky_run(
        [(0.0, 0.0), (20.0, 20.0), (0.0, 0.1), (0.0, -0.1)], start="2025-05-14T05:00:00", acquisition=(50.0, 40.0)
    )

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match="3 ON frames and 1 OFF frames.*2025-05-14T04:00:00"):
        _classify(log, threeOnOneOff + oneOnThreeOff)


def test_a_run_with_no_off_frame_fails_the_sof_and_suggests_the_stare_recipe(log: Any) -> None:
    # ARRANGE
    headers = fixed_sky_run([(0.0, 0.0), (0.0, 0.1)])

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"2025-05-14T04:00:00.*no OFF frame.*soxs-stare"):
        _classify(log, headers)


def test_a_balanced_sof_announces_its_cycles_and_runs(log: Any) -> None:
    # ARRANGE
    headers = fixed_sky_run(OSSO_STEPS)

    # ACT
    _classify(log, headers)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# PROCESSING 2 ON-OFF OFFSET CYCLES IN 1 TEMPLATE RUN"]


def test_a_generic_run_takes_on_and_off_from_dpr_type_and_needs_no_fixoff_cards(log: Any) -> None:
    # ARRANGE
    headers = generic_offset_run([(5.4, -7.9), (3.4, 12.1)], ["OBJECT", "SKY"])
    names = ["object.fits", "sky.fits"]

    # ACT
    frames, (frameOn, frameOff, namesOn, namesOff) = _classify(log, headers, names)

    # ASSERT
    assert frameOn == [frames[0]]
    assert frameOff == [frames[1]]
    assert (namesOn, namesOff) == (["object.fits"], ["sky.fits"])


def test_a_generic_sky_frame_is_off_wherever_the_telescope_pointed(log: Any) -> None:
    # ARRANGE
    headers = generic_offset_run([(0.0, 0.0), (0.0, 0.5), (0.4, 0.0), (3.0, 4.0)], ["OBJECT", "SKY", "OBJECT", "SKY"])

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, headers)

    # ASSERT
    assert frameOn == [frames[0], frames[2]]
    assert frameOff == [frames[1], frames[3]]


def test_a_generic_object_frame_beyond_one_arcsec_of_the_first_object_frame_fails_the_sof_and_names_it(
    log: Any,
) -> None:
    # ARRANGE
    headers = generic_offset_run(
        [(2.8, -16.3), (4.8, -14.3), (8.8, -10.3), (14.8, -4.2)], ["OBJECT", "SKY", "OBJECT", "SKY"]
    )
    names = ["pos1.fits", "pos2.fits", "pos3.fits", "pos4.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"pos3\.fits.*pos1\.fits") as error:
        _classify(log, headers, names)

    assert "pos2.fits" not in str(error.value)
    assert "pos4.fits" not in str(error.value)


def test_a_generic_run_with_no_sky_frame_fails_the_sof_and_suggests_the_stare_recipe(log: Any) -> None:
    # ARRANGE
    headers = generic_offset_run([(-4.1, -4.8)] * 4, ["OBJECT"] * 4)

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"2025-12-01T02:43:49.*no OFF frame.*soxs-stare"):
        _classify(log, headers)


def test_a_generic_run_of_object_frames_that_moved_off_the_star_fails_for_each_frame_and_for_the_missing_sky(
    log: Any,
) -> None:
    # ARRANGE
    headers = generic_offset_run([(-31.1, 28.1), (-29.1, 30.1), (-25.1, 34.1)], ["OBJECT"] * 3)
    names = ["pos1.fits", "pos2.fits", "pos3.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError) as error:
        _classify(log, headers, names)

    message = str(error.value)
    assert "pos2.fits" in message
    assert "pos3.fits" in message
    assert "no OFF frame" in message


def test_the_reference_of_a_generic_run_is_its_first_object_frame_not_its_first_exposure(log: Any) -> None:
    # ARRANGE
    headers = generic_offset_run(
        [(20.0, 20.0), (0.0, 0.0), (0.3, 0.0), (20.0, 20.0)], ["SKY", "OBJECT", "OBJECT", "SKY"]
    )

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, headers)

    # ASSERT
    assert frameOn == [frames[1], frames[2]]
    assert frameOff == [frames[0], frames[3]]


def test_a_generic_frame_that_is_neither_object_nor_sky_fails_the_sof_and_names_the_frame(log: Any) -> None:
    # ARRANGE
    headers = generic_offset_run([(0.0, 0.0), (0.0, 0.0), (20.0, 20.0)], ["OBJECT", "LAMP,FLAT", "SKY"])
    names = ["object.fits", "flat.fits", "sky.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"flat\.fits.*LAMP,FLAT"):
        _classify(log, headers, names)


def test_a_generic_standard_star_run_is_classified_like_an_object_run(log: Any) -> None:
    # ARRANGE
    headers = generic_offset_run([(0.0, 0.0), (20.0, 20.0)], ["STD,FLUX", "SKY"])

    # ACT
    frames, (frameOn, frameOff, _, _) = _classify(log, headers)

    # ASSERT
    assert (frameOn, frameOff) == ([frames[0]], [frames[1]])


def test_a_run_that_mixes_the_two_templates_fails_the_sof_and_names_the_frames(log: Any) -> None:
    # ARRANGE
    fixed = fixed_sky_run([(0.0, 0.0), (20.0, 20.0)], start="2025-12-01T02:43:49")
    generic = generic_offset_run(
        [(0.0, 0.0), (20.0, 20.0)], ["OBJECT", "SKY"], start="2025-12-01T02:43:49", firstExpno=3
    )
    names = ["fixed_1.fits", "fixed_2.fits", "generic_1.fits", "generic_2.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"generic_1\.fits.*generic_2\.fits.*fixed_1\.fits"):
        _classify(log, fixed + generic, names)


def test_a_frame_with_an_unreadable_card_fails_the_sof_and_names_the_frame(log: Any) -> None:
    # ARRANGE
    headers = fixed_sky_run(OSSO_STEPS)
    headers[1]["ESO TPL EXPNO"] = "second"
    names = ["obj_1.fits", "sky_1.fits", "sky_2.fits", "obj_2.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"sky_1\.fits.*unreadable.*ESO TPL EXPNO"):
        _classify(log, headers, names)


def test_a_frame_taken_with_a_technique_that_is_not_an_offset_template_fails_the_sof_and_names_the_frame(
    log: Any,
) -> None:
    # ARRANGE
    headers = fixed_sky_run(OSSO_STEPS)
    headers[3]["ESO DPR TECH"] = "ECHELLE,SLIT,STARE"
    names = ["obj_1.fits", "sky_1.fits", "sky_2.fits", "stare.fits"]

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=r"stare\.fits.*ECHELLE,SLIT,STARE"):
        _classify(log, headers, names)
