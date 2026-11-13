"""Flagging host emission that the nod subtraction has removed from the object spectrum (DY-2193).

An emission line that covers more of the slit than the nod throw is present at both nod positions, so the A-B
difference subtracts part of it from the object. The detector looks for narrow features at the negative trace that
the object's own spectrum does not explain.
"""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.nod_host_contamination import (
    find_host_contamination,
    host_contamination_message,
    host_feature,
    merge_host_features,
)
from soxspipe.commonutils.slit_alignment import build_two_d_map_images
from tests.factories.slit_alignment import synthetic_alignment_map, synthetic_nod_frame

pytestmark = pytest.mark.unit

# THE SLIT POSITIONS OF THE OBJECT AT THE A AND B NODS (ARCSEC): A 3 ARCSEC THROW
NOD_A = 1.5
NOD_B = -1.5
# 11 PIXELS OF 0.2 ARCSEC: A 2.2 ARCSEC APERTURE, SHORTER THAN THE THROW
APERTURE_PIXELS = 11
# THE HOST LINE: CENTRED 2 ARCSEC FROM THE OBJECT, WITH A 2 ARCSEC GAUSSIAN SIGMA ALONG THE SLIT, SO IT COVERS MORE
# OF THE SLIT THAN THE THROW
HOST_LINE_NM = 502.0
HOST_OFFSET_ARCSEC = 2.0
HOST_SLIT_SIGMA_ARCSEC = 2.0
HOST_LINE_SIGMA_NM = 0.04
HOST_LINE_PEAK = 150.0
# A SKY LINE THAT IS BRIGHTER IN THE A FRAME THAN IN THE B FRAME, SO IT LEAVES A RESIDUAL ALL ALONG THE SLIT
SKY_LINE_NM = 503.2
SKY_RESIDUAL_PEAK = 150.0
SKY_LINE_SIGMA_NM = 0.04


@pytest.fixture
def scene(tmp_path):
    return synthetic_alignment_map(tmp_path)


@pytest.fixture
def map_images(scene):
    return build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)


def _add_counts(frame, counts):
    """Return the frame with counts added to its data and their photon noise added to its error."""
    counts = np.nan_to_num(counts)
    frame.data = frame.data + counts
    frame.uncertainty.array = np.sqrt(frame.uncertainty.array**2 + np.clip(counts, 0.0, None))
    return frame


def _host_line(scene, objectSlit):
    """Return the counts of a narrow host line that is extended along the slit and moves with the object."""
    slitProfile = np.exp(-0.5 * ((scene.slit - objectSlit - HOST_OFFSET_ARCSEC) / HOST_SLIT_SIGMA_ARCSEC) ** 2)
    lineProfile = np.exp(-0.5 * ((scene.wavelength - HOST_LINE_NM) / HOST_LINE_SIGMA_NM) ** 2)
    return HOST_LINE_PEAK * slitProfile * lineProfile


def _sky_residual(scene, isFlat=True):
    """Return the counts of a sky line that is flat along the slit, or brighter towards both slit ends."""
    lineProfile = np.exp(-0.5 * ((scene.wavelength - SKY_LINE_NM) / SKY_LINE_SIGMA_NM) ** 2)
    slitProfile = 1.0 if isFlat else 0.2 + (scene.slit / 3.0) ** 2
    return SKY_RESIDUAL_PEAK * slitProfile * lineProfile


def _difference(scene, *, hasHost=False, hasSkyResidual=False, isSkyFlat=True):
    """Return the A-B difference of a nod pair of a continuum object, with an optional host line and sky residual."""
    frameA = synthetic_nod_frame(scene, objectSlit=NOD_A, mjd=60000.0, offset=1.5, hasSky=False, seed=1)
    frameB = synthetic_nod_frame(scene, objectSlit=NOD_B, mjd=60000.1, offset=-1.5, hasSky=False, seed=2)
    if hasHost:
        frameA = _add_counts(frameA, _host_line(scene, NOD_A))
        frameB = _add_counts(frameB, _host_line(scene, NOD_B))
    if hasSkyResidual:
        frameA = _add_counts(frameA, _sky_residual(scene, isFlat=isSkyFlat))
    difference = frameA.subtract(frameB)
    difference.header = frameA.header
    return difference


def test_a_host_line_longer_than_the_throw_is_flagged_at_its_wavelength(scene, map_images) -> None:
    # ARRANGE
    difference = _difference(scene, hasHost=True)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT: ONE FEATURE, IN THE ORDER THAT HOLDS THE LINE, AND ITS RANGE COVERS THE LINE
    assert len(features) == 1
    feature = features[0]
    assert feature.order == 1
    assert feature.wavelengthMinNm <= HOST_LINE_NM <= feature.wavelengthMaxNm
    assert feature.wavelengthMaxNm - feature.wavelengthMinNm < 0.5
    assert feature.significance > 5


def test_the_same_pair_without_the_host_line_has_no_feature(scene, map_images) -> None:
    # ARRANGE
    difference = _difference(scene)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT
    assert features == ()


def test_a_narrow_line_of_the_object_itself_is_explained_by_the_object(scene, map_images) -> None:
    # ARRANGE: A LINE THAT FOLLOWS THE OBJECT PROFILE ALONG THE SLIT APPEARS IN BOTH TRACES WITH THE SAME RATIO
    frameA = synthetic_nod_frame(scene, objectSlit=NOD_A, mjd=60000.0, offset=1.5, hasSky=False, seed=1)
    frameB = synthetic_nod_frame(scene, objectSlit=NOD_B, mjd=60000.1, offset=-1.5, hasSky=False, seed=2)
    lineProfile = np.exp(-0.5 * ((scene.wavelength - HOST_LINE_NM) / HOST_LINE_SIGMA_NM) ** 2)
    for frame, position in ((frameA, NOD_A), (frameB, NOD_B)):
        objectProfile = np.exp(-0.5 * ((scene.slit - position) / 0.3) ** 2)
        _add_counts(frame, 2000.0 * objectProfile * lineProfile)
    difference = frameA.subtract(frameB)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT
    assert features == ()


def test_a_sky_residual_flat_along_the_slit_is_removed_with_the_slit_background(scene, map_images) -> None:
    # ARRANGE: THE SKY CHANGED BETWEEN THE NODS. THE SKYLINE IS NOT IN THE CATALOGUE
    difference = _difference(scene, hasSkyResidual=True)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT
    assert features == ()


def test_a_sky_residual_that_varies_along_the_slit_does_not_count_at_a_catalogue_skyline(scene, map_images) -> None:
    # ARRANGE
    difference = _difference(scene, hasSkyResidual=True, isSkyFlat=False)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS, skylineWavelengthsNm=(SKY_LINE_NM,))

    # ASSERT
    assert features == ()


def test_the_same_varying_sky_residual_counts_when_the_skyline_is_not_in_the_catalogue(scene, map_images) -> None:
    # ARRANGE: THE SLIT BACKGROUND DOES NOT REMOVE IT, SO THE MASK IS WHAT REMOVES IT ABOVE
    difference = _difference(scene, hasSkyResidual=True, isSkyFlat=False)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT
    assert [feature.wavelengthMinNm <= SKY_LINE_NM <= feature.wavelengthMaxNm for feature in features] == [True]


def test_a_higher_sigma_limit_drops_a_weaker_feature(scene, map_images) -> None:
    # ARRANGE
    difference = _difference(scene, hasHost=True)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS, sigmaLimit=1000.0)

    # ASSERT
    assert features == ()


def test_a_masked_pixel_is_left_out_of_the_search(scene, map_images) -> None:
    # ARRANGE: A HOT PIXEL ON THE NEGATIVE TRACE WOULD BE A FEATURE IF IT WERE READ
    difference = _difference(scene)
    rows, columns = np.nonzero((scene.orderMap == 2) & (np.abs(scene.slit - NOD_B) < 0.3))
    hot = (rows >= 300) & (rows <= 302)
    difference.data[rows[hot], columns[hot]] += 5000.0
    difference.mask = np.zeros(difference.shape, dtype=bool)
    difference.mask[rows[hot], columns[hot]] = True

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT
    assert features == ()


def test_an_order_whose_traces_are_closer_than_the_aperture_is_skipped(scene, map_images) -> None:
    # ARRANGE: AN APERTURE LONGER THAN THE 3 ARCSEC THROW
    difference = _difference(scene, hasHost=True)

    # ACT
    features = find_host_contamination(difference, map_images, aperturePixels=20)

    # ASSERT
    assert features == ()


def test_overlapping_features_from_two_orders_or_cycles_merge_into_one() -> None:
    # ARRANGE
    features = [
        host_feature(order=3, wavelengthMinNm=667.10, wavelengthMaxNm=667.25, significance=9.0),
        host_feature(order=4, wavelengthMinNm=667.20, wavelengthMaxNm=667.30, significance=12.0),
        host_feature(order=2, wavelengthMinNm=656.20, wavelengthMaxNm=656.30, significance=7.0),
    ]

    # ACT
    merged = merge_host_features(features)

    # ASSERT
    assert merged == (
        host_feature(order=None, wavelengthMinNm=656.20, wavelengthMaxNm=656.30, significance=7.0),
        host_feature(order=None, wavelengthMinNm=667.10, wavelengthMaxNm=667.30, significance=12.0),
    )


def test_the_message_names_each_wavelength_range_and_is_none_without_features() -> None:
    # ARRANGE
    merged = (
        host_feature(order=None, wavelengthMinNm=656.2, wavelengthMaxNm=656.3, significance=7.0),
        host_feature(order=None, wavelengthMinNm=667.1, wavelengthMaxNm=667.3, significance=12.0),
    )

    # ACT
    message = host_contamination_message(merged)

    # ASSERT
    assert "656.20-656.30 nm, 667.10-667.30 nm" in message
    assert host_contamination_message(()) is None


def test_a_spike_in_the_positive_trace_alone_is_not_at_the_opposite_nod_position(scene, map_images) -> None:
    # ARRANGE: A COSMIC-RAY-LIKE SPIKE ON THE OBJECT IN THE A FRAME. THE NEGATIVE TRACE IS UNCHANGED THERE
    difference = _difference(scene)
    spike = np.exp(-0.5 * ((scene.slit - NOD_A) / 0.3) ** 2) * np.exp(
        -0.5 * ((scene.wavelength - HOST_LINE_NM) / HOST_LINE_SIGMA_NM) ** 2
    )
    difference.data = difference.data + 2000.0 * np.nan_to_num(spike)

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT
    assert features == ()


def test_a_compact_unmasked_hit_on_the_negative_trace_is_not_extended_emission(scene, map_images) -> None:
    # ARRANGE: A COSMIC RAY ONE PIXEL LONG ALONG THE SLIT, OVER THREE DETECTOR LINES, LEFT UNMASKED IN THE B FRAME
    difference = _difference(scene)
    rows, columns = np.nonzero((scene.orderMap == 2) & (np.abs(scene.slit - NOD_B) < 0.1))
    hit = (rows >= 300) & (rows <= 302)
    difference.data[rows[hit], columns[hit]] -= 5000.0

    # ACT
    features = find_host_contamination(difference, map_images, APERTURE_PIXELS)

    # ASSERT
    assert features == ()


def test_each_skipped_order_is_logged_with_its_reason(scene, map_images, log) -> None:
    # ARRANGE: AN APERTURE LONGER THAN THE THROW, SO NO ORDER CAN BE CHECKED
    difference = _difference(scene, hasHost=True)

    # ACT
    find_host_contamination(difference, map_images, aperturePixels=20, log=log)

    # ASSERT
    messages = [message for level, message in log.messages if level == "debug"]
    assert [message for message in messages if "not checked" in message and "closer than" in message] == [
        f"Host contamination: order {order} not checked: the traces are missing or closer than the aperture"
        for order in scene.orders
    ]
