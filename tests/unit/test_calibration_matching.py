"""Unit tests for the nearest-k calibration matcher in `soxspipe.commonutils.calibration_matching`.

Each test hand-builds small upstream and downstream product tables and checks which upstream product rows
`match_calibrations` pairs with each downstream SOF, for one row of the rule table in `CALIBRATION_RULES`.
"""

from __future__ import annotations

import pandas as pd
import pytest

from soxspipe.commonutils.calibration_matching import CALIBRATION_RULES, match_calibrations, sql_like

pytestmark = pytest.mark.unit

# RECIPE ORDER AS soxs_sof_map.yaml ENUMERATES IT
RECIPE_ORDER = {
    "mbias": 0,
    "mdark": 1,
    "disp_solution": 2,
    "order_centres": 3,
    "mflat": 4,
    "spat_solution": 5,
    "stare_std": 6,
    "nod_std": 7,
    "offset_std": 8,
    "stare_obj": 9,
    "nod_obj": 10,
    "offset_obj": 11,
}


def product(file, recipe, mjd, **overrides):
    """Return one `product_frames`-shaped row with VIS defaults."""
    row = {
        "file": file,
        "sof": f"{file}.sof",
        "eso pro catg": f"{recipe.upper()}_VIS",
        "filepath": f"./product/{file}",
        "recipe": recipe,
        "recipe_order": RECIPE_ORDER[recipe],
        "status": "pass",
        "complete": 1,
        "eso seq arm": "VIS",
        "mjd-obs": mjd,
        "exptime": 0.0,
        "binning": "1x1",
        "gain": 1.2,
        "rospeed": "1",
        "slit": "SLIT1.0",
    }
    row.update(overrides)
    return row


def downstream_sof(sof, recipe, mjd, **overrides):
    """Return one incomplete downstream SOF row (as `z_unique_product_sof_files` gives it)."""
    row = product(f"{sof}-product", recipe, mjd, **{"status": None, "complete": 0, **overrides})
    row["sof"] = sof
    return row


def matched_files(matches):
    """Return `{sof: sorted list of matched upstream files}`."""
    return {sof: sorted(group["file"]) for sof, group in matches.groupby("sof")}


def test_mbias_matches_the_nearest_bias_in_mjd():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("bias-early", "mbias", 100.0),
            product("bias-near", "mbias", 104.9),
            product("bias-late", "mbias", 110.0),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["bias-near"]}


def test_mbias_requires_matching_rospeed_binning_and_gain_on_vis():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("bias-wrong-rospeed", "mbias", 105.0, rospeed="2"),
            product("bias-wrong-binning", "mbias", 105.0, binning="2x2"),
            product("bias-wrong-gain", "mbias", 105.0, gain=2.4),
            product("bias-far-but-right", "mbias", 90.0),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["bias-far-but-right"]}


def test_nir_waives_the_binning_and_gain_match():
    # ARRANGE
    upstream = pd.DataFrame(
        [product("nir-bias", "mbias", 105.0, **{"eso seq arm": "NIR", "binning": None, "gain": None})]
    )
    downstream = pd.DataFrame(
        [downstream_sof("dsol.sof", "disp_solution", 105.0, **{"eso seq arm": "NIR", "binning": "1x1", "gain": None})]
    )

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["nir-bias"]}


def test_a_null_hard_match_value_never_matches():
    # ARRANGE
    upstream = pd.DataFrame([product("bias-null-rospeed", "mbias", 105.0, rospeed=None)])
    downstream = pd.DataFrame(
        [
            downstream_sof("null.sof", "disp_solution", 105.0, rospeed=None),
            downstream_sof("set.sof", "disp_solution", 105.0),
        ]
    )

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matches.empty


def test_upstream_must_be_complete_passing_same_arm_and_earlier_in_the_cascade():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("bias-incomplete", "mbias", 105.0, complete=0),
            product("bias-failed", "mbias", 105.0, status="fail"),
            product("bias-other-arm", "mbias", 105.0, **{"eso seq arm": "UVB"}),
            product("bias-null-status", "mbias", 80.0, status=None),
        ]
    )
    downstream = pd.DataFrame(
        [
            downstream_sof("dsol.sof", "disp_solution", 105.0),
            downstream_sof("done.sof", "disp_solution", 105.0, complete=1),
            downstream_sof("bias.sof", "mbias", 105.0),
        ]
    )

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["bias-null-status"]}


@pytest.mark.parametrize(
    ("target", "expected"),
    [(50.0, "bias-a"), (101.4, "bias-a"), (101.6, "bias-b"), (102.9, "bias-c"), (500.0, "bias-c")],
)
def test_nearest_is_found_on_either_side_and_beyond_both_array_edges(target, expected):
    # ARRANGE
    upstream = pd.DataFrame(
        [product("bias-c", "mbias", 103.0), product("bias-a", "mbias", 101.0), product("bias-b", "mbias", 102.0)]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", target)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": [expected]}


def test_ties_break_by_earlier_mjd_then_filepath():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("bias-after", "mbias", 106.0),
            product("bias-before-z", "mbias", 104.0, filepath="./product/z"),
            product("bias-before-a", "mbias", 104.0, filepath="./product/a"),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["bias-before-a"]}


def test_tie_break_considers_tied_rows_outside_the_k_window():
    # ARRANGE -- FOUR PRODUCTS SHARE ONE MJD; THE SMALLEST FILEPATH SORTS FURTHEST FROM THE INSERTION POINT
    upstream = pd.DataFrame(
        [product(f"bias-{name}", "mbias", 104.0, filepath=f"./product/{name}") for name in ("a", "b", "c", "d")]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["bias-a"]}


def test_empty_upstream_or_downstream_gives_no_matches():
    # ARRANGE
    upstream = pd.DataFrame([product("bias", "mbias", 105.0)])
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0)])

    # ACT
    noUpstream = match_calibrations(upstream.iloc[0:0], downstream, CALIBRATION_RULES["mbias"])
    noDownstream = match_calibrations(upstream, downstream.iloc[0:0], CALIBRATION_RULES["mbias"])

    # ASSERT
    assert noUpstream.empty
    assert noDownstream.empty
    assert list(noUpstream.columns) == list(noDownstream.columns)


def test_a_match_row_carries_the_view_columns():
    # ARRANGE
    upstream = pd.DataFrame([product("bias", "mbias", 104.5)])
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0, complete=-1)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matches.to_dict("records") == [
        {
            "sof": "dsol.sof",
            "file": "bias",
            "tag": "MBIAS_VIS",
            "filepath": "./product/bias",
            "complete": -1,
            "upstream_status": "pass",
            "upstream_recipe": "mbias",
            "upstream_recipe_order": 0,
            "downstream_recipe": "disp_solution",
            "downstream_recipe_order": 2,
            "mjd_diff": 0.5,
            "rank": 1,
        }
    ]


def nir(**overrides):
    """Return NIR overrides: NIR products carry no binning or gain."""
    return {"eso seq arm": "NIR", "binning": None, "gain": None, "rospeed": "3", **overrides}


def test_mdark_ranks_by_exptime_before_mjd_and_requires_rospeed():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("dark-near-wrong-exptime", "mdark", 105.0, **nir(exptime=10.0)),
            product("dark-far-right-exptime", "mdark", 90.0, **nir(exptime=300.0)),
            product("dark-farther-right-exptime", "mdark", 80.0, **nir(exptime=300.0)),
            product("dark-wrong-rospeed", "mdark", 105.0, **nir(exptime=300.0, rospeed="1")),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0, **nir(exptime=300.0))])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mdark"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["dark-far-right-exptime"]}


def test_mdark_keeps_both_exptimes_tied_in_distance_then_picks_by_mjd():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("dark-shorter", "mdark", 90.0, **nir(exptime=200.0)),
            product("dark-longer", "mdark", 104.0, **nir(exptime=400.0)),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0, **nir(exptime=300.0))])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mdark"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["dark-longer"]}


def test_mdark_uses_only_nir_darks():
    # ARRANGE
    upstream = pd.DataFrame([product("vis-dark", "mdark", 105.0, exptime=300.0)])
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0, exptime=300.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mdark"])

    # ASSERT
    assert matches.empty


def test_disp_solution_serves_only_order_centres_and_spat_solution():
    # ARRANGE
    upstream = pd.DataFrame([product("dsol", "disp_solution", 105.0, rospeed="2", slit="SLIT5.0")])
    downstream = pd.DataFrame(
        [
            downstream_sof("ocen.sof", "order_centres", 105.0),
            downstream_sof("ssol.sof", "spat_solution", 105.0),
            downstream_sof("mflat.sof", "mflat", 105.0),
        ]
    )

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["disp_solution"])

    # ASSERT
    assert matched_files(matches) == {"ocen.sof": ["dsol"], "ssol.sof": ["dsol"]}


def test_order_centres_serves_only_mflat():
    # ARRANGE
    upstream = pd.DataFrame([product("ocen", "order_centres", 105.0)])
    downstream = pd.DataFrame(
        [downstream_sof("mflat.sof", "mflat", 105.0), downstream_sof("ssol.sof", "spat_solution", 105.0)]
    )

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["order_centres"])

    # ASSERT
    assert matched_files(matches) == {"mflat.sof": ["ocen"]}


def test_mflat_keeps_the_two_nearest_rows_with_matching_slit_and_rospeed():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("flat-near-mflat", "mflat", 104.0, **{"eso pro catg": "MFLAT_VIS"}),
            product("flat-near-order", "mflat", 104.0, **{"eso pro catg": "ORDER_TAB_VIS"}),
            product("flat-far", "mflat", 90.0),
            product("flat-wrong-slit", "mflat", 105.0, slit="SLIT5.0"),
            product("flat-wrong-rospeed", "mflat", 105.0, rospeed="2"),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("stare.sof", "stare_obj", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mflat"])

    # ASSERT
    assert matched_files(matches) == {"stare.sof": ["flat-near-mflat", "flat-near-order"]}
    assert sorted(matches["rank"]) == [1, 2]


def test_mflat_for_spat_solution_waives_slit_and_rospeed_and_keeps_only_order_tables():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("flat-mflat", "mflat", 104.0, slit="SLIT5.0", rospeed="2", **{"eso pro catg": "MFLAT_VIS"}),
            product("flat-order", "mflat", 104.0, slit="SLIT5.0", rospeed="2", **{"eso pro catg": "order_tab_vis"}),
            product("flat-far-order", "mflat", 90.0, **{"eso pro catg": "ORDER_TAB_VIS"}),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("ssol.sof", "spat_solution", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mflat"])

    # ASSERT -- RANKED FIRST, THEN FILTERED, SO THE FAR ORDER TABLE NEVER RE-ENTERS
    assert matched_files(matches) == {"ssol.sof": ["flat-order"]}


def test_mflat_on_nir_lets_a_slit1p5_flat_serve_a_slit5p0_set():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("nir-flat-1.5", "mflat", 104.0, **nir(slit="SLIT1.5")),
            product("nir-flat-0.9", "mflat", 105.0, **nir(slit="SLIT0.9")),
        ]
    )
    downstream = pd.DataFrame(
        [
            downstream_sof("wide.sof", "stare_obj", 105.0, **nir(slit="SLIT5.0")),
            downstream_sof("narrow.sof", "stare_obj", 105.0, **nir(slit="SLIT1.0")),
        ]
    )
    visUpstream = pd.DataFrame([product("vis-flat-1.5", "mflat", 104.0, slit="SLIT1.5")])
    visDownstream = pd.DataFrame([downstream_sof("vis-wide.sof", "stare_obj", 105.0, slit="SLIT5.0")])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mflat"])
    visMatches = match_calibrations(visUpstream, visDownstream, CALIBRATION_RULES["mflat"])

    # ASSERT
    assert matched_files(matches) == {"wide.sof": ["nir-flat-1.5"]}
    assert visMatches.empty


def test_spat_solution_keeps_the_two_nearest_rows_with_no_mode_match():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("ssol-a", "spat_solution", 104.0, binning="2x2", rospeed="4", slit="SLIT5.0"),
            product("ssol-b", "spat_solution", 106.5),
            product("ssol-c", "spat_solution", 90.0),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("stare.sof", "stare_obj", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["spat_solution"])

    # ASSERT
    assert matched_files(matches) == {"stare.sof": ["ssol-a", "ssol-b"]}


def test_std_flux_uses_the_nearest_wide_slit_response_from_any_std_recipe():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("resp-nod", "nod_std", 104.0, slit="SLIT5.0", **{"eso pro catg": "RESP_TAB_VIS"}),
            product("resp-stare-far", "stare_std", 90.0, slit="SLIT5.0", **{"eso pro catg": "RESP_TAB_VIS"}),
            product("resp-narrow", "stare_std", 105.0, slit="SLIT1.0", **{"eso pro catg": "RESP_TAB_VIS"}),
            product("spectrum-wide", "stare_std", 105.0, slit="SLIT5.0", **{"eso pro catg": "SCI_SLIT_FLUX_VIS"}),
        ]
    )
    downstream = pd.DataFrame(
        [
            downstream_sof("stare.sof", "stare_obj", 105.0),
            downstream_sof("std.sof", "offset_std", 105.0),
        ]
    )

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["std_flux"])

    # ASSERT
    assert matched_files(matches) == {"stare.sof": ["resp-nod"]}


def test_a_null_filepath_among_tied_candidates_does_not_break_the_tie_break():
    # ARRANGE
    upstream = pd.DataFrame(
        [
            product("bias-no-path", "mbias", 104.0, filepath=None),
            product("bias-with-path", "mbias", 104.0, filepath="./product/a"),
        ]
    )
    downstream = pd.DataFrame([downstream_sof("dsol.sof", "disp_solution", 105.0)])

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT -- A NULL PATH SORTS FIRST, AS IN SQLITE
    assert matched_files(matches) == {"dsol.sof": ["bias-no-path"]}


def test_rows_with_no_mjd_obs_are_never_matched():
    # ARRANGE -- THE OLD VIEWS SORTED A NULL DISTANCE FIRST; THE MATCHER LEAVES SUCH ROWS OUT
    upstream = pd.DataFrame([product("bias-no-mjd", "mbias", None), product("bias", "mbias", 90.0)])
    downstream = pd.DataFrame(
        [
            downstream_sof("dsol.sof", "disp_solution", 105.0),
            downstream_sof("no-mjd.sof", "disp_solution", None),
        ]
    )

    # ACT
    matches = match_calibrations(upstream, downstream, CALIBRATION_RULES["mbias"])

    # ASSERT
    assert matched_files(matches) == {"dsol.sof": ["bias"]}


def test_sql_like_folds_case_only_for_ascii_letters_as_sqlite_does():
    # ARRANGE
    tags = pd.Series(["order_tab_vis", "İRDER_TAB", None])

    # ACT
    matched = sql_like(tags, "ORDER%").tolist() + sql_like(pd.Series(["İ"]), "i").tolist()

    # ASSERT
    assert matched == [True, False, False, False]
