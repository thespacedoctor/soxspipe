"""Contracts for the packaged SOXS NIR multi-pinhole predicted-line table.

DY-1989, DY-2084, DY-1970, DY-2173 and DY-2171.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.table import Table

import soxspipe
from tests.unit._nir_line_tables import (
    EDITED_ORDER,
    HGI_NM,
    MPH_O10_ROWS_BEFORE,
    MPH_TABLE,
    NEI_NM,
    PINHOLES_PER_LINE,
    POSITION_ATOL_PX,
    SLIT_INDEX_AT_CENTRE,
    SPH_TABLE,
    WAVELENGTH_ATOL_NM,
    XEI_NM,
    ions,
    order_checksums,
    positions_by_slit,
    read_table,
    set_mask,
    set_rows,
    untouched_rows_checksum,
)

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
RESOURCES_DIR = PACKAGE_DIR / "resources"

# DY-1989: O10 AND O11 ARE REBUILT FROM LINES CHECKED ON THE ARC FRAME
CURATED_ORDERS = [10, 11]
# A PINHOLE FURTHER THAN THIS FROM A SMOOTH CURVE ALONG THE SLIT IS A CORRUPT POSITION
MAX_GEOMETRY_RESIDUAL_PX = 0.5
SLIT_GEOMETRY_POLY_DEG = 2
# THE O11 RED HALF HAD NO FULL SET BEFORE DY-1989
O11_RED_HALF_NM = 1800.0
MIN_O11_RED_HALF_LINES = 8
# WITH THE O11 RED HALF NOW CONSTRAINED, THE X FIT NEEDS ONE MORE ORDER AND WAVELENGTH DEGREE
NIR_X_ORDER_DEG = 5
NIR_X_WAVELENGTH_DEG = 6
# DY-2084: XeI AND HgI WAVELENGTHS ARE NIST RITZ VALUES, NOT NIST OBSERVED VALUES; NO LINE MOVES MORE THAN 0.13 nm
# DY-2173: SO ARE THE Ar I AND Ar II WAVELENGTHS
RITZ_IONS = ["XeI", "HgI", "ArI", "ArII"]
MAX_RITZ_SHIFT_NM = 0.13
# (ORDER, ION, NIST OBSERVED nm, NIST RITZ nm). THE ARC MEASURES THE RITZ VALUE FOR THE FIRST SIX
RITZ_LINES = [
    (12, "XeI", 1605.203, 1605.329),
    (13, "XeI", 1465.985, 1466.080),
    (14, "XeI", 1365.649, 1365.706),
    (16, "XeI", 1227.289, 1227.1902),
    (16, "HgI", 1215.9, 1216.008),
    (14, "HgI", 1342.658, 1342.635),
    # THE ARC PREFERS NIST OBSERVED FOR 1719.868 AND 1732.942 (BY 0.3 A), BUT LEFT THERE THEY MAKE THE SET-OFFSET GUARD
    # CLIP 1692.017, AND THE O11 SKY GETS WORSE, SO ALL O11 HgI LINES TAKE RITZ
    (11, "HgI", 1692.017, 1692.027),
    (11, "HgI", 1719.868, 1719.828),
    (11, "HgI", 1732.942, 1732.896),
]
# (ORDER, ION, nm) SETS REMOVED IN DY-2084: THE ARC FEATURE IS NOT THE LISTED LINE, WHOSE NIST WAVELENGTH IS RIGHT
DROPPED_SETS = [
    (13, "NeI", 1534.82),
    (13, "Ar", 1534.88),
    (16, "ArII", 1167.79),
    (17, "ArII", 1167.79),
    (16, "ArII", 1172.05),
    (17, "ArII", 1172.05),
    (16, "ArII", 1202.71),
]
DROPPED_SET_ATOL_NM = 0.006
# DY-2173: (ORDER, ION, nm) SETS REMOVED BECAUSE NO NIST Ar I OR Ar II LINE MATCHES THEM. EACH WAS EITHER MEASURED ON
# BOTH NIR ARCS WITH NO NIST LINE WITHIN 0.25 A OF ITS ARC POSITION, OR MEASURED ON NEITHER ARC WITH NO NIST LINE NEAR
# ITS CATALOGUE WAVELENGTH. THE TABLE BUILDER RECORDS THE REASON FOR EACH SET
NO_NIST_AR_SETS = [
    (22, "ArI", 912.29839),
    (20, "ArII", 965.59884),
    (19, "ArII", 1016.38798),
    (18, "ArI", 1070.098),
    (19, "ArI", 1070.098),
    (18, "ArI", 1088.097),
    (17, "Ar", 1188.45),
    (16, "Ar", 1189.67),
    (17, "Ar", 1189.67),
    (16, "Ar", 1215.14),
    (15, "ArI", 1259.623),
    (16, "ArI", 1259.623),
    (15, "ArI", 1263.849),
    (16, "ArI", 1263.849),
    (15, "ArI", 1302.844),
    (15, "ArI", 1327.265),
    (14, "Ar", 1330.23),
    (15, "Ar", 1330.23),
    (14, "Ar", 1382.83),
    (14, "ArI", 1391.061),
    (14, "Ar", 1425.69),
    (13, "Ar", 1465.03),
    (13, "Ar", 1469.34),
    (13, "ArI", 1544.674),
    (12, "ArI", 1579.316),
    (11, "ArI", 1740.193),
    (11, "ArI", 1744.491),
    (11, "ArI", 1791.473),
    (10, "ArI", 1863.229),
    (10, "ArI", 1874.501),
    (9, "Ar", 2073.36),
]
NO_NIST_AR_ORDER_10_SET_COUNT = sum(1 for order, _, _ in NO_NIST_AR_SETS if order == EDITED_ORDER)
AR_IONS = ["ArI", "ArII"]
# (ORDER, CATALOGUE nm, NIST ION, NIST RITZ nm) OF THE BARE 'Ar' SETS THAT STAY, RELABELLED
RELABELLED_AR_SETS = [
    (17, 1139.37, "ArI", 1139.3703),
    (17, 1158.03, "ArI", 1158.0389),
    (16, 1173.32, "ArI", 1173.3236),
    (17, 1173.32, "ArI", 1173.3236),
    (13, 1463.44, "ArI", 1463.4414),
    (13, 1535.31, "ArI", 1535.3128),
    (12, 1654.93, "ArI", 1654.9307),
]
# EVERY Ar WAVELENGTH IN THE TABLE (nm), BY ION: THE NIST ASD AIR RITZ VALUE AS NIST PRINTS IT, OR THE OBSERVED VALUE
# FOR THE EIGHT Ar I LINES NIST LISTS WITHOUT A RITZ VALUE (1005.206, 1068.1773, 1073.387, 1075.916, 1110.646,
# 1322.8107, 1336.7111 AND 1390.7478)
NIST_AR_NM = {
    ion: [float(text) for text in values.split()]
    for ion, values in {
        "ArI": (
            "794.81759 800.61562 805.33084 810.36921 811.53109 826.45215 838.473 840.82085 842.46473 849.03066 "
            "852.14411 862.04603 866.7943 876.16864 879.90874 896.21468 907.5395 919.4638 922.4498 935.42179 "
            "965.77858 978.4502 1005.206 1047.00528 1047.8034 1050.65 1067.3566 1068.1773 1073.387 1075.916 1095.0726 "
            "1107.8868 1110.646 1139.3703 1144.1832 1148.81076 1158.0389 1166.8709 1173.3236 1211.2324 1213.9737 "
            "1234.3392 1235.6297 1240.2828 1243.9321 1245.6114 1248.7663 1255.4324 1270.228 1280.2737 1293.3196 "
            "1295.6658 1300.8264 1321.3991 1322.8107 1331.3209 1336.7111 1340.6587 1349.9406 1354.4205 1357.3618 "
            "1359.9333 1362.2659 1367.8549 1382.5716 1390.7478 1409.364 1463.4414 1473.9139 1487.6537 1504.6503 "
            "1530.1881 1532.9345 1535.3128 1540.264 1573.4909 1581.6777 1589.9687 1598.9491 1612.2656 1618.0023 "
            "1643.6575 1651.9867 1654.9307 1674.0078 1782.399"
        ),
        "ArII": (
            "801.752791 841.18823 885.07024 906.80229 919.25927 922.2264 933.59814 948.08797 965.59715 978.3081 "
            "1046.71793 1068.30352 1095.42866 1240.70844 1254.5326 1273.86035 1336.47865 1372.41907 1615.23873"
        ),
    }.items()
}
# DY-2171: NIST ASD LISTS NO Ar I LINE BETWEEN 1842.8 AND 2031.7 nm, SO THE TWO O10 Ar I SETS CARRY THE RITZ AIR
# WAVELENGTH OF THEIR NIST Ar I ENERGY LEVELS (PECK AND REEDER AIR, WHICH REPRODUCES THE NIST RITZ 2031.701 nm)
NIST_LEVEL_RITZ_AR_NM = {"ArI": [1981.7507, 1996.5729]}
# A CATALOGUE Ar WAVELENGTH MOVES AT MOST THIS FAR TO ITS NIST VALUE (THE LARGEST MOVE IS 1158.03 -> 1158.0389)
MAX_AR_SHIFT_NM = 0.009

# DY-1970: THREE O10 ArII SETS ARE OTHER LINES, 7 TO 9 A BLUER; THE TABLE SWAPS THEM FOR NeI, HgI AND XeI SETS THAT
# CARRY THE NIST WAVELENGTH OF THE LINE SEEN ON THE ARC
O10_REMOVED_ARII_NM = [1948.118, 1982.643, 1995.242]
# (ION, nm, [(X, Y) BY SLIT INDEX 0 TO 8]) OF THE O10 SETS ADDED IN DY-1970
O10_ADDED_SETS = [
    (
        "XeI",
        XEI_NM,
        [
            (238.94428400821, 384.2126726446532),
            (238.0291672711235, 388.3097409456561),
            (237.01387119679896, 392.4341127104791),
            (236.24280510734806, 396.55301744552594),
            (235.24139291562463, 400.64647062444766),
            (234.42978633556734, 404.72281130535845),
            (233.54475099718334, 408.9019194476989),
            (232.7501887397893, 413.02986753969475),
            (231.87786571163622, 417.10651614419385),
        ],
    ),
    (
        "NeI",
        NEI_NM,
        [
            (1133.1216187967545, 377.3162562167707),
            (1131.9770113938387, 381.3074396418454),
            (1130.9558039850574, 385.1798617035703),
            (1129.9770475807022, 389.10916048943017),
            (1128.92364794532, 393.0473836343341),
            (1127.9042316532855, 396.9720631194816),
            (1126.7902833534977, 400.9025266902693),
            (1125.8810051541302, 404.8512953425118),
            (1124.8812355416576, 408.7863302559909),
        ],
    ),
    (
        "HgI",
        HGI_NM,
        [
            (1278.8037203547938, 391.1359111152702),
            (1277.6575418458676, 394.996090481074),
            (1276.6034385248267, 398.85481309895357),
            (1275.5471443022102, 402.7272151588729),
            (1274.551644963742, 406.6691489833562),
            (1273.4727593716086, 410.4867301338521),
            (1272.3909330669321, 414.5750878728325),
            (1271.3225078958224, 418.3978550030975),
            (1270.2998868208856, 422.28452520525985),
        ],
    ),
]
# DY-2171: (ION, nm, [(X, Y) BY SLIT INDEX 0 TO 8]) OF THE STRONG O10 ARC LINES ADDED AS FULL SETS. THE NeI WAVELENGTHS
# ARE NIST ASD RITZ VALUES; 1862.5158 IS THE PEAK THE CLEAN LIST CALLS NeII 1862.734. POSITIONS COME FROM THE DY-2084
# METHOD (DEVELOP 095625 SSOL MAP PLUS THE MEDIAN TABLE-MINUS-MAP OFFSET OF THE O10 SETS)
O10_STRONG_ARC_LINE_SETS = [
    (
        "NeI",
        1859.1540,
        [
            (24.830957412719727, 409.0667724609375),
            (23.995861053466797, 413.2154541015625),
            (23.160764694213867, 417.36627197265625),
            (22.32566261291504, 421.5193176269531),
            (21.49056625366211, 425.6744384765625),
            (20.65546989440918, 429.8317565917969),
            (19.82036781311035, 433.9912109375),
            (18.985271453857422, 438.15283203125),
            (18.150175094604492, 442.3165588378906),
        ],
    ),
    (
        "NeI",
        1862.5158,
        [
            (61.23123359680176, 404.1734619140625),
            (60.38059043884277, 408.3154296875),
            (59.529951095581055, 412.459716796875),
            (58.67930030822754, 416.60626220703125),
            (57.82866096496582, 420.75506591796875),
            (56.978017807006836, 424.9061279296875),
            (56.127370834350586, 429.0594787597656),
            (55.2767276763916, 433.2150573730469),
            (54.42608451843262, 437.3729248046875),
        ],
    ),
    (
        "ArI",
        1981.7507,
        [
            (1419.8057308197021, 408.1239013671875),
            (1418.7250423431396, 411.96826171875),
            (1417.6443538665771, 415.8207092285156),
            (1416.5635433197021, 419.6813049316406),
            (1415.4828548431396, 423.54998779296875),
            (1414.4021663665771, 427.4267578125),
            (1413.3214778900146, 431.3116455078125),
            (1412.2407894134521, 435.2046203613281),
            (1411.1601009368896, 439.105712890625),
        ],
    ),
    (
        "ArI",
        1996.5729,
        [
            (1600.7295589447021, 435.2734375),
            (1599.6057796478271, 439.0909118652344),
            (1598.4820003509521, 442.9169616699219),
            (1597.3580989837646, 446.7515869140625),
            (1596.2343196868896, 450.59478759765625),
            (1595.1105403900146, 454.4465637207031),
            (1593.9867610931396, 458.3069152832031),
            (1592.8628597259521, 462.17584228515625),
            (1591.7390804290771, 466.0533447265625),
        ],
    ),
]
# DY-2171: THE PEAK AT 1862.5 nm IS NeI 1862.5158, NOT THE NeII 1862.734 OF THE CLEAN LIST
O10_MISIDENTIFIED_NEII_NM = 1862.734
# THE ORDER 10 ROWS THAT ARE NOT IN AN ADDED SET WERE NOT TOUCHED
TOUCHED_SETS = [(ion, wavelength) for ion, wavelength, _ in O10_ADDED_SETS + O10_STRONG_ARC_LINE_SETS]
# SHA-256 OF THE UNTOUCHED ORDER 10 ROWS (EVERY ROW BUT THE ADDED SETS), IN FILE ORDER, RE-RECORDED IN DY-2173
# REGENERATE WITH untouched_rows_checksum(read_table(MPH_TABLE), EDITED_ORDER, TOUCHED_SETS);
# A LATER LEGITIMATE EDIT TO ORDER 10 MUST RE-RECORD IT
UNTOUCHED_O10_CHECKSUM = "d45813e7c8681c928001a6c97b436547ab43ab38561116ef4cf60b0cc742d462"
# SHA-256 OF EVERY ROW OF EACH ORDER OTHER THAN 10, IN FILE ORDER, RE-RECORDED IN DY-2173 (ORDER 9 HELD ONLY THE BARE
# 'Ar' 2073.36 SET, SO IT IS GONE). REGENERATE WITH order_checksums(read_table(MPH_TABLE));
# A LATER LEGITIMATE EDIT TO ANY OF THESE ORDERS MUST RE-RECORD IT
UNEDITED_ORDER_CHECKSUMS = {
    11: "448f47fce8b45d467d500d9633311fcedcf7394191cd594b565b71bd1ecd5eca",
    12: "03ecde37d579c0f24d44cc0f099cd90114742bd9f238d3b056adb9361728bf8f",
    13: "ce3eb3217399b0aed9684425d99cd665730c7d0a08c56b0749b028fc00316ed8",
    14: "f6d049bc5f51ec1ef12fe854692b02d8da64726464255c7743bea10314849860",
    15: "17e8ace2b92e0e001cbbf1adf7827d56eda60059f788928bbb31ad9d6f332aa1",
    16: "394565decf9c6d809d9c0fde303845246ac58a22552f3a4fbcc38118809a61e9",
    17: "4344e57b9167fcb663d7977ec1acf47db537ae0490f4693509e96fb288b74d70",
    18: "3d507871812335184a9c0a1983e30647b2c9bf913c4ba9f9cc2b44d95c770d2e",
    19: "6e49e2dc0077810d789424d8e15a18e270eeb0e130abbac162fcb45877ccdd94",
    20: "8f46580808eae99f302d938d9380f67b398ed1ed9c3ee78d617af1be7914fe8d",
    21: "54402b0661dbe732e5394c5321ccc56f15baf1d7d4e0fb048659a47c3cca2071",
    22: "6497b2ebcc38080a51db3127ea1b41031c3b430b3babcb810d2d89dffe5247db",
    23: "ed364e57809d2596408af2ea1d7afd1f31ad6ac19f1996b963a2d90273d3edeb",
    24: "866fe14d5f2bf67fcc0bcee0c37a8d3f408c55170d4524996bccb1dde572a2fb",
}


def _uncurated_rows(name: str) -> Table:
    table = read_table(name)
    return table[~np.isin(table["order"], CURATED_ORDERS)]


def _stable_sorted(table: Table) -> Table:
    # THE SHIPPED TABLE REPEATS SOME (ORDER, WAVELENGTH, SLIT_INDEX) KEYS, SO POSITIONS BREAK THE TIES
    rows = table.copy()
    rows.sort(["order", "slit_index", "detector_x", "detector_y"], kind="stable")
    return rows


def _dropped_set_rows(table: Table) -> np.ndarray:
    rows = np.zeros(len(table), dtype=bool)
    for order, ion, wavelength in DROPPED_SETS + NO_NIST_AR_SETS:
        rows |= set_mask(table, order, ion, wavelength, DROPPED_SET_ATOL_NM)
    return rows


def _ar_rows(table: Table) -> np.ndarray:
    return np.char.startswith(ions(table), "Ar")


def _clean_list_rows_kept() -> Table:
    table = _uncurated_rows(SPH_TABLE)
    return table[~_dropped_set_rows(table)]


def _nir_spatial_solution_settings() -> dict:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return settings["soxs-spat-solution"]["nir"]


def _nir_line_tables() -> dict:
    detectorParams = yaml.safe_load((RESOURCES_DIR / "soxs_detector_parameters.yaml").read_text())
    return detectorParams["NIR"]["predicted pinhole lines"]


def test_nir_multi_pinhole_frames_read_the_curated_table() -> None:
    # ACT
    tables = _nir_line_tables()

    # ASSERT
    assert tables["multi"]["1x1"] == MPH_TABLE
    assert tables["single"]["1x1"] == SPH_TABLE


def test_nir_spatial_solution_x_fit_uses_the_raised_degrees() -> None:
    # ACT
    settings = _nir_spatial_solution_settings()

    # ASSERT
    assert settings["order-deg"][0] == NIR_X_ORDER_DEG
    assert settings["wavelength-deg"][0] == NIR_X_WAVELENGTH_DEG


@pytest.mark.parametrize("order", CURATED_ORDERS)
def test_every_curated_line_is_a_full_pinhole_set(order: int) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    rows = table[table["order"] == order]

    # ACT
    wavelengths, counts = np.unique(np.asarray(rows["wavelength"]), return_counts=True)

    # ASSERT
    assert len(wavelengths) > 0
    assert np.all(counts == PINHOLES_PER_LINE)
    for wavelength in wavelengths:
        line = rows[rows["wavelength"] == wavelength]
        assert sorted(line["slit_index"]) == list(range(PINHOLES_PER_LINE))


@pytest.mark.parametrize("order", CURATED_ORDERS)
def test_curated_pinhole_positions_follow_a_smooth_slit_geometry(order: int) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    rows = table[table["order"] == order]

    # ACT
    offenders = []
    for wavelength in np.unique(np.asarray(rows["wavelength"])):
        line = rows[rows["wavelength"] == wavelength]
        slit = np.asarray(line["slit_position"])
        for axis in ("detector_x", "detector_y"):
            position = np.asarray(line[axis])
            residual = position - np.polyval(np.polyfit(slit, position, SLIT_GEOMETRY_POLY_DEG), slit)
            if np.max(np.abs(residual)) >= MAX_GEOMETRY_RESIDUAL_PX:
                offenders.append(f"O{order} {wavelength} {axis}")

    # ASSERT
    assert offenders == []


def test_o11_red_half_has_enough_lines() -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    rows = table[(table["order"] == 11) & (table["wavelength"] > O11_RED_HALF_NM)]

    # ACT
    redLines = np.unique(np.asarray(rows["wavelength"]))

    # ASSERT
    assert len(redLines) >= MIN_O11_RED_HALF_LINES


@pytest.mark.parametrize(("order", "ion", "wavelength"), DROPPED_SETS)
def test_misidentified_set_is_removed(order: int, ion: str, wavelength: float) -> None:
    # ARRANGE
    clean = _uncurated_rows(SPH_TABLE)
    curated = read_table(MPH_TABLE)

    # ACT
    inClean = set_mask(clean, order, ion, wavelength, DROPPED_SET_ATOL_NM).sum()
    inCurated = set_mask(curated, order, ion, wavelength, DROPPED_SET_ATOL_NM).sum()

    # ASSERT
    assert inClean == PINHOLES_PER_LINE
    assert inCurated == 0


@pytest.mark.parametrize(("order", "ion", "wavelength"), NO_NIST_AR_SETS)
def test_ar_set_with_no_nist_line_is_removed(order: int, ion: str, wavelength: float) -> None:
    # ARRANGE
    clean = read_table(SPH_TABLE)
    curated = read_table(MPH_TABLE)
    nearby = (np.asarray(curated["order"]) == order) & np.isclose(
        np.asarray(curated["wavelength"], dtype=float), wavelength, rtol=0, atol=DROPPED_SET_ATOL_NM
    )

    # ACT
    inClean = set_mask(clean, order, ion, wavelength, DROPPED_SET_ATOL_NM).sum()
    arInCurated = int(np.sum(nearby & _ar_rows(curated)))

    # ASSERT
    assert inClean > 0
    assert arInCurated == 0


def test_no_set_keeps_the_bare_ar_label() -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    labels = set(ions(table[_ar_rows(table)]))

    # ASSERT
    assert labels == set(AR_IONS)


@pytest.mark.parametrize("ion", AR_IONS)
def test_every_ar_wavelength_is_a_nist_value(ion: str) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    wavelengths = sorted(set(np.asarray(table["wavelength"], dtype=float)[ions(table) == ion].tolist()))

    # ASSERT
    # THE TABLE STORES FLOAT64 PARSED FROM THE NIST TEXT, SO EACH VALUE MATCHES EXACTLY;
    # A SORTED LIST NAMES THE LINE THAT DIFFERS
    assert wavelengths == sorted(NIST_AR_NM[ion] + NIST_LEVEL_RITZ_AR_NM.get(ion, []))


@pytest.mark.parametrize(("order", "catalogue", "ion", "ritz"), RELABELLED_AR_SETS)
def test_bare_ar_set_takes_its_nist_ion_and_ritz_wavelength(
    order: int, catalogue: float, ion: str, ritz: float
) -> None:
    # ARRANGE
    clean = read_table(SPH_TABLE)
    curated = read_table(MPH_TABLE)

    # ACT
    before = set_rows(clean, order, "Ar", catalogue)
    after = set_rows(curated, order, ion, ritz, atol=0)

    # ASSERT
    assert len(before) > 0
    assert len(after) == len(before)
    assert len(set_rows(curated, order, "Ar", catalogue)) == 0


def test_ar_wavelengths_move_no_more_than_the_catalogue_rounding() -> None:
    # ARRANGE
    curated = _stable_sorted(_uncurated_rows(MPH_TABLE))
    shipped = _stable_sorted(_clean_list_rows_kept())
    ar = _ar_rows(curated)

    # ACT
    shift = np.asarray(curated["wavelength"], dtype=float) - np.asarray(shipped["wavelength"], dtype=float)

    # ASSERT
    assert np.any(shift[ar] != 0)
    assert np.max(np.abs(shift[ar])) <= MAX_AR_SHIFT_NM


def test_only_bare_ar_labels_differ_from_the_clean_list() -> None:
    # ARRANGE
    curated = _stable_sorted(_uncurated_rows(MPH_TABLE))
    shipped = _stable_sorted(_clean_list_rows_kept())

    # ACT
    relabelled = ions(curated) != ions(shipped)

    # ASSERT
    assert relabelled.any()
    assert set(ions(shipped[relabelled])) == {"Ar"}
    assert set(ions(curated[relabelled])) <= set(AR_IONS)


def test_orders_outside_the_curation_keep_the_clean_list_geometry() -> None:
    # ARRANGE
    columns = [name for name in read_table(MPH_TABLE).colnames if name not in ("wavelength", "ion")]

    # ACT
    curatedRest = _stable_sorted(_uncurated_rows(MPH_TABLE))
    shippedRest = _stable_sorted(_clean_list_rows_kept())

    # ASSERT
    assert curatedRest.colnames == shippedRest.colnames
    assert len(curatedRest) == len(shippedRest)
    for column in columns:
        assert np.array_equal(np.asarray(curatedRest[column]), np.asarray(shippedRest[column])), column


def test_only_ritz_ion_wavelengths_differ_from_the_clean_list() -> None:
    # ARRANGE
    curated = _stable_sorted(_uncurated_rows(MPH_TABLE))
    shipped = _stable_sorted(_clean_list_rows_kept())

    # ACT
    shift = np.asarray(curated["wavelength"], dtype=float) - np.asarray(shipped["wavelength"], dtype=float)
    moved = shift != 0

    # ASSERT
    assert moved.any()
    assert set(ions(curated[moved])) <= set(RITZ_IONS)
    assert np.max(np.abs(shift)) <= MAX_RITZ_SHIFT_NM


@pytest.mark.parametrize(("order", "ion", "observed", "ritz"), RITZ_LINES)
def test_line_carries_its_nist_ritz_wavelength(order: int, ion: str, observed: float, ritz: float) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    rows = table[(table["order"] == order) & (ions(table) == ion)]

    # ACT
    wavelengths = np.unique(np.asarray(rows["wavelength"], dtype=float))

    # ASSERT
    assert np.any(np.isclose(wavelengths, ritz, rtol=0, atol=WAVELENGTH_ATOL_NM))
    assert not np.any(np.isclose(wavelengths, observed, rtol=0, atol=WAVELENGTH_ATOL_NM))


@pytest.mark.parametrize("wavelength", O10_REMOVED_ARII_NM)
def test_mislabelled_arii_set_is_absent_from_order_10(wavelength: float) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    rows = set_rows(table, EDITED_ORDER, "ArII", wavelength)

    # ASSERT
    assert len(rows) == 0


@pytest.mark.parametrize(("ion", "wavelength", "positions"), O10_ADDED_SETS + O10_STRONG_ARC_LINE_SETS)
def test_added_order_10_line_is_a_full_set_with_its_measured_centre(
    ion: str, wavelength: float, positions: list[tuple[float, float]]
) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    centre = positions[SLIT_INDEX_AT_CENTRE]

    # ACT
    rows = set_rows(table, EDITED_ORDER, ion, wavelength)
    middle = rows[rows["slit_index"] == SLIT_INDEX_AT_CENTRE]

    # ASSERT
    assert sorted(rows["slit_index"]) == list(range(PINHOLES_PER_LINE))
    assert len(middle) == 1
    assert (float(middle["detector_x"][0]), float(middle["detector_y"][0])) == pytest.approx(
        centre, rel=0, abs=POSITION_ATOL_PX
    )


@pytest.mark.parametrize(("ion", "wavelength", "positions"), O10_ADDED_SETS + O10_STRONG_ARC_LINE_SETS)
def test_added_order_10_line_has_every_pinhole_at_its_recorded_position(
    ion: str, wavelength: float, positions: list[tuple[float, float]]
) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    rows = set_rows(table, EDITED_ORDER, ion, wavelength)
    actual = positions_by_slit(rows)

    # ASSERT
    assert len(rows) == PINHOLES_PER_LINE
    assert actual.keys() == set(range(PINHOLES_PER_LINE))
    for slit, position in enumerate(positions):
        assert actual[slit] == pytest.approx(position, rel=0, abs=POSITION_ATOL_PX)


def test_order_10_row_count_changes_only_by_the_removed_and_added_sets() -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)
    # DY-1970 SWAPPED THREE SETS OUT AND THREE IN; DY-2173 REMOVED TWO Ar I SETS; DY-2171 ADDED FOUR SETS
    expected = MPH_O10_ROWS_BEFORE + (len(O10_STRONG_ARC_LINE_SETS) - NO_NIST_AR_ORDER_10_SET_COUNT) * PINHOLES_PER_LINE

    # ACT
    rowCount = int(np.sum(np.asarray(table["order"]) == EDITED_ORDER))

    # ASSERT
    assert rowCount == expected


def test_order_10_rows_the_edit_did_not_touch_match_the_recorded_checksum() -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    checksum = untouched_rows_checksum(table, EDITED_ORDER, TOUCHED_SETS)

    # ASSERT
    assert checksum == UNTOUCHED_O10_CHECKSUM


@pytest.mark.parametrize("order", sorted(UNEDITED_ORDER_CHECKSUMS))
def test_orders_other_than_10_match_develop_row_for_row(order: int) -> None:
    # ARRANGE
    table = read_table(MPH_TABLE)

    # ACT
    checksums = order_checksums(table)

    # ASSERT
    assert checksums[order] == UNEDITED_ORDER_CHECKSUMS[order]


def test_order_10_does_not_carry_the_misidentified_neii_1862_734_line() -> None:
    # GUARDS A REBUILD FROM THE CLEAN LIST, WHICH LABELS THE NeI 1862.5158 PEAK AS NeII 1862.734
    # ARRANGE
    table = read_table(MPH_TABLE)
    order10 = table[np.asarray(table["order"]) == EDITED_ORDER]

    # ACT
    neii = np.isclose(
        np.asarray(order10["wavelength"], dtype=float), O10_MISIDENTIFIED_NEII_NM, rtol=0, atol=WAVELENGTH_ATOL_NM
    )

    # ASSERT
    assert not neii.any()
