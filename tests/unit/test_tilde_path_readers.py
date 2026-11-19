"""Readers of workspace and product paths expand only a leading `~`, as the writers do (DY-2287)."""

from __future__ import annotations

from pathlib import Path

import pytest

from soxspipe.commonutils import toolkit
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import read_dispersion_map_axes
from soxspipe.commonutils.scratch import absolute_path
from soxspipe.commonutils.set_of_files import _fits_path_from_sof_line, _supplementary_path_from_sof_line

pytestmark = pytest.mark.unit


@pytest.fixture(name="homeDir")
def fake_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A home directory that is not the real one, with HOME pointing at it."""
    home = (tmp_path / "tilde-home").resolve()
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    return home


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        ("/data/a~b/x.fits", "/data/a~b/x.fits"),
        ("/data/a~/x.fits", "/data/a~/x.fits"),
        ("~/x.fits", "{home}/x.fits"),
        ("~/a~b/x.fits", "{home}/a~b/x.fits"),
        ("~nosuchuser-zz/x.fits", "~nosuchuser-zz/x.fits"),
    ],
)
def test_absolute_path_expands_only_a_leading_tilde(homeDir: Path, path: str, expected: str) -> None:
    assert absolute_path(path) == expected.format(home=homeDir)


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        ("/data/a~/x.fits BIAS_VIS", "/data/a~/x.fits"),
        ("~/x.fits BIAS_VIS", "{home}/x.fits"),
        ("~/a~/x.fits BIAS_VIS", "{home}/a~/x.fits"),
    ],
)
def test_sof_frame_path_expands_only_a_leading_tilde(homeDir: Path, line: str, expected: str) -> None:
    assert _fits_path_from_sof_line(line) == expected.format(home=homeDir)


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        ("/data/a~/x.csv", "/data/a~/x.csv"),
        ("~/x.csv ORDER_TAB_VIS", "{home}/x.csv"),
        ("~/a~/x.csv", "{home}/a~/x.csv"),
    ],
)
def test_sof_supplementary_path_expands_only_a_leading_tilde(homeDir: Path, line: str, expected: str) -> None:
    assert _supplementary_path_from_sof_line(line) == expected.format(home=homeDir)


@pytest.mark.parametrize(
    ("givenPath", "expectedPath"),
    [("{tmp}/a~b/missing.fits", "{tmp}/a~b/missing.fits"), ("~/a~b/missing.fits", "{home}/a~b/missing.fits")],
)
@pytest.mark.parametrize("reader", ["order_table", "disp_map_dataframe", "disp_map_axes"])
def test_file_readers_expand_only_a_leading_tilde(
    homeDir: Path, tmp_path: Path, log: object, reader: str, givenPath: str, expectedPath: str
) -> None:
    given = givenPath.format(tmp=tmp_path.resolve())
    expected = expectedPath.format(tmp=tmp_path.resolve(), home=homeDir)
    readers = {
        "order_table": lambda: toolkit.unpack_order_table(log=log, orderTablePath=given),
        "disp_map_dataframe": lambda: toolkit.twoD_disp_map_image_to_dataframe(
            log=log, slit_length=11, twoDMapPath=given, associatedFrame=False, kw=lambda key: key
        ),
        "disp_map_axes": lambda: read_dispersion_map_axes(given),
    }

    with pytest.raises(FileNotFoundError) as raised:
        readers[reader]()

    assert expected in str(raised.value)
