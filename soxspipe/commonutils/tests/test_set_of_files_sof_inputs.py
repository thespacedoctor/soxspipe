from unittest.mock import Mock, patch

import pytest

from soxspipe.commonutils.set_of_files import set_of_files


@pytest.mark.parametrize("supplementary_in_frame_directory", [True, False])
def test_sof_collection_contains_only_fits_files(
    tmp_path, supplementary_in_frame_directory
):
    """Supplementary SOF members must not be passed to ImageFileCollection."""
    frame_directory = tmp_path / "frames"
    frame_directory.mkdir()
    fits_path = frame_directory / "a.fits"
    fits_path.touch()

    supplementary_directory = (
        frame_directory if supplementary_in_frame_directory else tmp_path / "products"
    )
    supplementary_directory.mkdir(exist_ok=True)
    supplementary_path = supplementary_directory / "VIS_ORDER_LOCATIONS.csv"
    supplementary_path.touch()

    sof_path = tmp_path / "input.sof"
    sof_path.write_text(f"{fits_path} BIAS_VIS\n{supplementary_path}\n")

    worker = set_of_files.__new__(set_of_files)
    worker.log = Mock()
    worker.inputFrames = str(sof_path)
    worker.currentSession = None
    worker.ext = 0
    worker.keys = ["file"]

    collection = Mock()
    with patch(
        "soxspipe.commonutils.set_of_files.ImageFileCollection",
        return_value=collection,
    ) as image_file_collection:
        result, supplementary = worker.get()

    assert result is collection
    image_file_collection.assert_called_once_with(
        filenames=[fits_path.name],
        keywords=worker.keys,
        location=str(frame_directory),
        ext=0,
    )
    assert supplementary == {
        "VIS": {"ORDER_LOCATIONS": str(supplementary_path)}
    }
