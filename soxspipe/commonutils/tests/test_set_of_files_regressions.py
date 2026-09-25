from unittest.mock import Mock

from astropy.io import fits

from soxspipe.commonutils.set_of_files import set_of_files


def test_sof_supplementary_files_are_not_added_to_image_collection(tmp_path):
    """Only FITS members of an SOF should determine the image collection."""
    frames_dir = tmp_path / "frames"
    supplementary_dir = tmp_path / "calibrations"
    frames_dir.mkdir()
    supplementary_dir.mkdir()

    frame_path = frames_dir / "a.fits"
    fits.PrimaryHDU().writeto(frame_path)
    order_table_path = supplementary_dir / "VIS_ORDER_LOCATIONS.csv"
    order_table_path.write_text("order,centre\n1,10\n")

    sof_path = tmp_path / "input.sof"
    sof_path.write_text(f"{frame_path} BIAS_VIS\n{order_table_path}\n")

    sof_reader = set_of_files.__new__(set_of_files)
    sof_reader.log = Mock()
    sof_reader.inputFrames = str(sof_path)
    sof_reader.currentSession = None
    sof_reader.ext = 0
    sof_reader.keys = ["file"]

    collection, supplementary = sof_reader.get()

    assert collection.location == str(frames_dir)
    assert collection.files == [frame_path.name]
    assert list(collection.summary["file"]) == [frame_path.name]
    assert supplementary == {
        "VIS": {"ORDER_LOCATIONS": str(order_table_path)},
    }
