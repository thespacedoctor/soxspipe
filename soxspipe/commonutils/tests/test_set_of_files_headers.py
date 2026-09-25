from astropy.io import fits

from soxspipe.commonutils.set_of_files import ImageFileCollection


def test_dict_from_fits_header_joins_comment_and_history_cards(tmp_path):
    """Repeated commentary cards are represented by one summary value."""
    header = fits.Header()
    header["OBJECT"] = "test target"
    header.add_comment("first comment")
    header.add_comment("second comment")
    header.add_history("first event")
    header.add_history("second event")
    fits_path = tmp_path / "repeated-commentary.fits"
    fits.PrimaryHDU(header=header).writeto(fits_path)

    collection = ImageFileCollection(
        location=tmp_path, filenames=[fits_path.name]
    )

    summary = collection._dict_from_fits_header(str(fits_path))

    assert summary["comment"] == ["first comment,second comment"]
    assert summary["history"] == ["first event,second event"]
    assert summary["OBJECT"] == ["test target"]
