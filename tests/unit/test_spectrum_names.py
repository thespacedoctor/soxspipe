"""The suffix and labels of the merged spectrum are the names the sof map and the recipes share (DY-2346)."""

from __future__ import annotations

import pytest

from soxspipe.commonutils.spectrum_names import SPECTRUM_ASCII_LABEL, SPECTRUM_SUFFIX, SPECTRUM_TABLE_LABEL

pytestmark = pytest.mark.unit


def test_the_suffix_is_the_one_the_sof_map_predicts() -> None:
    from pathlib import Path

    import yaml

    import soxspipe

    sofMap = yaml.safe_load((Path(soxspipe.__file__).parent / "resources" / "soxs_sof_map.yaml").read_text())
    predictedSuffixes = {
        table["replace"][0]["to"]
        for entry in sofMap.values()
        if isinstance(entry, dict)
        for table in (product.get("fits table", {}) for product in entry.get("products") or [])
        if table.get("eso pro catg") == "SCI_SLIT_FLUX"
    }

    assert len(predictedSuffixes) > 0
    assert predictedSuffixes == {f"{SPECTRUM_SUFFIX}.fits"}


def test_the_labels_name_the_spectrum_table_and_its_ascii_file() -> None:
    assert (SPECTRUM_TABLE_LABEL, SPECTRUM_ASCII_LABEL) == ("SPECTRUM_TABLE", "SPECTRUM_ASCII")
