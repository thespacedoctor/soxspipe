"""A bare SOXS extractor that runs the real order merge on the log-wavelength grid (DY-2348)."""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

from astropy.io import fits

from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.keyword_lookup import keyword_lookup


def soxs_merge_extractor(
    log: Any,
    arm: str = "NIR",
    cards: dict[str, Any] | None = None,
    notFlattened: str = "",
    carrySkyCounts: bool = True,
) -> horne_extraction:
    """Return an extractor with only what ``merge_extracted_orders`` reads for a SOXS arm.

    **Key Arguments:**

    - ``log`` -- the logger
    - ``arm`` -- ``VIS`` or ``NIR``
    - ``cards`` -- the header cards of the sky-subtracted frame, such as ``ESO DET BINY``
    - ``notFlattened`` -- ``_NOTFLAT`` for the not-flattened spectrum, which carries a flux density column
    - ``carrySkyCounts`` -- whether the merge carries the sky column

    **Return:**

    - ``extractor`` -- a ``horne_extraction`` object that has no other state than ``merge_extracted_orders`` needs

    **Usage:**

    ```python
    merged, joins = soxs_merge_extractor(log, "NIR").merge_extracted_orders(extractedOrdersDF)
    ```
    """
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.arm = arm
    extractor.inst = "SOXS"
    extractor.kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get
    extractor.detectorParams = detector_lookup(log=log, settings={"instrument": "soxs"}).get(arm)
    extractor.skySubtractedFrame = SimpleNamespace(header=fits.Header(cards or {}))
    extractor.notFlattened = notFlattened
    extractor.carrySkyCounts = carrySkyCounts
    extractor.recipeSettings = {}
    return extractor
