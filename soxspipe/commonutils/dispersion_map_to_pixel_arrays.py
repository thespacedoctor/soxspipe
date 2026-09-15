#!/usr/bin/env python
"""
*use a first-guess dispersion map to convert wavelengths to pixels*

Author
: David Young

Date Created
: April 15, 2021
"""

import os
from functools import lru_cache
from os.path import expanduser

from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials

os.environ["TERM"] = "vt100"


@lru_cache(maxsize=32)
def _read_dispersion_map_axes(resolvedPath, mtime):
    """*parse a dispersion-map FITS file's per-axis polynomial degrees and flat coefficient list, cached per (resolved path, mtime)*

    **Key Arguments:**

    - ``resolvedPath`` -- the fully resolved (realpath) path to the dispersion-map FITS file
    - ``mtime`` -- the file's last-modified time (cache is keyed on this so it self-invalidates if the file changes)

    **Return:**

    - ``axesData`` -- dictionary keyed by axis ("x"/"y") of {"orderDeg", "wavelengthDeg", "slitDeg", "coeff"}

    **Usage:**

    ```python
    axesData = _read_dispersion_map_axes(resolvedPath, mtime)
    ```
    """
    import math

    from astropy.table import Table

    dat = Table.read(resolvedPath, format="fits")
    tableData = dat.to_pandas()

    axesData = {}
    for index, row in tableData.iterrows():
        axis = row["axis"].decode("utf-8")
        axesData[axis] = {
            "orderDeg": int(row["order_deg"]),
            "wavelengthDeg": int(row["wavelength_deg"]),
            "slitDeg": int(row["slit_deg"]),
            "coeff": tuple(
                float(v)
                for k, v in row.items()
                if k not in ["axis", "order_deg", "wavelength_deg", "slit_deg"]
                and not math.isnan(v)
            ),
        }
    return axesData


def dispersion_map_to_pixel_arrays(
    log,
    dispersionMapPath,
    orderPixelTable,
    removeOffDetectorLocation=True,
    trimColumns=False,
):
    """*Use a dispersion solution to convert wavelength, slit-position and echelle order numbers to X,Y pixel positions.*

    Return a line-list with x,y fits given a first guess dispersion map.*

    **Key Arguments:**

    - `log` -- logger
    - `dispersionMapPath` -- path to the dispersion map
    - `orderPixelTable` -- a data-frame including 'order', 'wavelength' and 'slit_pos' columns
    - `removeOffDetectorLocation` -- if data points are found to lie off the detector plane then remove them from the results. Default *True*
    - `trimColumns` -- if True, only return the specified columns. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils import dispersion_map_to_pixel_arrays
    myDict = {
        "order": [11, 11, 11],
        "wavelength": [850.3, 894.3, 983.2],
        "slit_position": [0, 0, 0]
    }
    orderPixelTable = pd.DataFrame(myDict)
    orderPixelTable = dispersion_map_to_pixel_arrays(
        log=log,
        dispersionMapPath="/path/to/map.csv",
        orderPixelTable=orderPixelTable
    )
    ```
    """
    log.debug("starting the ``dispersion_map_to_pixel_arrays`` function")

    import numpy as np

    # READ THE FILE (CACHED PER RESOLVED PATH + MTIME SO REPEATED CALLS
    # AGAINST THE SAME MAP SKIP DISK I/O AND RE-PARSING)
    home = expanduser("~")
    dispersion_map = dispersionMapPath.replace("~", home)
    resolvedPath = os.path.realpath(dispersion_map)
    mtime = os.path.getmtime(resolvedPath)

    axesData = _read_dispersion_map_axes(resolvedPath, mtime)

    # READ IN THE X- AND Y- GENERATING POLYNOMIALS FROM DISPERSION MAP FILE
    coeff = {}
    poly = {}
    for axis in ("x", "y"):
        d = axesData[axis]
        coeff[axis] = d["coeff"]
        poly[axis] = chebyshev_order_wavelength_polynomials(
            log=log,
            orderDeg=d["orderDeg"],
            wavelengthDeg=d["wavelengthDeg"],
            slitDeg=d["slitDeg"],
            exponentsIncluded=False,
            axis=axis,
        ).poly

    # CONVERT THE ORDER-SORTED WAVELENGTH ARRAYS INTO ARRAYS OF PIXEL TUPLES
    orderPixelTable["fit_x"] = poly["x"](orderPixelTable, *coeff["x"])
    orderPixelTable["fit_y"] = poly["y"](orderPixelTable, *coeff["y"])

    orderPixelTable["fit_x"] = orderPixelTable["fit_x"].astype(np.float32)
    orderPixelTable["fit_y"] = orderPixelTable["fit_y"].astype(np.float32)

    if removeOffDetectorLocation:
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable["fit_x"] > 0) & (orderPixelTable["fit_y"] > 0)
        orderPixelTable = orderPixelTable.loc[mask]

    if trimColumns:
        orderPixelTable = orderPixelTable[
            ["order", "wavelength", "slit_position", "fit_x", "fit_y"]
        ]

    log.debug("completed the ``dispersion_map_to_pixel_arrays`` function")
    return orderPixelTable


def get_cached_coeffs(
    log, arm, settings, recipeName, orderDeg, wavelengthDeg, slitDeg, reset=False
):
    """*find cached coefficients (if they exist)*

    Return a line-list with x,y fits given a first guess dispersion map.*

    **Key Arguments:**

    - ``log`` -- logger
    - ``arm`` -- the spectrograph arm.
    - ``settings`` pipeline settings dictionary
    - ``recipeName`` -- the name of the recipe.
    - ``orderDeg`` -- the order deg
    - ``wavelengthDeg`` -- wavelength degree
    - ``slitDeg`` -- slit degree
    - ``reset`` -- always reset the coeffs. Don't use cached. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils import get_cached_coeffs
    xcoeff, ycoeff = get_cached_coeffs(
        log=log,
        arm=arm,
        settings=settings,
        recipeName=recipeName,
        orderDeg=orderDeg,
        wavelengthDeg=wavelengthDeg,
        slitDeg=slitDeg
    )
    ```
    """
    log.debug("starting the ``get_cached_coeffs`` function")

    import math

    import numpy as np
    from astropy.table import Table

    # READ THE FILE
    home = expanduser("~")
    cache = settings["workspace-root-dir"].replace("~", home) + "/.cache"
    polyOrders = [orderDeg, wavelengthDeg, slitDeg]
    if isinstance(orderDeg, list):
        merged_list = []
        for sublist in polyOrders:
            merged_list.extend(sublist)
        polyOrders = merged_list
    polyOrders[:] = [str(l) for l in polyOrders]
    polyOrders = "".join(polyOrders)
    filename = f"{recipeName}_{arm}_{polyOrders}.fits"
    filePath = f"{cache}/{filename}"

    coeff = {}

    if os.path.exists(filePath) and reset == False:
        dispersion_map = filePath
        # SPEC FORMAT TO PANDAS DATAFRAME
        dat = Table.read(dispersion_map, format="fits")
        tableData = dat.to_pandas()

        # READ IN THE X- AND Y- COEFF FROM DISPERSION MAP FILE
        for index, row in tableData.iterrows():
            axis = row["axis"].decode("utf-8")
            orderDeg = int(row["order_deg"])
            wavelengthDeg = int(row["wavelength_deg"])
            slitDeg = int(row["slit_deg"])
            coeff[axis] = [
                float(v)
                for k, v in row.items()
                if k not in ["axis", "order_deg", "wavelength_deg", "slit_deg"]
                and not math.isnan(v)
            ]
    else:
        if isinstance(orderDeg, list):
            coeff["x"] = np.ones(
                (orderDeg[0] + 1) * (wavelengthDeg[0] + 1) * (slitDeg[0] + 1),
                dtype=float,
            )
            coeff["y"] = np.ones(
                (orderDeg[1] + 1) * (wavelengthDeg[1] + 1) * (slitDeg[1] + 1),
                dtype=float,
            )
        else:
            coeff["x"] = np.ones(
                (orderDeg + 1) * (wavelengthDeg + 1) * (slitDeg + 1), dtype=float
            )
            coeff["y"] = np.ones(
                (orderDeg + 1) * (wavelengthDeg + 1) * (slitDeg + 1), dtype=float
            )

    log.debug("completed the ``get_cached_coeffs`` function")
    return coeff["x"], coeff["y"]
