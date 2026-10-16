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
                if k not in ["axis", "order_deg", "wavelength_deg", "slit_deg"] and not math.isnan(v)
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
    """*Use a dispersion solution to convert wavelength, slit-position and echelle order
    numbers to X,Y pixel positions.*

    Return a line-list with x,y fits given a first guess dispersion map.*

    **Key Arguments:**

    - ``log`` -- logger
    - ``dispersionMapPath`` -- path to the dispersion map
    - ``orderPixelTable`` -- a data-frame including 'order', 'wavelength' and 'slit_pos' columns
    - ``removeOffDetectorLocation`` -- if data points are found to lie off the detector plane then remove
      them from the results. Default *True*
    - ``trimColumns`` -- if True, only return the specified columns. Default *False*

    **Return:**

    - ``orderPixelTable`` -- the input data-frame with ``fit_x``/``fit_y`` pixel columns added, filtered
      and/or trimmed per the arguments above

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
    coeff, poly = _load_axis_polynomials(log, axesData)

    # CONVERT THE ORDER-SORTED WAVELENGTH ARRAYS INTO ARRAYS OF PIXEL TUPLES
    orderPixelTable["fit_x"] = poly["x"](orderPixelTable, *coeff["x"])
    orderPixelTable["fit_y"] = poly["y"](orderPixelTable, *coeff["y"])

    orderPixelTable["fit_x"] = orderPixelTable["fit_x"].astype(np.float32)
    orderPixelTable["fit_y"] = orderPixelTable["fit_y"].astype(np.float32)

    orderPixelTable = _filter_and_trim_pixel_table(orderPixelTable, removeOffDetectorLocation, trimColumns)

    log.debug("completed the ``dispersion_map_to_pixel_arrays`` function")
    return orderPixelTable


def _load_axis_polynomials(log, axesData):
    """*build per-axis coefficient tuples and fitted polynomial callables from parsed dispersion-map axis data*

    **Key Arguments:**

    - ``log`` -- logger
    - ``axesData`` -- dictionary keyed by axis ("x"/"y") of {"orderDeg", "wavelengthDeg", "slitDeg", "coeff"},
      as returned by ``_read_dispersion_map_axes``

    **Return:**

    - ``coeff`` -- dictionary keyed by axis of the flat coefficient tuple
    - ``poly`` -- dictionary keyed by axis of the fitted polynomial callable

    **Usage:**

    ```python
    coeff, poly = _load_axis_polynomials(log, axesData)
    ```
    """
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
    return coeff, poly


def _filter_and_trim_pixel_table(orderPixelTable, removeOffDetectorLocation, trimColumns):
    """*remove off-detector rows and/or trim to the reporting columns*

    **Key Arguments:**

    - ``orderPixelTable`` -- the data-frame with fitted ``fit_x``/``fit_y`` columns already added
    - ``removeOffDetectorLocation`` -- if data points are found to lie off the detector plane then remove
      them from the results
    - ``trimColumns`` -- if True, only return the specified columns

    **Return:**

    - ``orderPixelTable`` -- the filtered and/or trimmed data-frame

    **Usage:**

    ```python
    orderPixelTable = _filter_and_trim_pixel_table(orderPixelTable, removeOffDetectorLocation, trimColumns)
    ```
    """
    if removeOffDetectorLocation:
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable["fit_x"] > 0) & (orderPixelTable["fit_y"] > 0)
        orderPixelTable = orderPixelTable.loc[mask]

    if trimColumns:
        orderPixelTable = orderPixelTable[["order", "wavelength", "slit_position", "fit_x", "fit_y"]]

    return orderPixelTable


def default_coefficients(orderDeg, wavelengthDeg, slitDeg):
    """*build unit-valued fallback coefficient arrays sized for the given polynomial degrees*

    **Key Arguments:**

    - ``orderDeg`` -- the order deg, or a list of per-axis order degs
    - ``wavelengthDeg`` -- wavelength degree, or a list of per-axis wavelength degrees
    - ``slitDeg`` -- slit degree, or a list of per-axis slit degrees

    **Return:**

    - ``coeff`` -- dictionary keyed by axis ("x"/"y") of an array of ones sized for that axis's polynomial

    **Usage:**

    ```python
    coeff = default_coefficients(orderDeg, wavelengthDeg, slitDeg)
    ```
    """
    import numpy as np

    coeff = {}
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
        coeff["x"] = np.ones((orderDeg + 1) * (wavelengthDeg + 1) * (slitDeg + 1), dtype=float)
        coeff["y"] = np.ones((orderDeg + 1) * (wavelengthDeg + 1) * (slitDeg + 1), dtype=float)
    return coeff
