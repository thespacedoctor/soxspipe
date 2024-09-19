#!/usr/bin/env python
# encoding: utf-8
"""
*use a first-guess dispersion map to convert wavelengths to pixels*

Author
: David Young

Date Created
: April 15, 2021
"""

from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
from os.path import expanduser
from fundamentals import tools
from builtins import object
import sys
import os

os.environ['TERM'] = 'vt100'


def dispersion_map_to_pixel_arrays(
        log,
        dispersionMapPath,
        orderPixelTable,
        removeOffDetectorLocation=True):
    """*Use a dispersion solution to convert wavelength, slit-position and echelle order numbers to X,Y pixel positions.* 

    Return a line-list with x,y fits given a first guess dispersion map.*

    **Key Arguments:**

    - `log` -- logger
    - `dispersionMapPath` -- path to the dispersion map
    - `orderPixelTable` -- a data-frame including 'order', 'wavelength' and 'slit_pos' columns
    - `removeOffDetectorLocation` -- if data points are found to lie off the detector plane then remove them from the resutls. Default *True*

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
    log.debug('starting the ``dispersion_map_to_pixel_arrays`` function')

    from astropy.table import Table
    import math

    # READ THE FILE
    home = expanduser("~")
    dispersion_map = dispersionMapPath.replace("~", home)

    # SPEC FORMAT TO PANDAS DATAFRAME
    dat = Table.read(dispersion_map, format='fits')
    tableData = dat.to_pandas()

    # READ IN THE X- AND Y- GENERATING POLYNOMIALS FROM DISPERSION MAP FILE
    coeff = {}
    poly = {}
    check = 1

    for index, row in tableData.iterrows():
        axis = row["axis"].decode("utf-8")
        orderDeg = int(row["order_deg"])
        wavelengthDeg = int(row["wavelength_deg"])
        slitDeg = int(row["slit_deg"])

        # print(axis, orderDeg, wavelengthDeg, slitDeg)

        if check:
            for i in range(0, orderDeg + 1):
                orderPixelTable[f"order_pow_{axis}_{i}"] = orderPixelTable["order"].pow(i)
            for j in range(0, wavelengthDeg + 1):
                orderPixelTable[f"wavelength_pow_{axis}_{j}"] = orderPixelTable["wavelength"].pow(j)
            for k in range(0, slitDeg + 1):
                orderPixelTable[f"slit_position_pow_{axis}_{k}"] = orderPixelTable["slit_position"].pow(k)
            # check = 0

        coeff[axis] = [float(v) for k, v in row.items() if k not in [
            "axis", "order_deg", "wavelength_deg", "slit_deg"] and not math.isnan(v)]
        poly[axis] = chebyshev_order_wavelength_polynomials(
            log=log, orderDeg=orderDeg, wavelengthDeg=wavelengthDeg, slitDeg=slitDeg, exponentsIncluded=True, axis=axis).poly

    # CONVERT THE ORDER-SORTED WAVELENGTH ARRAYS INTO ARRAYS OF PIXEL TUPLES
    orderPixelTable["fit_x"] = poly['x'](orderPixelTable, *coeff['x'])
    orderPixelTable["fit_y"] = poly['y'](orderPixelTable, *coeff['y'])

    if removeOffDetectorLocation:
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable["fit_x"] > 0) & (orderPixelTable["fit_y"] > 0)
        orderPixelTable = orderPixelTable.loc[mask]

    log.debug('completed the ``dispersion_map_to_pixel_arrays`` function')
    return orderPixelTable


def get_cached_coeffs(
        log,
        arm,
        settings,
        recipeName,
        orderDeg,
        wavelengthDeg,
        slitDeg,
        reset=False):
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
    log.debug('starting the ``get_cached_coeffs`` function')

    from astropy.table import Table
    import math
    import numpy as np

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
        dat = Table.read(dispersion_map, format='fits')
        tableData = dat.to_pandas()

        # READ IN THE X- AND Y- COEFF FROM DISPERSION MAP FILE
        for index, row in tableData.iterrows():
            axis = row["axis"].decode("utf-8")
            orderDeg = int(row["order_deg"])
            wavelengthDeg = int(row["wavelength_deg"])
            slitDeg = int(row["slit_deg"])
            coeff[axis] = [float(v) for k, v in row.items() if k not in [
                "axis", "order_deg", "wavelength_deg", "slit_deg"] and not math.isnan(v)]
    else:
        if isinstance(orderDeg, list):
            coeff['x'] = np.ones((orderDeg[0] + 1) *
                                 (wavelengthDeg[0] + 1) * (slitDeg[0] + 1))
            coeff['y'] = np.ones((orderDeg[1] + 1) *
                                 (wavelengthDeg[1] + 1) * (slitDeg[1] + 1))
        else:
            coeff['x'] = np.ones((orderDeg + 1) *
                                 (wavelengthDeg + 1) * (slitDeg + 1))
            coeff['y'] = np.ones((orderDeg + 1) *
                                 (wavelengthDeg + 1) * (slitDeg + 1))

    log.debug('completed the ``get_cached_coeffs`` function')
    return coeff['x'], coeff['y']
