#!/usr/bin/env python
# encoding: utf-8
"""
*use a first-guess dispersion map to convert wavelengths to pixels*

:Author:
    David Young

:Date Created:
    April 15, 2021
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
    """*use a first-guess dispersion map to append x,y fits to line-list data frame.* 

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

        if check:
            for i in range(0, orderDeg + 1):
                orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)
            for j in range(0, wavelengthDeg + 1):
                orderPixelTable[f"wavelength_pow_{j}"] = orderPixelTable["wavelength"].pow(j)
            for k in range(0, slitDeg + 1):
                orderPixelTable[f"slit_position_pow_{k}"] = orderPixelTable["slit_position"].pow(k)
            check = 0

        coeff[axis] = [float(v) for k, v in row.items() if k not in [
            "axis", "order_deg", "wavelength_deg", "slit_deg"]]
        poly[axis] = chebyshev_order_wavelength_polynomials(
            log=log, orderDeg=orderDeg, wavelengthDeg=wavelengthDeg, slitDeg=slitDeg, exponentsIncluded=True).poly

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
