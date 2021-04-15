#!/usr/bin/env python
# encoding: utf-8
"""
*use a first-guess dispersion map to convert wavelengths to pixels*

:Author:
    David Young

:Date Created:
    April 15, 2021
"""
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from os.path import expanduser
import unicodecsv as csv
from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
import numpy as np


def dispersion_map_to_pixel_arrays(
        log,
        dispersionMapPath,
        orderWavelengthDict):
    """*use a first-guess dispersion map to convert wavelengths to pixels.* 

    Return a dictionary of {order-number:pixel-array} given a first guess dispersion map and a list of order-based wavelengths.*

    **Key Arguments:**

    - `log` -- logger
    - `dispersionMapPath` -- path to the dispersion map
    - `orderWavelengthDict` -- a dictionary with order number as keys and wavelength arrays for those orders are values; `{order-number:pixel-array}`

    **Usage:**

    ```python
    orderWavelengthDict = {11: [850.3, 894.3, 983.2]}
    pixelArrays = dispersion_map_to_pixel_arrays(
        log=log,
        dispersionMapPath="/path/to/map.csv",
        orderWavelengthDict=orderWavelengthDict
    )
    ```           
    """
    log.debug('starting the ``dispersion_map_to_pixel_arrays`` function')

    # READ THE FILE
    home = expanduser("~")
    dispersion_map = dispersionMapPath.replace("~", home)

    # READ IN THE X- AND Y- GENERATING POLYNOMIALS FROM DISPERSION MAP FILE
    coeff = {}
    poly = {}
    with open(dispersion_map, 'rb') as csvFile:
        csvReader = csv.DictReader(
            csvFile, dialect='excel', delimiter=',', quotechar='"')
        for row in csvReader:
            axis = row["axis"]
            order_deg = int(row["order-deg"])
            wavelength_deg = int(row["wavelength-deg"])
            coeff[axis] = [float(v) for k, v in row.items() if k not in [
                "axis", "order-deg", "wavelength-deg"]]
            poly[axis] = chebyshev_order_wavelength_polynomials(
                log=log, order_deg=order_deg, wavelength_deg=wavelength_deg).poly
    csvFile.close()

    # CONVERT THE ORDER-SORTED WAVELENGTH ARRAYS INTO ARRAYS OF PIXEL TUPLES
    pixelArrays = {}
    for order, wlArray in orderWavelengthDict.items():
        orderArray = np.ones(len(wlArray)) * order
        order_wave = (orderArray, wlArray)
        xArray = poly['x'](order_wave, *coeff['x'])
        yArray = poly['y'](order_wave, *coeff['y'])
        pixelArrays[order] = [(x, y) for x, y in zip(
            xArray, yArray) if (x > 0 and y > 0)]

    log.debug('completed the ``dispersion_map_to_pixel_arrays`` function')
    return pixelArrays
