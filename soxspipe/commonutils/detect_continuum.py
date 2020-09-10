#!/usr/bin/env python
# encoding: utf-8
"""
*find and fit the continuum in a pinhole flat frame with low-order polynomials. These polynominals are the central loctions of the orders.*

:Author:
    David Young

:Date Created:
    September 10, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils import detector_lookup
from os.path import expanduser
import numpy as np


class detect_continuum(object):
    """
    *The worker class for the detect_continuum module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``pinholeFlat`` -- calibrationed pinhole flat frame (CCDObject)
        - ``dispersion_map`` -- path to dispersion map csv file containing polynomial fits of the dispersion solution for the frame
        - ``settings`` -- the settings dictionary

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a detect_continuum object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``detect_continuum`` to documentation
        - create a blog post about what ``detect_continuum`` does
    ```

    ```python
    usage code
    ```

    """
    # Initialisation
    # 1. @flagged: what are the unique attrributes for each object? Add them
    # to __init__

    def __init__(
            self,
            log,
            pinholeFlat,
            dispersion_map,
            settings=False
    ):
        self.log = log
        log.debug("instansiating a new 'detect_continuum' object")
        self.settings = settings
        self.pinholeFlat = pinholeFlat
        self.dispersion_map = dispersion_map

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        self.arm = pinholeFlat.header[self.kw("SEQ_ARM")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        return None

    def get(self):
        """
        *get the detect_continuum object*

        **Return:**
            - ``detect_continuum``

        **Usage:**

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed
        ```

        ```python
        usage code
        ```
        """
        self.log.debug('starting the ``get`` method')

        orderNums, waveLengthMin, waveLengthMax = self.read_spectral_format()
        self.create_pixel_arrays(
            orderNums,
            waveLengthMin,
            waveLengthMax)

        detect_continuum = None

        self.log.debug('completed the ``get`` method')
        return detect_continuum

    def read_spectral_format(
            self):
        """*read the spectral format table to get some key parameters*

        **Return:**
            - ``orderNums`` -- a list of the order numbers
            - ``waveLengthMin`` -- a list of the maximum wavelengths reached by each order
            - ``waveLengthMax`` -- a list of the minimum wavelengths reached by each order
        """
        self.log.debug('starting the ``read_spectral_format`` method')

        kw = self.kw
        pinholeFlat = self.pinholeFlat
        dp = self.detectorParams

        # READ THE SPECTRAL FORMAT TABLE FILE
        home = expanduser("~")
        calibrationRootPath = self.settings[
            "calibration-data-root"].replace("~", home)
        spectralFormatFile = calibrationRootPath + \
            "/" + dp["spectral format table"]
        spectralFormat = np.genfromtxt(
            spectralFormatFile, delimiter=',', names=True)

        # print(spectralFormat.dtype.names)

        # EXTRACT REQUIRED PARAMETERS
        orderNums = spectralFormat["ORDER"]
        waveLengthMin = spectralFormat["WLMIN"]
        waveLengthMax = spectralFormat["WLMAX"]

        self.log.debug('completed the ``read_spectral_format`` method')
        return orderNums, waveLengthMin, waveLengthMax

    def create_pixel_arrays(
            self,
            orderNums,
            waveLengthMin,
            waveLengthMax):
        """*create a pixel array for the approximate centre of each order*

        **Key Arguments:**
            - ``orderNums`` -- a list of the order numbers
            - ``waveLengthMin`` -- a list of the maximum wavelengths reached by each order
            - ``waveLengthMax`` -- a list of the minimum wavelengths reached by each order

        **Return:**
            - None

        **Usage:**

        ```python
        usage code 
        ```

        ---

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - write a command-line tool for this method
            - update package tutorial with command-line tool info if needed
        ```
        """
        self.log.debug('starting the ``create_pixel_arrays`` method')

        # READ THE FILE
        home = expanduser("~")
        dispersion_map = self.dispersion_map.replace("~", home)
        dispersion_map = np.genfromtxt(
            dispersion_map, delimiter=',', names=True)

        print(dispersion_map.dtype.names)

        self.log.debug('completed the ``create_pixel_arrays`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
