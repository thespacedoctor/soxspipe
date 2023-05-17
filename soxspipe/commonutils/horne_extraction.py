#!/usr/bin/env python
# encoding: utf-8
"""
*perform optimal source extraction using the Horne method (Horne 1986)*

:Author:
    Marco Landoni & David Young

:Date Created:
    May 17, 2023
"""
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class horne_extraction(object):
    """
    *The worker class for the horne_extraction module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``skyModelFrame`` -- sky model frame
        - ``skySubFrame`` -- sky subtracted frame
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
        - ``twoDMapPath`` -- 2D dispersion map image path
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a horne_extraction object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``horne_extraction`` to documentation
        - create a blog post about what ``horne_extraction`` does
    ```

    ```python
    usage code 
    ```

    """

    def __init__(
            self,
            log,
            settings,
            skyModelFrame,
            skySubtractedFrame,
            twoDMapPath,
            recipeName=False,
            qcTable=False,
            productsTable=False,

    ):
        self.log = log
        log.debug("instansiating a new 'horne_extraction' object")
        self.settings = settings
        # xt-self-arg-tmpx

        return None

    def extract(self):
        """
        *get the horne_extraction object*

        **Return:**
            - ``horne_extraction``

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
        self.log.debug('starting the ``extract`` method')

        self.log.debug('completed the ``extract`` method')
        return None

    # xt-class-method
