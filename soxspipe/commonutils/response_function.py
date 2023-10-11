#!/usr/bin/env python
# encoding: utf-8
"""
*Given a standard star extracted spectrum, generate the instrument response function needed to flux calibrate science spectra*

:Author:
    David Young

:Date Created:
    July 28, 2023
"""
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class response_function(object):
    """
    *The worker class for the response_function module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``stdExtractionPath`` -- fits binary table containing the extracted standard spectrum

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a response_function object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``response_function`` to documentation
        - create a blog post about what ``response_function`` does
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
            stdExtractionPath,
            settings=False,
    ):
        self.log = log
        log.debug("instansiating a new 'response_function' object")
        self.settings = settings
        self.stdExtractionPath = stdExtractionPath
        from astropy.table import Table

        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # Initial Actions
        # OPEN EXTRACTED SPECTRUM
        # SPEC FORMAT TO PANDAS DATAFRAME
        stdExtractionDF = Table.read(self.stdExtractionPath, format='fits')
        stdExtractionDF = stdExtractionDF.to_pandas()

        from tabulate import tabulate
        self.log.print(tabulate(stdExtractionDF.head(100), headers='keys', tablefmt='psql'))

        return None

    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def get(self):
        """
        *get the response_function object*

        **Return:**
            - ``response_function``

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

        response_function = None

        self.log.debug('completed the ``get`` method')
        return response_function

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
