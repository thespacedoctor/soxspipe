#!/usr/bin/env python
# encoding: utf-8
"""
*Subtract the sky background using the Kelson Method*

:Author:
    David Young

:Date Created:
    April 14, 2022
"""
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


# OR YOU CAN REMOVE THE CLASS BELOW AND ADD A WORKER FUNCTION ... SNIPPET TRIGGER BELOW
# xt-worker-def

class subtract_sky(object):
    """
    *The worker class for the subtract_sky module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``objectFrame`` -- the image frame in need of sky subtraction
        - ``twoDMap`` -- 2D dispersion map image path
        - ``dispMap`` -- path to dispersion map. Default * False*

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a subtract_sky object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``subtract_sky`` to documentation
        - create a blog post about what ``subtract_sky`` does
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
            twoDMap,
            objectFrame,
            dispMap=False,
            settings=False,

    ):
        self.log = log
        log.debug("instansiating a new 'subtract_sky' object")
        self.settings = settings
        self.objectFrame = objectFrame
        self.twoDMap = twoDMap
        self.dispMap = dispMap

        # xt-self-arg-tmpx

        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # Initial Actions
        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
        self.mapDF = twoD_disp_map_image_to_dataframe(log=self.log, twoDMapPath=twoDMap, assosiatedFrame=self.objectFrame)

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=objectFrame, show=True, ext=False, stdWindow=5, title="Object Frame with dispersion solution grid", surfacePlot=True, dispMap=dispMap, dispMapDF=self.mapDF)

        from tabulate import tabulate
        print(tabulate(self.mapDF.head(100), headers='keys', tablefmt='psql'))

        return None

    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def get(self):
        """
        *get the subtract_sky object*

        **Return:**
            - ``subtract_sky``

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

        subtract_sky = None

        self.log.debug('completed the ``get`` method')
        return subtract_sky

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
