#!/usr/bin/env python
# encoding: utf-8
"""
*Fit and subtract background flux from scattered light from frame*

:Author:
    David Young

:Date Created:
    June  3, 2021
"""
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools


# OR YOU CAN REMOVE THE CLASS BELOW AND ADD A WORKER FUNCTION ... SNIPPET TRIGGER BELOW
# xt-worker-def

class subtract_background(object):
    """
    *The worker class for the subtract_background module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``frame`` -- the frame to subtract background light from
        - ``orderTable`` -- the order geometry table

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a subtract_background object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``subtract_background`` to documentation
        - create a blog post about what ``subtract_background`` does
    ```

    ```python
    usage code 
    ```

    """
    # Initialisation
    # 1. @flagged: what are the unique attributes for each object? Add them
    # to __init__

    def __init__(
            self,
            log,
            frame,
            orderTable,
            settings=False
    ):
        self.log = log
        log.debug("instansiating a new 'subtract_background' object")
        self.settings = settings
        # xt-self-arg-tmpx

        # 2. @flagged: what are the default attributes each object could have? Add them to variable attribute set here
        # Variable Data Attributes

        # 3. @flagged: what variable attributes need overridden in any baseclass(es) used
        # Override Variable Data Attributes

        # Initial Actions

        return None

    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def get(self):
        """
        *get the subtract_background object*

        **Return:**
            - ``subtract_background``

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

        subtract_background = None

        self.log.debug('completed the ``get`` method')
        return subtract_background

    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need amending? amend them here
    # Override Method Attributes
    # method-override-tmpx
