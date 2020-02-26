#!/usr/bin/env python
# encoding: utf-8
"""
*Given a keyword token and instrument name return the exact FITS Header keyword*

:Author:
    David Young

:Date Created:
    February 26, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools


class keyword_lookup(object):
    """
    *The worker class for the keyword_lookup module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary

    **Usage:**

        To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

        To initiate a keyword_lookup object, use the following:

        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this class
            - add a tutorial about ``keyword_lookup`` to documentation
            - create a blog post about what ``keyword_lookup`` does

        .. code-block:: python 

            usage code   
    """
    # Initialisation
    # 1. @flagged: what are the unique attrributes for each object? Add them
    # to __init__

    def __init__(
            self,
            log,
            settings=False,

    ):
        self.log = log
        log.debug("instansiating a new 'keyword_lookup' object")
        self.settings = settings
        # xt-self-arg-tmpx

        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes
        self.instrument = settings["instrument"]

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # Initial Actions

        return None

    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def get(self):
        """
        *get the keyword_lookup object*

        **Return:**
            - ``keyword_lookup``

        **Usage:**
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed

        .. code-block:: python 

            usage code 
        """
        self.log.debug('starting the ``get`` method')

        keyword_lookup = None

        self.log.debug('completed the ``get`` method')
        return keyword_lookup

    def _select_dictionary(
            self):
        """*select the keyword dictionary based on the instrument passed via the settings*

        **Return:**
            - ``kwDict`` -- the python dictionary of keywords (key = tag, value = fits keyword)

        **Usage:**

        ```python
        from soxspipe.commonutils import keyword_lookup
        this = keyword_lookup(
            log=log,
            settings=settings
        )
        kw = this._select_dictionary()
        ```
        """
        self.log.debug('starting the ``_select_dictionary`` method')

        # GENERATE PATH TO YAML DICTIONARY
        yamlFilePath = os.path.dirname(os.path.dirname(
            __file__)) + "/resources/" + self.instrument + "_keywords.yaml"

        # YAML CONTENT TO DICTIONARY
        import yaml
        with open(yamlFilePath, 'r') as stream:
            kwDict = yaml.load(stream)

        self.log.debug('completed the ``_select_dictionary`` method')
        return kwDict

    # use the tab-trigger below for new method
    # xt-class-method
