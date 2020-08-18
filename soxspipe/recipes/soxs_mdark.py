#!/usr/bin/env python
# encoding: utf-8
"""
*The recipe to generate a master dark frame*

:Author:
    David Young & Marco Landoni

:Date Created:
    January 27, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils import set_of_files
from ._base_recipe_ import _base_recipe_
import numpy as np
from astropy.nddata import CCDData
import ccdproc
from soxspipe.commonutils import keyword_lookup


class soxs_mdark(_base_recipe_):
    """
    *The soxs_mdark recipe*

    **Key Arguments**

        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths. Default []

    **Usage**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    See `produce_product` method for usage.

    ```eval_rst
    .. todo::

        - add a tutorial about ``soxs_mdark`` to documentation

    ```
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[]

    ):
        # INHERIT INITIALISATION FROM  _base_recipe_
        super(soxs_mdark, self).__init__(log=log, settings=settings)
        self.log = log
        log.debug("instansiating a new 'soxs_mdark' object")
        self.settings = settings
        self.inputFrames = inputFrames
        # xt-self-arg-tmpx

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames
        )
        self.inputFrames = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MDARK - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.verify_input_frames()
        print(self.inputFrames.summary)

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])

        return None

    def verify_input_frames(
            self):
        """*verify the input frame match those required by the soxs_mdark recipe*

        **Return:**
            - ``None``

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception shall be raised.
        """
        self.log.debug('starting the ``verify_input_frames`` method')

        self._verify_input_frames_basics()

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get

        imageTypes = self.inputFrames.values(
            keyword=kw("DPR_TYPE").lower(), unique=True)
        # MIXED INPUT IMAGE TYPES ARE BAD
        if len(imageTypes) > 1:
            imageTypes = " and ".join(imageTypes)
            print(self.inputFrames.summary)
            raise TypeError(
                "Input frames are a mix of %(imageTypes)s" % locals())

        self.log.debug('completed the ``verify_input_frames`` method')
        return None

    def produce_product(
            self):
        """*The code to generate the product of the soxs_mdark recipe*

        **Return:**
            - ``productPath`` -- the path to the final product

        **Usage**

        ```python
        from soxspipe.recipes import soxs_mdark
        recipe = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        mdarkFrame = recipe.produce_product()
        ```
        """
        self.log.debug('starting the ``produce_product`` method')

        # IMAGECOLLECTION FILEPATHS
        filepaths = self.inputFrames.files_filtered(include_path=True)

        productPath = None

        self.clean_up()

        self.log.debug('completed the ``produce_product`` method')
        return productPath

    # use the tab-trigger below for new method
    # xt-class-method

    # Override Method Attributes
    # method-override-tmpx
