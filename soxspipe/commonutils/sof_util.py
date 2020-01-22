#!/usr/bin/env python
# encoding: utf-8
"""
*Tools for working with 'set-of-files' (sof) files*

:Author:
    David Young

:Date Created:
    January 22, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from astropy.io import fits


class sof_util(object):
    """
    *The worker class for the sof module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a sof object, use the following:

    ```python
    usage code 
    ```

    ---

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``sof`` to documentation
        - create a blog post about what ``sof`` does
    ```
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,

    ):
        self.log = log
        log.debug("instansiating a new 'sof' object")
        self.settings = settings
        # xt-self-arg-tmpx

        # Initial Actions

        return None

    def generate_sof_file_from_directory(
            self,
            directory,
            sofPath):
        """*generate an sof file from a directory of FITS frames*

        **Key Arguments:**
            - ``directory`` -- the path to the directory to containing the FITS files.
            - ``sofPath`` -- the path to generate the sof file to

        **Return:**
            - ``sofPath`` -- the path to the sof file

        **Usage:**

        ```python
        from soxspipe.commonutils import sof_util
        sof = sof_util(
            log=log,
            settings=settings
        )
        sofFile = sof.generate_sof_file_from_directory(
            directory="path/to/directory", sofPath="/path/to/myFile.sof")
        ```

        ---

        ```eval_rst
        ..  todo::

            - write a command-line tool for this method
        ```

        """
        self.log.debug(
            'starting the ``generate_sof_file_from_directory`` method')

        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser
        home = expanduser("~")
        if directory[0] == "~":
            directory = directory.replace("~", home)
        if sofPath[0] == "~":
            sofPath = sofPath.replace("~", home)

        content = ""
        for d in sorted(os.listdir(directory)):
            if os.path.isfile(os.path.join(directory, d)) and (os.path.splitext(d)[-1].lower() == ".fits"):
                fitsPath = os.path.abspath(os.path.join(directory, d))
                # OPEN FITS FILE AT HDULIST - HDU (HEADER DATA UNIT) CONTAINS A HEADER AND A DATA ARRAY (IMAGE) OR
                # TABLE.
                with fits.open(fitsPath) as hdul:
                    # READ HEADER INTO MEMORY
                    hdr = hdul[0].header
                    # PRINT FULL FITS HEADER TO STDOUT
                    # print(repr(hdr).strip())
                    dpr_type = hdr['HIERARCH ESO DPR TYPE'].strip()
                    # CHECK ARM
                    arm = hdr['HIERARCH ESO SEQ ARM']
                    # CHECK BINNING
                    if 'CDELT1' in hdr:
                        xbin = str(int(hdr['CDELT1']))
                        ybin = str(int(hdr['CDELT2']))
                    catagory = dpr_type + "_" + arm.strip()
                    if 'CDELT1' in hdr:
                        catagory  += "_" + \
                            xbin.strip() + "x" + ybin.strip()

                    content += "%(fitsPath)s %(catagory)s\n" % locals()

        # Recursively create missing directories
        moduleDirectory = os.path.dirname(sofPath)
        if not os.path.exists(moduleDirectory):
            os.makedirs(moduleDirectory)

        # WRITE TO FILE
        with open(sofPath, 'w') as myFile:
            myFile.write(content)

        self.log.debug(
            'completed the ``generate_sof_file_from_directory`` method')
        return sofPath

    # use the tab-trigger below for new method
    # xt-class-method
