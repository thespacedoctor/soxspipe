#!/usr/bin/env python
# encoding: utf-8
"""
*Tools for working with 'set-of-files' (sof) files*

:Author:
    David Young & Marco Landoni

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
from ccdproc import ImageFileCollection
import codecs
from soxspipe.commonutils.keyword_lookup import keyword_lookup


class set_of_files(object):
    """
    *The worker class for the sof module used to homogenise various frame input formats (sof file, directory of fits fits, list of fits file paths) into a CCDProc ImageFileCollection*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``inputFrames`` -- can be a directory, a set-of-files (SOF) file or a list of fits frame paths. Default []
        - ``keys`` -- keys to report in the ImageFileCollection. Default ['cdelt1', 'cdelt2', 'eso det chip1 pszx',
            'eso dpr type', 'eso seq arm', 'exptime', 'naxis1', 'naxis2']

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a sof object, use the following:

    ```python
    # inputFrames = "/path/to/a/directory"
    # inputFrames = ['/path/to/one.fits','/path/to/two.fits','/path/to/three.fits']
    inputFrames = '/path/to/myfiles.sof'
    from soxspipe.commonutils import set_of_files
    sof = set_of_files(
        log=log,
        settings=settings,
        inputFrames=inputFrames
    )
    ```

    `inputFrames` can be a directory, a list of fits filepaths or a set-of-files (SOF) file
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[],
            keys=['MJDOBS', 'CDELT1', 'CDELT2', 'PSZX',
                  'DPR_TYPE', 'SEQ_ARM', 'EXPTIME', 'NAXIS1', 'NAXIS2']
    ):
        self.log = log
        log.debug("instansiating a new 'sof' object")
        self.settings = settings
        self.inputFrames = inputFrames

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        keys = kw(keys)
        self.keys = []
        self.keys[:] = [k.lower() for k in keys]

        # xt-self-arg-tmpx

        # Initial Actions
        # FIX RELATIVE HOME PATHS
        from os.path import expanduser
        home = expanduser("~")
        if isinstance(self.inputFrames, str) and self.inputFrames[0] == "~":
            self.inputFrames = home + "/" + self.inputFrames[1:]

        return None

    def _generate_sof_file_from_directory(
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
        from soxspipe.commonutils import set_of_files
        sof = set_of_files(
            log=log,
            settings=settings
        )
        sofFile = sof._generate_sof_file_from_directory(
            directory="path/to/directory", sofPath="/path/to/myFile.sof")
        ```
        """
        self.log.debug(
            'starting the ``_generate_sof_file_from_directory`` method')

        from soxspipe.commonutils import keyword_lookup
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get

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
                    dpr_type = hdr[kw("DPR_TYPE")].strip()
                    # CHECK ARM
                    arm = hdr[kw("SEQ_ARM")]
                    # CHECK BINNING
                    if kw('CDELT1') in hdr:
                        xbin = str(int(hdr[kw('CDELT1')]))
                        ybin = str(int(hdr[kw('CDELT2')]))
                    catagory = dpr_type + "_" + arm.strip()
                    if kw('CDELT1') in hdr:
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
            'completed the ``_generate_sof_file_from_directory`` method')
        return sofPath

    def get(
            self):
        """*return the set-of-files as a CCDProc [ImageFileCollection](https://ccdproc.readthedocs.io/en/latest/api/ccdproc.ImageFileCollection.html#ccdproc.ImageFileCollection)*

        **Return:**
            - ``sof`` -- a ccdproc ImageFileCollection of the frames

        **Usage:**

        To generate a ImageFileCollection from a directory, a list of fits filepaths or a set-of-files (SOF) file try the following:

        ```python
        # inputFrames = "/path/to/a/directory"
        # inputFrames = ['/path/to/one.fits','/path/to/two.fits','/path/to/three.fits']
        inputFrames = '/path/to/myfiles.sof'
        from soxspipe.commonutils import set_of_files
        sof = set_of_files(
            log=log,
            settings=settings,
            inputFrames=inputFrames
        )
        sofFile = sof.get()
        print(sofFile.summary)
        ```

        `inputFrames` can be a directory, a list of fits filepaths or a set-of-files (SOF) file.
        """
        self.log.debug('starting the ``get`` method')

        from os.path import expanduser
        home = expanduser("~")

        # DIRECTORY OF FRAMES
        if isinstance(self.inputFrames, str) and os.path.isdir(self.inputFrames):
            sof = ImageFileCollection(self.inputFrames, keywords=self.keys)
        elif isinstance(self.inputFrames, str) and os.path.isfile(self.inputFrames) and '.sof' in self.inputFrames:
            readFile = codecs.open(
                self.inputFrames, encoding='utf-8', mode='r')
            thisData = readFile.read()
            readFile.close()
            lines = thisData.split("\n")
            fitsFiles = []
            fitsFiles[:] = [l.split(".fits")[0].replace("~/", home + "/") +
                            ".fits" for l in lines if ".fits" in l]
            sof = ImageFileCollection(filenames=fitsFiles, keywords=self.keys)
        elif isinstance(self.inputFrames, list):
            sof = ImageFileCollection(
                filenames=self.inputFrames, keywords=self.keys)
        else:
            raise TypeError(
                "'inputFrames' should be the path to a directory of files, an SOF file or a list of FITS frame paths")

        self.log.debug('completed the ``get`` method')
        return sof

    # use the tab-trigger below for new method
    # xt-class-method
