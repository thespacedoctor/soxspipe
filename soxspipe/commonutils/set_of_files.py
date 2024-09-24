#!/usr/bin/env python
# encoding: utf-8
"""
*Tools for working with 'set-of-files' (sof) files*

Author
: David Young & Marco Landoni

Date Created
: January 22, 2020
"""
from os import listdir, path
from ccdproc import ImageFileCollection
import os
import sys
from builtins import object
from fundamentals import tools
from soxspipe.commonutils.keyword_lookup import keyword_lookup


class ImageFileCollection(ImageFileCollection):
    def _dict_from_fits_header(self, file_name, input_summary=None,
                               missing_marker=None):
        """

        """
        from astropy.io import fits
        from collections import OrderedDict

        def _add_val_to_dict(key, value, tbl_dict, n_previous, missing_marker):
            try:
                tbl_dict[key].append(value)
            except KeyError:
                tbl_dict[key] = [missing_marker] * n_previous
                tbl_dict[key].append(value)

        if input_summary is None:
            summary = OrderedDict()
            n_previous = 0
        else:
            summary = input_summary
            n_previous = len(summary['file'])

        try:
            h = fits.getheader(file_name, self.ext)
        except:
            h = fits.getheader(file_name, 0)

        assert 'file' not in h

        if self.location:
            # We have a location and can reconstruct the path using it
            name_for_file_column = path.basename(file_name)
        else:
            # No location, so use whatever path the user passed in
            name_for_file_column = file_name

        # Try opening header before this so that file name is only added if
        # file is valid FITS
        try:
            summary['file'].append(name_for_file_column)
        except KeyError:
            summary['file'] = [name_for_file_column]

        missing_in_this_file = [k for k in summary if (k not in h and
                                                       k != 'file')]

        multi_entry_keys = {'comment': [],
                            'history': []}

        alreadyencountered = set()
        for k, v in h.items():
            if k == '':
                continue

            if k in ['comment', 'history']:
                multi_entry_keys[k].append(str(v))
                # Accumulate these in a separate dictionary until the
                # end to avoid adding multiple entries to summary.
                continue
            elif k in alreadyencountered:
                # The "normal" multi-entries HISTORY, COMMENT and BLANK are
                # already processed so any further duplication is probably
                # a mistake. It would lead to problems in ImageFileCollection
                # to add it as well, so simply ignore those.
                import warnings
                warnings.warn(
                    'Header from file "{f}" contains multiple entries for '
                    '"{k}", the pair "{k}={v}" will be ignored.'
                    ''.format(k=k, v=v, f=file_name),
                    UserWarning)
                continue
            else:
                # Add the key to the already encountered keys so we don't add
                # it more than once.
                alreadyencountered.add(k)

            _add_val_to_dict(k, v, summary, n_previous, missing_marker)

        for k, v in multi_entry_keys.items():
            if v:
                joined = ','.join(v)
                _add_val_to_dict(k, joined, summary, n_previous,
                                 missing_marker)

        for missing in missing_in_this_file:
            summary[missing].append(missing_marker)

        return summary

    def _set_column_name_case_to_match_keywords(self, header_keys, summary_table):
        for k in header_keys:
            k_lower = k.lower()
            if k_lower != k:
                try:
                    summary_table.rename_column(k_lower, k)
                except KeyError:
                    pass


os.environ['TERM'] = 'vt100'


class set_of_files(object):
    """
    *The worker class for the sof module used to homogenize various frame input formats (sof file, directory of fits fits, list of fits file paths) into a CCDProc ImageFileCollection*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- can be a directory, a set-of-files (SOF) file or a list of fits frame paths. Default []
    - ``verbose`` -- verbose. True or False. Default *True*
    - ``recipeName`` -- the name of the recipe. Default *False*
    - ``ext`` -- the data extension for the frame. Default 0.

    **Usage**

    To initiate a sof object, use the following:

    ```python
    # inputFrames = "/path/to/a/directory"
    # inputFrames = ['/path/to/one.fits','/path/to/two.fits','/path/to/three.fits']
    inputFrames = '/path/to/myfiles.sof'
    from soxspipe.commonutils.set_of_files import set_of_files
    sof = set_of_files(
        log=log,
        settings=settings,
        inputFrames=inputFrames,
        ext=0
    )
    ```

    `inputFrames` can be a directory, a list of fits filepaths or a set-of-files (SOF) file
    """
    # Initialization

    def __init__(
            self,
            log,
            settings=False,
            inputFrames=[],
            verbose=True,
            recipeName=False,
            ext=0,
            session=None
    ):
        self.log = log
        log.debug("instantiating a new 'sof' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = recipeName
        self.ext = ext

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get

        if self.verbose:
            keys = self.settings['summary-keys']['verbose']
        else:
            keys = self.settings['summary-keys']['default']

        if recipeName and recipeName == "soxs-nod":
            keys += self.settings['summary-keys']['nodding_extras']

        keys = kw(keys)
        self.keys = []
        self.keys[:] = [k for k in keys]
        self.keys.append("file")
        # Initial Actions
        # FIX RELATIVE HOME PATHS
        from os.path import expanduser
        home = expanduser("~")
        if isinstance(self.inputFrames, str) and self.inputFrames[0] == "~":
            self.inputFrames = home + "/" + self.inputFrames[1:]

        # GRAB THE WORKSPACE SESSION
        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=self.log,
            rootDir="."
        )
        self.currentSession, allSessions = do.session_list(silent=True)

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

        **Usage**

        ```python
        from soxspipe.commonutils.set_of_files import set_of_files
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

        from astropy.io import fits

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
                        catagory += "_" + \
                            xbin.strip() + "x" + ybin.strip()

                    content += "%(fitsPath)s %(catagory)s\n" % locals()

        # RECURSIVELY CREATE MISSING DIRECTORIES
        moduleDirectory = os.path.dirname(sofPath)
        if not os.path.exists(moduleDirectory):
            os.makedirs(moduleDirectory)

        # WRITE TO FILE
        with open(sofPath, 'w') as myFile:
            myFile.write(content)

        self.log.debug(
            'completed the ``_generate_sof_file_from_directory`` method')
        return sofPath

    def get(self):
        """*return the set-of-files as a CCDProc ImageFileCollection*

        **Return:**

        - ``sof`` -- a ccdproc ImageFileCollection of the frames

        **Usage**

        To generate a ImageFileCollection from a directory, a list of fits filepaths or a set-of-files (SOF) file try the following:

        ```python
        # inputFrames = "/path/to/a/directory"
        # inputFrames = ['/path/to/one.fits','/path/to/two.fits','/path/to/three.fits']
        inputFrames = '/path/to/myfiles.sof'
        from soxspipe.commonutils.set_of_files import set_of_files
        sof = set_of_files(
            log=log,
            settings=settings,
            inputFrames=inputFrames
        )
        sofFile, supplementarySof = sof.get()
        print(sofFile.summary)
        ```

        `inputFrames` can be a directory, a list of fits filepaths or a set-of-files (SOF) file.
        """
        self.log.debug('starting the ``get`` method')

        from astropy.table import join
        import codecs
        from os.path import expanduser
        home = expanduser("~")

        if isinstance(self.inputFrames, str) and self.inputFrames[0] == "~":
            self.inputFrames = home + "/" + self.inputFrames[1:]

        # DIRECTORY OF FRAMES
        if isinstance(self.inputFrames, str) and os.path.isdir(self.inputFrames):
            if self.ext > 0:
                sofSeed = ImageFileCollection(location=self.inputFrames, ext=self.ext)
                foundKeys = [
                    k for k in self.keys if (k.lower() in sofSeed.summary.colnames or k in sofSeed.summary.colnames)]
                sof = ImageFileCollection(
                    keywords=foundKeys, location=self.inputFrames, ext=self.ext)
                missingKeys = [
                    k for k in self.keys if (k.lower() not in sofSeed.summary.colnames and k not in sofSeed.summary.colnames)]
                if len(missingKeys):
                    primExt = ImageFileCollection(
                        keywords=missingKeys, location=self.inputFrames, ext=0)
                    sof._summary = join(
                        primExt._summary, sof._summary, keys="file")
            else:
                sof = ImageFileCollection(location=self.inputFrames, keywords=self.keys, ext=self.ext)

            supplementaryFilepaths = []
            for d in os.listdir(self.inputFrames):
                filepath = os.path.join(self.inputFrames, d)
                if os.path.isfile(filepath) and ".fits" not in d.lower() and d[0] != ".":
                    supplementaryFilepaths.append(filepath)

        elif isinstance(self.inputFrames, str) and os.path.isfile(self.inputFrames) and '.sof' in self.inputFrames:

            readFile = codecs.open(
                self.inputFrames, encoding='utf-8', mode='r')
            thisData = readFile.read()
            readFile.close()
            lines = thisData.split("\n")

            # REMOVE COMMENTED LINES
            lines = [l for l in lines if len(l) and l[0] != "#"]

            fitsFiles = []
            fitsFiles[:] = [l.split(".fits")[0].replace("~/", home + "/") +
                            ".fits" for l in lines if ".fits" in l]

            supplementaryFilepaths = [
                l.replace("~/", home + "/") for l in lines if ".fits" not in l.lower() and len(l) > 3]

            # PREPEND SESSION PATHS
            if self.currentSession:
                fitsFiles[:] = [f.replace("./product", f"./sessions/{self.currentSession}/product") for f in fitsFiles]
                supplementaryFilepaths[:] = [f.replace("./product", f"./sessions/{self.currentSession}/product") for f in supplementaryFilepaths]

            # MAKE SURE FILES EXIST
            allFiles = fitsFiles.extend(supplementaryFilepaths)
            for f in fitsFiles + supplementaryFilepaths:
                exists = os.path.exists(f)
                if not exists:
                    raise FileNotFoundError(f"the input file `{f}` does not appear to exist")

            locations = [os.path.dirname(f) for f in fitsFiles]
            if len(set(locations)) == 1:
                location = locations[0]
                fitsFiles = [os.path.basename(
                    f) for f in fitsFiles]
            else:
                location = None

            if self.ext > 0:
                sofSeed = ImageFileCollection(
                    filenames=fitsFiles, location=location, ext=self.ext)
                foundKeys = [
                    k for k in self.keys if (k.lower() in sofSeed.summary.colnames or k in sofSeed.summary.colnames)]
                sof = ImageFileCollection(
                    filenames=fitsFiles, keywords=foundKeys, location=location, ext=self.ext)
                missingKeys = [
                    k for k in self.keys if (k.lower() not in sofSeed.summary.colnames and k not in sofSeed.summary.colnames)]
                if len(missingKeys):
                    primExt = ImageFileCollection(
                        filenames=fitsFiles, keywords=missingKeys, location=location, ext=0)
                    sof._summary = join(
                        primExt._summary, sof._summary, keys="file")
            else:
                sof = ImageFileCollection(
                    filenames=fitsFiles, keywords=self.keys, location=location, ext=self.ext)

        elif isinstance(self.inputFrames, list):
            fitsFiles = [f for f in self.inputFrames if ".fits" in f.lower()]
            # FIND UNIQUE FILE LOCATIONS
            locations = [os.path.dirname(f) for f in fitsFiles]
            if len(set(locations)) == 1:
                location = locations[0]
                fitsFiles = [os.path.basename(
                    f) for f in fitsFiles]
            else:
                location = None

            sof = ImageFileCollection(
                filenames=fitsFiles, keywords=self.keys, location=location, ext=self.ext)

            if self.ext > 0:
                foundKeys = [
                    k for k in self.keys if (k.lower() in sof.summary.colnames or k in sof.summary.colnames)]
                sof = ImageFileCollection(
                    filenames=fitsFiles, keywords=foundKeys, location=location, ext=self.ext)
                missingKeys = [
                    k for k in self.keys if (k.lower() not in sof.summary.colnames and k not in sof.summary.colnames)]
                if len(missingKeys):
                    primExt = ImageFileCollection(
                        filenames=fitsFiles, keywords=missingKeys, location=location, ext=0)
                    sof._summary = join(
                        primExt._summary, sof._summary, keys="file")
            fitsFiles = [os.path.basename(
                f) for f in fitsFiles]
            sof._summary["filename"] = fitsFiles
            self.keys = ["filename"] + self.keys
            supplementaryFilepaths = [
                f for f in self.inputFrames if ".fits" not in f.lower() and f[0] != "."]

        else:
            raise TypeError(
                "'inputFrames' should be the path to a directory of files, an SOF file or a list of FITS frame paths")

        supplementary_sof = self.create_supplementary_file_dictionary(
            supplementaryFilepaths)

        self.log.debug('completed the ``get`` method')
        return sof, supplementary_sof

    def create_supplementary_file_dictionary(
            self,
            supplementaryFilepaths):
        """*create supplementary file dictionary*

        **Key Arguments:**

        - ``supplementaryFilepaths`` -- the list of filepaths to generate the dictionary for

        **Return:**

        - ``supplementary_sof`` -- a dictionary of non-fits files needed for recipe
        """
        self.log.debug(
            'starting the ``create_supplementary_file_dictionary`` method')

        supplementary_sof = {}
        for f in supplementaryFilepaths:
            for a in ["NIR", "UVB", "VIS"]:
                if a.lower() in f.lower() and a not in supplementary_sof.keys():
                    supplementary_sof[a] = {}

        for f in supplementaryFilepaths:
            if "disp_map" in f.lower():
                for a in ["NIR", "UVB", "VIS"]:
                    if a.lower() in f.lower():
                        supplementary_sof[a]["DISP_MAP"] = f
            if "order_locations" in f.lower() or "order_centres" in f.lower():
                for a in ["NIR", "UVB", "VIS"]:
                    if a.lower() in f.lower():
                        supplementary_sof[a]["ORDER_LOCATIONS"] = f
            if "2d_map" in f.lower():
                for a in ["NIR", "UVB", "VIS"]:
                    if a.lower() in f.lower():
                        supplementary_sof[a]["2D_MAP"] = f

        self.log.debug(
            'completed the ``create_supplementary_file_dictionary`` method')
        return supplementary_sof

    # use the tab-trigger below for new method
    # xt-class-method
