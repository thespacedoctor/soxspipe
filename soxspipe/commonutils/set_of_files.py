#!/usr/bin/env python
"""
*Tools for working with 'set-of-files' (sof) files*

Author
: David Young & Marco Landoni

Date Created
: January 22, 2020
"""

import logging
import os
from os import path

from ccdproc import ImageFileCollection

from soxspipe.commonutils.keyword_lookup import keyword_lookup


class ImageFileCollection(ImageFileCollection):
    def _dict_from_fits_header(
        self, file_name, input_summary=None, missing_marker=None
    ):
        """*summarise one file's header, keeping keyword case (overrides the ccdproc method)*

        **Key Arguments:**

        - ``file_name`` -- the path of the FITS file
        - ``input_summary`` -- the summary built so far. Default *None*
        - ``missing_marker`` -- the value recorded for a keyword the file lacks. Default *None*

        **Return:**

        - ``summary`` -- the summary, with this file's row appended
        """
        from collections import OrderedDict

        from astropy.io import fits

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
            n_previous = len(summary["file"])

        try:
            h = fits.getheader(file_name, self.ext)
        except (OSError, IndexError, KeyError) as e:
            # THIS CLASS SUBCLASSES CCDPROC'S ImageFileCollection, WHICH CARRIES NO SOXSPIPE LOGGER
            logging.getLogger(__name__).debug(f"_dict_from_fits_header: `h = fits.getheade...` failed, continuing: {e}")
            h = fits.getheader(file_name, 0)

        # KEEP THE BARE ASSERT: AN EXPLICIT RAISE WOULD ALSO FIRE UNDER `python -O`, WHERE THIS CHECK IS SKIPPED TODAY.
        # IT MIRRORS THE UPSTREAM CCDPROC METHOD THIS OVERRIDES
        assert "file" not in h  # noqa: S101

        # WITH A LOCATION WE CAN RECONSTRUCT THE PATH USING IT; WITHOUT ONE, USE WHATEVER PATH THE USER PASSED IN
        name_for_file_column = path.basename(file_name) if self.location else file_name

        # TRY OPENING HEADER BEFORE THIS SO THAT FILE NAME IS ONLY ADDED IF
        # FILE IS VALID FITS
        try:
            summary["file"].append(name_for_file_column)
        except KeyError:
            summary["file"] = [name_for_file_column]

        missing_in_this_file = [k for k in summary if (k not in h and k != "file")]

        multi_entry_keys = {"comment": [], "history": []}

        alreadyencountered = set()
        for k, v in h.items():
            if k == "":
                continue

            if k in ["comment", "history"]:
                multi_entry_keys[k].append(str(v))
                # ACCUMULATE THESE IN A SEPARATE DICTIONARY UNTIL THE
                # END TO AVOID ADDING MULTIPLE ENTRIES TO SUMMARY.
                continue
            if k in alreadyencountered:
                # THE "NORMAL" MULTI-ENTRIES HISTORY, COMMENT AND BLANK ARE
                # ALREADY PROCESSED SO ANY FURTHER DUPLICATION IS PROBABLY
                # A MISTAKE. IT WOULD LEAD TO PROBLEMS IN IMAGEFILECOLLECTION
                # TO ADD IT AS WELL, SO SIMPLY IGNORE THOSE.
                import warnings

                warnings.warn(
                    f'Header from file "{file_name}" contains multiple entries for '
                    f'"{k}", the pair "{k}={v}" will be ignored.',
                    UserWarning,
                    stacklevel=1,
                )
                continue
            # ADD THE KEY TO THE ALREADY ENCOUNTERED KEYS SO WE DON'T ADD
            # IT MORE THAN ONCE.
            alreadyencountered.add(k)

            _add_val_to_dict(k, v, summary, n_previous, missing_marker)

        for k, v in multi_entry_keys.items():
            if v:
                joined = ",".join(v)
                _add_val_to_dict(k, joined, summary, n_previous, missing_marker)

        for missing in missing_in_this_file:
            summary[missing].append(missing_marker)

        return summary

    def _set_column_name_case_to_match_keywords(self, header_keys, summary_table):
        for k in header_keys:
            k_lower = k.lower()
            if k_lower != k:
                try:
                    summary_table.rename_column(k_lower, k)
                except KeyError as e:
                    # THIS CLASS SUBCLASSES CCDPROC'S ImageFileCollection, WHICH CARRIES NO SOXSPIPE LOGGER
                    logging.getLogger(__name__).debug(
                        f"_set_column_name_case_to_match_keywords: `rename_column(k_lower, k)` failed, continuing: {e}"
                    )


os.environ["TERM"] = "vt100"


def _supplementary_path_from_sof_line(line, home):
    """*return the path column from a supplementary SOF row*

    **Key Arguments:**

    - ``line`` -- one line of the SOF file
    - ``home`` -- the user's home directory, substituted for a leading ``~/``

    **Return:**

    - ``path`` -- the supplementary file path, with a trailing arm tag removed when the line is not itself an
      existing path
    """
    expandedLine = line.replace("~/", home + "/")
    if os.path.exists(expandedLine):
        return expandedLine
    columns = expandedLine.rsplit(maxsplit=1)
    if (
        len(columns) == 2
        and columns[1].isupper()
        and columns[1].replace("_", "").isalnum()
        and columns[1].endswith(("_NIR", "_UVB", "_VIS"))
    ):
        return columns[0]
    return expandedLine


def _join_fits_summaries_in_input_order(primarySummary, extensionSummary):
    """*join primary and extension headers without reordering input frames*

    **Key Arguments:**

    - ``primarySummary`` -- the summary table of primary-header keywords
    - ``extensionSummary`` -- the summary table of data-extension keywords, in input order

    **Return:**

    - ``joinedSummary`` -- the joined summary table, in the row order of ``extensionSummary``
    """
    from astropy.table import join

    inputOrderColumn = "_soxspipe_input_order"
    orderedExtension = extensionSummary.copy()
    orderedExtension[inputOrderColumn] = range(len(orderedExtension))
    joinedSummary = join(primarySummary, orderedExtension, keys="file")
    joinedSummary.sort(inputOrderColumn)
    joinedSummary.remove_column(inputOrderColumn)
    return joinedSummary


def _supplementary_files_in_directory(directory):
    """*return the non-FITS, non-hidden files in a directory of frames*

    **Key Arguments:**

    - ``directory`` -- the directory of frames

    **Return:**

    - ``supplementaryFilepaths`` -- the supplementary file paths, in directory-listing order
    """
    supplementaryFilepaths = []
    for d in os.listdir(directory):
        filepath = os.path.join(directory, d)
        if (
            os.path.isfile(filepath)
            and ".fits" not in d.lower()
            and d[0] != "."
        ):
            supplementaryFilepaths.append(filepath)
    return supplementaryFilepaths


def _common_location(fitsFiles):
    """*find the directory shared by all the frames*

    **Key Arguments:**

    - ``fitsFiles`` -- the frame paths

    **Return:**

    - ``location`` -- the shared directory, or None when the frames span several directories
    - ``fitsFiles`` -- the frame base names when a shared directory is found, otherwise the paths unchanged
    """
    locations = [os.path.dirname(f) for f in fitsFiles]
    if len(set(locations)) == 1:
        return locations[0], [os.path.basename(f) for f in fitsFiles]
    return None, fitsFiles


class set_of_files:
    """
    *The worker class for the sof module used to homogenize various frame input formats (sof file, directory of fits
      fits, list of fits file paths) into a CCDProc ImageFileCollection*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- can be a directory, a set-of-files (SOF) file or a list of fits frame paths. Default []
    - ``verbose`` -- verbose. True or False. Default *True*
    - ``recipeName`` -- the name of the recipe. Default *False*
    - ``ext`` -- the data extension for the frame. Default 0.
    - ``session`` -- unused; the workspace session is read from the data organiser. Default *None*

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

    # INITIALIZATION

    def __init__(
        self,
        log,
        settings=False,
        # THE DEFAULT LIST IS NEVER MUTATED. A None SENTINEL WOULD TURN AN EXPLICIT `inputFrames=None` FROM A
        # TypeError IN get() INTO AN EMPTY COLLECTION
        inputFrames=[],  # noqa: B006
        verbose=True,
        recipeName=False,
        ext=0,
        session=None,
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
        kw = keyword_lookup(log=self.log, settings=self.settings).get

        keys = self.settings["summary-keys"]["verbose"] if self.verbose else self.settings["summary-keys"]["default"]

        if recipeName and recipeName == "soxs-nod":
            # BUILD A NEW LIST -- NEVER MUTATE THE SETTINGS DICT'S OWN LIST (DY-73)
            keys = keys + self.settings["summary-keys"]["nodding_extras"]

        keys = kw(keys)
        self.keys = []
        self.keys[:] = list(keys)
        self.keys.append("file")
        # INITIAL ACTIONS
        # FIX RELATIVE HOME PATHS
        from os.path import expanduser

        home = expanduser("~")
        if isinstance(self.inputFrames, str) and self.inputFrames.startswith("~"):
            self.inputFrames = home + "/" + self.inputFrames[1:]

        # GRAB THE WORKSPACE SESSION
        from soxspipe.commonutils import data_organiser

        do = data_organiser(log=self.log, rootDir=".", dbConnect=False)
        self.currentSession, allSessions = do.session_list(silent=True)
        do.close()

        return

    def _generate_sof_file_from_directory(self, directory, sofPath):
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
        self.log.debug("starting the ``_generate_sof_file_from_directory`` method")

        from astropy.io import fits

        from soxspipe.commonutils import keyword_lookup

        kw = keyword_lookup(log=self.log, settings=self.settings).get

        # MAKE RELATIVE HOME PATH ABSOLUTE
        from os.path import expanduser

        home = expanduser("~")
        if directory[0] == "~":
            directory = directory.replace("~", home)
        if sofPath[0] == "~":
            sofPath = sofPath.replace("~", home)

        content = ""
        for d in sorted(os.listdir(directory)):
            if os.path.isfile(os.path.join(directory, d)) and (
                os.path.splitext(d)[-1].lower() == ".fits"
            ):
                fitsPath = os.path.abspath(os.path.join(directory, d))
                # OPEN FITS FILE AT HDULIST - HDU (HEADER DATA UNIT) CONTAINS A HEADER AND A DATA ARRAY (IMAGE) OR
                # TABLE.
                with fits.open(fitsPath) as hdul:
                    # READ HEADER INTO MEMORY
                    hdr = hdul[0].header
                    dpr_type = hdr[kw("DPR_TYPE")].strip()
                    # CHECK ARM
                    arm = hdr[kw("SEQ_ARM")]
                    # CHECK BINNING
                    if kw("CDELT1") in hdr:
                        xbin = str(int(hdr[kw("CDELT1")]))
                        ybin = str(int(hdr[kw("CDELT2")]))
                    catagory = dpr_type + "_" + arm.strip()
                    if kw("CDELT1") in hdr:
                        catagory += "_" + xbin.strip() + "x" + ybin.strip()

                    content += f"{fitsPath} {catagory}\n"

        # RECURSIVELY CREATE MISSING DIRECTORIES
        moduleDirectory = os.path.dirname(sofPath)
        if not os.path.exists(moduleDirectory):
            os.makedirs(moduleDirectory)

        # WRITE TO FILE
        with open(sofPath, "w") as myFile:
            myFile.write(content)

        self.log.debug("completed the ``_generate_sof_file_from_directory`` method")
        return sofPath

    def get(self):
        """*return the set-of-files as a CCDProc ImageFileCollection*

        **Return:**

        - ``sof`` -- a ccdproc ImageFileCollection of the frames

        **Usage**

        To generate a ImageFileCollection from a directory, a list of fits filepaths or a set-of-files (SOF) file try
        the following:

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
        self.log.debug("starting the ``get`` method")

        from os.path import expanduser

        home = expanduser("~")

        if isinstance(self.inputFrames, str) and self.inputFrames.startswith("~"):
            self.inputFrames = home + "/" + self.inputFrames[1:]

        # DIRECTORY OF FRAMES
        if isinstance(self.inputFrames, str) and os.path.isdir(self.inputFrames):
            sof = self._collection_from_frames(location=self.inputFrames)
            supplementaryFilepaths = _supplementary_files_in_directory(self.inputFrames)

        elif (
            isinstance(self.inputFrames, str)
            and os.path.isfile(self.inputFrames)
            and ".sof" in self.inputFrames
        ):
            fitsFiles, supplementaryFilepaths = self._frames_from_sof_file(home)
            location, fitsFiles = _common_location(fitsFiles)
            sof = self._collection_from_frames(location=location, filenames=fitsFiles)

        elif isinstance(self.inputFrames, list):
            fitsFiles = [f for f in self.inputFrames if ".fits" in f.lower()]
            # FIND UNIQUE FILE LOCATIONS
            location, fitsFiles = _common_location(fitsFiles)
            sof = self._collection_from_frames(location=location, filenames=fitsFiles)
            fitsFiles = [os.path.basename(f) for f in fitsFiles]
            sof._summary["filename"] = fitsFiles
            self.keys = ["filename"] + self.keys
            supplementaryFilepaths = [
                f for f in self.inputFrames if ".fits" not in f.lower() and f[0] != "."
            ]

        else:
            raise TypeError(
                "'inputFrames' should be the path to a directory of files, an SOF file or a list of FITS frame paths"
            )

        supplementary_sof = self.create_supplementary_file_dictionary(
            supplementaryFilepaths
        )

        self.log.debug("completed the ``get`` method")
        return sof, supplementary_sof

    def _frames_from_sof_file(self, home):
        """*read the FITS and supplementary file paths listed in the SOF file*

        **Key Arguments:**

        - ``home`` -- the user's home directory, substituted for a leading ``~/``

        **Return:**

        - ``fitsFiles`` -- the FITS frame paths, followed by the supplementary file paths
        - ``supplementaryFilepaths`` -- the supplementary (non-FITS) file paths
        """
        import codecs

        # KEEP codecs.open: UNLIKE open() IT DOES NOT TRANSLATE CRLF LINE ENDINGS
        with codecs.open(self.inputFrames, encoding="utf-8", mode="r") as readFile:
            thisData = readFile.read()
        lines = thisData.split("\n")

        # REMOVE COMMENTED LINES
        lines = [sofLine for sofLine in lines if len(sofLine) and sofLine[0] != "#"]

        fitsFiles = []
        fitsFiles[:] = [
            sofLine.split(".fits")[0].replace("~/", home + "/") + ".fits"
            for sofLine in lines
            if ".fits" in sofLine
        ]

        supplementaryFilepaths = [
            _supplementary_path_from_sof_line(sofLine, home)
            for sofLine in lines
            if ".fits" not in sofLine.lower() and len(sofLine) > 3
        ]

        # PREPEND SESSION PATHS
        if self.currentSession:
            fitsFiles[:] = [
                f.replace("./reduced", f"./sessions/{self.currentSession}/reduced")
                for f in fitsFiles
            ]
            supplementaryFilepaths[:] = [
                f.replace("./reduced", f"./sessions/{self.currentSession}/reduced")
                for f in supplementaryFilepaths
            ]

        # MAKE SURE FILES EXIST
        fitsFiles.extend(supplementaryFilepaths)
        for f in fitsFiles + supplementaryFilepaths:
            exists = os.path.exists(f)
            if not exists:
                raise FileNotFoundError(
                    f"the input file `{f}` does not appear to exist"
                )

        return fitsFiles, supplementaryFilepaths

    def _collection_from_frames(self, location, filenames=None):
        """*build the ImageFileCollection, filling keys missing from a data extension from the primary header*

        **Key Arguments:**

        - ``location`` -- the directory holding the frames, or None when the frames span several directories
        - ``filenames`` -- the frame file names. Default *None*, which collects every FITS file in ``location``

        **Return:**

        - ``sof`` -- a ccdproc ImageFileCollection of the frames
        """
        if self.ext > 0:
            sofSeed = ImageFileCollection(
                filenames=filenames, location=location, ext=self.ext
            )
            foundKeys = [
                k
                for k in self.keys
                if (
                    k.lower() in sofSeed.summary.colnames
                    or k in sofSeed.summary.colnames
                )
            ]
            sof = ImageFileCollection(
                filenames=filenames,
                keywords=foundKeys,
                location=location,
                ext=self.ext,
            )
            missingKeys = [
                k
                for k in self.keys
                if (
                    k.lower() not in sofSeed.summary.colnames
                    and k not in sofSeed.summary.colnames
                )
            ]
            if len(missingKeys):
                primExt = ImageFileCollection(
                    filenames=filenames,
                    keywords=missingKeys,
                    location=location,
                    ext=0,
                )
                sof._summary = _join_fits_summaries_in_input_order(
                    primExt._summary,
                    sof._summary,
                )
        else:
            sof = ImageFileCollection(
                filenames=filenames,
                keywords=self.keys,
                location=location,
                ext=self.ext,
            )
        return sof

    def create_supplementary_file_dictionary(self, supplementaryFilepaths):
        """*create supplementary file dictionary*

        **Key Arguments:**

        - ``supplementaryFilepaths`` -- the list of filepaths to generate the dictionary for

        **Return:**

        - ``supplementary_sof`` -- a dictionary of non-fits files needed for recipe
        """
        self.log.debug("starting the ``create_supplementary_file_dictionary`` method")

        supplementary_sof = {}
        for f in supplementaryFilepaths:
            for a in ["NIR", "UVB", "VIS"]:
                if a.lower() in f.lower() and a not in supplementary_sof:
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

        self.log.debug("completed the ``create_supplementary_file_dictionary`` method")
        return supplementary_sof

    # USE THE TAB-TRIGGER BELOW FOR NEW METHOD
    # xt-class-method
