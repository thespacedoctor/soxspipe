# set_of_files

In ESO parlance,  a set-of-files (`.sof`) file contains a list of the input data (or ingredients) needed for a pipeline recipe. This is a plain-text file where each input file path is specified with an associated classification label (one per line). Here is an example of the contents of a sof file needed by the `soxs-mdark` recipe:

```text
./raw_frames/SOXS_GEN_DARK_NIR_122_0001.fits  DARK_NIR
./raw_frames/SOXS_GEN_DARK_NIR_122_0002.fits  DARK_NIR
./raw_frames/SOXS_GEN_DARK_NIR_122_0003.fits  DARK_NIR
./raw_frames/SOXS_GEN_DARK_NIR_122_0004.fits  DARK_NIR
./raw_frames/SOXS_GEN_DARK_NIR_122_0005.fits  DARK_NIR
```

The primary purpose of the [`set_of_files`](#soxspipe.commonutils.set_of_files) utility is to read and translate these sof files into a [CCDProc ImageFileCollection](https://ccdproc.readthedocs.io/en/latest/api/ccdproc.ImageFileCollection.html), an in-memory database of the files and their header keywords. Alongside sof files, the utility also accepts a Python list of file paths or a path to a directory of FITS files as acceptable inputs.

Lines in a sof file beginning with a `#` are considered comments and, therefore, ignored by the pipeline. This is helpful to quickly remove a file from the recipe input by commenting out its line in the sof file or for adding user notes to the sof file.

### Utility API

:::{autodoc2-object} soxspipe.commonutils.set_of_files.set_of_files
:::

