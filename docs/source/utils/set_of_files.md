# `set_of_files`

The [`set_of_files`](#soxspipe.commonutils.set_of_files) utility helps to translate and homogenise various recipe input-frame lists. This allows recipes to accept any of the following inputs:

* an ESORex-like sof file,
* a directory of FITS files
* a list of fits file paths

Behind the scenes [`set_of_files`](#soxspipe.commonutils.set_of_files) converts the lists into a [CCDProc ImageFileCollection](https://ccdproc.readthedocs.io/en/latest/api/ccdproc.ImageFileCollection.html).

Lines in a sof file beginning with a `#` are considered as comments and therefore ignored by the pipeline.

:::{autodoc2-object} soxspipe.commonutils.set_of_files.set_of_files
:::
