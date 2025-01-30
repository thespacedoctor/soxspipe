
The [`soxs_stare`](../../../recipes/soxs_stare.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] stare <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
```

To rerun a previously executed `soxs_stare` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_stare` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe stare sof/2019.10.31T02.29.50.606_UVB_2X2_SLOW_STARE_1200.0S_XSHOOTER_NGC___985.sof -s ./sessions/base/soxspipe.yaml   -x
```

To adjust the default settings for the `soxs_stare` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_stare` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_stare_parameters`.


Product files are written in the `products/soxs_stare`, and QC plots are in the `qc/soxs_stare` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_stare` are found in {numref}`soxs_stare_qc`.


### Reduction Tips


Typically, this recipe fails because the object trace is not found or fitted correctly. This may simply be because the object is too faint in the data, but a few settings can be adjusted to try and secure a fit. First, increase `slice-length` by a few pixels and try the recipe again. If this does not work, adjust `poly-fitting-residual-clipping-sigma` to both lower and higher levels (by 1 or 2 sigma each time) and rerun the recipe each time.
