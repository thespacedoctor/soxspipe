The [`soxs_nod_std`](../../../recipes/soxs_nod.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] nod_std <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
```

Select `nod_std` when reducing spectroscopic standard stars taken in nodding mode. For spectroscopic standard stars, the pipeline produces a response curve that can be used to flux-calibrate the science sources.

To rerun a previously executed `nod_std` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_nod_std ` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe nod_std sof/20260111T041914_NIR_3_NOD_STD_FLUX_SLIT1_0_150_0S_SOXS.sof -s ./sessions/base/soxspipe.yaml  -x
```

To adjust the default settings for the `soxs_nod_std` and recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs-nod-std` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_nod_parameters`.


Product files are written in the `products/soxs_nod_std ` and QC plots are in the qc/soxs_nod_std` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_nod_std` are found in {numref}`soxs_nod_qc`.


### Reduction Tips


Typically, this recipe fails because the object trace is not found or fitted correctly. This may simply be because the object is too faint in the data, but a few settings can be adjusted to try and secure a fit. First, increase `slice-length` by a few pixels and try the recipe again. If this does not work, adjust `poly-fitting-residual-clipping-sigma` to both lower and higher levels (by 1 or 2 sigma each time) and rerun the recipe each time.
