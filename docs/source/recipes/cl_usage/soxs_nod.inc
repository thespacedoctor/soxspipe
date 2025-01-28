The [`soxs_nod`](../../../recipes/soxs_nod.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] nod <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
```

To rerun a previously executed `soxs_nod` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_nod` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe nod sof/2019.10.30T23.34.48.0961_NIR_NOD_200.0S_XSHOOTER_STD_FLUX.sof -s ./sessions/base/soxspipe.yaml  -x
```

To adjust the default settings for the `soxs_nod` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_nod` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_nod_parameters`.


Product files are written in the `products/soxs_nod`, and QC plots are in the `qc/soxs_nod` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_nod` are found in {numref}`soxs_nod_qc`.


### Reduction Tips


