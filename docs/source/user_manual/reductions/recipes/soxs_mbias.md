## soxs_mbias

The [`soxs_mbias`](../../../recipes/soxs_mbias.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] mbias <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
```

To rerun a previously executed `soxs_mbias` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_mbias` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe mbias sof/2019.10.30T09.39.11.127_VIS_1X1_SLOW_MBIAS_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml -x
```

To adjust the default settings for the `soxs_mbias` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_mbias` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`table_soxs_mbias_parameters`.

:::{include} ../../../recipes/parameters/soxs_mbias.md
:::

Product files are written in the `products/soxs_mbias`, and QC plots are in the `qc/soxs_mbias` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_mbias` are found in {numref}`table_soxs_mbias_qc`. 

:::{include} ../../../recipes/qcs/soxs_mbias.md
:::

## Input

:::{include} ../../../recipes/inputs/soxs_mbias.md
:::
