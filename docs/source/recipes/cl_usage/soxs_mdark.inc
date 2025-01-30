The [`soxs_mdark`](../../../recipes/soxs_mdark.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] mdark <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
```

To rerun a previously executed `soxs_mdark` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_mdark` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe mdark sof/2019.10.30T10.27.28.6175_NIR_MDARK_60.0S_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml  -x
```

To adjust the default settings for the `soxs_mdark` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_mdark` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_mdark_parameters`.

Product files are written in the `products/soxs_mdark`, and QC plots are in the `qc/soxs_mdark` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_mdark` are found in {numref}`soxs_mdark_qc`. 

### Reduction Tips

As with `soxs_mbias`, the `soxs_mdark` recipe is one of the simpler recipes, and therefore, it is typically very robust against failures. To reduce the time it takes to run the recipe (at the expense of allowing a few more outlying pixels into your master bias frame), you can experiment with reducing both `frame-clipping-iterations` and `stacked-clipping-iterations` (although speed gains will be very moderate).

