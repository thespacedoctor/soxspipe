The [`soxs_mflat`](../../../recipes/soxs_mflat.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] mflat <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
```

To rerun a previously executed `soxs_mflat` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_mflat` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe mflat sof/2020.10.22T12.59.53.634_VIS_1X1_FAST_MFLAT_23.6953S_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml  -x
```

To adjust the default settings for the `soxs_mflat` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_mflat` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_mflat_parameters`.


Product files are written in the `products/soxs_mflat`, and QC plots are in the `qc/soxs_mflat` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_mflat` are found in {numref}`soxs_mflat_qc` and a typical QC plot in {numref}`soxs_mflat_qc_fig`.

### Reduction Tips

If this recipe fails during the fitting of the order edge, the first parameter to adjust is `poly-fitting-residual-clipping-sigma`. Try reducing this to 3-5 sigma and rerun the recipe to see if a fit is found. If the fit still fails, next try and increase the `slice-width-for-edge-detection` to 5-9 pixels to give the code a better chance of detecting the order-edges. 

If the previous adjustments have not helped, you can slowly reduce `min-percentage-threshold-for-edge-detection` and/or increase `max-percentage-threshold-for-edge-detection` by 2-3 percent and attempt a refit.

You can also try adjusting the polynomial fitting orders in `order-deg` and `disp-axis-deg`. However, please be advised the pipeline itself will dynamically adjust these values if they fail to fit the default set. It will slowly reduce the orders and refit until it finds a fit or decides a fit can not be found (after five iterations of decreasing the orders).
