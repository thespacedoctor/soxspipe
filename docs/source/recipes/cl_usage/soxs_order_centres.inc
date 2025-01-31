
The [`soxs_order_centres`](../../../recipes/soxs_order_centres.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] order_centres <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<ooww>]
```

To rerun a previously executed `soxs_order_centres` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_order_centres` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe order_centres sof/2019.10.31T10.02.46.516_UVB_1X1_FAST_ORDER_LOCATIONS_QLAMP_UVB_HIGH_30.0S_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml 
```

While the default polynomial fitting orders have been carefully tuned to robustly reduce most data, it is possible to execute this recipe while providing the order-trace polynomial fitting orders through the command line, with command line settings overriding those in the YAML settings file. For instance, to attempt a order trace fit using 3rd (x) and 6th-order (y) spectral-order (`oo`) component and a 4th-order (x) and 3rd-order (y) wavelength (`ww`) component:

```bash
soxspipe order_centres sof/2019.10.31T10.02.46.516_UVB_1X1_FAST_ORDER_LOCATIONS_QLAMP_UVB_HIGH_30.0S_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml  --poly=3643
```

To adjust the default settings for the `soxs_order_centres` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_order_centres` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_order_centres_parameters`.

Product files are written in the `products/soxs_order_centres`, and QC plots are in the `qc/soxs_order_centres` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_order_centres` are found in {numref}`soxs_order_centres_qc` and a typical QC plot in {numref}`soxs_order_centres_qc_fig`.

### Reduction Tips

If this recipe fails the order centre traces, the first parameter to adjust is `poly-fitting-residual-clipping-sigma`. Try reducing this to 3-5 sigma and rerun the recipe to see if a fit is found. If the fit still fails, next try and increase the `slice-width` to 5-9 pixels to give the code a better chance of detecting the trace. 

You can also try adjusting the polynomial fitting orders in `order-deg` and `disp-axis-deg`. However, please be advised the pipeline itself will dynamically adjust these values if they fail to fit the default set. It will slowly reduce the orders and refit until it finds a fit or decides a fit can not be found (after five iterations of decreasing the orders).
