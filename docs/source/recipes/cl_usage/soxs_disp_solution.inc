The [`soxs_disp_solution`](../../../recipes/soxs_disp_solution.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] disp_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<od>]
```

To rerun a previously executed `soxs_disp_solution` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_disp_solution` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe disp_sol sof/2019.10.31T09.56.05.151_UVB_1X1_FAST_DISP_SOL_THAR_40.0S_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml  -x
```

While the default polynomial fitting orders have been carefully tuned to robustly reduce most data, it is possible to execute this recipe while providing the dispersion-solution polynomial fitting orders through the command line, with command line settings overriding those in the YAML settings file. For instance, to attempt a dispersion solution fit using 3rd (x) and 4th-order (y) spectral-order (`oo`) component and a 5th-order (for both x and y) wavelength (`ww`) component:

```bash
soxspipe disp_sol sof/2019.10.31T09.56.05.151_UVB_1X1_FAST_DISP_SOL_THAR_40.0S_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml --poly=3455
```

To adjust the default settings for the `soxs_disp_solution` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_disp_solution` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_disp_solution_parameters`.

Product files are written in the `products/soxs_disp_solution`, and QC plots are in the `qc/soxs_disp_solution` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_disp_solution` are found in {numref}`soxs_disp_solution_qc` and a typical QC plot in {numref}`soxs_disp_solution_qc_fig`.

### Reduction Tips

If this recipe fails because a fit for the dispersion solution is not found, the first parameter to adjust is `poly-fitting-residual-clipping-sigma`. Try setting a value lower than the default value by 0.5 and running the recipe again. If the recipe still fails, repeat the reduction in steps of 0.5 down to 3 sigma.

You can also try adjusting the polynomial fitting orders in `order-deg`, `wavelength-deg` and `slit-deg`. However, please be advised the pipeline itself will dynamically adjust these values if they fail to fit the default set. It will slowly reduce the orders and refit until it finds a fit or decides a fit can not be found (after five iterations of decreasing the orders).
