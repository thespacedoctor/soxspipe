:::{table} Quality Control metrics calculated in the `soxs_disp_solution` recipe.
:name: soxs_disp_solution_qc
:class: longtable
:widths: 15, 45, 10, 30

| Label | Description                                                          | Unit  | Acceptable Range |
|:----|:-------------------------------------------------------------------|:-----|:-----|
| `DETLINES CLIP NUM` | Total number of detected lines clipped during solution fitting       | lines | - |
| `DETLINES NUM` | Number of lines detected in single pinhole frame                     | lines | - |
| `DETLINES FRAC` | Proportion of input line-list lines detected on single pinhole frame |       | VIS: [0.9,1.0], NIR: [0.9,1.0] |
| `DETLINES TOT` | Total number of line in single line-list                             | lines | - |
| `GOODLINES FRAC` | Proportion of good, unclipped lines in single pinhole frame | lines | - |
| `X RES MAX`              | Maximum residual in fitting of pinholes along x-axis (global)  | pixels    | - |
| `X RES MIN`              | Minimum residual in fitting of pinholes along x-axis (global)  | pixels    | - |
| `X RES MEDIAN` | Median of absolute residuals in fitting of pinholes along x-axis (global) | pixels | VIS: [0.0,0.15], NIR: [0.0,0.35] |
| `X RES SD`              | Std-dev of residuals in fitting of pinholes along x-axis (global)  | pixels    | VIS: [0,0.4], NIR: [0,1.0] |
| `Y RES MAX`              | Maximum residual in fitting of pinholes along y-axis (global)  | pixels    | - |
| `Y RES MIN`              | Minimum residual in fitting of pinholes along y-axis (global)  | pixels    | - |
| `Y RES MEDIAN` | Median of absolute residuals in fitting of pinholes along y-axis (global) | pixels | VIS: [0.0,0.9], NIR: [0.0,0.13] |
| `Y RES SD`              | Std-dev of residuals in fitting of pinholes along y-axis (global)  | pixels    | VIS: [0,1.2], NIR: [0,0.25] |
| `XY RES MAX`              | Maximum residual in fitting of pinholes (global)  | pixels    | - |
| `XY RES MIN`              | Minimum residual in fitting of pinholes (global)  | pixels    | - |
| `XY RES MEDIAN` | Median residual in fitting of pinholes (global) | pixels | NIR: [0,1.0] |
| `XY RES SD`              | Std-dev of residuals in fitting of pinholes (global)  | pixels    | VIS: [0,1.5], NIR: [0,1.0] |
| `X DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations along the x-axis (global) | pixels    | - |
| `X DIFF SD` | Std-dev in difference between observed and predicted pinhole detector locations along the x-axis (global) | pixels | VIS: [0,3.0], NIR: [0,1.3] |
| `Y DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations along the y-axis (global)  | pixels    | - |
| `Y DIFF SD` | Std-dev in difference between observed and predicted pinhole detector locations along the y-axis (global) | pixels | VIS: [0,10.0], NIR: [0,1.5] |
| `XY DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations (global)  | pixels    | NIR: [0,6.0] |
| `XY DIFF SD` | Std-dev of difference between observed and predicted pinhole detector locations (global) | pixels | NIR: [0,1.9] |
| `FWHM PIN MEDIAN` |  Median FWHM of detected lines in pinhole frames (global) | pixels    | VIS: [1.60,1.75], NIR: [1.75,1.89] |
| `FWHM PIN SD` | Std-dev in FWHM of detected lines in pinhole frames (global) | pixels | VIS: [0.1,0.2], NIR: [0,0.211] |
| `R PIN MEDIAN` | Median spectral resolution measured from detected lines in pinhole frames  (global) | | - |
| `X RES MAX O<N>`              | Maximum residual in fitting of pinholes along x-axis (order N)  | pixels    | - |
| `X RES MIN O<N>`              | Minimum residual in fitting of pinholes along x-axis (order N)  | pixels    | - |
| `X RES MEDIAN O<N>` | Median of absolute residuals in fitting of pinholes along x-axis (order N) | pixels | - |
| `X RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along x-axis (order N)  | pixels    | - |
| `Y RES MAX O<N>`              | Maximum residual in fitting of pinholes along y-axis (order N)  | pixels    | - |
| `Y RES MIN O<N>`              | Minimum residual in fitting of pinholes along y-axis (order N)  | pixels    | - |
| `Y RES MEDIAN O<N>` | Median of absolute residuals in fitting of pinholes along y-axis (order N) | pixels | - |
| `Y RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along y-axis (order N)  | pixels    | - |
| `XY RES MAX O<N>`              | Maximum residual in fitting of pinholes (order N)  | pixels    | - |
| `XY RES MIN O<N>`              | Minimum residual in fitting of pinholes (order N)  | pixels    | - |
| `XY RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along (order N)  | pixels    | - |
| `X DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations along the x-axis (order N)  | pixels    | - |
| `X DIFF SD O<N>` | Std-dev in difference between observed and predicted pinhole detector locations along the x-axis (order N) | pixels | - |
| `Y DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations along the y-axis (order N)  | pixels    | - |
| `Y DIFF SD O<N>` | Std-dev in difference between observed and predicted pinhole detector locations along the y-axis (order N) | pixels | - |
| `XY DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations (order N)  | pixels    | - |
| `FWHM PIN MEDIAN O<N>` |  Median FWHM of detected lines in pinhole frames (order N) | pixels    | - |
| `FWHM PIN SD O<N>` | Std-dev in FWHM of detected lines in pinhole frames (order N) | pixels | - |
| `R PIN MEDIAN O<N>` | Median spectral resolution measured from detected lines in pinhole frames  (order N) | | - |

:::


