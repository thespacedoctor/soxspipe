:::{table} Quality Control metrics calculated in the `soxs_spatial_solution` recipe.
:name: soxs_spatial_solution_qc
:class: longtable
:widths: 15, 45, 10, 30

| Label | Description                                                          | Unit  | Acceptable Range |
|:----|:-------------------------------------------------------------------|:-----|:-----|
| `DETLINES CLIP NUM` | Total number of detected lines clipped during solution fitting       | lines | - |
| `DETLINES NUM` | Number of lines detected in single pinhole frame                     | lines | - |
| `DETLINES FRAC` | Proportion of input line-list lines detected on single pinhole frame |       | VIS: [0.95,1.0]| 
| `DETLINES TOT` | Total number of line in single line-list                             | lines | - |
| `GOODLINES FRAC` | Proportion of good, unclipped lines in single pinhole frame | lines | VIS: [0.55,1.0] NIR: [0.4,1.0] |
| `PINHOLE COUNT MIN` | The minimum number of pinholes found across all orders (should be 9) | pinholes | VIS: 9, NIR: 9 |
| `X RES MAX`              | Maximum residual in fitting of pinholes along x-axis (global)  | pixels    | - |
| `X RES MIN`              | Minimum residual in fitting of pinholes along x-axis (global)  | pixels    | - |
| `X RES MEDIAN` | Median of absolute residuals in fitting of pinholes along x-axis (global) | pixels | - |
| `X RES SD`              | Std-dev of residuals in fitting of pinholes along x-axis (global)  | pixels    | VIS: [0,0.25], NIR: [0,0.42] |
| `Y RES MAX`              | Maximum residual in fitting of pinholes along y-axis (global)  | pixels    | - |
| `Y RES MIN`              | Minimum residual in fitting of pinholes along y-axis (global)  | pixels    | - |
| `Y RES MEDIAN` | Median of absolute residuals in fitting of pinholes along y-axis (global) | pixels | - |
| `Y RES SD`              | Std-dev of residuals in fitting of pinholes along y-axis (global)  | pixels    | VIS: [0,0.48], NIR: [0,0.12] |
| `XY RES MAX`              | Maximum residual in fitting of pinholes (global)  | pixels    | - |
| `XY RES MIN`              | Minimum residual in fitting of pinholes (global)  | pixels    | - |
| `XY RES MEDIAN` | Median residual in fitting of pinholes (global) | pixels | VIS: [0,0.7], NIR: [0,0.3] |
| `XY RES SD`              | Std-dev of residuals in fitting of pinholes (global)  | pixels    | VIS: [0,0.5], NIR: [0,0.4] |
| `X DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations along the x-axis (global)  | pixels    | - |
| `X DIFF SD` | Standard deviation in difference between observed and predicted pinhole detector locations along the x-axis (global) | pixels | - |
| `Y DIFF MEDIAN` | SMedian difference between observed and predicted pinhole detector locations along the y-axis (global)  | pixels    | - |
| `Y DIFF SD` | Standard deviation in difference between observed and predicted pinhole detector locations along the y-axis (global) | pixels | - |
| `XY DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations (global)  | pixels    | - |
| `XY DIFF SD` | Standard deviation of difference between observed and predicted pinhole detector locations (global) | pixels | - |
| `FWHM PIN MEDIAN` |  Median FWHM of detected lines in pinhole frames (global) | pixels    | VIS: [1.60,1.72], NIR: [1.55,1.85] |
| `FWHM PIN SD` | Standard deviation in FWHM of detected lines in pinhole frames (global) | pixels | VIS: [0.0,0.16], NIR: [0.0,0.145] |
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
| `X DIFF SD O<N>` | Standard deviation in difference between observed and predicted pinhole detector locations along the x-axis (order N) | pixels | - |
| `Y DIFF MEDIAN O<N>` | SMedian difference between observed and predicted pinhole detector locations along the y-axis (order N)  | pixels    | - |
| `Y DIFF SD O<N>` | Standard deviation in difference between observed and predicted pinhole detector locations along the y-axis (order N) | pixels | - |
| `XY DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations (order N)  | pixels    | - |
| `FWHM PIN MEDIAN O<N>` |  Median FWHM of detected lines in pinhole frames (order N) | pixels    | - |
| `FWHM PIN SD O<N>` | Standard-deviation in FWHM of detected lines in pinhole frames (order N) | pixels | - |
| `R PIN MEDIAN O<N>` | Median spectral resolution measured from detected lines in pinhole frames (order N) | | - |

:::


