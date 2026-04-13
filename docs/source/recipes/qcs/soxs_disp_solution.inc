:::{table} Quality Control metrics calculated in the `soxs_disp_solution` recipe.
:name: soxs_disp_solution_qc

| Label | Description                                                          | Unit  | Acceptable Range |
| :---- | :------------------------------------------------------------------- | ----- | ----- |
| `DETLINES CLIP NUM` | Total number of detected lines clipped during solution fitting       | lines | VIS: , NIR: |
| `DETLINES NUM` | Number of lines detected in single pinhole frame                     | lines | VIS: , NIR: |
| `DETLINES FRAC` | Proportion of input line-list lines detected on single pinhole frame |       | VIS: [0.9,1.0], NIR: [0.9,1.0] |
| `DETLINES TOT` | Total number of line in single line-list                             | lines | VIS: , NIR: |
| `GOODLINES FRAC` | Proportion of good, unclipped lines in single pinhole frame | lines | VIS: , NIR: |
| `X RES MAX`              | Maximum residual in fitting of pinholes along x-axis (global)  | pixels    | VIS: , NIR: |
| `X RES MIN`              | Minimum residual in fitting of pinholes along x-axis (global)  | pixels    | VIS: , NIR: |
| `X RES MEDIAN` | Median of absolute residuals in fitting of pinholes along x-axis (global) | pixels | VIS: [0.0,0.15], NIR: [0.0,0.35] |
| `X RES SD`              | Std-dev of residuals in fitting of pinholes along x-axis (global)  | pixels    | VIS: [0,0.4], NIR: [0,1.0] |
| `Y RES MAX`              | Maximum residual in fitting of pinholes along y-axis (global)  | pixels    | VIS: , NIR: |
| `Y RES MIN`              | Minimum residual in fitting of pinholes along y-axis (global)  | pixels    | VIS: , NIR: |
| `Y RES MEDIAN` | Median of absolute residuals in fitting of pinholes along y-axis (global) | pixels | VIS: [0.0,0.9], NIR: [0.0,0.13] |
| `Y RES SD`              | Std-dev of residuals in fitting of pinholes along y-axis (global)  | pixels    | VIS: [0,1.2], NIR: [0,0.25] |
| `XY RES MAX`              | Maximum residual in fitting of pinholes (global)  | pixels    | VIS: , NIR: |
| `XY RES MIN`              | Minimum residual in fitting of pinholes (global)  | pixels    | VIS: , NIR: |
| `XY RES MEDIAN` | Median residual in fitting of pinholes (global) | pixels | VIS: , NIR: [0,1.0] |
| `XY RES SD`              | Std-dev of residuals in fitting of pinholes (global)  | pixels    | VIS: [0,1.5], NIR: [0,1.0] |
| `X DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations along the x-axis (global) | pixels    | VIS: , NIR: |
| `X DIFF SD` | Standard deviation in difference between observed and predicted pinhole detector locations along the x-axis (global) | pixels | VIS: [0,3.0], NIR: [0,1.3] |
| `Y DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations along the y-axis (global)  | pixels    | VIS: , NIR: |
| `Y DIFF SD` | Standard deviation in difference between observed and predicted pinhole detector locations along the y-axis (global) | pixels | VIS: [0,10.0], NIR: [0,1.5] |
| `XY DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations (global)  | pixels    | VIS: , NIR: [0,6.0] |
| `XY DIFF SD` | Standard deviation of difference between observed and predicted pinhole detector locations (global) | pixels | VIS: , NIR: [0,1.9] |
| `FWHM PIN MEDIAN` |  Median FWHM of detected lines in pinhole frames (global) | pixels    | VIS: [1.60,1.75], NIR: [1.75,1.89] |
| `FWHM PIN SD` | Standard deviation in FWHM of detected lines in pinhole frames (global) | pixels | VIS: [0.1,0.2], NIR: [0,0.211] |
| `R PIN MEDIAN` | Median spectral resolution measured from detected lines in pinhole frames  (global) | | VIS: , NIR: |
| `X RES MAX O<N>`              | Maximum residual in fitting of pinholes along x-axis (order N)  | pixels    | VIS: , NIR: |
| `X RES MIN O<N>`              | Minimum residual in fitting of pinholes along x-axis (order N)  | pixels    | VIS: , NIR: |
| `X RES MEDIAN O<N>` | Median of absolute residuals in fitting of pinholes along x-axis (order N) | pixels | VIS: , NIR: |
| `X RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along x-axis (order N)  | pixels    | VIS: , NIR: |
| `Y RES MAX O<N>`              | Maximum residual in fitting of pinholes along y-axis (order N)  | pixels    | VIS: , NIR: |
| `Y RES MIN O<N>`              | Minimum residual in fitting of pinholes along y-axis (order N)  | pixels    | VIS: , NIR: |
| `Y RES MEDIAN O<N>` | Median of absolute residuals in fitting of pinholes along y-axis (order N) | pixels | VIS: , NIR: |
| `Y RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along y-axis (order N)  | pixels    | VIS: , NIR: |
| `XY RES MAX O<N>`              | Maximum residual in fitting of pinholes (order N)  | pixels    | VIS: , NIR: |
| `XY RES MIN O<N>`              | Minimum residual in fitting of pinholes (order N)  | pixels    | VIS: , NIR: |
| `XY RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along (order N)  | pixels    | VIS: , NIR: |
| `X DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations along the x-axis (order N)  | pixels    | VIS: , NIR: |
| `X DIFF SD O<N>` | Standard deviation in difference between observed and predicted pinhole detector locations along the x-axis (order N) | pixels | VIS: , NIR: |
| `Y DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations along the y-axis (order N)  | pixels    | VIS: , NIR: |
| `Y DIFF SD O<N>` | Standard deviation in difference between observed and predicted pinhole detector locations along the y-axis (order N) | pixels | VIS: , NIR: |
| `XY DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations (order N)  | pixels    | VIS: , NIR: |
| `FWHM PIN MEDIAN O<N>` |  Median FWHM of detected lines in pinhole frames (order N) | pixels    | VIS: , NIR: |
| `FWHM PIN SD O<N>` | Standard-deviation in FWHM of detected lines in pinhole frames (order N) | pixels | VIS: , NIR: |
| `R PIN MEDIAN O<N>` | Median spectral resolution measured from detected lines in pinhole frames  (order N) | | VIS: , NIR: |

:::




==UPDATE QCs==
