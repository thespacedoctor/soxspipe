:::{table} Quality Control metrics calculated in the `soxs_disp_solution` recipe.
:name: soxs_disp_solution_qc

| Label | Description                                                          | Unit  | Acceptable Range |
| :---- | :------------------------------------------------------------------- | ----- | ----- |
| `DETLINES CLIP NUM` | Total number of detected lines clipped during solution fitting       | lines | VIS: , NIR: |
| `DETLINES NUM` | Number of lines detected in single pinhole frame                     | lines | VIS: , NIR: |
| `DETLINES FRAC` | Proportion of input line-list lines detected on single pinhole frame |       | VIS: [0.9,1.01], NIR: |
| `DETLINES TOT` | Total number of line in single line-list                             | lines | VIS: , NIR: |
| `GOODLINES FRAC` | Proportion of good, unclipped lines in single pinhole frame | lines | VIS: , NIR: |
| `X RES MAX`              | Maximum residual in fitting of pinholes along x-axis (global)  | pixels    | VIS: , NIR: |
| `X RES MIN`              | Minimum residual in fitting of pinholes along x-axis (global)  | pixels    | VIS: , NIR: |
| `X RES SD`              | Std-dev of residuals in fitting of pinholes along x-axis (global)  | pixels    | VIS: , NIR: |
| `Y RES MAX`              | Maximum residual in fitting of pinholes along y-axis (global)  | pixels    | VIS: , NIR: |
| `Y RES MIN`              | Minimum residual in fitting of pinholes along y-axis (global)  | pixels    | VIS: , NIR: |
| `Y RES SD`              | Std-dev of residuals in fitting of pinholes along y-axis (global)  | pixels    | VIS: , NIR: |
| `XY RES MAX`              | Maximum residual in fitting of pinholes (global)  | pixels    | VIS: , NIR: |
| `XY RES MIN`              | Minimum residual in fitting of pinholes (global)  | pixels    | VIS: , NIR: |
| `XY RES SD`              | Std-dev of residuals in fitting of pinholes (global)  | pixels    | VIS: , NIR: |
| `X DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations along the x-axis (global) | pixels    | VIS: , NIR: |
| `Y DIFF MEDIAN` | SMedian difference between observed and predicted pinhole detector locations along the y-axis (global)  | pixels    | VIS: , NIR: |
| `XY DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations (global)  | pixels    | VIS: , NIR: |
| `FWHM PIN MEDIAN` |  Median FWHM of detected lines in pinhole frames (global) | pixels    | VIS: [1.60,1.75], NIR: [1.75,1.89] |
| `FWHM PIN SD` | Standard-deviation in FWHM of detected lines in pinhole frames (global) | pixels | VIS: [0.1,0.2], NIR: [0,0.211] |
| `R PIN MEDIAN` | Median spectral resolution measured from detected lines in pinhole frames  (global) | | VIS: , NIR: |
| `X RES MAX O<N>`              | Maximum residual in fitting of pinholes along x-axis (order N)  | pixels    | VIS: , NIR: |
| `X RES MIN O<N>`              | Minimum residual in fitting of pinholes along x-axis (order N)  | pixels    | VIS: , NIR: |
| `X RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along x-axis (order N)  | pixels    | VIS: , NIR: |
| `Y RES MAX O<N>`              | Maximum residual in fitting of pinholes along y-axis (order N)  | pixels    | VIS: , NIR: |
| `Y RES MIN O<N>`              | Minimum residual in fitting of pinholes along y-axis (order N)  | pixels    | VIS: , NIR: |
| `Y RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along y-axis (order N)  | pixels    | VIS: , NIR: |
| `XY RES MAX O<N>`              | Maximum residual in fitting of pinholes (order N)  | pixels    | VIS: , NIR: |
| `XY RES MIN O<N>`              | Minimum residual in fitting of pinholes (order N)  | pixels    | VIS: , NIR: |
| `XY RES SD O<N>`              | Std-dev of residuals in fitting of pinholes along (order N)  | pixels    | VIS: , NIR: |
| `X DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations along the x-axis (order N)  | pixels    | VIS: , NIR: |
| `Y DIFF MEDIAN O<N>` | SMedian difference between observed and predicted pinhole detector locations along the y-axis (order N)  | pixels    | VIS: , NIR: |
| `XY DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations (order N)  | pixels    | VIS: , NIR: |
| `FWHM PIN MEDIAN O<N>` |  Median FWHM of detected lines in pinhole frames (order N) | pixels    | VIS: , NIR: |
| `FWHM PIN SD O<N>` | Standard-deviation in FWHM of detected lines in pinhole frames (order N) | pixels | VIS: , NIR: |
| `R PIN MEDIAN O<N>` | Median spectral resolution measured from detected lines in pinhole frames  (order N) | | VIS: , NIR: |

:::




==UPDATE QCs==
