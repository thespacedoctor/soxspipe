:::{table} Quality Control metrics calculated in the `soxs_spatial_solution` recipe.
:name: soxs_spatial_solution_qc

| Label | Description                                                          | Unit  |
| :---- | :------------------------------------------------------------------- | ----- |
| `DETLINES CLIP NUM` | Total number of detected lines clipped during solution fitting       | lines |
| `DETLINES NUM` | Number of lines detected in single pinhole frame                     | lines |
| `DETLINES FRAC` | Proportion of input line-list lines detected on single pinhole frame |       |
| `DETLINES TOT` | Total number of line in single line-list                             | lines |
| `X RES MAX`              | Maximum residual in fitting of pinholes along x-axis (global)  | pixels    |
| `X RES MIN`              | Minimum residual in fitting of pinholes along x-axis (global)  | pixels    |
| `X RES RMS`              | Std-dev of residuals in fitting of pinholes along x-axis (global)  | pixels    |
| `Y RES MAX`              | Maximum residual in fitting of pinholes along y-axis (global)  | pixels    |
| `Y RES MIN`              | Minimum residual in fitting of pinholes along y-axis (global)  | pixels    |
| `Y RES RMS`              | Std-dev of residuals in fitting of pinholes along y-axis (global)  | pixels    |
| `XY RES MAX`              | Maximum residual in fitting of pinholes (global)  | pixels    |
| `XY RES MIN`              | Minimum residual in fitting of pinholes (global)  | pixels    |
| `XY RES RMS`              | Std-dev of residuals in fitting of pinholes (global)  | pixels    |
| `X DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations along the x-axis (global)  | pixels    |
| `Y DIFF MEDIAN` | SMedian difference between observed and predicted pinhole detector locations along the y-axis (global)  | pixels    |
| `XY DIFF MEDIAN` | Median difference between observed and predicted pinhole detector locations (global)  | pixels    |
| `FWHM PIN MEDIAN` |  Median FWHM of detected lines in pinhole frames (global) | pixels    |
| `R PIN MEDIAN` | Median spectral resolution measured from detected lines in pinhole frames  (global) | |
| `X RES MAX O<N>`              | Maximum residual in fitting of pinholes along x-axis (order N)  | pixels    |
| `X RES MIN O<N>`              | Minimum residual in fitting of pinholes along x-axis (order N)  | pixels    |
| `X RES RMS O<N>`              | Std-dev of residuals in fitting of pinholes along x-axis (order N)  | pixels    |
| `Y RES MAX O<N>`              | Maximum residual in fitting of pinholes along y-axis (order N)  | pixels    |
| `Y RES MIN O<N>`              | Minimum residual in fitting of pinholes along y-axis (order N)  | pixels    |
| `Y RES RMS O<N>`              | Std-dev of residuals in fitting of pinholes along y-axis (order N)  | pixels    |
| `XY RES MAX O<N>`              | Maximum residual in fitting of pinholes (order N)  | pixels    |
| `XY RES MIN O<N>`              | Minimum residual in fitting of pinholes (order N)  | pixels    |
| `XY RES RMS O<N>`              | Std-dev of residuals in fitting of pinholes along (order N)  | pixels    |
| `X DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations along the x-axis (order N)  | pixels    |
| `Y DIFF MEDIAN O<N>` | SMedian difference between observed and predicted pinhole detector locations along the y-axis (order N)  | pixels    |
| `XY DIFF MEDIAN O<N>` | Median difference between observed and predicted pinhole detector locations (order N)  | pixels    |
| `FWHM PIN MEDIAN O<N>` |  Median FWHM of detected lines in pinhole frames (order N) | pixels    |
| `R PIN MEDIAN O<N>` | Median spectral resolution measured from detected lines in pinhole frames  (order N) | |

:::






