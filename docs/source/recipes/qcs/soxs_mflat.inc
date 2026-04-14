:::{table} Quality Control metrics calculated in the `soxs_mflat` recipe.
:name: soxs_mflat_qc

| Label                | Description                                           | Unit      | Acceptable Range |
| -------------------- | ----------------------------------------------------- | --------- | --------- |
| `COLDPIX FRAC`      | Fraction of colf pixels                                |           | - |
| `INNER ORDER PIX MEAN` | Mean inner-order pixel value                     | electrons | VIS: [0.9,1.1], NIR: [0.9,1.1] |
| `INNER ORDER PIX SUM`  | Sum of all inner-order pixel values              | electrons | VIS: [100000,700000], NIR: [1000000,1200000] |
| `COLDPIX NUM `        | Number of cold pixels                                  |           | VIS: [0.0,0.01], NIR: [0.0,0.01] |
| `N LOW SENS`             | Number of low-sensitivity pixels found in master flat | pixels    | - |
| `ORDEXP10`             | 10th percentile inter-order flux                 | electrons | - |
| `ORDEXP50`             | 50th percentile inter-order flux                 | electrons | - |
| `ORDEXP90`             | 90th percentile inter-order flux                 | electrons | VIS: [10000,150000], NIR: [30000,60000] |
| `X RES MAX`              | Maximum residual in order edge fit along x-axis  | pixels    | - |
| `X RES MIN`              | Minimum residual in order edge fit along x-axis  | pixels    | - |
| `X RES SD`              | Std-dev of residual order edge fit along x-axis  | pixels    | VIS: [0,0.37], NIR: [0,0.2] |



:::

==UPDATE QCs==
