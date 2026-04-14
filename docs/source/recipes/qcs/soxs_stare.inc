:::{table} Quality Control metrics calculated in the `soxs_stare` recipe. Note that the acceptable ranges are enforced only for standard stars, not for science frames.
:name: soxs_stare_qc

| Label                | Description                                                | Unit      | Acceptable Range |
| -------------------- | ---------------------------------------------------------- | --------- | --------- |
| `BADPIX FRAC`      | Fraction of bad pixels                                     |           | VIS: [0.0,0.01], NIR: [0.00,0.01] |
| `INNER ORDER PIX MEAN` | Mean inner-order pixel value                          | electrons | VIS: [-40,60], NIR: [-40,60] |
| `INNER ORDER PIX SUM`  | Sum of all inner-order pixel values                   | electrons | - |
| `BADPIX NUM`         | Number of bad pixels                                       |           | - |
| `N ORDERS`              | Number of orders containing an object trace                |           | VIS: 4 , NIR: 15 |
| `SAMPLES DET NUM`                | Number of samples where a continuum is detected            |           | - |
| `SAMPLES DET FRAC`                | Fraction of samples where a continuum is detected        |           | VIS: [0.93,1.0] NIR: [0.7,1.0] |
| `SAMPLES TOT NUM`                | Total number of samples along orders                       |           | - |
| `SAMPLES CLIP NUM`                | Number of continuum sample clipped during solution fitting |           | - |
| `SAMPLES CLIP FRAC`                | Fraction of detected continuum samples clipped during solution fitting |           | VIS: [0.0,0.2], NIR: [0.0,0.3] |
| `X RES MAX`              | Maximum residual in continuum fit along x-axis        | px        | - |
| `X RES MIN`              | Minimum residual in continuum fit along x-axis        | px        | - |
| `X RES SD`              | Std-dev of residual continuum fit along x-axis        | px        | VIS: [0.0,0.2] | 
| `Y RES MAX`              | Maximum residual in continuum fit along y-axis        | px        | - |
| `Y RES MIN`              | Minimum residual in continuum fit along y-axis        | px        | - |
| `Y RES SD`              | Std-dev of residual continuum fit along y-axis        | px        | NIR: [0.0,0.5] |
| `EFF MEDIAN` | The median efficiency (global) |  | VIS: [0.04,0.30], NIR: [0.02,0.2] |
| `SNR MEDIAN` | The median signal-to-noise ratio across all orders (global) |  | VIS: [50,1000], NIR: [7,200] |
| `EFF MEDIAN O<N>` | The median efficiency (order N) |  | - |
| `SNR MEDIAN O<N>` | The median signal-to-noise ratio across all orders (order N) |  | - |


:::


==UPDATE QCs==
