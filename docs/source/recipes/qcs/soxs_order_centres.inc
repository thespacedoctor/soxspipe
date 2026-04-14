:::{table} Quality Control metrics calculated in the `soxs_order_centres` recipe.
:name: soxs_order_centres_qc

| Label   | Description                                                | Unit | Acceptable Range |
| ------- | ---------------------------------------------------------- | ---- | ---- |
| `N ORDERS` | Number of order centre traces found                        |      | - |
| `SAMPLES DET NUM`                | Number of samples where a continuum is detected            |           | - |
| `SAMPLES DET FRAC`                | Fraction of samples where a continuum is detected        |           | VIS: [0.85,1.0], NIR: [0.80,1.0] |
| `SAMPLES TOT NUM`                | Total number of samples along orders                       |           | - |
| `SAMPLES CLIP NUM`                | Number of continuum sample clipped during solution fitting |           | - |
| `SAMPLES CLIP FRAC `               | Fraction of detected continuum samples clipped during solution fitting |           | - |
| `X RES MAX` | Maximum residual in order centre fit along x-axis     | px   | - |
| `X RES MIN` | Minimum residual in order centre fit along x-axis     | px   | - |
| `X RES SD` | Std-dev of residual order centre fit along x-axis     |      | VIS: [0,0.1]| 
| `Y RES MAX` | Maximum residual in order centre fit along y-axis    | px   | - |
| `Y RES MIN` | Minimum residual in order centre fit along y-axis    | px   | - |
| `Y RES SD` | Std-dev of residual order centre fit along y-axis    |      | NIR: [0,0.1] |

:::


==UPDATE QCs==
