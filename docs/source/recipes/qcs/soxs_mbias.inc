:::{table} Quality Control metrics calculated in the `soxs_mbias` recipe.
:name: table_soxs_mbias_qc
:widths: 15, 45, 10, 30

| Label | Description | Unit | Acceptable Range |
| :------------ | :----------- | ------------- | ------------- |
| `BADPIX FRAC` |                Fraction of bad pixels |   -       | - |
| `MBIAS MEDIAN` | Median flux level of master bias |electrons | VIS: [0, 2000] |
| `BADPIX NUM` | Number of bad pixels |  -        | - |
| `FPN FRACMAX` | Maximum periodic pattern noise ratio in raw bias frames | - | VIS: [0,3] |
| `RAW RON` | RON in single BIAS |electrons | VIS: [0,15] |
| `MASTER RON` | Combined RON in MBIAS |electrons | VIS: [0,10] |
| `STRUCTX` | Slope of BIAS in X direction | - | VIS: [-3,1] |
| `STRUCTY` | Slope of BIAS in Y direction | - | VIS: [0.0,0.6] |



:::
