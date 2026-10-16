:::{table} Quality Control metrics calculated in the `soxs_stare` recipe. 
:name: soxs_stare_qc
:class: longtable
:widths: 15, 45, 10, 30

| Label                | Description                                                | Unit      | Acceptable Range |
| -------------------- | ---------------------------------------------------------- | --------- | --------- |
| `BADPIX FRAC`      | Fraction of bad pixels                                     |           | - |
| `INNER ORDER PIX MEAN` | Mean inner-order pixel value                          | electrons | - |
| `INNER ORDER PIX SUM`  | Sum of all inner-order pixel values                   | electrons | - |
| `BADPIX NUM`         | Number of bad pixels                                       |           | - |
| `N ORDERS`              | Number of orders containing an object trace                |           | - |
| `SAMPLES DET NUM`                | Number of samples where a continuum is detected            |           | - |
| `SAMPLES DET FRAC`                | Fraction of samples where a continuum is detected        |           | - |
| `SAMPLES TOT NUM`                | Total number of samples along orders                       |           | - |
| `SAMPLES CLIP NUM`                | Number of continuum sample clipped during solution fitting |           | - |
| `SAMPLES CLIP FRAC`                | Fraction of detected continuum samples clipped during solution fitting |           | - |
| `X RES MAX`              | Maximum residual in continuum fit along x-axis        | px        | - |
| `X RES MIN`              | Minimum residual in continuum fit along x-axis        | px        | - |
| `X RES SD`              | Std-dev of residual continuum fit along x-axis        | px        | - | 
| `Y RES MAX`              | Maximum residual in continuum fit along y-axis        | px        | - |
| `Y RES MIN`              | Minimum residual in continuum fit along y-axis        | px        | - |
| `Y RES SD`              | Std-dev of residual continuum fit along y-axis        | px        | - |
| `SNR MEDIAN` | The median signal-to-noise ratio across all orders (global) |  | - |
| `SNR MEDIAN O<N>` | The median signal-to-noise ratio in order N |  | - |
| `SLIT DRIFT` | Fitted slit position of the object trace 5% in from the red end of the reddest order, minus that 5% in from the blue end of the bluest order. Signed. Fallback orders are skipped. Not written if every order fell back | arcsec | - |
| `SLIT ZERO WLEN` | Wavelength where a fitted order curve crosses slit position 0. A sign change between two orders does not count. Written only if such a crossing exists. If there are several, the bluest is recorded | nm | - |
| `N ORDERS SLIT EDGE` | Number of orders where the object spills over the slit edge. Flux beyond the slit edge is lost and the extraction in these orders is unreliable. An order spills when more than 20% of its wavelength slices are cut off at the slit edge | | - |
| `SLIT ILLUM OFFSET O<N>` | Along-slit offset of the sky image from the master flat (order N). Written for each order with a usable offset, and only if `slit-illumination: measure-offset` is on. Positive when the sky image sits toward +SLIT | arcsec | - |
| `SLIT ILLUM APPLIED` | 1 if the slit offsets were used to mask the sky model, 0 if no order had a usable offset | | - |
| `SLIT ILLUM N ORDERS` | Number of orders with a usable slit offset | | - |


:::

In the labels that end in `O<N>`, `<N>` is the order number. For the SOXS UV-VIS arm the order numbers are 1 (r), 2 (u), 3 (g), and 4 (i). The same numbers are used in the products, plots, and logs. The UV-VIS spectrograph page lists the wavelength range of each order.
