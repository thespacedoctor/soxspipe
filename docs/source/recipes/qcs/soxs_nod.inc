:::{table} Quality Control metrics calculated in the `soxs_nod` recipe. 
:name: soxs_nod_qc
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
| `SKY TUNING` | `skylines` if the wavelength solution was tuned to the sky lines, `skipped: attached arc` if a SOXS arc taken within 7 minutes and 5 degrees of rotator angle of the science frame made the tuning unnecessary | | - |
| `SKY SHIFT O<N>` | Wavelength shift that the skyline tuning added to the wavelength solution of order N. Written for every order. For SOXS VIS orders 2 (u) and 3 (g) it is the modelled shift if the model applied, and 0 if not | nm | - |
| `SKY SHIFT O<N> PX` | The same shift as detector motion along the dispersion axis, using the signed dispersion of the order. Positive means the features moved to larger pixel numbers. Not written for an order with no usable dispersion | px | - |
| `SKY LINE SNR O<N>` | Combined signal-to-noise of the catalogue sky lines in order N. Written for each SOXS VIS order measured from its own sky lines (orders 1 and 4). The skyline tuning does not apply the shift of an order below 5.4 | | - |
| `SLIT ALIGN FALLBACK NUM` | Number of slit shifts that took the median shift of the other orders or frames, because no reliable object trace was found. Written to the header | | - |
| `UG SHIFT SOURCE` | How the SOXS VIS u and g shifts were found: `model` (predicted from the r and i shifts and the rotator angle) or `none` (u and g left unshifted). The comment gives the reason. Written only for SOXS VIS | | - |

:::

In the labels that end in `O<N>`, `<N>` is the order number. For the SOXS UV-VIS arm the order numbers are 1 (r), 2 (u), 3 (g), and 4 (i). The same numbers are used in the products, plots, and logs. The UV-VIS spectrograph page lists the wavelength range of each order.
