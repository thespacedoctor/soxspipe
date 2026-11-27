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
| `COSMIC RAY PIX NUM` | Total number of pixels L.A.Cosmic flagged as cosmic rays, over all raw frames. Written only if `use_lacosmic` is on. The count in each frame is in the header of its extension of the `LACOSMIC_MASK` product | pixel | - |
| `N ORDERS`              | Number of orders containing an object trace                |           | - |
| `N ORDERS FAINT TRACE` | Number of faint SOXS VIS orders whose trace was transferred from the neighbouring order. Written to the header only when above 0. See [`detect_continuum`](../utils/detect_continuum.md) | | - |
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
| `SKY TUNING` | `skylines` if the wavelength solution was tuned to the sky lines, `skipped: attached arc` if a SOXS arc taken within 7 minutes of the science exposures (end of the earlier exposure to start of the later) and within 5 degrees of rotator angle made the tuning unnecessary | | - |
| `SKY SHIFT O<N>` | Wavelength shift that the skyline tuning added to the wavelength solution of order N. Written for every order. For SOXS VIS orders 2 (u) and 3 (g) it is the modelled shift if the model applied, and 0 if not | nm | - |
| `SKY SHIFT O<N> PX` | The same shift as detector motion along the dispersion axis, using the signed dispersion of the order. Positive means the features moved to larger pixel numbers. Not written for an order with no usable dispersion | px | - |
| `SKY SHIFT O<N> SRC` | Where the shift of order N came from: `measured` (from the sky lines of the order), `model` (predicted: the SOXS VIS r, u and g orders, and the NIR orders whose shift comes from the cross-order global fit) or `none` (no shift was applied, so the 0 in `SKY SHIFT O<N>` is a placeholder). The comment gives the reason. Written for every order that has a `SKY SHIFT O<N>` card | | - |
| `SKY SHIFT O<N> ERR` | The uncertainty of the shift of order N, set by where the shift came from (`SKY SHIFT O<N> SRC`): for `measured`, the standard error of the median of the kept sky-line offsets; for the SOXS VIS r order predicted from i, or combined with that prediction, the r-from-i model rms times the dispersion (the combined value is the inverse-variance uncertainty); for the u and g orders from the transfer model, the leave-one-night-out rms of the model times the dispersion; for `none`, and for every order when the tuning is switched off, the per-arm `untuned` value of `sky-shift-uncertainty-nm`; for an attached arc, its `attached-arc` value. A nod or offset recipe keeps the largest value over its extractions. `SPEC_SYE` of the merged spectrum is the largest of these values. Written for every order of a SOXS VIS frame, and not yet for NIR | nm | - |
| `SKY LINE SNR O<N>` | Combined signal-to-noise of the catalogue sky lines in order N. Written for each SOXS VIS order measured from its own sky lines (orders 1 and 4). The skyline tuning does not apply the shift of an order below 5.4 | | - |
| `TELL SHIFT` | Median shift of the star's telluric absorption lines (fitted minus catalogue wavelength, air) over the lines of the measured bands, after the sky-line tuning. For a star in a wide slit it is mostly the centring offset of the star across the slit width. It is informational: it changes no wavelength and has no acceptable range. Reads `not measured` when no band is measured. See the [horne_extraction](../utils/horne_extraction.md) utility | nm | - |
| `TELL SHIFT ERR` | Robust scatter of the line shifts (1.4826 times the median absolute deviation) divided by √N. Reads `not measured` when no band is measured | nm | - |
| `TELL NLINES` | Number of telluric lines used for `TELL SHIFT`. 0 when no band is measured | | - |
| `TELL SHIFT <BAND>`, `TELL SHIFT ERR <BAND>`, `TELL NLINES <BAND>` | The same three values for one telluric band. A band counts as measured when it has at least 5 lines with a centroid error below 0.3 Å, and reads `not measured` otherwise. Written to the QC table only, not to the header | nm | - |
| `SLIT ALIGN FALLBACK NUM` | Number of slit shifts that took the median shift of the other orders or frames, because no reliable object trace was found. Written to the header | | - |
| `NOD COMBINED` | 1 if the stacked A-B and B-A differences were combined and extracted once, 0 if they were extracted separately. Written to the header | | - |
| `NOD COMBINE FALLBACK NUM` | Number of orders whose shift between the nod positions took the median of the other orders, because a nod position had no reliable object trace. 0 when `NOD COMBINED` is 0. Written to the header | | - |
| `NOD HOST CONTAM` | Number of narrow features at the negative trace of the A-B difference that the object does not explain. A feature is a run of at least 2 detector lines where both the residual and the negative trace are above `host-contamination-sigma` times the local noise. Features within 0.2 nm are counted once. Written to the header. It has no acceptable range and never fails the recipe. See [`soxs_nod`](../soxs_nod.md) | | - |
| `UG SHIFT SOURCE` | How the SOXS VIS u and g shifts were found: `model` (predicted from the r and i shifts and the rotator angle) or `none` (u and g left unshifted). The comment gives the reason. Written only for SOXS VIS | | - |

:::

In the labels that end in `O<N>`, `<N>` is the order number. For the SOXS UV-VIS arm the order numbers are 1 (r), 2 (u), 3 (g), and 4 (i). The same numbers are used in the products, plots, and logs. The UV-VIS spectrograph page lists the wavelength range of each order.
