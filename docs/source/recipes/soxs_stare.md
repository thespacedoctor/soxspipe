# soxs_stare

:::{include} ./descriptions/soxs_stare.inc
:::


## Input



:::{include} ./inputs/soxs_stare.inc
:::

:::{include} ./static_files/soxs_stare.inc
:::



## Parameters

:::{include} parameters/soxs_stare.inc
:::

## Method

The algorithm used in the `soxs_stare` and `soxs_stare_std` recipes is shown in {numref}`soxs_stare_diagram`.


:::{figure-md} soxs_stare_diagram
![](soxs_stare.png){width=600px}

The `soxs_stare`/`soxs_stare_std` recipe algorithm. At the top of the diagram, NIR input data is found on the right and VIS on the left. 
:::

If more than one stare mode frame is passed to the `soxs_stare` recipe, there is a call to [`clip_and_stack`](../utils/clip_and_stack.md) to combine the data into a single frame. The single stare-mode frame is detrended using the [`detrend`](../utils/detrend.md), optionally dividing by a master flat field and fitting and removing the background scattered light. The sky-flux is modelled and removed using the `subtract_sky` util, and finally, the object is optimally extracted using the [`horne_extraction`](../utils/horne_extraction.md) utility.

Before the sky is modelled, SOXS stare frames get a slit-illumination check. Spatial flexure can shift the sky image along the slit, relative to the master flat and the 2D dispersion map, by an amount that differs from order to order. For each order, the recipe measures this offset from the skylines (`measure_slit_illumination`), finds the slit range that is lit in both the sky image and the master flat, and leaves the pixels outside that range out of the sky model. The offsets are written as `SLIT ILLUM OFFSET O<N>` QCs, with a QC plot, but they have no acceptable range, so flexure never fails the recipe. The step runs only if `slit-illumination: measure-offset` is `True` (the SOXS default; X-Shooter and a missing `slit-illumination` block leave it off). The offset follows the rotator angle, so the recipe logs `ABSROT` and warns if the stacked input frames differ in `ABSROT` by more than 10 degrees. The SOXS VIS u and g orders have no skylines below 557.7 nm, so their offset cannot be measured and they are not masked.

Note a boxcar extraction is also preformed alongside the Horne extraction. Use the `horne-extraction-slit-length` to control the size of the aperture used to preform the boxcar extraction on the object trace.

## Output

:::{include} output/soxs_stare.inc
:::


## QC Metrics


:::{include} qcs/soxs_stare.inc
:::

:::{figure-md} soxs_stare_qc
![image-20260414115256301](../_images/image-20260414115256301.png)

A QC plot resulting from the `soxs_stare` recipe. This is a SOXS NIR wavelength and flux calibrated spectrum of the standard star CD-325613. The top- and middle-panels show the flux and wavelength calibrated spectrum, the top in linear-flux and the middle in log-flux scale. The bottom panel shows the signal-to-noise ratio across the entire wavelength range covered by the spectrum.

:::



## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_stare.soxs_stare
:::
