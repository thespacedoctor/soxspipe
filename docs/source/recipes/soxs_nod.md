# soxs_nod

:::{include} ./descriptions/soxs_nod.inc
:::


## Input

:::{include} ./inputs/soxs_nod.inc
:::

:::{include} ./static_files/soxs_nod.inc
:::


## Parameters


:::{include} parameters/soxs_nod.inc
:::


## Method

The algorithm used in the `soxs_nod` recipe is shown in {numref}`soxs_nod_diagram`.

:::{figure-md} soxs_nod_diagram
![](soxs_nod.png){width=600px}

The `soxs_nod` recipe reduced SOXS data acquired in nodding mode with one or multiple ABBA (jitter offsets) sequences. 
:::

During nodding observations, as described in the SOXS Observing Mode sections, 2 locations on the slit, A and B, separated by a few arcseconds, are considered. Observations are performed in pairs of equal-length exposures called cycles. First, a spectrum is taken with the target located at slit-position A. The telescope is then shifted, or 'nodded', so the target is at slit-position B and a second spectrum is taken. This process is then repeated in reverse to acquire the third and fourth exposures at slit-positions B and A, a BA cycle.

Two or more cycles are then combined into a sequence. The shortest sequence possible is AB and BA or ABBA, but sequences can be compiled by many cycles, e.g., an 8-cycle sequence ABBAABBAABBAABBA. 

Within this schema, jittering can also be considered. In this case, the target is shifted by a slight random offset each time it is placed at slit positions A and B, avoiding putting the targets always on a bad row or column of the CCD. When including jitter an ABBA sequence becomes A $^{1}$ B $^{1}$ B $^{2}$ A $^{2}$, where the subscript indicates a location at the slit-positions A or B but now including a small random offset (a jitter).

Adopting this assumption, the `soxs_nod` recipe proceeds as follows. First, it identifies how many jitter positions are in the provided frames. If only one is present (in other words, no jitter offsets are considered), the recipe detrends the images using bias, flat field and subtracts the scattered background light. Then, all the frames in the A position are stacked together. The same is also done for frames in the B position. Then, two images, $A - B$ and $B - A$, are computed where A and B are the stacked images described before. Those images have the sky background and dark current subtracted. At this stage, the standard Horne extraction utility is run, and orders are extracted from the first and second subtracted images and merged. The two spectra obtained from the $A - B$ and $B - A$ images are stacked together for a unique spectrum.


If the jitter is present, the `soxs_nod` recipe determines how many different offsets are present in the provided frames, which are grouped according to their offsets before stacking. The procedure described above is then repeated for each offset detected in the input data. At the end of the procedure, there will be a number of spectra equal to the number of detected offsets in the data. Those spectra are then stacked and merged to obtain a unique spectrum.


## Output
 
:::{include} output/soxs_nod.inc
:::

## QC Metrics


:::{include} qcs/soxs_nod.inc
:::

## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_nod.soxs_nod
:::
