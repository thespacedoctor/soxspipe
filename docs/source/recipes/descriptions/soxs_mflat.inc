The [`soxs_mflat`](#soxspipe.recipes.soxs_mflat) recipe creates a single normalised master-flat frame used to correct for non-uniformity in response to light across the detector plane. Hot and dead pixels are also detected and added to a bad-pixel mask. Finally, the echelle order edges are detected and fitted with a polynomial model.

Sources of this non-uniformity include:
Varying pixel sensitivities.
Obstructions in the optical path (e.g., dust or pollen grains).
Vignetting at the edges of the detector.
A flat frame is ideally an image taken with uniform illumination across the detector's light-collecting pixels. This evenly exposed image can be used to identify irregularities in the detector's response.
