# Recipes

SOXSPIPE borrows the informative concept of `recipes' employed by ESO's data reduction pipelines to define the modular components of the pipeline. These recipes can be strung together to create an end-to-end workflow that takes as input the raw and calibration frames from the instrument and telescope and processes them all the way through to fully reduced, calibrated, ESO Phase III compliant science products.

## Standard Calibrations

```eval_rst
.. toctree::
   :maxdepth: 1

   soxs_mbias
   soxs_mdark
   soxs_mflat
```

## Dispersion and Spatial Solutions

There is a strong curvature in the traces of the NIR orders and spectral-lines do not run perpendicular to the dispersion direction, but are highly tilted. Therefore wavelength cannot be expressed as simply a function of pixel position, but instead detector pixel positions ($X, Y$) much be mapped as a function of:

1. wavelength $\lambda$
2. order number $n$, and 
3. slit position $s$

This 2D mapping function is determined incrementally via the `soxs_disp_solution`, `soxs_order_centres` and `soxs_spatial_solution` recipes. The `soxs_straighten` recipe can then be used to transform spectral images from detector pixel-space to wavelength and slit-position space 

```eval_rst
.. toctree::
   :maxdepth: 1

   soxs_disp_solution
   soxs_order_centres
   soxs_spatial_solution
   soxs_straighten
```
