# Recipes

## Standard Calibrations

```eval_rst
.. toctree::
   :maxdepth: 1

   soxs_mbias
   soxs_mdark
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
