## `soxs_straighten ` - PLANNED

This recipe takes the full dispersion map given by `soxs_spatial_solution` and uses it to map images from their original representation of the detector surface to one that presents the signal in a wavelength by slit-position coordinate system.

### Input

| Data Type | Content | Related OB |
|:----|:----|:---|
| CSV File | Coefficients of polynomials providing a full dispersion-spatial solution | |
| FITS Image | An associated spectral image requiring rectification | Many |

### Parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| `straighten_grid_res_wavelength`  |  size of the grid cell in wavelength dimension (nm)  |  float  |  settings file  | | 
| `straighten_grid_res_split` |  size of the grid cell in slit dimension (arcsec)  |  float  |  settings file  | | 

### Method

We now have a pair of polynomials that can be used to give the exact pixel on the detector containing flux resulting from a specific order, with a given wavelength and slit position.

$$
X = \sum\limits_{ijk} c_{ijk} \times n^i \times \lambda^j \times s^k \\
$$

$$
Y = \sum\limits_{ijk} c_{ijk} \times n^i \times \lambda^j \times s^k \\
$$

To begin we want to create a full wavelength and slit-position map; a 2D grid of wavelengths along one axis and slit-position along the other. Using the polynomial solutions above, we populate each cell in the grid with its corresponding detector pixel coordinate. The 2D grid is of fine enough resolution so that many cells in the grid are mapped to each individual detector pixel. With this map in hand we can now assign flux recorded in each detector pixel to the corresponding cells in the 2D wavelength and slit-position grid. The flux from each detector is evenly distributed between all cells found to be associate with that pixel; so if 9 cells are associated then each cell gets 1/9th of the pixel flux.

The error and bad-pixel extensions go through the same mapping process.

![](soxs_straighten.png)

### Output
 
| Data Type | Content |
|:----|:----|
| FITS Images | The straightened images containing flux represented in wavelength space; one for each order | 

### QC Metrics

| Metric  | Description |
| :------------ | :----------- |
| TBC     | ...  |

### Recipe API

```eval_rst
.. autoclass:: soxspipe.recipes. soxs_straighten
    :members:
```
