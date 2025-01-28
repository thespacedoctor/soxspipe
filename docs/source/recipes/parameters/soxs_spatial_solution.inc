:::{table} The `soxs_spatial_solution` recipe parameters.
:name: soxs_spatial_solution_parameters
:widths: 25, 30, 10, 15, 20

| Parameter                            | Description                                                                                                                     | Type     | Entry Point              | Related Util                                                 |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------- | -------- | ------------------------ | ------------------------------------------------------------ |
| use_flat                             | divide image by master flat frame                                                                                               | bool     | settings                 | -                                                            |
| subtract_background                  | fit and subtract the intra-order background light                                                                               | bool     | settings                 | [`subtract_background`](../utils/subtract_background.md)     |
| pixel-window-size                    | the size of the square window used to search for an arc-lamp emission line, centred on the predicted pixel position of the line | int      | settings                 | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| pinhole-detection-thres-sigma | minimum significance required for arc-line to be considered 'detected' | float  | settings                 | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| order-deg                            | degree of echelle order number component of global polynomial fit to the dispersion solution [x, y]                             | int/list | settings file or command-line | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| wavelength-deg                       | degree of wavelength component of global polynomial fit to the dispersion solution [x, y]                                       | int/list | settings file or command-line | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| slit-deg                             | degree of slit position component of global polynomial fit to the dispersion solution [x, y]                                    | int/list | settings file or command-line | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| poly-clipping-iteration-limit        | number of sigma-clipping iterations to perform before settings on a polynomial fit for the dispersion solution                  | int      | settings                 | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| poly-fitting-residual-clipping-sigma | sigma clipping limit when fitting global polynomial to the dispersion solution                                                  | float    | settings                 | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| poly-clipping-pinhole-sets           | clipping performed on multi-pinhole sets (true) or individual pinholes (false)                                                  | bool     | settings                 | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| map-to-image-displacement-threshold  | maximum distance allowed from the pixel centre when calculating wavelength, order and slit-position for 2d disp-sol image       | float    | settings                 | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| mph_line_set_min                     | full multi-pinholes sets (same arc line) with fewer than mph_line_set_min lines detected get clipped                            | int      | settings                 | [`create_dispersion_map`](../utils/create_dispersion_map.md) |

:::




