:::{table} The `soxs_disp_solution` recipe parameters.
:name: soxs_disp_solution_parameters
:widths: 25, 30, 10, 15, 20


| Parameter | Description | Type | Entry Point | Related Util |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| `pixel_window_size` | the size of the square window used to search for an arc-lamp emission line, centred on the predicted pixel position of the line | int | settings file | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| `pinhole_detection_thres_sigma` | minimum significance required for arc-line to be considered 'detected' | float  |settings file                | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| `order_deg` | degree of echelle order number component of global polynomial fit to the dispersion solution [x, y] | int/list | settings file or command-line | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| `wavelength_deg` | degree of wavelength component of global polynomial fit to the dispersion solution [x, y] | int/list | settings file or command-line | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| `poly_clipping_iteration_limit` | number of sigma-clipping iterations to perform before settling on a polynomial fit for the dispersion solution | int | settings file | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| `poly_fitting_residual_clipping_sigma` | sigma clipping limit when fitting global polynomial to the dispersion solution | float | settings file | [`create_dispersion_map`](../utils/create_dispersion_map.md) |

:::



