:::{table} The `soxs_disp_solution` recipe parameters.
:name: soxs_disp_solution_parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| pixel-window-size | the size of the square window used to search for an arc-lamp emission line, centred on the predicted pixel position of the line | int | settings file | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| order-deg | degree of echelle order number component of global polynomial fit to the dispersion solution | int | settings file or cmd-line  | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| wavelength-deg | degree of wavelength component of global polynomial fit to the dispersion solution | int | settings file or cmd-line | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| poly-clipping-iteration-limit | number of sigma-clipping iterations to perform before settling on a polynomial fit for the dispersion solution | int | settings file | [`create_dispersion_map`](../utils/create_dispersion_map.md) |
| poly-fitting-residual-clipping-sigma | sigma clipping limit when fitting global polynomial to the dispersion solution | float | settings file | [`create_dispersion_map`](../utils/create_dispersion_map.md) |

:::



