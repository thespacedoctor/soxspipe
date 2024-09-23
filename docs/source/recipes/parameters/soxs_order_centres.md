:::{table} The `soxs_order_centres` recipe parameters.
:name: soxs_order_centres_parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| order-sample-count  | number of times along the order in the dispersion direction to measure the order-centre trace  |  int | settings file |  [detect_continuum utility](../utils/detect_continuum.md) |
| peak-sigma-limit  |  minimum value a peak must be above the median value of pixel to be considered for order-trace fitting  | int | settings file  |  [detect_continuum utility](../utils/detect_continuum.md) |
| disp-axis-deg | degree of dispersion axis component of polynomal fit to order-centre traces |   | int | settings file  |  [detect_continuum utility](../utils/detect_continuum.md) |
| order-deg | degree of order component of polynomal fit to order-centre traces |   | int | settings file  |  [detect_continuum utility](../utils/detect_continuum.md) |
| poly-fitting-residual-clipping-sigma  | sigma distance limit, where distance is the difference between the detected and polynomial fitted positions of the order-trace, outside of which to remove lines from the fit   | float   | settings file |  [detect_continuum utility](../utils/detect_continuum.md) | 
|  poly-clipping-iteration-limit  |  number of sigma-clipping iterations to perform before settings on a polynomial fit for the order-centre traces  |  int   | settings file | [detect_continuum utility](../utils/detect_continuum.md) |


:::
