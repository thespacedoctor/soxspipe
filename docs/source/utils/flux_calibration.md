# `flux_calibration` 

The purpose of the [`flux_calibration`](#soxspipe.commonutils.flux_calibration) utility is to flux calibrate a scientific spectrum applying the response function computed by [`response_function`](#soxspipe.commonutils.response_function) utility.


## Input

| Frame.                   | Description                                   | 
| ------------------------ | --------------------------------------------- |
| Extracted  1D scientific spectrum | FITS table containing the 1D spectrum (any observing mode) of the object to be flux calibrated |  
| Extinction curve | FITS table containing the tabulated value of the extinction curve (wavelenght, mag per airmass) of the observing site.|
| Response function | FITS table containing the fit parameters modelling the response function calculed by [`response_function`](#soxspipe.commonutils.response_function) utility |

## Parameters

N/A

## Method
The general algorithm and steps performed by [`flux_calibration`](#soxspipe.commonutils.flux_calibration) are the one reported in the flow chart below:

![](flux_calibration.png).

In details, 


## Output


| Data Type | Content |
| ------------------------ | --------------------------------------------- |
|FITS table |FITS table containing the flux calibrated input spectrum|