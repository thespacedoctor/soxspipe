# Release Notes

## v0.5.0 - June 10, 2021

* **FEATURE** Added a new `filenamer` module that implements a strict intermediate and reduced file-naming scheme
* **FEATURE:** `soxs_mflat` recipe now included
* **FEATURE:** `soxs_spatial_solution` recipe is now included
* **FEATURE:** `subtract_background` utility added
* **FEATURE:** added a `detect_order_edges` object
* **FEATURE:** Added a `dispersion_map_to_pixel_arrays` function to convert from order-based and wavelength arrays to pixel arrays (first guess dispersion map only so far)
* **FEATURE:** added a quicklook function in toolkit to quickly visualise a frame
* **FEATURE:** added a toolkit module for small functions used throughout soxspipe 
* **FEATURE:** added function in toolkit to unpack an order table into lists of coordinates, one list per order.
* **FEATURE:** added image slice tool to toolkit
* **ENHANCEMENT** Added a `-o <outputDirectory>` switch to the command-line to optionally override the 'intermediate-data-root' setting in the settings file.
* **ENHANCEMENT:** added a fraction of a second tolerance when matching exptimes between darks and science/calibration frames 
* **ENHANCEMENT:** y limits now added to the order table to show limits of order locations on detector
* **REFACTOR:** Change the "SOXSPIPE PRE" date stamp keyword to "SXSPRE" to future-proof for phase III (8 character keyword limit)
* **REFACTOR:** Pandas tables are now used through-out code to pass line-lists between methods
* **REFACTOR:** refactoring of polynomial fitting has made creation of dispersion maps ~50 times faster
* **REFACTOR:** removed OBID from file names and added readout mode. This information is more helpful at the glance.
* **FIXED:** correct binning reported in product file names
* **FIXED:** lines in a sof file beginning with a `#` are considered as comments and therefore ignored by the pipeline.

## v0.4.1 - September 15, 2020

* **FEATURE:** add command-line util for soxs order_centres recipe
* **FEATURE** added the `detect_continuum` utility to fit order centre locations in single pinhole flat frames.
* **ENHANCEMENT:** added a supplementary file list for non-fits input files in set-of-file util
* **ENHANCEMENT:** adding more information residual plots & visualisation of fitting for disp solution
* **ENHANCEMENT:** check that files in the sof files exist before proceeding.
* **ENHANCEMENT:** added spectral format table lookup to detector settings file
* **REFACTOR:** moved chebyshev order/wavelength polynomials into its own class - decoupled from create_dispersion_map class

## v0.4.0 - September 3, 2020

* **FEATURE:** added create_dispersion_map class to be used in `soxs_disp_solution` and `soxs_spatial_solution`
* **FEATURE:** added a `subtract_calibrations` method to subtract calibration frames (bias and dark) from an input frame
* **FEATURE:** added the dispersion solution recipe and unit tests
* **FEATURE:** added the disp_solution command-line tool
* **DOCS:** major docs overhaul
* **ENHANCEMENT:** added predicted lines lists to detector parameter file
* **ENHANCEMENT:** DPR CATG and DPR TECH added to metadata of sof imagefilecollection objects
* **ENHANCEMENT:** wcs copied from a single frame into the combined frames during clip and stack
* **REFACTOR:** bad-pixel map paths abstracted to detector settings files
* **REFACTOR:** renaming of unit-testing test data directories
* **REFACTOR:** only filenames reported by sof summaries when files are found in the same directory (easier to read on terminal) 
* **FIXED:** fixed detector science pixels for UVB

## v0.3.1 - August 25, 2020

* **FEATURE:** recipe & command-line tool for master dark creation (`mdark`)
* **ENHANCEMENT:** default binning add to detector settings file
* **ENHANCEMENT:** added mixed exposure time unit test for dark-frames
* **ENHANCEMENT:** added default values for gain and ron in the detector settings files. Default values can be overwritten if correct GAIN and RON are found in fits-headers (overwritten for UVB and VIS but not NIR for XShooter)
* **ENHANCEMENT:** can now interrupt "~" as home directory in sof file path
* **FIXED:** binning factor used when trimming frames
* **FIXED:** the trimming dimensions of NIR frames - bad-pixel map now aligns correctly with data frame
* **FIXED:** science pixels for all 3 xshooter detectors in parameters file

## v0.3.0 - August 18, 2020

* **FEATURE:** added a `write` method to the `_base_recipe` to write frames to disk (renames extensions to ESO preferred naming scheme)
* **FEATURE:** detector lookup class added alongside yaml files to host detector specific parameters (rotation, science-pixels etc). Code has been updated to remove hard-wired detector values.
* **FEATURE:** added a cleanup method to remove intermediate file once recipe completes
* **ENHANCEMENT:** parameters for clip and stack method added to the settings files
* **ENHANCEMENT:** added strict typing of data and variables with astropy units to avoid silent mistakes in frame arithmetic 
* **ENHANCEMENT:** added mixing of readout speeds to input frame verification checks
* **ENHANCEMENT:** added readnoise and gain to the list of keyword values to check during frame verification
* **ENHANCEMENT:** inject a 'SOXSPIPE PRE' keyword with timestamp value into prepared frames
* **ENHANCEMENT:** check frames for 'SOXSPIPE PRE' keyword before preparing - raises exception if found
* **ENHANCEMENT:** ron and gain are added to the recipe's detector lookup dictionary during frame verification (so they don't need read again later)
* **REFACTOR:** moved stacking code to it own `clip_and_stack` method hosted in the `_base_recipe`
* **REFACTOR:** moved basic input frame verifications to the `_base_recipe` - so not to repeat code
* **REFACTOR:** removed python 2.7 support - not feasible with CCDProc
* **DOCS:** added workflow diagrams to the documentation for many of the methods implemented (`prepare_frames()`, ``clip_and_stack``)

## v0.2.0 - February 27, 2020

* **FEATURE** added keyword lookups - abstracting exact keyword names from code
