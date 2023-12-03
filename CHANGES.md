
# Release Notes  

## v0.9.2 - November 29, 2023

* **ENHANCEMENT:** intra-order background (scattered light) fits are now being written to FITS image files in the QC directories and reported at the end of a recipe run.
* **ENHANCEMENT:** added a `create_dispersion_solution_grid_lines_for_plot` function to allow adding dispersion solution grid to QC plots.  This is extremely useful for quickly diagnosing problems with the fits.
* **REFACTOR:** All product FITS files now pass fitverify without error or warnings. All issues were due to using '-' instead of underscores in FITS binary table column names.
* **REFACTOR:** bad-pixel values set to 0 in data extensions of products
* **REFACTOR:** nans have been replaced by zero in FITS image product
* **FIXED**: a mismatch between daofind results and the original input pixel table was causing dispersion solution to break (a recent bug introduced during code optimisations)
* **FIXED**: the internal soxspipe logger was being interfered with by astropy so that logs were somtimes getting redirected to the wrong place

## v0.9.0 - October 11, 2023

* **FEATURE:** added a `predict_product_path` function to determine the product path from a recipe's input sof file
* **FEATURE:** Merging of individual order extracted spectra from object frame into a single spectrum for each arm
* **FEATURE:** Object spectra are now extracted from the sky-subtracted frames using the Horne 86 method
* **FEATURE:** Real SOXS data is now included in the unit-test suite (starting to replace simulated data unit-tests). `soxs-disp-solu` recipe so far.
* **FEATURE:** SOXS NIR Xe line-lists added to static-calibration suite (single and multi pinhole).
* **FEATURE:** when running a recipe, `soxspipe` writes informative logs to stdoutAND to a log file adjacent to the recipe's product file(s). Error logs are also written if a recipe fails (see docs).
* **ENHANCEMENT:** recipe timing added to the end of the logs
* **ENHANCEMENT:** fitted lines from the dispersion solution are written out to file as a QC product
* **ENHANCEMENT:** flux (and other daostarfinder metrics) are now recorded in the detected line-list QC file. This will help measure degradation of arc-lamps over time.
* **ENHANCEMENT:** FWHM and pixel-scale added to fitted lines from the dispersion solution
* **ENHANCEMENT:** legends added to many of the QC plots
* **ENHANCEMENT:** OB ids now getting add to the data-organiser database tables.
* **ENHANCEMENT:** object trace FITS binary table added to stare-mode products (alongside complimentary QC plot)
* **ENHANCEMENT:** products and QC outputs are differentiated in the table reported upon recipe completion (see label column).
* **ENHANCEMENT:** verifying the master flat used to calibrate object/std spectra has the same slit-witdh as used to take the science frames
* **REFACTOR:** `init` command has been subsumed into the prep command. The `prep` command will generate a settings file to live within the prepared workspace.
* **REFACTOR:** `misc/` directory created by data-organiser even if empty
* **REFACTOR:** close matplotlib plot after writing plots to file
* **REFACTOR:** command-line startup speeds improved
* **REFACTOR:** continuum fitting code made more robust against edge cases (orders of the fit are automatically reduced if fit does not converge)
* **REFACTOR:** soxspipe now has a 'full' and a 'lite' test-suite. Using the lite suite will speed up deploying of new releases.
* **DOCS:** updated docs with a more robust SOXSPIPE upgrade path (users having issue with `conda update ...`)
* **FIXED**: sky-subtraction code and data-organiser fixed to work with binned data


## v0.8.0 - May 18, 2023

* **FEATURE:** we now have a data-organiser to sort data, prepare the required SOF files and generate reduction scripts.
* **ENHANCEMENT:** '.db', '.yaml', '.sh' and '.log' extensions skipped when moving items to the misc folder
* **ENHANCEMENT:** move information printed to STDOUT when preparing a workspace to inform the user of how the data is organised
* **ENHANCEMENT:** code can automatically adjust polynomial fitting parameters to find a dispersion solution if those provided in the settings file fail.
* **ENHANCEMENT:** uncompression of fits.Z files (if any) occurs before data-organising
* **REFACTOR:** speed & robustness improvements to dispersion solution to 2D image map conversion.
* **REFACTOR:** much fast check for product existence so recipes are quickly skipped if they have already run.
* **REFACTOR:** removed the `intermediate-data-root` setting renamed to a more accurate `workspace-root-dir`
* **REFACTOR:** removed the `reduced-data-root` setting.
* **REFACTOR:** updating all depreciated pandas commands so pipeline is now compatible with 1.X and 2.X versions of pandas 
* **FIXED** pandas 1.X and pandas 2.X were doing different things when renaming columns in data-frames. Both 1.X and 2.X now work in the pipeline.

## v0.7.2 - March 3, 2023

* **REFACTOR:** Big improvements on sky-subtraction  
* **REFACTOR:** UV order-edge detection more robust  
* **REFACTOR:** changed quickstart guide compress to gzipped tar  
* **REFACTOR:** updated default settings to be more robust  

## v0.7.1 - November 4, 2022  

* **FEATURE:** UV D-Lamp and QTH-Lamp master flats now being stitched together  
* **FEATURE:** errors in error maps now being treated correctly and propagating to combined images  
* **FEATURE:** Pipeline can now 'remember' where it left off in the reduction cascade. If it has run a recipe before it will exit with a message to the user informing them how to force the recipe to rerun.  
* **FEATURE:** added a `twoD_disp_map_image_to_dataframe` function to toolkit  
* **ENHANCEMENT:** PRO CATG now written to product FITS header  
* **ENHANCEMENT:** Handling of binned images when generating flats and order-locations  
* **ENHANCEMENT:** Where possible, product files are given the same name as the SOF file used to generate them (replacing `.sof` extension with `.fits`)  
* **ENHANCEMENT:** SOF files can now contain a file 'tag' to allow users to read the SOF file contents and know exactly which files are being passed to the recipe (e.g. `MASTER_BIAS_UVB`, `LAMP,DORDERDEF_UVB` ... )  
* **ENHANCEMENT:** dispersion solution now working with simulated NIR SOXS data  
* **ENHANCEMENT:** quicklook now renders dispersion solution grid  
* **ENHANCEMENT:** \~40% speed gain in combining images.  
* **REFACTOR:** 2D Map generation now ~6-8 times faster (seeding solutions with nearest neighbour with cubic spline method)  
* **REFACTOR:** SOF filenames reworked to contain the UTC observation date instead of MJD (more in-line with ESO ecosystems)  
* **REFACTOR:** updated workflow for master bias combination  
* **REFACTOR:** updated workflow for master dark combination  
* **REFACTOR:** QC PDF plots now added to their own directory separate from the products    
* **REFACTOR:** products now sub-divided into recipe directories (e.g. `./products/soxs-mbias/`)    
* **DOCS:** mflat docs brought up-to-date    
* **DOCS:** mflat docs brought up-to-date    
* **FIXED:** mflat recipe now exits if flat frames are not of a consistent exptime.    


## v0.6.2 - April 13, 2022

* **ENHANCEMENT:** quickstart guide added for calibration recipes  
* **FEATURE:** QCs added for dispersion solution and order centre recipes  
* **REFACTOR:** clean up of stdout information  

## v0.6.1 - April 11, 2022

* **FEATURE:** shipping static calibration files with the code (one less thing for end-users to install and set-up)

## v0.6.0 - April 10, 2022

This is only a summary of some of the updates included in this release:

* **ENHANCEMENT:** All CSV files moved to FITS binary tables - metadata very useful for developing data organiser
* **FEATURE:** 2D image map now created by create_dispersion_solution
`subtract_calibrations` util renamed to `detrend` and added ability to flat correct
* **FEATURE:** 2D image map of wavelength values, slit-position values and order values written alongside polynomial solutions of full dispersion solution
* **FEATURE:** soxspipe now on conda
* **FEATURE:** QCs now being written to FITS header
* **FEATURE:** adding QC and product collection in mbias recipe
* **ENHANCEMENT** RON and bias structure QCs now reported by mbias
* **ENHANCEMENT** nan ignored when scaling quicklook images
* **ENHANCEMENT** RON and bias structure QCs now reported by mdark
* **ENHANCEMENT:** QCs have an option to *NOT* (`to_header`) write to FITS header (default is to write)
* **REFACTOR:** better treatment of masked pixels when stacking images (e.g. in mbias and mdark)
* **REFACTOR:** removed raw frame reports and neater QC table
* **REFACTOR:** fits header keywords neatly sorted before writing to file
* **FIX:** Correct management of mask when determining RON on bias and darks

## v0.5.1 - September 29, 2021

* **FEATURE:** recipes now have a `qc` and `products` attribute. These are pandas data frames used to collect QCs and generated products throughout the life-time of the recipe. They are printed to STDOUT at the end of the recipe (can be used in the future to send post request to health monitor API with JSON content in request body).
* **ENHANCEMENT** added code-base to conda-forge
* **ENHANCEMENT** added bottleneck to the install requirement (makes image combination more efficient)
* **ENHANCEMENT** masked pixel now coloured red in quicklook plots (easier to differentiate from good pixels)
* **ENHANCEMENT** low-sensitivity pixels in lamp-flats now identified and added to bad-pixel mask
* **ENHANCEMENT** add a verbosity flag to the command-line and a verbose parameter to each recipe
* **REFACTOR** inter-order pixel value in flats now set to unity (instead of running background fitting and subtraction)
* **REFACTOR:** recipes now have their recipe name as a `recipeName` attribute

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
* **FEATURE:** added a `detrend` method to subtract calibration frames (bias and dark) from an input frame
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
