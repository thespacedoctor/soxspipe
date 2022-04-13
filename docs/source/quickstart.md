# Quickstart Guide

```eval_rst
.. warning::
    This quickstart guide is subject to (much) change during development of the pipeline. New features and ways of operating the pipeline are still being added.
```

## INSTALL

The best way to install soxspipe is to use `conda` and install the package in its own isolated environment (using either [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [Minicoda](https://docs.conda.io/en/latest/miniconda.html)), as shown here:

``` bash
conda create -n soxspipe python=3.8 soxspipe -c conda-forge
conda activate soxspipe
```

To check installation was successful run `soxspipe -v`. This should return the version number of the install.

To upgrade to the latest version of soxspipe use the command:

``` bash
conda upgrade soxspipe
```

## Initialising SOXSPIPE for the first time

If this is the first time you have installed soxspipe on the machine you are using, you must first initialise it by running:

```bash
soxspipe init
```

This adds a default soxspipe settings file and log file within `~/.config/soxspipe/`. The default parameters that are set within this settings file are good for typical reductions, but as a user, you can adjust these settings to your needs. 

```eval_rst
.. note::
    It is also possible to copy this settings file, tailor it to your needs and use the amended settings file with the `--settings` option. Run `soxspipe -h` to see how.
```

## Demo Data

```eval_rst
.. note::
    All of the Set-of-Files (`sof`) files have been pre-generated for this demo and in a future release of the pipeline these files will probably be completely obfuscated from the typical user. Instead, a built-in data organiser will intelligently select the correct files to process.
```


```eval_rst
.. note::
    Although demo data here is XShooter data, the pipeline can also run on simulated SOXS data (at least for dark frames)
```

The demo XShooter data (stare-mode) is of the X-ray binary SAX J1808.4-3658 taken during a 2019 outburst. You can download and unpack the data with the following commands:

```bash
curl -L "https://www.dropbox.com/sh/6zlxu80rfy5svjc/AAD0ze4r2pdtKTOWmmac6r_Ca?dl=1" > soxspipe-quickstart-demo.zip
unzip soxspipe-quickstart-demo.zip -d soxspipe-quickstart-demo
```

You may also retrieve the raw data directly from the [ESO archive](http://archive.eso.org/eso/eso_archive_main.html) with the following parameters:

```text
RA = 18 08 27.54
Dec = -36 58 44.3
Night = 2019 08 30
Spectroscopy = XSHOOTER/VLT
Science
```

To help you get familiar with running the pipeline, one or two example commands will be shown for each pipeline recipe. Within the demo-data package, there is also a `run_pipeline.sh` script that can be used to a) view all the commands needed to reduce these data and b) to reduce the in one go by running `sh run_pipeline.sh` within your conda environment.


## MBIAS

Change into the `soxspipe-quickstart-demo` direwctory and create a master bias frames with a command like:

```bash
cd soxspipe-quickstart-demo
soxspipe  mbias sof/58725_VIS_BIAS_1x1_slow.sof -o ./
```

Upon completion, you should be presented with some quality control (QC) metrics for the data and a list of products produced by the recipe. QCs are also written to the FITS headers of the products.

[![](https://live.staticflickr.com/65535/51999455194_dede3217a4_b.jpg)](https://live.staticflickr.com/65535/51999455194_dede3217a4_b.jpg)

```eval_rst
.. note::
    QC checks have not been added to all recipes yet
```

## MDARK

Create the master dark frames with commands like:

```bash
soxspipe  mdark sof/58724_NIR_DARK_300pt0.sof -o ./
```

## DISP_SOL

We will now create the first guess dispersion solution using the single pinhole frames. Here's an example command to create the initial dispersion solution:

```bash
soxspipe disp_sol ./sof/58725_NIR_LAMP,FMTCHK_10.sof -o ./
```

Here, for the first time, you will find a PDF file included in the products table:

[![](https://live.staticflickr.com/65535/51999630094_f97cb55f7f_b.jpg)](https://live.staticflickr.com/65535/51999630094_f97cb55f7f_b.jpg)

The PDF includes a visualisation of the dispersion map fit and its associated residuals.

[![](https://live.staticflickr.com/65535/51999627639_9b1c73e26a_z.png)](https://live.staticflickr.com/65535/51999627639_9b1c73e26a_o.png)


## ORDER_CENTRES

To create an initial order table, fitting a global polynomial to the order centre traces from a single pinhole lamp-flat frame, run:

```bash
soxspipe  order_centres ./sof/58725_NIR_LAMP,ORDERDEF_1.sof -o ./
```

## MFLAT

```eval_rst
.. note::
    We still need to do a little work to stitch/combine the master flats generated from the 2 UVB flat-lamps.
```

Now create the master flats and order tables with order edges defined:

```bash
soxspipe mflat ./sof/58725_NIR_LAMP,FLAT_8pt320.sof -o ./
```

## SPAT_SOL

```eval_rst
.. note::
    This recipe takes some time to run as it generates the 2D wavelength/slit/order map image required in later recipes.
```

We will now use the multi-pinhole frames to generate the full dispersion solution.


Now run the `spat_sol` recipe like:

```bash
soxspipe spat_sol sof/58725_NIR_LAMP,WAVE_0pt6651.sof -o ./
```


... more to come!

