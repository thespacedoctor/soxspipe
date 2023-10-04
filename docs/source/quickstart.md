# Quickstart Guide

```eval_rst
.. warning::
    This quickstart guide is subject to (much) change during the development of the pipeline. New features and ways of operating the pipeline are still being added. Current data taken in stare mode can be reduced to the point of sky subtraction. 
```

## Install

The best way to install soxspipe is to use `conda` and install the package in its own isolated environment (using either [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [Minicoda](https://docs.conda.io/en/latest/miniconda.html)), as shown here:

``` bash
conda create -n soxspipe python=3.9 soxspipe -c conda-forge
conda activate soxspipe
```

To check installation was successful run `soxspipe -v`. This should return the version number of the installation.

With each new release of soxspipe, you can upgrade to the latest version using the command:

``` bash
conda upgrade soxspipe -c conda-forge
```

## Demo Data

The demo XShooter data (stare-mode) is of the X-ray binary SAX J1808.4-3658 taken during a 2019 outburst. You can download and unpack the data with the following commands:

```bash
curl -L "https://www.dropbox.com/s/t3adwc86bcwonkj/soxspipe-quickstart-demo-lite.tgz?dl=1" > soxspipe-quickstart-demo.tgz
tar -xzvf soxspipe-quickstart-demo.tgz
```

You may also retrieve the raw data directly from the [ESO archive](http://archive.eso.org/eso/eso_archive_main.html) with the following parameters:

```text
RA = 18 08 27.54
Dec = -36 58 44.3
Night = 2019 08 30
Spectroscopy = XSHOOTER/VLT
Science
```

## Preparing the Data-Reduction Workspace

Now you have a sample data set to work with, it is time to prepare the `soxspipe-quickstart-demo` workspace. Change into the `soxspipe-quickstart-demo` directory and run the `soxspipe prep` command:

```bash
cd soxspipe-quickstart-demo
soxspipe prep .
```

Once the workspace has been prepared, you should find it contains the following files and folders:

   - `raw_frames/`: all raw-frames to be reduced
   - `misc/`: an archive of other files that may have been found at the root of the workspace when running the `prep` command
   - `sof/`: the set-of-files (sof) files required for each reduction step
   - `soxspipe.db`: a sqlite database needed by the data-organiser, please do not delete
   - `_reduce_all.sh`: a single script to reduce all the data in the workspace

## Reduce the Data

All of the `soxspipe` recipe commands needed to reduce workspace data can be found in the `_reduce_all.sh` script. To reduce the data, simple execute that script:

```bash
sh _reduce_all.sh
```
