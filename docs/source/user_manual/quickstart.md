# Quickstart Guide

:::{warning}
This quickstart guide is subject to (much) change during the development of the pipeline. New features and ways of operating the pipeline are still being added. Current data taken in stare mode can be reduced to the point of sky subtraction. 
:::

## Install

The best way to install soxspipe is to use `conda` and install the package in its own isolated environment (preferably using [Miniforge](https://github.com/conda-forge/miniforge); a minimal installation for the conda package and environment manager), as shown here:

``` bash
conda create -n soxspipe python=3.9 soxspipe -c conda-forge
conda activate soxspipe
```

If you have previously installed soxspipe, a warning will be issued stating that a conda environment already exists; select `y` when asked to remove the existing environment. This has proven to be the cleanest way to upgrade soxspipe.

To check installation was successful run `soxspipe -v`. This should return the version number of the installation.

For alternative methods of installing `soxspipe`, please refer to the [installation section](installation.md).

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

- `misc/`: a lost-and-found archive of non-fits files
- `raw_frames/`: all raw-frames to be reduced
- `sessions/`: directory of data-reduction sessions
- `sof/`: the set-of-files (sof) files required for each reduction step
- `soxspipe.db`: a sqlite database needed by the data-organiser, please do not delete
- `soxspipe.yaml`: file containing the default settings for each pipeline recipe

soxspipe reduces data within a [`reduction session`](./sessions.md) and an initial `base` session is automatically created when running the `prep` command. For a more detailed guide to preparing a workspace, please refer to the [Preparing a Data-Reduction Workspace](preparing_a_workspace.md) section.

## Reduce the Data

In most use case, you will want to reduce all of the raw frames contained within your workspace. To do this run the command:

```bash
soxspipe reduce all .
```

The `reduce` command stops when all data within the workspace have been reduced. For a more in-depth guide to the using the `reduce` command, please refer to the [Reducing Data](reductions/index.md) section.
