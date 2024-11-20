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

The demo XShooter data is of AT2020xnd, a super-luminous supernova taken in stare mode, and HD 168076, a bright O-type star taken in nodding mode. You can download and unpack the data with the following commands:

```bash
curl -L "https://www.dropbox.com/scl/fi/g7ie2i4ijh0w2xrjq67xo/soxspipe-quickstart-demo-lite.tgz?rlkey=eow6ujhyyo1drmzv2yt2qpo8i&dl=1" > soxspipe-quickstart-demo.tgz
tar -xzvf soxspipe-quickstart-demo.tgz
```

These are exactly the same data you can download from the [ESO archive](http://archive.eso.org/eso/eso_archive_main.html).

## Preparing the Data-Reduction Workspace

Now you have a sample data set to work with, it is time to prepare the `soxspipe-quickstart-demo` workspace. Change into the `soxspipe-quickstart-demo` directory and run the `soxspipe prep` command:

```bash
cd soxspipe-quickstart-demo
soxspipe prep .
```

Once the workspace has been prepared, you should find it contains the following files and folders:

- `misc/`: a lost-and-found archive of non-fits files
- `raw/`: all raw-frames to be reduced
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
