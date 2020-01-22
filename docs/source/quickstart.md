# Quickstart

This quickstart guide assumes you have already installed `soxspipe` in a conda environment on your machine; if you haven't please go ahead and follow the [install instructions](index.html#installation).

## Downloading the SOXS Calibration Data and Demo Data-Suite

Before we begin you will need to download the SOXS calibration data and also some raw-data to reduce.

{{_includes/downloading_data_suite.md}}

## Initialising soxspipe

Before you begin using soxspipe you need to use the `init` command to generate a user settings file. Running the following:

```bash
> soxspipe init

default settings have been added to '~/.config/soxspipe/soxspipe.yaml'. Tailor these settings before proceeding to run soxspipe
```

This command creates a yaml settings file in your home foler under `~/.config/soxspipe/soxspipe.yaml`. This is where most of the `soxspipe` setting can be adjusted to liking.

<!-- Once created, open the settings file in any text editor and follow the in-file instructions to populate the missing settings values (usually given an ``XXX`` placeholder).  -->



## Command-Line Usage

{{usage.md}}

