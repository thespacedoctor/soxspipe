# Installation

The best way to install soxspipe is to use `conda` and install the package in its own isolated environment, as shown here:

``` bash
conda create -n soxspipe python=3.11 soxspipe -c conda-forge
conda activate soxspipe
```

To check installation was successful run `soxspipe -v`. This should return the version number of the install.

To upgrade to the latest version of soxspipe use the command:

``` bash
conda upgrade soxspipe -c conda-forge
```
