# soxspipe



[![](https://zenodo.org/badge/DOI/10.5281/zenodo.8038264.svg)](https://zenodo.org/doi/10.5281/zenodo.8038264) 




<!-- INFO BADGES -->  

[![](https://img.shields.io/pypi/pyversions/soxspipe)](https://pypi.org/project/soxspipe/)
[![](https://img.shields.io/pypi/v/soxspipe)](https://pypi.org/project/soxspipe/)
[![](https://img.shields.io/conda/vn/conda-forge/soxspipe)](https://anaconda.org/conda-forge/soxspipe)
[![](https://static.pepy.tech/badge/soxspipe/month)](https://pepy.tech/project/soxspipe)
[![](https://img.shields.io/github/license/thespacedoctor/soxspipe)](https://github.com/thespacedoctor/soxspipe)

<!-- STATUS BADGES -->  

[![](https://soxs-eso-data.org/ci/buildStatus/icon?job=soxspipe%2Fmaster&subject=build%20master)](https://soxs-eso-data.org/ci/blue/organizations/jenkins/soxspipe/activity?branch=master)
[![](https://soxs-eso-data.org/ci/buildStatus/icon?job=soxspipe%2Fdevelop&subject=build%20dev)](https://soxs-eso-data.org/ci/blue/organizations/jenkins/soxspipe/activity?branch=develop)
[![](https://cdn.jsdelivr.net/gh/thespacedoctor/soxspipe@master/coverage.svg)](https://raw.githack.com/thespacedoctor/soxspipe/master/htmlcov/index.html)
[![](https://readthedocs.org/projects/soxspipe/badge/?version=master)](https://soxspipe.readthedocs.io/en/master/)
[![](https://img.shields.io/github/issues/thespacedoctor/soxspipe/type:%20bug?label=bug%20issues)](https://github.com/thespacedoctor/soxspipe/issues?q=is%3Aissue+is%3Aopen+label%3A%22type%3A+bug%22+)

*The data-reduction pipeline for the SOXS instrument* (a python package with command-line tools).

Documentation for soxspipe is hosted by [Read the Docs](https://soxspipe.readthedocs.io/en/master/) ([development version](https://soxspipe.readthedocs.io/en/develop/) and [master version](https://soxspipe.readthedocs.io/en/master/)). The code lives on [github](https://github.com/thespacedoctor/soxspipe). Please report any issues you find [here](https://github.com/thespacedoctor/soxspipe/issues).

## Installation

The best way to install or upgrade soxspipe is to use `conda` to install the package in its own isolated environment, as shown here:

``` bash
conda create -n soxspipe python=3.11 soxspipe -c conda-forge
conda activate soxspipe
```

If you have previously installed soxspipe, a warning will be issued stating that a conda environment already exists; select 'y' when asked to remove the existing environment.

To check installation was successful run `soxspipe -v`. This should return the version number of the install.

## How to cite soxspipe

If you use `soxspipe` in your work, please cite using the following BibTeX entry: 

```bibtex
@software{Young_soxspipe,
    author = {Young, David R. & Landoni, Marco},
    doi = {10.5281/zenodo.8038264},
    license = {GPL-3.0-only},
    title = {{soxspipe. The SOXS data-reduction pipeline}},
    url = {https://zenodo.org/doi/10.5281/zenodo.8038264}
}
```

