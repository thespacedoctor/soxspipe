# soxspipe

[![](https://readthedocs.org/projects/soxspipe/badge/)](http://soxspipe.readthedocs.io/en/latest/?badge)   [![](https://cdn.jsdelivr.net/gh/thespacedoctor/soxspipe@master/coverage.svg)](https://cdn.jsdelivr.net/gh/thespacedoctor/soxspipe@master/htmlcov/index.html)

*The data-reduction pipeline for the SOXS instrument* (a python package with command-line tools).

Documentation for soxspipe is hosted by [Read the Docs](http://soxspipe.readthedocs.org/en/stable/) (last
[stable version](http://soxspipe.readthedocs.org/en/stable/) and [latest version](http://soxspipe.readthedocs.org/en/latest/)). The code lives on [github](https://github.com/thespacedoctor/soxspipe). Please report any issues you find [here](https://github.com/thespacedoctor/soxspipe/issues).

## Installation

The easiest way to install soxspipe is to use `pip` (here we show the install inside of a conda environment):

``` bash
conda create -n soxspipe python=3.7 pip
conda activate soxspipe
pip install soxspipe
```

Or you can clone the [github repo](https://github.com/thespacedoctor/soxspipe) and install from a local version of the code:

``` bash
git clone git@github.com:thespacedoctor/soxspipe.git
cd soxspipe
python setup.py install
```

To upgrade to the latest version of soxspipe use the command:

``` bash
pip install soxspipe --upgrade
```

## Development

If you want to tinker with the code, then install in development mode. This means you can modify the code from your cloned repo:

``` bash
git clone git@github.com:thespacedoctor/soxspipe.git
cd soxspipe
python setup.py develop
```

[Pull requests](https://github.com/thespacedoctor/soxspipe/pulls) are welcomed! 


Note the data-sets are fairly bulky so make sure you have plenty of space on the drive you are downloading this data to.

<!-- ### Sublime Snippets

If you use [Sublime Text](https://www.sublimetext.com/) as your code editor, and you're planning to develop your own python code with soxspipe, you might find [my Sublime Snippets](https://github.com/thespacedoctor/soxspipe-Sublime-Snippets) useful. -->





