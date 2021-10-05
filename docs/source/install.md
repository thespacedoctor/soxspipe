# Installation

The best way to install soxspipe is to use `conda` and install the package in its own isolated environment, as shown here:

``` bash
conda create -n soxspipe python=3.8 soxspipe
conda activate soxspipe
```

To check installation was successful run `soxspipe -v`. This should return the version number of the install.

To upgrade to the latest version of soxspipe use the command:

``` bash
conda upgrade soxspipe
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


