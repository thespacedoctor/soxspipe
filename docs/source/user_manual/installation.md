# Installation

There are multiple methods for installing soxspipe, and you can choose the method that best suits your needs.

::::{tab-set}

:::{tab-item} Miniforge

[Miniforge](https://github.com/conda-forge/miniforge) provides a minimal installation for the conda package and environment manager. It uses the community-maintained [conda-forge](https://conda-forge.org) distribution channel by default and, therefore, completely avoids any commercial licensing issues. To install soxspipe via Miniforge, you will first need an up-to-date installation on your machine. Please refer to [the Miniforge documentation](https://github.com/conda-forge/miniforge) for installation instructions. 

With Miniforge initialised, generate an isolated environment for soxspipe and install it from the conda-forge distribution channel with the following command:

```bash
conda create -n soxspipe python=3.9 soxspipe -c conda-forge
```

If you have previously installed soxspipe, a warning will be issued stating that a conda environment already exists; select `y` when asked to remove the existing environment. This has proven to be the cleanest way to upgrade soxspipe.

Once the installation completes, activate the soxspipe conda environment and check soxspipe is installed.

```bash
conda activate soxspipe
soxspipe -v
```
:::


:::{tab-item} Anaconda

[Anaconda](https://docs.anaconda.com/anaconda/install/index.html) and its lighter-weight alternative [Minicoda](https://docs.conda.io/en/latest/miniconda.html), are distributions for Python and other programming languages. To install soxspipe via Anaconda or Miniconda, you will first need an up-to-date installation on your machine. Please refer to their documentation for installation instructions. 

With Anaconda/Mininconda initialised, generate an isolated environment for soxspipe and install it from the conda-forge distribution channel with the following command:

```bash
conda create -n soxspipe python=3.9 soxspipe -c conda-forge
```

If you have previously installed soxspipe, a warning will be issued stating that a conda environment already exists; select `y` when asked to remove the existing environment. This has proven to be the cleanest way to upgrade soxspipe.

Once the installation completes, activate the soxspipe conda environment and check soxspipe is installed.

```bash
conda activate soxspipe
soxspipe -v
```
:::

:::{tab-item} Pip

You can install soxspipe via pip within a virtual environment (recommended) or directly into your machine's native Python environment.

To create and activate a virtual environment in your home folder, run the following:

```bash
cd ~/
mkdir soxspipe
cd soxspipe
python -m venv venv
source venv/bin/activate
```

You can now install soxspipe via pip:

```bash
pip install soxspipe
```

Finally, check soxspipe is installed correctly by checking the version number:

```bash
soxspipe -v
```

Each time you wish to use soxspipe, remember to first activate the virtual environment with:

```bash
source ~/soxspipe/venv/bin/activate
```
:::


:::{tab-item} Github

You can download the source code for soxspipe directly from Github. This is the best installation method if you plan to change/bug-fix the source code (pull requests are welcome!). First, clone the repository somewhere on your machine:

```bash
cd ~/
git clone https://github.com/thespacedoctor/soxspipe.git
cd soxspipe
```

You can install source code via pip within a virtual environment (recommended) or directly into your machine's native Python environment.

To create and activate a virtual environment in your home folder, run the following:

```bash
python -m venv venv
source venv/bin/activate
```

You can now install soxspipe via pip:

```bash
pip install -e .
```

Finally, check soxspipe is installed correctly by checking the version number:

```bash
soxspipe -v
```

Each time you wish to use soxspipe, remember to first activate the virtual environment with:

```bash
source ~/soxspipe/venv/bin/activate
```
:::


::::
