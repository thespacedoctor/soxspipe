# soxspipe



[![](https://zenodo.org/badge/DOI/10.5281/zenodo.8038264.svg)](https://zenodo.org/doi/10.5281/zenodo.8038264) 




<!-- INFO BADGES -->  

[![](https://img.shields.io/pypi/pyversions/soxspipe)](https://pypi.org/project/soxspipe/)
[![](https://img.shields.io/pypi/v/soxspipe)](https://pypi.org/project/soxspipe/)
[![](https://img.shields.io/conda/vn/conda-forge/soxspipe)](https://anaconda.org/conda-forge/soxspipe)
[![](https://static.pepy.tech/badge/soxspipe/month)](https://pepy.tech/project/soxspipe)
[![](https://img.shields.io/github/license/thespacedoctor/soxspipe)](https://github.com/thespacedoctor/soxspipe)

<!-- STATUS BADGES -->  

[![Required tests](https://github.com/thespacedoctor/soxspipe/actions/workflows/tests.yml/badge.svg?branch=develop)](https://github.com/thespacedoctor/soxspipe/actions/workflows/tests.yml)
[![](https://readthedocs.org/projects/soxspipe/badge/?version=main)](https://soxspipe.readthedocs.io/en/main/)
[![](https://img.shields.io/github/issues/thespacedoctor/soxspipe/type:%20bug?label=bug%20issues)](https://github.com/thespacedoctor/soxspipe/issues?q=is%3Aissue+is%3Aopen+label%3A%22type%3A+bug%22+)

*The data-reduction pipeline for the SOXS instrument* (a python package with command-line tools).

Documentation for soxspipe is hosted by [Read the Docs](https://soxspipe.readthedocs.io/en/main/) ([development version](https://soxspipe.readthedocs.io/en/develop/) and [main version](https://soxspipe.readthedocs.io/en/main/)). The code lives on [github](https://github.com/thespacedoctor/soxspipe). Please report any issues you find [here](https://github.com/thespacedoctor/soxspipe/issues).

## Installation

The best way to install or upgrade soxspipe is to use `conda` to install the package in its own isolated environment, as shown here:

``` bash
conda create -n soxspipe python=3.12 soxspipe -c conda-forge
conda activate soxspipe
```

If you have previously installed soxspipe, a warning will be issued stating that a conda environment already exists; select 'y' when asked to remove the existing environment.

To check installation was successful run `soxspipe -v`. This should return the version number of the install.

## Testing

Install the test dependencies and run the required offline suite with Python 3.12:

```bash
python -m pip install -e ".[tests]"
python -m pytest tests/unit tests/integration -m "not slow"
```

The real-data acceptance test is opt-in. The scheduled GitHub Actions workflow downloads the immutable archive named by `tests/real_data/manifest.json`, verifies it, creates a disposable workspace, and runs the representative NIR-offset reduction. To run the acceptance assertions locally after preparing that workspace, set `SOXSPIPE_REAL_DATA_DIR` to its absolute path:

```bash
SOXSPIPE_REAL_DATA_DIR=/absolute/path/to/workspace python -m pytest tests/real_data
```

## Changed-line gates

Two gates run on every pull request into `develop`, and both judge only the lines the pull request changed. Neither one asks the whole package to be clean, so existing debt never blocks an unrelated change:

- **Coverage.** `diff-cover reports/coverage.xml --compare-branch=origin/develop --fail-under=80` requires 80% coverage of changed lines.
- **Lint.** `python tools/lint_ratchet.py --compare-branch origin/develop` fails when a ruff finding lands on a changed line. Ruff has no baseline feature, so the tool intersects `ruff check --output-format json` with the diff hunks itself. Pre-existing findings in the files you touched are counted and printed as context, never gated; the package carries roughly 1,380 of them in total, and `ruff check soxspipe` reports that whole-package figure.

A finding is matched on the line ruff anchors it to. A finding that covers many lines, such as `PLR0915` for an over-long function, is anchored at the `def` line: writing a new over-long function fails the gate, while adding a statement to one that is already over-long does not.

Run the lint gate locally against your staged changes before committing:

```bash
python tools/lint_ratchet.py --staged
```

To run it automatically on every commit, install the hook. It is a convenience, not the gate — `git commit --no-verify` skips it, and the CI step does not:

```bash
python -m pip install -e ".[dev]"
pre-commit install
```

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
