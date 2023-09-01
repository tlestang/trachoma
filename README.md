# NTDMC trachoma model

This repository holds a Python package for the simulation of the
spread of trachoma over a population of individuals.  This is an
individuals-based model, meaning that the spread of trachoma is
simulated by tracking the state of every single individuals in a
finite population, over time. For a description of the underlying
mathematical model and transmission rules, see [Pinsent A,
Hollingsworth TD
(2018)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0006531).

## Installation

The package is currently in early development. No official package
distribution is available on either PyPI or Anaconda cloud yet, but it
shouldn't be too hard to install the package from source.

```
git clone https://github.com:NTDMC-Modelling-consortium/trachoma
```

In addition to Python (3.9+), you'll need a C compiler, for instance
gcc or clang. You can then install the model package using `pip`

``` shell
$ cd trachoma
$ python -m pip install .
```

The package only depends on NumPy, but it is still recommended to use
a virtual environment to isolate this verion on NumPy from others you
might have installed on your machine.

### Updating

The package is currently changing on a daily basis. To update, your
local copy of the repository, run

```
git pull origin main
```

from the repo root directory.

In you installed the package in editable mode (see next section), then
no more actions are required. If not, you'll need to reinstall the
Python package into your environment:

```
$ python -m pip install .
```

## Contributing

Start by install the package in [editable
mode](https://packaging.python.org/en/latest/guides/distributing-packages-using-setuptools/#working-in-development-mode):

``` shell
$ python -m pip install -e .
```

### Running the tests

Currently there are only tests covering the C core functions in
`tests/C`. You can run these tests using
[`pytest`](https://docs.pytest.org/):

```
$ python -m pip install pytest
$ pytest tests/C
```

### Building the documentation website

You'll need [Sphinx](https://www.sphinx-doc.org/en/master/) and
[Hawkmoth](https://hawkmoth.readthedocs.io/en/stable/). The latter is
used to automatically populate the documentation source with
documentation comments written directly in the C source files. it uses
the [Clang](https://clang.llvm.org/) python interface to parse C
source code and extract documentaion comments.

```
$ python -m pip install sphinx hawkmoth libclang
```
