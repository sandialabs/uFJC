###################
uFJC Python Package
###################

|docs| |build| |coverage| |pyversions| |pypi| |conda| |docker| |license| |zenodo|

The is the Python package ``ufjc`` developed for the uFJC single-chain model, a freely-joined chain model with arbitrary link potential u.

************
Installation
************

The package can be installed using ``pip`` via the `Python Package Index (PyPI) <https://pypi.org/project/ufjc/>`_,

::

    pip install ufjc

or using ``conda`` via the ``mrbuche`` channel on `Anaconda <https://anaconda.org/mrbuche/ufjc>`_,

::

    conda install --channel mrbuche ufjc
    
Alternatively, a branch can be directly installed using

::

    pip install git+https://github.com/sandialabs/ufjc.git@<branch-name>

or after cloning a branch and executing ``python setup.py install``.

***********
Information
***********

- `Contributing <https://sandialabs.github.io/ufjc/CONTRIBUTING.html>`__
- `Documentation <https://sandialabs.github.io/ufjc>`__
- `Examples <https://sandialabs.github.io/ufjc/examples>`__
- `License <https://github.com/sandialabs/ufjc/blob/main/LICENSE>`__
- `Releases <https://github.com/sandialabs/ufjc/releases>`__
- `Tutorial <https://sandialabs.github.io/ufjc/TUTORIAL.html>`__

********
Citation
********

\M. R. Buche and S. J. Grutzik, ``ufjc``: the Python package for the uFJC single-chain model, `Zenodo (2022) <https://doi.org/10.5281/zenodo.6114263>`_.

*********
Copyright
*********

Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

..
    Badges ========================================================================

.. |docs| image:: https://github.com/sandialabs/ufjc/actions/workflows/docs.yml/badge.svg
    :target: https://sandialabs.github.io/ufjc

.. |build| image:: https://github.com/sandialabs/ufjc/workflows/main/badge.svg
    :target: https://github.com/sandialabs/ufjc

.. |coverage| image:: https://coveralls.io/repos/github/sandialabs/ufjc/badge.svg?branch=main
    :target: https://coveralls.io/github/sandialabs/ufjc?branch=main

.. |pyversions| image:: https://img.shields.io/pypi/pyversions/ufjc.svg?logo=python&logoColor=FBE072&color=4B8BBE&label=Python
    :target: https://pypi.org/project/ufjc/

.. |pypi| image:: https://img.shields.io/pypi/v/ufjc?logo=pypi&logoColor=FBE072&label=PyPI&color=4B8BBE
    :target: https://pypi.org/project/ufjc/

.. |conda| image:: https://img.shields.io/conda/v/mrbuche/ufjc.svg?logo=anaconda&color=3EB049&label=Anaconda
    :target: https://anaconda.org/mrbuche/ufjc/
    :alt: conda

.. |docker| image:: https://img.shields.io/docker/v/mrbuche/ufjc?color=0db7ed&label=Docker%20Hub&logo=docker&logoColor=0db7ed
    :target: https://hub.docker.com/r/mrbuche/ufjc
    :alt: docker

.. |license| image:: https://img.shields.io/github/license/sandialabs/ufjc
    :target: https://github.com/sandialabs/ufjc/blob/main/LICENSE

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6114263.svg
    :target: https://doi.org/10.5281/zenodo.6114263
