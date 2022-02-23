###################
uFJC Python Package
###################

|docs| |build| |coverage| |pyversions| |pypi| |anaconda| |license| |zenodo|

The is the Python package ``ufjc`` developed for the uFJC single-chain model, a freely-joined chain model with arbitrary link potential u.

************
Installation
************

The package can be installed using ``pip`` via the `Python Package Index (PyPI) <https://pypi.org/project/ufjc/>`_, which is recommended,

::

    pip install ufjc

or using ``conda`` via the ``mrbuche`` channel on `Anaconda <https://anaconda.org/mrbuche/ufjc>`_,

::

    conda install --channel mrbuche ufjc
    
Alternatively, a branch can be directly installed using ``pip``,

::

    pip install git+https://github.com/sandialabs/ufjc.git@<branch-name>

or after cloning a branch and executing ``python setup.py install``.

***********
Information
***********

- `Documentation <https://sandialabs.github.io/ufjc>`__
- `Examples <https://sandialabs.github.io/ufjc/examples>`__
- `Release History <https://github.com/sandialabs/ufjc/releases>`__
- `Tutorial <https://sandialabs.github.io/ufjc/Tutorial.html>`__

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

.. |pyversions| image:: https://img.shields.io/pypi/pyversions/ufjc.svg?logo=python&logoColor=FBE072
    :target: https://pypi.org/project/ufjc/

.. |pypi| image:: https://img.shields.io/pypi/v/ufjc?logo=pypi&logoColor=FBE072
    :target: https://pypi.org/project/ufjc/

.. |anaconda| image:: https://img.shields.io/conda/v/mrbuche/ufjc.svg?logo=anaconda
    :target: https://anaconda.org/mrbuche/ufjc/
    :alt: conda

.. |license| image:: https://img.shields.io/github/license/sandialabs/ufjc
    :target: https://github.com/sandialabs/ufjc/blob/main/LICENSE

.. |zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.6114263.svg
    :target: https://doi.org/10.5281/zenodo.6114263
